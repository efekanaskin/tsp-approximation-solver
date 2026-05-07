
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <random>
#include <stdexcept>
#include <limits>
#include <algorithm>
#include <chrono>
#include <iomanip>   // nicer printing

// -------------------- DATA STRUCTURES --------------------
struct City { double x, y; };

using DistMatrix  = std::vector<std::vector<double>>;
using PheroMatrix = std::vector<std::vector<double>>;

// -------------------- LOAD INSTANCE --------------------
std::vector<City> loadInstance(const std::string& filename) {
    std::ifstream in(filename);
    if (!in) throw std::runtime_error("Cannot open file: " + filename);

    int n;
    in >> n;
    if (n <= 0) throw std::runtime_error("Invalid number of cities in file.");

    std::vector<City> cities(n);
    for (int i = 0; i < n; i++) in >> cities[i].x >> cities[i].y;
    return cities;
}

// -------------------- DISTANCE MATRIX --------------------
double euclidean(const City& a, const City& b) {
    double dx = a.x - b.x, dy = a.y - b.y;
    return std::sqrt(dx*dx + dy*dy);
}

DistMatrix computeDistanceMatrix(const std::vector<City>& cities) {
    int n = (int)cities.size();
    DistMatrix D(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double d = euclidean(cities[i], cities[j]);
            D[i][j] = d;
            D[j][i] = d;
        }
    }
    return D;
}

// -------------------- VALIDATION --------------------
bool isValidTour(const std::vector<int>& tour, int n) {
    if ((int)tour.size() != n) return false;
    std::vector<int> seen(n, 0);
    for (int v : tour) {
        if (v < 0 || v >= n) return false;
        seen[v]++;
        if (seen[v] > 1) return false;
    }
    return true;
}

// -------------------- TOUR LENGTH --------------------
double tourLength(const std::vector<int>& tour, const DistMatrix& D) {
    double length = 0.0;
    int n = (int)tour.size();
    for (int i = 0; i < n - 1; i++) length += D[tour[i]][tour[i + 1]];
    length += D[tour[n - 1]][tour[0]];
    return length;
}

// =======================================================
//                    ACO  (Algorithm #1)
// =======================================================
PheroMatrix initializePheromones(int n, double tau0) {
    PheroMatrix tau(n, std::vector<double>(n, tau0));
    for (int i = 0; i < n; i++) tau[i][i] = 0.0;
    return tau;
}

int chooseNextCity(
    int current,
    const std::vector<bool>& visited,
    const DistMatrix& D,
    const PheroMatrix& tau,
    double alpha,
    double beta,
    std::mt19937& rng
) {
    int n = (int)visited.size();
    std::vector<double> weight(n, 0.0);
    double sum = 0.0;

    for (int j = 0; j < n; j++) {
        if (!visited[j]) {
            double dist = D[current][j];
            double eta = (dist > 0.0) ? (1.0 / dist) : 1e9;
            double w = std::pow(tau[current][j], alpha) * std::pow(eta, beta);
            if (!std::isfinite(w) || w < 0.0) w = 0.0;
            weight[j] = w;
            sum += w;
        }
    }

    if (sum <= 0.0) {
        for (int j = 0; j < n; j++) if (!visited[j]) return j;
        return -1;
    }

    std::uniform_real_distribution<double> dist(0.0, sum);
    double r = dist(rng);

    double cumulative = 0.0;
    for (int j = 0; j < n; j++) {
        if (!visited[j]) {
            cumulative += weight[j];
            if (cumulative >= r) return j;
        }
    }

    for (int j = 0; j < n; j++) if (!visited[j]) return j;
    return -1;
}

std::vector<int> buildAntTour(
    int startCity,
    const DistMatrix& D,
    const PheroMatrix& tau,
    double alpha,
    double beta,
    std::mt19937& rng
) {
    int n = (int)D.size();
    std::vector<bool> visited(n, false);
    std::vector<int> tour;
    tour.reserve(n);

    int current = startCity;
    visited[current] = true;
    tour.push_back(current);

    while ((int)tour.size() < n) {
        int next = chooseNextCity(current, visited, D, tau, alpha, beta, rng);
        if (next == -1) break;
        visited[next] = true;
        tour.push_back(next);
        current = next;
    }

    if ((int)tour.size() < n) {
        std::vector<int> remaining;
        for (int j = 0; j < n; j++) if (!visited[j]) remaining.push_back(j);
        std::shuffle(remaining.begin(), remaining.end(), rng);
        for (int v : remaining) tour.push_back(v);
    }

    return tour;
}

void evaporate(PheroMatrix& tau, double rho) {
    double factor = 1.0 - rho;
    const double minTau = 1e-12;
    int n = (int)tau.size();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            tau[i][j] *= factor;
            if (tau[i][j] < minTau) tau[i][j] = minTau;
        }
        tau[i][i] = 0.0;
    }
}

void deposit(PheroMatrix& tau, const std::vector<int>& tour, double L, double Q) {
    if (L <= 0.0 || !std::isfinite(L)) return;
    double add = Q / L;
    int n = (int)tour.size();

    for (int i = 0; i < n - 1; i++) {
        int a = tour[i], b = tour[i + 1];
        tau[a][b] += add; tau[b][a] += add;
    }
    int last = tour[n - 1], first = tour[0];
    tau[last][first] += add; tau[first][last] += add;

    for (int i = 0; i < n; i++) tau[i][i] = 0.0;
}

struct ACOResult {
    double bestLen;
    std::vector<int> bestTour;
};

ACOResult runACO(const DistMatrix& D, int iterations, int numAnts,
                 double tau0, double alpha, double beta, double rho, double Q,
                 std::mt19937& rng) {
    int n = (int)D.size();
    auto pheromone = initializePheromones(n, tau0);

    double globalBestLen = std::numeric_limits<double>::infinity();
    std::vector<int> globalBestTour;

    std::uniform_int_distribution<int> startDist(0, n - 1);

    for (int it = 1; it <= iterations; it++) {
        double bestLenIter = std::numeric_limits<double>::infinity();
        std::vector<int> bestTourIter;

        for (int k = 0; k < numAnts; k++) {
            int startCity = startDist(rng);
            auto tour = buildAntTour(startCity, D, pheromone, alpha, beta, rng);
            if (!isValidTour(tour, n)) continue;

            double len = tourLength(tour, D);
            if (len < bestLenIter) { bestLenIter = len; bestTourIter = tour; }
        }

        if (!bestTourIter.empty() && bestLenIter < globalBestLen) {
            globalBestLen = bestLenIter;
            globalBestTour = bestTourIter;
        }

        evaporate(pheromone, rho);
        if (!bestTourIter.empty() && std::isfinite(bestLenIter)) {
            deposit(pheromone, bestTourIter, bestLenIter, Q);
        }
    }

    return {globalBestLen, globalBestTour};
}

// =======================================================
//             Simulated Annealing (Algorithm #2)
// =======================================================
std::vector<int> createInitialTour(int n) {
    std::vector<int> tour(n);
    for (int i = 0; i < n; i++) tour[i] = i;
    return tour;
}

void applyTwoOpt(std::vector<int>& tour, int i, int j) {
    std::reverse(tour.begin() + i, tour.begin() + j + 1);
}

struct SAResult {
    double bestLen;
    std::vector<int> bestTour;
};

SAResult runSA(const DistMatrix& D,
               int iterations,
               double T0,
               double Tmin,
               double coolingRate,
               std::mt19937& rng) {
    int n = (int)D.size();

    std::vector<int> current = createInitialTour(n);
    std::shuffle(current.begin(), current.end(), rng);

    double currentLen = tourLength(current, D);
    std::vector<int> best = current;
    double bestLen = currentLen;

    std::uniform_int_distribution<int> idxDist(0, n - 1);
    std::uniform_real_distribution<double> uni(0.0, 1.0);

    double T = T0;

    for (int it = 0; it < iterations && T > Tmin; it++) {
        int i = idxDist(rng);
        int j = idxDist(rng);
        if (i > j) std::swap(i, j);

        if (i == j) {
            T *= coolingRate;
            continue;
        }

        std::vector<int> candidate = current;
        applyTwoOpt(candidate, i, j);

        double candLen = tourLength(candidate, D);
        double delta = candLen - currentLen;

        if (delta <= 0.0 || uni(rng) < std::exp(-delta / T)) {
            current = std::move(candidate);
            currentLen = candLen;

            if (currentLen < bestLen) {
                bestLen = currentLen;
                best = current;
            }
        }

        T *= coolingRate;
    }

    return {bestLen, best};
}

// -------------------- STATS HELPERS --------------------
struct Stats {
    double best = std::numeric_limits<double>::infinity();
    double mean = 0.0;
    double stddev = 0.0;
};

Stats computeStats(const std::vector<double>& values) {
    Stats s;
    if (values.empty()) return s;

    double sum = 0.0;
    for (double v : values) {
        if (v < s.best) s.best = v;
        sum += v;
    }
    s.mean = sum / (double)values.size();

    double var = 0.0;
    for (double v : values) {
        double diff = v - s.mean;
        var += diff * diff;
    }
    var /= (double)values.size();
    s.stddev = std::sqrt(var);
    return s;
}

// =======================================================
//     EXACT OPTIMAL (ONLY FOR n <= 12)  - Branch & Bound
// =======================================================
struct ExactResult {
    double optLen = std::numeric_limits<double>::infinity();
    std::vector<int> optTour;
};

void dfsExactTSP(
    int start,
    int current,
    int depth,
    double currentLen,
    const DistMatrix& D,
    std::vector<int>& path,
    std::vector<bool>& used,
    ExactResult& best
) {
    int n = (int)D.size();
    if (currentLen >= best.optLen) return; // prune

    if (depth == n) {
        double total = currentLen + D[current][start];
        if (total < best.optLen) {
            best.optLen = total;
            best.optTour = path;
        }
        return;
    }

    for (int nxt = 0; nxt < n; nxt++) {
        if (!used[nxt]) {
            used[nxt] = true;
            path.push_back(nxt);

            dfsExactTSP(start, nxt, depth + 1, currentLen + D[current][nxt], D, path, used, best);

            path.pop_back();
            used[nxt] = false;
        }
    }
}

ExactResult computeExactOptimalIfSmall(const DistMatrix& D) {
    int n = (int)D.size();
    ExactResult res;
    if (n > 12) return res; // INF = "not computed"

    int start = 0;
    std::vector<bool> used(n, false);
    std::vector<int> path;
    path.reserve(n);

    used[start] = true;
    path.push_back(start);

    dfsExactTSP(start, start, 1, 0.0, D, path, used, res);
    return res;
}

// -------------------- WINNER HELPER --------------------
std::string winnerLabel(double acoLen, double saLen) {
    const double eps = 1e-9;
    if (acoLen + eps < saLen) return "ACO";
    if (saLen + eps < acoLen) return "SA";
    return "Tie";
}

// -------------------- CSV WRITER --------------------
void writeCSV(
    const std::string& outName,
    int n,
    bool hasOpt,
    double optLen,
    const std::vector<double>& acoLen,
    const std::vector<double>& saLen,
    const std::vector<double>& acoMs,
    const std::vector<double>& saMs,
    const std::vector<double>& acoRatio,
    const std::vector<double>& saRatio
) {
    std::ofstream out(outName);
    if (!out) {
        std::cerr << "Warning: could not write CSV: " << outName << "\n";
        return;
    }

    out << "n,has_opt,opt_len,run,aco_len,sa_len,aco_ms,sa_ms,aco_ratio,sa_ratio\n";
    for (size_t i = 0; i < acoLen.size(); i++) {
        out << n << "," << (hasOpt ? 1 : 0) << ",";
        if (hasOpt) out << optLen; else out << "";
        out << "," << (i + 1) << ","
            << acoLen[i] << "," << saLen[i] << ","
            << acoMs[i] << "," << saMs[i] << ",";
        if (hasOpt) {
            out << acoRatio[i] << "," << saRatio[i];
        } else {
            out << ","; // aco_ratio empty
            out << "";  // sa_ratio empty
        }
        out << "\n";
    }
}

// -------------------- MAIN --------------------
int main(int argc, char* argv[]) {
    try {
        if (argc < 2) {
            std::cout << "Usage: ./tsp data/instance_small_12.txt\n";
            return 1;
        }

        std::string filename = argv[1];
        auto cities = loadInstance(filename);
        auto D = computeDistanceMatrix(cities);
        int n = (int)cities.size();

        // ---- ACO parameters ----
        int acoIterations = 50;
        int numAnts = 50;
        double tau0 = 1.0;
        double alpha = 1.0;
        double beta  = 2.0;
        double rho   = 0.5;
        double Q     = 100.0;

        // ---- SA parameters ----
        int saIterations = 20000;
        double T0 = 1000.0;
        double Tmin = 1e-6;
        double coolingRate = 0.9995;

        // ---- EXPERIMENT SETTINGS ----
        int runs = 10;
        bool printEachRun = true;

        std::random_device rd;

        std::vector<double> acoLens, saLens;
        std::vector<double> acoTimesMs, saTimesMs;
        std::vector<double> acoRatios, saRatios; // only if OPT computed

        acoLens.reserve(runs);
        saLens.reserve(runs);
        acoTimesMs.reserve(runs);
        saTimesMs.reserve(runs);
        acoRatios.reserve(runs);
        saRatios.reserve(runs);

        std::cout << std::fixed << std::setprecision(6);

        std::cout << "Instance: " << filename << "\n";
        std::cout << "Cities: " << n << "\n";
        std::cout << "Runs: " << runs << "\n\n";

        // Exact optimal only if n <= 12
        bool hasOpt = false;
        double optLen = std::numeric_limits<double>::infinity();
        if (n <= 12) {
            std::cout << "Computing EXACT optimal for n <= 12 ...\n";
            auto exact = computeExactOptimalIfSmall(D);
            optLen = exact.optLen;
            if (!std::isfinite(optLen)) {
                throw std::runtime_error("Exact optimal failed (optLen is INF).");
            }
            hasOpt = true;
            std::cout << "OPT length (exact): " << optLen << "\n\n";
        } else {
            std::cout << "Exact optimal is skipped (n > 12).\n\n";
        }

        using Clock = std::chrono::high_resolution_clock;

        for (int r = 1; r <= runs; r++) {
            unsigned seed = (unsigned)rd();
            std::mt19937 rngACO(seed);
            std::mt19937 rngSA(seed + 1u);

            auto t1 = Clock::now();
            auto aco = runACO(D, acoIterations, numAnts, tau0, alpha, beta, rho, Q, rngACO);
            auto t2 = Clock::now();
            double acoMs = std::chrono::duration<double, std::milli>(t2 - t1).count();

            auto t3 = Clock::now();
            auto sa  = runSA(D, saIterations, T0, Tmin, coolingRate, rngSA);
            auto t4 = Clock::now();
            double saMs = std::chrono::duration<double, std::milli>(t4 - t3).count();

            if (aco.bestTour.empty() || !isValidTour(aco.bestTour, n)) {
                throw std::runtime_error("ACO produced an invalid tour in run " + std::to_string(r));
            }
            if (sa.bestTour.empty() || !isValidTour(sa.bestTour, n)) {
                throw std::runtime_error("SA produced an invalid tour in run " + std::to_string(r));
            }

            acoLens.push_back(aco.bestLen);
            saLens.push_back(sa.bestLen);
            acoTimesMs.push_back(acoMs);
            saTimesMs.push_back(saMs);

            double acoRatio = 0.0, saRatio = 0.0;
            if (hasOpt && optLen > 0.0) {
                acoRatio = aco.bestLen / optLen;
                saRatio  = sa.bestLen  / optLen;
                acoRatios.push_back(acoRatio);
                saRatios.push_back(saRatio);
            }

            if (printEachRun) {
                std::cout << "Run " << r << ":\n";
                std::cout << "  ACO best length: " << aco.bestLen << "\n";
                std::cout << "  SA  best length: " << sa.bestLen  << "\n";
                std::cout << "  ACO time (ms):   " << acoMs << "\n";
                std::cout << "  SA  time (ms):   " << saMs << "\n";
                if (hasOpt) {
                    std::cout << "  ACO ratio vs OPT: " << acoRatio << "\n";
                    std::cout << "  SA  ratio vs OPT: " << saRatio  << "\n";
                }
                std::cout << "  Winner (length): " << winnerLabel(aco.bestLen, sa.bestLen) << "\n\n";
            }
        }

        auto acoLenStats = computeStats(acoLens);
        auto saLenStats  = computeStats(saLens);
        auto acoTimeStats = computeStats(acoTimesMs);
        auto saTimeStats  = computeStats(saTimesMs);

        std::cout << "==================== SUMMARY (Report Friendly) ====================\n";
        std::cout << "Length:\n";
        std::cout << "  ACO  best=" << acoLenStats.best
                  << "  avg=" << acoLenStats.mean
                  << "  std=" << acoLenStats.stddev << "\n";
        std::cout << "  SA   best=" << saLenStats.best
                  << "  avg=" << saLenStats.mean
                  << "  std=" << saLenStats.stddev << "\n";

        std::cout << "Time (ms):\n";
        std::cout << "  ACO  best=" << acoTimeStats.best
                  << "  avg=" << acoTimeStats.mean
                  << "  std=" << acoTimeStats.stddev << "\n";
        std::cout << "  SA   best=" << saTimeStats.best
                  << "  avg=" << saTimeStats.mean
                  << "  std=" << saTimeStats.stddev << "\n";

        if (hasOpt) {
            auto acoRatioStats = computeStats(acoRatios);
            auto saRatioStats  = computeStats(saRatios);
            std::cout << "Optimal:\n";
            std::cout << "  OPT length=" << optLen << "\n";
            std::cout << "Ratio (len/OPT):\n";
            std::cout << "  ACO  best=" << acoRatioStats.best
                      << "  avg=" << acoRatioStats.mean
                      << "  std=" << acoRatioStats.stddev << "\n";
            std::cout << "  SA   best=" << saRatioStats.best
                      << "  avg=" << saRatioStats.mean
                      << "  std=" << saRatioStats.stddev << "\n";
        } else {
            std::cout << "Ratios vs OPT: skipped (n > 12)\n";
        }

        // Write CSV
        std::string csvName = "results_" + std::to_string(n) + ".csv";
        writeCSV(csvName, n, hasOpt, hasOpt ? optLen : 0.0,
                 acoLens, saLens, acoTimesMs, saTimesMs, acoRatios, saRatios);
        std::cout << "\nSaved CSV: " << csvName << "\n";

        return 0;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
