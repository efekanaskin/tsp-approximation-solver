# TSP Metaheuristics: ACO vs. Simulated Annealing 🗺️

## Project Overview
This repository contains a comparative study of two metaheuristic algorithms—Ant Colony Optimization (ACO) and Simulated Annealing (SA)—implemented to solve the metric Traveling Salesman Problem (TSP). This project explores the practical trade-offs between exploration, solution quality, and scalability in heuristic optimization.

The TSP is a classical NP-hard problem. Exact algorithms become computationally infeasible as the number of locations increases. To address real-world routing challenges in logistics, we implemented these nature-inspired algorithms to find high-quality approximations within limited time constraints.

## Algorithms Implemented

### 1. Ant Colony Optimization (ACO)
* Inspired by the foraging behavior of ants communicating via pheromone trails.
* Utilizes a pheromone matrix that updates over iterations, favoring shorter paths alongside a heuristic value (the inverse of Euclidean distance).
* Employs pheromone evaporation to prevent premature convergence.

### 2. Simulated Annealing (SA)
* Inspired by the physical process of annealing in metallurgy.
* Starts with a random feasible tour and explores neighboring solutions using modifications like the 2-opt swap.
* Uses a cooling schedule and an acceptance probability function to occasionally accept worse solutions early on, effectively escaping local minima.

## Experimental Setup
The algorithms were evaluated on synthetic 2D-plane instances using Euclidean distance:
* **Small Instance (12 cities):** Used to compare heuristic solutions against the exact optimal tour found via exhaustive search.
* **Medium (50 cities) & Large (100 cities) Instances:** Used to analyze scalability, runtime behavior, and convergence without absolute optimal baselines.

## Key Results & Findings
Both algorithms successfully produced high-quality solutions, but their scalability varied significantly:
* **Small Instances:** Both ACO and SA consistently found the optimal solution (approximation ratio of 1) across all runs.
* **Scalability & Runtime:** As instance size increased to 100 cities, Simulated Annealing maintained competitive tour lengths while achieving an average runtime of **44 ms**, compared to ACO's **1215 ms**.
* **Conclusion:** Simulated Annealing demonstrated superior computational efficiency and scalability, whereas Ant Colony Optimization incurred much higher overhead due to its population-based design and pheromone update mechanisms.

