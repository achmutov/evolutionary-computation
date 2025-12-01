#include <evolutionary_computation/solver/lns.h>
#include <evolutionary_computation/solver/random.h>
#include <algorithm>
#include <chrono>
#include <random>
#include <limits>
#include <numeric>

LNSSolver::LNSSolver(LocalSearchSolver& baseLS,
                     Solver::Duration targetTime,
                     DestroyType destroyType,
                     bool useLocalSearchAfterRepair,
                     double destroyFraction)
    : baseLS{baseLS},
      targetTime{targetTime},
      destroyType{destroyType},
      useLocalSearchAfterRepair{useLocalSearchAfterRepair},
      destroyFraction{destroyFraction},
      iterations{0} {}

std::string LNSSolver::name() const {
    std::string destroyName;
    switch (destroyType) {
        case DestroyType::Random:
            destroyName = "random";
            break;
        case DestroyType::SingleSubpath:
            destroyName = "single_subpath";
            break;
        case DestroyType::MultipleSubpaths:
            destroyName = "multiple_subpaths";
            break;
        case DestroyType::Heuristic:
            destroyName = "heuristic";
            break;
    }
    
    std::string lsSuffix = useLocalSearchAfterRepair ? "" : "-noLS";
    return std::string("lns-") + destroyName + lsSuffix;
}

void LNSSolver::init(Data const& data) {
    Solver::init(data);
    baseLS.init(data);
    iterations = 0;
}

Solver::Indices LNSSolver::_solve(int i) {
    iterations = 0;
    
    // 1. Generate random initial solution
    RandomSolver randomSolver;
    randomSolver.init(data);
    Indices x = randomSolver._solve(i);
    
    // 2. Apply LS to initial solution (always)
    x = baseLS._solveFromSolution(x);
    
    int xCost = this->cost(x);
    auto startTime = std::chrono::high_resolution_clock::now();
    
    // 3. Main LNS loop
    while (true) {
        auto currentTime = std::chrono::high_resolution_clock::now();
        auto elapsed = currentTime - startTime;
        
        if (elapsed >= targetTime) {
            break;  // Stop when we've used target time
        }
        
        // Destroy: remove nodes from solution
        Indices partialSolution = destroy(x);
        
        // Repair: insert nodes from all available nodes to reach target count
        Indices y = repair(partialSolution);
        
        // Optional local search after repair
        if (useLocalSearchAfterRepair) {
            y = baseLS._solveFromSolution(y);
        }
        
        int yCost = this->cost(y);
        
        // Accept if better (greedy acceptance)
        if (yCost < xCost) {
            x = y;
            xCost = yCost;
        }
        
        iterations++;
    }
    
    return x;
}

Solver::Indices LNSSolver::destroy(Indices const& solution) {
    switch (destroyType) {
        case DestroyType::Random:
            return destroyRandom(solution);
        case DestroyType::SingleSubpath:
            return destroySingleSubpath(solution);
        case DestroyType::MultipleSubpaths:
            return destroyMultipleSubpaths(solution);
        case DestroyType::Heuristic:
            return destroyHeuristic(solution);
    }
    return solution;
}

Solver::Indices LNSSolver::destroyRandom(Indices const& solution) {
    int n = solution.size();
    int numToRemove = static_cast<int>(std::round(n * destroyFraction));
    
    // Create a copy and shuffle indices
    std::vector<int> indices(n);
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), this->mt);
    
    // Select nodes to remove
    std::set<int> removedSet;
    for (int i = 0; i < numToRemove; i++) {
        removedSet.insert(solution[indices[i]]);
    }
    
    // Build partial solution (nodes not removed)
    Indices partialSolution;
    for (int node : solution) {
        if (removedSet.find(node) == removedSet.end()) {
            partialSolution.push_back(node);
        }
    }
    
    return partialSolution;
}

Solver::Indices LNSSolver::destroySingleSubpath(Indices const& solution) {
    int n = solution.size();
    int numToRemove = static_cast<int>(std::round(n * destroyFraction));
    
    // Pick random starting position
    int startIdx = this->mt() % n;
    
    std::set<int> removedSet;
    
    // Remove contiguous subpath
    for (int i = 0; i < numToRemove; i++) {
        int idx = (startIdx + i) % n;
        removedSet.insert(solution[idx]);
    }
    
    // Build partial solution (nodes not removed)
    Indices partialSolution;
    for (int node : solution) {
        if (removedSet.find(node) == removedSet.end()) {
            partialSolution.push_back(node);
        }
    }
    
    return partialSolution;
}

Solver::Indices LNSSolver::destroyMultipleSubpaths(Indices const& solution) {
    int n = solution.size();
    int numToRemove = static_cast<int>(std::round(n * destroyFraction));
    
    // Random number of subpaths: 5-10
    int numSubpaths = 5 + (this->mt() % 6);  // 5 to 10
    int nodesPerSubpath = numToRemove / numSubpaths;
    int remainder = numToRemove % numSubpaths;
    
    std::set<int> removedSet;
    std::vector<bool> removed(n, false);
    
    // Remove multiple subpaths
    for (int s = 0; s < numSubpaths; s++) {
        // Pick random starting position (not already removed)
        int startIdx;
        int attempts = 0;
        do {
            startIdx = this->mt() % n;
            attempts++;
        } while (removed[startIdx] && attempts < 100);
        
        if (attempts >= 100) {
            // Fallback: find any non-removed position
            for (int i = 0; i < n; i++) {
                if (!removed[i]) {
                    startIdx = i;
                    break;
                }
            }
        }
        
        // Determine length of this subpath
        int subpathLength = nodesPerSubpath;
        if (s < remainder) {
            subpathLength++;  // Distribute remainder
        }
        
        // Remove contiguous subpath (checking for wrap-around)
        for (int i = 0; i < subpathLength; i++) {
            int idx = (startIdx + i) % n;
            if (!removed[idx]) {
                removed[idx] = true;
                removedSet.insert(solution[idx]);
            }
        }
    }
    
    // Build partial solution (nodes not removed)
    Indices partialSolution;
    for (int node : solution) {
        if (removedSet.find(node) == removedSet.end()) {
            partialSolution.push_back(node);
        }
    }
    
    return partialSolution;
}

Solver::Indices LNSSolver::destroyHeuristic(Indices const& solution) {
    int n = solution.size();
    int numToRemove = static_cast<int>(std::round(n * destroyFraction));
    
    // Compute "badness" score for each node in solution
    std::vector<double> scores(n);
    std::vector<int> nodeIndices(n);  // Map from solution index to actual node
    
    for (int i = 0; i < n; i++) {
        int node = solution[i];
        nodeIndices[i] = node;
        
        // Get incident edges
        int prevNode = solution[(i - 1 + n) % n];
        int nextNode = solution[(i + 1) % n];
        
        // Compute detour cost: if we remove this node, how much longer does the path become?
        // Current: prev -> node -> next (with node cost)
        int currentCost = this->distances[prevNode][node] 
                         + this->distances[node][nextNode] 
                         + this->data.entries[node].cost;
        
        // After removal: prev -> next (direct, no node cost)
        int afterCost = this->distances[prevNode][nextNode];
        
        // Detour cost = how much worse it gets if we remove this node
        // Higher detour cost = worse node (more expensive to remove)
        int detourCost = currentCost - afterCost;
        
        scores[i] = static_cast<double>(detourCost);
    }
    
    // Find min and max for normalization
    double minScore = *std::min_element(scores.begin(), scores.end());
    double maxScore = *std::max_element(scores.begin(), scores.end());
    double range = maxScore - minScore;
    
    // Normalize and compute probabilities (higher score = higher probability)
    // Add small epsilon to avoid zero probabilities
    std::vector<double> probabilities(n);
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        // Normalize to [0, 1] and add small base probability for randomization
        double normalized = range > 0 ? (scores[i] - minScore) / range : 0.5;
        probabilities[i] = normalized + 0.1;  // Base 0.1 + normalized [0,1] = [0.1, 1.1]
        sum += probabilities[i];
    }
    
    // Normalize to probabilities (sum = 1.0)
    for (int i = 0; i < n; i++) {
        probabilities[i] /= sum;
    }
    
    // Build cumulative distribution for roulette wheel
    std::vector<double> cumulative(n);
    cumulative[0] = probabilities[0];
    for (int i = 1; i < n; i++) {
        cumulative[i] = cumulative[i - 1] + probabilities[i];
    }
    
    // Roulette wheel selection
    std::set<int> removedSet;
    std::vector<bool> selected(n, false);
    
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    
    int removed = 0;
    while (removed < numToRemove) {
        double r = dist(this->mt);
        
        // Find which node to select based on cumulative distribution
        int selectedIdx = 0;
        for (int i = 0; i < n; i++) {
            if (r <= cumulative[i] && !selected[i]) {
                selectedIdx = i;
                break;
            }
        }
        
        // If already selected, find next available
        if (selected[selectedIdx]) {
            // Find first non-selected node
            for (int i = 0; i < n; i++) {
                if (!selected[i]) {
                    selectedIdx = i;
                    break;
                }
            }
        }
        
        selected[selectedIdx] = true;
        removedSet.insert(nodeIndices[selectedIdx]);
        removed++;
    }
    
    // Build partial solution (nodes not removed)
    Indices partialSolution;
    for (int node : solution) {
        if (removedSet.find(node) == removedSet.end()) {
            partialSolution.push_back(node);
        }
    }
    
    return partialSolution;
}

Solver::Indices LNSSolver::repair(Indices const& partialSolution) {
    // Calculate target count: 50% of all nodes (rounded)
    int targetCount = static_cast<int>(std::round(static_cast<float>(this->data.entries.size()) / 2));
    int numToInsert = targetCount - partialSolution.size();
    
    // Start with partial solution
    Indices solution = partialSolution;
    std::vector<bool> visited(this->data.entries.size(), false);
    
    // Mark nodes already in solution as visited
    for (int node : solution) {
        visited[node] = true;
    }
    
    // Build list of ALL available nodes (not in partial solution)
    std::vector<int> availableNodes;
    for (size_t i = 0; i < this->data.entries.size(); i++) {
        if (!visited[i]) {
            availableNodes.push_back(i);
        }
    }
    
    // Use weighted regret to insert nodes (same logic as NearestNeighborPosSolver)
    // repairSolver has alpha=0.5, beta=0.5 for weighted regret
    double alpha = 0.5;
    double beta = 0.5;
    
    // Insert nodes until we reach targetCount
    for (int insertCount = 0; insertCount < numToInsert && !availableNodes.empty(); insertCount++) {
        int bestCity = -1;
        int bestPos = -1;
        double maxScore = std::numeric_limits<double>::lowest();
        int minDeltaForMaxScore = std::numeric_limits<int>::max();
        
        // For each available node, find best position and compute regret
        for (int city : availableNodes) {
            if (visited[city]) continue;
            
            // Find best and second-best positions for this city
            int minDelta1 = std::numeric_limits<int>::max();
            int minDelta2 = std::numeric_limits<int>::max();
            int bestPosForCity = -1;
            
            // Try all positions
            for (size_t pos = 0; pos <= solution.size(); pos++) {
                int delta = getDelta(city, pos, solution);
                
                if (delta < minDelta1) {
                    minDelta2 = minDelta1;
                    minDelta1 = delta;
                    bestPosForCity = pos;
                } else if (delta < minDelta2) {
                    minDelta2 = delta;
                }
            }
            
            int regret = minDelta2 - minDelta1;
            // Weighted regret: α * regret + β * (-best_cost)
            double score = alpha * regret + beta * (-minDelta1);
            
            if (score > maxScore || (score == maxScore && minDelta1 < minDeltaForMaxScore)) {
                maxScore = score;
                bestCity = city;
                bestPos = bestPosForCity;
                minDeltaForMaxScore = minDelta1;
            }
        }
        
        // Insert best city at best position
        if (bestCity != -1) {
            solution.insert(solution.begin() + bestPos, bestCity);
            visited[bestCity] = true;
            
            // Remove from availableNodes
            availableNodes.erase(std::remove(availableNodes.begin(), availableNodes.end(), bestCity), availableNodes.end());
        } else {
            // Fallback: insert remaining greedily
            break;
        }
    }
    
    // Fallback: greedy insertion for any remaining nodes
    for (int city : availableNodes) {
        if (visited[city] || solution.size() >= targetCount) continue;
        
        int bestPos = -1;
        int minDelta = std::numeric_limits<int>::max();
        
        for (size_t pos = 0; pos <= solution.size(); pos++) {
            int delta = getDelta(city, pos, solution);
            if (delta < minDelta) {
                minDelta = delta;
                bestPos = pos;
            }
        }
        
        if (bestPos != -1) {
            solution.insert(solution.begin() + bestPos, city);
            visited[city] = true;
        }
    }
    
    return solution;
}

int LNSSolver::getDelta(int city, int pos, Indices const& solution) const {
    if (solution.empty()) {
        return this->data.entries[city].cost;
    }
    
    if (pos == 0) {
        auto after = solution[pos];
        return this->data.entries[city].cost + this->distances[city][after];
    } else if (pos == solution.size()) {
        auto before = solution[pos - 1];
        return this->data.entries[city].cost + this->distances[before][city];
    }
    auto before = solution[pos - 1];
    auto after = solution[pos];
    return this->data.entries[city].cost
        + this->distances[before][city]
        + this->distances[city][after]
        - this->distances[before][after];
}

