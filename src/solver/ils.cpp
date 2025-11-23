#include <evolutionary_computation/solver/ils.h>
#include <algorithm>
#include <limits>
#include <chrono>

ILSSolver::ILSSolver(LocalSearchSolver& baseLS, 
                     Solver::Duration targetTime,
                     int perturbationSize,
                     PerturbationType pertType)
    : baseLS{baseLS}, 
      targetTime{targetTime}, 
      perturbationSize{perturbationSize},
      pertType{pertType},
      lsRuns{0} {}

std::string ILSSolver::name() const {
    std::string pertName = (pertType == PerturbationType::TwoOpt) 
        ? "2opt" + std::to_string(perturbationSize)
        : "swap" + std::to_string(perturbationSize);
    return std::string("ils-") + baseLS.name() + "-pert" + pertName;
}

void ILSSolver::init(Data const& data) {
    Solver::init(data);
    baseLS.init(data);
    lsRuns = 0;
}

Solver::Indices ILSSolver::_solve(int i) {
    lsRuns = 0;
    
    // Start with random solution and run local search
    auto x = baseLS._solve(i);  // Random init + LS
    int xCost = this->cost(x);
    lsRuns = 1;
    
    auto startTime = std::chrono::high_resolution_clock::now();
    
    while (true) {
        auto currentTime = std::chrono::high_resolution_clock::now();
        auto elapsed = currentTime - startTime;
        
        if (elapsed >= targetTime) {
            break;  // Stop when we've used target time
        }
        
        // Perturb current best solution
        auto y = perturb(x);
        
        // Run local search on perturbed solution
        y = baseLS._solveFromSolution(y);
        lsRuns++;
        
        int yCost = this->cost(y);
        
        // Accept if better (greedy acceptance)
        if (yCost < xCost) {
            x = y;
            xCost = yCost;
        }
    }
    
    return x;
}

Solver::Indices ILSSolver::perturb(Indices const& solution) {
    switch (pertType) {
        case PerturbationType::TwoOpt:
            return perturbTwoOpt(solution);
        case PerturbationType::NodeSwap:
            return perturbNodeSwap(solution);
    }
    return solution;
}

Solver::Indices ILSSolver::perturbTwoOpt(Indices const& solution) {
    Indices perturbed = solution;
    int n = solution.size();
    
    // Apply k random 2-opt moves (without checking improvement)
    for (int k = 0; k < perturbationSize; k++) {
        // Pick two random indices
        int idx1 = this->mt() % (n - 1);
        int idx2 = (idx1 + 1 + (this->mt() % (n - idx1 - 1))) % n;
        
        // Ensure idx1 < idx2
        if (idx1 > idx2) {
            std::swap(idx1, idx2);
        }
        
        // Apply 2-opt: reverse segment from idx1+1 to idx2
        if (idx1 + 1 < idx2) {
            std::reverse(perturbed.begin() + idx1 + 1, 
                        perturbed.begin() + idx2 + 1);
        }
    }
    
    return perturbed;
}

Solver::Indices ILSSolver::perturbNodeSwap(Indices const& solution) {
    Indices perturbed = solution;
    int n = solution.size();
    
    // Apply k random node swaps (without checking improvement)
    for (int k = 0; k < perturbationSize; k++) {
        // Pick two random indices
        int idx1 = this->mt() % n;
        int idx2 = this->mt() % n;
        
        // Swap if different
        if (idx1 != idx2) {
            std::swap(perturbed[idx1], perturbed[idx2]);
        }
    }
    
    return perturbed;
}

