#include <evolutionary_computation/solver/msls.h>
#include <limits>

MSLSSolver::MSLSSolver(LocalSearchSolver& baseLS, int numIterations)
    : baseLS{baseLS}, numIterations{numIterations} {}

std::string MSLSSolver::name() const {
    return std::string("msls-") + baseLS.name() + "-" + std::to_string(numIterations);
}

void MSLSSolver::init(Data const& data) {
    Solver::init(data);
    baseLS.init(data);
}

Solver::Indices MSLSSolver::_solve(int i) {
    Indices bestSolution;
    int bestCost = std::numeric_limits<int>::max();
    
    // Run local search numIterations times, each from a random start
    for (int iter = 0; iter < numIterations; iter++) {
        auto solution = baseLS._solve(iter);  // Each iteration gets different random seed
        int cost = this->cost(solution);
        
        if (cost < bestCost) {
            bestCost = cost;
            bestSolution = solution;
        }
    }
    
    return bestSolution;
}

