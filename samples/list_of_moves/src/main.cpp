#include <evolutionary_computation/solver/random.h>
#include <evolutionary_computation/solver/local_search.h>
#include <evolutionary_computation/cli.h>
#include <memory>
#include <vector>

int main(int argc, char* argv[]) {
    auto uniqueSolvers = std::vector<std::unique_ptr<Solver>>();
    auto randomSolver = RandomSolver();
    
    // 1. Steepest, edge-swap, WITH move list (LM)
    uniqueSolvers.push_back(std::make_unique<LocalSearchSolver>(
        LocalSearchType::Steep,
        IntraNeighborhoodType::Edges,
        randomSolver,
        -1,  // no candidate moves
        true // useMoveList = true
    ));
    
    // 2. Steepest, edge-swap, WITHOUT move list (baseline)
    uniqueSolvers.push_back(std::make_unique<LocalSearchSolver>(
        LocalSearchType::Steep,
        IntraNeighborhoodType::Edges,
        randomSolver,
        -1,  // no candidate moves
        false // useMoveList = false
    ));

    auto solvers = std::vector<Solver*>();
    for (auto& solver : uniqueSolvers) {
        solvers.push_back(solver.get());
    }
    cli(argc, argv, solvers);

    return 0;
}

