#include <evolutionary_computation/solver/random.h>
#include <evolutionary_computation/solver/local_search.h>
#include <evolutionary_computation/solver/lns.h>
#include <evolutionary_computation/cli.h>
#include <memory>
#include <vector>
#include <chrono>

int main(int argc, char* argv[]) {
    // Create base local search solver: Steep, 2-edges, Random init, with candidate moves
    auto randomSolver = RandomSolver();  // Random initializer
    auto baseLS = LocalSearchSolver(
        LocalSearchType::Steep,
        IntraNeighborhoodType::Edges,
        randomSolver,  // init = Random
        -1,            // candidate moves = 10
        false          // useMoveList = false (no LM)
    );
    
    // Target time: 1050ms (measured from MSLS)
    auto targetTime = std::chrono::milliseconds(1050);
    
    auto uniqueSolvers = std::vector<std::unique_ptr<Solver>>();
    
    // LNS variants – Random destroy only
    uniqueSolvers.push_back(std::make_unique<LNSSolver>(
        baseLS, targetTime, DestroyType::Random, true, 0.3));  // with LS
    uniqueSolvers.push_back(std::make_unique<LNSSolver>(
        baseLS, targetTime, DestroyType::Random, false, 0.3)); // without LS
    
    // Single subpath destroy (disabled)
    // uniqueSolvers.push_back(std::make_unique<LNSSolver>(
    //     baseLS, targetTime, DestroyType::SingleSubpath, true, 0.5));  // with LS
    // uniqueSolvers.push_back(std::make_unique<LNSSolver>(
    //     baseLS, targetTime, DestroyType::SingleSubpath, false, 0.5)); // without LS
    
    // Multiple subpaths destroy (ENABLED, 50% destroy)
    // uniqueSolvers.push_back(std::make_unique<LNSSolver>(
    //     baseLS, targetTime, DestroyType::MultipleSubpaths, true, 0.3));   // with LS
    // uniqueSolvers.push_back(std::make_unique<LNSSolver>(
    //     baseLS, targetTime, DestroyType::MultipleSubpaths, false, 0.3));  // without LS
    
    // Heuristic destroy (disabled)
    // uniqueSolvers.push_back(std::make_unique<LNSSolver>(
    //     baseLS, targetTime, DestroyType::Heuristic, true, 0.5));  // with LS
    // uniqueSolvers.push_back(std::make_unique<LNSSolver>(
    //     baseLS, targetTime, DestroyType::Heuristic, false, 0.5)); // without LS
    
    auto solvers = std::vector<Solver*>();
    for (auto& solver : uniqueSolvers) {
        solvers.push_back(solver.get());
    }
    
    cli(argc, argv, solvers);
    
    return 0;
}

