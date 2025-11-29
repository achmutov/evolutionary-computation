#include <evolutionary_computation/solver/random.h>
#include <evolutionary_computation/solver/local_search.h>
#include <evolutionary_computation/solver/nearest_neighbor_pos.h>
#include <evolutionary_computation/solver/msls.h>
#include <evolutionary_computation/solver/lns.h>
#include <evolutionary_computation/cli.h>
#include <evolutionary_computation/solver/mode.h>
#include <memory>
#include <vector>
#include <chrono>

int main(int argc, char* argv[]) {
    // Create base local search solver (steepest, edges, random init, with LM)
    auto randomSolver = RandomSolver();
    auto baseLS = LocalSearchSolver(
        LocalSearchType::Steep,
        IntraNeighborhoodType::Edges,
        randomSolver,
        -1,  // no candidate moves
        true // useMoveList = true (LM)
    );
    
    // Create repair solver: NearestNeighborPos with weighted regret (50/50)
    auto repairSolver = NearestNeighborPosSolver(
        Mode::WeightedRegret,
        0.5,  // alpha
        0.5   // beta
    );
    
    // Create MSLS solver for time measurement (200 iterations)
    auto mslsSolver = std::make_unique<MSLSSolver>(baseLS, 200);
    
    // Placeholder target time (will be updated dynamically by CLI based on MSLS measurement)
    auto targetTime = std::chrono::milliseconds(1050);
    
    auto uniqueSolvers = std::vector<std::unique_ptr<Solver>>();
    
    // Add MSLS solver for time measurement
    uniqueSolvers.push_back(std::move(mslsSolver));
    
    // Create LNS variants with 30% destroy fraction:
    // Only heuristic destroy for testing
    
    // Random destroy - COMMENTED OUT FOR TESTING
    uniqueSolvers.push_back(std::make_unique<LNSSolver>(
        baseLS, repairSolver, targetTime, DestroyType::Random, true, 0.3));  // with LS
    uniqueSolvers.push_back(std::make_unique<LNSSolver>(
        baseLS, repairSolver, targetTime, DestroyType::Random, false, 0.3)); // without LS
    
    // Single subpath destroy - COMMENTED OUT FOR TESTING
    uniqueSolvers.push_back(std::make_unique<LNSSolver>(
        baseLS, repairSolver, targetTime, DestroyType::SingleSubpath, true, 0.3));  // with LS
    uniqueSolvers.push_back(std::make_unique<LNSSolver>(
        baseLS, repairSolver, targetTime, DestroyType::SingleSubpath, false, 0.3)); // without LS
    
    // Multiple subpaths destroy - COMMENTED OUT FOR TESTING
    uniqueSolvers.push_back(std::make_unique<LNSSolver>(
        baseLS, repairSolver, targetTime, DestroyType::MultipleSubpaths, true, 0.3));  // with LS
    uniqueSolvers.push_back(std::make_unique<LNSSolver>(
        baseLS, repairSolver, targetTime, DestroyType::MultipleSubpaths, false, 0.3)); // without LS
    
    // Heuristic destroy (30% destroy, using detour cost)
    uniqueSolvers.push_back(std::make_unique<LNSSolver>(
        baseLS, repairSolver, targetTime, DestroyType::Heuristic, true, 0.3));  // with LS
    uniqueSolvers.push_back(std::make_unique<LNSSolver>(
        baseLS, repairSolver, targetTime, DestroyType::Heuristic, false, 0.3)); // without LS
    
    auto solvers = std::vector<Solver*>();
    for (auto& solver : uniqueSolvers) {
        solvers.push_back(solver.get());
    }
    
    cli(argc, argv, solvers);
    
    return 0;
}

