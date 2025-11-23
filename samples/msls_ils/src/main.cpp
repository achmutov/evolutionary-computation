#include <evolutionary_computation/solver/random.h>
#include <evolutionary_computation/solver/local_search.h>
#include <evolutionary_computation/solver/msls.h>
#include <evolutionary_computation/solver/ils.h>
#include <evolutionary_computation/cli.h>
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
    
    auto uniqueSolvers = std::vector<std::unique_ptr<Solver>>();
    
    // MSLS solver - 200 iterations
    uniqueSolvers.push_back(std::make_unique<MSLSSolver>(baseLS, 200));
    
    // ILS solvers with different perturbation types
    // Target time will be set dynamically per instance by CLI based on MSLS measurement
    // Using a placeholder initial value (will be updated by CLI)
    auto targetTime = std::chrono::milliseconds(1050);
    
    // ILS with 2-opt perturbation (pert4)
    uniqueSolvers.push_back(std::make_unique<ILSSolver>(
        baseLS, targetTime, 4, PerturbationType::TwoOpt));  // pert4 (2-opt)
    
    // ILS with node swap perturbation (pert4)
    uniqueSolvers.push_back(std::make_unique<ILSSolver>(
        baseLS, targetTime, 4, PerturbationType::NodeSwap));  // pert4 (node swap)
    
    auto solvers = std::vector<Solver*>();
    for (auto& solver : uniqueSolvers) {
        solvers.push_back(solver.get());
    }
    
    cli(argc, argv, solvers);
    
    return 0;
}

