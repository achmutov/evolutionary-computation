#include <evolutionary_computation/loader/csv.h>
#include <evolutionary_computation/solver/random.h>
#include <evolutionary_computation/solver/nearest_neighbor.h>
#include <evolutionary_computation/solver/nearest_neighbor_pos.h>
#include <evolutionary_computation/solver/greedy_cycle.h>
#include <evolutionary_computation/cli.h>
#include <memory>
#include <limits.h>
#include <vector>

int main(int argc, char* argv[]) {
    auto randomSolver= std::make_unique<RandomSolver>();
    auto nearestNeighborSolver = std::make_unique<NearestNeighborSolver>();
    auto nearestNeighborPosSolver = std::make_unique<NearestNeighborPosSolver>();
    auto greedyCycle = std::make_unique<GreedyCycle>();
    auto solvers = std::vector<Solver*> {
        randomSolver.get(),
        nearestNeighborSolver.get(),
        nearestNeighborPosSolver.get(),
        greedyCycle.get(),
    };

    cli(argc, argv, solvers);
    return 0;
}
