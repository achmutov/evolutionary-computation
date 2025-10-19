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
    auto nearestNeighborPosSolver1 = std::make_unique<NearestNeighborPosSolver>(Mode::Regret);
    auto nearestNeighborPosSolver2 = std::make_unique<NearestNeighborPosSolver>(Mode::WeightedRegret);
    auto greedyCycle1 = std::make_unique<GreedyCycle>(Mode::Regret);
    auto greedyCycle2 = std::make_unique<GreedyCycle>(Mode::WeightedRegret);
    auto solvers = std::vector<Solver*> {
        nearestNeighborPosSolver1.get(),
        nearestNeighborPosSolver2.get(),
        greedyCycle1.get(),
        greedyCycle2.get(),
    };

    cli(argc, argv, solvers);
    return 0;
}
