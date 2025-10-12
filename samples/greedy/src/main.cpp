#include <fstream>
#include <iostream>
#include <evolutionary_computation/loader/csv.h>
#include <evolutionary_computation/solver/random.h>
#include <evolutionary_computation/solver/nearest_neighbor.h>
#include <evolutionary_computation/solver/nearest_neighbor_pos.h>
#include <evolutionary_computation/solver/greedy_cycle.h>
#include <limits>
#include <memory>
#include <stdexcept>
#include <ranges>
#include <limits.h>
#include <vector>

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "usage: " << argv[0] << " N_SAMPLES  SOURCES..." << std::endl;
        std::exit(64);
    }

    int n{};
    try {
        n = std::stoi(argv[1]);
    } catch (std::invalid_argument) {
        std::cerr << "Error: N_SAMPLES must be a valid integer" << std::endl;
        std::exit(64);
    }

    std::cout << "instance,method,min,max,mean,best_solution\n";

    auto solvers = std::vector<Solver*>();

    auto randomSolver= std::make_unique<RandomSolver>();
    solvers.push_back(randomSolver.get());

    auto nearestNeighborSolver = std::make_unique<NearestNeighborSolver>();
    solvers.push_back(nearestNeighborSolver.get());

    auto nearestNeighborPosSolver = std::make_unique<NearestNeighborPosSolver>();
    solvers.push_back(nearestNeighborPosSolver.get());

    auto greedyCycle = std::make_unique<GreedyCycle>();
    solvers.push_back(greedyCycle.get());

    for (int i = 2; i < argc; i++) {
        auto const instancePath = argv[i];
        auto stream = std::ifstream(instancePath);
        auto loader = CSVLoader(stream);
        auto data = loader.load();

        for (auto solver : solvers) {
            solver->init(data);
            auto minObj = std::numeric_limits<int>::max();
            auto maxObj = std::numeric_limits<int>::min();
            float meanObj = 0;
            Solver::Indices best;

            for (int j = 0; j < n; j++) {
                auto const [indices, candidate] = solver->solve(j);
                maxObj = std::max(maxObj, candidate);
                if (candidate < minObj) {
                    minObj = candidate;
                    best = indices;
                }
                minObj = std::min(minObj, candidate);
                meanObj += candidate;
            }
            meanObj /= n;

            // Rest of the stats
            std::cout << instancePath
                << ',' << solver->name()
                << ',' << minObj
                << ',' << maxObj
                << ',' << meanObj
                << ',';

            // Best indices
            auto const allButLast = best | std::ranges::views::take(best.size() - 1);
            for (auto const& index : allButLast)
                std::cout << index << ' ';
            std::cout << best.back() << "\n";
        }
    }
    return 0;
}
