#include <evolutionary_computation/cli.h>
#include <iostream>
#include <fstream>
#include <evolutionary_computation/loader/csv.h>
#include <ranges>

void cli(int argc, char* argv[], std::vector<Solver*>& solvers) {
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
}


