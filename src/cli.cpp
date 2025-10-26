#include <chrono>
#include <evolutionary_computation/cli.h>
#include <iostream>
#include <fstream>
#include <evolutionary_computation/loader/csv.h>
#include <ranges>
#include <ratio>

void cli(int argc, char* argv[], std::vector<Solver*>& solvers) {
    if (argc < 3) {
        std::cerr << "usage: " << argv[0] << " N_SAMPLES  SOURCES..." << std::endl;
        std::exit(64);
    }

    int n{};
#ifdef EXCEPTIONS
    try {
#endif
        n = std::stoi(argv[1]);
#ifdef EXCEPTIONS
    } catch (std::invalid_argument) {
        std::cerr << "Error: N_SAMPLES must be a valid integer" << std::endl;
        std::exit(64);
    }
#endif

    std::cout << "instance,method,min,max,mean,duration_min,duration_max,duration_mean,best_solution\n";

    for (int i = 2; i < argc; i++) {
        auto const instancePath = argv[i];
        auto stream = std::ifstream(instancePath);
        auto loader = CSVLoader(stream, instancePath);
        auto data = loader.load();

        for (auto solver : solvers) {
            solver->init(data);
            auto minObj = std::numeric_limits<int>::max();
            auto maxObj = std::numeric_limits<int>::min();
            float meanObj = 0;
            Solver::Indices best;
            Solver::Duration minDuration = std::chrono::nanoseconds::max();
            Solver::Duration meanDuration = std::chrono::nanoseconds(0);
            Solver::Duration maxDuration = std::chrono::nanoseconds::min();

            for (int j = 0; j < n; j++) {
                auto const [indices, candidate, duration] = solver->solve(j);
                maxObj = std::max(maxObj, candidate);
                if (candidate < minObj) {
                    minObj = candidate;
                    best = indices;
                }
                minObj = std::min(minObj, candidate);
                meanObj += candidate;

                minDuration = std::min(minDuration, duration);
                maxDuration = std::max(maxDuration, duration);
                meanDuration += duration;

            }
            meanObj /= n;
            meanDuration /= n;

            // Rest of the stats
            std::cout << instancePath
                << ',' << solver->name()
                << ',' << minObj
                << ',' << maxObj
                << ',' << meanObj
                << ',' << std::chrono::duration<double, std::milli>(minDuration).count()
                << ',' << std::chrono::duration<double, std::milli>(maxDuration).count()
                << ',' << std::chrono::duration<double, std::milli>(meanDuration).count()
                << ',';

            // Best indices
            auto const allButLast = best | std::ranges::views::take(best.size() - 1);
            for (auto const& index : allButLast)
                std::cout << index << ' ';
            std::cout << best.back() << "\n";
        }
    }
}


