#include <chrono>
#include <evolutionary_computation/cli.h>
#include <evolutionary_computation/solver/ils.h>
#include <iostream>
#include <fstream>
#include <evolutionary_computation/loader/csv.h>
#include <ranges>
#include <ratio>
#include <algorithm>
#include <string>

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

    std::cout << "instance,method,min,max,mean,duration_min,duration_max,duration_mean,ls_runs_min,ls_runs_max,ls_runs_mean,best_solution\n";

    for (int i = 2; i < argc; i++) {
        auto const instancePath = argv[i];
        auto stream = std::ifstream(instancePath);
        auto loader = CSVLoader(stream, instancePath);
        auto data = loader.load();

        // Find MSLS solver and ILS solvers to set dynamic target time
        // We identify them by name prefix to avoid RTTI
        Solver* mslsSolver = nullptr;
        std::vector<Solver*> ilsSolvers;
        
        for (auto solver : solvers) {
            std::string name = solver->name();
            if (name.find("msls-") == 0) {
                mslsSolver = solver;
            } else if (name.find("ils-") == 0) {
                ilsSolvers.push_back(solver);
            }
        }
        
        // If we have MSLS and ILS, measure MSLS time and update ILS target
        if (mslsSolver && !ilsSolvers.empty()) {
            mslsSolver->init(data);
            
            // Measure MSLS average time on a few runs
            Solver::Duration totalTime = std::chrono::nanoseconds(0);
            int measureRuns = std::min(n, 5);  // Measure on first 5 runs (or n if n < 5)
            
            for (int j = 0; j < measureRuns; j++) {
                auto [indices, cost, duration] = mslsSolver->solve(j);
                totalTime += duration;
            }
            
            Solver::Duration avgMSLSTime = totalTime / measureRuns;
            
            // Update ILS target time (cast to ILSSolver* - safe because we checked name)
            for (auto* solver : ilsSolvers) {
                // We know it's an ILSSolver based on name check
                static_cast<ILSSolver*>(solver)->setTargetTime(avgMSLSTime);
            }
        }

        for (auto solver : solvers) {
            solver->init(data);
            auto minObj = std::numeric_limits<int>::max();
            auto maxObj = std::numeric_limits<int>::min();
            float meanObj = 0;
            Solver::Indices best;
            Solver::Duration minDuration = std::chrono::nanoseconds::max();
            Solver::Duration meanDuration = std::chrono::nanoseconds(0);
            Solver::Duration maxDuration = std::chrono::nanoseconds::min();
            
            // Track LS runs
            int minLSRuns = std::numeric_limits<int>::max();
            int maxLSRuns = std::numeric_limits<int>::min();
            float meanLSRuns = 0;
            bool hasLSRuns = false;

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
                
                // Get LS runs if available (after solve) - uses virtual method
                int lsRuns = solver->getLSRuns();
                if (lsRuns >= 0) {
                    minLSRuns = std::min(minLSRuns, lsRuns);
                    maxLSRuns = std::max(maxLSRuns, lsRuns);
                    meanLSRuns += lsRuns;
                    hasLSRuns = true;
                }

            }
            meanObj /= n;
            meanDuration /= n;
            if (hasLSRuns) {
                meanLSRuns /= n;
            }

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
            
            // LS runs stats
            if (hasLSRuns) {
                std::cout << minLSRuns << ',' << maxLSRuns << ',' << meanLSRuns << ',';
            } else {
                std::cout << ",,";  // Empty if not available
            }
            
            // Best indices
            auto const allButLast = best | std::ranges::views::take(best.size() - 1);
            for (auto const& index : allButLast)
                std::cout << index << ' ';
            std::cout << best.back() << "\n";
        }
    }
}


