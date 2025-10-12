#include <fstream>
#include <iostream>
#include <evolutionary_computation/loader/csv.h>
#include <evolutionary_computation/solver/random.h>
#include <stdexcept>
#include <ranges>
#include <limits.h>

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

    for (int i = 2; i < argc; i++) {
        auto const instancePath = argv[i];
        auto stream = std::ifstream(instancePath);
        auto loader = CSVLoader(stream);
        auto data = loader.load();
        auto solver = RandomSolver(data);
        float minObj = INT_MAX;
        float maxObj = -1;
        float meanObj = 0;
        Solver::Indices best;
        for (int j = 0; j < n; j++) {
            auto const [indices, candidate] = solver.solve();
            auto const fcandidate = static_cast<float>(candidate);
            maxObj = std::max(maxObj, fcandidate);
            if (minObj < fcandidate) {
                minObj = fcandidate;
                best = indices;
            }
            minObj = std::min(minObj, static_cast<float>(candidate));
            meanObj += candidate;
        }
        meanObj /= n;

        // Rest of the stats
        std::cout << instancePath << ",random," << minObj << ',' << maxObj << ',' << meanObj << ',';

        // Best indices
        auto const allButLast = best | std::ranges::views::take(best.size() - 1);
        for (auto const& index : allButLast)
            std::cout << index << ';';
        std::cout << best.back() << "\n";
    }
    return 0;
}
