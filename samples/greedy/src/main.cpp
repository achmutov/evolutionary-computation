#include <fstream>
#include <iostream>
#include <evolutionary_computation/loader/csv.h>
#include <evolutionary_computation/solver/random.h>
#include <stdexcept>

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

    for (int i = 2; i < argc; i++) {
        auto stream = std::ifstream(argv[i]);
        auto loader = CSVLoader(stream);
        auto data = loader.load();
        auto solver = RandomSolver(data);
        for (int j = 0; j < n; j++) {
            auto [indices, result] = solver.solve();
            for (auto i : indices) {
                std::cout << i << ' ';
            }
            std::cout << '\n' << result << '\n';
        }
    }
    return 0;
}
