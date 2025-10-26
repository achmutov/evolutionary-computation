#include <evolutionary_computation/solver/random.h>
#include <evolutionary_computation/solver/local_search.h>
#include <evolutionary_computation/cli.h>
#include <memory>
#include <vector>

int main(int argc, char* argv[]) {
    auto uniqueSolvers = std::vector<std::unique_ptr<Solver>>();
    auto initSolver = RandomSolver();
    for (int localSearchType_ = 0; localSearchType_ < 2; localSearchType_++) {
        auto localSearchType = static_cast<LocalSearchType>(localSearchType_);
        for (int interNeighborhoodType_ = 0; interNeighborhoodType_ < 2; interNeighborhoodType_++) {
            auto interNeighborhoodType = static_cast<IntraNeighborhoodType>(interNeighborhoodType_);
            uniqueSolvers.push_back(std::make_unique<LocalSearchSolver>(localSearchType, interNeighborhoodType, initSolver));
        }
    }

    auto solvers = std::vector<Solver*>();
    for (auto& solver : uniqueSolvers) {
        solvers.push_back(solver.get());
    }
    cli(argc, argv, solvers);

    return 0;
}
