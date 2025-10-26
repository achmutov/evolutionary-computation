#include <evolutionary_computation/solver/random.h>
#include <evolutionary_computation/solver/local_search.h>
#include <evolutionary_computation/solver/nearest_neighbor_pos.h>
#include <evolutionary_computation/cli.h>
#include <memory>
#include <vector>
#include <iostream>

class InstanceSolver : public Solver {
    std::unique_ptr<Solver> solver;
public:
    void init(Data const& data) override {
        if (data.filename.has_value()) {
            auto value = data.filename.value();
            auto lowered = std::string(value);
            std::transform(value.begin(), value.end(), lowered.begin(), [](char c) { return std::tolower(c); });
            std::cout << value << ' ' << lowered << '\n';
            if (lowered.contains("tspa")) {
                // this->solver = std::make_unique<NearestNeighborPosSolver>(Mode::WeightedRegret);
                this->solver = std::make_unique<NearestNeighborPosSolver>();
                this->solver->init(data);
            } else if (lowered.contains("tspb")) {
                this->solver = std::make_unique<NearestNeighborPosSolver>();
                this->solver->init(data);
            } else {
#ifdef EXCEPTIONS
                // throw std::runtime_error("Instance unknown: " + value);
#endif
            }
        } else {
#ifdef EXCEPTIONS
            throw std::runtime_error("Filename required");
#endif
        }
    }

    std::string name() const override {
        return this->solver.get()->name();
    }

    Indices _solve(int i) override {
        return this->solver.get()->_solve(i);
    }
};

int main(int argc, char* argv[]) {
    auto uniqueSolvers = std::vector<std::unique_ptr<Solver>>();
    auto randomSolver = RandomSolver();
    auto instanceSolver = InstanceSolver();
    for (int localSearchType_ = 0; localSearchType_ < 2; localSearchType_++) {
        auto localSearchType = static_cast<LocalSearchType>(localSearchType_);
        for (int interNeighborhoodType_ = 0; interNeighborhoodType_ < 2; interNeighborhoodType_++) {
            auto interNeighborhoodType = static_cast<IntraNeighborhoodType>(interNeighborhoodType_);
            uniqueSolvers.push_back(std::make_unique<LocalSearchSolver>(localSearchType, interNeighborhoodType, randomSolver));
            uniqueSolvers.push_back(std::make_unique<LocalSearchSolver>(localSearchType, interNeighborhoodType, instanceSolver));
        }
    }

    auto solvers = std::vector<Solver*>();
    for (auto& solver : uniqueSolvers) {
        solvers.push_back(solver.get());
    }
    cli(argc, argv, solvers);

    return 0;
}
