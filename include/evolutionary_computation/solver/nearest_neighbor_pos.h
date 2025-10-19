#pragma once

#include <evolutionary_computation/solver/base.h>
#include <evolutionary_computation/solver/mode.h>

class NearestNeighborPosSolver : public Solver {
public:
    using Solver::Solver;

    NearestNeighborPosSolver(Mode mode = Mode::None, double alpha = 1.0, double beta = 1.0);

    virtual std::string name() const override;
protected:
    Mode mode;
    double alpha, beta;

    virtual Indices _solve(int i) override;

    int get_delta(int city, int pos, Indices& indices);
    std::tuple<int, int> next_greedy(Indices& indices, std::vector<bool>& visited);

    std::tuple<int, int, int> get_regret(int city, Indices& indices);
    std::tuple<int, int> next_regret(Indices& indices, std::vector<bool>& visited);
    std::tuple<int, int> next_weighted_regret(Indices& indices, std::vector<bool>& visited);
};
