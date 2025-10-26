#pragma once

#include <evolutionary_computation/solver/base.h>

class NearestNeighborSolver : public Solver {
public:
    using Solver::Solver;

    virtual std::string name() const override;
    virtual Indices _solve(int i) override;
};
