#pragma once

#include <evolutionary_computation/solver/base.h>

class RandomSolver : public Solver {
public:
    using Solver::Solver;

    virtual std::string name() const override;
    virtual Indices _solve(int i) override;
};
