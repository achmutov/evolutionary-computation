#pragma once

#include <evolutionary_computation/solver/base.h>

class RandomSolver : public Solver {
public:
    using Solver::Solver;

protected:
    virtual Indices _solve() override;
};
