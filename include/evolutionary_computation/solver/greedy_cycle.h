#pragma once

#include <evolutionary_computation/solver/base.h>

class GreedyCycle : public Solver {
public:
    using Solver::Solver;

    virtual std::string name() const override;
protected:
    virtual Indices _solve(int i) override;
};
