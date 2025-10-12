#pragma once

#include <evolutionary_computation/solver/base.h>

class NearestNeighborPosSolver : public Solver {
public:
    using Solver::Solver;

    virtual std::string name() const override;
protected:
    virtual Indices _solve(int i) override;
};
