#pragma once

#include <evolutionary_computation/solver/base.h>
#include <evolutionary_computation/solver/local_search.h>
#include <limits>

class MSLSSolver : public Solver {
public:
    MSLSSolver(LocalSearchSolver& baseLS, int numIterations = 200);
    
    virtual std::string name() const override;
    virtual void init(Data const& data) override;
    virtual Indices _solve(int i) override;
    
    // Get number of local search runs (for reporting) - always numIterations
    int getLSRuns() const override { return numIterations; }
    
protected:
    LocalSearchSolver& baseLS;
    const int numIterations;
};

