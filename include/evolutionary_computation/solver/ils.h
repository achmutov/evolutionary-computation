#pragma once

#include <evolutionary_computation/solver/base.h>
#include <evolutionary_computation/solver/local_search.h>
#include <chrono>
#include <limits>

enum class PerturbationType {
    TwoOpt,      // Random k-opt moves (2-opt)
    NodeSwap     // Random node swaps
};

class ILSSolver : public Solver {
public:
    ILSSolver(LocalSearchSolver& baseLS, 
              Solver::Duration targetTime,
              int perturbationSize = 4,
              PerturbationType pertType = PerturbationType::TwoOpt);
    
    virtual std::string name() const override;
    virtual void init(Data const& data) override;
    virtual Indices _solve(int i) override;
    
    // Get number of local search runs (for reporting)
    int getLSRuns() const override { return lsRuns; }
    
    // Set target time (for dynamic adjustment per instance)
    void setTargetTime(Solver::Duration newTargetTime) { targetTime = newTargetTime; }
    
protected:
    LocalSearchSolver& baseLS;
    Solver::Duration targetTime;  // Non-const to allow per-instance adjustment
    const int perturbationSize;
    const PerturbationType pertType;
    mutable int lsRuns;  // Track number of LS runs
    
    Indices perturb(Indices const& solution);
    Indices perturbTwoOpt(Indices const& solution);
    Indices perturbNodeSwap(Indices const& solution);
};

