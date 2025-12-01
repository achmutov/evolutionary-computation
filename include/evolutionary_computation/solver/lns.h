#pragma once

#include <evolutionary_computation/solver/base.h>
#include <evolutionary_computation/solver/local_search.h>
#include <chrono>
#include <limits>
#include <set>

enum class DestroyType {
    Random,           // Remove 30% of nodes randomly
    SingleSubpath,    // Remove single contiguous subpath of 30% nodes
    MultipleSubpaths, // Remove multiple subpaths (5-10 subpaths, proportional lengths)
    Heuristic         // Remove nodes using roulette wheel based on edge length + node cost
};

class LNSSolver : public Solver {
public:
    LNSSolver(LocalSearchSolver& baseLS,
              Solver::Duration targetTime,
              DestroyType destroyType,
              bool useLocalSearchAfterRepair = true,
              double destroyFraction = 0.3);
    
    virtual std::string name() const override;
    virtual void init(Data const& data) override;
    virtual Indices _solve(int i) override;
    
    // Get number of iterations (destroy-repair cycles) for reporting
    int getLSRuns() const override { return iterations; }
    
    // Set target time (for dynamic adjustment per instance)
    void setTargetTime(Solver::Duration newTargetTime) { targetTime = newTargetTime; }
    
protected:
    LocalSearchSolver& baseLS;
    Solver::Duration targetTime;  // Non-const to allow per-instance adjustment
    const DestroyType destroyType;
    const bool useLocalSearchAfterRepair;
    const double destroyFraction;
    mutable int iterations;  // Track number of destroy-repair iterations
    
    // Destroy operators
    Indices destroy(Indices const& solution);
    Indices destroyRandom(Indices const& solution);
    Indices destroySingleSubpath(Indices const& solution);
    Indices destroyMultipleSubpaths(Indices const& solution);
    Indices destroyHeuristic(Indices const& solution);
    
    // Repair operator: insert nodes from all available nodes to reach target count
    Indices repair(Indices const& partialSolution);
    
    // Helper method for repair: compute delta of inserting city at position
    int getDelta(int city, int pos, Indices const& solution) const;
};

