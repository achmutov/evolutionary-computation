#pragma once

#include <evolutionary_computation/solver/base.h>
#include <evolutionary_computation/solver/local_search.h>
#include <evolutionary_computation/solver/nearest_neighbor_pos.h>
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
              NearestNeighborPosSolver& repairSolver,
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
    NearestNeighborPosSolver& repairSolver;
    Solver::Duration targetTime;  // Non-const to allow per-instance adjustment
    const DestroyType destroyType;
    const bool useLocalSearchAfterRepair;
    const double destroyFraction;
    mutable int iterations;  // Track number of destroy-repair iterations
    
    // Destroy operators
    std::pair<Indices, std::set<int>> destroy(Indices const& solution);
    std::pair<Indices, std::set<int>> destroyRandom(Indices const& solution);
    std::pair<Indices, std::set<int>> destroySingleSubpath(Indices const& solution);
    std::pair<Indices, std::set<int>> destroyMultipleSubpaths(Indices const& solution);
    std::pair<Indices, std::set<int>> destroyHeuristic(Indices const& solution);
    
    // Repair operator: insert removed nodes back using greedy heuristic
    Indices repair(Indices const& partialSolution, std::set<int> const& removedNodes);
    
    // Helper method for repair: compute delta of inserting city at position
    int getDelta(int city, int pos, Indices const& solution) const;
};

