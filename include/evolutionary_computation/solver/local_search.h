#pragma once

#include <evolutionary_computation/solver/base.h>

enum class LocalSearchType {
    Greedy,
    Steep,
};

enum class IntraNeighborhoodType {
    Nodes,
    Edges,
};

enum class MoveType {
    Selection,
    Nodes,
    Edges,
};

class LocalSearchSolver : public Solver {
public:
    using Solver::Solver;
    typedef std::tuple<int, int, MoveType> Move;

    LocalSearchSolver(LocalSearchType localSearchType,
                      IntraNeighborhoodType InterNeighborhoodType,
                      Solver& initSolver);


    virtual std::string name() const override;
    virtual void init(Data const& data) override;
protected:
    const LocalSearchType localSearchType;
    const IntraNeighborhoodType intraNeighborhoodType;
    Solver& initSolver;

    virtual Indices _solve(int i) override;
    void doGreedy(Indices& solution);
    void doSteep(Indices& solution);

    std::vector<Move> getMoves(Indices const& solution);
    int evalMove(Indices const& solution, Move& move);
    void applyMove(Indices& solution, Move& move);

    int getDelta(int target, int city1, int city2);
};
