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
                      Solver& initSolver,
                      int nCandidateMoves = -1);


    virtual std::string name() const override;
    virtual void init(Data const& data) override;
    virtual Indices _solve(int i) override;
protected:
    const LocalSearchType localSearchType;
    const IntraNeighborhoodType intraNeighborhoodType;
    const int nCandidateMoves;
    std::vector<Indices> closestCities;
    std::vector<std::vector<bool>> isClosestCity;
    Solver& initSolver;

    void learnCandidateMoves();

    void doGreedy(Indices& solution);
    void doSteep(Indices& solution);

    std::vector<Move> getMoves(Indices const& solution);
    std::vector<Move> buildMoves(Indices const& solution);
    int evalMove(Indices const& solution, Move& move);
    void applyMove(Indices& solution, Move& move);

    int getDelta(int target, int city1, int city2);
};
