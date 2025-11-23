#pragma once

#include <evolutionary_computation/solver/base.h>
#include <unordered_map>
#include <unordered_set>
#include <set>

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

// Move type for LM (List of Moves) - stores nodes and delta
enum class LMMoveType {
    TWO_OPT,                    // Intra-edge (2-opt)
    EXCHANGE_SELECTED_UNSELECTED // Inter-node
};

struct MoveWithDelta {
    LMMoveType type;
    int delta;
    int nodeA, nodeB, nodeC, nodeD; // Nodes defining the move
    
    static MoveWithDelta forIntraEdge(int delta, int nodeA, int nodeB, int nodeC, int nodeD) {
        MoveWithDelta move;
        move.type = LMMoveType::TWO_OPT;
        move.delta = delta;
        move.nodeA = nodeA;
        move.nodeB = nodeB;
        move.nodeC = nodeC;
        move.nodeD = nodeD;
        return move;
    }
    
    static MoveWithDelta forInterRoute(int delta, int prev, int next, int nodeInCycle, int nodeOutOfCycle) {
        MoveWithDelta move;
        move.type = LMMoveType::EXCHANGE_SELECTED_UNSELECTED;
        move.delta = delta;
        move.nodeA = prev;
        move.nodeB = next;
        move.nodeC = nodeInCycle;
        move.nodeD = nodeOutOfCycle;
        return move;
    }
};

struct ValidationResult {
    bool applyMove;
    bool keepMove;
    bool applyReversed;
    
    ValidationResult(bool apply, bool keep, bool reversed)
        : applyMove(apply), keepMove(keep), applyReversed(reversed) {}
};

class LocalSearchSolver : public Solver {
public:
    using Solver::Solver;
    typedef std::tuple<int, int, MoveType> Move;

    LocalSearchSolver(LocalSearchType localSearchType,
                      IntraNeighborhoodType InterNeighborhoodType,
                      Solver& initSolver,
                      int nCandidateMoves = -1,
                      bool useMoveList = false);


    virtual std::string name() const override;
    virtual void init(Data const& data) override;
    virtual Indices _solve(int i) override;
    Indices _solveFromSolution(Indices const& initialSolution);  // For ILS
protected:
    const LocalSearchType localSearchType;
    const IntraNeighborhoodType intraNeighborhoodType;
    const int nCandidateMoves;
    const bool useMoveList;
    std::vector<Indices> closestCities;
    std::vector<std::vector<bool>> isClosestCity;
    Solver& initSolver;
    
    // LM (List of Moves) data structures
    std::vector<MoveWithDelta> improvingMoveList;

    void learnCandidateMoves();

    void doGreedy(Indices& solution);
    void doSteep(Indices& solution);
    void doSteepLM(Indices& solution, std::set<int>& remainingNodes);

    std::vector<Move> getMoves(Indices const& solution);
    std::vector<Move> buildMoves(Indices const& solution);
    int evalMove(Indices const& solution, Move& move);
    void applyMove(Indices& solution, Move& move);

    int getDelta(int target, int city1, int city2);
    
    // LM-specific methods
    bool performSteepestStepLM(Indices& solution, std::set<int>& remainingNodes);
    void populateMoveList(Indices const& solution, std::set<int> const& remainingNodes);
    ValidationResult validateMove(MoveWithDelta const& move, 
                                 std::unordered_map<int, int> const& succMap,
                                 std::unordered_map<int, int> const& predMap,
                                 std::set<int> const& remaining) const;
    std::set<int> applyMoveLM(Indices& solution, std::set<int>& remaining, 
                              MoveWithDelta const& move, ValidationResult const& val);
    void updateLocalMoves(Indices const& solution, std::set<int> const& remaining,
                          std::set<int> const& changedNodes, LMMoveType lastMoveType);
    
    // Helper methods for LM
    std::unordered_map<int, int> buildSuccMap(Indices const& solution) const;
    std::unordered_map<int, int> buildPredMap(Indices const& solution) const;
    bool hasEdge(std::unordered_map<int, int> const& succMap, int u, int v) const;
    int deltaInter(Indices const& solution, int selectedIndex, int unselectedNode) const;
    int deltaTwoOpt(Indices const& solution, int i, int j) const;
    int findNodeIndex(Indices const& solution, int node) const;
    void reverseSublist(Indices& solution, int start, int end);
};
