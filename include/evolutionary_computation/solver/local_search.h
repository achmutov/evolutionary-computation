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
    void doSteepLM(Indices& solution, std::set<int>& unselectedNodes);

    std::vector<Move> getMoves(Indices const& solution);
    std::vector<Move> buildMoves(Indices const& solution);
    int evalMove(Indices const& solution, Move& move);
    void applyMove(Indices& solution, Move& move);

    int getDelta(int target, int city1, int city2);
    
    // LM-specific methods
    bool executeSteepestStep(Indices& solution, std::set<int>& unselectedNodes);
    void initializeMoveList(Indices const& solution, std::set<int> const& unselectedNodes);
    ValidationResult checkMoveValidity(MoveWithDelta const& move, 
                                 std::unordered_map<int, int> const& nextMap,
                                 std::unordered_map<int, int> const& prevMap,
                                 std::set<int> const& unselected) const;
    std::set<int> executeMove(Indices& solution, std::set<int>& unselected, 
                              MoveWithDelta const& move, ValidationResult const& val);
    void refreshMoveList(Indices const& solution, std::set<int> const& unselected,
                          std::set<int> const& affectedNodes, LMMoveType lastMoveType);
    
    // Helper methods for LM
    std::unordered_map<int, int> createNextNodeMap(Indices const& solution) const;
    std::unordered_map<int, int> createPrevNodeMap(Indices const& solution) const;
    bool edgeExists(std::unordered_map<int, int> const& nextMap, int u, int v) const;
    int computeInterRouteDelta(Indices const& solution, int selectedIndex, int unselectedNode) const;
    int computeTwoOptDelta(Indices const& solution, int i, int j) const;
    int getNodePosition(Indices const& solution, int node) const;
    void reverseSegment(Indices& solution, int start, int end);
};
