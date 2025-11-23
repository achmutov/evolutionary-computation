#include <algorithm>
#include <evolutionary_computation/solver/local_search.h>
#include <limits>
#include <unordered_set>

LocalSearchSolver::LocalSearchSolver(
    LocalSearchType localSearchType,
    IntraNeighborhoodType intraNeighborhoodType,
    Solver &initSolver,
    int nCandidateMoves,
    bool useMoveList)
    : localSearchType{localSearchType},
      intraNeighborhoodType{intraNeighborhoodType},
      initSolver{initSolver},
      nCandidateMoves{nCandidateMoves},
      useMoveList{useMoveList} {}

std::string LocalSearchSolver::name() const {
    return std::string("local_search-")
        + (this->localSearchType == LocalSearchType::Steep ? "steep-" : "greedy-")
        + (this->intraNeighborhoodType == IntraNeighborhoodType::Edges ? "two_edges-" : "two_nodes-")
        + "init_" + initSolver.name()
        + (this->nCandidateMoves == -1 ? "" : "-candidate_moves")
        + (this->useMoveList ? "-LM" : "");
}

void LocalSearchSolver::init(Data const& data) {
    Solver::init(data);
    this->learnCandidateMoves();
    this->initSolver.init(data);
}

void LocalSearchSolver::learnCandidateMoves() {
    if (this->nCandidateMoves == -1) return;

    auto const n = this->data.entries.size();
    auto indices = Indices(n);
    std::iota(indices.begin(), indices.end(), 0);

    this->closestCities = std::vector<Indices>(n);
    this->isClosestCity = std::vector<std::vector<bool>>(n);

    Indices currentIndices;
    std::vector<Move>* interMoves;
    std::vector<Move>* intraMoves;
    for (auto i : indices) {
        currentIndices = indices;
        std::nth_element(
            currentIndices.begin(),
            currentIndices.begin() + this->nCandidateMoves,
            currentIndices.end(),
            [&](auto a, auto b) {
                auto aScore = this->data.entries[a].cost + this->distances[i][a];
                auto bScore = this->data.entries[b].cost + this->distances[i][b];
                return aScore < bScore;
            }
        );
        currentIndices.resize(this->nCandidateMoves);
        this->closestCities[i] = currentIndices;
        this->isClosestCity[i] = std::vector<bool>(n, false);
        for (auto j : currentIndices) this->isClosestCity[i][j] = true;
    }
}

LocalSearchSolver::Indices LocalSearchSolver::_solve(int i) {
    auto currentSolution = this->initSolver._solve(i);

    if (this->useMoveList && this->localSearchType == LocalSearchType::Steep) {
        // Build set of unselected nodes (not in solution)
        std::set<int> unselectedNodes;
        std::set<int> selected(currentSolution.begin(), currentSolution.end());
        for (size_t j = 0; j < this->data.entries.size(); j++) {
            if (selected.find(j) == selected.end()) {
                unselectedNodes.insert(j);
            }
        }
        
        // Clear move list for new run
        improvingMoveList.clear();
        
        bool improved = true;
        while (improved) {
            improved = executeSteepestStep(currentSolution, unselectedNodes);
        }
        
        return currentSolution;
    }

    int currentCost = std::numeric_limits<int>::max();
    int lastCost;

    do {
        lastCost = currentCost;

        switch (this->localSearchType) {
            case LocalSearchType::Greedy:
                this->doGreedy(currentSolution);
                break;
            case LocalSearchType::Steep:
                this->doSteep(currentSolution);
                break;
        }

        currentCost = this->cost(currentSolution);
    } while (lastCost > currentCost);

    return currentSolution;
}

void LocalSearchSolver::doGreedy(Indices& solution) {
    for (auto& move : this->getMoves(solution)) {
        if (this->evalMove(solution, move) < 0) {
            this->applyMove(solution, move);
            return;
        }
    }
}

void LocalSearchSolver::doSteep(Indices& solution) {
    // This is the baseline method - always use it when useMoveList is false
    auto moves = this->getMoves(solution);
    if (!moves.size()) return;

    Move bestMove;
    int bestDelta = std::numeric_limits<int>::max();

    for (auto& move : moves) {
        auto delta = this->evalMove(solution, move);
        if (bestDelta > delta) {
            bestDelta = delta;
            bestMove = move;
        }
    }
    if (bestDelta < 0) this->applyMove(solution, bestMove);
}

void LocalSearchSolver::doSteepLM(Indices& solution, std::set<int>& unselectedNodes) {
    // This method is not used directly - executeSteepestStep is called from _solve
    // Keeping for potential future use
}

std::vector<LocalSearchSolver::Move> LocalSearchSolver::getMoves(Indices const& solution) {
    auto n = this->data.entries.size();
    auto sn = solution.size();
    auto moves = std::vector<Move>();
    moves.reserve(
        sn * (n - sn) + // Inter
        sn * sn         // Intra
    );

    auto visited = std::vector(n, false);
    for (auto i : solution) visited[i] = true;

    if (nCandidateMoves == -1) {
        // Inter
        for (int i = 0; i < sn; i++)
            for (int j = 0; j < n; j++)
                if (!visited[j])
                    moves.emplace_back(i, j, MoveType::Selection);

        // Intra
        auto intraNeighborhoodType = this->intraNeighborhoodType == IntraNeighborhoodType::Edges
            ? MoveType::Edges
            : MoveType::Nodes;
        for (int i = 0; i < sn; i++)
            for (int j = i + 1; j < sn; j++)
                moves.emplace_back(i, j, intraNeighborhoodType);
    } else {
        // Inter
        for (int i = 0; i < sn; i++)
            for (auto j : this->closestCities[solution[i]])
                if (!visited[j])
                    moves.emplace_back((i + sn - 1) % sn, j, MoveType::Selection);

        // Infra
        auto intraNeighborhoodType = this->intraNeighborhoodType == IntraNeighborhoodType::Edges
            ? MoveType::Edges
            : MoveType::Nodes;

        for (int i = 0; i < sn; i++) {
            for (int j = i + 1; j < sn; j++)
                if (this->isClosestCity[solution[i]][solution[i != 0 ? i - 1 : sn - 1]] ||
                        this->isClosestCity[solution[j]][solution[j != sn - 1 ? j + 1 : 0]])
                    moves.emplace_back(i, j, intraNeighborhoodType);
        }
    }

    std::shuffle(moves.begin(), moves.end(), this->mt);
    return moves;
}

int LocalSearchSolver::evalMove(Indices const& solution, Move& move) {
    auto [a, b, moveType] = move;
    switch (moveType) {
        case MoveType::Selection: {
            auto before = solution[(a + solution.size() - 1) % solution.size()];
            auto after = solution[(a + 1) % solution.size()];
            return this->getDelta(b, before, after) - this->getDelta(solution[a], before, after);
        }
        case MoveType::Nodes: {
            auto first = solution[a];
            auto beforeFirst = solution[(a + solution.size() - 1) % solution.size()];
            auto afterFirst = solution[(a + 1) % solution.size()];

            auto second = solution[b];
            auto beforeSecond = solution[(b + solution.size() - 1) % solution.size()];
            auto afterSecond = solution[(b + 1) % solution.size()];

            if (afterFirst == second) {
                auto newEdges = this->distances[beforeFirst][second] + this->distances[first][afterSecond];
                auto oldEdges = this->distances[beforeFirst][first] + this->distances[second][afterSecond];
                return newEdges - oldEdges;
            }
            else if (second == beforeFirst) {
                auto oldEdges = this->distances[beforeFirst][second] + this->distances[first][afterSecond];
                auto newEdges = this->distances[beforeFirst][first] + this->distances[second][afterSecond];
                return newEdges - oldEdges;
            } else {
                auto firstOld = this->distances[beforeFirst][first] + this->distances[first][afterFirst];
                auto secondOld = this->distances[beforeSecond][second] + this->distances[second][afterSecond];

                auto firstNew = this->distances[beforeSecond][first] + this->distances[first][afterSecond];
                auto secondNew = this->distances[beforeFirst][second] + this->distances[second][afterFirst];

                return firstNew + secondNew - firstOld - secondOld;
            }
        }
        case MoveType::Edges: {
            if (a == 0 && b == solution.size() - 1) return 0;
            auto first = solution[a];
            auto second = solution[b];
            auto beforeFirst = solution[(a + solution.size() - 1) % solution.size()];
            auto afterSecond = solution[(b + 1) % solution.size()];
            return this->distances[first][afterSecond] - this->distances[first][beforeFirst]
                + this->distances[second][beforeFirst] - this->distances[second][afterSecond];
        }
    }
    return 0;
}

int LocalSearchSolver::getDelta(int target, int city1, int city2) {
    return this->data.entries[target].cost
        + this->distances[city1][target]
        + this->distances[target][city2];
}

void LocalSearchSolver::applyMove(Indices& solution, Move& move) {
    auto [a, b, moveType] = move;
    switch (moveType) {
        case MoveType::Selection:
            solution[a] = b;
            break;
        case MoveType::Nodes:
            std::swap(solution[a], solution[b]);
            break;
        case MoveType::Edges:
            std::reverse(solution.begin() + a, solution.begin() + b + 1);
            break;
    }
}

// ========== LM (List of Moves) Implementation ==========

bool LocalSearchSolver::executeSteepestStep(Indices& solution, std::set<int>& unselectedNodes) {
    // 1. Populate LM if it's empty (first iteration or after a local optimum)
    if (improvingMoveList.empty()) {
        initializeMoveList(solution, unselectedNodes);
        if (improvingMoveList.empty()) {
            return false; // No improving moves found at all
        }
    }

    // Build successor/predecessor maps for fast edge validation
    auto nextMap = createNextNodeMap(solution);
    auto prevMap = createPrevNodeMap(solution);

    // 2. Recheck moves in LM, best first
    std::sort(improvingMoveList.begin(), improvingMoveList.end(),
              [](const MoveWithDelta& a, const MoveWithDelta& b) {
                  return a.delta < b.delta;
              });

    for (auto it = improvingMoveList.begin(); it != improvingMoveList.end();) {
        const MoveWithDelta& move = *it;

        // 3. Validate the move against the current solution
        ValidationResult val = checkMoveValidity(move, nextMap, prevMap, unselectedNodes);

        if (!val.keepMove) {
            it = improvingMoveList.erase(it); // Remove move, edges no longer exist
            continue;
        }

        if (!val.applyMove) {
            // Keep move, but don't apply (e.g., direction is wrong)
            ++it;
            continue;
        }

        // 4. Apply the move
        std::set<int> affectedNodes = executeMove(solution, unselectedNodes, move, val);
        it = improvingMoveList.erase(it); // Remove the move we just applied

        refreshMoveList(solution, unselectedNodes, affectedNodes, move.type);

        return true; // Found and applied the best move
    }

    return false; // No valid improving move found in the list
}

void LocalSearchSolver::initializeMoveList(Indices const& solution, std::set<int> const& unselectedNodes) {
    int n = solution.size();
    improvingMoveList.clear();

    // Inter-route moves
    for (int i = 0; i < n; i++) {
        int nodeInCycle = solution[i];
        for (int nodeOutOfCycle : unselectedNodes) {
            int delta = computeInterRouteDelta(solution, i, nodeOutOfCycle);
            if (delta < 0) {
                int prev = solution[(i - 1 + n) % n];
                int next = solution[(i + 1) % n];
                improvingMoveList.push_back(MoveWithDelta::forInterRoute(
                    delta, prev, next, nodeInCycle, nodeOutOfCycle
                ));
            }
        }
    }

    // Intra-route moves (2-opt)
    for (int i = 0; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++) {
            int delta = computeTwoOptDelta(solution, i, j);

            if (delta < 0) {
                int a = solution[i];
                int b = solution[(i + 1) % n];
                int c = solution[j];
                int d = solution[(j + 1) % n];

                if (c == b) continue; // Avoid adjacent/overlapping

                improvingMoveList.push_back(MoveWithDelta::forIntraEdge(delta, a, b, c, d));
            }
        }
    }

    std::sort(improvingMoveList.begin(), improvingMoveList.end(),
              [](const MoveWithDelta& a, const MoveWithDelta& b) {
                  return a.delta < b.delta;
              });
}

ValidationResult LocalSearchSolver::checkMoveValidity(MoveWithDelta const& move,
                                                  std::unordered_map<int, int> const& nextMap,
                                                  std::unordered_map<int, int> const& prevMap,
                                                  std::set<int> const& unselected) const {
    switch (move.type) {
        case LMMoveType::EXCHANGE_SELECTED_UNSELECTED: {
            // A=prev, B=next, C=inCycle, D=outOfCycle
            int prevNode = move.nodeA;
            int nextNode = move.nodeB;
            int selectedNode = move.nodeC;
            int unselectedNode = move.nodeD;

            bool edgesExist = (edgeExists(nextMap, prevNode, selectedNode) || edgeExists(nextMap, selectedNode, prevNode)) &&
                            (edgeExists(nextMap, selectedNode, nextNode) || edgeExists(nextMap, nextNode, selectedNode));
            bool d_exists = unselected.find(unselectedNode) != unselected.end();

            if (!edgesExist || !d_exists) {
                return ValidationResult(false, false, false); // Case 1: Remove
            }

            // Check for forward direction: (A -> C -> B)
            if (edgeExists(nextMap, prevNode, selectedNode) && edgeExists(nextMap, selectedNode, nextNode)) {
                return ValidationResult(true, true, false); // Case 3: Apply forward
            }
            // Check for reversed direction: (B -> C -> A)
            if (edgeExists(nextMap, nextNode, selectedNode) && edgeExists(nextMap, selectedNode, prevNode)) {
                return ValidationResult(true, true, true); // Case 3: Apply reversed
            }

            return ValidationResult(false, true, false); // Case 2: Keep
        }

        case LMMoveType::TWO_OPT: {
            // A=a, B=b, C=c, D=d
            int firstNode = move.nodeA;
            int secondNode = move.nodeB;
            int thirdNode = move.nodeC;
            int fourthNode = move.nodeD;

            bool ab_exists = edgeExists(nextMap, firstNode, secondNode) || edgeExists(nextMap, secondNode, firstNode);
            bool cd_exists = edgeExists(nextMap, thirdNode, fourthNode) || edgeExists(nextMap, fourthNode, thirdNode);

            if (!ab_exists || !cd_exists) {
                return ValidationResult(false, false, false); // Case 1: Remove
            }

            // Check for forward direction: (A -> B) and (C -> D)
            if (edgeExists(nextMap, firstNode, secondNode) && edgeExists(nextMap, thirdNode, fourthNode)) {
                return ValidationResult(true, true, false); // Case 3: Apply forward
            }
            // Check for reversed direction: (B -> A) and (D -> C)
            if (edgeExists(nextMap, secondNode, firstNode) && edgeExists(nextMap, fourthNode, thirdNode)) {
                return ValidationResult(true, true, true); // Case 3: Apply reversed
            }

            return ValidationResult(false, true, false); // Case 2: Keep
        }
    }
    return ValidationResult(false, false, false);
}

std::set<int> LocalSearchSolver::executeMove(Indices& solution, std::set<int>& unselected,
                                             MoveWithDelta const& move, ValidationResult const& val) {
    std::set<int> affectedNodes;

    switch (move.type) {
        case LMMoveType::EXCHANGE_SELECTED_UNSELECTED: {
            // A=prev, B=next, C=inCycle, D=outOfCycle
            int oldNode = move.nodeC;
            int newNode = move.nodeD;
            int replaceIndex = getNodePosition(solution, oldNode);

            if (replaceIndex != -1) {
                solution[replaceIndex] = newNode;
                unselected.erase(newNode);
                unselected.insert(oldNode);

                affectedNodes.insert(move.nodeA); // prev
                affectedNodes.insert(move.nodeB); // next
                affectedNodes.insert(oldNode);      // old node (now in remaining)
                affectedNodes.insert(newNode);      // new node (now in route)
            }
            break;
        }

        case LMMoveType::TWO_OPT: {
            // A=a, B=b, C=c, D=d
            int nodeA = move.nodeA;
            int nodeB = move.nodeB;
            int nodeC = move.nodeC;
            int nodeD = move.nodeD;

            int startIdx, endIdx;

            if (!val.applyReversed) {
                // Forward: A->B, C->D. Reverse (B...C)
                startIdx = getNodePosition(solution, nodeB);
                endIdx = getNodePosition(solution, nodeC);
            } else {
                // Reversed: B->A, D->C. Reverse (A...D)
                startIdx = getNodePosition(solution, nodeA);
                endIdx = getNodePosition(solution, nodeD);
            }

            // Handle wrap-around
            if (startIdx > endIdx) {
                Indices sublist;
                for (int i = startIdx; i < solution.size(); i++) {
                    sublist.push_back(solution[i]);
                }
                for (int i = 0; i <= endIdx; i++) {
                    sublist.push_back(solution[i]);
                }

                std::reverse(sublist.begin(), sublist.end());

                int k = 0;
                for (int i = startIdx; i < solution.size(); i++) {
                    solution[i] = sublist[k++];
                }
                for (int i = 0; i <= endIdx; i++) {
                    solution[i] = sublist[k++];
                }
            } else {
                reverseSegment(solution, startIdx, endIdx);
            }

            affectedNodes.insert(nodeA);
            affectedNodes.insert(nodeB);
            affectedNodes.insert(nodeC);
            affectedNodes.insert(nodeD);
            break;
        }
    }
    return affectedNodes;
}

void LocalSearchSolver::refreshMoveList(Indices const& solution, std::set<int> const& unselected,
                                         std::set<int> const& affectedNodes, LMMoveType lastMoveType) {
    int n = solution.size();

    // Remove moves involving affected nodes
    improvingMoveList.erase(
        std::remove_if(improvingMoveList.begin(), improvingMoveList.end(),
            [&affectedNodes](const MoveWithDelta& m) {
                return affectedNodes.find(m.nodeA) != affectedNodes.end() ||
                       affectedNodes.find(m.nodeB) != affectedNodes.end() ||
                       affectedNodes.find(m.nodeC) != affectedNodes.end() ||
                       affectedNodes.find(m.nodeD) != affectedNodes.end();
            }),
        improvingMoveList.end()
    );

    std::set<int> indicesToCheck;
    for (int node : affectedNodes) {
        int idx = getNodePosition(solution, node); // Check nodes now in the route
        if (idx != -1) {
            indicesToCheck.insert((idx - 1 + n) % n); // neighbor
            indicesToCheck.insert(idx);               // self
            indicesToCheck.insert((idx + 1) % n);    // neighbor
        }
    }

    // A. Re-evaluate neighbors of affected positions with all unselected nodes
    for (int i : indicesToCheck) {
        int nodeInCycle = solution[i];
        for (int nodeOutOfCycle : unselected) {
            int delta = computeInterRouteDelta(solution, i, nodeOutOfCycle);
            if (delta < 0) {
                int prev = solution[(i - 1 + n) % n];
                int next = solution[(i + 1) % n];
                improvingMoveList.push_back(MoveWithDelta::forInterRoute(
                    delta, prev, next, nodeInCycle, nodeOutOfCycle
                ));
            }
        }
    }

    // B. Re-evaluate all cycle nodes with the newly available unselected node
    if (lastMoveType == LMMoveType::EXCHANGE_SELECTED_UNSELECTED) {
        // Find the node that was just moved from the route to 'unselected'
        int newlyAvailableNode = -1;
        for (int node : affectedNodes) {
            if (unselected.find(node) != unselected.end()) {
                newlyAvailableNode = node;
                break;
            }
        }

        if (newlyAvailableNode != -1) {
            for (int i = 0; i < n; i++) {
                int nodeInCycle = solution[i];
                int delta = computeInterRouteDelta(solution, i, newlyAvailableNode);
                if (delta < 0) {
                    int prev = solution[(i - 1 + n) % n];
                    int next = solution[(i + 1) % n];
                    improvingMoveList.push_back(MoveWithDelta::forInterRoute(
                        delta, prev, next, nodeInCycle, newlyAvailableNode
                    ));
                }
            }
        }
    }

    // Find all nodes in the cycle that were affected (A, B, C, D)
    std::set<int> affectedCycleNodes;
    for (int node : affectedNodes) {
        if (getNodePosition(solution, node) != -1) {
            affectedCycleNodes.insert(node);
        }
    }

    // Check all new pairs (A', C') where A' is an affected node
    for (int changedCycleNode : affectedCycleNodes) {
        int i = getNodePosition(solution, changedCycleNode);
        if (i == -1) continue;

        for (int j = 0; j < n; j++) {
            if (i == j) continue;

            // Ensure i < j for consistent delta calculation
            int idx_i = i, idx_j = j;
            if (idx_i > idx_j) {
                idx_i = j;
                idx_j = i;
            }

            // This delta calculation is for 2-opt
            int delta = computeTwoOptDelta(solution, idx_i, idx_j);

            if (delta < 0) {
                int a = solution[idx_i];
                int b = solution[(idx_i + 1) % n];
                int c = solution[idx_j];
                int d = solution[(idx_j + 1) % n];

                if (c == b) continue;

                improvingMoveList.push_back(MoveWithDelta::forIntraEdge(delta, a, b, c, d));
            }
        }
    }
}

// Helper methods for LM
std::unordered_map<int, int> LocalSearchSolver::createNextNodeMap(Indices const& solution) const {
    std::unordered_map<int, int> map;
    for (size_t i = 0; i < solution.size(); i++) {
        map[solution[i]] = solution[(i + 1) % solution.size()];
    }
    return map;
}

std::unordered_map<int, int> LocalSearchSolver::createPrevNodeMap(Indices const& solution) const {
    std::unordered_map<int, int> map;
    int n = solution.size();
    for (int i = 0; i < n; i++) {
        map[solution[i]] = solution[(i - 1 + n) % n];
    }
    return map;
}

bool LocalSearchSolver::edgeExists(std::unordered_map<int, int> const& nextMap, int u, int v) const {
    auto it = nextMap.find(u);
    return it != nextMap.end() && it->second == v;
}

int LocalSearchSolver::computeInterRouteDelta(Indices const& solution, int selectedIndex, int unselectedNode) const {
    int selectedNode = solution[selectedIndex];
    int prev = solution[(selectedIndex - 1 + solution.size()) % solution.size()];
    int next = solution[(selectedIndex + 1) % solution.size()];

    int before = distances[prev][selectedNode] + distances[selectedNode][next];
    int after = distances[prev][unselectedNode] + distances[unselectedNode][next];

    before += data.entries[selectedNode].cost;
    after += data.entries[unselectedNode].cost;

    return after - before;
}

int LocalSearchSolver::computeTwoOptDelta(Indices const& solution, int i, int j) const {
    int n = solution.size();
    // Ensure i < j
    if (i > j) {
        int temp = i;
        i = j;
        j = temp;
    }

    int a = solution[i];
    int b = solution[(i + 1) % n];
    int c = solution[j];
    int d = solution[(j + 1) % n];

    int before = distances[a][b] + distances[c][d];
    int after = distances[a][c] + distances[b][d];
    return after - before;
}

int LocalSearchSolver::getNodePosition(Indices const& solution, int node) const {
    for (size_t i = 0; i < solution.size(); i++) {
        if (solution[i] == node) {
            return i;
        }
    }
    return -1;
}

void LocalSearchSolver::reverseSegment(Indices& solution, int start, int end) {
    while (start < end) {
        std::swap(solution[start++], solution[end--]);
    }
}

// Method to run local search from a given initial solution (for ILS)
LocalSearchSolver::Indices LocalSearchSolver::_solveFromSolution(Indices const& initialSolution) {
    auto currentSolution = initialSolution;

    if (this->useMoveList && this->localSearchType == LocalSearchType::Steep) {
        // Build set of unselected nodes (not in solution)
        std::set<int> unselectedNodes;
        std::set<int> selected(currentSolution.begin(), currentSolution.end());
        for (size_t j = 0; j < this->data.entries.size(); j++) {
            if (selected.find(j) == selected.end()) {
                unselectedNodes.insert(j);
            }
        }
        
        // Clear move list for new run
        improvingMoveList.clear();
        
        bool improved = true;
        while (improved) {
            improved = executeSteepestStep(currentSolution, unselectedNodes);
        }
        
        return currentSolution;
    }

    int currentCost = std::numeric_limits<int>::max();
    int lastCost;

    do {
        lastCost = currentCost;

        switch (this->localSearchType) {
            case LocalSearchType::Greedy:
                this->doGreedy(currentSolution);
                break;
            case LocalSearchType::Steep:
                this->doSteep(currentSolution);
                break;
        }

        currentCost = this->cost(currentSolution);
    } while (lastCost > currentCost);

    return currentSolution;
}
