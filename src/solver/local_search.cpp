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
        // Build set of remaining nodes (not in solution)
        std::set<int> remainingNodes;
        std::set<int> selected(currentSolution.begin(), currentSolution.end());
        for (size_t j = 0; j < this->data.entries.size(); j++) {
            if (selected.find(j) == selected.end()) {
                remainingNodes.insert(j);
            }
        }
        
        // Clear move list for new run
        improvingMoveList.clear();
        
        bool improved = true;
        while (improved) {
            improved = performSteepestStepLM(currentSolution, remainingNodes);
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

void LocalSearchSolver::doSteepLM(Indices& solution, std::set<int>& remainingNodes) {
    // This method is not used directly - performSteepestStepLM is called from _solve
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

bool LocalSearchSolver::performSteepestStepLM(Indices& solution, std::set<int>& remainingNodes) {
    // 1. Populate LM if it's empty (first iteration or after a local optimum)
    if (improvingMoveList.empty()) {
        populateMoveList(solution, remainingNodes);
        if (improvingMoveList.empty()) {
            return false; // No improving moves found at all
        }
    }

    // Build successor/predecessor maps for fast edge validation
    auto succMap = buildSuccMap(solution);
    auto predMap = buildPredMap(solution);

    // 2. Recheck moves in LM, best first
    std::sort(improvingMoveList.begin(), improvingMoveList.end(),
              [](const MoveWithDelta& a, const MoveWithDelta& b) {
                  return a.delta < b.delta;
              });

    for (auto it = improvingMoveList.begin(); it != improvingMoveList.end();) {
        const MoveWithDelta& move = *it;

        // 3. Validate the move against the current solution
        ValidationResult val = validateMove(move, succMap, predMap, remainingNodes);

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
        std::set<int> changedNodes = applyMoveLM(solution, remainingNodes, move, val);
        it = improvingMoveList.erase(it); // Remove the move we just applied

        updateLocalMoves(solution, remainingNodes, changedNodes, move.type);

        return true; // Found and applied the best move
    }

    return false; // No valid improving move found in the list
}

void LocalSearchSolver::populateMoveList(Indices const& solution, std::set<int> const& remainingNodes) {
    int n = solution.size();
    improvingMoveList.clear();

    // Inter-route moves
    for (int i = 0; i < n; i++) {
        int nodeInCycle = solution[i];
        for (int nodeOutOfCycle : remainingNodes) {
            int delta = deltaInter(solution, i, nodeOutOfCycle);
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
            int delta = deltaTwoOpt(solution, i, j);

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

ValidationResult LocalSearchSolver::validateMove(MoveWithDelta const& move,
                                                  std::unordered_map<int, int> const& succMap,
                                                  std::unordered_map<int, int> const& predMap,
                                                  std::set<int> const& remaining) const {
    switch (move.type) {
        case LMMoveType::EXCHANGE_SELECTED_UNSELECTED: {
            // A=prev, B=next, C=inCycle, D=outOfCycle
            int A_inter = move.nodeA;
            int B_inter = move.nodeB;
            int C_inter = move.nodeC;
            int D_inter = move.nodeD;

            bool edgesExist = (hasEdge(succMap, A_inter, C_inter) || hasEdge(succMap, C_inter, A_inter)) &&
                            (hasEdge(succMap, C_inter, B_inter) || hasEdge(succMap, B_inter, C_inter));
            bool d_exists = remaining.find(D_inter) != remaining.end();

            if (!edgesExist || !d_exists) {
                return ValidationResult(false, false, false); // Case 1: Remove
            }

            // Check for forward direction: (A -> C -> B)
            if (hasEdge(succMap, A_inter, C_inter) && hasEdge(succMap, C_inter, B_inter)) {
                return ValidationResult(true, true, false); // Case 3: Apply forward
            }
            // Check for reversed direction: (B -> C -> A)
            if (hasEdge(succMap, B_inter, C_inter) && hasEdge(succMap, C_inter, A_inter)) {
                return ValidationResult(true, true, true); // Case 3: Apply reversed
            }

            return ValidationResult(false, true, false); // Case 2: Keep
        }

        case LMMoveType::TWO_OPT: {
            // A=a, B=b, C=c, D=d
            int A_intra = move.nodeA;
            int B_intra = move.nodeB;
            int C_intra = move.nodeC;
            int D_intra = move.nodeD;

            bool ab_exists = hasEdge(succMap, A_intra, B_intra) || hasEdge(succMap, B_intra, A_intra);
            bool cd_exists = hasEdge(succMap, C_intra, D_intra) || hasEdge(succMap, D_intra, C_intra);

            if (!ab_exists || !cd_exists) {
                return ValidationResult(false, false, false); // Case 1: Remove
            }

            // Check for forward direction: (A -> B) and (C -> D)
            if (hasEdge(succMap, A_intra, B_intra) && hasEdge(succMap, C_intra, D_intra)) {
                return ValidationResult(true, true, false); // Case 3: Apply forward
            }
            // Check for reversed direction: (B -> A) and (D -> C)
            if (hasEdge(succMap, B_intra, A_intra) && hasEdge(succMap, D_intra, C_intra)) {
                return ValidationResult(true, true, true); // Case 3: Apply reversed
            }

            return ValidationResult(false, true, false); // Case 2: Keep
        }
    }
    return ValidationResult(false, false, false);
}

std::set<int> LocalSearchSolver::applyMoveLM(Indices& solution, std::set<int>& remaining,
                                             MoveWithDelta const& move, ValidationResult const& val) {
    std::set<int> changedNodes;

    switch (move.type) {
        case LMMoveType::EXCHANGE_SELECTED_UNSELECTED: {
            // A=prev, B=next, C=inCycle, D=outOfCycle
            int nodeC = move.nodeC;
            int nodeD = move.nodeD;
            int idxC = findNodeIndex(solution, nodeC);

            if (idxC != -1) {
                solution[idxC] = nodeD;
                remaining.erase(nodeD);
                remaining.insert(nodeC);

                changedNodes.insert(move.nodeA); // prev
                changedNodes.insert(move.nodeB); // next
                changedNodes.insert(nodeC);      // old node (now in remaining)
                changedNodes.insert(nodeD);      // new node (now in route)
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
                startIdx = findNodeIndex(solution, nodeB);
                endIdx = findNodeIndex(solution, nodeC);
            } else {
                // Reversed: B->A, D->C. Reverse (A...D)
                startIdx = findNodeIndex(solution, nodeA);
                endIdx = findNodeIndex(solution, nodeD);
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
                reverseSublist(solution, startIdx, endIdx);
            }

            changedNodes.insert(nodeA);
            changedNodes.insert(nodeB);
            changedNodes.insert(nodeC);
            changedNodes.insert(nodeD);
            break;
        }
    }
    return changedNodes;
}

void LocalSearchSolver::updateLocalMoves(Indices const& solution, std::set<int> const& remaining,
                                         std::set<int> const& changedNodes, LMMoveType lastMoveType) {
    int n = solution.size();

    // Remove moves involving changed nodes
    improvingMoveList.erase(
        std::remove_if(improvingMoveList.begin(), improvingMoveList.end(),
            [&changedNodes](const MoveWithDelta& m) {
                return changedNodes.find(m.nodeA) != changedNodes.end() ||
                       changedNodes.find(m.nodeB) != changedNodes.end() ||
                       changedNodes.find(m.nodeC) != changedNodes.end() ||
                       changedNodes.find(m.nodeD) != changedNodes.end();
            }),
        improvingMoveList.end()
    );

    std::set<int> indicesToCheck;
    for (int node : changedNodes) {
        int idx = findNodeIndex(solution, node); // Check nodes now in the route
        if (idx != -1) {
            indicesToCheck.insert((idx - 1 + n) % n); // neighbor
            indicesToCheck.insert(idx);               // self
            indicesToCheck.insert((idx + 1) % n);    // neighbor
        }
    }

    // A. Re-evaluate neighbors of affected positions with all remaining nodes
    for (int i : indicesToCheck) {
        int nodeInCycle = solution[i];
        for (int nodeOutOfCycle : remaining) {
            int delta = deltaInter(solution, i, nodeOutOfCycle);
            if (delta < 0) {
                int prev = solution[(i - 1 + n) % n];
                int next = solution[(i + 1) % n];
                improvingMoveList.push_back(MoveWithDelta::forInterRoute(
                    delta, prev, next, nodeInCycle, nodeOutOfCycle
                ));
            }
        }
    }

    // B. Re-evaluate all cycle nodes with the newly available remaining node
    if (lastMoveType == LMMoveType::EXCHANGE_SELECTED_UNSELECTED) {
        // Find the node that was just moved from the route to 'remaining'
        int newNodeInRemaining = -1;
        for (int node : changedNodes) {
            if (remaining.find(node) != remaining.end()) {
                newNodeInRemaining = node;
                break;
            }
        }

        if (newNodeInRemaining != -1) {
            for (int i = 0; i < n; i++) {
                int nodeInCycle = solution[i];
                int delta = deltaInter(solution, i, newNodeInRemaining);
                if (delta < 0) {
                    int prev = solution[(i - 1 + n) % n];
                    int next = solution[(i + 1) % n];
                    improvingMoveList.push_back(MoveWithDelta::forInterRoute(
                        delta, prev, next, nodeInCycle, newNodeInRemaining
                    ));
                }
            }
        }
    }

    // Find all nodes in the cycle that were changed (A, B, C, D)
    std::set<int> changedCycleNodes;
    for (int node : changedNodes) {
        if (findNodeIndex(solution, node) != -1) {
            changedCycleNodes.insert(node);
        }
    }

    // Check all new pairs (A', C') where A' is a changed node
    for (int nodeA_prime : changedCycleNodes) {
        int i = findNodeIndex(solution, nodeA_prime);
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
            int delta = deltaTwoOpt(solution, idx_i, idx_j);

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
std::unordered_map<int, int> LocalSearchSolver::buildSuccMap(Indices const& solution) const {
    std::unordered_map<int, int> map;
    for (size_t i = 0; i < solution.size(); i++) {
        map[solution[i]] = solution[(i + 1) % solution.size()];
    }
    return map;
}

std::unordered_map<int, int> LocalSearchSolver::buildPredMap(Indices const& solution) const {
    std::unordered_map<int, int> map;
    int n = solution.size();
    for (int i = 0; i < n; i++) {
        map[solution[i]] = solution[(i - 1 + n) % n];
    }
    return map;
}

bool LocalSearchSolver::hasEdge(std::unordered_map<int, int> const& succMap, int u, int v) const {
    auto it = succMap.find(u);
    return it != succMap.end() && it->second == v;
}

int LocalSearchSolver::deltaInter(Indices const& solution, int selectedIndex, int unselectedNode) const {
    int selectedNode = solution[selectedIndex];
    int prev = solution[(selectedIndex - 1 + solution.size()) % solution.size()];
    int next = solution[(selectedIndex + 1) % solution.size()];

    int before = distances[prev][selectedNode] + distances[selectedNode][next];
    int after = distances[prev][unselectedNode] + distances[unselectedNode][next];

    before += data.entries[selectedNode].cost;
    after += data.entries[unselectedNode].cost;

    return after - before;
}

int LocalSearchSolver::deltaTwoOpt(Indices const& solution, int i, int j) const {
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

int LocalSearchSolver::findNodeIndex(Indices const& solution, int node) const {
    for (size_t i = 0; i < solution.size(); i++) {
        if (solution[i] == node) {
            return i;
        }
    }
    return -1;
}

void LocalSearchSolver::reverseSublist(Indices& solution, int start, int end) {
    while (start < end) {
        std::swap(solution[start++], solution[end--]);
    }
}
