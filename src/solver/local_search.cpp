#include <algorithm>
#include <evolutionary_computation/solver/local_search.h>
#include <limits>

LocalSearchSolver::LocalSearchSolver(
    LocalSearchType localSearchType,
    IntraNeighborhoodType intraNeighborhoodType,
    Solver &initSolver)
    : localSearchType{localSearchType},
      intraNeighborhoodType{intraNeighborhoodType},
      initSolver{initSolver} {}

std::string LocalSearchSolver::name() const {
    return std::string("local_search-")
        + ((localSearchType == LocalSearchType::Steep) ? "steep-" : "greedy-")
        + ((intraNeighborhoodType == IntraNeighborhoodType::Edges) ? "two_edges-" : "two_nodes-")
        + "init_" + initSolver.name();
}

void LocalSearchSolver::init(Data const& data) {
    Solver::init(data);
    this->initSolver.init(data);
}

LocalSearchSolver::Indices LocalSearchSolver::_solve(int i) {
    auto currentSolution = this->initSolver._solve(i);

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

    // Inter
    for (int i = 0; i < solution.size(); i++)
        for (int j = 0; j < n; j++)
            if (!visited[j])
                moves.emplace_back(i, j, MoveType::Selection);

    // Intra
    auto intraNeighborhoodType = this->intraNeighborhoodType == IntraNeighborhoodType::Edges
        ? MoveType::Edges
        : MoveType::Nodes;
    for (int i = 0; i < solution.size(); i++)
        for (int j = i + 1; j < solution.size(); j++)
            moves.emplace_back(i, j, intraNeighborhoodType);

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
