#include <evolutionary_computation/solver/greedy_cycle.h>

GreedyCycle::GreedyCycle(Mode mode, double alpha, double beta) 
    : mode {mode}, alpha {alpha}, beta {beta} {}

std::string GreedyCycle::name() const {
    auto result = std::string("greedy_cycle");
    if (this->mode == Mode::Regret)
        result += "_regret";
    else if (this->mode == Mode::WeightedRegret)
        result += "_weighted_regret";
    return result;
}

GreedyCycle::Indices GreedyCycle::_solve(int i) {
    auto const n = this->data.size();
    auto const half = static_cast<Indices::size_type>(std::round(static_cast<float>(n) / 2));

    auto visited = std::vector<bool>(n, false);
    auto indices = std::vector<int>();
    indices.reserve(half);

    indices.push_back(i);
    visited[i] = true;

    auto next = &GreedyCycle::next_greedy;
    if (this->mode == Mode::Regret)
        next = &GreedyCycle::next_regret;
    else if (this->mode == Mode::WeightedRegret)
        next = &GreedyCycle::next_weighted_regret;

    for (int k = 1; k < half; k++) {
        auto [bestCity, bestPos] = (*this.*next)(indices, visited);
        indices.insert(indices.begin() + bestPos, bestCity);
        visited[bestCity] = true;
    }

    return indices;
}

int GreedyCycle::get_delta(int city, int pos, Indices& indices) {
    auto before = indices[pos - 1];
    auto after = indices[pos % indices.size()];
    return this->data[city].cost
        + this->distances[before][city]
        + this->distances[city][after]
        - this->distances[before][after];
}


std::tuple<int, int> GreedyCycle::next_greedy(Indices& indices, std::vector<bool>& visited) {
    int bestCity = -1;
    int bestPos = -1;
    auto minDelta = std::numeric_limits<int>::max();

    for (int j = 0; j < visited.size(); ++j) {
        if (visited[j]) continue;

        for (size_t pos = 1; pos <= indices.size(); ++pos) {
            auto delta = this->get_delta(j, pos, indices);
            if (minDelta > delta) {
                minDelta = delta;
                bestCity = j;
                bestPos = pos;
            }
        }
    }

    return { bestCity, bestPos };
}

std::tuple<int, int, int> GreedyCycle::get_regret(int city, Indices& indices) {
    int minDelta1 = std::numeric_limits<int>::max();
    int minDelta2 = std::numeric_limits<int>::max();
    int bestPos = -1;

    for (size_t pos = 1; pos <= indices.size(); ++pos) {
        auto delta = this->get_delta(city, pos, indices);
        if (minDelta1 > delta) {
            minDelta2 = minDelta1;
            minDelta1 = delta;
            bestPos = pos;
        } else if (minDelta2 > delta) {
            minDelta2 = delta;
        }
    }
    auto regret = minDelta2 - minDelta1;
    return { minDelta1, bestPos, regret };
}

std::tuple<int, int> GreedyCycle::next_regret(Indices& indices, std::vector<bool>& visited) {
    int bestCity = -1;
    int bestPos = -1;
    int maxRegret = -1;
    int minCostForMaxRegretCity = std::numeric_limits<int>::max();

    for (int j = 0; j < visited.size(); ++j) {
        if (visited[j]) continue;
        auto [currentMinDelta, currentBestPos, currentRegret] = this->get_regret(j, indices);
        if (currentRegret > maxRegret || (currentRegret == maxRegret && currentMinDelta < minCostForMaxRegretCity)) {
            maxRegret = currentRegret;
            bestCity = j;
            bestPos = currentBestPos;
            minCostForMaxRegretCity = currentMinDelta;
        }
    }

    return { bestCity, bestPos };
}

std::tuple<int, int> GreedyCycle::next_weighted_regret(Indices& indices, std::vector<bool>& visited) {
    int bestCity = -1;
    int bestPos = -1;
    int maxScore = std::numeric_limits<int>::min();
    auto minDeltaForMaxRegret = std::numeric_limits<int>::max();

    for (int j = 0; j < visited.size(); ++j) {
        if (visited[j]) continue;
        auto [currentMinDelta, currentBestPos, currentRegret] = this->get_regret(j, indices);
        // Weighted regret: α * regret + β * (-best_cost)
        auto score = this->alpha * currentRegret + this->beta * (-currentMinDelta);
        if (score > maxScore || (score == maxScore && currentMinDelta < minDeltaForMaxRegret)) {
            maxScore = score;
            bestCity = j;
            bestPos = currentBestPos;
            minDeltaForMaxRegret = currentMinDelta;
        }
    }

    return { bestCity, bestPos };
}
