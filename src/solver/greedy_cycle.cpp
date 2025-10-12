#include <evolutionary_computation/solver/greedy_cycle.h>

std::string GreedyCycle::name() const {
    return std::string("greedy_cycle");
}

GreedyCycle::Indices GreedyCycle::_solve(int i) {
    auto const n = this->data.size();
    auto const half = static_cast<Indices::size_type>(std::round(static_cast<float>(n) / 2));

    auto visited = std::vector<bool>(n, false);
    auto indices = std::vector<int>();
    indices.reserve(half);

    indices.push_back(i);
    visited[i] = true;

    for (int k = 1; k < half; k++) {
        int current = -1;
        int currentPos = -1;
        auto minDelta = std::numeric_limits<int>::max();

        for (int j = 0; j < n; ++j) {
            if (visited[j]) continue;

            for (size_t pos = 1; pos <= indices.size(); ++pos) {
                auto before = indices[pos - 1];
                auto after = indices[pos % indices.size()];
                auto delta =
                    this->data[j].cost
                    + this->distances[before][j]
                    + this->distances[j][after]
                    - this->distances[before][after];

                if (minDelta > delta) {
                    minDelta = delta;
                    current = j;
                    currentPos = pos;
                }
            }
        }

        indices.insert(indices.begin() + currentPos, current);
        visited[current] = true;
    }

    return indices;
}
