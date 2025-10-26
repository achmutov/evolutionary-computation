#include <evolutionary_computation/solver/nearest_neighbor.h>

std::string NearestNeighborSolver::name() const {
    return std::string("nearest_neighbor");
}

NearestNeighborSolver::Indices NearestNeighborSolver::_solve(int i) {
    auto const n = this->data.entries.size();
    auto const half = static_cast<Indices::size_type>(std::round(static_cast<float>(n) / 2));

    auto visited = std::vector<bool>(n, false);
    auto indices = std::vector<int>();
    indices.reserve(half);

    auto& current = i;
    indices.push_back(current);
    visited[current] = true;

    for (int j = 1; j < half; j++) {
        int next = -1;
        auto minObj = std::numeric_limits<double>::max();

        for (int candidate = 0; candidate < n; candidate++) {
            auto const currentObj = this->distances[current][candidate] + this->data.entries[candidate].cost;
            if (!visited[candidate] && currentObj < minObj) {
                minObj = currentObj;
                next = candidate;
            }
        }

        current = next;
        indices.push_back(current);
        visited[current] = true;
    }

    return indices;
}
