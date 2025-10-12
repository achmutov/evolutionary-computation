#include <evolutionary_computation/solver/nearest_neighbor_pos.h>
#include <limits>

std::string NearestNeighborPosSolver::name() const {
    return std::string("nearest_neighbor_pos");
}

NearestNeighborPosSolver::Indices NearestNeighborPosSolver::_solve(int i) {
    auto const n = this->data.size();
    auto const half = static_cast<Indices::size_type>(std::round(static_cast<float>(n) / 2));

    auto visited = std::vector<bool>(n, false);
    auto indices = std::vector<int>();
    indices.reserve(half);

    indices.push_back(i);
    visited[i] = true;

    for (int j = 1; j < half; j++) {
        // Identify closest node to the current path
        int closest = -1;
        int closestVisitedPos = -1;
        auto minObj = std::numeric_limits<double>::max();
        for (int iToAdd = 0; iToAdd < n; iToAdd++) {
            if (visited[iToAdd]) continue;
            for (size_t pos = 0; pos < indices.size(); pos++) {
                auto const iInPath = indices[pos];
                auto const currentObj = this->distances[iInPath][iToAdd] + this->data[iToAdd].cost;
                if (minObj > currentObj) {
                    minObj = currentObj;
                    closest = iToAdd;
                    closestVisitedPos = pos;
                }
            }
        }


        // Decide whether to insert before or after the closest visited node
        auto prevPos = (closestVisitedPos - 1 + indices.size()) % indices.size();
        auto nextPos = (closestVisitedPos + 1) % indices.size();
        auto candidatePos = std::vector<std::pair<int, int>> {
            {prevPos, closestVisitedPos},
            {closestVisitedPos, nextPos}
        };
        int closestPos = -1;
        auto minDelta = std::numeric_limits<int>::max();
        for (auto const [prevPos, pos] : candidatePos) {
            auto before = indices[prevPos];
            auto after = indices[pos];
            auto delta =
                this->data[j].cost
                + this->distances[before][j]
                + this->distances[j][after]
                - this->distances[before][after];

            if (minDelta > delta) {
                minDelta = delta;
                closestPos = pos;
            }
        }

        indices.insert(indices.begin() + closestPos, closest);
        visited[closest] = true;
    }

    return indices;
}
