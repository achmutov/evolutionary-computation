#include <algorithm>
#include <numeric>
#include <evolutionary_computation/solver/random.h>

std::string RandomSolver::name() const {
    return std::string("random");
}

RandomSolver::Indices RandomSolver::_solve(int i) {
    auto const n = this->data.entries.size();
    auto indices = std::vector<int>(n);

    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), this->mt);

    auto const half = static_cast<Indices::size_type>(std::round(static_cast<float>(n) / 2));
    indices.resize(half);
    return indices;
}
