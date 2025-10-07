#include <algorithm>
#include <evolutionary_computation/solver/random.h>

RandomSolver::Indices RandomSolver::_solve() {
    auto const n = this->data.size();
    auto indices = std::vector<int>(n);

    std::ranges::iota(indices, 0);
    std::ranges::shuffle(indices, this->mt);

    auto const half = static_cast<Indices::size_type>(std::round(static_cast<float>(n) / 2));
    indices.resize(half);
    return indices;
}
