#include <evolutionary_computation/solver/base.h>
#include <random>
#include <ranges>
#include <tuple>
#include <vector>

Solver::Solver(std::vector<Data> data) : data{data}, distances{toMatrix(data)} {}

std::tuple<Solver::Indices, int> Solver::solve() {
    auto indices = this->_solve();
    auto cost = this->cost(indices);
    return std::make_tuple(indices, cost);
}

std::vector<std::vector<int>> Solver::toMatrix(std::vector<Data> const& data) {
    auto result = std::vector<std::vector<int>>(data.size());

    for (auto const [i, a] : data | std::ranges::views::enumerate) {
        result[i] = std::vector<int>(data.size());
        for (auto const [j, b] : data | std::ranges::views::enumerate) {
            result[i][j] = euclideanDistance(a.x, a.y, b.x, b.y);
        }
    }

    return result;
}

int Solver::cost(Indices const& indices) {
    int result = 0;

    for (int i = 1; i < indices.size(); i++)
        result += distances[indices[i - 1]][indices[i]];
    result += distances[indices.front()][indices.back()];

    for (auto const i : indices)
        result += data[i].cost;

    return result;
}

std::mt19937 Solver::mt {std::random_device{}()};

int Solver::euclideanDistance(int x1, int y1, int x2, int y2) {
    return std::round(std::sqrt(std::pow(x1 - x2, 2) + std::pow(y1 - y2, 2)));
}
