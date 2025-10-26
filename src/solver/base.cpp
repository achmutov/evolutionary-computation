#include <evolutionary_computation/solver/base.h>
#include <random>
#include <ranges>
#include <tuple>
#include <vector>
#include <chrono>

void Solver::init(Data const& data) {
    this->data = data;
    this->distances = toMatrix(data);
}

std::tuple<Solver::Indices, int, Solver::Duration> Solver::solve(int i) {
    auto clock = std::chrono::high_resolution_clock();
    auto start = clock.now();
    auto indices = this->_solve(i);
    auto duration  = clock.now() - start;
    auto cost = this->cost(indices);
    return std::make_tuple(indices, cost, duration);
}

std::vector<std::vector<int>> Solver::toMatrix(Data const& data) {
    auto result = std::vector<std::vector<int>>(data.entries.size());

    for (size_t i = 0; i < data.entries.size(); ++i) {
        auto const& a = data.entries[i];
        result[i] = std::vector<int>(data.entries.size());
        for (size_t j = 0; j < data.entries.size(); ++j) {
            auto const& b = data.entries[j];
            result[i][j] = euclideanDistance(a.x, a.y, b.x, b.y);
        }
    }

    return result;
}

int Solver::cost(Indices const& indices) const {
    int result = 0;

    for (int i = 1; i < indices.size(); i++)
        result += distances[indices[i - 1]][indices[i]];
    result += distances[indices.front()][indices.back()];

    for (auto const i : indices)
        result += data.entries[i].cost;

    return result;
}

std::mt19937 Solver::mt { 17 };

int Solver::euclideanDistance(int x1, int y1, int x2, int y2) {
    return std::round(std::sqrt(std::pow(x1 - x2, 2) + std::pow(y1 - y2, 2)));
}
