#pragma once

#include <cmath>
#include <random>
#include <tuple>
#include <vector>
#include <evolutionary_computation/data.h>

class Solver {
public:
    Solver(std::vector<Data> data);

    typedef std::vector<int> Indices;

    virtual std::tuple<Indices, int> solve();
protected:
    std::vector<Data> data;
    std::vector<std::vector<int>> distances;

    static std::mt19937 mt;

    static std::vector<std::vector<int>> toMatrix(std::vector<Data> const& data);
    static int euclideanDistance(int x1, int y1, int x2, int y2);

    int cost(Indices const& indices);

    virtual Indices _solve() = 0;
};
