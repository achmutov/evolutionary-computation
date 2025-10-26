#pragma once

#include <chrono>
#include <cmath>
#include <random>
#include <tuple>
#include <vector>
#include <evolutionary_computation/data.h>

class Solver {
public:
    typedef std::vector<int> Indices;
    typedef std::vector<std::vector<int>> Distances;
    typedef std::chrono::high_resolution_clock::duration Duration;

    virtual void init(Data const& data);
    virtual std::tuple<Indices, int, Duration> solve(int i);
    virtual Indices _solve(int i) = 0;
    virtual std::string name() const = 0;

    virtual ~Solver() = default;
protected:
    Data data;
    std::vector<std::vector<int>> distances;

    static std::mt19937 mt;

    static std::vector<std::vector<int>> toMatrix(Data const& data);
    static int euclideanDistance(int x1, int y1, int x2, int y2);

    int cost(Indices const& indices) const;
};
