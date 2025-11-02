#include <evolutionary_computation/solver/random.h>
#include <evolutionary_computation/solver/local_search.h>
#include <evolutionary_computation/solver/nearest_neighbor_pos.h>
#include <evolutionary_computation/cli.h>
#include <vector>

int main(int argc, char* argv[]) {
    auto targetLSType = LocalSearchType::Steep;
    auto targetIntraNeighborhoodType = IntraNeighborhoodType::Edges;

    auto initSolver1 = RandomSolver();
    auto ls = LocalSearchSolver(
        targetLSType,
        targetIntraNeighborhoodType,
        initSolver1
    );

    auto initSolver2 = RandomSolver();
    auto ls_candidate_moves = LocalSearchSolver(
        targetLSType,
        targetIntraNeighborhoodType,
        initSolver2,
        10
    );

    auto solvers = std::vector<Solver*> { &ls_candidate_moves, &ls };

    cli(argc, argv, solvers);

    return 0;
}
