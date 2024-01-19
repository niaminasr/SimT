#include <mpi.h>
#include "geometry/Grid.hpp"
#include "io/input/parameters.hpp"
#include "io/output/Writers.hpp"
#include "solver/PastixSolver.hpp"
#include "solver/SolverFabric.hpp"
#include "solver/System.hpp"
#include "solver/Solver.hpp"
#include <chrono>

// EIT_DECL_PARALLEL_DEFINITIONS;
/*
!!!!!! We are in the [-2 2] square
*/
void eitGetOptions(int argc, char* argv[], float* , int*, unsigned int*,  int*);

int
main( int argc, char **argv )
{

    unsigned int gridDimension;
    float formRadius;
    int nbElectrodes, solverType;
    eitGetOptions(argc, argv, &formRadius, &nbElectrodes, &gridDimension, &solverType);

    GridDimension N = { gridDimension, gridDimension };
    SolverFabric solverFabric;
    switch( solverType)
    {
        case (int)SolverType::PASTIX_SOLVER:
            solverFabric.setSolverType( SolverType::PASTIX_SOLVER );
        default:
            std::cout << "bad solver type" << std::endl;
    }
    Solver        *mySolver = solverFabric.makeSolver();

      char* salut ;
      salut = "SOL.txt" ;
     mySolver->solve(N, 1, 4, 1);
     mySolver->exportResults(salut);
    
    return 0;
}
