/*!
 * @brief Header file for the SolverFabric class. This class allows the user
 * to create an instance of a solver class that implements the Solver interface.
 * The solver will only be instantiated when the makeSolver method is called.
 * This allows for efficient parameters personalisation without creating a new
 * instance every time.
 *
 * @author Dimitri Walther
 * @date 02/03/2023
 */

#ifndef SIMT_SOLVERFABRIC_HPP
#define SIMT_SOLVERFABRIC_HPP

#include "PastixSolver.hpp"
#include "PastixSolverParams.hpp"
#include "Solver.hpp"

/// Allows the user to easily choose what type of solver the SolverFabric shall return
enum class SolverType
{
    BICONJUGATE_GRADIENT_SOLVER,
    PASTIX_SOLVER,
    UNDEFINED_SOLVER = -1
};

/**
 *
 */
class SolverFabric
{
   private:
    struct PastixSolverParams *_pastixParams;
    SolverType                 _solverType;

   public:
    SolverFabric();
    SolverFabric( SolverType t );
    // TODO : add SolverFabric( SolverType t, int solverArgc, char ** solverArgv, int mpiArgc, char
    // ** mpiArgv)
    ~SolverFabric();
    void    setSolverType( SolverType t );
    void    setDefaultParam( void );
    Solver *makeSolver( void );
};

#endif // SIMT_SOLVERFABRIC_HPP
