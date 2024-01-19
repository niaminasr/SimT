/**
 * @brief Class file for the SolverFabric class. Implements methods defined in the SolverFabric.hpp
 * header file.
 *
 * @author Dimitri Walther
 * @date 02/03/2023
 */

#include "SolverFabric.hpp"
#include <pastix.h>
#include <stdlib.h>
#include <stdexcept>
#include "PastixSolver.hpp"
#include "SolverParams.hpp"

SolverFabric::SolverFabric() : _pastixParams( NULL ), _solverType( SolverType::UNDEFINED_SOLVER ) {}

SolverFabric::SolverFabric( SolverType t ) : _pastixParams( NULL ), _solverType( t )
{
    setSolverType( _solverType );
}

SolverFabric::~SolverFabric() {}

void
SolverFabric::setSolverType( SolverType t )
{
    _solverType = t;
    switch ( _solverType )
    {
    case SolverType::PASTIX_SOLVER:
        if ( _pastixParams == NULL )
        {
            _pastixParams = new PastixSolverParams();
            _pastixParams->setDefault();
        }
        break;
    case SolverType::BICONJUGATE_GRADIENT_SOLVER:
    default:
        throw std::runtime_error( "Error : Undefined solver type" );
    }
}

void
SolverFabric::setDefaultParam()
{
    switch ( _solverType )
    {
    case SolverType::PASTIX_SOLVER:
        _pastixParams->setDefault();
    default:
        throw std::runtime_error( "Error : Undefined solver type" );
    }
}

Solver *
SolverFabric::makeSolver()
{
    switch ( _solverType )
    {
    case SolverType::PASTIX_SOLVER: {
        int     solverArgc = _pastixParams->getSolverArgc();
        char  **solverArgv = _pastixParams->getSolverArgv();
        int    *mpiArgc    = _pastixParams->getMPIArgc();
        char ***mpiArgv    = _pastixParams->getMPIArgv();
        return new PastixSolver( solverArgc, solverArgv, mpiArgc, mpiArgv );
    }
    default:
        throw std::runtime_error( "Error : Undefined solver type" );
    }
}