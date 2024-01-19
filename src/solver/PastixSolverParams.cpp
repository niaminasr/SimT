//
// Created by dwalther on 07/04/23.
//

#include <cstdlib>
#include <cstring>
#include "PastixSolverParams.hpp"

const char *default_argv[18] = { "eit", "--mm",   "matrix.mtx", "--threads", "-1", "--gpus",  "0",
                                 "--sched",      "1",      "fact",      "2",  "--check", "1",
                                 "--ord",        "scotch", "--verbose", "1",  NULL };
PastixSolverParams::PastixSolverParams()
: _solver_argc(0), _solver_argv(NULL), _mpi_argc(NULL), _mpi_argv(NULL)
{
}

PastixSolverParams::~PastixSolverParams()
{
    if (_solver_argv != NULL){
        delete [] _solver_argv;
    }
}

void
PastixSolverParams::setDefault()
{
    if ( _solver_argv != NULL )
    {
        free( _solver_argv );
        _solver_argv = NULL;
    }
    _solver_argc           = 17;
    _solver_argv           = new char *[_solver_argc + 1];

    for (int i = 0; i < _solver_argc; i++){
        _solver_argv[i] = new char[64]; // 64 is an arbitrary value (it is enough for most cases)
        std::strcpy(_solver_argv[i], default_argv[i]);
    }
    _solver_argv[_solver_argc] = NULL;

    _mpi_argc = NULL;
    _mpi_argv = NULL;
}

int
PastixSolverParams::getSolverArgc() const
{
    return _solver_argc;
}

char **
PastixSolverParams::getSolverArgv() const
{
    return _solver_argv;
}

int*
PastixSolverParams::getMPIArgc() const
{
    return _mpi_argc;
}

char ***
PastixSolverParams::getMPIArgv() const
{
    return _mpi_argv;
}

void
PastixSolverParams::setIthDParm( int IPARM_IDX, double value )
{
    char **newSolverArgv = new char *[_solver_argc + 4];
    memcpy(newSolverArgv, _solver_argv, _solver_argc);
    // TODO : add "--i", IPARM_IDX, value, NULL to newSolverArgv


    if (_solver_argv != NULL){
        delete [] _solver_argv;
    }
}