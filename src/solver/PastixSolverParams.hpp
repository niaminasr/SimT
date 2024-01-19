//
// Created by dwalther on 07/04/23.
//

#ifndef SIMT_PASTIXSOLVERPARAMS_HPP
#define SIMT_PASTIXSOLVERPARAMS_HPP

#include "SolverParams.hpp"

class PastixSolverParams : public SolverParams
{
   private:
    int     _solver_argc;
    char  **_solver_argv;
    int    *_mpi_argc;
    char ***_mpi_argv;

   public:
    /**
     * @brief Default constructor for PastixSolverParams /!\ HAS TO BE INITIALIZED WITH setDefault()
     * BEFORE MANUALLY SETTING PARAMETERS
     */
    PastixSolverParams();
    PastixSolverParams( int solver_argc, char *solver_argv[] );
    ~PastixSolverParams();
    void    setDefault() override;
    int     getSolverArgc() const;
    char  **getSolverArgv() const;
    int    *getMPIArgc() const;
    char ***getMPIArgv() const;
    void    setIthIParm( int IPARM_IDX, int value );
    void    setIthDParm( int IPARM_IDX, double value );
};

#endif // SIMT_PASTIXSOLVERPARAMS_HPP
