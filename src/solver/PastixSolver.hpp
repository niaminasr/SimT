//
// Created by dwalther on 23/03/23.
//

#ifndef SIMT_PASTIXSOLVER_HPP
#define SIMT_PASTIXSOLVER_HPP

#include "Solver.hpp"
#include <pastix.h>
#include <spm.h>
#include "geometry/Point.hpp"

class PastixSolver : public Solver
{
   private:
    pastix_data_t *pastix_data;       /*< Pointer to the storage structure required by pastix */
    pastix_int_t   iparm[IPARM_SIZE]; /*< Integer in/out parameters for pastix                */
    double         dparm[DPARM_SIZE]; /*< Floating in/out parameters for pastix               */
    spm_driver_t   driver;
    char          *filename;
    spmatrix_t    *spm, spm2;
    void          *x, *b, *x0;
    size_t         size;
    int            check;
    int            scatter;
    int            nrhs;

   public:
    PastixSolver( int solver_argc, char *solver_argv[], int *mpi_argc, char **mpi_argv[] );
    virtual ~PastixSolver();
    virtual void solve( GridDimension N, unsigned int formRadius, unsigned int nbElectrodes, real_t coeff ) override;
    virtual void exportMetrics( char *filename ) const override;
    virtual void exportResults( char *filename ) const override;
};

#endif // SIMT_PASTIXSOLVER_HPP
