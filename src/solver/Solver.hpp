/**
 * @brief Header file for the Solver interface.
 *
 * @author Dimitri Walther
 * @date 02/03/2023
 */

#ifndef SIMT_SOLVER_HPP
#define SIMT_SOLVER_HPP

#include "geometry/Point.hpp"

class Solver
{
   public:
    virtual ~Solver() {}
    virtual void solve( GridDimension N, unsigned int formRadius, unsigned int nbElectrodes, real_t coeff ) = 0;
    virtual void exportMetrics( char *fileName ) const                                        = 0;
    virtual void exportResults( char *fileName ) const                                        = 0;
};

#endif // SIMT_SOLVER_HPP
