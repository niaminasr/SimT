/**
 * @file Solver.hpp
 * @brief Contains the implementation of iterative linear solvers for the reconstruction process.
 */
 
#ifndef SRC_SOLVER_SOLVER
#define SRC_SOLVER_SOLVER

#include "ParEITConfig.hpp"
#include "solver/Matrix.hpp"
#include "tools/EitBiVector.hpp"

/**
 * @brief Solves a linear system using the BiConjugate Gradient Method.
 * @param A The matrix representing the linear system.
 * @param b The right-hand side vector of the linear system.
 * @return The norm of the error between the solution and the right-hand side.
 */
real_t
BiConjuguateGradient( Matrix A, const EitBiVector<real_t> &b )
{
    // uint_t pSize = b.posSize();
    // uint_t nSize = b.negSize();

    size_t Tot  = 2 * 500 * 500;
    double eps0 = 1.0E-8, eps = 1.0E-8;

    EitBiVector<real_t> x( b ), r0( b ), r( b ), p( b ), s( b ), r_bis( b ), Ap( b ), As( b );
    x.ones();

    real_t alpha, omega, beta;

    r  = b - A * x;
    r0 = r;
    p  = r0;

    for ( uint_t i = 0; i < Tot; ++i )
    {
        Ap = A * p;

        alpha = specialDot( r, r0 ) / specialDot( Ap, r0 );
        s     = r - alpha * Ap;

        if ( norm( s ) < eps0 )
        {
            x += alpha * p;
            break;
        }

        As    = A * s;
        omega = specialDot( As, s ) / specialDot( As, As );
        x += alpha * p + omega * s;

        r_bis = s - omega * As;

        if ( norm( r_bis ) < eps )
        {
            r = r_bis;
            break;
        }

        beta = ( alpha / omega ) * specialDot( r_bis, r0 ) / specialDot( r, r0 );
        p    = r_bis + beta * ( p - omega * Ap );

        if ( specialDot( r_bis, r0 ) < 1.0E-6 )
        {
            r0 = r_bis;
            p  = r_bis;
        }

        r = r_bis;
    }
    std::cout << "norm error" << norm( A * x - b ) << std::endl;

    return 1;
}

/**
 * @brief Solves a linear system using the Conjugate Gradient Method.
 * @param A The matrix representing the linear system.
 * @param b The right-hand side vector of the linear system.
 * @return The norm of the error between the solution and the right-hand side.
 */
real_t
ConjuguateGradient( Matrix A, const EitBiVector<real_t> &b )
{
    size_t Tot = 2 * 300;
    double eps = 1e-8;

    EitBiVector<real_t> x( b ), r( b ), p( b ), Ap( b ), As( b ), r_bis( b );
    x.ones();

    r = b - A * x;
    p = r;

    size_t k;
    for ( k = 0; k < Tot; k++ )
    {
        Ap = A * p;

        double alpha = specialDot( r, r ) / specialDot( Ap, p );

        x += alpha * p;

        r_bis = r;
        r_bis -= alpha * Ap;

        double value = std::sqrt( specialDot( r_bis, r_bis ) );
        std::cout << value << "  " << k << std::endl;

        if ( value < eps )
        {
            r = r_bis;
            break;
        }

        double beta = specialDot( r_bis, r_bis ) / specialDot( r, r );

        // p *= beta;
        // p += r_bis;
        r = r_bis;
        p = beta * p + r_bis;
    }

    return 1;
}
#endif /* SRC_SOLVER_SOLVER */
