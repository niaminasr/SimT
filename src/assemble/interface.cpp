#include "solver/System.hpp"

/**
 * @brief If no comprendo, refer to Matrix.hpp
 *
 * @param grid The grid that the matrix is built from
 * @param sigma the sigma vector that defines the problem
 * @param u the right hand term of the problem
 * @param i current interface point
 *
 */
void
System::assembleInterfacePoint( const Grid &grid, const EitBiVector<real_t> &sigma,
                                const EitBiVector<real_t> &u, uint_t i )
{
    const std::vector<Point *> &pts  = grid.getPoints();
    const std::vector<Point *> &itfs = grid.getInterfaces();

    const uint_t &n_pts  = pts.size();
    const uint_t &n_itfs = itfs.size();
    (void) n_pts;

    real_t alphaj, betaj, alphai, betai, alphak, betak;
    real_t denom;

    // The Interface points lines of the matricial vector product. <- weird comment
    uint_t sizeStencil_itfs = 3;
    int    idx_itfs[sizeStencil_itfs];
    real_t s_itfs[sizeStencil_itfs], uij_itfs[sizeStencil_itfs], kij_itfs[sizeStencil_itfs];
    real_t x[sizeStencil_itfs], y[sizeStencil_itfs];

    idx_itfs[2] = ( -(int_t)i - 1 );
    s_itfs[2]   = sigma.at( idx_itfs[2] );
    uij_itfs[2] = u.at( idx_itfs[2] );
    x[2]        = itfs[i]->getCoordinate( 0 );
    y[2]        = itfs[i]->getCoordinate( 1 );

    /* Get the stencil needed to descretizise the normal
    derivative on the interface pointes (the negative part of Bivector)
    */
    idx_itfs[0] = itfs[i]->getNormalStencil()[0];
    idx_itfs[1] = itfs[i]->getNormalStencil()[1];

    uij_itfs[0] = u.at( idx_itfs[0] );
    uij_itfs[1] = u.at( idx_itfs[1] );
    (void) uij_itfs;
    s_itfs[0]   = sigma.at( idx_itfs[0] );
    s_itfs[1]   = sigma.at( idx_itfs[0] );

    // set coordinates
    x[0] = pts[idx_itfs[0]]->getCoordinate( 0 );
    y[0] = pts[idx_itfs[0]]->getCoordinate( 1 );

    x[1] = pts[idx_itfs[1]]->getCoordinate( 0 );
    y[1] = pts[idx_itfs[1]]->getCoordinate( 1 );

    denom  = ( x[2] - x[0] ) * ( y[2] - y[1] ) - ( x[2] - x[1] ) * ( y[2] - y[0] );
    alphaj = ( y[0] - y[1] ) / denom;
    betaj  = ( x[1] - x[0] ) / denom;

    denom  = ( x[1] - x[0] ) * ( y[1] - y[2] ) - ( x[1] - x[2] ) * ( y[1] - y[0] );
    alphai = ( y[0] - y[2] ) / denom;
    betai  = ( x[2] - x[0] ) / denom;

    denom  = ( x[0] - x[2] ) * ( y[0] - y[1] ) - ( x[0] - x[1] ) * ( y[0] - y[2] );
    alphak = ( y[2] - y[1] ) / denom;
    betak  = ( x[1] - x[2] ) / denom;

    kij_itfs[0] =
        s_itfs[0] * ( alphak * itfs[i]->getNormal()[0] + betak * itfs[i]->getNormal()[1] );
    kij_itfs[1] =
        s_itfs[1] * ( alphai * itfs[i]->getNormal()[0] + betai * itfs[i]->getNormal()[1] );
    kij_itfs[2] =
        s_itfs[2] * ( alphaj * itfs[i]->getNormal()[0] + betaj * itfs[i]->getNormal()[1] );

    for ( uint_t k = 0; k < sizeStencil_itfs - 1; k++ )
    {
        x_coords->push_back( i );
        y_coords->push_back( idx_itfs[k] + n_itfs + EIT_ELEC );
        coefs->push_back( kij_itfs[k] );
    }

    if ( itfs[-idx_itfs[2] - 1]->getisOnElectrode() == -1 )
    {
        x_coords->push_back( i );
        y_coords->push_back( i );
        coefs->push_back( kij_itfs[2] );
    }
    else
    {
        x_coords->push_back( i );
        y_coords->push_back( i );
        coefs->push_back( kij_itfs[2] );

        /* adding the term for the electrode line */
        x_coords->push_back( i );
        y_coords->push_back( n_itfs + itfs[-idx_itfs[2] - 1]->getisOnElectrode() );
        coefs->push_back( -1 );
    }

    b->push_back( u.at( -(uint_t)i - 1 ) );
}