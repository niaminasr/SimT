#include "solver/System.hpp"
#include <fstream>
#include <sstream>
#include <string>

/**
 * @brief If no comprendo, refer to Matrix.hpp
 *
 * @param grid The grid that the matrix is built from
 * @param sigma the sigma vector that defines the problem
 * @param u the right hand term of the problem
 * @param i current interface point
 *
 * @return exit code
 */
uint_t
System::assembleRegularPoint( const Grid &grid, const EitBiVector<real_t> &sigma,
                              const EitBiVector<real_t> &u, uint_t i )
{
    
    const std::vector<Point *> &pts  = grid.getPoints();
    const std::vector<Point *> &itfs = grid.getInterfaces();

    const uint_t &n_itfs = grid.getNumberOfInterfaces();

    // Get the step size
    real_t hx = grid.getPads()[0];
    real_t hy = grid.getPads()[1];

    // Initialize the stencil
    /*
    idx (or idx_itfs) is going to represent the local index of each point used in the stencil for
    the laplacien, s (or s_itfs) represents the value of the conductivity sigma at each one of these
    nodes, uij (or uij_itfs) the value of the potential "u" on ecah stencil node, kij(or kij_itfs)
    the laplacien facteur, and finally h  the distance between the center stencil node uij[4] and
    the other ones.
    */
    uint_t sizeStencil;
    switch ( grid.getPoints()[i]->getPlace() )
    {
    case GRID:
        sizeStencil = 5;
        break;
    case TOP:
    case BOTTOM:
    case LEFT:
    case RIGHT:
        sizeStencil = 4;
        break;
    case CORNER:
        sizeStencil = 3;
        break;
    default:
        fprintf( stderr, "error with point type\n" );
        return -1;
    }

    int    idx[sizeStencil];
    real_t s[sizeStencil], uij[sizeStencil], kij[sizeStencil], h[sizeStencil];

    /* Regular Points */
    if ( grid.getPoints()[i]->getPlace() == GRID )
    {
        idx[4]    = i;
        s[4]      = sigma.at( idx[4] );
        uij[4]    = u.at( idx[4] );
        real_t xi = pts[idx[4]]->getCoordinate( 0 );
        real_t yj = pts[idx[4]]->getCoordinate( 1 );

        for ( uint_t k = 0; k < sizeStencil - 1; ++k )
        {
            idx[k] = pts[idx[4]]->getNeighs()[k];
        }

        // Defining Pads :: distance between the central node and the other nodes (for each node).
        // In the case of regular nodes:
        for ( uint_t k = 0; k < sizeStencil - 3; ++k )
        {
            h[k]     = hx;
            h[k + 2] = hy;
        }

        // Modifying pads in the case of an interface point as a neighbor
        // In the case of Irregular nodes: checked in this order West, Est, South, North.
        if ( idx[0] < 0 )
            h[0] = xi - itfs[-idx[0] - 1]->getCoordinate( 0 );
        if ( idx[1] < 0 )
            h[1] = itfs[-idx[1] - 1]->getCoordinate( 0 ) - xi;
        if ( idx[2] < 0 )
            h[2] = yj - itfs[-idx[2] - 1]->getCoordinate( 1 );
        if ( idx[3] < 0 )
            h[3] = itfs[-idx[3] - 1]->getCoordinate( 1 ) - yj;

        // getting values for building coefficients
        for ( uint_t k = 0; k < sizeStencil - 1; ++k )
        {
            uij[k] = u.at( idx[k] );
            s[k]   = sigma.at( idx[k] );
        }

        // building coefficients
        kij[4] = 0.0;
        for ( uint_t k = 0; k < sizeStencil - 3; ++k )
        {
            kij[k]     = ( ( ( s[4] + s[k] ) / 2.0 ) / h[k] ) / hx;
            kij[k + 2] = ( ( ( s[4] + s[k + 2] ) / 2.0 ) / h[k + 2] ) / hx;
            kij[4] -= ( kij[k] + kij[k + 2] );
        }

        for ( uint_t k = 0; k < sizeStencil; k++ )
        {
            // change index from neighbour's idx on the grid to neighbour's line number on the
            // matrix;
            int neigh_idx = idx[k] < 0 ? -idx[k] - 1 : idx[k] + n_itfs + EIT_ELEC;

            x_coords->push_back( idx[4] + n_itfs + EIT_ELEC );
            y_coords->push_back( neigh_idx );
            coefs->push_back( -kij[k] );
        }

        b->push_back( u.at( idx[4] ) );
    }
    (void) uij;

    /* Corner points */
    if ( grid.getPoints()[i]->getPlace() == CORNER )
    {
        x_coords->push_back( i + n_itfs + EIT_ELEC );
        y_coords->push_back( i + n_itfs + EIT_ELEC );
        coefs->push_back( -4 * ( 1. / ( hy * hx ) ) );

        for ( int k = 0; k < 4; k++ )
        {
            if ( pts[i]->getNeighs()[k] != EIT_INVALID_INT )
            {
                x_coords->push_back( i + n_itfs + EIT_ELEC );
                y_coords->push_back( pts[i]->getNeighs()[k] + n_itfs + EIT_ELEC );
                coefs->push_back( 1. / ( hy * hx ) );
            }
        }
    }

    /* Top Points */
    if ( grid.getPoints()[i]->getPlace() == TOP || grid.getPoints()[i]->getPlace() == BOTTOM ||
         grid.getPoints()[i]->getPlace() == LEFT || grid.getPoints()[i]->getPlace() == RIGHT )
    {
        idx[3] = i;
        s[3]   = sigma.at( idx[3] );
        kij[3] = 0.0;

        /* Initializing the non-nul neighbors */
        uint_t         current = 0;
        PointNeighbors neighs  = pts[idx[3]]->getNeighs();
        for ( uint_t i = 0; ( i < sizeStencil ) && ( current < 3 ); i++ )
        {
            if ( neighs[i] != EIT_INVALID_INT )
            {
                idx[current++] = neighs[i];
            }
        }

        for ( uint_t k = 0; k < sizeStencil - 1; ++k )
        {
            s[k] = sigma.at( idx[k] );

            kij[k] = ( ( ( s[3] + s[k] ) / 2.0 ) / hx ) / hx;
            kij[3] -= ( kij[k] );

            // change index from neighbour's idx on the grid to neighbour's line number on the
            // matrix;
            int neigh_idx = idx[k] < 0 ? -idx[k] - 1 : idx[k] + n_itfs + EIT_ELEC;

            x_coords->push_back( idx[3] + n_itfs + EIT_ELEC );
            y_coords->push_back( neigh_idx );
            coefs->push_back( kij[k] );
        }

        x_coords->push_back( idx[3] + n_itfs + EIT_ELEC );
        y_coords->push_back( idx[3] + n_itfs + EIT_ELEC );
        coefs->push_back( kij[3] );
        b->push_back( u.at( idx[3] ) );
    }

    return 0;
}