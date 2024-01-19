#include "solver/System.hpp"

/**
 * @brief If no comprendo, refer to Matrix.hpp
 *
 * @param grid The grid that the matrix is built from
 * @param sigma the sigma vector that defines the problem
 * @param u the right hand term of the problem
 * @param i current interface point
 */
void
System::assembleElectrodePoint( const Grid &grid, const EitBiVector<real_t> &sigma,
                                const EitBiVector<real_t> &u, uint_t i )
{
   real_t hx = grid.getPads()[0];
    real_t hy = grid.getPads()[1];

    int sizeStencil = 4;

    const std::vector<Point *> &pts  = grid.getPoints();
    const std::vector<Point *> &itfs = grid.getInterfaces();

    // const uint_t &n_pts  = pts.size();
    const uint_t &n_itfs = itfs.size();
    // (void) n_pts;

    std::vector<uint_t> idx;
    std::vector<real_t> kij;
    real_t              Um = 0.0;

    PointNeighbors neighs;

    for ( unsigned int k = 0; k < n_itfs; k++ )
    {
        if ( itfs[k]->getisOnElectrode() == (int) i )
        {
            neighs = itfs[k]->getNeighs();
            for ( int l = 0; l < sizeStencil; l++ )
            {
                if ( neighs[l] != EIT_INVALID_INT )
                {
                    Um += hx * hy * grid.getForm()->approximateIntegral( *( pts[neighs[l]] ), hy );

                    kij.push_back( -hx * hy *
                                   grid.getForm()->approximateIntegral( *( pts[neighs[l]] ), hy ) );
                    idx.push_back( neighs[l] );
                }
            }
        }
    }

    x_coords->push_back( n_itfs + i );
    y_coords->push_back( n_itfs + i );
    coefs->push_back( i == 0 ? Um + 0.0000000001 : Um );

    for ( uint_t k = 0; k < idx.size(); k++ )
    {
        x_coords->push_back( n_itfs + i );
        y_coords->push_back( idx.at( k ) + n_itfs + EIT_ELEC );
        coefs->push_back( kij.at( k ) );
    }

    return;
}