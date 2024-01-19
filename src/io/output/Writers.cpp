#include "io/output/Writers.hpp"

#include <fstream>
#include <iomanip>
EIT_DECL_PARALLEL_DEFINITIONS

void
writeGnuplot( std::string filename, const Grid &grid, const EitBiVector<real_t> &vec, bool onlyItf )
{
    std::ofstream file( filename + std::to_string( parallel::rank ) + ".dat" );

    file << "#" << std::setw( 30 ) << "X"
#if EIT_DIM >= 2
         << std::setw( 30 ) << "Y"
#if EIT_DIM >= 3
         << std::setw( 30 ) << "Z"
#endif
#endif
         << std::setw( 30 ) << "V" << std::endl;

    if ( onlyItf )
    {
        const std::vector<Point *> &itfs   = grid.getInterfaces();
        const int_t                &n_itfs = itfs.size();

        for ( int_t i = 0; i < n_itfs; ++i )
            file << std::setw( 30 ) << itfs[i]->getCoordinate( 0 )
#if EIT_DIM >= 2
                 << std::setw( 30 ) << itfs[i]->getCoordinate( 1 )
#if EIT_DIM >= 3
                 << std::setw( 30 ) << itfs[i]->getCoordinate( 1 )
#endif
#endif
                 << std::setw( 30 ) << vec.at( -i ) << std::endl;
    }
    else
    {
        const std::vector<Point *> &pts   = grid.getPoints();
        const uint_t               &n_pts = pts.size();

        for ( uint_t i = 0U; i < n_pts; ++i )
            file << std::setw( 30 ) << pts[i]->getCoordinate( 0 )
#if EIT_DIM >= 2
                 << std::setw( 30 ) << pts[i]->getCoordinate( 1 )
#if EIT_DIM >= 3
                 << std::setw( 30 ) << pts[i]->getCoordinate( 1 )
#endif
#endif
                 << std::setw( 30 ) << vec.at( i ) << std::endl;
    }

    file.close();
    return;
}
