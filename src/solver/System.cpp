#include "System.hpp"

// Not 3D yet
void
System::buildFullSystem( const Grid &grid, const EitBiVector<real_t> &sigma,
                         const EitBiVector<real_t> &u )
{

    this->assembleRegularPoint(grid, sigma, u, 0 );
    // for ( int i = 0; i < grid.getNumberOfInterfaces(); i++ )
    // {
    //     this->assembleInterfacePoint( grid, sigma, u, i );
    // }

    //  for ( int i = 0; i < EIT_ELEC; i++ )
    //  {
    //    this->assembleElectrodePoint( grid, sigma, u, i );
    //  }

    // for ( int i = 0; i < grid.getNumberOfPoints(); i++ )
    // {
    //     ;
    //    uint_t rc = this->assembleRegularPoint( grid, sigma, u, i );
    //     if ( rc != 0 )
    //         return;
    // }

    /* update nnz number */
    nnz  = coefs->size();
    ncol = grid.getNumberOfInterfaces() + grid.getNumberOfPoints() + EIT_ELEC;
    nrow = grid.getNumberOfInterfaces() + grid.getNumberOfPoints() + EIT_ELEC;
}