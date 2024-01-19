#include <spm.h>
#include <fstream>
#include "solver/System.hpp"

void
System::writeMatrixMarket()
{
    std::ofstream file( "matrix.mtx" );

    // Write header
    file << "%%MatrixMarket matrix coordinate real general\n\n";
    file << nrow << " " << ncol << " " << nnz << "\n";

    // Write matrix elements
    for ( int i = 0; i < nnz; i++ )
    {
        file << x_coords->at( i ) << " " << y_coords->at( i ) << " " << coefs->at( i ) << "\n";
    }
    file << "%%WHS\n\n";
    for ( int i = 0; i < b->size(); i++ )
    {
        file << b->at( i ) << " " << i << "\n";
    }

    file.close();
}

spmatrix_t *
System::buildSpm()
{
    spmatrix_t *spm = (spmatrix_t *)malloc( sizeof( spmatrix_t ) );

    spmInit( spm );

    spm->mtxtype = SpmGeneral;
    spm->flttype = SpmDouble;
    spm->fmttype = SpmIJV;
    spm->baseval = 0;

    spm->n   = this->ncol;
    spm->nnz = this->nnz;
    spm->dof = 1; /* for now */

    /*
     * Update the computed fields.
     * If dof < 0, please see example_user2.c
     **/
    spmUpdateComputedFields( spm );

    /*  spm->colptr = (spm_int_t*)this->x_coords->data();
        spm->rowptr = (spm_int_t*)this->y_coords->data();
        spm->values = this->coefs->data(); */

    spm->colptr  = (spm_int_t *)malloc( sizeof( int64_t ) * this->nnz );
    spm->rowptr  = (spm_int_t *)malloc( sizeof( int64_t ) * this->nnz );
    spm->values  = (real_t *)malloc( sizeof( real_t ) * this->nnz );
    real_t *data = (real_t *)spm->values;

    for ( int i = 0; i < this->nnz; i++ )
    {
        spm->colptr[i] = (spm_int_t)this->x_coords->at( i );
        spm->rowptr[i] = (spm_int_t)this->y_coords->at( i );
        data[i]        = (real_t)this->coefs->at( i );
    }

    return spm;
}