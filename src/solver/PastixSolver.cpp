//
// Created by dwalther on 23/03/23.
//

#include "PastixSolver.hpp"
#include <mpi.h>
#include <pastix.h>
#include <spm.h>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include "System.hpp"
#include <chrono>
//#include "solver/Matrix.hpp"

PastixSolver::PastixSolver( int solver_argc, char *solver_argv[], int *mpi_argc, char **mpi_argv[] )
    : filename( NULL ), x( NULL ), b( NULL ), x0( NULL ), check( 1 ), scatter( 0 ), nrhs( 1 )
{
    MPI_Init( mpi_argc, mpi_argv );
    pastixInitParam( iparm, dparm );
    pastixGetOptions( solver_argc, solver_argv, iparm, dparm, &check, &scatter, &driver,
                      &filename );
                 iparm[IPARM_THREAD_NBR] = 6; 
    pastixInit( &pastix_data, (PASTIX_Comm)MPI_COMM_WORLD, iparm, dparm );

}

PastixSolver::~PastixSolver()
{
    spmExit( spm );
    free( spm );
    free( x );
    free( b );

    pastixFinalize( &pastix_data );
    MPI_Finalize();

    spm = NULL;
    x   = NULL;
    b   = NULL;
}

void
PastixSolver::solve( GridDimension N, unsigned int formRadius, unsigned int nbElectrodes,
                     real_t coeff )
{
    int  rc;
    Form myForm( formRadius, nbElectrodes );
    myForm.coeff( 0 ) = coeff;
      
    for ( uint_t i = 0; i < 4; ++i )
        myForm.begAngle( i ) = (i)*2. * EIT_PI / 4 - EIT_PI;
        //myForm.begAngle( i ) = (i)*2. * EIT_PI / 4 - EIT_PI;

    // Set the length of the elctrodes to 0.4 :: this will set the end angles via computation
    myForm.buildThetas(1.5);
    // Print the beginning angle and the end angle of each electrode for ( uint_t i = 0; i < 4; ++i
    for ( uint_t i = 0; i < 4; ++i )
    {
        std::cout << myForm.getThetas()[i][0] << "  " << myForm.getThetas()[i][1] << std::endl;
    }

    // Create the object grid of type Grid related to the object myForm of type form.
    Grid grid( N, &myForm );
    // Build the Grid Points
    grid.buildPoints();

    grid.buildInterfaces();
    // Modify neighbors
    grid.correctConnections();

    EitBiVector<real_t> sigma = grid.createBiVector<real_t>();
    EitBiVector<real_t> u     = grid.createBiVector<real_t>();
    EitBiVector<real_t> v     = grid.createBiVector<real_t>();



    for ( uint_t i = 0U; i < grid.getNumberOfPoints(); ++i )
    {
        sigma.at( i ) = 1;
    }


 auto start = std::chrono::high_resolution_clock::now();
    System system( grid );
    system.buildFullSystem( grid, sigma, u );
    system.writeMatrixMarket();
    spm = system.buildSpm();
    auto stop = std::chrono::high_resolution_clock::now();

    spmPrintInfo( spm, stdout );

    rc = spmCheckAndCorrect( spm, &spm2 );
    if ( rc != 0 )
    {
        spmExit( spm );
        *spm = spm2;
        rc   = 0;
    }

    //Generate a Fake values array if needed for the numerical part
    if ( spm->flttype == SpmPattern )
    {
        spmGenFakeValues( spm );
    }
 
    // Perform ordering, symbolic factorization, and analyze steps
    pastix_task_analyze( pastix_data, spm );

    //Normalize A matrix (optional, but recommended for low-rank functionality)
    double normA = spmNorm( SpmFrobeniusNorm, spm );
    spmScal( 1. / normA, spm );

   // Perform the numerical factorization
    pastix_task_numfact( pastix_data, spm );

    // Generates the b and x vector such that A * x = b
    size = pastix_size_of( spm->flttype ) * spm->nexp * nrhs;
    x    = malloc( size );
    b    = malloc( size );

    if ( check > 1 )
    {
        x0 = malloc( size );
    }


    pastix_subtask_spm2bcsc( pastix_data, spm );
    pastix_subtask_bcsc2ctab( pastix_data );
    pastix_subtask_sopalin( pastix_data );

    if ( check )
    {
        spmGenRHS( SpmRhsRndX, nrhs, spm, x0, spm->nexp, b, spm->nexp );
        memcpy( b, x, size );
         
    } else {
        spmGenRHS( SpmRhsRndB, nrhs, spm, NULL, spm->nexp, x, spm->nexp );
 
         /* Apply also normalization to b vectors */
        spmScalMat( 1./normA, spm, nrhs, b, spm->nexp );
 
        /* Save b for refinement */
        memcpy( b, x, size );
    }

        /* Solve the linear system
             */
        pastix_task_solve( pastix_data, spm->nexp, nrhs, x, spm->nexp );
        pastix_task_refine( pastix_data, spm->nexp, nrhs, b, spm->nexp, x, spm->nexp );
 
        if ( check ) {
            rc |= spmCheckAxb( dparm[DPARM_EPSILON_REFINEMENT], nrhs, spm, x0, spm->nexp, b, spm->nexp, x, spm->nexp );
        }
 


    //Solve the linear system (and perform the optional refinement)
    // auto start = std::chrono::high_resolution_clock::now();
    // pastix_task_solve( pastix_data, spm->nexp, nrhs, x, spm->nexp );
    // pastix_task_refine( pastix_data, spm->nexp, nrhs, b, spm->nexp, x, spm->nexp );
    //auto stop = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    std::cout << "Time taken: " << duration.count() << " milliseconds" <<"  "<<grid.getPads()<< std::endl;


    // if ( check )
    // {
    //     rc = spmCheckAxb( dparm[DPARM_EPSILON_REFINEMENT], nrhs, spm, x0, spm->nexp, b,
    //     spm->nexp,
    //                       x, spm->nexp );

    //     if ( x0 )
    //     {
    //         free( x0 );
    //     }
    // }

    // if ( rc != SPM_SUCCESS )
    // {
    //     throw std::runtime_error( "spmCheckAxb failed" );
    // }

//     // double   *data = (double *)b;
//     // spm_int_t i;

//     //     double errmax = 0;
//     //     for ( i = 0; i <grid.getNumberOfPoints(); i++ )
//     //     {
//     //         if(myForm.levelSet(*grid.getPoint(i))<0){

//     //            double x =  (grid.getPoint(i))->getCoordinate(0);
//     //            double y =  (grid.getPoint(i))->getCoordinate(1);
//     //            double exactsol = exp(y*y + x*x);
//     //            std::cout<<abs(data[i])<<"   "<<exactsol<<std::endl;

//     //         //    if (errmax < abs((data[i] - exactsol)))
//     //         //    {
//     //         //      errmax = abs((data[i] - exactsol));
//     //         //    }
//     //         }
              
//     //     }

//     // std::cout << '\n';
 
//     // std::cout << "errmax:"
//     //     << " " << std::setprecision(17) << errmax << std::endl;

}

void
PastixSolver::exportMetrics( char *filename ) const
{
    fprintf( stderr, "%s\n", filename ); // shut up compiler
}

void
PastixSolver::exportResults( char *filename ) const
{
    // Get the file extension
    char *ext = strrchr( filename, '.' );
    if ( !ext )
    {
        printf( "Invalid filename\n" );
        return;
    }

    double   *data = (double *)x;
    spm_int_t i;

    ofstream file;
    file.open( filename );

    // Write the data based on the file extension
    if ( strcmp( ext, ".txt" ) == 0 )
    {
        for ( i = 0; i < spm->n; i++ )
        {
            file << data[i] << std::endl;
        }
    }
    else if ( strcmp( ext, ".csv" ) == 0 )
    {
        file << "Index,Value" << std::endl;
        for ( i = 0; i < spm->n; i++ )
        {
            file  << data[i] << std::endl;
        }
    }
    else
    {
        printf( "Unsupported file format\n" );
        file.close();
        return;
    }

    file.close();
}