#include "tests.hpp"

// EIT_DECL_PARALLEL_DEFINITIONS;


int test_assemble(){


    float epsilon = 0.001;
    MPI_Init( 0, NULL);
    MPI_Comm_rank( MPI_COMM_WORLD, &parallel::rank );
    MPI_Comm_size( MPI_COMM_WORLD, &parallel::size );
    GridDimension N = { 3, 3 };

    // Create the form :: a disk of radius 0.5 with 4 electrodes
    Form myform( 1, 4 );
    myform.coeff( 0 ) = 1;
    // Setting the beginning angle of an electrode
    for ( uint_t i = 0; i < 4; ++i )
    {
        myform.begAngle( i ) = (i)*2. * EIT_PI / 4;
    }
    // Set the length of the elctrodes to 0.4 :: this will set the end angles via computation
    myform.buildThetas( 0.4 );
    Grid grid( N, &myform );
    // Build the Grid Points
    grid.buildPoints();

    grid.buildInterfaces();
    // Build the halos for communications
    grid.buildHalos();
    // Modify neighbors
    grid.correctConnections();
    MPI_Barrier( MPI_COMM_WORLD );
    EitBiVector<real_t> sigma = grid.createBiVector<real_t>();

    EitBiVector<real_t> u = grid.createBiVector<real_t>();

    for ( uint_t i = 0U; i < grid.getNumberOfPoints(); ++i )
    {
        u.at( i )     = 1;
        //int pair      = rand() % 1;
        sigma.at( i ) = 1;
    }

    for ( int_t i = 1; i <=  EIT_ELEC; ++i )
    {
        u.at( -i ) = 1;
    }

    System system( grid );

    system.buildFullSystem(grid, sigma, u);

// Itfs points
    TEST(system.x_coords->at(2) == 0);
    TEST(system.y_coords->at(2) == 0);
    TEST((system.coefs->at(2) - 1 ) < epsilon);
    
    TEST(system.x_coords->at(3) == 0);
    TEST(system.y_coords->at(3) == 6);
    TEST((system.coefs->at(3) - (-1) ) < epsilon);
    
// Electrode points
    TEST(system.x_coords->at(7) == 1);
    TEST(system.y_coords->at(7) == 4);
    TEST((system.coefs->at(7) - (-1) ) < epsilon);
    
    TEST(system.x_coords->at(4) == 1);
    TEST(system.y_coords->at(4) == 12);
    TEST((system.coefs->at(4) - (-1) ) < epsilon);

//Regular point
    TEST(system.x_coords->at(16) == 4);
    TEST(system.y_coords->at(16) == 12);
    TEST((system.coefs->at(16) - 0 ) < epsilon);
    
    TEST(system.x_coords->at(15) == 4);
    TEST(system.y_coords->at(15) == 4);
    TEST((system.coefs->at(15) - 1 ) < epsilon);

    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Finalize();

    RETURN_TEST;
}


int ( *test_list[] )() = { test_assemble };

int
main( int argc, char **argv )
{
    int test_number = 1;
    if ( argc > 1 )
        test_number = atoi( argv[1] );

    return test_list[test_number - 1]();
}