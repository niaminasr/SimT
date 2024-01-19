#include "tests.hpp"

// EIT_DECL_PARALLEL_DEFINITIONS;

int
test_getRaduis()
{
    Form myform( 1, 4 );
    myform.coeff( 0 ) = 1;

    TEST( myform.radius( 0.0 ) == 1.0 );
    TEST( myform.radius( EIT_PI / 4 ) == 1.0 );
    TEST( myform.radius( EIT_PI ) == 1.0 );

    Form   myform1( 7, 4 );
    real_t coeffs[7] = { 1.51, 0.01, 0.05, 0.2, 0.035, 0.01, 0.1 };
    for ( int i = 0; i < 7; i++ )
    {
        myform1.coeff( i ) = coeffs[i];
    }

    TEST( myform1.radius( 0.0 ) > 0 );
    TEST( myform1.radius( EIT_PI / 2 ) > 0 );
    TEST( abs( 1.76 - myform1.radius( 0.0 ) ) < 0.1 );
    TEST( abs( 1.4 - myform1.radius( EIT_PI / 2 ) ) < 0.1 );

    Form   myform2( 7, 4 );
    real_t coeffs_[7] = { 1.6, 0.002, 0.01, 0.003, 0.035, 0.2, 0.15 };
    for ( int i = 0; i < 7; i++ )
    {
        myform2.coeff( i ) = coeffs_[i];
    }

    TEST( myform2.radius( EIT_PI ) > 0 );
    TEST( abs( -1.4 - myform1.radius( EIT_PI ) < 0.1 ) );

    RETURN_TEST;
}

int
test_getFirstDerivativeRaduis()
{
    Form myform( 1, 4 );
    myform.coeff( 0 ) = 1;

    TEST( myform.radius( 0.0 ) == 1.0 );
    TEST( myform.radius( EIT_PI / 4 ) == 1.0 );
    TEST( myform.radius( EIT_PI ) == 1.0 );

    Form   myform1( 7, 4 );
    real_t coeffs[7] = { 1.51, 0.01, 0.05, 0.2, 0.035, 0.01, 0.1 };
    for ( int i = 0; i < 7; i++ )
    {
        myform1.coeff( i ) = coeffs[i];
    }

    TEST( myform1.firstDerivativeofRadius( 0.0 ) > 0 );
    TEST( abs( 1.865 - myform1.firstDerivativeofRadius( 0.0 ) ) < 0.1 );

    TEST( myform1.firstDerivativeofRadius( EIT_PI / 2 ) > 0 );
    //TEST( abs( 1.52 - sqrt( 2 ) * 0.15 - myform1.firstDerivativeofRadius( EIT_PI / 2 ) ) < 0.1 );

    RETURN_TEST;
}

int
test_getSecondDerivativeRaduis()
{
    Form myform( 1, 4 );
    myform.coeff( 0 ) = 1;

    TEST( myform.radius( 0.0 ) == 1.0 );
    TEST( myform.radius( EIT_PI / 4 ) == 1.0 );
    TEST( myform.radius( EIT_PI ) == 1.0 );

    Form   myform1( 7, 4 );
    real_t coeffs[7] = { 1.51, 0.01, 0.05, 0.2, 0.035, 0.01, 0.1 };
    for ( int i = 0; i < 7; i++ )
    {
        myform1.coeff( i ) = coeffs[i];
    }

    TEST( myform1.secondDerivativeofRadius( 0.0 ) < 0 );
    TEST( abs( -0.5 - myform1.secondDerivativeofRadius( 0.0 ) ) < 0.01 );
    TEST( myform1.secondDerivativeofRadius( EIT_PI / 2 ) > 0 );
    //TEST( abs( 1.675 + sqrt( 2 ) * 0.45 - myform1.secondDerivativeofRadius( EIT_PI / 2 ) ) < 0.01 );

    RETURN_TEST;
}

int
test_getEndAngleForLength()
{
    real_t epsilon = 0.01;
    Form   myform( 1, 1 );
    myform.coeff( 0 )    = 1.5;
    myform.begAngle( 0 ) = 0;
    TEST( abs( myform.getEndAngleForLength( 2 * EIT_PI * 1.5 - epsilon, 0 ) ) -
          2 * EIT_PI < epsilon); // One electrode of length 2*pi*r, his beg angle is endAngle

    Form myform1( 1, 4 );
    myform1.begAngle( 0 ) = EIT_PI / 6;
    myform1.coeff(0) = 1;

    TEST( myform1.radius( 0.0 ) * abs( myform1.getEndAngleForLength( 0.5, 0 ) - EIT_PI / 6 ) - 0.5 <
          epsilon ); // L = r * (O2 - O1)

    Form   myform2( 7, 16 );
    real_t coeffs[7] = { 1.51, 0.01, 0.05, 0.2, 0.035, 0.01, 0.1 };
    for ( int i = 0; i < 7; i++ )
    {
        myform2.coeff( i ) = coeffs[i];
    }
    for ( int i = 0; i < 16; i++ )
    {
        myform2.begAngle( i ) = -EIT_PI + ( i - 1 ) * EIT_PI / 8;
    }

    TEST( myform2.radius( 0.0 ) * abs( myform2.getEndAngleForLength( 0.35, 9 ) - 0 ) - 0.35 <
          epsilon ); // L = r * (O2 - O1)

    RETURN_TEST;
}

int
test_buildThetas()
{
    real_t epsilon = 0.1;

    Form myform( 1, 4 );
    myform.coeff( 0 ) = 1.5;

    for ( int i = 0; i < 4; i++ )
    {
        myform.begAngle( i ) = EIT_PI*i/2;
    }

    myform.buildThetas( 0.3 );

    std::vector<ElectrodeThetas> thetas = myform.getThetas();
    TEST( abs( 0.2 - 1.5 * abs((abs(thetas[0][1]) - abs(thetas[0][0])))) < epsilon );
    TEST( abs( 0.2 - 1.5 * abs( thetas[1][1] - thetas[1][0] ) ) < epsilon );
    TEST( abs( 0.2 - 1.5 * abs( thetas[3][1] - thetas[3][0] ) ) < epsilon );

    RETURN_TEST;
}

int test_levelSet(){
    
    MPI_Init( 0, NULL);
    MPI_Comm_rank( MPI_COMM_WORLD, &parallel::rank );
    MPI_Comm_size( MPI_COMM_WORLD, &parallel::size );
    Form myform( 1, 4 );

    GridDimension N = { 3, 3};
    Point p0;
    Point p1;
    Point p2;
    Point p3;
    Grid grid( N, &myform );

    
    myform.coeff( 0 ) = 1;
    for ( uint_t i = 0; i < 4; ++i )
    {
        myform.begAngle( i ) = (i)*2. * EIT_PI / 4 - EIT_PI;
    }
    myform.buildThetas( 0.4 );
    // Build the Grid Points
    grid.buildPoints();

    grid.buildInterfaces();
    // Build the halos for communications
    grid.buildHalos();
    // Modify neighbors
    grid.correctConnections();
    p0.setCoordinates( { 0, 0 } );   // center
    p1.setCoordinates( { 1, 1 } ); // outside form
    p2.setCoordinates( { 1, 2 } );   // outside form
    p3.setCoordinates( { 1, 0 } ); // in form

    real_t p_0;
    real_t p_1;
    real_t p_2;
    p_0 = myform.levelSet( p0 ); // levelset of center is -r
    TEST( p_0 - (-1) < 0.01);
    p_1 = myform.levelSet( p1 );
    TEST(p_1 - (sqrt(2) - 1) < 0.01); 
    p_2 = myform.levelSet(p2);
    TEST(p_2 - 1.23 < 0.01);

    
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Finalize();
    
    RETURN_TEST;
}

int
test_derivativeLevelSet()
{
    
    MPI_Init( 0, NULL);
    MPI_Comm_rank( MPI_COMM_WORLD, &parallel::rank );
    MPI_Comm_size( MPI_COMM_WORLD, &parallel::size );
    real_t epsilon = 0.1;
    Form myform( 1, 4 );
    myform.coeff( 0 ) = 1;

    GridDimension N = { 3, 3 };
    Point p0;
    Point p1;
    Point p2;
    Grid grid( N, &myform );
    // Build the Grid Points
    grid.buildPoints();

    grid.buildInterfaces();
    // Build the halos for communications
    grid.buildHalos();
    // Modify neighbors
    grid.correctConnections();
    p0.setCoordinates( { 0, 0 } );   // center
    p1.setCoordinates( { 1, 0 } ); // in form
    p2.setCoordinates( { 1, 1 } );   // outside form

    real_t p_x;
    real_t p_y;
    p_x = myform.derivativeLevelSet( p0, 0, grid.getPads()[0] );
    p_y = myform.derivativeLevelSet( p0, 1, grid.getPads()[0] );
    TEST(abs(p_x - 0) < epsilon);
    TEST(abs(p_y - 0) < epsilon);
    //TEST( abs( p_x * p_x + p_y * p_y ) == 1 ); //||delta(levelset(p))|| = 1

    p_x = myform.derivativeLevelSet( p1, 0, grid.getPads()[0] ); 
    p_y = myform.derivativeLevelSet( p1, 1, grid.getPads()[1] );
    
    TEST(abs(p_x - 2) < epsilon);
    TEST(abs(p_y - 0) < epsilon);

    p_x = myform.derivativeLevelSet( p2, 0, grid.getPads()[0]);
    p_y = myform.derivativeLevelSet( p2, 1, grid.getPads()[1] );
    
    TEST(abs(p_x - 1.7) < epsilon);
    TEST(abs(p_y - 1.7) < epsilon);

    Form myform2( 7, 16 );

    real_t coeffs[7] = { 1.51, 0.01, 0.05, 0.2, 0.035, 0.01, 0.1 };
    for ( int i = 0; i < 7; i++ )
    {
        myform2.coeff( i ) = coeffs[i];
    }
    for ( int i = 0; i < 16; i++ )
    {
        myform2.begAngle( i ) = -EIT_PI + ( i - 1 ) * EIT_PI / 8;
    }

    p_x = myform2.derivativeLevelSet( p0, 0, grid.getPads()[0]);
    p_y = myform2.derivativeLevelSet( p0, 1, grid.getPads()[1]);
    
    TEST(abs(p_x - (-0.4)) < epsilon);
    TEST(abs(p_y - (0.13)) < epsilon);

    p_x = myform2.derivativeLevelSet( p1, 0, grid.getPads()[0]);
    p_y = myform2.derivativeLevelSet( p1, 1, grid.getPads()[1]);
    
    TEST(abs(p_x - 1.5) < epsilon);
    TEST(abs(p_y - -(0.04)) < epsilon);

    p_x = myform2.derivativeLevelSet( p2, 0, grid.getPads()[0] );
    p_y = myform2.derivativeLevelSet( p1, 1, grid.getPads()[1] );
    
    TEST(abs(p_x - 1.7) < epsilon);
    TEST(abs(p_y - (-0.04)) < epsilon);
    
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Finalize();
    
    RETURN_TEST;
}

int
test_computeInterfacePoint()
{
    real_t epsilon = 0.1;

    Form myform( 1, 4 );
    myform.coeff( 0 ) = 1.5;
   
    Point p;
    p.setCoordinates( { 1, 0 } );

    Point p_res = myform.computeInterfacePoint( p, 0, 1 );
    TEST( abs( 2 - p_res.getCoordinate( 0 ) ) < epsilon );
    TEST( abs( p.getCoordinate( 1 ) - p_res.getCoordinate( 1 ) ) < epsilon );

    Point p2;
    p2.setIJK( { 1, 2 } );

    Point p_res2 = myform.computeInterfacePoint( p2, 1, 1 );
    TEST( abs( p2.getCoordinate( 0 ) - p_res2.getCoordinate( 0 ) ) < epsilon );
    TEST( abs( 1.05 - p_res2.getCoordinate( 1 ) ) < epsilon );

    RETURN_TEST;
}

int ( *test_list[] )() = { test_getRaduis,
                           test_getFirstDerivativeRaduis,
                           test_getSecondDerivativeRaduis,
                           test_getEndAngleForLength,
                           test_buildThetas,
                           test_levelSet,
                           test_derivativeLevelSet,
                           test_computeInterfacePoint };

int
main( int argc, char **argv )
{
    
    int test_number = 1;
    if ( argc > 1 )
        test_number = atoi( argv[1] );

    return test_list[test_number - 1]();
}