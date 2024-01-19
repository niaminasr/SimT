#include "tests.hpp"

// Don't compile with MPI if you forget this line
// EIT_DECL_PARALLEL_DEFINITIONS;s

int
Point_getElectrodeIdForAngle()
{
    // Compute the index of corresponding electrode with point's angle
    // The result is obtained with Point::getisOnElectrode
    // And in 3D ?

    // Dependencies :
    //- Point::setCoordinates
    //- Point::getisOnElectrode

    double                       epsilon = 0.1;
    std::vector<ElectrodeThetas> thetas  = {
        { 0.0, EIT_PI / 4.0 },                                      //[0, pi/4]
        { EIT_PI / 2.0, EIT_PI / 2.0 },                             //[pi/2, pi/2]
        { -( EIT_PI / 2.0 ) - epsilon, 0.0 - epsilon },             //[(-pi/2)-eps, 0-eps]
        { EIT_PI - ( 2 * EIT_PI ), -( EIT_PI / 2.0 ) + epsilon } }; //[pi-2pi, (-pi/2)+eps]
    Point p;

    p.setCoordinates( { 0.0, 0.0 } ); // angle ???
    p.getElectrodeIdForAngle( thetas );
    TEST( p.getisOnElectrode() == -1 );

    p.setCoordinates( { 1.0, 0.0 } ); // angle = 0
    p.getElectrodeIdForAngle( thetas );
    TEST( p.getisOnElectrode() == 0 );

    p.setCoordinates( { 1.0, 1.0 } ); // angle = pi/4
    p.getElectrodeIdForAngle( thetas );
    TEST( p.getisOnElectrode() == 0 );

    p.setCoordinates( { 0.0, 1.0 } ); // angle = pi/2
    p.getElectrodeIdForAngle( thetas );
    TEST( p.getisOnElectrode() == 1 );

    p.setCoordinates( { -1.0, 1.0 } ); // angle = 3pi/4
    p.getElectrodeIdForAngle( thetas );
    TEST( p.getisOnElectrode() == -1 );

    p.setCoordinates( { -1.0, 0.0 } ); // angle = pi
    p.getElectrodeIdForAngle( thetas );
    TEST( p.getisOnElectrode() == 3 );

    p.setCoordinates( { 1.0, -1.0 } ); // angle = -pi/4
    p.getElectrodeIdForAngle( thetas );
    TEST( p.getisOnElectrode() == 2 );

    p.setCoordinates( { 0.0, -1.0 } ); // angle = -pi/2
    p.getElectrodeIdForAngle( thetas );
    TEST( p.getisOnElectrode() == 2 || p.getisOnElectrode() == 3 );

    p.setCoordinates( { -1.0, -1.0 } ); // angle = -3pi/4
    p.getElectrodeIdForAngle( thetas );
    TEST( p.getisOnElectrode() == 3 );

    RETURN_TEST;
}

int
Point_computeGlobalIndex()
{
    // Return an unique index for a grid element (starting from 0 on top left)
    // With integers grid coordinates ijk in grid boundaries (if not return EIT_INVALID_INT)

    GridDimension N = { 10, 3 }; // 2D : lines columns

    TEST( computeGlobalIndex( N, { 0, 0 } ) == 0 );  // origin (north and west boundary)
    TEST( computeGlobalIndex( N, { 0, 1 } ) == 1 );  // regular
    TEST( computeGlobalIndex( N, { 0, 2 } ) == 2 );  // east boundary
    TEST( computeGlobalIndex( N, { 1, 1 } ) == 4 );  // regular 2
    TEST( computeGlobalIndex( N, { 5, 0 } ) == 15 ); // west center boundary
    TEST( computeGlobalIndex( N, { 9, 2 } ) == 29 ); // last (south and east boundary)

    TEST( computeGlobalIndex( N, { 0, 3 } ) == EIT_INVALID_INT );   // east out
    TEST( computeGlobalIndex( N, { -1, 0 } ) == EIT_INVALID_INT );  // north out
    TEST( computeGlobalIndex( N, { 10, -1 } ) == EIT_INVALID_INT ); // south and west out
    TEST( computeGlobalIndex( N, { 0.5, 0 } ) == EIT_INVALID_INT ); // not integer

    RETURN_TEST;
}

int
Point_setAsRegular()
{
    // Set regular bool
    // Set interface bool
    // Compute neighbors
    // And set type

    // For the moment we only test in 2D

    // Dependencies :
    //- Point::setIJK
    //- Point::getNeighs
    //- Point::getType
    //- Point::isIrregular
    //- Point::isInterface

    GridDimension  N = { 3, 5 };
    Point          p;
    PointNeighbors neighs;

    p.setIJK( { 1, 2 } ); // type = GRID
    p.setAsRegular( N );
    neighs = p.getNeighs();
    TEST( p.isIrregular() == false ); // test booleans : regular and interface
    TEST( p.isInterface() == false );
    TEST( p.getPlace() == GRID );
    TEST( neighs[0] == 0 * 5 + 2 );
    TEST( neighs[1] == 2 * 5 + 2 );
    TEST( neighs[2] == 1 * 5 + 1 );
    TEST( neighs[3] == 1 * 5 + 3 );

    p.setIJK( { 2, 4 } ); // type = CORNER
    p.setAsRegular( N );
    neighs = p.getNeighs();
    TEST( p.getPlace() == CORNER );
    TEST( neighs[0] == 1 * 5 + 4 );
    TEST( neighs[2] == 2 * 5 + 3 );
    TEST( neighs[1] == EIT_INVALID_INT && neighs[3] == EIT_INVALID_INT );

    p.setIJK( { 0, 2 } ); // type = TOP (north)
    p.setAsRegular( N );
    neighs = p.getNeighs();
    TEST( p.getPlace() == TOP );
    TEST( neighs[0] == EIT_INVALID_INT );
    TEST( neighs[1] == 1 * 5 + 2 );
    TEST( neighs[2] == 0 * 5 + 1 );
    TEST( neighs[3] == 0 * 5 + 3 );

    p.setIJK( { 2, 2 } ); // type = BOTTOM (south)
    p.setAsRegular( N );
    neighs = p.getNeighs();
    TEST( p.getPlace() == BOTTOM );
    TEST( neighs[0] == 1 * 5 + 2 );
    TEST( neighs[1] == EIT_INVALID_INT );
    TEST( neighs[2] == 2 * 5 + 1 );
    TEST( neighs[3] == 2 * 5 + 3 );

    p.setIJK( { 1, 0 } ); // type = LEFT (west)
    p.setAsRegular( N );
    neighs = p.getNeighs();
    TEST( p.getPlace() == LEFT );
    TEST( neighs[0] == 0 * 5 + 0 );
    TEST( neighs[1] == 2 * 5 + 0 );
    TEST( neighs[2] == EIT_INVALID_INT );
    TEST( neighs[3] == 1 * 5 + 1 );

    p.setIJK( { 1, 4 } ); // type = RIGHT (east)
    p.setAsRegular( N );
    neighs = p.getNeighs();
    TEST( p.getPlace() == RIGHT );
    TEST( neighs[0] == 0 * 5 + 4 );
    TEST( neighs[1] == 2 * 5 + 4 );
    TEST( neighs[2] == 1 * 5 + 3 );
    TEST( neighs[3] == EIT_INVALID_INT );

    p.setIJK( { 3, -1 } ); // invalid point (out of grid)
    p.setAsRegular( N );
    neighs = p.getNeighs();
    // TEST(???)

    p.setIJK( { 1, 0.5 } ); // invalid point (not integer grid coords -> interface)
    p.setAsRegular( N );
    neighs = p.getNeighs();
    // TEST(???)

    RETURN_TEST;
}

int
Point_setAsInterface()
{
    // Set interface bool
    // Compute neighbors
    // Set type to interface
    // 3D ?

    // Dependencies :
    //- Point::setIJK
    //- Point::getNeighs
    //- Point::getType
    //- Point::isInterface

    GridDimension  N = { 3, 4 };
    Point          p;
    PointNeighbors neighs;

    p.setIJK( { 1, 1.5 } ); // regular (center)
    p.setAsInterface( N );
    neighs = p.getNeighs();
    TEST( p.isInterface() == true );  // test boolean : interface,
    TEST( p.getPlace() == INTERFACE ); // and type
    TEST( neighs[0] == EIT_INVALID_INT );
    TEST( neighs[1] == EIT_INVALID_INT );
    TEST( neighs[2] == 1 * 4 + 1 );
    TEST( neighs[3] == 1 * 4 + 2 );

    p.setIJK( { 1.7, 3 } ); // regular (east boundary)
    p.setAsInterface( N );
    neighs = p.getNeighs();
    TEST( neighs[0] == 1 * 4 + 3 );
    TEST( neighs[1] == 2 * 4 + 3 );
    TEST( neighs[2] == EIT_INVALID_INT );
    TEST( neighs[3] == EIT_INVALID_INT );

    p.setIJK( { 0, 0.2 } ); // regular (north west corner)
    p.setAsInterface( N );
    neighs = p.getNeighs();
    TEST( neighs[0] == EIT_INVALID_INT );
    TEST( neighs[1] == EIT_INVALID_INT );
    TEST( neighs[2] == 0 * 4 + 0 );
    TEST( neighs[3] == 0 * 4 + 1 );

    p.setIJK( { 1, 2 } ); // invalid (regular point / not interface)
    p.setAsInterface( N );
    TEST( isErr() );

    p.setIJK( { 0.6, -1 } ); // invalid (out of grid < 0)
    p.setAsInterface( N );
    TEST( isErr() );

    p.setIJK( { 3, 0.2 } ); // invalid (out of grid > N)
    p.setAsInterface( N );
    TEST( isErr() );

    p.setIJK( { 1.1, 2.9 } ); // invalid (more than one grid coord not integer / cannot be an interface)
    p.setAsInterface( N );
    TEST( isErr() );

    RETURN_TEST;
}

int ( *test_list[] )() = { Point_getElectrodeIdForAngle, Point_computeGlobalIndex,
                           Point_setAsRegular, Point_setAsInterface };

// Check LastTest.log for sub-testing results

int
main( int argc, char **argv )
{
    int test_number = 0;
    if ( argc > 1 )
        test_number = atoi( argv[1] );

    return test_list[test_number - 1]();
}