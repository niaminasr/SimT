/**
 * @file Point.hpp
 * @brief Implementation of the Point class
 */

#ifndef __POINT_HPP__
#define __POINT_HPP__

#include <array>
#include <map>

#include "ParEITConfig.hpp"
#include "tools/EitArray.hpp"

typedef std::array<real_t, EIT_DIM> Vertex;
using GridDimension  = std::array<uint_t, EIT_DIM>;
using PointNeighbors = std::array<int_t, 2 * EIT_DIM>;

typedef enum pointPlace
{
    GRID,
    TOP,
    BOTTOM,
    LEFT,
    RIGHT,
    CORNER,
    INTERFACE
} pointPlace;


/**
 * @brief Computes the global index of a point given its IJK coordinates and the grid dimensions.
 * @param N the grid dimensions.
 * @param ijk the IJK coordinates of the point.
 * @return the global index of the point.
 * @retval EIT_INVALID_INT if one of the IJK coordinates is greater than or equal to its corresponding 
 *         dimension in N, or if the resulting index is negative.
 */
int_t computeGlobalIndex( const GridDimension &N, const Vertex &ijk );


/**
 * @class Point
 * @brief Creates of an object point that represents the point of a cartezian
 *        grid (x_i,y_j,z_k) up to the third dimension.
 * These points can be regular cartesian grid points. They also can be the intersection between 
 * the boundary and the geometrical form.
 * Let's call them "interface points".
 */
class Point
{
   public:
    /**
     * @brief Default constructor of a new Point object with default values.
     */
    Point();
    /**
     * @brief Constructs a new Point object with the given IJK coordinates.
     * @param ijk the IJK coordinates of the point.
     */
    Point( Vertex ijk );
    /**
     * @brief Copy constructor of a new Point object that is a copy of the given Point object.
     * @param that the Point object to copy.
     */
    Point( const Point &that );
    /**
     * @brief Copy assignment operator (Equalizing operator) for the Point class.
     * @param that the Point object to copy.
     * @return a reference to this Point object.
     */
    Point &operator=( const Point &that );
    /**
     * @brief Destructor for the Point class.
     */
    ~Point();

    /* getting information on the object Point */

    /**
     * @brief Returns the neighbors of an object Point.
     * @return The neighbors of the point.
     */
    const PointNeighbors &getNeighs() const;
    /**
     * @brief Returns the placement of a Point.
     * @return The type of the point.
     */
    const pointPlace &getPlace() const;
    /**
     * @brief Returns the coordinate array of an object Point.
     * @return The coordinates of the point.
     */
    const Vertex &getCoordinates() const;
    /**
     * @brief Returns the i_th coordinate of an object Point.
     * @param i The dimension of the coordinate to retrieve.
     * @return The coordinate value of the point for the specified dimension.
     */
    const real_t &getCoordinate( uint i ) const;
    /**
     * @brief Retuns the component of the normal vector at the object Point.
     * @return The normal vector of the point.
     */
    const Vertex &getNormal() const;
    /**
     * @brief Retuns the position (i,j,k) in cartesian grid of an object Point.
     * @return The IJK value of the point.
     */
    const Vertex &getIJK() const;
    /**
     * @brief Returns the normal stencil used in the discretization of the normal on the boundary. 
     * This is used for the interface Points. 
     * @return The normal vector stencil of the point.
     */
    const Vertex &getNormalStencil() const;
    /**
     * @brief Is this an interface Point ? 
     * @return True if the point is an interface point, false otherwise.
     */
    bool isInterface() const;
    /**
     * @brief Is my Point irregular ?
     * @return True if the point is an irregular point, false otherwise.
     */
    bool isIrregular() const;
    /**
     * @brief Is my point placed on an electrode ?
     * @return The electrode status of the point.
     */
    const int &getisOnElectrode() const;

    /* Setting information for an object Point */

    /**
     * @brief Sets the point's coordinates to the given vertex.
     * @param val The vertex to set the coordinates to.
     */
    void setCoordinates( const Vertex &val );
    /**
     * @brief Sets the i-th coordinate of the point's internal coordinate vector.
     * @param i The index of the coordinate to be set.
     * @param val The value to set the coordinate to.
     */
    void setCoordinate( uint_t i, const real_t &val );
    /**
     * @brief Sets the i-th coordinate (component) of the point's internal normal vector.
     * @param i The index of the coordinate to be set.
     * @param val The value to set the coordinate to.
     */
    void setNormal( uint_t i, const real_t &val );
    /**
     * @brief Sets the i-th coordinate (component) of the point's internal normal stencil vector.
     * @param i The index of the coordinate to be set.
     * @param val The value to set the coordinate to.
     */
    void setNormalStencil( uint_t i, const real_t &val );
    /**
     * @brief Get the number of the electrode on wich we are placed.
     * Determines the index of the electrode for which the angle of the point falls within 
     * the range of angles for that electrode.
     * @param thetas A vector containing the range of angles for each electrode
     */
    void getElectrodeIdForAngle( std::vector<ElectrodeThetas> &thetas );
    /**
     * @brief Sets the i-th coordinate of the point's internal 3D index.
     * @param val A vertex containing the i-th coordinate value.
     */
    void setIJK( const Vertex &val );
    /**
     * @brief Set the value of (i,j,k) individually (meaning component by component).
     * @param i The index of the coordinate to be set.
     * @param val The value to set the coordinate to.
     */
    void setIJK( uint_t i, const real_t &val );
    /**
     * @brief Sets the point's coordinates to zero (as the center).
     */
    void setCenter();
    /**
     * @brief Calculates the distance between this point and a given point.
     * @param center The point to calculate the distance to (the center Point).
     * @return The Euclidean distance between this point and the given point.
     */
    real_t distance( const Point &center ) const;
    /**
     * @brief Calculates the angle of the point in radians.
     * @return The angle associated to the point in radians.
     */
    real_t angle() const;
    /**
     * @brief Sets  m_isIrregular to false.
     * Sets the point's type to regular and calculates its neighbors based on its
     * internal 3D index and the given grid dimensions.
     * @param N The dimensions of the grid containing the point.
     */
    void setAsRegular( const GridDimension &N );
    /**
     * @brief Sets  m_isIrregular to true.
     * Sets the point's type to irregular and sets its neighbors to the given set of neighbors.
     * @param neighs The set of neighbors for the point
     */
    void setAsIrregular( const PointNeighbors &neighs );
    /**
     * @brief Sets m_isInterface to true.
     * Sets the point's type to interface and calculates its neighbors based on 
     * its internal 3D index and the given grid dimensions.
     * @param N The dimensions of the grid containing the point
     */
    void setAsInterface( const GridDimension &N );

    /* Operations on object Point */

    /**
     * @brief Add another point to this point.
     * @param b The point to add.
     * @return A reference to this point.
     */
    Point &operator+=( const Point &b );
    /**
     * @brief Subtract another point from this point.
     * @param b The point to subtract.
     * @return A reference to this point.
     */
    Point &operator-=( const Point &b );
    /**
     * @brief Multiply this point by a scalar value.
     * @param k The scalar value.
     * @return A reference to this point.
     */
    Point &operator*=( const real_t &k );
    /**
     * @brief Print an object Point.
     */
    friend std::ostream &operator<<( std::ostream &os, const Point &p );

   protected:
    Vertex         m_ijk;           // store (i,j,k)
    Vertex         m_coor;          // stores the coordinates of a Point
    PointNeighbors m_neighs;        // stores the neighbors of a Point, [North, South, East, West]
    bool           m_isInterface;   // a bool that indicates if our Point is on the interface or not
    bool           m_isIrregular;   // a bool that indicates the type of a grid Point: is it regular grid point or an irregular one
    int            m_isOnElectrode; // electrode number
    Vertex         m_normal, m_norStencil; // stores normal values and normal stencils of a Point
    pointPlace     m_place; // stores the position in the grid of a Point
};

/* Operations on object Point with two entries */

/**
 * @brief Add two points together.
 * @param a The first point.
 * @param b The second point.
 * @return The sum of the two points.
 */
Point operator+( const Point &a, const Point &b );
/**
 * @brief Subtract one point from another.
 * @param a The first point.
 * @param b The second point.
 * @return The difference between the two points.
 */
Point operator-( const Point &a, const Point &b );
/**
 * @brief Multiply a point by a scalar value.
 * @param a The point to multiply.
 * @param k The scalar value.
 * @return The scaled point.
 */
Point operator*( const Point &a, const real_t &k );
/**
 * @brief Multiply a point by a scalar value.
 * @param k The scalar value.
 * @param a The point to multiply.
 * @return The scaled point.
 */
Point operator*( const real_t &k, const Point &a );
/**
 * @brief Write a point to an output stream.
 * @param os The output stream.
 * @param p The point to write.
 * @return The output stream.
 */
std::ostream &operator<<( std::ostream &os, const Point &p );

#endif /* __POINT_HPP__ */
