#include "geometry/Point.hpp"

#include <set>

#include "tools/EitMath.hpp"

int_t
computeGlobalIndex( const GridDimension &N, const Vertex &ijk )
{
    // i
    // i*Ny + j
    // (i*Ny + j)*Nz + k

    int_t index = 0;
    for ( uint_t id = 0U; id < EIT_DIM; ++id )
    {
        if ( !isInteger(ijk[id]) ) return EIT_INVALID_INT;

        if ( ijk[id] >= N[id] )
            return EIT_INVALID_INT;
        index = index * N[id] + ijk[id];
    }

    if ( index < 0 )
        return EIT_INVALID_INT;
    return index;
}

Point::Point()
    : m_ijk( {} )
    , m_coor( {} )
    , m_neighs( {} )
    , m_isInterface( false )
    , m_isIrregular( false )
    , m_isOnElectrode( -1 )
    , m_normal( {} )
    , m_norStencil( {} )
    , m_place( GRID )
{
}

Point::Point( Vertex ijk )
    : m_ijk( ijk )
    , m_coor( {} )
    , m_neighs( {} )
    , m_isInterface( false )
    , m_isIrregular( false )
    , m_isOnElectrode( -1 )
    , m_normal( {} )
    , m_norStencil( {} )
    , m_place( GRID )
{
}

Point::Point( const Point &that )
    : m_ijk( that.m_ijk )
    , m_coor( that.m_coor )
    , m_neighs( that.m_neighs )
    , m_isInterface( that.m_isInterface )
    , m_isIrregular( that.m_isIrregular )
    , m_isOnElectrode( that.m_isOnElectrode )
    , m_normal( that.m_normal )
    , m_norStencil( that.m_norStencil )
    , m_place( that.m_place )
{
}

Point &
Point::operator=( const Point &that )
{
    if ( this != &that )
    {
        m_ijk           = that.m_ijk;
        m_coor          = that.m_coor;
        m_neighs        = that.m_neighs;
        m_isInterface   = that.m_isInterface;
        m_isIrregular   = that.m_isIrregular;
        m_isOnElectrode = that.m_isOnElectrode;
        m_normal        = that.m_normal;
        m_norStencil    = that.m_norStencil;
    }
    return *this;
}

Point::~Point() {}

// getting information on the object Point
const PointNeighbors &
Point::getNeighs() const
{
    return m_neighs;
}

const pointPlace &
Point::getPlace() const
{
    return m_place;
}

const Vertex &
Point::getCoordinates() const
{
    return m_coor;
}

const real_t &
Point::getCoordinate( uint i ) const
{
    return m_coor[i];
}

const Vertex &
Point::getNormal() const
{
    return m_normal;
}

const Vertex &
Point::getNormalStencil() const
{
    return m_norStencil;
}

const Vertex &
Point::getIJK() const
{
    return m_ijk;
}

bool
Point::isInterface() const
{
    return m_isInterface;
}

const int &
Point::getisOnElectrode() const
{
    return m_isOnElectrode;
}

bool
Point::isIrregular() const
{
    return m_isIrregular;
}

void
Point::setIJK( const Vertex &val )
{
    for ( uint_t i = 0; i < EIT_DIM; ++i )
        m_ijk[i] = val[i];
}

void
Point::setIJK( uint_t i, const real_t &val )
{
    m_ijk[i] = val;
}

void
Point::setCoordinates( const Vertex &val )
{
    for ( uint_t i = 0; i < EIT_DIM; ++i )
        m_coor[i] = val[i];
}

void
Point::setNormal( uint_t i, const real_t &val )
{
    m_normal[i] = val;
}

void
Point::setCoordinate( uint_t i, const real_t &val )
{
    m_coor[i] = val;
}

void
Point::setNormalStencil( uint_t i, const real_t &val )
{
    m_norStencil[i] = val;
}

void
Point::setCenter()
{
    for ( uint_t i = 0; i < EIT_DIM; ++i )
        m_coor[i] = 0;
}

real_t
Point::distance( const Point &center ) const
{
    return norm( m_coor - center.m_coor );
}

real_t
Point::angle() const
{
    return atan2( m_coor[1], m_coor[0] );
}

void
Point::getElectrodeIdForAngle( std::vector<ElectrodeThetas> &thetas )
{
    m_isOnElectrode = -1;
    if(norm(m_coor) <= EPSILON) return;

    for ( uint_t i = 0; i < EIT_ELEC; ++i )
    {
        if ( isInIntervalMod2PI(this->angle(), thetas[i]) )
        {
            m_isOnElectrode = i;
        }
    }

    return;
}

void
Point::setAsRegular( const GridDimension &N )
{
    m_place = GRID;
    std::array<real_t, EIT_DIM> previous, next;
    for ( uint_t id = 0U; id < EIT_DIM; ++id )
    {
        previous = next = m_ijk;
        next[id]++;
        previous[id]--;

        if ( previous[0] >= 0 && previous[1] >= 0 )
            m_neighs[2 * id] = computeGlobalIndex( N, previous );
        else
            m_neighs[2 * id] = EIT_INVALID_INT;

        if ( next[0] <= N[0] && next[1] <= N[1] )
            m_neighs[2 * id + 1] = computeGlobalIndex( N, next );
        else
            m_neighs[2 * id + 1] = EIT_INVALID_INT;
    }

    if ( m_neighs[0] == EIT_INVALID_INT || m_neighs[1] == EIT_INVALID_INT )
    {
        if ( m_neighs[2] == EIT_INVALID_INT || m_neighs[3] == EIT_INVALID_INT )
        {
            m_place = CORNER;
        }
        else
        {
            m_place = m_neighs[0] == EIT_INVALID_INT ? TOP : BOTTOM;
        }
    }
    else if ( m_neighs[2] == EIT_INVALID_INT )
    {
        m_place = LEFT;
    }
    else if ( m_neighs[3] == EIT_INVALID_INT )
    {
        m_place = RIGHT;
    }
    /* Looking for opinions, will be deleted on merge.
     * Is this clearer for you ? Does the same thing.

     switch(true) {
    case (m_neighs[0] == EIT_INVALID_INT || m_neighs[1] == EIT_INVALID_INT) &&
         (m_neighs[2] == EIT_INVALID_INT || m_neighs[3] == EIT_INVALID_INT):
      m_place = CORNER;
      break;
    case m_neighs[0] == EIT_INVALID_INT:
      m_place = TOP;
      break;
    case m_neighs[1] == EIT_INVALID_INT:
      m_place = BOTTOM;
      break;
    case m_neighs[2] == EIT_INVALID_INT:
      m_place = LEFT;
      break;
    case m_neighs[3] == EIT_INVALID_INT:
      m_place = RIGHT;
      break;
    default:
      // handle unexpected case
      break;
  }
     */

    m_isIrregular = false;
    m_isInterface = false;
    return;
}

void
Point::setAsIrregular( const PointNeighbors &neighs )
{
    m_neighs      = neighs;
    m_isIrregular = true;
    return;
}

void
Point::setAsInterface( const GridDimension &N )
{
    for ( uint_t id = 0U; id < EIT_DIM; ++id )
        m_neighs[2 * id] = m_neighs[2 * id + 1] = EIT_INVALID_INT; // no neighs

    m_place = INTERFACE;

    std::set<Vertex> unique;
    Vertex           previous, next;

    for ( uint_t dim = 0U; dim < EIT_DIM; ++dim )
    {
        if(m_ijk[dim] < 0 || m_ijk[dim] >= N[dim]){
            ERROR(
                std::cerr << "invalid coordinates (out of the grid) " << m_ijk << std::endl;
                MPI_Abort( MPI_COMM_WORLD, -1 );
            );
        }
        previous = next = m_ijk;
        
        previous[dim]   = std::floor( myAround( previous[dim] ) );
        next[dim]       = std::ceil( myAround( next[dim] ) );

        if ( norm( previous - m_ijk ) > EPSILON )
            unique.insert( previous );

        if ( norm( next - m_ijk ) > EPSILON )
            unique.insert( next );
    }

    uint_t nbToFind = std::pow( 2, EIT_DIM - 1 );

    if ( unique.size() != nbToFind )
    {
        ERROR(
        std::cerr << "error too much extremities detected on an interface. (Found " << unique.size()
                  << ")" << std::endl;
        for ( auto &p : unique )
            std::cout << p << std::endl;
        MPI_Abort( MPI_COMM_WORLD, -1 );
        );
    }

    for ( Vertex p : unique )
    {
        Vertex m   = p - m_ijk;
        uint_t dim = indexOfMaxAbs( m );

        if ( m[dim] < 0 )
            m_neighs[2 * dim] = computeGlobalIndex( N, p );
        else
            m_neighs[2 * dim + 1] = computeGlobalIndex( N, p );
    }

    m_isInterface = true;
}

Point &
Point::operator+=( const Point &b )
{
    this->m_coor = this->m_coor + b.m_coor;
    return *this;
}

Point &
Point::operator-=( const Point &b )
{
    m_coor = m_coor - b.m_coor;
    return *this;
}

Point &
Point::operator*=( const real_t &k )
{
    m_coor = m_coor * k;
    return *this;
}

Point
operator+( const Point &a, const Point &b )
{
    Point out = a;
    out += b;
    return out;
}

Point
operator-( const Point &a, const Point &b )
{
    Point out = a;
    out -= b;
    return out;
}

Point
operator*( const Point &a, const real_t &k )
{
    Point out = a;
    out *= k;
    return out;
}

Point
operator*( const real_t &k, const Point &a )
{
    return a * k;
}

std::ostream &
operator<<( std::ostream &os, const Point &p )
{
    os <<
#if EIT_DIM >= 1
        "[i = " << p.m_ijk[0] <<
#if EIT_DIM >= 2
        ", j = " << p.m_ijk[1] <<
#if EIT_DIM >= 3
        ", k = " << p.m_ijk[2] <<
#endif
#endif
#endif
        "] " << p.m_coor << std::flush;
    return os;
}
