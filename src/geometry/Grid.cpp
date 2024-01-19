#include "geometry/Grid.hpp"


Grid::Grid( GridDimension N, Form *form )
    : m_N( N )
    , m_NxBeg( 0U )
    , m_NxEnd( N[0] )
    , m_h( {} )
    , m_points( {} )
    , m_halo_left( {} )
    , m_halo_right( {} )
    , m_interfaces( {} )
    , m_form( form )
    , m_sizeproc( {} )
{
    for ( uint_t i = 0; i < EIT_DIM; i++ )
        m_h[i] = 4.0 / (real_t)( m_N[i] - 1); // the square -2-->2

    //int value = m_N[1] % ( parallel::size );
    //m_NxBeg = parallel::rank * ( m_N[1] / ( parallel::size ) ) + std::min( parallel::rank, value );
    //m_NxEnd = m_NxBeg + ( m_N[1] / ( parallel::size ) ) - static_cast<int>( parallel::rank >= value );
    int value = m_N[1] ;
    m_NxBeg = parallel::rank * ( m_N[1] ) + std::min( parallel::rank, value );
    m_NxEnd = m_NxBeg + ( m_N[1] ) - static_cast<int>( parallel::rank >= value );
}

// !+=========================================================
// MEMORY LEAK HERE, EACH POINT* IS NOW SHARED BY TWO DIFFERENT GRID
// !+=========================================================
Grid::Grid( const Grid &that )
    : m_N( that.m_N )
    , m_NxBeg( that.m_NxBeg )
    , m_NxEnd( that.m_NxEnd )
    , m_h( that.m_h )
    , m_points( that.m_points )
    , m_halo_left( that.m_halo_left )
    , m_halo_right( that.m_halo_right )
    , m_interfaces( that.m_interfaces )
    , m_form( that.m_form )
    , m_sizeproc( that.m_sizeproc )
{
}

Grid &
Grid::operator=( const Grid &that )
{
    if ( this != &that )
    {
        m_N          = that.m_N;
        m_NxBeg      = that.m_NxBeg;
        m_NxEnd      = that.m_NxEnd;
        m_h          = that.m_h;
        m_points     = that.m_points;
        m_interfaces = that.m_interfaces;
        m_form       = that.m_form;
        m_sizeproc   = that.m_sizeproc;
    }
    return *this;
}

Grid::~Grid()
{
    for ( auto &l : { m_points, m_interfaces, m_halo_left, m_halo_right } )
        for ( Point *p : l )
            delete p;
}

const GridDimension &
Grid::getSizes() const
{
    return m_N;
}

const std::array<real_t, EIT_DIM> &
Grid::getPads() const
{
    return m_h;
}

uint_t
Grid::getLocalBeg() const
{
    return m_NxBeg;
}

uint_t
Grid::getLocalEnd() const
{
    return m_NxEnd;
}

uint_t
Grid::getNumberOfPoints() const
{
    return m_points.size();
}

uint_t
Grid::getNumberOfInterfaces() const
{
    return m_interfaces.size();
}

const std::vector<Point *> &
Grid::getInterfaces() const
{
    return m_interfaces;
}

Point *
Grid::getPoint( uint_t i )
{
    return m_points[i];
}

Point *
Grid::getInterface( uint_t i )
{
    return m_interfaces[i];
}

const std::vector<Point *> &
Grid::getPoints() const
{
    return m_points;
}

const std::vector<Point *> &
Grid::getHaloL() const
{
    return m_halo_left;
}

const std::vector<Point *> &
Grid::getHaloR() const
{
    return m_halo_right;
}

std::vector<uint_t>
Grid::getSizeAllProc()
{
    return m_sizeproc;
}

void
Grid::buildPoints()
{
    uint   NxLoc = m_NxEnd - m_NxBeg + 1;
    uint_t size  = NxLoc;
    for ( uint_t i = 1U; i < EIT_DIM; ++i )
        size *= m_N[i];

    // Intialize a vector of pointers :: i have to fill it with nullpts
    m_points.reserve( size );

    Vertex ijk, coor;
   //faire gaffe quand c'est parallèle a m_NxEnd doit etre inclu dans le truc
    for ( ijk[0] = m_NxBeg; ijk[0] < m_NxEnd; ++ijk[0] )
    {
#if EIT_DIM >= 2
        for ( ijk[1] = 0; ijk[1] < m_N[1]; ++ijk[1] )
#endif
        {
#if EIT_DIM >= 3
            for ( ijk[2] = 0; ijk[2] < m_N[2]; ++ijk[2] )
#endif
            {
                Point *p = new Point( ijk );
                toCoordinate( ijk, coor );
                p->setCoordinates( coor );
                p->setAsRegular( m_N );
                m_points.push_back( p );
            }
        }
    }


    return;
}

void
Grid::buildHalos()
{
    uint_t size = 1;
    for ( uint_t i = 1U; i < EIT_DIM; ++i )
        size *= m_N[i];

    m_halo_left.reserve( size );
    m_halo_right.reserve( size );

    Vertex ijk, coor;

    if ( parallel::rank > 0 )
    {
        ijk[0] = m_NxBeg - 1;

#if EIT_DIM >= 2
        for ( ijk[1] = 0; ijk[1] < m_N[1]; ++ijk[1] )
#endif
        {
#if EIT_DIM >= 3
            for ( ijk[2] = 0; ijk[2] < m_N[2]; ++ijk[2] )
#endif
            {
                Point *p = new Point( ijk );
                toCoordinate( ijk, coor );
                p->setCoordinates( coor );
                p->setAsRegular( m_N );
                m_halo_left.push_back( p );
            }
        }
    }
    std::cout << "halo left size : " << m_halo_left.size() << "  " << EIT_PROC_RANK << std::endl;

    if ( parallel::rank < parallel::size - 1 )
    {
        ijk[0] = m_NxEnd + 1;

#if EIT_DIM >= 2
        for ( ijk[1] = 0; ijk[1] < m_N[1]; ++ijk[1] )
#endif
        {
#if EIT_DIM >= 3
            for ( ijk[2] = 0; ijk[2] < m_N[2]; ++ijk[2] )
#endif
            {
                Point *p = new Point( ijk );
                toCoordinate( ijk, coor );
                p->setCoordinates( coor );
                p->setAsRegular( m_N );
                m_halo_right.push_back( p );
            }
        }
    }

    std::cout << "halo right size : " << m_halo_right.size() << "  " << EIT_PROC_RANK << std::endl;

    return;
}

void
Grid::buildInterfaces()
{
    for ( Point *p : m_interfaces )
        delete p;
    m_interfaces.clear();
    /**
     * @brief An enum to keep track of the point indices in the arrays `p`, `phi`,
     * `ijk`, and `coor`. `cur` represents the current point, `neigh` represents
     * the neighboring point, and `itf` represents the interface point.
     */
    enum // c'est une structure ? non
    {
        cur   = 0,
        neigh = 1,
        itf   = 2
    };

    // 0: current point, 1: p+h, 2: itf
    Point  p[3];
    real_t phi[3];
    Vertex ijk[3];
    Vertex coor[3];

    for ( uint_t dim = 0U; dim < EIT_DIM; ++dim )
    {
        ijk[cur] = Vertex();

        // all grid must be studied ? Yes
        for ( ijk[cur][0] = 0; ijk[cur][0] < m_N[0]; ++ijk[cur][0] )
        {
#if EIT_DIM >= 2
            for ( ijk[cur][1] = 0; ijk[cur][1] < m_N[1]; ++ijk[cur][1] )
#endif
            {
#if EIT_DIM >= 3
                for ( ijk[cur][2] = 0; ijk[cur][2] < m_N[2]; ++ijk[cur][2] )
#endif
                {
                    toCoordinate( ijk[cur], coor[cur] );
                    p[cur].setCoordinates( coor[cur] );
                    phi[cur] = m_form->levelSet( p[cur] );

                    ijk[neigh] = ijk[cur];
                    ijk[neigh][dim]++;

                    if ( ijk[neigh][dim] >= 0 )
                    {
                        toCoordinate( ijk[neigh], coor[neigh] );
                        p[neigh].setCoordinates( coor[neigh] );
                        phi[neigh] = m_form->levelSet( p[neigh] );

                        if ( phi[cur] * phi[neigh] < 0 )
                        {
                            // Is this enaugh or i have to do a dichotomie ?
                            // phi = a t + b
                            // gamma(t) = cur * t + (1-t)* neigh
                            real_t a = ( phi[cur] - phi[neigh] ) / ( 0.0 - 1.0 );
                            real_t b = phi[cur] - a * 0.0;
                            real_t t = -b / a;

                            coor[itf] = t * coor[cur] + ( 1.0 - t ) * coor[neigh];
                            toIJK( coor[itf], ijk[itf] );
                            p[itf].setIJK( ijk[itf] );
                            p[itf].setCoordinates( coor[itf] );
                            p[itf].setAsInterface( m_N );
                            p[itf].getElectrodeIdForAngle( m_form->getThetas() );

                            if ( ijk[cur][1] ==
                                 ijk[neigh][1] ) // interface point in the x direction.
                            {
                                real_t levelx = m_form->derivativeLevelSet( p[cur], 0, m_h[0] ) *
                                                    ( coor[neigh][0] - coor[itf][0] ) -
                                                m_form->derivativeLevelSet( p[neigh], 0, m_h[0] ) *
                                                    ( coor[cur][0] - coor[itf][0] );

                                real_t levely = m_form->derivativeLevelSet( p[cur], 1, m_h[1] ) *
                                                    ( coor[neigh][1] - coor[itf][0] ) -
                                                m_form->derivativeLevelSet( p[neigh], 1, m_h[1] ) *
                                                    ( coor[cur][1] - coor[itf][0] );

                                real_t norm = sqrt( levelx * levelx + levely * levely );

                                levelx = levelx / norm;
                                levely = levely / norm;

                                p[itf].setNormal( 0, levelx );
                                p[itf].setNormal( 1, levely );

                                real_t stencil_option[6][2] = {
                                    { coor[cur][0], coor[cur][1] },
                                    { coor[cur][0], coor[cur][1] + m_h[1] },
                                    { coor[neigh][0], coor[neigh][1] + m_h[1] },
                                    { coor[neigh][0], coor[neigh][1] },
                                    { coor[neigh][0], coor[neigh][1] - m_h[1] },
                                    { coor[cur][0], coor[cur][1] - m_h[1] } };
                                uint_t id = 0;
                                for ( int_t i = 0; i < 6; ++i )
                                {
                                    real_t x = ( coor[itf][0] - stencil_option[i][0] );
                                    real_t y = ( coor[itf][1] - stencil_option[i][1] );

                                    real_t n_x = x * levelx;
                                    real_t n_y = y * levely;

                                    Vertex v;

                                    if ( n_x >= 0 && n_y >= 0 )
                                    {
                                        v[0] = ( stencil_option[i][0] + 2.0 ) / m_h[0];
                                        v[1] = ( stencil_option[i][1] + 2.0 ) / m_h[1];

                                        p[itf].setNormalStencil( id, v[0] * m_N[0] + v[1] );
                                        id += 1;
                                    }
                                }
                            }
                            else
                            {
                                real_t levelx = m_form->derivativeLevelSet( p[cur], 0, m_h[0] ) *
                                                    ( coor[neigh][1] - coor[itf][1] ) -
                                                m_form->derivativeLevelSet( p[neigh], 0, m_h[0] ) *
                                                    ( coor[cur][1] - coor[itf][1] );

                                real_t levely = m_form->derivativeLevelSet( p[cur], 1, m_h[1] ) *
                                                    ( coor[neigh][1] - coor[itf][1] ) -
                                                m_form->derivativeLevelSet( p[neigh], 1, m_h[1] ) *
                                                    ( coor[cur][1] - coor[itf][1] );

                                real_t norm = sqrt( levelx * levelx + levely * levely );
                                levelx      = levelx / norm;
                                levely      = levely / norm;

                                p[itf].setNormal( 0, levelx );
                                p[itf].setNormal( 1, levely );

                                real_t stencil_option[6][2] = {
                                    { coor[cur][0], coor[cur][1] },
                                    { coor[cur][0] - m_h[0], coor[cur][1] },
                                    { coor[neigh][0] - m_h[0], coor[neigh][1] },
                                    { coor[neigh][0], coor[neigh][1] },
                                    { coor[neigh][0] + m_h[0], coor[neigh][1] },
                                    { coor[cur][0] + m_h[1], coor[cur][1] } };
                                uint_t id = 0;
                                for ( int_t i = 0; i < 6; ++i )
                                {
                                    real_t x = ( coor[itf][0] - stencil_option[i][0] );
                                    real_t y = ( coor[itf][1] - stencil_option[i][1] );

                                    real_t n_x = x * levelx;
                                    real_t n_y = y * levely;

                                    Vertex v;

                                    if ( n_x >= 0 && n_y >= 0 )
                                    {
                                        v[0] = ( stencil_option[i][0] + 2.0 ) / m_h[0];
                                        v[1] = ( stencil_option[i][1] + 2.0 ) / m_h[1];

                                        p[itf].setNormalStencil( id, v[0] * m_N[0] + v[1] );
                                        id += 1;
                                    }
                                }
                            }

                            m_interfaces.push_back( new Point( p[itf] ) );
                        }
                    }
                }
            }
        }
    }

    return;
}

void
Grid::correctConnections()
{
    // Firstly, we need to reset as regular all grid points
    for ( auto &l : { m_points, m_halo_left, m_halo_right } )
        for ( Point *p : l )
            p->setAsRegular( m_N );

    Vertex Nmin, Nmax;
    Nmin[0] = m_NxBeg;
    Nmax[0] = m_NxEnd;
    for ( uint_t dim = 1U; dim < EIT_DIM; ++dim )
    {
        Nmin[dim] = 0;
        Nmax[dim] = m_N[dim] - 1;
    }

    int_t min_gID = computeGlobalIndex( m_N, Nmin );
    int_t max_gID = computeGlobalIndex( m_N, Nmax );

    // std::cout << min_gID << "  " << max_gID << "  " << EIT_PROC_RANK << std::endl;

    uint_t n_itfs = getNumberOfInterfaces();
    for ( uint_t itf = 0U; itf < n_itfs; ++itf )
    {
        Point                *p      = m_interfaces[itf];
        const PointNeighbors &neighs = p->getNeighs();

        for ( uint_t dim = 0U; dim < EIT_DIM; ++dim )
        {
            for ( int_t id : { neighs[2 * dim], neighs[2 * dim + 1] } )
            {
                if ( ( id != -1 ) and ( min_gID <= id ) and ( id <= max_gID ) )
                {
                    id -= min_gID;

                    Point *p_neigh = m_points[id];

                    PointNeighbors n_neigh = p_neigh->getNeighs();

                    if ( p_neigh->getCoordinate( dim ) < p->getCoordinate( dim ) )
                        n_neigh[2 * dim + 1] = -itf - 1; // negative number for interfaces
                    else
                        n_neigh[2 * dim] = -itf - 1; // negative number for interfaces

                    p_neigh->setAsIrregular( n_neigh );
                }
            }
        }
    }
};

void
Grid::setForm( Form *form )
{
    m_form = form;
    return;
}

Form *
Grid::getForm() const
{
    return m_form;
}

void
Grid::toCoordinate( const Vertex &ijk, Vertex &coor )
{
    coor = (ijk)*m_h - 2.0;
    return;
}

void
Grid::toIJK( const Vertex &coor, Vertex &ijk )
{
    ijk = ( coor + 2.0 ) / m_h;
    return;
}
// Bradcast to all processors the number of nodes managed by each proc and stoch it into a array of
// size ::: parallel::size
void
Grid::bcastSizes()
{
    int nproc = parallel::size;
    m_sizeproc.resize( nproc );

    m_sizeproc[parallel::rank] = this->getNumberOfPoints();
    for ( uint_t i = 0; i < parallel::size; ++i )
    {
        MPI_Bcast( &m_sizeproc[i], 1, EIT_MPI_INT, i, MPI_COMM_WORLD );
    }
}