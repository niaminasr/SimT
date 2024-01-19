#include "geometry/Form.hpp"
EIT_DECL_PARALLEL_DEFINITIONS

Form::Form() : m_center( new Point() ), m_alphas( {} ), m_thetas( {} ) {}

Form::Form( uint_t n, uint_t ne ) : m_center( new Point() ), m_alphas( {} ), m_thetas( {} )
{
    m_alphas.resize( n, 0. );
    m_center->setCenter();
    m_thetas.resize( ne );
}

Form::Form( const Form &that )
    : m_center( new Point( *that.m_center ) ), m_alphas( that.m_alphas ), m_thetas( that.m_thetas )
{
}

real_t &
Form::coeff( uint_t i )
{
    return m_alphas[i];
}

real_t &
Form::begAngle( uint_t i )
{
    return m_thetas[i][0];
}

uint_t
Form::size() const
{
    return m_alphas.size();
}

real_t
Form::radius( real_t theta )
{
    uint_t n   = ( size() - 1 ) / 2;
    real_t out = coeff( 0 );

    for ( uint_t i = 1; i <= n; i++ )
    {
        out += coeff( i ) * cos( i * theta );
        out += coeff( n + i ) * sin( i * theta );
    }

    return out;
}

real_t
Form::firstDerivativeofRadius( real_t theta )
{
    uint_t n   = ( size() - 1 ) / 2;
    real_t out = coeff( 0 );

    for ( uint_t i = 1; i <= n; i++ )
    {
        out -= i * coeff( i ) * sin( i * theta );
        out += i * coeff( n + i ) * cos( i * theta );
    }

    return out;
}

real_t
Form::secondDerivativeofRadius( real_t theta )
{
    uint_t n   = ( size() - 1 ) / 2;
    real_t out = coeff( 0 );

    for ( uint_t i = 1; i <= n; i++ )
    {
        out -= i * i * coeff( i ) * cos( i * theta );
        out -= i * i * coeff( n + i ) * sin( i * theta );
    }

    return out;
}

real_t
Form::rho( real_t theta )
{
    double r  = radius( theta );
    double dr = firstDerivativeofRadius( theta );
    return sqrt( r * r + dr * dr );
}

real_t
Form::getEndAngleForLength( real_t L, uint_t i )
{
    // number of time subdivision
    uint_t n = 100;

    real_t F  = L;
    real_t t  = m_thetas[i][0] - L / 0.5;
    real_t dt = 0.0;

    for ( uint_t k = 0; ( k < EIT_ITMAX ) && ( F > EIT_EPSILON ); k++ )
    {
        // Firstly, we compute the approximation of (L-integral_{tm-->t}) using the midpoints
        // algorithm
        dt = ( t - m_thetas[i][0] ) / n;
        F  = L;
        for ( uint_t j = 0; j < n; j++ )
            F -= dt * rho( m_thetas[i][0] + ( j + 1 / 2 ) * dt );

        // Secondly, we compute second end of the electrode t
        t += F / rho( t );
    }
    return t;
}

void
Form::buildThetas( real_t L )
{
    for ( uint_t i = 0; i < m_thetas.size(); ++i )
    {
        m_thetas[i][1] = this->getEndAngleForLength( L, i );
    }
}

std::vector<ElectrodeThetas> &
Form::getThetas()
{
    return m_thetas;
}

real_t
Form::levelSet( const Point &pt )
{
    real_t theta = pt.angle();
    real_t dist  = pt.distance( *m_center );

    return ( dist - radius( theta ) );
}

real_t
Form::operator()( const Point &pt )
{
    return levelSet( pt );
}

real_t
Form::derivativeLevelSet( const Point &pt, real_t dim, real_t step )
{
    real_t Dxlevel, Dylevel, norm = 0;
    Point  pt_minusx, pt_plusx, pt_minusy, pt_plusy;

    Vertex coor, coor_minusx, coor_plusx, coor_minusy, coor_plusy;

    coor = pt.getCoordinates();

    coor_minusx    = coor;
    coor_minusx[0] = coor[0] - step;

    coor_plusx    = coor;
    coor_plusx[0] = coor[0] + step;

    coor_minusy    = coor;
    coor_minusy[1] = coor[1] - step;

    coor_plusy    = coor;
    coor_plusy[1] = coor[1] + step;

    pt_minusx.setCoordinates( coor_minusx );
    pt_plusx.setCoordinates( coor_plusx );
    pt_minusy.setCoordinates( coor_minusy );
    pt_plusy.setCoordinates( coor_plusy );

    Dxlevel = ( levelSet( pt_plusx ) - levelSet( pt_minusx ) ) / 2 * step;
    Dylevel = ( levelSet( pt_plusy ) - levelSet( pt_minusy ) ) / 2 * step;

    norm = sqrt( Dxlevel * Dxlevel + Dylevel * Dylevel );

    if ( dim == 0 )
    {
        return Dxlevel;
    }
    else if ( dim == 1 )
    {
        return Dylevel;
    }
}

Point
Form::computeInterfacePoint( const Point &a, uint_t dim, real_t dir )
{
    real_t Fbeg, Fmid;
    Point  beg, mid, end;
    beg = a;
    end = a;
    end.setCoordinate( dim, end.getCoordinate( dim ) + dir );

    for ( uint_t n = 0; n < 1000; ++n )
    {
        mid = ( beg + end ) * ( 1. / 2. );

        Fbeg = beg.distance( *m_center ) - radius( beg.angle() );
        Fmid = mid.distance( *m_center ) - radius( mid.angle() );

        if ( Fbeg * Fmid < 0. )
            end.setCoordinate( dim, mid.getCoordinate( dim ) );
        else
            beg.setCoordinate( dim, mid.getCoordinate( dim ) );
    }

    return mid;
}

real_t
Form::approximateIntegral( const Point &pt, real_t &step )
{
    Point  pt_minusx, pt_plusx, pt_minusy, pt_plusy;
    real_t term;

    Vertex coor, coor_minusx, coor_plusx, coor_minusy, coor_plusy;

    coor = pt.getCoordinates();

    coor_minusx    = coor;
    coor_minusx[0] = coor[0] - step;

    coor_plusx    = coor;
    coor_plusx[0] = coor[0] + step;

    coor_minusy    = coor;
    coor_minusy[1] = coor[1] - step;

    coor_plusy    = coor;
    coor_plusy[1] = coor[1] + step;

    pt_minusx.setCoordinates( coor_minusx );
    pt_plusx.setCoordinates( coor_plusx );
    pt_minusy.setCoordinates( coor_minusy );
    pt_plusy.setCoordinates( coor_plusy );

    real_t dxphiplus  = ( levelSet( pt_plusx ) - levelSet( pt ) ) / step;
    real_t dxphimoins = ( levelSet( pt ) - levelSet( pt_minusx ) ) / step;
    real_t dxphizero  = ( levelSet( pt_plusx ) - levelSet( pt_minusx ) ) / ( 2. * step );

    real_t dyphiplus  = ( levelSet( pt_plusy ) - levelSet( pt ) ) / step;
    real_t dyphimoins = ( levelSet( pt ) - levelSet( pt_minusy ) ) / step;
    real_t dyphizero  = ( levelSet( pt_plusy ) - levelSet( pt_minusy ) ) / ( 2. * step );

    real_t absgrad = sqrt( dxphizero * dxphizero + dyphizero * dyphizero + EIT_EPSILON );

    real_t delta = 0.;
    if ( levelSet( pt ) * levelSet( pt_plusx ) <= 0. )
    {
        term = std::abs( levelSet( pt_plusx ) * dxphizero ) /
               ( step * step * std::abs( dxphiplus * absgrad ) );
        delta = delta + term;
    }

    if ( levelSet( pt ) * levelSet( pt_minusx ) < 0. )
    {
        term = std::abs( levelSet( pt_minusx ) * dxphizero ) /
               ( step * step * std::abs( dxphimoins * absgrad ) );
        delta = delta + term;
    }

    if ( levelSet( pt ) * levelSet( pt_plusy ) <= 0. )
    {
        term = std::abs( levelSet( pt_plusy ) * dyphizero ) /
               ( step * step * std::abs( dyphiplus * absgrad ) );
        delta = delta + term;
    }

    if ( levelSet( pt ) * levelSet( pt_minusy ) < 0. )
    {
        term = std::abs( levelSet( pt_minusy ) * dyphizero ) /
               ( step * step * std::abs( dyphimoins * absgrad ) );
        delta = delta + term;
    }

    return delta;
}

std::ostream &
operator<<( std::ostream &os, const Form &f )
{
    os << "Form: " << std::endl;
    os << " order : " << f.size() << std::endl;
    for ( uint_t i = 0; i < f.size(); i++ )
    {
        os << " coefs : {" << f.m_alphas[i] << "}" << std::endl; // transpose a vector.
    }
    return os;
}
