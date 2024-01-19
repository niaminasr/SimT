#include "tools/EitVector.hpp"

#include <cmath>

EitVector::EitVector() : m_data( new std::vector<real_t>() ) {}

EitVector::EitVector( uint_t n ) : m_data( new std::vector<real_t>() )
{
    resize( n );
    fill( 0.0 );
}

EitVector::EitVector( const EitVector &that ) : m_data( new std::vector<real_t>( *that.m_data ) ) {}

EitVector::~EitVector()
{
    if ( m_data )
        delete m_data;
}

EitVector &
EitVector::operator=( const EitVector &that )
{
    if ( this != &that )
        *m_data = *that.m_data;
    return *this;
}

uint_t
EitVector::size() const
{
    return m_data->size();
}

void
EitVector::resize( uint_t n )
{
    m_data->resize( n );
    return;
}

real_t &
EitVector::operator[]( uint_t index )
{
    return ( *m_data )[index];
}

const real_t &
EitVector::operator[]( uint_t index ) const
{
    return ( *m_data )[index];
}

real_t &
EitVector::at( uint_t index )
{
    return m_data->at( index );
}

const real_t &
EitVector::at( uint_t index ) const
{
    return m_data->at( index );
}

void
EitVector::fill( real_t value )
{
    for ( real_t &v : *this )
        v = value;
    return;
}

EitVector
EitVector::reciproqual() const
{
    EitVector temp = *this;
    for ( real_t &v : temp )
        v = 1.0 / v;
    return temp;
}

EitVector
EitVector::zero() const
{
    EitVector temp = *this;
    for ( real_t &v : temp )
        v = 0;
    return temp;
}

real_t
EitVector::dotProduct( const EitVector &that ) const
{
    checkSizes( __func__, size(), that.size() );

    uint_t n = size();
    real_t v = 0.0;

    for ( uint_t i = 0; i < n; ++i )
        v += ( *this )[i] * that[i];

    return v;
}

real_t
EitVector::norml2() const
{
    return sqrt( this->dotProduct( *this ) );
}

typename EitVector::iterator
EitVector::begin()
{
    return m_data->begin();
}

typename EitVector::const_reverse_iterator
EitVector::rbegin() const
{
    return m_data->rbegin();
}

typename EitVector::const_iterator
EitVector::begin() const
{
    return m_data->begin();
}

typename EitVector::const_iterator
EitVector::cbegin() const
{
    return m_data->cbegin();
}

typename EitVector::const_reverse_iterator
EitVector::crbegin() const
{
    return m_data->crbegin();
}

typename EitVector::iterator
EitVector::end()
{
    return m_data->end();
}

typename EitVector::const_reverse_iterator
EitVector::rend() const
{
    return m_data->rend();
}

typename EitVector::const_iterator
EitVector::end() const
{
    return m_data->end();
}

typename EitVector::const_iterator
EitVector::cend() const
{
    return m_data->cend();
}

typename EitVector::const_reverse_iterator
EitVector::crend() const
{
    return m_data->crend();
}

void
checkSizes( std::string func, uint_t n, uint_t m )
{
    if ( n != m )
    {
        std::cerr << "Assert the sizes are different in call of " << func << "(" << m << "," << n
                  << ")" << std::endl;
        abort();
    }
}

real_t *
EitVector::data()
{
    return m_data->data();
}

const real_t *
EitVector::data() const
{
    return m_data->data();
}

EitVector &
EitVector::operator+=( const EitVector &that )
{
    checkSizes( __func__, size(), that.size() );

    uint_t n = size();
    for ( uint_t i = 0; i < n; ++i )
        ( *this )[i] += that[i];
    return *this;
}

EitVector &
EitVector::operator+=( const real_t &coeff )
{
    uint_t n = size();
    for ( uint_t i = 0; i < n; ++i )
        ( *this )[i] += coeff;
    return *this;
}

EitVector &
EitVector::operator-=( const EitVector &that )
{
    checkSizes( __func__, size(), that.size() );

    uint_t n = size();
    for ( uint_t i = 0; i < n; ++i )
        ( *this )[i] -= that[i];
    return *this;
}

EitVector &
EitVector::operator-=( const real_t &coeff )
{
    uint_t n = size();
    for ( uint_t i = 0; i < n; ++i )
        ( *this )[i] -= coeff;
    return *this;
}

EitVector &
EitVector::operator*=( const EitVector &that )
{
    checkSizes( __func__, size(), that.size() );

    uint_t n = size();
    for ( uint_t i = 0; i < n; ++i )
        ( *this )[i] *= that[i];
    return *this;
}

EitVector &
EitVector::operator*=( const real_t &coeff )
{
    uint_t n = size();
    for ( uint_t i = 0; i < n; ++i )
        ( *this )[i] *= coeff;
    return *this;
}

EitVector &
EitVector::operator/=( const EitVector &that )
{
    checkSizes( __func__, size(), that.size() );

    uint_t n = size();
    for ( uint_t i = 0; i < n; ++i )
        if ( abs( that[i] ) > EIT_EPSILON )
            ( *this )[i] /= that[i];

    return *this;
}

EitVector &
EitVector::operator/=( const real_t &coeff )
{
    if ( abs( coeff ) < EIT_EPSILON )
        return *this;

    uint_t n = size();
    for ( uint_t i = 0; i < n; ++i )
        ( *this )[i] /= coeff;
    return *this;
}

EitVector
operator+( const EitVector &a, const EitVector &b )
{
    EitVector temp = a;
    temp += b;
    return temp;
}

EitVector
operator+( const EitVector &a, const real_t &coeff )
{
    EitVector temp = a;
    temp += coeff;
    return temp;
}

EitVector
operator+( const real_t &coeff, const EitVector &a )
{
    EitVector temp = a;
    temp += coeff;
    return temp;
}

EitVector
operator-( const EitVector &a, const EitVector &b )
{
    EitVector temp = a;
    temp -= b;
    return temp;
}

EitVector
operator-( const EitVector &a, const real_t &coeff )
{
    EitVector temp = a;
    temp -= coeff;
    return temp;
}

EitVector
operator-( const real_t &coeff, const EitVector &a )
{
    EitVector temp = a * -1.0;
    temp += coeff;
    return temp;
}

EitVector
operator*( const EitVector &a, const EitVector &b )
{
    EitVector temp = a;
    temp *= b;
    return temp;
}

EitVector
operator*( const EitVector &a, const real_t &coeff )
{
    EitVector temp = a;
    temp *= coeff;
    return temp;
}

EitVector
operator*( const real_t &coeff, const EitVector &a )
{
    EitVector temp = a;
    temp *= coeff;
    return temp;
}

EitVector
operator/( const EitVector &a, const EitVector &b )
{
    EitVector temp = a;
    temp /= b;
    return temp;
}

EitVector
operator/( const EitVector &a, const real_t &coeff )
{
    EitVector temp = a;
    temp /= coeff;
    return temp;
}

EitVector
operator/( const real_t &coeff, const EitVector &a )
{
    EitVector temp = a.reciproqual();
    temp *= coeff;
    return temp;
}