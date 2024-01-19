/**
 * @file EitBiVector.hpp
 * @brief Header file for the EitBiVector class.
 */

#ifndef __EITBIVECTOR_HPP__
#define __EITBIVECTOR_HPP__

#include <algorithm>

#include "ParEITConfig.hpp"

/**
 * @class EitBiVector
 * @brief Bi-vector class.
 *
 * The EitBiVector class is a bi-vector class that contains two contiguous arrays
 * of different sizes, one for positive indices and one for negative indices.
 * It is used to store solution and data arrays in the EIT solver.
 *
 * @tparam _Type The type of data stored in the bi-vector.
 */
template <class _Type>
class EitBiVector
{
   public:
    /**
     * @brief Constructor for the EitBiVector class.
     * @param Nneg The size of the negative index array.
     * @param Npos The size of the positive index array.
     */
    EitBiVector( uint_t Nneg, uint_t Npos )
        : m_pos( Npos ), m_neg( Nneg ), m_size( Nneg + Npos ), m_data( nullptr )
    {
        resize( Nneg, Npos );
    }
    /**
     * @brief Constructor for the EitBiVector class.
     * @param Nneg The size of the negative index array.
     * @param Npos The size of the positive index array.
     * @param value The initial value for all entries in the bi-vector.
     */
    EitBiVector( uint_t Nneg, uint_t Npos, const _Type &value )
        : m_pos( Npos ), m_neg( Nneg ), m_size( Nneg + Npos ), m_data( nullptr )
    {
        resize( Nneg, Npos, value );
    }
    /**
     * @brief Copy constructor for the EitBiVector class.
     * @param that The EitBiVector object to be copied.
     */
    EitBiVector( const EitBiVector &that )
        : m_pos( that.m_pos ), m_neg( that.m_neg ), m_size( that.m_size ), m_data( nullptr )
    {
        m_data = new _Type[m_size];
        std::copy_n( that.m_data, m_size, m_data );
    }
    /**
     * @brief Assignment operator for the EitBiVector class.
     * @param that The EitBiVector object to be assigned.
     * @return A reference to the assigned object.
     */
    EitBiVector &
    operator=( const EitBiVector &that )
    {
        if ( this != &that )
        {
            m_pos  = that.m_pos;
            m_neg  = that.m_neg;
            m_size = m_neg + m_pos;

            delete[] m_data;
            m_data = new _Type[m_size];
            std::copy_n( that.m_data, m_size, m_data );
        }
        return *this;
    }
    /**
     * @brief Destructor for the EitBiVector class.
     */
    ~EitBiVector() { delete[] m_data; }
    /**
     * @brief Overload of the [] operator for the EitBiVector class.
     * @param i The index of the entry to be returned.
     * @return A reference to the requested entry.
     */
    _Type &
    operator[]( int_t i )
    {
        return m_data[m_neg + i];
    }
    /**
     * @brief Overload of the [] operator for the EitBiVector class.
     * @param i The index of the entry to be returned.
     * @return A const reference to the requested entry.
     */
    const _Type &
    operator[]( int_t i ) const
    {
        return m_data[m_neg + i];
    }
    /**
     * @brief Method to access a specific entry in the bi-vector.
     * @param i The index of the entry to be accessed.
     * @return A reference to the requested entry.
     */
    _Type &
    at( int_t i )
    {
        assert( ( -i <= (int_t)m_neg ) and ( i < (int_t)m_pos ) );
        return m_data[m_neg + i];
    }
    /**
     * @brief Method to access a specific entry in the bi-vector.
     * @param i The index of the entry to be accessed.
     * @return A const reference to the requested entry.
     */
    const _Type &
    at( int_t i ) const
    {
        // assert ((-i <= m_neg) and (i < m_pos));
        return m_data[m_neg + i];
    }
    /**
     * @brief Method to get the size of the bi-vector.
     * @return A const reference to the size of the bi-vector.
     */
    const uint_t &
    size() const
    {
        return m_size;
    }
    /**
     * @brief Method to get the size of the positive index array.
     * @return A const reference to the size of the positive index array.
     */
    const uint_t &
    posSize() const
    {
        return m_pos;
    }
    /**
     * @brief Method to get the size of the negative index array.
     * @return A const reference to the size of the negative index array.
     */
    const uint_t &
    negSize() const
    {
        return m_neg;
    }
    /**
     * @brief Method to resize the vector, preserving its content.
     * @param Nneg The new size of the negative index array.
     * @param Npos The new size of the positive index array.
     */
    void
    resize( const uint_t &Nneg, const uint_t &Npos )
    {
        _Type *data = new _Type[Nneg + Npos];

        if ( m_data != nullptr )
        {
            _Type *middle = m_data + m_neg;
            _Type *beg    = middle - std::min( m_neg, Nneg );
            _Type *end    = middle + std::min( m_pos, Npos );
            std::copy( beg, end, data );

            delete[] m_data;
        }

        m_data = data;
        return;
    }
    /**
     * @brief Method to resize the vector and fill it with a given value.
     * @param Nneg The new size of the negative index array.
     * @param Npos The new size of the positive index array.
     * @param value The value to fill the vector with.
     */
    void
    resize( const uint_t &Nneg, const uint_t &Npos, const _Type &value )
    {
        resize( Nneg, Npos );
        fill( value );
        return;
    }
    /**
     * @brief Method to get a pointer to the underlying data array.
     * @return A pointer to the underlying data array.
     */
    _Type *
    data()
    {
        return m_data;
    }
    /**
     * @brief Method to get a const pointer to the underlying data array.
     * @return A const pointer to the underlying data array.
     */
    const _Type *
    data() const
    {
        return m_data;
    }
    /**
     * @brief Method to get a pointer to the positive index part of the underlying data array.
     * @return A pointer to the positive index part of the underlying data array.
     */
    _Type *
    dataPos()
    {
        return m_data + m_neg;
    }
    /**
     * @brief Method to get a const pointer to the positive index part of the underlying data array.
     * @return A const pointer to the positive index part of the underlying data array.
     */
    const _Type *
    dataPos() const
    {
        return m_data + m_neg;
    }
    /**
     * @brief Method to get a pointer to the negative index part of the underlying data array.
     * @return A pointer to the negative index part of the underlying data array.
     */
    _Type *
    dataNeg()
    {
        return m_data;
    }
    /**
     * @brief Method to get a const pointer to the negative index part of the underlying data array.
     * @return A const pointer to the negative index part of the underlying data array.
     */
    const _Type *
    dataNeg() const
    {
        return m_data;
    }
    /**
     * @brief Method to add another vector to this one.
     * @param a The vector to be added.
     * @return A reference to this vector after the addition operation.
     */
    EitBiVector &
    operator+=( const EitBiVector &a )
    {
        assert( size() == a.size() );
        for ( uint_t i = 0U; i < size(); ++i )
            m_data[i] += a.m_data[i];
        return *this;
    }
    /**
     * @brief Method to subtract another vector from this one.
     * @param a The vector to be subtracted.
     * @return A reference to this vector after the subtraction operation.
     */
    EitBiVector &
    operator-=( const EitBiVector &a )
    {
        assert( size() == a.size() );
        for ( uint_t i = 0U; i < size(); ++i )
            m_data[i] -= a.m_data[i];
        return *this;
    }
    /**
     * @brief Method to multiply this vector by a scalar factor.
     * @param factor The scalar factor to multiply the vector by.
     * @return A reference to this vector after the multiplication operation.
     */
    EitBiVector &
    operator*=( const _Type &factor )
    {
        for ( uint_t i = 0U; i < size(); ++i )
            m_data[i] *= factor;
        return *this;
    }
    /**
     * @brief Method to fill the vector with a given value.
     * @param v The value to fill the vector with.
     */
    void
    fill( const _Type &v )
    {
        for ( uint_t i = 0U; i < size(); ++i )
            m_data[i] = v;
        return;
    }
    /**
     * @brief Method to set all elements of the vector to zero.
     */
    void
    zero()
    {
        return fill( (_Type)0 );
    }
    /**
     * @brief Method to set all elements of the vector to zero.
     */
    void
    ones()
    {
        return fill( (_Type)1 );
    }

   protected:
    uint_t m_pos; /**< The size of the positive index array. */
    uint_t m_neg; /**< The size of the negative index array. */
    uint_t m_size; /**< The total size of the vector. */
    _Type *m_data; /**< A pointer to the underlying data array. */
};

/**
 * @brief Calculates the dot product of two arrays of type `_Type`.
 * @tparam _Type The type of elements in the arrays.
 * @param a The first array.
 * @param b The second array.
 * @param n The number of elements in the arrays.
 * @param reduce If true, the result of the dot product is reduced (summed) across all processes
 * in MPI_COMM_WORLD. Default is false.
 * @return The dot product of `a` and `b`.
 */
template <class _Type>
_Type
dot( const _Type *a, const _Type *b, const uint_t &n, bool reduce )
{
    _Type value = _Type();
    for ( uint_t i = 0U; i < n; ++i )
        value += a[i] * b[i];
    if ( reduce )
        MPI_Allreduce( MPI_IN_PLACE, &value, 1, EitDataType<_Type>::get(), MPI_SUM,
                       MPI_COMM_WORLD );
    return value;
}

/**
 * @brief Calculates the dot product of two `EitBiVector` objects.
 * @tparam _Type The type of elements in the `EitBiVector` objects.
 * @param a The first `EitBiVector` object.
 * @param b The second `EitBiVector` object.
 * @return The dot product of `a` and `b`.
 */
template <class _Type>
_Type
specialDot( const EitBiVector<_Type> &a, const EitBiVector<_Type> &b )
{
    assert( a.posSize() == b.posSize() );
    assert( a.negSize() == b.negSize() );

    _Type value = dot( a.dataPos(), b.dataPos(), a.posSize(), true );
    value += dot( a.dataNeg(), b.dataNeg(), a.negSize(), false );

    return value;
}

/**
 * @brief Calculates the Euclidean norm of an `EitBiVector` object.
 * @tparam _Type The type of elements in the `EitBiVector` object.
 * @param a The `EitBiVector` object.
 * @return The Euclidean norm of `a`.
 */
template <class _Type>
_Type
norm( const EitBiVector<_Type> &a )
{
    return std::sqrt( specialDot( a, a ) );
}

/**
 * @brief Adds two `EitBiVector` objects and returns the result.
 * @tparam _Type The type of elements in the `EitBiVector` objects.
 * @param a The first `EitBiVector` object.
 * @param b The second `EitBiVector` object.
 * @return A new `EitBiVector` object that is the sum of `a` and `b`.
 */
template <class _Type>
EitBiVector<_Type>
operator+( const EitBiVector<_Type> &a, const EitBiVector<_Type> &b )
{
    EitBiVector<_Type> out = a;
    out += b;
    return out;
}

/**
 * @brief Subtracts one `EitBiVector` object from another and returns the result.
 * @tparam _Type The type of elements in the `EitBiVector` objects.
 * @param a The `EitBiVector` object to subtract from.
 * @param b The `EitBiVector` object to subtract.
 * @return A new `EitBiVector` object that is the difference between `a` and `b`.
 */
template <class _Type>
EitBiVector<_Type>
operator-( const EitBiVector<_Type> &a, const EitBiVector<_Type> &b )
{
    EitBiVector<_Type> out = a;
    out -= b;
    return out;
}

/**
 * @brief Multiplies an `EitBiVector` object by a scalar factor and returns the result.
 * @tparam _Type The type of elements in the `EitBiVector` object.
 * @param factor The scalar factor to multiply with.
 * @param b The `EitBiVector` object to multiply.
 * @return A new `EitBiVector` object that is the result of multiplying `b` by `factor`.
 */
template <class _Type>
EitBiVector<_Type>
operator*( const _Type &factor, const EitBiVector<_Type> &b )
{
    EitBiVector<_Type> out = b;
    out *= factor;
    return out;
}

/**
 * @brief Multiplies each element of an EitBiVector by a scalar value.
 * @tparam _Type The type of elements in the EitBiVector.
 * @param a The EitBiVector to multiply.
 * @param factor The scalar value to multiply each element by.
 * @return EitBiVector<_Type> The resulting EitBiVector.
 */
template <class _Type>
EitBiVector<_Type>
operator*( const EitBiVector<_Type> &a, const _Type &factor )
{
    return factor * a;
}

#endif /* __EITBIVECTOR_HPP__ */
