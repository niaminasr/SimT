/**
 * @file EitArray.hpp
 * @brief This file defines the EitArray class and several helper functions for std::array manipulation.
 */

#ifndef __EITARRAY_HPP__
#define __EITARRAY_HPP__

#include <array>
#include <cmath>
#include <iostream>

#include "ParEITConfig.hpp"

/**
 * @brief The EitArray class represents a std::array of any data type with a fixed size.
 * @tparam _Type The data type of the std::array.
 * @tparam _Dim The size of the std::array.
 */
template <class _Type, std::size_t _Dim>
using EitArray = std::array<_Type, _Dim>;

/**
 * @brief The EitArrayInt class represents a std::array of signed integers with a fixed size.
 * @tparam _Dim The size of the std::array.
 */
template <std::size_t _Dim>
using EitArrayInt = EitArray<int_t, _Dim>;

/**
 * @brief The EitArrayUInt class represents a std::array of unsigned integers with a fixed size.
 * @tparam _Dim The size of the std::array.
 */
template <std::size_t _Dim>
using EitArrayUInt = EitArray<std::size_t, _Dim>;

/**
 * @brief The EitArrayReal class represents a std::array of real numbers with a fixed size.
 * @tparam _Dim The size of the std::array.
 */
template <std::size_t _Dim>
using EitArrayReal = EitArray<real_t, _Dim>;

/**
 * @brief Adds two std::arrays element-wise.
 * @tparam _Type The data type of the std::arrays.
 * @tparam _Dim The size of the std::arrays.
 * @param[in] a The first std::array.
 * @param[in] b The second std::array.
 * @return The result of the addition.
 */

template <class _Type, std::size_t _Dim>
std::array<_Type, _Dim>
operator+( const std::array<_Type, _Dim> &a, const std::array<_Type, _Dim> &b )
{
    std::array<_Type, _Dim> out = a;
    for ( std::size_t i = 0; i < _Dim; ++i )
        out[i] += b[i];
    return out;
}

/**
 * @brief Subtracts two std::arrays element-wise.
 * @tparam _Type The data type of the std::arrays.
 * @tparam _Dim The size of the std::arrays.
 * @param[in] a The first std::array.
 * @param[in] b The second std::array.
 * @return The result of the subtraction.
 */
template <class _Type, std::size_t _Dim>
std::array<_Type, _Dim>
operator-( const std::array<_Type, _Dim> &a, const std::array<_Type, _Dim> &b )
{
    std::array<_Type, _Dim> out = a;
    for ( std::size_t i = 0; i < _Dim; ++i )
        out[i] -= b[i];
    return out;
}

/**
 * @brief Adds a scalar to a std::array element-wise.
 * @tparam _Type The data type of the std::array.
 * @tparam _Dim The size of the std::array.
 * @param[in] a The std::array.
 * @param[in] k The scalar.
 * @return The result of the addition.
 */
template <class _Type, std::size_t _Dim>
std::array<_Type, _Dim>
operator+( const std::array<_Type, _Dim> &a, const _Type &k )
{
    std::array<_Type, _Dim> out = a;
    for ( std::size_t i = 0; i < _Dim; ++i )
        out[i] += k;
    return out;
}

/**
 * @brief Adds the scalar value to each element of the array.
 * @tparam _Type The data type of the elements in the array.
 * @tparam _Dim The size of the array.
 * @param k The scalar value to add to the array.
 * @param a The array to which the scalar value is to be added.
 * @return std::array<_Type, _Dim> The resulting array after adding the scalar value.
 */
template <class _Type, std::size_t _Dim>
std::array<_Type, _Dim>
operator+( const _Type &k, const std::array<_Type, _Dim> &a )
{
    std::array<_Type, _Dim> out = a;
    for ( std::size_t i = 0; i < _Dim; ++i )
        out[i] += k;
    return out;
}

/**
 * @brief Subtracts the scalar value from each element of the array.
 * @tparam _Type The data type of the elements in the array.
 * @tparam _Dim The size of the array.
 * @param a The array from which the scalar value is to be subtracted.
 * @param k The scalar value to subtract from the array.
 * @return std::array<_Type, _Dim> The resulting array after subtracting the scalar value.
 */
template <class _Type, std::size_t _Dim>
std::array<_Type, _Dim>
operator-( const std::array<_Type, _Dim> &a, const _Type &k )
{
    std::array<_Type, _Dim> out = a;
    for ( std::size_t i = 0; i < _Dim; ++i )
        out[i] -= k;
    return out;
}

/**
 * @brief Subtracts each element of the array from the scalar value.
 * @tparam _Type The data type of the elements in the array.
 * @tparam _Dim The size of the array.
 * @param k The scalar value from which the array is to be subtracted.
 * @param a The array to be subtracted from the scalar value. 
 * @return std::array<_Type, _Dim> The resulting array after subtracting each element of the array from the scalar value.
 */
template <class _Type, std::size_t _Dim>
std::array<_Type, _Dim>
operator-( const _Type &k, const std::array<_Type, _Dim> &a )
{
    std::array<_Type, _Dim> out = a;
    for ( std::size_t i = 0; i < _Dim; ++i )
        out[i] = k - out[i];
    return out;
}

/**
 * @brief Computes the dot product of two arrays.
 *
 * @tparam _Type The type of the array elements.
 * @tparam _Dim The number of dimensions of the arrays.
 * @param a The first array.
 * @param b The second array.
 * @return The dot product of the two input arrays.
 */
template <class _Type, std::size_t _Dim>
_Type
dot( const std::array<_Type, _Dim> &a, const std::array<_Type, _Dim> &b )
{
    _Type out = _Type();
    for ( std::size_t i = 0; i < _Dim; ++i )
        out += a[i] * b[i];
    return out;
}

/**
 * @brief Computes the Euclidean norm (magnitude) of an array.
 *
 * @tparam _Type The type of the array elements.
 * @tparam _Dim The number of dimensions of the array.
 * @param a The array to compute the norm of.
 * @return The Euclidean norm of the input array.
 */
template <class _Type, std::size_t _Dim>
auto
norm( const std::array<_Type, _Dim> &a ) -> decltype( std::sqrt( dot( a, a ) ) )
{
    return std::sqrt( dot( a, a ) );
}

/**
 * @brief Computes the product of a scalar value and an array.
 *
 * @tparam _Type The type of the array elements.
 * @tparam _Dim The number of dimensions of the array.
 * @param k The scalar value to multiply the array by.
 * @param a The array to multiply by the scalar value.
 * @return A new array that is the product of the scalar value and the input array.
 */
template <class _Type, std::size_t _Dim>
std::array<_Type, _Dim>
operator*( const _Type &k, const std::array<_Type, _Dim> &a )
{
    std::array<_Type, _Dim> out = a;
    for ( std::size_t i = 0; i < _Dim; ++i )
        out[i] *= k;
    return out;
}

/**
 * @brief Multiplies an array to a scalar element-wise
 * @tparam _Type The data type of the array
 * @tparam _Dim The dimension of the array
 * @param a The array
 * @param k The scalar
 * @return The resulting array after scalar multiplication
 */
template <class _Type, std::size_t _Dim>
std::array<_Type, _Dim>
operator*( const std::array<_Type, _Dim> &a, const _Type &k )
{
    return k * a;
}

/**
 * @brief Divides a scalar to an array element-wise
 * @tparam _Type The data type of the array
 * @tparam _Dim The dimension of the array
 * @param k The scalar
 * @param a The array
 * @return The resulting array after scalar division
 */
template <class _Type, std::size_t _Dim>
std::array<_Type, _Dim>
operator/( const _Type &k, const std::array<_Type, _Dim> &a )
{
    std::array<_Type, _Dim> out = a;
    for ( std::size_t i = 0; i < _Dim; ++i )
        out[i] = k / out[i];
    return out;
}

/**
 * @brief Divide each element of an array by a scalar value.
 * @tparam _Type The type of the array elements.
 * @tparam _Dim The size of the array.
 * @param a The array to divide.
 * @param k The scalar value to divide by.
 * @return The resulting array.
 */
template <class _Type, std::size_t _Dim>
std::array<_Type, _Dim>
operator/( const std::array<_Type, _Dim> &a, const _Type &k )
{
    return a * ( 1.0 / k );
}

/**
 * @brief Divide two arrays element-wise.
 *
 * @tparam _Type The type of the array elements.
 * @tparam _Dim The size of the arrays.
 * @param a The first array.
 * @param b The second array.
 * @return The resulting array.
 */
template <class _Type, std::size_t _Dim>
std::array<_Type, _Dim>
operator/( const std::array<_Type, _Dim> &a, const std::array<_Type, _Dim> &b )
{
    std::array<_Type, _Dim> out;
    for ( std::size_t i = 0; i < _Dim; ++i )
        out[i] = a[i] / b[i];
    return out;
}

/**
 * @brief Multiply two arrays element-wise.
 *
 * @tparam _Type The type of the array elements.
 * @tparam _Dim The size of the arrays.
 * @param a The first array.
 * @param b The second array.
 * @return The resulting array.
 */
template <class _Type, std::size_t _Dim>
std::array<_Type, _Dim>
operator*( const std::array<_Type, _Dim> &a, const std::array<_Type, _Dim> &b )
{
    std::array<_Type, _Dim> out = a;
    for ( std::size_t i = 0; i < _Dim; ++i )
        out[i] *= b[i];
    return out;
}

/**
 * @brief Compute the product of all elements in an array.
 *
 * @tparam _Type The type of the array elements.
 * @tparam _Dim The size of the array.
 * @param a The array.
 * @return The product of all elements in the array.
 */
template <class _Type, std::size_t _Dim>
_Type
product( const std::array<_Type, _Dim> &a )
{
    if ( _Dim == 0U )
        return _Type( 0 );

    _Type out = a[0];

    for ( std::size_t i = 1U; i < _Dim; ++i )
        out *= a[i];
    return out;
}

/**
 * @brief Find the index of the element with maximum absolute value in an array.
 *
 * @tparam _Type The type of the array elements.
 * @tparam _Dim The size of the array.
 * @param a The array.
 * @return The index of the element with maximum absolute value in the array.
 */
template <class _Type, std::size_t _Dim>
std::size_t
indexOfMaxAbs( const std::array<_Type, _Dim> &a )
{
    if ( _Dim == 0 )
    {
        std::cerr << "can not find the maximum pos" << std::endl;
        MPI_Abort( MPI_COMM_WORLD, -1 );
    }

    std::size_t imax = 0U;
    _Type       vmax = std::abs( a[imax] );

    for ( uint_t i = 1U; i < _Dim; ++i )
    {
        if ( vmax < std::abs( a[i] ) )
        {
            imax = i;
            vmax = a[i];
        }
    }
    return imax;
}

/**
 * @brief Converts each element of an array to a different type.
 * 
 * This function takes an array of type `_Type` and converts each element to type `_Target`.
 * The resulting array will have the same size as the input array.
 * 
 * @tparam _Target The target type to convert each element to.
 * @tparam _Type The type of the input array.
 * @tparam _Dim The size of the input and output arrays.
 * @param a The input array to convert.
 * @return std::array<_Target, _Dim> The resulting array with each element converted to type `_Target`.
 */
template <class _Target, class _Type, std::size_t _Dim>
std::array<_Target, _Dim>
convert( const std::array<_Type, _Dim> &a )
{
    std::array<_Target, _Dim> out;
    for ( std::size_t i = 0; i < _Dim; ++i )
        out[i] = static_cast<_Target>( a[i] );
    return out;
}

/**
 * @brief Prints an array to an output stream.
 * 
 * This operator takes an array of type `_Type` and prints it to an output stream.
 * The elements of the array are separated by spaces, enclosed by parentheses and flushed at the end.
 * 
 * @tparam _Type The type of the input array.
 * @tparam _Dim The size of the input array.
 * @param os The output stream to print to.
 * @param a The input array to print.
 * @return std::ostream& A reference to the output stream.
 */
template <class _Type, std::size_t _Dim>
std::ostream &
operator<<( std::ostream &os, const std::array<_Type, _Dim> &a )
{
    os << "(";
    for ( uint_t i = 0; i < _Dim - 1; ++i )
        os << a[i] << " ";
    os << a[_Dim - 1] << ")" << std::flush;
    return os;
}

#endif /* __EITARRAY_HPP__ */
