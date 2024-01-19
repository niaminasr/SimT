/**
 * @file EitMath.cpp
 * @brief Implementation file for math helper functions used by EIT simulations.
 */

#ifndef __EITMATH_HPP__
#define __EITMATH_HPP__

#include <cmath>
#include <functional>

#include "ParEITConfig.hpp"

/**
 * @brief Rounds a real number to the nearest integer.
 * @param x The real number to round.
 * @return The nearest integer to x.
 *
 * This function rounds a real number to the nearest integer using the "round half up" method.
 * If the fractional part of x is exactly 0.5, this function rounds up to the nearest even integer.
 */
real_t myAround( real_t x );

/**
 * @brief Computes the fractional part of a real number.
 * @param x The real number to compute the fractional part of.
 * @return The fractional part of x.
 */
real_t myDecimal( real_t x );

bool isInteger( real_t x );

real_t modulo2PI( real_t theta );

bool isInIntervalMod2PI( real_t x, std::array<real_t, 2> interval);

/**
 * @brief Returns the sign of a value.
 * @tparam T Type of the value.
 * @param[in] val The value.
 * @returns -1 if the value is negative, 1 if it is positive, and 0 if it is zero.
 * @note This function uses the fact that boolean expressions can be implicitly converted to integers
 *       (1 for true, 0 for false).
 */
template <typename T>
int
sgn( T val )
{
    return ( T( 0 ) < val ) - ( val < T( 0 ) );
}

#endif /* __EITMATH_HPP__ */
