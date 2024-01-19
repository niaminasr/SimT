#include "tools/EitMath.hpp"

real_t
myAround( real_t x )
{
    real_t c = ceil( x );
    real_t f = floor( x );

    if ( std::abs( x - c ) < EPSILON )
        return c;
    else if ( std::abs( x - f ) < EPSILON )
        return f;
    return x;
}

real_t
myDecimal( real_t x )
{
    real_t f = floor( x );
    return x - f;
}

bool isInteger( real_t x ){
    return (( std::abs( x - ceil(x) ) < EPSILON ) || (std::abs( x - floor(x) ) < EPSILON ));
}

real_t
modulo2PI( real_t theta )
{
    real_t factor;
    if(theta > 0.0) factor = floor(theta / (2*EIT_PI));
    else            factor = ceil(theta / (2*EIT_PI));
    real_t r = theta - (factor*2*EIT_PI);
    return r;
}

bool isInIntervalMod2PI( real_t x, std::array<real_t, 2> interval ){
    interval = {modulo2PI(interval[0]), modulo2PI(interval[1])};
    
    if ((std::abs(modulo2PI(x) - EIT_PI) < EPSILON) || (std::abs(modulo2PI(x) + EIT_PI) < EPSILON))
        return ((interval[0] <= -x && -x <= interval[1]) || (interval[0] <= x && x <= interval[1]));  
    return (interval[0] <= x && x <= interval[1]);
}