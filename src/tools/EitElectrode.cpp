#include "tools/EitElectrode.hpp"

// Finds on which electrode we are located, knowing the angle "theta", the length L and the angles
// at which each elctrode begin
int_t
getElectrodeIdForAngle( real_t theta, std::vector<ElectrodeThetas> &thetas )
{
    uint_t n = thetas.size();
    for ( uint_t i = 0; i < n; ++i )
        if ( ( thetas[i][0] <= theta ) and ( theta < thetas[i][1] ) )
            return i;
    return -1;
}
