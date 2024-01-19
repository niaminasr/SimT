/**
 * @file EitElectrode.hpp
 * @brief Contains functions related to Electrodes for EIT (Electrical Impedance Tomography)
*/ 

#ifndef SRC_TOOLS_EITELECTRODE
#define SRC_TOOLS_EITELECTRODE

#include <array>
#include <vector>
#include "ParEITConfig.hpp"

using ElectrodeThetas = std::array<real_t, 2>;

/**
 * @brief Finds on which electrode we are located, knowing the angle "theta", the length L and the angles
 *        at which each electrode begin.
 * @param theta The angle at which we want to find the corresponding electrode.
 * @param thetas A vector of ElectrodeThetas objects containing the beginning and ending angles 
                 of each electrode.
 *@return The index of the electrode corresponding to the given angle, or -1 if no electrode is found.
*/
int_t getElectrodeIdForAngle( real_t theta, std::vector<ElectrodeThetas> &thetas );

// revoir ici
#endif /* SRC_TOOLS_EITELECTRODE */
