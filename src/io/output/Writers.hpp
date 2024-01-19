/**
* @file Writers.hpp
* @brief Declaration of output writers for various formats
*/

#ifndef __GNUPLOT_HPP__
#define __GNUPLOT_HPP__

#include "geometry/Grid.hpp"

/**
 * @brief Writes a GNUplot data file with the grid coordinates and the vector values.
 * @param filename Name of the file to be created.
 * @param grid The grid object from which the coordinates will be extracted.
 * @param vec The vector object containing the values to be written.
 * @param onlyItf Flag indicating whether only the grid interface points should be written.
 */
void writeGnuplot( std::string filename, const Grid &grid, const EitBiVector<real_t> &vec,
                   bool onlyItf );

#endif /* __GNUPLOT_HPP__ */
