/**
* @file System.hpp
* @brief Header file for the System class, which represents a system of equations in a 
*        finite element method context.
*/

#ifndef __SYSTEM_HPP_
#define __SYSTEM_HPP_

#include "ParEITConfig.hpp"
#include "geometry/Grid.hpp"
#include "tools/EitBiVector.hpp"
#include <spm.h>

using namespace std;


/**
 * @class System
 * @brief Represents a system of equations in a IBM finite difference method context.
 */
class System
{
   private:
   public:
    vector<real_t> *b; /* The right-hand side vector.*/
    vector<real_t> *coefs; /* Coefficients vector. */
    vector<int64_t> *x_coords; /* The x coordinates vector. */
    vector<int64_t> *y_coords; /* The y coordinates vector. */
    uint_t          ncol; /* The number of columns in the system matrix. */
    /* not updated until system is built, if constructed without grid. */
    uint_t          nrow; /* The number of rows in the system matrix. */
    /* The number of non-zero elements in the system matrix. */
    uint_t          nnz;

    /**
     * @brief Default constructor.
     * Creates an empty system with zero size.
     */
    System()
    {
        b        = new vector<real_t>();
        coefs    = new vector<real_t>();
        x_coords = new vector<int64_t>();
        y_coords = new vector<int64_t>();
        ncol     = 0;
        nrow     = 0;
        nnz      = 0;
    }
    /**
     * @brief Constructor that creates a system with the same size as the grid.
     * @param grid The grid object.
     */
    System( const Grid &grid )
    {
        b        = new vector<real_t>();
        coefs    = new vector<real_t>();
        x_coords = new vector<int64_t>();
        y_coords = new vector<int64_t>();
        ncol     = grid.getNumberOfInterfaces() + grid.getNumberOfPoints() + EIT_ELEC;
        nrow     = grid.getNumberOfInterfaces() + grid.getNumberOfPoints() + EIT_ELEC;
    }
    /**
     * @brief Destructor that frees the memory used by the system..
     */
    ~System()
    {
        delete b;
        delete coefs;
        delete x_coords;
        delete y_coords;
    }
    /**
     * @brief Assembles the full system of equations.
     * This function iterates over all interfaces, electrodes, and regular points 
     * in the grid to assemble the full system of equations.
     * @param grid The grid object.
     * @param sigma The conductivity vector.
     * @param u The voltage vector.
     */
    void buildFullSystem( const Grid &grid, const EitBiVector<real_t> &sigma,
                          const EitBiVector<real_t> &u );
    /**
     * @brief Assembles the system matrix for an interface point.
     *
     * @param grid The grid object.
     * @param sigma The conductivity vector.
     * @param u The voltage vector.
     * @param interfaceIndex The index of the interface.
     */
    void assembleInterfacePoint( const Grid &grid, const EitBiVector<real_t> &sigma,
                                 const EitBiVector<real_t> &u, uint_t i );
    /**
     * @brief Assembles the system matrix for an electrode point.
     *
     * @param grid The grid object.
     * @param sigma The conductivity vector.
     * @param u The voltage vector.
     * @param electrodeIndex The index of the electrode.
     */
    void assembleElectrodePoint( const Grid &grid, const EitBiVector<real_t> &sigma,
                                 const EitBiVector<real_t> &u, uint_t i );
    /**
     * @brief Assembles the system matrix for a regular point.
     *
     * @param grid The grid object.
     * @param sigma The conductivity vector.
     * @param u The voltage vector.
     * @param pointIndex The index of the regular point.
     * @return 0 if successful, -1 otherwise.
     */
    uint_t assembleRegularPoint( const Grid &grid, const EitBiVector<real_t> &sigma,
                                 const EitBiVector<real_t> &u, uint_t i );
    /**
     * @brief Writes the system matrix to a Matrix Market file.
     * The Matrix Market file format is a standard format for sparse matrices.
     * @param filename The name of the file to write to.
     */
    void writeMatrixMarket();

    /**
     * @brief builds an spm structure
     *
     */
    spmatrix_t* buildSpm();
};

#endif /* __SYSTEM_HPP_ */