/**
 * @file Grid.hpp
 * @author 
 * @brief Implementation of the Grid class
 */
 
#ifndef __GRID_HPP__
#define __GRID_HPP__

#include "Form.hpp"
#include "Point.hpp"

/**
 * @class Grid
 * @brief Class Grid creates an object grid that stores all the information needed 
 *        to create the linear system Ax=b*.
 */
class Grid
{
   public:
    /**
     * @brief Default Constructor of a new Grid object.
     * 
     * @param N GridDimension object containing the grid dimensions.
     * @param form Pointer to a Form object.
     */
    Grid( GridDimension N, Form *form );
    /**
     * @brief Copy constructor for Grid objects.
     * 
     * @param that Reference to the Grid object to be copied.
     */
    Grid( const Grid &that );
    /**
     * @brief Assignment operator overload for Grid objects.
     * 
     * @param that Reference to the Grid object to be assigned.
     * @return Grid& Reference to the assigned Grid object.
     */
    Grid &operator=( const Grid &that );
    /**
     * @brief Destroy the Grid object.
     * 
     * Deletes all Point objects in the grid.
     */
    ~Grid();

    /**
     * @brief Creates a new vector of type bivector with the specified sizes.
     *
     * The vector has a size of n_neg = the number of interface points + nb of electrodes
     * (m_interfaces.size() + EIT_ELEC), and n_pos = the number of grid points Nx * Ny
     * (m_points.size()).
     *
     * @tparam _Type The type of elements in the vector.
     * @return A pointer to a newly created EitBiVector object.
    */
    template <class _Type> 
    EitBiVector<_Type> * 
    createVector()
    {
        return new EitBiVector<_Type>( m_interfaces.size() + EIT_ELEC, m_points.size() );
    }

    /**
     * @brief Creates a new vector of type bivector with the specified sizes.
     *
     * The vector has a size of n_neg = the number of interface points + nb of electrodes
     * (m_interfaces.size() + EIT_ELEC), and n_pos = the number of grid points Nx * Ny
     * (m_points.size()).
     *
     * @tparam _Type The type of elements in the vector.
     * @return A new EitBiVector object.
     */
    template <class _Type>
    EitBiVector<_Type>
    createBiVector()
    {
        return EitBiVector<_Type>( m_interfaces.size() + EIT_ELEC, m_points.size() );
    }
    /**
     * @brief Get the grid dimensions (the size of the cartesian grid).
     *  
     * @return const GridDimension& Reference to a const GridDimension object 
     *         containing the grid dimensions.
     */
    const GridDimension &getSizes() const;
    /**
     * @brief Get the grid padding (pads of a cartesian grid).
     * 
     * @return const std::array<real_t, EIT_DIM>& Reference to a 
     *         const std::array<real_t, EIT_DIM> object containing the grid padding.
     */
    const std::array<real_t, EIT_DIM> &getPads() const;
    /**
     * @brief Returns the local beginning index.
     * @return The local beginning index.
     */
    uint_t getLocalBeg() const;
    /**
     * @brief Returns the local end index.
     * @return The local end index.
     */
    uint_t getLocalEnd() const;
    /**
     * @brief Returns the total number of points in the grid (so it returns the size of m_points).
     * @return The total number of points in the grid.
     */
    uint_t getNumberOfPoints() const;
    /**
     * @brief Returns the total number of interfaces in the grid (m_interface belongs to all procs).
     * @return The total number of interfaces in the grid.
     */
    uint_t getNumberOfInterfaces() const;
    /**
     * @brief Returns a reference to the vector containing the interfaces.
     * @return A reference to the vector containing the interfaces.
     */
    const std::vector<Point *> &getInterfaces() const;
    /**
     * @brief Returns a pointer to the point at the specified index. 
     * ith component of a vector of type Point that stores the grid points.
     * @param i The index of the point.
     * @return A pointer to the point at the specified index.
     */
    Point *getPoint( uint_t i );
    /**
     * @brief Returns a pointer to the interface at the specified index.
     * ith component of a vector of type Point that stores the grid interfaces
     * @param i The index of the interface.
     * @return A pointer to the interface at the specified index.
     */
    Point *getInterface( uint_t i );
    /**
     * @brief Returns a reference to the vector containing all points (locally).
     * @return A reference to the vector containing all points.
     */
    const std::vector<Point *> &getPoints() const;
    /**
     * @brief Returns a reference to the vector containing the right halo points.
     * Get the value of the communicated object of type Point from the right.
     * @return A reference to the vector containing the right halo points.
     */
    const std::vector<Point *> &getHaloR() const;
    /**
     * @brief Returns a reference to the vector containing the left halo points.
     * @return A reference to the vector containing the left halo points.
     */
    const std::vector<Point *> &getHaloL() const;
    /**
     * @brief Returns a vector containing the number of points in each process.
     * @return A vector containing the number of points in each process.
     */
    std::vector<uint_t> getSizeAllProc();
    /**
     * @brief Builds m_points locally.
     */
    void buildPoints();
    /**
     * @brief Builds the left and right halo vectors for ecah proc.
     */
    void buildHalos();
    /**
     * @brief Build m_intefaces. 
     *        Constructs the interface points between two different materials using 
     *        the level-set function to identify the zero-level surface.
     */
    void buildInterfaces();
    /**
     * @brief Resets grid points to be regular and corrects connections between
     *        neighboring grid points and interfaces.
     * Modifies the values of the neighbors of some objects Point so we can build 
     *        later the correct linear system.
     * @details The function resets all grid points to be regular by calling the
     *          `setAsRegular()` function for each point. Then, it computes the
     *          minimum and maximum global indices of grid points that are owned by
     *          the current process. Next, it iterates over each interface point
     *          and its neighboring points to correct the connections between
     *          them. If a neighboring point is within the range of global indices
     *          owned by the current process, it is marked as irregular by calling
     *          the `setAsIrregular()` function with the new neighboring information.
     */
    void correctConnections();
    /**
     * @brief Set the Form object associated with this Grid.
     *
     * @param[in] form  Pointer to the Form object to associate with this Grid.
     */
    void setForm( Form *form );
    /**
     * @brief Broadcast to all processors the number of nodes managed by each processor
     *        and store it in an array of size parallel::size.
     */
    void bcastSizes();
    /**
     * @brief Get the Form object associated with this Grid.
     *
     * @return Pointer to the Form object associated with this Grid.
     */
    Form *getForm() const;

   protected:
    GridDimension        m_N;
    uint_t               m_NxBeg, m_NxEnd;
    Vertex               m_h;
    std::vector<Point *> m_points;                  // points :: stores the grid Points
    std::vector<Point *> m_halo_left, m_halo_right; // points on other proc [me-1, me+1]
    std::vector<Point *> m_interfaces; // interfaces points :: stores the interface Points
    Form                *m_form;
    std::vector<uint_t>  m_sizeproc;

    /**
     * @brief Convert a vertex in index space (i,j,k) to coordinate space (x,y,z).
     *
     * @param[in]  ijk   Vertex in index space.
     * @param[out] coor  Vertex in coordinate space.
     */
    void toCoordinate( const Vertex &ijk, Vertex &coor );
    /**
     * @brief Convert a vertex in coordinate space (x,y,z) to index space (i,j,k).
     *
     * @param[in]  coor  Vertex in coordinate space.
     * @param[out] ijk   Vertex in index space.
     */
    void toIJK( const Vertex &coor, Vertex &ijk );
};

#endif /* __GRID_HPP__ */
