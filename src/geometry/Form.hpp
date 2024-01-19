/**
 * @file Form.hpp
 * @author 
 * @brief Implementation of the Form class
 */

#ifndef __FORM_HPP__
#define __FORM_HPP__

#include "geometry/Point.hpp"
#include "tools/EitBiVector.hpp"
#include "tools/EitMath.hpp"
#include "tools/EitVector.hpp"

using FormCoeffs = std::vector<real_t>;

/**
 * @class Form
 * @brief Creates an object form that represents the geometry in wich we resolve the EIT problem.
 * Represents the form of an object using Fourier descriptors.
 */
class Form
{
   public:
    /**
     * @brief Default constructor.
     */
    Form();
     /**
     * @brief Constructor that creates a Form object with a given number of Fourier descriptors.
     *
     * @param n Number of Fourier descriptors.
     * @param ne Number of electrodes.
     */
    Form( uint_t n, uint ne );
    /**
     * @brief Copy constructor.
     *
     * @param that Another Form object.
     */
    Form( const Form &that );
    /**
     * @brief Returns the ith component of the shape vector m_alpha(shape of the domain).
     *
     * @param i Index of the Fourier descriptor coefficient.
     * @return Reference to the i-th Fourier descriptor coefficient.
     */
    real_t &coeff( uint_t i );
    /**
     * @brief Returns the ith component of the shape vector m_thetas(positions of the electrodes).
     *
     * @param i Index of the electrode.
     * @return Reference to the first angle of the i-th electrode.
     */
    real_t &begAngle( uint_t i );
    /**
     * @brief Returns the shape vector m_thetas.
     *
     * @return Vector of ElectrodeThetas.
     */
    std::vector<ElectrodeThetas> &getThetas();
    /**
     * @brief Returns the size of the shape vector alpha.
     *
     * @return Number of Fourier descriptors.
     */
    uint_t size() const;
    /**
     * @brief Returns the radius of the form at each angle theta.
     *
     * @param theta Angle in radians.
     * @return Radius of the object at the given angle.
     */
    real_t radius( real_t theta );
    /**
     * @brief Returns the first derivative of the radius.
     *
     * @param theta Angle in radians.
     * @return First derivative of the radius of the object at the given angle.
     */
    real_t firstDerivativeofRadius( real_t theta );
    /**
     * @brief Returns the segond derivative of the radius.
     *
     * @param theta Angle in radians.
     * @return Second derivative of the radius of the object at the given angle.
     */
    real_t secondDerivativeofRadius( real_t theta );
    /**
     * @brief Returns  Rho=sqrt(radius*2 +radiusdarivative*2) at angle theta.
     *
     * @param theta Angle in radians.
     * @return Rho value of the object at the given angle.
     */
    real_t rho( real_t theta );
    /**
     * @brief Computes the end of the elctrode, knowing it's length "L" and it's begening angle "tm".
     *
     * @param L Length of the electrode.
     * @param i Index of the electrode.
     * @return End angle of the electrode given a desired length.
     */
    real_t getEndAngleForLength( real_t L, uint_t i );
    /**
     * @brief Computes the level set function of a given point with respect to the object.
     *
     * @param pt Point.
     * @return Level set function of the point with respect to the object.
     */
    real_t levelSet( const Point &pt );
    /**
     * @brief Parathesis operator that returnsthe value of the Levelset at Point pt. 
     * Evaluates the Form object at a given point.
     *
     * @param pt the point to evaluate the Form at.
     * @return the value of the level set function at the given point.
     */
    real_t operator()( const Point &pt );
    /**
     * @brief Computes the first derivative of the level set function at a given point.
     *
     * @param pt The point at which to compute the derivative.
     * @param dim The dimension for which to compute the derivative.
     * @param step The step size to use for the finite difference approximation.
     * @return The first derivative of the level set function in the specified 
     *         dimension (in a certain direction dir) at the given point.
     */
    real_t derivativeLevelSet( const Point &pt, real_t dir, real_t step );
    /**
     * @brief Sets the electrodes plaaced on the Form. 
     * Builds the theta angles for each electrode given a desired length.
     *
     * @param L Length of the electrode.
     */
    void buildThetas( real_t L );
    /**
     * @brief Computes the point on the interface of a circle defined by a point and 
     *        direction using a dichotomie method.
     *
     * @param a The initial point on the circle.
     * @param dim The dimension along which to change the point.
     * @param dir The direction in which to change the point.
     * @return The point on the interface of the circle.
     */
    Point computeInterfacePoint( const Point &a, uint_t dim, real_t dir );
    /**
     * @brief Approximates the integral of the level set function at a given point 
     *        using the Eikonal equation.
     * Computes the value of the delta function of Smereka
     *
     * @param pt The point at which to approximate the integral.
     * @param step The step size to use for the finite difference approximation.
     * @return The approximate integral of the level set function at the given point.
     */
    real_t approximateIntegral( const Point &pt, real_t &step );
    /**
     * @brief Overloads the stream insertion operator to print the contents of a Form object 
     *        (print operator).
     *
     * @param os the output stream to write to.
     * @param f the Form object to print.
     * @return the output stream.
     */
    friend std::ostream &operator<<( std::ostream &os, const Form &f );

   private:
    Point     *m_center; // The center Point
    FormCoeffs m_alphas; // The values of the vector that determines shape are stored here
    std::vector<ElectrodeThetas> m_thetas; // The angles of the ectrodes are stored here
};

/**
 * @brief Overloads the stream insertion operator to print 
 *        the contents of a Form object (print operator).
 *
 * @param os the output stream to write to.
 * @param f the Form object to print.
 * @return the output stream.
 */
std::ostream &operator<<( std::ostream &os, const Form &f );

#endif /* __FORM_HPP__ */
