/**
 * @file EitVector.hpp
 * @brief This file provides the declaration for the EitVector class.
 */

#ifndef SRC_TOOLS_EITVECTOR
#define SRC_TOOLS_EITVECTOR

#include <vector>
#include "ParEITConfig.hpp"

/**
 * @class EitVector
 * @brief A class representing a vector of real numbers.
 */
class EitVector
{
   public:
    using iterator               = typename std::vector<real_t>::iterator;
    using const_iterator         = typename std::vector<real_t>::const_iterator;
    using const_reverse_iterator = typename std::vector<real_t>::const_reverse_iterator;
    using reverse_iterator       = typename std::vector<real_t>::reverse_iterator;

   public:
    /**
     * @brief Default constructor.
     * Constructs an empty vector.
     */
    EitVector();
    /**
     * @brief Constructor with size.
     * Constructs a vector with @p n elements and initializes all elements to 0.
     * @param n The size of the vector.
    */
    EitVector( uint_t n );
    /**
    * @brief Copy constructor.
    * Constructs a vector by copying another vector.
    * @param that The vector to copy.
    */
    EitVector( const EitVector &that );
    /**
    * @brief Destructor.
    * Destroys the vector and releases the memory allocated for the vector.
    */
    ~EitVector();
    /**
     * @brief Copy assignment operator.
     * Assigns the contents of one vector to another vector.
     * @param that The vector to copy.
     * @return A reference to the vector after the assignment.
     */
    EitVector &operator=( const EitVector &that );
    /**
     * @brief Returns the size of the vector.
     * @return The size of the vector.
     */
    uint_t size() const;
    /**
     * @brief Resizes the vector.
     * Changes the size of the vector and initializes new elements to 0.
     * @param n The new size of the vector.
     */
    void resize( uint_t n );
    /**
     * @brief Accesses an element of the vector.
     * @param index The index of the element to be accessed.
     * @return A reference to the accessed element.
     */
    real_t &operator[]( uint_t index );
    /**
     * @brief Accesses a const element of the vector.
     * @param index The index of the element to be accessed.
     * @return A const reference to the accessed element.
     */
    const real_t &operator[]( uint_t index ) const;
    /**
     * @brief Accesses an element of the vector with bounds checking.
     * @param index The index of the element to be accessed.
     * @return A reference to the accessed element.
     */
    real_t &at( uint_t index );
    /**
     * @brief Accesses a const element of the vector with bounds checking.
     * @param index The index of the element to be accessed.
     * @return A const reference to the accessed element.
     */
    const real_t &at( uint_t index ) const;
    /**
     * @brief Sets all elements of the vector to a given value.
     * @param value The value to be set.
     */
    void fill( real_t value );
    /**
     * @brief Calculates the reciprocal vector of the current vector.
     * @return A new EitVector with the reciprocal values of the current vector.
     */
    EitVector reciproqual() const;
    /**
     * @brief Creates a new EitVector with all elements initialized to zero.
     * @return A new EitVector with all elements set to zero.
     */
    EitVector zero() const;
    /**
     * @brief Calculates the dot product of the current vector and the specified vector.
     * @param that The EitVector to calculate the dot product with.
     * @return The dot product of the two vectors.
     */
    real_t dotProduct( const EitVector &that ) const;
    /**
     * @brief Calculates the L2-norm (Euclidean length) of the current vector.
     * @return The L2-norm of the vector.
     */
    real_t norml2() const;
    /**
     * @brief Returns an iterator to the beginning of the vector.
     * @return An iterator pointing to the beginning of the vector.
     */
    iterator begin();
    /**
     * @brief Returns a const reverse iterator to the beginning of the vector.
     * @return A const reverse iterator pointing to the beginning of the vector.
     */
    const_reverse_iterator rbegin() const;
    /**
     * @brief Returns a const iterator to the beginning of the vector.
     * @return A const iterator pointing to the beginning of the vector.
     */
    const_iterator begin() const;
    /**
    * @brief Returns a const iterator to the beginning of the vector.
    * @return A const iterator pointing to the beginning of the vector.
    */
    const_iterator cbegin() const;
    /**
     * @brief Returns a const reverse iterator to the beginning of the vector.
     * @return A const reverse iterator pointing to the beginning of the vector.
     */
    const_reverse_iterator crbegin() const; 
    /**
     * @brief Returns an iterator to the end of the vector.
     * @return An iterator pointing to the end of the vector.
     */
    iterator end();
    /**
     * @brief Returns a const reverse iterator to the end of the vector.
     * @return A const reverse iterator pointing to the end of the vector.
     */
    const_reverse_iterator rend() const; // j'ai modifie
    /**
     * @brief Returns a const iterator to the end of the vector.
     * @return A const iterator pointing to the end of the vector.
     */
    const_iterator end() const;
    /**
     * @brief Returns a const iterator to the end of the vector.
     * @return A const iterator pointing to the end of the vector.
     */
    const_iterator cend() const;
    /**
     * @brief Returns a const reverse iterator to the end of the vector.
     * @return A const reverse iterator pointing to the end of the vector.
     */
    const_reverse_iterator crend() const;
    /**
     * @brief Accessor for the data array of the vector.
     * @return A pointer to the data array.
     */
    real_t *data();
    /**
     * @brief Const accessor for the data array of the vector.
     * @return A const pointer to the data array.
     */
    const real_t *data() const;
    /**
    * @brief Adds the elements of another vector to this vector.
    *
    * @param that The vector to add to this vector.
    * @return A reference to this vector.
    */
    EitVector &operator+=( const EitVector &that );
    /**
     * @brief Adds a scalar value to each element of this vector.
     *
     * @param coeff The scalar value to add.
     * @return A reference to this vector.
     */
    EitVector &operator+=( const real_t &coeff );
    /**
     * @brief Subtracts the elements of another vector from this vector.
     *
     * @param that The vector to subtract from this vector.
     * @return A reference to this vector.
     */
    EitVector &operator-=( const EitVector &that );
    /**
     * @brief Subtracts a scalar value from each element of this vector.
     *
     * @param coeff The scalar value to subtract.
     * @return A reference to this vector.
     */
    EitVector& operator-=(const real_t& coeff);
    /**
     * @brief Multiplies the elements of this vector by the corresponding elements of another vector.
     *
     * @param that The vector to multiply with this vector.
     * @return A reference to this vector.
     */
    EitVector& operator*=(const EitVector& that);
    /**
     * @brief Multiplies each element of this vector by a scalar value.
     *
     * @param coeff The scalar value to multiply.
     * @return A reference to this vector.
     */
    EitVector& operator*=(const real_t& coeff);
    /**
     * @brief Divides the elements of this vector by the corresponding elements 
     *        of another vector, ignoring elements of the other vector with an 
     *        absolute value less than a threshold.
     *
     * @param that The vector to divide this vector by.
     * @return A reference to this vector.
     */
    EitVector& operator/=(const EitVector& that);
    /**
     * @brief Divides each element of this vector by a scalar value, ignoring the operation if the scalar is too small.
     *
     * @param coeff The scalar value to divide by.
     * @return A reference to this vector.
     */
    EitVector& operator/=(const real_t& coeff);

   protected:
    std::vector<real_t> *m_data = nullptr;
};

/**
 * @brief Function to check the sizes of two vectors and raise an error if they are different
 * @param func Name of the calling function
 * @param n Size of first vector
 * @param m Size of second vector
 */
void checkSizes( std::string func, uint_t n, uint_t m );
/**
 * @brief Addition operator overload for two EitVectors.
 * @param a The first EitVector.
 * @param b The second EitVector.
 * @return The result of adding the two EitVectors.
 */
EitVector operator+( const EitVector &a, const EitVector &b );
/**
 * @brief Addition operator overload for an EitVector and a real number.
 * @param a The EitVector.
 * @param coeff The real number.
 * @return The result of adding the real number to the EitVector.
 */
EitVector operator+( const EitVector &a, const real_t &coeff );
/**
 * @brief Addition operator overload for a real number and an EitVector.
 * @param coeff The real number.
 * @param a The EitVector.
 * @return The result of adding the real number to the EitVector.
 */
EitVector operator+( const real_t &coeff, const EitVector &a );
/**
 * @brief Subtraction operator overload for two EitVectors.
 * @param a The first EitVector.
 * @param b The second EitVector.
 * @return The result of subtracting the second EitVector from the first EitVector.
 */
EitVector operator-( const EitVector &a, const EitVector &b );
/**
 * @brief Subtraction operator overload for an EitVector and a real number.
 * @param a The EitVector.
 * @param coeff The real number.
 * @return The result of subtracting the real number from the EitVector.
 */
EitVector operator-( const EitVector &a, const real_t &coeff );
/**
 * @brief Subtraction operator overload for a real number and an EitVector.
 * @param coeff The real number.
 * @param a The EitVector.
 * @return The result of subtracting the EitVector from the real number.
 */
EitVector operator-( const real_t &coeff, const EitVector &a );
/**
 * @brief Overloaded operator * for element-wise multiplication of two EitVector objects
 * @param a EitVector object to be multiplied
 * @param b EitVector object to be multiplied
 * @return EitVector object containing the result of the multiplication
 */
EitVector operator*( const EitVector &a, const EitVector &b );
/**
 * @brief Overloaded operator * for multiplying an EitVector object by a scalar value
 * @param a EitVector object to be multiplied
 * @param coeff Scalar value to be multiplied
 * @return EitVector object containing the result of the multiplication
 */
EitVector operator*( const EitVector &a, const real_t &coeff );
/**
 * @brief Overloaded operator * for multiplying a scalar value by an EitVector object
 * @param coeff Scalar value to be multiplied
 * @param a EitVector object to be multiplied
 * @return EitVector object containing the result of the multiplication
 */
EitVector operator*( const real_t &coeff, const EitVector &a );
/**
 * @brief Overloads the division operator "/" for two EitVector objects.
 * Returns a new EitVector object which is the result of dividing
 * each component of vector 'a' by the corresponding component in vector 'b'.
 * @param a The first EitVector object to be divided.
 * @param b The second EitVector object to be divided.
 * @return A new EitVector object that is the result of the division operation.
 */
EitVector operator/( const EitVector &a, const EitVector &b );
/**
 * @brief Overloads the division operator "/" for an EitVector object and a scalar value.
 * Returns a new EitVector object which is the result of dividing
 * each component of vector 'a' by the scalar value 'coeff'.
 * @param a The EitVector object to be divided.
 * @param coeff The scalar value to divide each component of the vector by.
 * @return A new EitVector object that is the result of the division operation.
 */
EitVector operator/( const EitVector &a, const real_t &coeff );
/**
 * @brief Overloads the division operator "/" for a scalar value and an EitVector object.
 * Returns a new EitVector object which is the result of dividing
 * each component of vector 'a' by the scalar value 'coeff'.
 * @param coeff The scalar value to divide each component of the vector by.
 * @param a The EitVector object to be divided.
 * @return A new EitVector object that is the result of the division operation.
 */
EitVector operator/( const real_t &coeff, const EitVector &a );

#endif /* SRC_TOOLS_EITVECTOR */
