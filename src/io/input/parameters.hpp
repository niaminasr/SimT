/**
 * @file parameters.hpp
 * @brief Header file for the Parameters class. Input file reader.
 */

#ifndef SRC_IO_INPUT_PARAMETERS
#define SRC_IO_INPUT_PARAMETERS

#include <iostream>
#include <map>
#include <vector>
#include "ParEITConfig.hpp"

 /** 
 * @class Parameters
 * @brief Reads parameters from input file and provides access to them.
 * Manages options for simulations.
 * @details Options are usually read from an input file, but can be set manually.
 * The input line format is the following:
 *
 * keyword = a list of options
 *
 * If the line ends with a backslash character (\),
 * options can be continued on the next line
 *
 */
class Parameters
{
   public:
    /**
     * @brief Default constructor.
     */
    Parameters();
    /**
     * @brief Constructor with a file name that reads the input file.
     * @param fileName Name of the input file, it can be absolute or relative.
     */
    Parameters( std::string fileName );
    /**
     * @brief Copy constructor.
     * @param cfg Parameters to copy from.
     * @details All options are copied from the configuration state,
     *          not the input file (in case some entries were added/modified).
     */
    Parameters( const Parameters &cfg ); // cfg=configuration.
    /**
     * @brief Assignment operator.
     * @param that Parameters to copy from.
     * @details All options are copied from the configuration state,
     *          not the input file (in case some entries were added/modified).
     * @return Reference to the Parameters object.
     */
    Parameters &operator=( Parameters &cfg );

     /**
     * @brief Destructor. Nothing special to do.
     */
    ~Parameters() {};

    // ===============================================================
    // Getters (Get values from config)
    // ===============================================================

    /**
     * @brief Get the value of a real parameter.
     * @param key Key of the parameter.
     * @return Value of the parameter associated to key. If key is missing from file, program stops.
     */
    double getReal( std::string key );
    /**
     * @brief Get the value of an integer parameter.
     * @param key Key of the parameter.
     * @return Value of the parameter.
     */
    int_t getInt( std::string key );
    /**
     * @brief Get the value of an unsigned integer parameter.
     * @param key Key of the parameter.
     * @return Value of the parameter.
     */
    uint_t getUInt( std::string key );
    /**
     *  @brief Get the value of a string parameter.
     * @param key Key of the parameter.
     * @return Value of the parameter.
     */
    std::string getString( std::string key );
   /**
     * @brief Get the value of a real parameter with default value.
     * @param key Key of the parameter.
     * @param defVal Default value.
     * @return Value of the parameter associated to key if it exists, default value otherwise.
     */
    double getReal( std::string key, double defVal );
    /**
     * @brief Get the value of an integer parameter with default value.
     * @param key Key of the parameter.
     * @param defVal Default value.
     * @return Value of the parameter if it exists, default value otherwise.
     */
    int_t getInt( std::string key, int_t defVal );
    /**
     * @brief Get the value of an unsigned integer parameter with default value.
     * @param key Key of the parameter.
     * @param defVal Default value.
     * @return Value of the parameter if it exists, default value otherwise.
     */
    uint_t getUInt( std::string key, uint_t defVal );
    /**
     * @brief Get the value of a string parameter with default value.
     * @param key Key of the parameter.
     * @param defVal Default value.
     * @return Value of the parameter if it exists, default value otherwise.
     */
    std::string getString( std::string key, std::string defVal );

    // ===============================================================
    // Modify configuration
    // ===============================================================

    /**
     * @brief Set the value of a parameter: Replaces value in key, or creates entry if not existing.
     * @param key Key of the parameter.
     * @param value Value of the parameter.
     */
    void set( std::string key, std::string value );

    // ===============================================================
    // Quick information
    // ===============================================================

    /**
     * @brief Check if a parameter exists in the Parameters object.
     * @param key The key to check for.
     * @return True if the key exists, false otherwise.
     */
    bool hasKey( std::string key );

    // ===============================================================
    // I/O
    // ===============================================================

    /**
    * @brief Save all parameters in the Parameters object to a file.
    * Write configuration in given file.
    * @param fileName The name of the file to save the parameters to.
    */
    void save( std::string fileName );
    /**
     * @brief Display all parameters in the Parameters object to stdout.
     * Display configuration on standard output.
     */
    void display();
    /**
     * @brief Overloaded operator to display all parameters to an output stream
     * @param os The output stream to write to
     * @param p The Parameters object to display
     * @return The output stream
     */
    friend std::ostream &operator<<( std::ostream &os, const Parameters &p );

   private:
     /**
     * @brief Read a parameter file and store the parameters in the Parameters object.
     * @param clear If true, clear any previously-stored parameters before reading.
     * @details  Parses input file and store configuration.
     * clear: If true(default), previously defined options are deleted.
     * When false, can be used to complete a configuration with another file.
     * Previously existing entries are overwritten though.
     *
     * The input line format is the following:
     *
     * keyword = a list of options*.
     *
     * If the line ends with a backslash character (\\), options can be continued on the next line.
     */
    void readInputFile( bool clear = true );
    /**
     * @brief Print an error message for a missing key.
     * @param key The missing key.
     */
    void errNoKey( std::string key ) const;

    std::string                        m_fileName;   // input file
    std::map<std::string, std::string> m_parameters; // configuration info
};

/**
 * @brief Overloaded operator to display all parameters to an output stream
 * @param os The output stream to write to
 * @param p The Parameters object to display
 * @return The output stream
 */
std::ostream &operator<<( std::ostream &os, const Parameters &p );

#endif /* SRC_IO_INPUT_PARAMETERS */
