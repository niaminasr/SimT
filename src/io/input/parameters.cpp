/* parameters.cpp
 * Parameters class. Input file reader.
 */
#include "parameters.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include "stringTools.hpp"

using std::string;
using std::vector;

// ======================================================================
// Class Parameters
// ======================================================================
// PUBLIC MEMBERS
// ======================================================================

// ======================================================================
// Construct, assign

Parameters::Parameters() {}

Parameters::Parameters( string fileName )
{
    m_fileName = fileName;
    readInputFile();
}

Parameters::Parameters( const Parameters &that )
{
    this->m_fileName = that.m_fileName;
    // We do not use readInputFile as the config
    // may have been changed manually...
    this->m_parameters = that.m_parameters;
}

Parameters &
Parameters::operator=( Parameters &that )
{
    if ( this != &that )
    {
        m_fileName   = that.m_fileName;
        m_parameters = that.m_parameters;
    }
    return *this;
}

// ======================================================================
// Get values from configuration

double
Parameters::getReal( string key )
{
    string KEY = toUpper( stringTrim( key ) );
    if ( !hasKey( KEY ) )
        errNoKey( key );
    else
        return std::stod( m_parameters.at( KEY ) );
    return 0.e0;
}

int_t
Parameters::getInt( string key )
{
    string KEY = toUpper( stringTrim( key ) );
    if ( !hasKey( KEY ) )
        errNoKey( key );
    else
        return (int_t)std::stoi( m_parameters.at( KEY ) );
    return 0;
}

uint_t
Parameters::getUInt( string key )
{
    string KEY = toUpper( stringTrim( key ) );
    if ( !hasKey( KEY ) )
        errNoKey( key );
    else
        return (uint_t)std::stoi( m_parameters.at( KEY ) );
    return 0;
}

string
Parameters::getString( string key )
{
    string KEY = toUpper( stringTrim( key ) );
    if ( !hasKey( KEY ) )
        errNoKey( key );
    else
        return m_parameters.at( KEY );
    return string( "" );
}

double
Parameters::getReal( string key, double defVal )
{
    string KEY = toUpper( stringTrim( key ) );
    return hasKey( KEY ) ? std::stod( m_parameters.at( KEY ) ) : defVal;
}

int
Parameters::getInt( string key, int defVal )
{
    string KEY = toUpper( stringTrim( key ) );
    return hasKey( KEY ) ? std::stoi( m_parameters.at( KEY ) ) : defVal;
}

uint_t
Parameters::getUInt( string key, uint_t defVal )
{
    string KEY = toUpper( stringTrim( key ) );
    return hasKey( KEY ) ? (uint_t)std::stoi( m_parameters.at( KEY ) ) : defVal;
}

string
Parameters::getString( string key, string defVal )
{
    string KEY = toUpper( stringTrim( key ) );
    return hasKey( KEY ) ? m_parameters.at( KEY ) : defVal;
}

// =====================================================================
// Modify configuration

void
Parameters::set( string key, string value )
{
    string KEY        = toUpper( stringTrim( key ) );
    string VALUE      = stringTrim( value );
    m_parameters[KEY] = VALUE;
}

// =====================================================================
// Quick info

bool
Parameters::hasKey( string key )
{
    return m_parameters.find( stringTrim( toUpper( key ) ) ) != m_parameters.end();
}

// =====================================================================
// I/O

void
Parameters::display()
{
    std::cout << "Parameters from input file " << m_fileName << std::endl;
    for ( auto it : m_parameters )
        std::cout << " - " << it.first << " : " << it.second << std::endl;
}

void
Parameters::save( string fileName )
{
    std::ofstream f( fileName.c_str() );
    for ( auto it : m_parameters )
        f << it.first << " = " << it.second << std::endl;
}

// ======================================================================
// PRIVATE MEMBERS
// ======================================================================

void
Parameters::readInputFile( bool clear )
{
    if ( clear )
        m_parameters.clear();

    std::ifstream pfile( m_fileName.c_str() );
    if ( !pfile.is_open() )
        std::cerr << "Couldn't open parameters file: " << m_fileName << std::endl
                  << "           Will use default values if available." << std::endl;

    vector<string> lines;
    string         l, key, value, KEY;
    // Store all the lines of the file
    while ( getline( pfile, l ) )
        lines.push_back( l );

    // Loop on lines
    for ( string &line : lines )
    {
        // The line must have an = sign and must not be a comment
        if ( line.find( '=' ) != line.npos and line.at( 0 ) != '#' )
        {
            std::stringstream ss( line );
            getline( ss, key, '=' );
            getline( ss, value );

            key   = stringTrim( key );
            KEY   = toUpper( key );
            value = stringTrim( value );

            if ( KEY.length() == 0 )
                std::cerr << "Warning: Empty key detected in configuration file" << std::endl;
            if ( value.length() == 0 )
                std::cerr << "Warning: Empty value for key \"" << key << "\"" << std::endl;
            if ( m_parameters.find( KEY ) != m_parameters.end() )
                std::cerr << "Warning: Key \"" << key
                          << "\" found several times in parameters file." << std::endl
                          << "         The value begin read last will be used." << std::endl;

            m_parameters.insert( std::make_pair( KEY, value ) );
        }
    }
}

void
Parameters::errNoKey( string key ) const
{
    std::cerr << "Error : no key \"" << key
              << "\" in configuration file and no default value provided." << std::endl;
    abort();
}

std::ostream &
operator<<( std::ostream &os, const Parameters &p )
{
    os << "\nParameters from input file: \"" << p.m_fileName << "\"" << std::endl;
    for ( auto pair : p.m_parameters )
        // os << " \033[0;31m➤\033[0m " <<pair.first << ": " << pair.second << std::endl;
        os << "😈 " << pair.first << ": " << pair.second << std::endl;
    return os;
}