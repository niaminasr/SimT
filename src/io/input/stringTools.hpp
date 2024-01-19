/**
 * @file stringTools.hpp
 * @brief A collection of string manipulation functions.
 */
 
#ifndef SRC_IO_INPUT_STRINGTOOLS
#define SRC_IO_INPUT_STRINGTOOLS

#include <string>
#include <vector>

 /**
 * @brief Splits a string into substrings based on delimiters.
 * @details All trailing and preceeding spaces are removed from the original string,
 *          prior to any operation.
 * Successive delimiters lead to empty strings being inserted.
 * Example: stringSplit("   foo ,,, bar ,,  ",",") returns the std::vector {"foo","","","bar","",""}.
 * 
 * @param s The input string.
 * @param delims A string containing the delimiters.
 * @return A vector of substrings.
 */
std::vector<std::string> stringSplit( std::string s, std::string delim = " \t" );

/**
 * @brief Trims whitespace from the beginning and end of a string.
 * Removes trailing and preceeding spaces from given string.
 * 
 * @param s The input string.
 * @return The trimmed string.
 */
std::string stringTrim( std::string );

/**
 * @brief Converts a string to uppercase.
 * Converts eligible characters of input std::string to upper case.
 * 
 * @param s The input string.
 * @return The uppercase version of the string.
 */
std::string toUpper( std::string );

#endif /* SRC_IO_INPUT_STRINGTOOLS */
