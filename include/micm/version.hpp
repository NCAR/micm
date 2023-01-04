#pragma once

#include <string>

namespace micm{

/**
 * @brief Retrive all componenets of this version of MICM
 *
 * @return std::string
 */
std::string getmicmVersion();

/**
 * @brief  Retrieve the major version number
 *
 * @return unsigned
 */
unsigned getmicmVersionMajor();

/**
 * @brief  Retrieve the minor version number
 *
 * @return unsigned
 */
unsigned getmicmVersionMinor();

/**
 * @brief  Retrieve the patch version number
 *
 * @return unsigned
 */
unsigned getmicmVersionPatch();

/**
 * @brief  Retrieve the tweak version number
 *
 * @return unsigned
 */
unsigned getmicmVersionTweak();

}