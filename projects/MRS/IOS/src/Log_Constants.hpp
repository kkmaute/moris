#ifndef MORIS_LOG_CONSTANTS_HPP_
#define MORIS_LOG_CONSTANTS_HPP_

#include "core.hpp"

/**
 * Macro constants used for formatting the output by
 * the logger and the queries
 *
 */

const moris::sint LOGGER_DEFAULT_SEVERITY_LEVEL = 2;
const moris::uint LOGGER_DEFAULT_DIRECT_OUTPUT  = 3;
const moris::uint LOGGER_FLOAT_PRECISION        = 14;
const moris::uint LOGGER_MAX_NUMEL_MATRIX_PRINT = 32;

const std::string LOGGER_NON_SPECIFIC_ENTITY_TYPE = "NoType";
const std::string LOGGER_ARBITRARY_DESCRIPTOR     = "Arbitrary";
const std::string LOGGER_HEADER_BEGIN             = ">>>>>";
const std::string LOGGER_HEADER_END               = "<<<<<";

const moris::uint QUERY_MAX_COLUMN_WIDTH = 18;

const bool PRINT_WALL_TIME = true;

#endif /* MORIS_LOG_CONSTANTS_HPP_ */
