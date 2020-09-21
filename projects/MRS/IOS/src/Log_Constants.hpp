#ifndef MORIS_LOG_CONSTANTS_HPP_
#define MORIS_LOG_CONSTANTS_HPP_

#include "core.hpp"

/**
 * Macro constants used for formatting the output by
 * the logger and the queries
 *
 */

const moris::sint LOGGER_DEFAULT_SEVERITY_LEVEL = 3;
const moris::uint LOGGER_DEFAULT_DIRECT_OUTPUT = 1;
const moris::uint LOGGER_FLOAT_PRECISION = 6;

const std::string LOGGER_HEADER_BEGIN = ">>>>>";
const std::string LOGGER_HEADER_END = "<<<<<";

const moris::uint QUERY_MAX_COLUMN_WIDTH = 18;

#endif /* MORIS_LOG_CONSTANTS_HPP_ */
