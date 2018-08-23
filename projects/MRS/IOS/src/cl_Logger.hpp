#ifndef MORIS_IOS_CL_LOGGER_HPP_
#define MORIS_IOS_CL_LOGGER_HPP_
#include <iostream>


/**
 * @brief Macro for function scope markup.
 * The scope name is constructed with help of compiler
 * and contains current function name.
 *
 * @include "IOS/src/cl_Logger/log_function.inc"
 */
#define MORIS_LOG_FUNCTION() std::cout

/**
 * Log an informational severity message.
 *
 * @include "IOS/src/cl_Logger/log.inc"
 */
#define MORIS_LOG            std::cout

/**
 * Log a trace severity message.
 *
 * @include "IOS/src/cl_Logger/log_trace.inc"
 */
#define MORIS_LOG_TRACE      std::cerr

/**
 * Log a debug severity message.
 *
 * @include "IOS/src/cl_Logger/log_debug.inc"
 */
#define MORIS_LOG_DEBUG      std::cerr

/**
 * Log an informational sevezrity message.
 *
 * @include "IOS/src/cl_Logger/log_info.inc"
 */
//#define MORIS_LOG_INFO       BOOST_LOG_SEV(moris_logger::get(), trivial::info)

#define MORIS_LOG_INFO       std::cout

/**
 * Log a warning severity message.
 *
 * @include "IOS/src/cl_Logger/log_warning.inc"
 */
#define MORIS_LOG_WARNING    std::cerr

/**
 * Log an error severity message.
 *
 * @include "IOS/src/cl_Logger/log_error.inc"
 */
#define MORIS_LOG_ERROR      std::cerr

/**
 * Log a fatal severity message.
 *
 * @include "IOS/src/cl_Logger/log_fatal.inc"
 */
#define MORIS_LOG_FATAL      std::cerr


#endif	/* MORIS_IOS_CL_LOGGER_HPP_ */
