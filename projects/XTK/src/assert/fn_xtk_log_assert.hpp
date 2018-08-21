/*
 * fn_log_assert.hpp
 *
 *  Created on: Jun 19, 2017
 *      Author: ktdoble
 */

#ifndef SRC_ASSERT_FN_LOG_ASSERT_HPP_
#define SRC_ASSERT_FN_LOG_ASSERT_HPP_

// C++ header files.
#include <sstream>
#include <iostream>
#include <stdexcept>

// MPI Header
#include "tools/cl_MPI_Tools.hpp"

#include "core/xtk_typedefs.hpp"
#include "ios/cl_Logger.hpp"

// ----------------------------------------------------------------------------

namespace xtk
{

namespace assert
{
  // Forward declarion of get comm
  MPI_Comm get_comm();

/**
 * Overloaded xtk::assert::error.
 *
 * @param[in] msg Error message.
 */
template<typename Exception = std::runtime_error>
void error(std::string const & msg)
{
    XTK_ERROR << "*** Error: " << msg;

    throw Exception(msg.c_str());
}

/**
 * Log error message with global logger.
 *
 * @param[in] location  Location where error was raised.
 * @param[in] task      Task where error was raised.
 * @param[in] check     Check that raised assertion.
 * @param[in] exception Exception raised by check.
 */
template<typename Exception>
void error(std::string const & location, std::string const & task, std::string const & check, Exception const & exception =
                   Exception())
{
    std::istringstream exception_msg(exception.what());
    std::string exception_line;

#ifdef XTK_HAVE_PARALLEL
    int tProcRank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &tProcRank);
    while (std::getline(exception_msg, exception_line))
    {
        XTK_ERROR<<"\n"
        << "*** ---------------------------------------------------------------------------\n"
        << "***\n"
        << "*** " << "Error:   Unable to " << task << ".\n"
        << "*** " << "Reason:  " << check<<"\n"
        << "***          " << exception_line <<"\n"
        << "*** " << "Where:   This error was encountered inside " << location << ".\n"
        << "*** " << "Version: " << "1.0\n"
        << "***\n"
        << "*** " << "Process: " << tProcRank<<"\n"
        << "*** ---------------------------------------------------------------------------\n";
    }
#else
    while (std::getline(exception_msg, exception_line))
    {
        XTK_ERROR<<"\n"
        << "*** ---------------------------------------------------------------------------\n"
        << "***\n"
        << "*** " << "Error:   Unable to " << task << ".\n"
        << "*** " << "Reason:  " << check<<"\n"
        << "***          " << exception_line <<"\n"
        << "*** " << "Where:   This error was encountered inside " << location << ".\n"
        << "*** " << "Version: " << "1.0\n"
        << "***\n"
        << "*** ---------------------------------------------------------------------------\n";
    }
#endif
    throw exception;
}

/**
 * Overloaded xtk::assert::xtk_assert.
 *
 * @param[in] file      File where assertion was raised.
 * @param[in] line      Line where assertion was raised.
 * @param[in] function  Function where assertion was raised.
 * @param[in] check     Check that raised assertion.
 * @param[in] exception Exception raised by check.
 */
template<typename Exception>
void xtk_assert(std::string const & file,
                xtk::size_t const & line,
                std::string const & function,
                std::string const & check,
                Exception const & exception)
{
    std::stringstream location;
    location << file << " (line " << line << ")";

    std::stringstream task;
    task << "complete call to function " << function << "()";

    std::stringstream reason;
    reason << "Assertion " << check << " failed.";

    xtk::assert::error(location.str(), task.str(), reason.str(), exception);
}

/**
 * Overloaded xtk::assert::xtk_assert.
 *
 * @param[in] file      File where assertion was raised.
 * @param[in] line      Line where assertion was raised.
 * @param[in] function  Function where assertion was raised.
 * @param[in] check     Check that raised assertion.
 * @param[in] msg       Error message to build exception.
 */
inline
void xtk_assert(std::string const & file,
                xtk::size_t const & line,
                std::string const & function,
                std::string const & check,
                char const * msg)
{
    xtk::assert::xtk_assert(file, line, function, check, std::runtime_error(msg));
}

}    // namespace assert
}    // namespace xtk
#endif
