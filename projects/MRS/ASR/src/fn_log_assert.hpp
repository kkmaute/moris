#ifndef MORIS_ASSERT_FN_LOG_ASSERT_HPP_
#define MORIS_ASSERT_FN_LOG_ASSERT_HPP_

// C++ header files.
#include <sstream>
#include <iostream>
#include <stdexcept>

// MORIS header files.
#include "core.hpp"
#include "cl_Logger.hpp"
#ifdef MORIS_HAVE_PARALLEL
#include <mpi.h>
#endif

// ----------------------------------------------------------------------------

namespace moris
{
namespace assert
{
    /**
     * Overloaded moris::assert::error.
     *
     * @param[in] msg Error message.
     */
    template<typename Exception = std::runtime_error>
    void
    error(
            std::string const & msg)
    {
        MORIS_LOG_ERROR << "*** Error: " << msg;

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
    void
    error(
            std::string const & location,
            std::string const & task,
            std::string const & check,
            Exception   const & exception = Exception())
    {
        MORIS_LOG_ERROR << "*** ---------------------------------------------------------------------------\n";
        MORIS_LOG_ERROR << "*** \n";
        MORIS_LOG_ERROR << "*** Moris encountered an error. If you are not able to resolve this issue\n";
        MORIS_LOG_ERROR << "*** using the information listed below, you can ask for help at\n";
        MORIS_LOG_ERROR << "***\n";
        MORIS_LOG_ERROR << "***     kurt.maute@colorado.edu\n";
        MORIS_LOG_ERROR << "***\n";
        MORIS_LOG_ERROR << "*** Remember to include the error message listed below and, if possible,\n";
        MORIS_LOG_ERROR << "*** include a *minimal* running example to reproduce the error.\n";
        MORIS_LOG_ERROR << "***\n";
        MORIS_LOG_ERROR << "*** ---------------------------------------------------------------------------\n";
        MORIS_LOG_ERROR << "***\n";
        MORIS_LOG_ERROR << "*** " << "Error:   Unable to " << task << ".\n";
        MORIS_LOG_ERROR << "*** " << "Reason:  " << check<<"\n";
        std::istringstream exception_msg(exception.what());
        std::string exception_line;
        while (std::getline(exception_msg, exception_line))
        {
            MORIS_LOG_ERROR << "***          " << exception_line<<"\n";
        }
        MORIS_LOG_ERROR << "*** " << "Where:   This error was encountered inside " << location << ".\n";
#ifdef MORIS_HAVE_PARALLEL
        int mpi_rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
        MORIS_LOG_ERROR << "*** " << "Process: " << mpi_rank<<"\n";
#endif
        MORIS_LOG_ERROR << "*** " << "Version: " << "1.0\n";
        MORIS_LOG_ERROR << "***\n";
        MORIS_LOG_ERROR << "*** ---------------------------------------------------------------------------\n";

        throw exception;
    }

    /**
     * Overloaded moris::assert::moris_assert.
     *
     * @param[in] file      File where assertion was raised.
     * @param[in] line      Line where assertion was raised.
     * @param[in] function  Function where assertion was raised.
     * @param[in] check     Check that raised assertion.
     * @param[in] exception Exception raised by check.
     */
    template<typename Exception>
    void
    moris_assert(
            std::string   const & file,
            moris::size_t const & line,
            std::string   const & function,
            std::string   const & check,
            Exception     const & exception)
    {
        std::stringstream location;
        location << file << " (line " << line << ")";

        std::stringstream task;
        task << "complete call to function " << function << "()";

        std::stringstream reason;
        reason << "Assertion " << check << " failed.";

        moris::assert::error(location.str(), task.str(), reason.str(), exception);
    }

    /**
     * Overloaded moris::assert::moris_assert.
     *
     * @param[in] file      File where assertion was raised.
     * @param[in] line      Line where assertion was raised.
     * @param[in] function  Function where assertion was raised.
     * @param[in] check     Check that raised assertion.
     * @param[in] msg       Error message to build exception.
     */
    inline
    void
    moris_assert(
            std::string   const & file,
            moris::size_t const & line,
            std::string   const & function,
            std::string   const & check,
            char          const * msg)
    {
        moris::assert::moris_assert(file, line, function, check, std::runtime_error(msg));
    }

}    // namespace assert
}    // namespace moris

#endif /* MORIS_ASSERT_FN_LOG_ASSERT_HPP_ */
