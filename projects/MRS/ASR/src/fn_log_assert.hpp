/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_log_assert.hpp
 *
 */

#ifndef MORIS_ASSERT_FN_LOG_ASSERT_HPP_
#define MORIS_ASSERT_FN_LOG_ASSERT_HPP_

// C++ header files.
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <cstdio>
#include <string>

// MORIS header files.
#include "core.hpp"
#include "cl_Logger.hpp"

// ----------------------------------------------------------------------------

namespace moris::assert
{
    /**
     * Overloaded moris::assert::error.
     *
     * @param[in] msg Error message.
     */
    // template< typename Exception = std::runtime_error >
    inline void
    error( std::string const & msg )
    {
        MORIS_LOG_ERROR( "*** Error: " );
        MORIS_LOG_ERROR( "%s", msg.c_str() );

        throw;
    }

    /**
     * Log error message with global logger.
     *
     * @param[in] location  Location where error was raised.
     * @param[in] task      Task where error was raised.
     * @param[in] check     Check that raised assertion.
     * @param[in] exception Exception raised by check.
     */
    // template< typename Exception >
    inline void
    error(
            std::string const &        location,
            std::string const &        task,
            std::string const &        check,
            std::runtime_error const & exception )
    {
        std::string tString;

        tString = "\n";
        tString = "*** ---------------------------------------------------------------------------\n";
        tString += "*** \n";
        tString += "*** Error - Moris encountered an error.";
        tString += "***\n";
        tString += "*** Unable to ";
        tString += task.c_str();
        tString += ".\n";
        tString += "*** Reason:  ";
        tString += check;
        tString += "\n";

        std::istringstream exception_msg( exception.what() );
        std::string        exception_line;

        while ( std::getline( exception_msg, exception_line ) )
        {
            tString += "***          ";
            tString += exception_line.c_str();
        }
        tString += "\n*** Where:   This error was encountered on processor ";
        tString += std::to_string( gLogger.logger_par_rank() );
        tString += " inside ";
        tString += location.c_str();
        tString += ".\n";
        tString += "***\n";
        tString += "*** ---------------------------------------------------------------------------";

        MORIS_LOG_ERROR( "%s", tString.c_str() );

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
    // template<typename Exception>
    inline void
    moris_assert(
            std::string const &        file,
            moris::size_t const &      line,
            std::string const &        function,
            std::string const &        check,
            std::runtime_error const & exception )
    {
        std::stringstream location;
        location << file << ":" << line << " (line " << line << ")";

        std::stringstream task;
        task << "complete call to function " << function << "()";

        std::stringstream reason;
        reason << "Assertion " << check << " failed.";

        moris::assert::error( location.str(), task.str(), reason.str(), exception );
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
    template< typename... Args >
    inline void
    moris_assert(
            const std::string&   file,
            const moris::size_t& line,
            const std::string&   function,
            const std::string&   check,
            const char*          format,
            const Args... aArgs )
    {
        // initialize string for assert message
        std::string tString;

        // check for print argument
        if constexpr ( sizeof...( Args ) != 0 )
        {
            // Determine size of string
            auto tSize = snprintf( nullptr, 0, format, aArgs... );

            // create char pointer with size of string length + 1 for \0
            std::unique_ptr< char[] > tMsg( new char[ tSize + 1 ] );

            // write string into buffered char pointer
            snprintf( tMsg.get(), tSize + 1, format, aArgs... );

            tString = std::string( tMsg.get(), tMsg.get() + tSize );
        }
        else
        {
            tString = format;
        }

        moris::assert::moris_assert(
                file,
                line,
                function,
                check,
                std::runtime_error( tString.c_str() ) );
    }

    /**
     * moris::assert::moris_warning.
     *
     * @param[in] file      File where assertion was raised.
     * @param[in] line      Line where assertion was raised.
     * @param[in] function  Function where assertion was raised.
     * @param[in] check     Check that raised assertion.
     * @param[in] msg       Error message to build exception.
     */
    template< typename... Args >
    inline void
    moris_warning(
            std::string const &   file,
            moris::size_t const & line,
            std::string const &   function,
            std::string const &   check,
            const char*           format,
            const Args... aArgs )
    {
        std::stringstream location;
        location << file << ":" << line << " (line " << line << ")";

        std::stringstream task;
        task << "call to function " << function << "()";

        std::stringstream reason;
        reason << "Assertion " << check << " may indicate an problem.";

        // initialize string for warning message
        std::string tString;

        // check for print argument
        if constexpr ( sizeof...( Args ) != 0 )
        {
            // Determine size of string
            auto tSize = snprintf( nullptr, 0, format, aArgs... );

            // create char pointer with size of string length + 1 for \0
            std::unique_ptr< char[] > tMsg( new char[ tSize + 1 ] );

            // write string into buffered char pointer
            snprintf( tMsg.get(), tSize + 1, format, aArgs... );

            tString = std::string( tMsg.get(), tMsg.get() + tSize );
        }
        else
        {
            tString = format;
        }

        const std::runtime_error exception( tString.c_str() );

        // build information for printing
        tString = "\n";
        tString = "*** ---------------------------------------------------------------------------\n";
        tString += "*** \n";
        tString += "*** Warning - Moris encountered a potential issue.";
        tString += "***\n";
        tString += "*** Issue:   Potential issue in ";
        tString += task.str();
        tString += ".\n";
        tString += "*** Reason:  ";
        tString += check;
        tString += "\n";

        std::istringstream exception_msg( exception.what() );
        std::string        exception_line;

        while ( std::getline( exception_msg, exception_line ) )
        {
            tString += "***          ";
            tString += exception_line.c_str();
        }

        tString += "\n*** Where:   This warning was encountered on processor ";
        tString += std::to_string( gLogger.logger_par_rank() );
        tString += " inside ";
        tString += location.str();
        tString += ".\n";
        tString += "***\n";
        tString += "*** ---------------------------------------------------------------------------";

        MORIS_LOG_WARNING( "%s", tString.c_str() );
    }
}    // namespace moris::assert

#endif /* MORIS_ASSERT_FN_LOG_ASSERT_HPP_ */
