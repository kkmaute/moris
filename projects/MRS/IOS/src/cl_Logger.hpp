/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Logger.hpp
 *
 */

#ifndef MORIS_IOS_CL_LOGGER_HPP_
#define MORIS_IOS_CL_LOGGER_HPP_
#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <cstring>

#include "typedefs.hpp"
#include "IO_Tools.hpp"

// for the global clock
#include "cl_GlobalClock.hpp"
#include "fn_stringify.hpp"
#include "Log_Constants.hpp"

// need to be revised later
#include <mpi.h>

namespace moris
{
    class Logger
    {
      public:
        // Output File
        std::ofstream mStream;

        /**
         * @brief Flag to control output of logger information to file
         */
        bool mWriteToAscii = false;

        /**
         * @brief Variable to control which processor writes logger information to screen
         */
        int mOutputRank = 0;

        /*
         * Severity Levels:
         * 1 - everything gets written
         * 2 - print general logs and output values but suppress verbose info logs
         * 3 - suppress everything except errors and sign-in/sign-out
         */
        // decide which outputs get written
        sint mSeverityLevel = LOGGER_DEFAULT_SEVERITY_LEVEL;

        /**
         * @brief Variable to control level of memory reporting
         */
        sint mMemoryOutput = 1;

        /**
         * Output formatting mode for console output:
         * 1 - legacy mode, everything is written to screen as handed to log macros/functions
         * 2 - output is written in a clean tree structure, general MORIS_LOG and MORIS_LOG_INFO is suppressed. Good for live tracing.
         * 3 - all output is written in a tree structure, no filtering
         */
        uint mDirectOutputFormat = LOGGER_DEFAULT_DIRECT_OUTPUT;

        // Clock for tracing and timing
        GlobalClock mGlobalClock;

        uint mIteration = 0;    // FIXME this is absolutely a hack, it doesn't even store the iteration correctly :)

        inline int
        logger_par_rank()
        {
            int tProcRank;
            MPI_Comm_rank( MPI_COMM_WORLD, &tProcRank );
            return tProcRank;
        }

        inline real
        logger_sum_all( real& aLocalInput )
        {
            real aGlobalMax;
            MPI_Allreduce( &aLocalInput, &aGlobalMax, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
            return aGlobalMax;
        }

        inline real
        logger_max_all( real& aLocalInput )
        {
            real aGlobalMax;
            MPI_Allreduce( &aLocalInput, &aGlobalMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
            return aGlobalMax;
        }

        inline real
        logger_min_all( real& aLocalInput )
        {
            real aGlobalMin;
            MPI_Allreduce( &aLocalInput, &aGlobalMin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
            return aGlobalMin;
        }

      public:
        Logger(){};

        //------------------------------------------------------------------------------

        ~Logger()
        {
            // report error if prematurely stopped
            if ( mGlobalClock.mIndentationLevel != 0 )
            {
                // set indentation level to zero
                mGlobalClock.mIndentationLevel = 0;

                // report error
                std::cout << "Global clock & logger prematurely stopped.";

                if ( mWriteToAscii )
                {
                    // log current position in code
                    this->log_to_file_error( "Global clock & logger prematurely stopped." );
                }
            }

            // Stop Global Clock Timer
            real tElapsedTime = ( (moris::real)std::clock() - mGlobalClock.mTimeStamps[ mGlobalClock.mIndentationLevel ] ) / CLOCKS_PER_SEC;

            if ( mWriteToAscii )
            {
                // log runtime of clock
                this->log_to_file( "ElapsedTime", tElapsedTime );

                // close file
                mStream.close();
            }

            // log end of Global Clock to console - only processor mOutputRank prints message
            std::cout << "Global Clock Stopped. ElapsedTime = " << tElapsedTime << " \n"
                      << std::flush;

            if ( PRINT_WALL_TIME )
            {
                std::chrono::duration< double > tChronoElapsedWallTime = ( std::chrono::system_clock::now() - mGlobalClock.mWallTimeStamps[ mGlobalClock.mIndentationLevel ] );
                real                            tElapsedWallTime       = tChronoElapsedWallTime.count();
                std::cout << "Global Clock Stopped. ElapsedWallTime = " << tElapsedWallTime << " \n"
                          << std::flush;
            }
        };

        //------------------------------------------------------------------------------

        /**
         * Initialize MORIS_LOGGER
         *
         * Severity Levels
         * 1 ... all outputs
         * 2 ... all outputs except info
         * 3 ... all outputs except info and warning
         *
         * @include "IOS/src/cl_Logger/log.inc"
         */

        void
        initialize( const moris::sint aSeverityLevel )
        {
            mSeverityLevel = aSeverityLevel;

            // determine memory usage
            std::string tMemUsage = this->memory_usage();

            // log start of Global Clock to console - only processor mOutputRank prints message
            if ( logger_par_rank() == mOutputRank )
            {
                std::cout << "Global Clock Initialized ... \n"
                          << std::flush;

                if ( mMemoryOutput > 0 )
                {
                    std::cout << tMemUsage << std::endl
                              << std::flush;
                }
            }
        };

        //------------------------------------------------------------------------------

        void initialize( int& argc, char* argv[] );

        //------------------------------------------------------------------------------

        void
        initialize(
                const std::string aPath,
                const moris::sint aSeverityLevel      = LOGGER_DEFAULT_SEVERITY_LEVEL,
                const moris::uint aDirectOutputFormat = LOGGER_DEFAULT_DIRECT_OUTPUT )
        {
            mDirectOutputFormat = aDirectOutputFormat;
            mSeverityLevel      = aSeverityLevel;

            mStream.open( aPath + "." + std::to_string( logger_par_rank() ), std::ofstream::out );

            mWriteToAscii = true;

            // print header
            this->print_header();

            // log start of Global Clock to file
            if ( mWriteToAscii )
            {
                // formatted output to log file
                this->log_to_file( "SignIn", 1.0 );
            }

            // determine memory usage
            std::string tMemUsage = this->memory_usage();

            // log start of Global Clock to console - only processor mOutputRank prints message
            if ( logger_par_rank() == mOutputRank )
            {
                std::cout << "Global Clock Initialized ... \n"
                          << std::flush;

                if ( mMemoryOutput > 0 )
                {
                    std::cout << tMemUsage << std::endl
                              << std::flush;
                }
            }
        };

        //------------------------------------------------------------------------------

        void
        set_severity_level( const moris::sint aSeverityLevel )
        {
            mSeverityLevel = aSeverityLevel;
        };

        //------------------------------------------------------------------------------

        moris::sint
        get_severity_level()
        {
            return mSeverityLevel;
        };

        //------------------------------------------------------------------------------

        void
        set_screen_output_rank( const moris::sint aOutputRank )
        {
            mOutputRank = aOutputRank;
        };

        //------------------------------------------------------------------------------

        template< typename... Args >
        void
        log_section( const char* tFormat, const Args... aArgs )
        {
            // build log message
            std::string tString = print_log( tFormat, aArgs... );

            // only processor mOutputRank prints message
            if ( logger_par_rank() == mOutputRank )
            {
                if ( mDirectOutputFormat == 1 )
                {
                    char tTempString[ 1000 ];
                    std::strcpy( tTempString, "===============================================================================\n" );
                    std::strcat( tTempString, "\n" );
                    std::strcat( tTempString, tString.c_str() );
                    std::strcat( tTempString, "\n" );
                    std::strcat( tTempString, "===============================================================================" );

                    std::cout << tTempString << "\n";
                }
            }

            if ( mWriteToAscii )
            {
                this->log_to_file( tString );
            }
        }

        //------------------------------------------------------------------------------

        template< typename... Args >
        void
        log( const char* tFormat, const Args... aArgs )
        {
            // build log message
            std::string tString = print_log( tFormat, aArgs... );

            // only processor mOutputRank prints message
            if ( logger_par_rank() == mOutputRank )
            {
                // switch based on OutputFormat provided
                if ( mDirectOutputFormat == 3 )
                {
                    std::cout << print_empty_line( mGlobalClock.mIndentationLevel ) << "_"    //
                              << "Log: " << tString << " \n";
                }
                else if ( mDirectOutputFormat == 1 )
                {
                    std::cout << tString << "\n";
                }
            }

            if ( mWriteToAscii )
            {
                this->log_to_file( tString );
            }
        }

        //------------------------------------------------------------------------------

        template< typename... Args >
        void
        log_info( const char* tFormat, const Args... aArgs )
        {
            // build log message
            std::string tString = print_log( tFormat, aArgs... );

            // only processor mOutputRank prints message
            if ( logger_par_rank() == mOutputRank )
            {
                // switch based on OutputFormat provided
                if ( mDirectOutputFormat == 3 )
                {
                    // check if Entity Type has been specified, if not use entity base for printing
                    if ( mGlobalClock.mCurrentType[ mGlobalClock.mIndentationLevel ] == LOGGER_NON_SPECIFIC_ENTITY_TYPE )
                    {
                        std::cout << print_empty_line( mGlobalClock.mIndentationLevel ) << "_"                 //
                                  << mGlobalClock.mCurrentEntity[ mGlobalClock.mIndentationLevel ] << " - "    //
                                  << tString << " \n"
                                  << std::flush;
                    }
                    else    // if yes use use entity type for printing
                    {
                        std::cout << print_empty_line( mGlobalClock.mIndentationLevel ) << "_"               //
                                  << mGlobalClock.mCurrentType[ mGlobalClock.mIndentationLevel ] << " - "    //
                                  << tString << " \n"
                                  << std::flush;
                    }
                }
                else if ( mDirectOutputFormat == 1 )
                {
                    std::cout << tString << "\n"
                              << std::flush;
                }
            }

            // write to file if requested
            if ( mWriteToAscii )
            {
                this->log_to_file_info( tString );
            }
        }

        //------------------------------------------------------------------------------

        // FIXME: design this function in a more clean manner
        template< typename... Args >
        void
        log_info_all_proc( const char* tFormat, const Args... aArgs )
        {
            int tOutputRank = mOutputRank;
            mOutputRank = logger_par_rank();
            this->log_info( tFormat, aArgs... );
            mOutputRank = tOutputRank;
        }

        //------------------------------------------------------------------------------

        template< typename... Args >
        void
        log_debug( const char* tFormat, const Args... aArgs )
        {
            // build log message
            std::string tString = print_log( tFormat, aArgs... );

            // only processor mOutputRank prints message
            if ( logger_par_rank() == mOutputRank )
            {
                // switch based on OutputFormat provided
                if ( ( mDirectOutputFormat == 3 ) || ( mDirectOutputFormat == 2 ) )
                {
                    std::cout << print_empty_line( mGlobalClock.mIndentationLevel ) << "_"
                              << "Debug: " << tString << " \n"
                              << std::flush;
                }
                else
                {
                    std::cout << tString << "\n"
                              << std::flush;
                }
            }

            // write to file if requested
            if ( mWriteToAscii )
            {
                this->log_to_file_debug( tString );
            }
        }

        //------------------------------------------------------------------------------

        template< typename... Args >
        void
        log_warning( const char* tFormat, const Args... aArgs )
        {
            // build log message
            std::string tString = print_log( tFormat, aArgs... );

            // only processor mOutputRank prints message
            if ( logger_par_rank() == mOutputRank )
            {
                // switch based on OutputFormat provided
                if ( ( mDirectOutputFormat == 3 ) || ( mDirectOutputFormat == 2 ) )
                {
                    std::cout << print_empty_line( mGlobalClock.mIndentationLevel ) << "_"
                              << mGlobalClock.mCurrentEntity[ mGlobalClock.mIndentationLevel ] << " - "
                              << "Proc #" << logger_par_rank() << " - WARNING: " << tString << " \n"
                              << std::flush;
                }
                else
                {
                    std::cout << "Proc #" << logger_par_rank() << " - WARNING: " << tString << "\n"
                              << std::flush;
                }
            }

            // write to file if requested
            if ( mWriteToAscii )
            {
                this->log_to_file_warning( tString );
            }
        }

        //------------------------------------------------------------------------------

        template< typename... Args >
        void
        log_error( const char* tFormat, const Args... aArgs )
        {
            // build log message
            std::string tString = print_log( tFormat, aArgs... );

            // switch based on OutputFormat provided
            if ( ( mDirectOutputFormat == 3 ) || ( mDirectOutputFormat == 2 ) )
            {
                std::cout << "Processor = " << logger_par_rank() << " : " << tString << "\n"
                          << std::flush;
            }
            else
            {
                std::cout << "Processor = " << logger_par_rank() << " : " << tString << "\n"
                          << std::flush;
            }

            // write to file if requested
            if ( mWriteToAscii )
            {
                this->log_to_file_error( tString );
            }
        }

        // ---------------------------------------------------------------------------------------------- //
        // FUNCTIONS ENABLING TRACING AND CLOCK LOGGING ------------------------------------------------- //
        // ---------------------------------------------------------------------------------------------- //

        // log with specified output type
        template< class T >
        void
        log_specific(
                std::string aOutputSpecifier,
                T           aOutputValue )
        {
            // only processor 0 prints message
            if ( logger_par_rank() == mOutputRank )
            {
                // switch based on OutputFormat provided
                if ( ( mDirectOutputFormat == 3 ) || ( mDirectOutputFormat == 2 ) )
                {
                    if ( mGlobalClock.mCurrentType[ mGlobalClock.mIndentationLevel ] == LOGGER_NON_SPECIFIC_ENTITY_TYPE )
                    {
                        std::cout << print_empty_line( mGlobalClock.mIndentationLevel ) << "_"
                                  << mGlobalClock.mCurrentEntity[ mGlobalClock.mIndentationLevel ] << " - "
                                  << aOutputSpecifier << ": "
                                  << ios::stringify( aOutputValue ) << " \n"
                                  << std::flush;
                    }
                    else
                    {
                        std::cout << print_empty_line( mGlobalClock.mIndentationLevel ) << "_"
                                  << mGlobalClock.mCurrentType[ mGlobalClock.mIndentationLevel ] << " - "
                                  << aOutputSpecifier << ": "
                                  << ios::stringify( aOutputValue ) << " \n"
                                  << std::flush;
                    }
                }
                else
                {
                    std::cout << aOutputSpecifier << ": " << ios::stringify( aOutputValue ) << " \n"
                              << std::flush;
                }
            }

            // write to file if requested
            if ( mWriteToAscii )
            {
                this->log_to_file( aOutputSpecifier, aOutputValue );
            }
        }

        //------------------------------------------------------------------------------

        /**
         * Sign in to the logger with an entity action.
         *
         * @param aEntityBase Entity base
         * @param aEntityType Entity type
         * @param aEntityAction Entity action
         */
        void sign_in(
                std::string aEntityBase,
                std::string aEntityType,
                std::string aEntityAction );

        //------------------------------------------------------------------------------

        // signing out
        void sign_out();

        //------------------------------------------------------------------------------

        // checks whether an instance exits
        bool exists(
                const std::string& aEntityBase,
                const std::string& aEntityType,
                const std::string& aEntityAction );

        //------------------------------------------------------------------------------

        // increment iteration count of current instance
        void iterate();

        //------------------------------------------------------------------------------

        /**
         * request/get the iteration id of a logged instance; by default iteration id
         * of first instance is returned.
         *
         * @param aEntityBase Entity base
         * @param aEntityType Entity type
         * @param aEntityAction Entity action
         * @param aGetLastInstance Flag (if true) iteration id of last instance is returned
         */

        uint get_iteration(
                const std::string& aEntityBase,
                const std::string& aEntityType,
                const std::string& aEntityAction,
                bool               aGetLastInstance = false );

        //------------------------------------------------------------------------------

        // set the iteration of a logged instance
        void set_iteration(
                const std::string& aEntityBase,
                const std::string& aEntityType,
                const std::string& aEntityAction,
                const uint         aIter );

        //------------------------------------------------------------------------------

        // get data defined by key for a logged instance
        real
        get_action_data(
                const std::string& aEntityBase,
                const std::string& aEntityType,
                const std::string& aEntityAction,
                const std::string& aEntityData );

        //------------------------------------------------------------------------------

        void
        set_action_data(
                const std::string& aEntityBase,
                const std::string& aEntityType,
                const std::string& aEntityAction,
                const std::string& aEntityData,
                const real&        aActionValue );

        //------------------------------------------------------------------------------

        // request/get the iteration of the optimization algorithm
        uint get_opt_iteration();

        //------------------------------------------------------------------------------

        // request/get the iteration of the optimization algorithm
        void set_opt_iteration( uint aIter );

        //------------------------------------------------------------------------------

        // for logging memory usage
        std::string memory_usage();

        //------------------------------------------------------------------------------

        // write logged info to formatted file
        template< class T >
        void
        log_to_file( std::string aOutputSpecifier, T aOutputValue )
        {
            std::string tLine =
                    ios::stringify( mGlobalClock.mIndentationLevel ) + ";"
                    + ios::stringify( mGlobalClock.mCurrentFunctionID[ mGlobalClock.mIndentationLevel ] ) + ";"
                    + mGlobalClock.mCurrentEntity[ mGlobalClock.mIndentationLevel ] + ";"
                    + mGlobalClock.mCurrentType[ mGlobalClock.mIndentationLevel ] + ";"
                    + mGlobalClock.mCurrentAction[ mGlobalClock.mIndentationLevel ] + ";"
                    + aOutputSpecifier + ";"
                    + ios::stringify( aOutputValue ) + "\n";

            mStream << tLine << std::flush;
        }

        //------------------------------------------------------------------------------

        template< class T1, class T2 >
        void log2_to_file(
                std::string aOutputSpecifier1,
                T1          aOutputValue1,
                std::string aOutputSpecifier2,
                T2          aOutputValue2 );

        //------------------------------------------------------------------------------

        template< class T1, class T2, class T3 >
        void log3_to_file(
                std::string aOutputSpecifier1,
                T1          aOutputValue1,
                std::string aOutputSpecifier2,
                T2          aOutputValue2,
                std::string aOutputSpecifier3,
                T3          aOutputValue3 );

        //------------------------------------------------------------------------------

        // logging operation using Clock Info and string
        void log_to_file( std::string aOutputString );
        void log_to_file_info( std::string aOutputString );
        void log_to_file_debug( std::string aOutputString );
        void log_to_file_warning( std::string aOutputString );
        void log_to_file_error( std::string aOutputString );

        // -----------------------------------------------------------------------------
        // FORMATTING TOOLS FOR OUTPUT -------------------------------------------------
        // -----------------------------------------------------------------------------

        // print header
        void print_header();

        //------------------------------------------------------------------------------

        // print empty line
        std::string print_empty_line( uint aIndentationLevel );

    };    // end class Logger

}    // end namespace moris

// -----------------------------------------------------------------------------
// ---------------------------------------- MACROS -----------------------------
// -----------------------------------------------------------------------------

extern moris::Logger gLogger;

// dummy macro used to check format and print argument
#ifdef __GNUC__
__attribute__( ( format( printf, 1, 2 ) ) )
#endif
void
check_args( const char* tFormat, ... );

/**
 * Log an section severity message.
 *
 * @include "IOS/src/cl_Logger/log.inc"
 */
#define MORIS_SECTION( ... )                    \
    do                                          \
    {                                           \
        if ( false ) check_args( __VA_ARGS__ ); \
        gLogger.log_section( __VA_ARGS__ );     \
    } while ( false )

// -----------------------------------------------------------------------------

/**
 * Log an informational severity message.
 *
 * @include "IOS/src/cl_Logger/log.inc"
 */
#define MORIS_LOG( ... )                            \
    do                                              \
    {                                               \
        if ( gLogger.get_severity_level() < 3 )     \
        {                                           \
            if ( false ) check_args( __VA_ARGS__ ); \
            gLogger.log_info( __VA_ARGS__ );        \
        }                                           \
    } while ( false )

// -----------------------------------------------------------------------------

/**
 * Log an informational severity message.
 *
 * @include "IOS/src/cl_Logger/log.inc"
 */
#define MORIS_LOG_INFO( ... )                       \
    do                                              \
    {                                               \
        if ( gLogger.get_severity_level() < 2 )     \
        {                                           \
            if ( false ) check_args( __VA_ARGS__ ); \
            gLogger.log_info( __VA_ARGS__ );        \
        }                                           \
    } while ( false )

// -----------------------------------------------------------------------------

/**
 * Log an informational severity message and let all processors output to console.
 *
 * @include "IOS/src/cl_Logger/log.inc"
 */
#define MORIS_LOG_INFO_ALL_PROCS( ... )                 \
    do                                                  \
    {                                                   \
        if ( gLogger.get_severity_level() < 2 )         \
        {                                               \
            if ( false ) check_args( __VA_ARGS__ );     \
            gLogger.log_info_all_proc( __VA_ARGS__ );   \
        }                                               \
    } while ( false )

// -----------------------------------------------------------------------------

/**
 * Log single output value with specified type
 *
 * @include "IOS/src/cl_Logger/log.inc"
 */
#define MORIS_LOG_SPEC( ... )                    \
    do                                           \
    {                                            \
        if ( gLogger.get_severity_level() < 3 )  \
        {                                        \
            gLogger.log_specific( __VA_ARGS__ ); \
        }                                        \
    } while ( false )

// -----------------------------------------------------------------------------

/**
 * Log an iteration increment
 *
 * @include "IOS/src/cl_Logger/log.inc"
 */
#define MORIS_LOG_ITERATION( ... )      \
    do                                  \
    {                                   \
        gLogger.iterate( __VA_ARGS__ ); \
    } while ( false )

// -----------------------------------------------------------------------------

/**
 * Log a debug severity message.
 *
 * @include "IOS/src/cl_Logger/log_debug.inc"
 */
#if defined( MORIS_HAVE_DEBUG )
#define MORIS_LOG_DEBUG( ... )                  \
    do                                          \
    {                                           \
        if ( false ) check_args( __VA_ARGS__ ); \
        gLogger.log_debug( __VA_ARGS__ );       \
    } while ( false )
#else
#define MORIS_LOG_DEBUG( ... )
#endif

// -----------------------------------------------------------------------------

/**
 * Log a warning severity message.
 *
 * @include "IOS/src/cl_Logger/log_warning.inc"
 */
#define MORIS_LOG_WARNING( ... )                    \
    do                                              \
    {                                               \
        if ( gLogger.get_severity_level() < 3 )     \
        {                                           \
            if ( false ) check_args( __VA_ARGS__ ); \
            gLogger.log_warning( __VA_ARGS__ );     \
        }                                           \
    } while ( false )

// -----------------------------------------------------------------------------

/**
 * Log an error severity message.
 *
 * @include "IOS/src/cl_Logger/log_error.inc"
 */
#define MORIS_LOG_ERROR( ... )                      \
    do                                              \
    {                                               \
        if ( gLogger.get_severity_level() < 4 )     \
        {                                           \
            if ( false ) check_args( __VA_ARGS__ ); \
            gLogger.log_error( __VA_ARGS__ );       \
        }                                           \
    } while ( false )

#endif /* MORIS_IOS_CL_LOGGER_HPP_ */
