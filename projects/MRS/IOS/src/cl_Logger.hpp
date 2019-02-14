#ifndef MORIS_IOS_CL_LOGGER_HPP_
#define MORIS_IOS_CL_LOGGER_HPP_
#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <cstring>

#include "typedefs.hpp"
#include "IO_Tools.hpp"



namespace moris
{
    class Logger
    {
    public:
        std::ofstream mStream;

        moris::sint mSeverityLevel = -1;

        bool mWriteToAscii = false;

    public:
        Logger(){};

        ~Logger()
        {
            if( mWriteToAscii )
            {
                //close file
                mStream.close();
            }
        };

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

        void initialize( const moris::sint aSverityLevel )
        {
            mSeverityLevel = aSverityLevel;
        };

        void initialize( const std::string aPath,
                         const moris::sint aSverityLevel = -1)
        {
            mSeverityLevel = aSverityLevel;

            mStream.open(aPath,std::ofstream::out ),

            mWriteToAscii = true;
       };

        void set_severity_level( const moris::sint aSverityLevel )
        {
            mSeverityLevel = aSverityLevel;
        };

        moris::sint get_severity_level()
        {
            return mSeverityLevel;
        };

        template< typename ... Args >
        void log_section( const Args ... aArgs )
        {
            // MORIS_ASSERT(mStream.is_open(),"Logger error, the output file stream ofstream is not open);
            std::string tString = print_log( aArgs ... );

            char tTempString[1000];
            std::strcpy( tTempString, "===============================================================================\n" );
            std::strcat( tTempString, "\n" );
            std::strcat( tTempString, tString.c_str() );
            std::strcat( tTempString, "\n" );
            std::strcat( tTempString, "===============================================================================" );

            std::cout << tTempString;

            if( mWriteToAscii )
            {
                mStream << tTempString;
            }
        }

        template< typename ... Args >
        void log( const Args ... aArgs )
        {
           // MORIS_ASSERT(mStream.is_open(),"Logger error, the output file stream ofstream is not open);

            std::string tString = print_log( aArgs ... );

            std::cout << tString;

            if( mWriteToAscii )
            {
                mStream << tString;
            }
        }

        template< typename ... Args >
        void log_info( const Args ... aArgs )
        {
           // MORIS_ASSERT(mStream.is_open(),"Logger error, the output file stream ofstream is not open);

            std::string tString = print_log( aArgs ... );

            std::cout << tString;

            if( mWriteToAscii )
            {
                mStream << tString;
            }
        }

        template< typename ... Args >
        void log_debug( const Args ... aArgs )
        {
           // MORIS_ASSERT(mStream.is_open(),"Logger error, the output file stream ofstream is not open);

            std::string tString = print_log( aArgs ... );

            std::cout << tString;

            if( mWriteToAscii )
            {
                mStream << tString;
            }
        }

        template< typename ... Args >
        void log_warning( const Args ... aArgs )
        {
           // MORIS_ASSERT(mStream.is_open(),"Logger error, the output file stream ofstream is not open);

            std::string tString = print_log( aArgs ... );

            std::cout << tString;

            if( mWriteToAscii )
            {
                mStream << tString;
            }
        }

        template< typename ... Args >
        void log_error( const Args ... aArgs )
        {
           // MORIS_ASSERT(mStream.is_open(),"Logger error, the output file stream ofstream is not open);

            std::string tString = print_log( aArgs ... );

            std::cout << tString;

            if( mWriteToAscii )
            {
                mStream << tString;
            }
        }
    };
}

extern moris::Logger gLogger;
/**
 * Log an section severity message.
 *
 * @include "IOS/src/cl_Logger/log.inc"
 */

#define MORIS_SECTION( ... ) \
    do \
    { \
        gLogger.log_section( __VA_ARGS__ ); \
    } while (false)

/**
 * Log an informational severity message.
 *
 * @include "IOS/src/cl_Logger/log.inc"
 */
#define MORIS_LOG( ... ) \
    do \
    { \
        gLogger.log( __VA_ARGS__ ); \
    } while (false)

/**
 * Log an informational severity message.
 *
 * @include "IOS/src/cl_Logger/log.inc"
 */
#define MORIS_LOG_INFO( ... ) \
    do \
    { \
        if ( gLogger.get_severity_level() < 1 )\
        {\
           gLogger.log_info( __VA_ARGS__ ); \
        }\
    } while (false)

/**
 * Log a debug severity message.
 *
 * @include "IOS/src/cl_Logger/log_debug.inc"
 */
#if !defined(NDEBUG) || defined(DEBUG)
#define MORIS_LOG_DEBUG( ... ) \
    do \
    { \
        gLogger.log_debug( __VA_ARGS__ ); \
    } while (false)
#else
#define MORIS_LOG_DEBUG( ... )
#endif

/**
 * Log a warning severity message.
 *
 * @include "IOS/src/cl_Logger/log_warning.inc"
 */
#define MORIS_LOG_WARNING( ... ) \
    do \
    { \
        if ( gLogger.get_severity_level() < 2 ) \
        {\
            gLogger.log_warning( __VA_ARGS__ ); \
        }\
    } while (false)

/**
 * Log an error severity message.
 *
 * @include "IOS/src/cl_Logger/log_error.inc"
 */
#define MORIS_LOG_ERROR( ... )\
    do \
    { \
        if ( gLogger.get_severity_level() < 2 ) \
        {\
            gLogger.log_warning( __VA_ARGS__ ); \
        }\
    } while (false)

#endif	/* MORIS_IOS_CL_LOGGER_HPP_ */
