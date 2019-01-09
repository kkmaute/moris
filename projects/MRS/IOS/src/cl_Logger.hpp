#ifndef MORIS_IOS_CL_LOGGER_HPP_
#define MORIS_IOS_CL_LOGGER_HPP_
#include <iostream>
#include <fstream>
#include <cstdio>

#include "IO_Tools.hpp"


namespace moris
{
    class Logger
    {
    public:
        std::ofstream mStream;

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

        void initialize( const std::string aPath )
        {
            mStream.open(aPath,std::ofstream::out ),

            mWriteToAscii=true;
       };

        template< typename ... Args >
        void log( const Args ... aArgs )
        {
           // MORIS_ASSERT(mStream.is_open(),"Logger error, the output file stream ofstream is not open);

            if( mWriteToAscii )
            {
                mStream << print_log( aArgs ... );
            }
        }

        template< typename ... Args >
        void log_function( const Args ... aArgs )
        {
           // MORIS_ASSERT(mStream.is_open(),"Logger error, the output file stream ofstream is not open);

            if( mWriteToAscii )
            {
                mStream << print_log( aArgs ... );
            }
        }

        template< typename ... Args >
        void log_info( const Args ... aArgs )
        {
           // MORIS_ASSERT(mStream.is_open(),"Logger error, the output file stream ofstream is not open);

            std::string tString = print_log( aArgs ... );

            std::cout << tString << std::endl;

            if( mWriteToAscii )
            {
                mStream << tString;
            }
        }

        template< typename ... Args >
        void log_debug( const Args ... aArgs )
        {
           // MORIS_ASSERT(mStream.is_open(),"Logger error, the output file stream ofstream is not open);

            if( mWriteToAscii )
            {
                mStream << print_log( aArgs ... );
            }
        }

        template< typename ... Args >
        void log_trace( const Args ... aArgs )
        {
           // MORIS_ASSERT(mStream.is_open(),"Logger error, the output file stream ofstream is not open);

            if( mWriteToAscii )
            {
                mStream << print_log( aArgs ... );
            }
        }

        template< typename ... Args >
        void log_warning( const Args ... aArgs )
        {
           // MORIS_ASSERT(mStream.is_open(),"Logger error, the output file stream ofstream is not open);

            if( mWriteToAscii )
            {
                mStream << print_log( aArgs ... );
            }
        }

        template< typename ... Args >
        void log_error( const Args ... aArgs )
        {
           // MORIS_ASSERT(mStream.is_open(),"Logger error, the output file stream ofstream is not open);

            if( mWriteToAscii )
            {
                mStream << print_log( aArgs ... );
            }
        }

        template< typename ... Args >
        void log_fatal( const Args ... aArgs )
        {
           // MORIS_ASSERT(mStream.is_open(),"Logger error, the output file stream ofstream is not open);

            if( mWriteToAscii )
            {
                mStream << print_log( aArgs ... );
            }
        }
    };
}
/**
 * @brief Macro for function scope markup.
 * The scope name is constructed with help of compiler
 * and contains current function name.
 *
 * @include "IOS/src/cl_Logger/log_function.inc"
 */
#define MORIS_LOG_FUNCTION( ... ) \
    do \
    { \
        gLogger.log_function( __VA_ARGS__ ); \
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
 * Log a trace severity message.
 *
 * @include "IOS/src/cl_Logger/log_trace.inc"
 */
#define MORIS_LOG_TRACE( ... ) \
    do \
    { \
        gLogger.log_trace( __VA_ARGS__ ); \
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
 * Log an informational severity message.
 *
 * @include "IOS/src/cl_Logger/log_info.inc"
 */
//#define MORIS_LOG_INFO       BOOST_LOG_SEV(moris_logger::get(), trivial::info)

//#define MORIS_LOG_INFO       std::cout

#define MORIS_LOG_INFO( ... ) \
    do \
    { \
        gLogger.log_info( __VA_ARGS__ ); \
    } while (false)

/**
 * Log a warning severity message.
 *
 * @include "IOS/src/cl_Logger/log_warning.inc"
 */
#define MORIS_LOG_WARNING( ... ) \
    do \
    { \
        gLogger.log_warning( __VA_ARGS__ ); \
    } while (false)

/**
 * Log an error severity message.
 *
 * @include "IOS/src/cl_Logger/log_error.inc"
 */
#define MORIS_LOG_ERROR    std::cerr

/**
 * Log a fatal severity message.
 *
 * @include "IOS/src/cl_Logger/log_fatal.inc"
 */
#define MORIS_LOG_FATAL( ... ) \
    do \
    { \
        gLogger.log_fatal( __VA_ARGS__ ); \
    } while (false)


#endif	/* MORIS_IOS_CL_LOGGER_HPP_ */
