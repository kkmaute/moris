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
//#include "cl_GlobalClock.hpp" // MRS/IOS/src
//#include "cl_Tracer_Enums.hpp"

//class GlobalClock;

namespace moris
{
    class Logger
    {
    public:
        std::ofstream mStream;

        moris::sint mSeverityLevel = -1;

        bool mWriteToAscii = false;

//        GlobalClock mGlobalClock;


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
        void log_section( const Args ... aArgs );

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

//        // FUNCTIONS ENABLING TRACING AND CLOCK LOGGING ------------------------------------------------- //
//
//        // sign in
//        void sign_in( enum EntityBase aEntityBase, enum EntityType aEntityType, enum EntityAction aEntityAction )
//        {
//            // pass save info to clock
//            mGlobalClock.sign_in(aEntityBase, aEntityType, aEntityAction);
//
//            // formated output to log file
//            this->log_clock( OutputSpecifier::SignIn, 1.0);
//        }
//
//        // signing out
//        void sign_out()
//        {
//            // log current position in code
//            this->log(OutputSpecifier::ElapsedTime, ( (moris::real) std::clock() - mTimeStamps(mIndentationLevel) ) / CLOCKS_PER_SEC );
//
//            // decrement clock
//            mGlobalClock.sign_out();
//        }
//
//        // logging operation using Clock Info, Output value is allowed to be any common number type
//        template <class T>
//        void log_clock(enum OutputSpecifier aOutputSpecifier, T aOutputValue)
//        {
//            std::string tLine = std::to_string(mGlobalClock.mIndentationLevel) + ";"
//                                + mCurrentFunctionID(mGlobalClock.mIndentationLevel) + ";"
//                                + get_enum_str(mGlobalClock.mCurrentEntity(mGlobalClock.mIndentationLevel)) + ";"
//                                + get_enum_str(mGlobalClock.mCurrentType(mGlobalClock.mIndentationLevel)) + ";"
//                                + get_enum_str(mGlobalClock.mCurrentAction(mGlobalClock.mIndentationLevel)) + ";"
//                                + get_enum_str(mGlobalClock.aOutputSpecifier) + ";"
//                                + std::to_string(aOutputValue) + "\n";
//
//            if( mWriteToAscii )
//            {
//                mStream << tLine;
//            }
//        }
//
//        // logging operation using Clock Info and string
//        void log_clock(std::string aOutputString)
//        {
//            std::string tLine = std::to_string(mGlobalClock.mIndentationLevel) + ";"
//                                + mCurrentFunctionID(mGlobalClock.mIndentationLevel) + ";"
//                                + get_enum_str(mGlobalClock.mCurrentEntity(mGlobalClock.mIndentationLevel)) + ";"
//                                + get_enum_str(mGlobalClock.mCurrentType(mGlobalClock.mIndentationLevel)) + ";"
//                                + get_enum_str(mGlobalClock.mCurrentAction(mGlobalClock.mIndentationLevel)) + ";"
//                                + get_enum_str(OutputSpecifier::FreeText) + ";"
//                                + aOutputString + "\n";
//
//            if( mWriteToAscii )
//            {
//                mStream << tLine;
//            }
//        }

        // -------------------------------------------------------------------------------- //
        // print header
        void print_header()
        {

            mStream << ">>>";

            // get date and time at runtime
            time_t tTimeStamp = time(NULL);

            // print info about test
            mStream << "\n--- RUN ---\n";
    //        mStream << "Executable: " << argv[ 0 ] << " \n";
    //        mStream << "Arguments: ";
    //        for (int ia = 1; ia < argc; ++ia)
    //        {
    //            mStream << argv[ ia ] << "; ";
    //        }
    //        mStream << "Number of Processors used: " << ( int ) par_size() << " \n";
            mStream << "Date of execution: " << ctime(&tTimeStamp) ;


            // print user information
            mStream << "\n--- HOST ---\n";
            mStream << "User: " << std::getenv( "USER" ) << "\n";
            mStream << "Host: " << std::getenv( "HOSTNAME" ) << "\n";
            mStream << "Operating System: " << std::getenv( "OSTYPE" ) << "\n";
            mStream << "Host Type: " << std::getenv( "HOSTTYPE" ) << "\n";


            // print build info
            mStream << "\n--- BUILD ---\n";
            mStream << "Build date: " << __DATE__ << "; " << __TIME__ << "\n";
            mStream << "DEBUG: ";

    #if !defined(NDEBUG) || defined(DEBUG)
            mStream << "ON\n";
    #endif
    #if defined(NDEBUG) || !defined(DEBUG)
            mStream << "OFF\n";
    #endif

            mStream << "Matrix Library: ";

    #ifdef MORIS_USE_ARMA
            mStream << "ARMA \n";
    #endif
    #ifdef MORIS_USE_EIGEN
            mStream << "Eigen \n";
    #endif

            mStream << "<<<\n";


    //        // Top row in table
    //        mStream << "_________________________________________________________________________________________________________\n";
    //        mStream << "Indentation Level ; Function ID ; Entity ; Entity Type ; Entity Action ; Output Specifier ; Output Value \n";
    //        mStream << "---------------------------------------------------------------------------------------------------------\n";

        }

    };

    template< typename ... Args >
    void Logger::log_section( const Args ... aArgs )
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
