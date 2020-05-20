#include "cl_Logger.hpp"

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <cstring>

#include <sstream>
#include <iomanip>
#include <limits>

#include "typedefs.hpp"
#include "IO_Tools.hpp"
#include "paths.hpp"

// for the global clock
#include "cl_GlobalClock.hpp" // MRS/IOS/src
#include "cl_Tracer_Enums.hpp"

#include "fn_stringify.hpp"
#include "Log_Constants.hpp"

namespace moris
{


// -------------------------------------------------------------------------------------------------- //
// log with specified output type
//template <class T>
//void Logger::log_specific(enum OutputSpecifier aOutputSpecifier, T aOutputValue)

// -------------------------------------------------------------------------------------------------- //
// initialize Logger based arguments inputed by user
void Logger::initialize( int  & argc, char * argv[] )
{
    // default parameters
    mSeverityLevel = LOGGER_DEFAULT_SEVERITY_LEVEL;
    mDirectOutputFormat = LOGGER_DEFAULT_DIRECT_OUTPUT;

    // go through user arguments and look for flags
    for( int k=0; k<argc; ++k )
    {
        // user requests output log file
        if ( std::string( argv[ k ] ) == "--outputlog" || std::string( argv[ k ] ) == "-ol" )
        {
            mStream.open( std::string( argv[ k+1 ] ) , std::ofstream::out );
            mWriteToAscii = true;
            //std::cout << "\n Logger: writing to: " << std::string( argv[ k+1 ] ) << "\n";
        }

        // user sets severity level
        if ( std::string( argv[ k ] ) == "--severity" || std::string( argv[ k ] ) == "-sl" )
        {
            mSeverityLevel = std::stoi( std::string( argv[ k+1 ] ) );
            //std::cout << "\n Logger: setting severity level to: " << std::stoi( std::string( argv[ k+1 ] ) ) << "\n";
        }

        // user sets format for output to console
        if ( std::string( argv[ k ] ) == "--directoutput" || std::string( argv[ k ] ) == "-do" )
        {
            uint Integer = std::stoi( std::string( argv[ k+1 ] ) );
            if ( (Integer == 1) || (Integer == 2) || (Integer == 3) )
                mDirectOutputFormat = Integer;
            else
                std::cout << "\n <MRS::IOS::cl_Logger::initialize()>: Unknown direct output format, using standard mode. \n";
        }

    } // end for each input argument

    // print header
    this->print_header();

    // log start of Global Clock to file
    if( mWriteToAscii )
    {
        // formated output to log file
        this->log_to_file( OutputSpecifier::SignIn, 1.0);
    }
}


// -------------------------------------------------------------------------------------------------- //
// FUNCTIONS ENABLING TRACING AND CLOCK LOGGING ----------------------------------------------------- //
// -------------------------------------------------------------------------------------------------- //

// -------------------------------------------------------------------------------------------------- //
// sign in
void Logger::sign_in( enum EntityBase aEntityBase, enum EntityType aEntityType, enum EntityAction aEntityAction )
{
    // pass save info to clock
    mGlobalClock.sign_in(aEntityBase, aEntityType, aEntityAction);

    // log to file
    if( mWriteToAscii )
    {
        // formated output to log file
        this->log_to_file( OutputSpecifier::SignIn, 1.0);
    }

    if (mSeverityLevel < 2)
    {
        // switch based on OutputFormat provided
        if ((mDirectOutputFormat == 3) || (mDirectOutputFormat == 2))
        {
            std::cout << print_empty_line(mGlobalClock.mIndentationLevel - 1)
                      << "__Start: " << get_enum_str(aEntityBase)
                      << " - " << get_enum_str(aEntityType)
                      << " - " << get_enum_str(aEntityAction)
                      << " \n";
            std::cout << print_empty_line(mGlobalClock.mIndentationLevel) << " \n";
        }
        else
            std::cout << "Signing in: " << get_enum_str(aEntityBase)
                      << " - " << get_enum_str(aEntityType)
                      << " - " << get_enum_str(aEntityAction)
                      << " \n";
    }
}


// -------------------------------------------------------------------------------------------------- //
// signing out
void Logger::sign_out()
{
    // log to file
    if( mWriteToAscii )
    {
        // log current position in code
        this->log_to_file(OutputSpecifier::ElapsedTime,
                ( (moris::real) std::clock() - mGlobalClock.mTimeStamps(mGlobalClock.mIndentationLevel) ) / CLOCKS_PER_SEC );
    }

    // log to console
    if (mSeverityLevel < 2)
    {
        // switch based on OutputFormat provided
        if ((mDirectOutputFormat == 3) || (mDirectOutputFormat == 2))
        {
            std::cout << print_empty_line(mGlobalClock.mIndentationLevel)
                      << "_ElapsedTime = "
                      << ( (moris::real) std::clock() - mGlobalClock.mTimeStamps(mGlobalClock.mIndentationLevel) ) / CLOCKS_PER_SEC
                      << " \n";
            std::cout << print_empty_line(mGlobalClock.mIndentationLevel - 1) << " \n";
        }
        else
            std::cout << "Signing out current instance. Elapsed time = "
                      << ( (moris::real) std::clock() - mGlobalClock.mTimeStamps(mGlobalClock.mIndentationLevel) ) / CLOCKS_PER_SEC
                      << " \n";
    }

    // decrement clock
    mGlobalClock.sign_out();
}

// -------------------------------------------------------------------------------------------------- //
// logging operation using Clock Info, Output value is allowed to be any common number type
//void Logger::log_to_file(enum OutputSpecifier aOutputSpecifier, T aOutputValue)


// -------------------------------------------------------------------------------------------------- //
// logging operation using Clock Info, using multiple output values
template <class T1, class T2>
void Logger::log2_to_file(enum OutputSpecifier aOutputSpecifier1, T1 aOutputValue1,
                          enum OutputSpecifier aOutputSpecifier2, T2 aOutputValue2)
{
    this->log_to_file(aOutputSpecifier1, aOutputValue1);
    this->log_to_file(aOutputSpecifier2, aOutputValue2);
}


template <class T1, class T2, class T3>
void Logger::log3_to_file(enum OutputSpecifier aOutputSpecifier1, T1 aOutputValue1,
                          enum OutputSpecifier aOutputSpecifier2, T2 aOutputValue2,
                          enum OutputSpecifier aOutputSpecifier3, T3 aOutputValue3)
{
    this->log_to_file(aOutputSpecifier1, aOutputValue1);
    this->log_to_file(aOutputSpecifier2, aOutputValue2);
    this->log_to_file(aOutputSpecifier3, aOutputValue3);
}


// -------------------------------------------------------------------------------------------------- //
// logging operations for specific types of output texts

void Logger::log_to_file(std::string aOutputString)
{
    this->log_to_file(OutputSpecifier::FreeText, aOutputString);
}
void Logger::log_to_file_info(std::string aOutputString)
{
    this->log_to_file(OutputSpecifier::InfoText, aOutputString);
}
void Logger::log_to_file_debug(std::string aOutputString)
{
    this->log_to_file(OutputSpecifier::DebugText, aOutputString);
}
void Logger::log_to_file_warning(std::string aOutputString)
{
    this->log_to_file(OutputSpecifier::Warning, aOutputString);
}
void Logger::log_to_file_error(std::string aOutputString)
{
    this->log_to_file(OutputSpecifier::Error, aOutputString);
}

// -------------------------------------------------------------------------------------------------- //
// FORMATTING TOOLS FOR OUTPUT ---------------------------------------------------------------------- //
// -------------------------------------------------------------------------------------------------- //

// -------------------------------------------------------------------------------------------------- //
// print header
void Logger::print_header()
{

    mStream << LOGGER_HEADER_BEGIN;

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
    mStream << "Severity Level: " << mSeverityLevel << " \n";


    // print user information
    mStream << "\n--- HOST ---\n";
    mStream << "User: " << std::getenv( "USER" ) << "\n";
    mStream << "Host: " << std::getenv( "HOSTNAME" ) << "\n";
    mStream << "Operating System: " << std::getenv( "OSTYPE" ) << "\n";
    mStream << "Host Type: " << std::getenv( "HOSTTYPE" ) << "\n";


    // print build info
    mStream << "\n--- BUILD ---\n";
    mStream << "Build date: " << __DATE__ << "; " << __TIME__ << "\n";
    mStream << "Git branch: " << get_moris_git_branch() << "\n";
    mStream << "Git hash: " << get_moris_git_hash() << "\n";
    mStream << "DEBUG: ";

#if defined(DEBUG)
    mStream << "ON\n";
#endif
#if !defined(DEBUG)
    mStream << "OFF\n";
#endif

    mStream << "Matrix Library: ";

#ifdef MORIS_USE_ARMA
    mStream << "ARMA \n";
#endif
#ifdef MORIS_USE_EIGEN
    mStream << "Eigen \n";
#endif

    mStream << LOGGER_HEADER_END << "\n";

    // Top row in table
    mStream << "_________________________________________________________________________________________________________\n";
    mStream << "Indentation Level ; Function ID ; Entity ; Entity Type ; Entity Action ; Output Specifier ; Output Value \n";
    mStream << "---------------------------------------------------------------------------------------------------------\n";
}


// -------------------------------------------------------------------------------------------------- //
// prints empty line with vertical markers according to indentation level
std::string Logger::print_empty_line(uint aIndentationLevel)
{
    // initialize
    std::string tEmptyLine = "|";

    // go through all lines and get max of cell array
    for(uint ii = 1; ii <= aIndentationLevel; ii++)
    {
        tEmptyLine.append("  |");
    }

    // return map that lists the first and last lines for logging for each instance by ID
    return tEmptyLine;
    //mStream << tEmptyLine;
}


} // end namespace moris

