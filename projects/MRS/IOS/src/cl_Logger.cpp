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
    // -----------------------------------------------------------------------------

    // log with specified output type
    //template <class T>
    //void Logger::log_specific( std::string aOutputSpecifier, T aOutputValue )

    // -----------------------------------------------------------------------------

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
                mStream.open( std::string( argv[ k+1 ] ) + "." + std::to_string( logger_par_rank() ) , std::ofstream::out );
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
            this->log_to_file( "SignIn" , 1.0);
        }

        // log start of Global Clock to console - only processor mOutputRank prints message
        if ( logger_par_rank() == mOutputRank )
        {
            std::cout << "Global Clock Initialized ... \n" << std::flush;
        }

    }

    // -----------------------------------------------------------------------------
    // FUNCTIONS ENABLING TRACING AND CLOCK LOGGING --------------------------------
    // -----------------------------------------------------------------------------

    void Logger::sign_in(
            std::string aEntityBase,
            std::string aEntityType,
            std::string aEntityAction)
    {
        // pass save info to clock
        mGlobalClock.sign_in(aEntityBase, aEntityType, aEntityAction);

        // log to file
        if( mWriteToAscii )
        {
            // formated output to log file
            this->log_to_file( "SignIn", 1.0 );
        }

        // log to console - only processor mOutputRank prints message
        if ( logger_par_rank() == mOutputRank )
        {
            if (mSeverityLevel < 1)
            {
                // switch based on OutputFormat provided
                if ((mDirectOutputFormat == 3) || (mDirectOutputFormat == 2))
                {
                    std::cout << print_empty_line(mGlobalClock.mIndentationLevel - 1) << " \n";
                    std::cout << print_empty_line(mGlobalClock.mIndentationLevel - 1)
                            << "__" << aEntityBase
                            << " - " << aEntityType
                            << " - " << aEntityAction
                            << " \n";
                }
                else
                {
                    std::cout << "Signing in: " << aEntityBase
                            << " - " << aEntityType
                            << " - " << aEntityAction
                            << " \n";
                }
            }
        }
    }

    // -----------------------------------------------------------------------------

    // signing out
    void Logger::sign_out()
    {
        // stop timer
        real tElapsedTime = ( (moris::real) std::clock() - mGlobalClock.mTimeStamps[mGlobalClock.mIndentationLevel] ) / CLOCKS_PER_SEC;

        // compute maximum and minimum time used by processors
        real tElapsedTimeMax = logger_max_all(tElapsedTime);
        real tElapsedTimeMin = logger_min_all(tElapsedTime);

        // debug - wall clock time measurements
        real tElapsedWallTimeMax = 0.0;
        real tElapsedWallTimeMin = 0.0;
        if ( PRINT_WALL_TIME )
        {
            std::chrono::duration<double> tChronoElapsedWallTime = ( std::chrono::system_clock::now() - mGlobalClock.mWallTimeStamps[mGlobalClock.mIndentationLevel] );
            real tElapsedWallTime = tChronoElapsedWallTime.count();
            //std::cout << "Proc #" << logger_par_rank() << ": Wall clock time at sign out: " << tElapsedWallTime << " seconds. \n" << std::endl;
            tElapsedWallTimeMax = logger_max_all(tElapsedWallTime);
            tElapsedWallTimeMin = logger_min_all(tElapsedWallTime);
        }

        // log to file
        if( mWriteToAscii )
        {
            // print timing for previous iteration if previous iterations are present
            if ( mGlobalClock.mCurrentIteration[ mGlobalClock.mIndentationLevel ] > 0 )
            {
                // compute iteration time on each proc
                real tIndividualIterationTime = ( (moris::real) std::clock() - mGlobalClock.mIterationTimeStamps[mGlobalClock.mIndentationLevel] ) / CLOCKS_PER_SEC;

                // log iteration time to file
                this->log_to_file( "IterationTime", tIndividualIterationTime );
            }

            // log current position in code
            this->log_to_file( "ElapsedTime", tElapsedTime );
        }
            
        // log to console - only processor mOutputRank prints message
        if ( logger_par_rank() == mOutputRank )
        {
            if (mSeverityLevel < 1)
            {
                // switch based on OutputFormat provided
                if ((mDirectOutputFormat == 3) || (mDirectOutputFormat == 2))
                {
                    // print timing for previous iteration if previous iterations are present
                    if ( mGlobalClock.mCurrentIteration[ mGlobalClock.mIndentationLevel ] > 0 )
                    {
                        // stop timer
                        real tIterationTime = ( (moris::real) std::clock() - mGlobalClock.mIterationTimeStamps[mGlobalClock.mIndentationLevel] ) / CLOCKS_PER_SEC;

                        std::cout << print_empty_line(mGlobalClock.mIndentationLevel) << "_" <<
                                mGlobalClock.mCurrentType[ mGlobalClock.mIndentationLevel ] << " - " <<
                                "IterationTime: " << tIterationTime << " \n" << std::flush;
                    }

                    // print elapsed time for whole entity
                    std::cout << print_empty_line(mGlobalClock.mIndentationLevel) << "_" <<
                            "ElapsedTime (max/min) = " <<
                            tElapsedTimeMax <<
                            " / " <<
                            tElapsedTimeMin <<
                            " \n" << std::flush;

                    if ( PRINT_WALL_TIME )
                        std::cout << print_empty_line(mGlobalClock.mIndentationLevel) << "_" << 
                        "ElapsedWallTime (max/min) = " << 
                        tElapsedWallTimeMax << 
                        " / " <<
                        tElapsedWallTimeMin <<
                        " \n" << std::flush;

                    std::cout << print_empty_line(mGlobalClock.mIndentationLevel - 1) << " \n";
                }
                else                
                {
                    std::cout << "Signing out " <<
                            mGlobalClock.mCurrentType[mGlobalClock.mIndentationLevel] <<
                            ". Elapsed time (max/min) = " <<
                            tElapsedTimeMax <<
                            " / " <<
                            tElapsedTimeMin <<
                           " \n" << std::flush;
                }
            }
        }

        // decrement clock
        mGlobalClock.sign_out();
    }

    //------------------------------------------------------------------------------

    void Logger::iterate()
    {
        // only processor 0 prints message
        if (logger_par_rank() == mOutputRank )
        {
            // print iteration to screen
            // switch based on OutputFormat provided
            if ((mDirectOutputFormat == 3) || (mDirectOutputFormat == 2))
            {
                // print timing for previous iteration if previous iterations are present
                if ( mGlobalClock.mCurrentIteration[ mGlobalClock.mIndentationLevel ] > 0 )
                {
                    // stop timer
                    real tIterationTime = ( (moris::real) std::clock() - mGlobalClock.mIterationTimeStamps[mGlobalClock.mIndentationLevel] ) / CLOCKS_PER_SEC;

                    std::cout << print_empty_line(mGlobalClock.mIndentationLevel) << "_" <<
                            mGlobalClock.mCurrentType[ mGlobalClock.mIndentationLevel ] << " - " <<
                            "IterationTime: " << tIterationTime << " \n" << std::flush;
                }

                std::cout << print_empty_line( mGlobalClock.mIndentationLevel ) <<
                        "===================================================================\n" << std::flush;
                std::cout << print_empty_line( mGlobalClock.mIndentationLevel ) << "_" <<
                        mGlobalClock.mCurrentType[ mGlobalClock.mIndentationLevel ] << " - " <<
                        "Iteration" << ": " <<
                        ios::stringify( mGlobalClock.mCurrentIteration[ mGlobalClock.mIndentationLevel ] + 1 ) <<
                        " \n" << std::flush;
            }
            else
            {
                std::cout << "Iteration" << ": " <<
                        mGlobalClock.mCurrentType[ mGlobalClock.mIndentationLevel ] << " - " <<
                        ios::stringify( mGlobalClock.mCurrentIteration[ mGlobalClock.mIndentationLevel ] + 1 ) <<
                        " \n" << std::flush;
            }

            // write to file if requested
            if( mWriteToAscii )
            {
                // print timing for previous iteration if previous iterations are present
                if ( mGlobalClock.mCurrentIteration[ mGlobalClock.mIndentationLevel ] > 0 )
                {
                    // compute iteration time on each proc
                    real tIndividualIterationTime =
                            ( (moris::real) std::clock() - mGlobalClock.mIterationTimeStamps[mGlobalClock.mIndentationLevel] ) / CLOCKS_PER_SEC;

                    // log iteration time to file
                    this->log_to_file( "IterationTime", tIndividualIterationTime );
                }

                // formated output to log file
                this->log_to_file(
                        "Iteration",
                        mGlobalClock.mCurrentIteration[ mGlobalClock.mIndentationLevel ] + 1 );
            }
        }

        // increment iteration count of current instance inside global clock
        mGlobalClock.iterate();
    }

    //------------------------------------------------------------------------------

    uint Logger::get_iteration(
            std::string aEntityBase,
            std::string aEntityType,
            std::string aEntityAction )
    {
        // initialize
        uint tIndentLevel = 0;
        bool tInstanceFound = false;

        // go through global clock stack from bottom and look for requested instance
        while( tIndentLevel <= mGlobalClock.mIndentationLevel && tInstanceFound == false )
        {
            // check if Instance matches the instance searched for
            if(
                    aEntityBase     == mGlobalClock.mCurrentEntity[ tIndentLevel ] and
                    ( aEntityType   == mGlobalClock.mCurrentType[ tIndentLevel ]   or aEntityType == "Arbitrary" ) and
                    ( aEntityAction == mGlobalClock.mCurrentAction[ tIndentLevel ] or aEntityAction == "Arbitrary" ) )
            {
                tInstanceFound = true;
            }

            // increment cursor
            tIndentLevel++;
        }

        // return iteration count if the instance was found
        if ( tInstanceFound )
        {
            return mGlobalClock.mCurrentIteration[ tIndentLevel - 1 ];
        }

        return 0;
    }

    //------------------------------------------------------------------------------

    void Logger::set_iteration(
            std::string aEntityBase,
            std::string aEntityType,
            std::string aEntityAction,
            uint              aIter)
    {
        // initialize
        uint tIndentLevel = 0;
        bool tInstanceFound = false;

        // go through global clock stack from bottom and look for requested instance
        while( tIndentLevel <= mGlobalClock.mIndentationLevel && tInstanceFound == false )
        {
            // check if Instance matches the instance searched for
            if(
                    aEntityBase     == mGlobalClock.mCurrentEntity[ tIndentLevel ] and
                    ( aEntityType   == mGlobalClock.mCurrentType  [ tIndentLevel ] or aEntityType   == "Arbitrary" ) and
                    ( aEntityAction == mGlobalClock.mCurrentAction[ tIndentLevel ] or aEntityAction == "Arbitrary" ) )
            {
                tInstanceFound = true;
            }

            // increment cursor
            tIndentLevel++;
        }

        // return iteration count if the instance was found
        if ( tInstanceFound )
        {
            mGlobalClock.mCurrentIteration[ tIndentLevel - 1 ] = aIter;
        }
        // throw error if instance was not found
        else
        {
            std::cout << "Logger::set_iteration: Iteration index was not found and could not be set.\n";
            throw;
        }
    }

    //------------------------------------------------------------------------------

     uint Logger::get_opt_iteration()
     {
         // iteration of arbitrary Optimization Algorithm
         // note: this can be done, as there's always only one active opt.-alg.
         return this->get_iteration(
                 "OptimizationAlgorithm" ,
                 "Arbitrary",
                 "Arbitrary" );
     }

     //------------------------------------------------------------------------------

      void Logger::set_opt_iteration( uint aIter)
      {
          // iteration of arbitrary Optimization Algorithm
          // note: this can be done, as there's always only one active opt.-alg.
          this->set_iteration(
                  "OptimizationAlgorithm" ,
                  "Arbitrary",
                  "Arbitrary",
                  aIter);
      }

    // -----------------------------------------------------------------------------

    // logging operation using Clock Info, using multiple output values
    template <class T1, class T2>
    void Logger::log2_to_file(
            std::string aOutputSpecifier1, T1 aOutputValue1,
            std::string aOutputSpecifier2, T2 aOutputValue2)
    {
            this->log_to_file(aOutputSpecifier1, aOutputValue1);
            this->log_to_file(aOutputSpecifier2, aOutputValue2);
    }

    template <class T1, class T2, class T3>
    void Logger::log3_to_file(
            std::string aOutputSpecifier1, T1 aOutputValue1,
            std::string aOutputSpecifier2, T2 aOutputValue2,
            std::string aOutputSpecifier3, T3 aOutputValue3)
    {
            this->log_to_file(aOutputSpecifier1, aOutputValue1);
            this->log_to_file(aOutputSpecifier2, aOutputValue2);
            this->log_to_file(aOutputSpecifier3, aOutputValue3);
    }

    // -----------------------------------------------------------------------------

    // logging operations for specific types of output texts

    void Logger::log_to_file(std::string aOutputString)
    {
        this->log_to_file( "FreeText", aOutputString);
    }
    void Logger::log_to_file_info(std::string aOutputString)
    {
        this->log_to_file( "InfoText", aOutputString);
    }
    void Logger::log_to_file_debug(std::string aOutputString)
    {
        this->log_to_file( "DebugText", aOutputString);
    }
    void Logger::log_to_file_warning(std::string aOutputString)
    {
        this->log_to_file( "Warning", aOutputString);
    }
    void Logger::log_to_file_error(std::string aOutputString)
    {
        this->log_to_file( "Error", aOutputString);
    }

    // -----------------------------------------------------------------------------
    // FORMATTING TOOLS FOR OUTPUT -------------------------------------------------
    // -----------------------------------------------------------------------------

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

    // -----------------------------------------------------------------------------

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

