#ifndef MORIS_IOS_CL_GLOBALCLOCK_HPP_
#define MORIS_IOS_CL_GLOBALCLOCK_HPP_

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <cstring>

// define time functions
#include <ctime>

// Define Cells
#include "cl_Cell.hpp"

// Define uint, real, etc.
#include "typedefs.hpp"

// needed?
#include "IO_Tools.hpp"

// Define enums used
#include "cl_Tracer_Enums.hpp"

// Get access to processor count, Ids, etc.
#include "cl_Communication_Tools.hpp" // COM/src


namespace moris
{

class GlobalClock
    {

    //------------------------------------ PRIVATE -------------------------------------//
    private:

    // indentation level indicating how far down the tree the trace is
    uint mIndentationLevel = 0;

    // lists of currently active entities in tracing tree
    Cell< enum EntityBase > mCurrentEntity;

    // lists of currently active entity types in tracing tree
    Cell< enum EntityType > mCurrentType;

    // list of currently active action
    Cell< enum EntityAction > mCurrentAction;

    // list of starting times for each active entity
    Cell< real > mTimeStamps;

    // text file for outputting info
    std::ofstream mGlobalLog;

    //------------------------------------ PUPLIC --------------------------------------//
    public:


    // -------------------------------------------------------------------------------- //
    // constructor
    GlobalClock()
    {
        // initialize list of entities
        mCurrentEntity.resize( 1 , EntityBase::GlobalClock );

        // initialize list of entity types
        mCurrentType.resize( 1 , EntityType::Base );

        // initialize list of Actions
        mCurrentAction.resize( 1 , EntityAction::Create );

        // record starting time
        mTimeStamps.resize( 1 , (real) std::clock() );

        // open text file
        mGlobalLog.open("gClock.log");

        // create header
        this->print_header();
        //mGlobalLog << "This is a header. \n Things like date of run, machine, code version, etc. go here. \n The following lines contain all information. \n \n";
        //mGlobalLog << "Indentation Level ; Entity Type ; Entity Action | Output Specifier : Output Value | Output Specifier : Output Value | ... \n";

        // sign in global clock
        this->log(OutputSpecifier::Signing, 1.0);

    };


    // -------------------------------------------------------------------------------- //
    // destructor
    ~GlobalClock()
    {

        // log last line
        if (mIndentationLevel != 0)
        {
            // set indentation level to zero
            mIndentationLevel = 0;

            // sign out clock and report error
            this->log3(OutputSpecifier::Signing, 0.0,
                       OutputSpecifier::Time, ( (moris::real) std::clock() - mTimeStamps(mIndentationLevel) ) / CLOCKS_PER_SEC ,
                       OutputSpecifier::Error, 1.0);

        }
        else
        {
            // sign out clock
            this->log2(OutputSpecifier::Signing, 0.0,
                       OutputSpecifier::Time, ( (moris::real) std::clock() - mTimeStamps(mIndentationLevel) ) / CLOCKS_PER_SEC );
        }

        // footer
        mGlobalLog << "CLOCK STOPPED.";

        // close text file
        mGlobalLog.close();
    };


    // -------------------------------------------------------------------------------- //
    // logging operation for one output value
    void log(enum OutputSpecifier aOutputSpecifier, real aOutputValue)
    {
        mGlobalLog << mIndentationLevel << "; "
                  << get_enum_str(mCurrentEntity(mIndentationLevel)) << "; "
                  << get_enum_str(mCurrentType(mIndentationLevel)) << "; "
                  << get_enum_str(mCurrentAction(mIndentationLevel)) << " | "
                  << get_enum_str(aOutputSpecifier) << "; "
                  << aOutputValue << " | \n";
    }


    // logging operation for two output values
    void log2(enum OutputSpecifier aOutputSpecifier1, real aOutputValue1,
              enum OutputSpecifier aOutputSpecifier2, real aOutputValue2)
    {
        mGlobalLog << mIndentationLevel << "; "
                  << get_enum_str(mCurrentEntity(mIndentationLevel)) << "; "
                  << get_enum_str(mCurrentType(mIndentationLevel)) << "; "
                  << get_enum_str(mCurrentAction(mIndentationLevel)) << " | "
                  << get_enum_str(aOutputSpecifier1) << "; "
                  << aOutputValue1 << "; "
                  << get_enum_str(aOutputSpecifier2) << "; "
                  << aOutputValue2 << "| \n";
    }


    // logging operation for three output values
    void log3(enum OutputSpecifier aOutputSpecifier1, real aOutputValue1,
              enum OutputSpecifier aOutputSpecifier2, real aOutputValue2,
              enum OutputSpecifier aOutputSpecifier3, real aOutputValue3)
    {
        mGlobalLog << mIndentationLevel << "; "
                  << get_enum_str(mCurrentEntity(mIndentationLevel)) << "; "
                  << get_enum_str(mCurrentType(mIndentationLevel)) << "; "
                  << get_enum_str(mCurrentAction(mIndentationLevel)) << " | "
                  << get_enum_str(aOutputSpecifier1) << "; "
                  << aOutputValue1 << "; "
                  << get_enum_str(aOutputSpecifier2) << "; "
                  << aOutputValue2 << "; "
                  << get_enum_str(aOutputSpecifier3) << "; "
                  << aOutputValue3 << " | \n";
    }

    // -------------------------------------------------------------------------------- //
    // operation to start tracing new entity
    void sign_in( enum EntityBase aEntityBase, enum EntityType aEntityType, enum EntityAction aEntityAction )
    {
        // increment indentation level
        mIndentationLevel++;

        // save entity to list of active entities
        mCurrentEntity.push_back(aEntityBase);

        // save entity to list of active entities
        mCurrentType.push_back(aEntityType);

        // save entity action to list of actions
        mCurrentAction.push_back(aEntityAction);

        // log current position in code
        this->log(OutputSpecifier::Signing, 1.0);

        // create time stamp for new entity
        mTimeStamps.push_back( (real) std::clock() );
    }


    // -------------------------------------------------------------------------------- //
    // operation to stop tracing an entity
    void sign_out()
    {
        // log current position in code
        this->log2(OutputSpecifier::Signing, 0.0,
                   OutputSpecifier::Time, ( (moris::real) std::clock() - mTimeStamps(mIndentationLevel) ) / CLOCKS_PER_SEC );

        // remove time stamp from list of active entities
        mTimeStamps.pop_back();

        // remove entity from list of active entities
        mCurrentEntity.pop_back();

        // remove type from list of active entity types
        mCurrentType.pop_back();

        // save entity action to list of actions
        mCurrentAction.pop_back();

        // decrement indentation level
        mIndentationLevel--;
    }

    // -------------------------------------------------------------------------------- //
    // operation to stop global clock
    void stop()
    {
        this->~GlobalClock();
    }

    // -------------------------------------------------------------------------------- //
    // print header
    void print_header()
    {

        mGlobalLog << "Output from global Clock. >>> \n";

        // get date and time at runtime
        time_t tTimeStamp = time(NULL);

        // print info about test
        mGlobalLog << "\n--- RUN ---\n";
//        mGlobalLog << "Executable: " << argv[ 0 ] << " \n";
//        mGlobalLog << "Arguments: ";
//        for (int ia = 1; ia < argc; ++ia)
//        {
//            mGlobalLog << argv[ ia ] << "; ";
//        }
//        mGlobalLog << "Number of Processors used: " << ( int ) par_size() << " \n";
        mGlobalLog << "Date of execution: " << ctime(&tTimeStamp) << " \n";


        // print user information
        mGlobalLog << "\n--- HOST ---\n";
        mGlobalLog << "User: " << std::getenv( "USER" ) << "\n";
        mGlobalLog << "Host: " << std::getenv( "HOSTNAME" ) << "\n";
        mGlobalLog << "Operating System: " << std::getenv( "OSTYPE" ) << "\n";
        mGlobalLog << "Host Type: " << std::getenv( "HOSTTYPE" ) << "\n";


        // print build info
        mGlobalLog << "\n--- BUILD ---\n";
        mGlobalLog << "Build date: " << __DATE__ << "; " << __TIME__ << "\n";
        mGlobalLog << "DEBUG: ";

#if !defined(NDEBUG) || defined(DEBUG)
        mGlobalLog << "ON\n";
#endif
#if defined(NDEBUG) || !defined(DEBUG)
        mGlobalLog << "OFF\n";
#endif

        mGlobalLog << "Matrix Library: ";

#ifdef MORIS_USE_ARMA
                mGlobalLog << "ARMA \n";
#endif
#ifdef MORIS_USE_EIGEN
        mGlobalLog << "Eigen \n";
#endif


        // Top row in table
        mGlobalLog << "\n____________________________________________________________________________________________________________________________________\n";
        mGlobalLog << "Indentation Level ; Entity ; Entity Type ; Entity Action | Output Specifier ; Output Value ; Output Specifier ; Output Value ; ... |\n";
        mGlobalLog << "------------------------------------------------------------------------------------------------------------------------------------\n";

    }

    // -------------------------------------------------------------------------------- //
    }; // class GlobalClock
} // namespace moris


#endif  /* MORIS_IOS_CL_GLOBALCLOCK_HPP_ */

