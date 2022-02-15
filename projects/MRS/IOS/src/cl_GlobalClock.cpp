#include "cl_GlobalClock.hpp"
#include "Log_Constants.hpp"

#include <cstdio>
#include <string>
#include <cstring>
#include <iostream>

// define time functions
#include <ctime>

// Define uint, real, etc.
#include "typedefs.hpp"

namespace moris
{

    // --------------------------------------------------------------------------------

    GlobalClock::GlobalClock()
    {
        // initialize list of function IDs
        mCurrentFunctionID.resize( 1, mMaxFunctionID );

        // initialize list of entities
        mCurrentEntity.resize( 1, "GlobalClock" );

        // initialize list of entity types
        mCurrentType.resize( 1, LOGGER_NON_SPECIFIC_ENTITY_TYPE );

        // initialize list of Actions
        mCurrentAction.resize( 1, "Create" );

        // initialize list of Actions
        mCurrentIteration.resize( 1, 0 );

        // record starting time of iterations
        mIterationTimeStamps.resize( 1, (real)std::clock() );

        // record starting time
        mTimeStamps.resize( 1, (real)std::clock() );

        // record action data for new entry
        mActionData.resize( 1, std::unordered_map< std::string, real >() );

        // record starting wall clock time
        if ( PRINT_WALL_TIME )
            mWallTimeStamps.resize( 1, std::chrono::system_clock::now() );
    }

    // --------------------------------------------------------------------------------

    void
    GlobalClock::sign_in(
            std::string aEntityBase,
            std::string aEntityType,
            std::string aEntityAction )
    {
        // increment indentation level
        mIndentationLevel++;

        // create new function ID
        mMaxFunctionID++;

        // assign function ID
        mCurrentFunctionID.push_back( mMaxFunctionID );

        // save entity to list of active entities
        mCurrentEntity.push_back( aEntityBase );

        // save entity to list of active entities
        mCurrentType.push_back( aEntityType );

        // save entity action to list of actions
        mCurrentAction.push_back( aEntityAction );

        // initiate iteration counter with zero
        mCurrentIteration.push_back( 0 );

        // create time stamp for entity
        mIterationTimeStamps.push_back( (real)std::clock() );

        // create time stamp for new entity
        mTimeStamps.push_back( (real)std::clock() );

        // create map for action data for new entry
        mActionData.push_back( std::unordered_map< std::string, real >() );

        // create wall clock time stamp for new entity
        if ( PRINT_WALL_TIME )
            mWallTimeStamps.push_back( std::chrono::system_clock::now() );

#ifdef DEBUG
        // check that indentation level and array size match
        if ( mIndentationLevel != mCurrentFunctionID.size() - 1 )
        {
            std::cout << "GlobalClock::sign_in - indentation level and array sizes do not match.\n";
            throw;
        }
#endif
    }

    // --------------------------------------------------------------------------------

    void
    GlobalClock::sign_out()
    {
        // remove time stamp from list of active entities
        mTimeStamps.pop_back();

        // remove iteration time stamp from list of active entities
        mIterationTimeStamps.pop_back();

        // remove function ID
        mCurrentFunctionID.pop_back();

        // remove entity from list of active entities
        mCurrentEntity.pop_back();

        // remove type from list of active entity types
        mCurrentType.pop_back();

        // remove entity action from list of actions
        mCurrentAction.pop_back();

        // remove iteration counter from list of actions
        mCurrentIteration.pop_back();

        // remove map for action data for new entry
        mActionData.pop_back();

        // decrement indentation level
        mIndentationLevel--;

        // remove wall time stamp from list of active entities
        if ( PRINT_WALL_TIME )
            mWallTimeStamps.pop_back();

#ifdef DEBUG
        // check that indentation level and array size match
        if ( mIndentationLevel != mCurrentFunctionID.size() - 1 )
        {
            std::cout << "GlobalClock::sign_in - indentation level and array sizes do not match.\n";
            throw;
        }
#endif
    }

    // --------------------------------------------------------------------------------

    void
    GlobalClock::iterate()
    {
        // increment iteration counter of currently active action
        mCurrentIteration[mIndentationLevel]++;

        // renew time stamp at beginning of an iteration
        mIterationTimeStamps[mIndentationLevel] = (real)std::clock();
    }

    // --------------------------------------------------------------------------------
}    // namespace moris
