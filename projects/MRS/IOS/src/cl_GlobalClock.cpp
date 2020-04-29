#include "cl_GlobalClock.hpp"

#include <cstdio>
#include <string>
#include <cstring>

// define time functions
#include <ctime>

// Define Cells
#include "cl_Cell.hpp"

// Define uint, real, etc.
#include "typedefs.hpp"

// Define enums used
#include "cl_Tracer_Enums.hpp"


namespace moris
{

    // -------------------------------------------------------------------------------- //
    // constructor initializing all members
    GlobalClock::GlobalClock()
    {
        // initialize list of function IDs
        mCurrentFunctionID.resize( 1 , mMaxFunctionID );
//        mCurrentFunctionID = new Cell<uint> ( 1 , mMaxFunctionID );
        // ^^ do this for all

        // initialize list of entities
        mCurrentEntity.resize( 1 , EntityBase::GlobalClock );

        // initialize list of entity types
        mCurrentType.resize( 1 , EntityType::Base );

        // initialize list of Actions
        mCurrentAction.resize( 1 , EntityAction::Create );

        // record starting time
        mTimeStamps.resize( 1 , (real) std::clock() );
    }


    // -------------------------------------------------------------------------------- //
    // destructor
    GlobalClock::~GlobalClock()
    {

//        // log last line
//        if (mIndentationLevel != 0)
//        {
//            // set indentation level to zero
//            mIndentationLevel = 0;
//
//            // report error
//            std::cout << "Global Clock prematurely stopped.";
//        }
//        else
//        {
//            // do nothing
//        }

    }

    // -------------------------------------------------------------------------------- //
    // operation to start tracing new entity, increment all lists
    void GlobalClock::sign_in( enum EntityBase aEntityBase, enum EntityType aEntityType, enum EntityAction aEntityAction )
    {
        // increment indentation level
        mIndentationLevel++;

        // create new function ID
        mMaxFunctionID++;

        // assign function ID
        mCurrentFunctionID.push_back(mMaxFunctionID);

        // save entity to list of active entities
        mCurrentEntity.push_back(aEntityBase);

        // save entity to list of active entities
        mCurrentType.push_back(aEntityType);

        // save entity action to list of actions
        mCurrentAction.push_back(aEntityAction);

        // create time stamp for new entity
        mTimeStamps.push_back( (real) std::clock() );
    }


    // -------------------------------------------------------------------------------- //
    // operation to stop tracing an entity, decrement all lists
    void GlobalClock::sign_out()
    {
        // remove time stamp from list of active entities
        mTimeStamps.pop_back();

        // remove function ID
        mCurrentFunctionID.pop_back();

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

} // namespace moris


