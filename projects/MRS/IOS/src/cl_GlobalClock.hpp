#ifndef MORIS_IOS_CL_GLOBALCLOCK_HPP_
#define MORIS_IOS_CL_GLOBALCLOCK_HPP_

#include <cstdio>
#include <string>
#include <cstring>

// define time functions
#include <ctime>

// Define uint, real, etc.
#include "typedefs.hpp"

namespace moris
{

// Define enums used
enum class EntityBase;
enum class EntityType;
enum class EntityAction;
enum class OutputSpecifier;


class GlobalClock
    {

    //------------------------------------ PRIVATE ------------------------------------

    private:

    friend class Logger;

    // indentation level indicating how far down the tree the trace is
    uint mIndentationLevel = 0;

    // lists of currently active processes
    std::vector< uint > mCurrentFunctionID;

    // lists of currently active entities in tracing tree
    std::vector< std::string > mCurrentEntity;

    // lists of currently active entity types in tracing tree
    std::vector< std::string > mCurrentType;

    // list of currently active action
    std::vector< std::string > mCurrentAction;

    // list of current iteration for each instance
    std::vector< uint > mCurrentIteration;

    // list of starting times for each active entity
    std::vector< real > mTimeStamps;

    // track number of function IDs
    uint mMaxFunctionID = 0;

    //------------------------------------ PUBLIC -------------------------------------

    public:

    // --------------------------------------------------------------------------------
    // constructor initializing all members
    GlobalClock();

    // --------------------------------------------------------------------------------
    // destructor
    ~GlobalClock(){};

    // --------------------------------------------------------------------------------
    // operation to start tracing new entity, increment all lists

    /**
     * Sign in to the clock with an entity action.
     *
     * @param aEntityBase Entity base
     * @param aEntityType Entity type
     * @param aEntityAction Entity action
     */
    void sign_in(
            std::string aEntityBase,
            std::string aEntityType,
            std::string aEntityAction);

    // --------------------------------------------------------------------------------
    // operation to stop tracing an entity, decrement all lists
    void sign_out();

    // --------------------------------------------------------------------------------
    // operation to increment iteration count of currently active instance
    void iterate();

    // --------------------------------------------------------------------------------
    }; // class GlobalClock
} // namespace moris

#endif  /* MORIS_IOS_CL_GLOBALCLOCK_HPP_ */
