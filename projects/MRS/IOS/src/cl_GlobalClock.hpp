/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GlobalClock.hpp
 *
 */

#ifndef MORIS_IOS_CL_GLOBALCLOCK_HPP_
#define MORIS_IOS_CL_GLOBALCLOCK_HPP_

#include <cstdio>
#include <string>
#include <cstring>
#include <unordered_map>

// define time functions
#include <ctime>
#include <chrono>

// Define uint, real, etc.
#include "moris_typedefs.hpp"

namespace moris
{

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

        // list of starting times for each iteration
        std::vector< real > mIterationTimeStamps;

        // list of starting times for each active entity
        std::vector< real > mTimeStamps;

        // list of maps for action data for each active entry
        std::vector< std::unordered_map< std::string, real > > mActionData;

        // track number of function IDs
        uint mMaxFunctionID = 0;

        // wall clock timer for debugging purposes
        std::vector< std::chrono::_V2::system_clock::time_point > mWallTimeStamps;

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
                const std::string& aEntityBase,
                const std::string& aEntityType,
                const std::string& aEntityAction );

        // --------------------------------------------------------------------------------
        // operation to stop tracing an entity, decrement all lists
        void sign_out();

        // --------------------------------------------------------------------------------
        // operation to increment iteration count of currently active instance
        void iterate();

        // --------------------------------------------------------------------------------
    };    // class GlobalClock
}    // namespace moris

#endif /* MORIS_IOS_CL_GLOBALCLOCK_HPP_ */

