/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Arguments.hpp
 *
 */

#pragma once
#include <string>

#include "cl_HMR_State.hpp"

namespace moris::hmr
{
//--------------------------------------------------------------------------------
    class Arguments
    {
        std::string mParameterPath  = "";
        State       mState;
        double      mTimestep = 0.0;
        bool        mMapWhileRefine = false;

//--------------------------------------------------------------------------------
    public:
//--------------------------------------------------------------------------------

        Arguments( int  & argc,
                   char * argv[] );

//---------------------------------------------------------------------------------

        void print_usage();

//---------------------------------------------------------------------------------

        void print_help();

//---------------------------------------------------------------------------------

        /**
         * return the run state of the executable
         */
        State get_state() const
        {
            return mState;
        }

//---------------------------------------------------------------------------------

        /**
         * return the parameter path
         */
        const std::string & get_parameter_path() const
        {
            return mParameterPath;
        }

//---------------------------------------------------------------------------------

        /**
         * return the timestep that is written into the exodus file
         */
        double get_timestep() const
        {
            return mTimestep;
        }

//---------------------------------------------------------------------------------

        /**
         * returns true if user wants to call refinement and mapping at the same time
         */
        bool map_while_refine() const
        {
            return mMapWhileRefine;
        }
    };
//---------------------------------------------------------------------------------
}
