/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SDF_Arguments.hpp
 *
 */

#ifndef PROJECTS_GEN_SDF_SRC_CL_SDF_ARGUMENTS_HPP_
#define PROJECTS_GEN_SDF_SRC_CL_SDF_ARGUMENTS_HPP_

#include <string>
#include "cl_SDF_State.hpp"

namespace moris::sdf
{
    class Arguments
    {
        std::string mParameterPath  = "";
        std::string mInputMeshPath  = "";
        std::string mOutputMeshPath = "";
        double      mTimestep       = 0.0;
        State       mState;
        //--------------------------------------------------------------------------------

      public:
        //--------------------------------------------------------------------------------

        Arguments(
                int  &argc,
                char *argv[] );

        //--------------------------------------------------------------------------------

        void
        print_usage();

        //--------------------------------------------------------------------------------

        void
        print_help();

        //--------------------------------------------------------------------------------

        const State &
        get_state() const
        {
            return mState;
        }

        //--------------------------------------------------------------------------------

        const std::string &
        get_parameter_path() const
        {
            return mParameterPath;
        }

        //--------------------------------------------------------------------------------

        const std::string &
        get_input_mesh_path() const
        {
            return mInputMeshPath;
        }

        //--------------------------------------------------------------------------------

        const std::string &
        get_output_mesh_path() const
        {
            return mOutputMeshPath;
        }

        //--------------------------------------------------------------------------------

        double
        get_timestep() const
        {
            return mTimestep;
        }
    };
    }

#endif /* PROJECTS_GEN_SDF_SRC_CL_SDF_ARGUMENTS_HPP_ */

