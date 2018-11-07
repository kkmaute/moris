/*
 * cl_MTK_Mapper_Arguments.hpp
 *
 *  Created on: Nov 5, 2018
 *      Author: messe
 */

#ifndef PROJECTS_MTK_MAP_SRC_CL_MTK_MAPPER_ARGUMENTS_HPP_
#define PROJECTS_MTK_MAP_SRC_CL_MTK_MAPPER_ARGUMENTS_HPP_

#include <string>

#include "cl_MTK_Mapper_State.hpp"


namespace moris
{
    namespace mapper
    {
// -----------------------------------------------------------------------------

        class Arguments
        {
            std::string mParameterPath  = "";
            State       mState;

//--------------------------------------------------------------------------------
        public:
//--------------------------------------------------------------------------------

            Arguments(
                        int  & argc,
                        char * argv[] );

//--------------------------------------------------------------------------------

            void
            print_usage();

//---------------------------------------------------------------------------------

            void
            print_help();

//---------------------------------------------------------------------------------

            /**
             * return the run state of the executable
             */
            State
            get_state() const
            {
                return mState;
            }
//---------------------------------------------------------------------------------

            /**
             * return the parameter path
             */
            const std::string &
            get_parameter_path() const
            {
                return mParameterPath;
            }

//---------------------------------------------------------------------------------
        };

// -----------------------------------------------------------------------------
    } /* namespace mapper */
} /* namespace moris */



#endif /* PROJECTS_MTK_MAP_SRC_CL_MTK_MAPPER_ARGUMENTS_HPP_ */
