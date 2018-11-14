/*
 * cl_SDF_State.hpp
 *
 *  Created on: Oct 15, 2018
 *      Author: messe
 */

#ifndef PROJECTS_GEN_SDF_SRC_CL_SDF_STATE_HPP_
#define PROJECTS_GEN_SDF_SRC_CL_SDF_STATE_HPP_

namespace moris
{
    namespace sdf
    {
//-------------------------------------------------------------------------------

        enum class State
        {
            PRINT_USAGE,
            PRINT_HELP,
            PRINT_VERSION,
            CALCULATE_RAYCAST,
            CALCULATE_RAYCAST_AND_SDF
        };

//-------------------------------------------------------------------------------
    }
}



#endif /* PROJECTS_GEN_SDF_SRC_CL_SDF_STATE_HPP_ */
