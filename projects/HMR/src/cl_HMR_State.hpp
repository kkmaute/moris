/*
 * cl_HMR_Runstate.hpp
 *
 *  Created on: Sep 12, 2018
 *      Author: messe
 */

#ifndef PROJECTS_HMR_SRC_CL_HMR_STATE_HPP_
#define PROJECTS_HMR_SRC_CL_HMR_STATE_HPP_

namespace moris
{
    namespace hmr
    {
// -----------------------------------------------------------------------------

        enum class State
        {
            PRINT_USAGE,
            PRINT_HELP,
            PRINT_VERSION,
            INITIALIZE_MESH,
            REFINE_MESH,
            MAP_FIELDS
        };
// -----------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */


#endif /* PROJECTS_HMR_SRC_CL_HMR_STATE_HPP_ */
