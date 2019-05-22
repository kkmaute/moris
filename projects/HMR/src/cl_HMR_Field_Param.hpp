/*
 * cl_HMR_Field_Param.hpp
 *
 *  Created on: Dec 5, 2018
 *      Author: messe
 */

#ifndef PROJECTS_HMR_SRC_CL_HMR_FIELD_PARAM_HPP_
#define PROJECTS_HMR_SRC_CL_HMR_FIELD_PARAM_HPP_

#include "typedefs.hpp"

namespace moris
{
    namespace hmr
    {
// -----------------------------------------------------------------------------

        struct Field_Param
        {
            std::string mLabel;
            moris_id    mID = gNoID;
            uint        mInputLagrangeOrder = 0;
            uint        mInputBSplineOrder = 0;
            uint        mOutputBSplineOrder = 0;
            std::string mSource;
            std::string mTarget;
            real        mL2alpha = 0.0;

            Field_Param(){};

            ~Field_Param(){};
        };

// -----------------------------------------------------------------------------
    }
}

#endif /* PROJECTS_HMR_SRC_CL_HMR_FIELD_PARAM_HPP_ */
