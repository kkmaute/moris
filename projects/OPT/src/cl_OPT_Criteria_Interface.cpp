/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_OPT_Criteria_Interface.cpp
 *
 */

#include "cl_OPT_Criteria_Interface.hpp"

namespace moris
{
    namespace opt
    {

        //--------------------------------------------------------------------------------------------------------------

        Vector< real > Criteria_Interface::get_criteria( Vector< real >& aNewADVs )
        {
            mEvaluated = false;

            return this->perform( aNewADVs );
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Criteria_Interface::get_dcriteria_dadv()
        {
            if( !mInitializeOptimizationRestart)
            {
                if (!mEvaluated)
                {
                    mSensitivities = this->compute_dcriteria_dadv();
                    mEvaluated = true;
                }
            }
            return mSensitivities;
        }

        //--------------------------------------------------------------------------------------------------------------
    }
}

