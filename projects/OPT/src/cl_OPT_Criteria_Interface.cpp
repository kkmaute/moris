#include "cl_OPT_Criteria_Interface.hpp"

namespace moris
{
    namespace opt
    {

        //--------------------------------------------------------------------------------------------------------------

        Matrix< DDRMat >
        Criteria_Interface::get_criteria( Matrix< DDRMat >& aNewADVs )
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
