#include "cl_GEN_Constant_Property.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        real Constant_Property::get_field_value(
                uint                  aNodeIndex,
                const Matrix<DDRMat>& aCoordinates)
        {
            return *mFieldVariables(0);
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix<DDRMat>& Constant_Property::get_field_sensitivities(
                uint                  aNodeIndex,
                const Matrix<DDRMat>& aCoordinates)
        {
            mSensitivities = {{1.0}};
            return mSensitivities;
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
