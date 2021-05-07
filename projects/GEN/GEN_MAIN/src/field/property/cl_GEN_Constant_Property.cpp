#include "cl_GEN_Constant_Property.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        real Constant_Property::get_field_value(const Matrix<DDRMat>& aCoordinates)
        {
            return *mFieldVariables(0);
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix<DDRMat>& Constant_Property::get_dfield_dadvs(const Matrix<DDRMat>& aCoordinates)
        {
            return mSensitivities;
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
