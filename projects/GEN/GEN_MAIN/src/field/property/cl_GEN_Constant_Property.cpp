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

        void Constant_Property::get_dfield_dcoordinates(
                const Matrix<DDRMat>& aCoordinates,
                Matrix<DDRMat>&       aSensitivities)
        {
            MORIS_ERROR(false, "get_dfield_dcoordinates not implemented for constant property.");
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
