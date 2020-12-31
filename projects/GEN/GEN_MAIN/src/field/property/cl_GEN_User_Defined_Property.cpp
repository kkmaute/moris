#include "cl_GEN_User_Defined_Property.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        real User_Defined_Property::get_field_value(const Matrix<DDRMat>& aCoordinates)
        {
            return this->get_field_value_user_defined(aCoordinates, mFieldVariables);
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix<DDRMat>& User_Defined_Property::get_field_sensitivities(const Matrix<DDRMat>& aCoordinates)
        {
            this->get_field_sensitivities_user_defined(aCoordinates, mFieldVariables, mSensitivities);
            return mSensitivities;
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
