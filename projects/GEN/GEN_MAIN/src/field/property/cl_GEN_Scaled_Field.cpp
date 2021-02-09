#include "cl_GEN_Scaled_Field.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        real Scaled_Field::get_field_value(
                uint aNodeIndex,
                const Matrix<DDRMat>& aCoordinates)
        {
            return *mFieldVariables(0) * mField->get_field_value(aNodeIndex, aCoordinates);
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix<DDRMat>& Scaled_Field::get_field_sensitivities(
                uint aNodeIndex,
                const Matrix<DDRMat>& aCoordinates)
        {
            mSensitivities = *mFieldVariables(0) *
                    mField->get_field_sensitivities(aNodeIndex, aCoordinates);
            return mSensitivities;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDSMat> Scaled_Field::get_determining_adv_ids(
                uint aNodeIndex,
                const Matrix<DDRMat>& aCoordinates)
        {
            return mField->get_determining_adv_ids(aNodeIndex, aCoordinates);
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
