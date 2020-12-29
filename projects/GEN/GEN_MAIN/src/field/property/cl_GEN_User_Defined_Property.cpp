#include "cl_GEN_User_Defined_Property.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        User_Defined_Property::User_Defined_Property(
                Matrix<DDRMat>&                aADVs,
                Matrix<DDUMat>                 aPropertyVariableIndices,
                Matrix<DDUMat>                 aADVIndices,
                Matrix<DDRMat>                 aConstants,
                Cell<std::shared_ptr<Field>>   aFieldDependencies,
                MORIS_GEN_FIELD_FUNCTION       aFieldEvaluationFunction,
                MORIS_GEN_SENSITIVITY_FUNCTION aSensitivityEvaluationFunction,
                Field_Parameters               aParameters)
                : Field(aADVs, aPropertyVariableIndices, aADVIndices, aConstants, aParameters)
        {
            get_field_value_user_defined = aFieldEvaluationFunction;
            get_field_sensitivities_user_defined = aSensitivityEvaluationFunction;
        }

        //--------------------------------------------------------------------------------------------------------------

        User_Defined_Property::User_Defined_Property(
                sol::Dist_Vector*              aOwnedADVs,
                Matrix<DDUMat>                 aPropertyVariableIndices,
                Matrix<DDUMat>                 aADVIndices,
                Matrix<DDRMat>                 aConstants,
                Cell<std::shared_ptr<Field>>   aFieldDependencies,
                MORIS_GEN_FIELD_FUNCTION       aFieldEvaluationFunction,
                MORIS_GEN_SENSITIVITY_FUNCTION aSensitivityEvaluationFunction,
                Field_Parameters               aParameters)
                : Field(aOwnedADVs, aPropertyVariableIndices, aADVIndices, aConstants, aParameters)
        {
            get_field_value_user_defined = aFieldEvaluationFunction;
            get_field_sensitivities_user_defined = aSensitivityEvaluationFunction;
        }

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
