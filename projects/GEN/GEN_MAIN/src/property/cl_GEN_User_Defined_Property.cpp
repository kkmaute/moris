//
// Created by christopherson on 5/18/20.
//

#include "cl_GEN_User_Defined_Property.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        User_Defined_Property::User_Defined_Property(Matrix<DDRMat>& aADVs,
                                                     Matrix<DDUMat> aPropertyVariableIndices,
                                                     Matrix<DDUMat> aADVIndices,
                                                     Matrix<DDRMat> aConstantParameters,
                                                     Cell<std::shared_ptr<Property>> aPropertyDependencies,
                                                     MORIS_GEN_FIELD_FUNCTION aFieldEvaluationFunction,
                                                     MORIS_GEN_SENSITIVITY_FUNCTION aSensitivityEvaluationFunction)
                                                     : Field(aADVs,
                                                       aPropertyVariableIndices,
                                                       aADVIndices,
                                                       aConstantParameters),
                                                     Property(aPropertyDependencies)
        {
            evaluate_field_value_user_defined = aFieldEvaluationFunction;
            evaluate_sensitivity_user_defined = aSensitivityEvaluationFunction;
        }

        //--------------------------------------------------------------------------------------------------------------

        real User_Defined_Property::evaluate_field_value(const Matrix<DDRMat>& aCoordinates)
        {
            return this->evaluate_field_value_user_defined(aCoordinates, mFieldVariables);
        }

        //--------------------------------------------------------------------------------------------------------------

        void User_Defined_Property::evaluate_all_sensitivities(const Matrix<DDRMat>& aCoordinates, Matrix<DDRMat>& aSensitivities)
        {
            this->evaluate_sensitivity_user_defined(aCoordinates, mFieldVariables, aSensitivities);
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
