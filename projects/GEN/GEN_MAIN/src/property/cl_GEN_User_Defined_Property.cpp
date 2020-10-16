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
                                                     MORIS_GEN_SENSITIVITY_FUNCTION aSensitivityEvaluationFunction,
                                                     sint aNumRefinements,
                                                     sint aRefinementFunctionIndex,
                                                     sint aBSplineMeshIndex,
                                                     real aBSplineLowerBound,
                                                     real aBSplineUpperBound)
                : Field(aADVs,
                        aPropertyVariableIndices,
                        aADVIndices,
                        aConstantParameters,
                        aNumRefinements,
                        aRefinementFunctionIndex,
                        aBSplineMeshIndex,
                        aBSplineLowerBound,
                        aBSplineUpperBound)
        {
            get_field_value_user_defined = aFieldEvaluationFunction;
            get_field_adv_sensitivities_user_defined = aSensitivityEvaluationFunction;
        }

        //--------------------------------------------------------------------------------------------------------------

        real User_Defined_Property::get_field_value(const Matrix<DDRMat>& aCoordinates)
        {
            return this->get_field_value_user_defined(aCoordinates, mFieldVariables);
        }

        //--------------------------------------------------------------------------------------------------------------

        void User_Defined_Property::evaluate_sensitivities(const Matrix<DDRMat>& aCoordinates, Matrix<DDRMat>& aSensitivities)
        {
            this->get_field_adv_sensitivities_user_defined(aCoordinates, mFieldVariables, aSensitivities);
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
