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
                Matrix<DDRMat>                 aConstantParameters,
                Cell<std::shared_ptr<Field>>   aFieldDependencies,
                MORIS_GEN_FIELD_FUNCTION       aFieldEvaluationFunction,
                MORIS_GEN_SENSITIVITY_FUNCTION aSensitivityEvaluationFunction,
                std::string                    aName,
                Matrix<DDSMat>  aNumRefinements,
                Matrix<DDSMat>  aRefinementMeshIndices,
                sint                           aRefinementFunctionIndex,
                sint                           aBSplineMeshIndex,
                real                           aBSplineLowerBound,
                real                           aBSplineUpperBound)
                : Field(aADVs,
                        aPropertyVariableIndices,
                        aADVIndices,
                        aConstantParameters,
                        aName,
                        aNumRefinements,
                        aRefinementMeshIndices,
                        aRefinementFunctionIndex,
                        aBSplineMeshIndex,
                        aBSplineLowerBound,
                        aBSplineUpperBound)
        {
            get_field_value_user_defined = aFieldEvaluationFunction;
            get_field_sensitivities_user_defined = aSensitivityEvaluationFunction;
        }

        //--------------------------------------------------------------------------------------------------------------

        User_Defined_Property::User_Defined_Property(
                sol::Dist_Vector*              aOwnedADVs,
                Matrix<DDUMat>                 aPropertyVariableIndices,
                Matrix<DDUMat>                 aADVIndices,
                Matrix<DDRMat>                 aConstantParameters,
                Cell<std::shared_ptr<Field>>   aFieldDependencies,
                MORIS_GEN_FIELD_FUNCTION       aFieldEvaluationFunction,
                MORIS_GEN_SENSITIVITY_FUNCTION aSensitivityEvaluationFunction,
                std::string                    aName,
                Matrix<DDSMat>  aNumRefinements,
                Matrix<DDSMat>  aRefinementMeshIndices,
                sint                           aRefinementFunctionIndex,
                sint                           aBSplineMeshIndex,
                real                           aBSplineLowerBound,
                real                           aBSplineUpperBound)
                : Field(aOwnedADVs,
                        aPropertyVariableIndices,
                        aADVIndices,
                        aConstantParameters,
                        aName,
                        aNumRefinements,
                        aRefinementMeshIndices,
                        aRefinementFunctionIndex,
                        aBSplineMeshIndex,
                        aBSplineLowerBound,
                        aBSplineUpperBound)
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

        Matrix<DDRMat> User_Defined_Property::get_field_sensitivities(const Matrix<DDRMat>& aCoordinates)
        {
            Matrix<DDRMat> tSensitivities;
            this->get_field_sensitivities_user_defined(aCoordinates, mFieldVariables, tSensitivities);
            return tSensitivities;
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
