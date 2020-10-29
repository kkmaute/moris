#include "cl_GEN_User_Defined_Geometry.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        User_Defined_Geometry::User_Defined_Geometry(
                Matrix<DDRMat>&                aADVs,
                Matrix<DDUMat>                 aGeometryVariableIndices,
                Matrix<DDUMat>                 aADVIndices,
                Matrix<DDRMat>                 aConstantParameters,
                MORIS_GEN_FIELD_FUNCTION       aFieldEvaluationFunction,
                MORIS_GEN_SENSITIVITY_FUNCTION aSensitivityEvaluationFunction,
                std::string                    aName,
                Matrix<DDSMat>  aNumRefinements,
                Matrix<DDSMat>  aNumPatterns,
                sint                           aRefinementFunctionIndex,
                sint                           aBSplineMeshIndex,
                real                           aBSplineLowerBound,
                real                           aBSplineUpperBound)
                : Field(aADVs,
                        aGeometryVariableIndices,
                        aADVIndices,
                        aConstantParameters,
                        aName,
                        aNumRefinements,
                        aNumPatterns,
                        aRefinementFunctionIndex,
                        aBSplineMeshIndex,
                        aBSplineLowerBound,
                        aBSplineUpperBound)
        {
            this->set_user_defined_functions(aFieldEvaluationFunction, aSensitivityEvaluationFunction);
        }

        //--------------------------------------------------------------------------------------------------------------

        User_Defined_Geometry::User_Defined_Geometry(
                sol::Dist_Vector* aOwnedADVs,
                Matrix<DDUMat>                 aGeometryVariableIndices,
                Matrix<DDUMat>                 aADVIndices,
                Matrix<DDRMat>                 aConstantParameters,
                MORIS_GEN_FIELD_FUNCTION       aFieldEvaluationFunction,
                MORIS_GEN_SENSITIVITY_FUNCTION aSensitivityEvaluationFunction,
                std::string                    aName,
                Matrix<DDSMat>  aNumRefinements,
                Matrix<DDSMat>  aNumPatterns,
                sint                           aRefinementFunctionIndex,
                sint                           aBSplineMeshIndex,
                real                           aBSplineLowerBound,
                real                           aBSplineUpperBound)
                : Field(aOwnedADVs,
                        aGeometryVariableIndices,
                        aADVIndices,
                        aConstantParameters,
                        aName,
                        aNumRefinements,
                        aNumPatterns,
                        aRefinementFunctionIndex,
                        aBSplineMeshIndex,
                        aBSplineLowerBound,
                        aBSplineUpperBound)
        {
            this->set_user_defined_functions(aFieldEvaluationFunction, aSensitivityEvaluationFunction);
        }

        //--------------------------------------------------------------------------------------------------------------

        User_Defined_Geometry::User_Defined_Geometry(
                Matrix<DDRMat>           aConstantParameters,
                MORIS_GEN_FIELD_FUNCTION aFieldEvaluationFunction,
                std::string              aName,
                Matrix<DDSMat>  aNumRefinements,
                Matrix<DDSMat>  aNumPatterns,
                sint                     aRefinementFunctionIndex,
                sint                     aBSplineMeshIndex,
                real                     aBSplineLowerBound,
                real                     aBSplineUpperBound)
                : Field(aConstantParameters,
                        aName,
                        aNumRefinements,
                        aNumPatterns,
                        aRefinementFunctionIndex,
                        aBSplineMeshIndex,
                        aBSplineLowerBound,
                        aBSplineUpperBound)
        {
            this->set_user_defined_functions(aFieldEvaluationFunction, nullptr);
        }

        //--------------------------------------------------------------------------------------------------------------

        real User_Defined_Geometry::get_field_value(const Matrix<DDRMat>& aCoordinates)
        {
            return this->get_field_value_user_defined(aCoordinates, mFieldVariables);
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> User_Defined_Geometry::get_field_sensitivities(const Matrix<DDRMat>& aCoordinates)
        {
            Matrix<DDRMat> tSensitivities(0, 0);
            this->get_field_sensitivities_user_defined(aCoordinates, mFieldVariables, tSensitivities);
            return tSensitivities;
        }

        //--------------------------------------------------------------------------------------------------------------

        void User_Defined_Geometry::set_user_defined_functions(
                MORIS_GEN_FIELD_FUNCTION aFieldEvaluationFunction,
                MORIS_GEN_SENSITIVITY_FUNCTION aSensitivityEvaluationFunction)
        {
            // Set field evaluation function
            get_field_value_user_defined = aFieldEvaluationFunction;

            // Check field evaluation function
            MORIS_ERROR(std::isfinite(this->get_field_value_user_defined({{0.0, 0.0, 0.0}}, mFieldVariables)),
                    "There is an error in a user-defined geometry field (field evaluates to nan/infinity).");

            // Set sensitivity evaluation function
            if (aSensitivityEvaluationFunction == nullptr)
            {
                get_field_sensitivities_user_defined = &(User_Defined_Geometry::no_sensitivities);
            }
            else
            {
                get_field_sensitivities_user_defined = aSensitivityEvaluationFunction;

                // Check sensitivity function
                Matrix<DDRMat> tSensitivities(0, 0);
                this->get_field_sensitivities_user_defined({{0.0, 0.0, 0.0}}, mFieldVariables, tSensitivities);

                // Check for row vector
                MORIS_ERROR(tSensitivities.n_rows() == 1,
                        "A user-defined geometry must provide a row vector for sensitivities.");

                // Check for size
                MORIS_ERROR(tSensitivities.n_cols() == mFieldVariables.size(),
                        "A user-defined geometry must have a sensitivity vector with a length equal to the total "
                        "number of geometry variables (ADVs + constant parameters).");

                // Check for values not nan/infinity
                for (uint tSensitivityIndex = 0; tSensitivityIndex < tSensitivities.n_cols(); tSensitivityIndex++)
                {
                    MORIS_ERROR(std::isfinite(tSensitivities(tSensitivityIndex)),
                            "There is an error in a user-defined geometry sensitivity (evaluates to nan/infinity).");
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void User_Defined_Geometry::no_sensitivities(
                const Matrix<DDRMat>& aCoordinates,
                const Cell<real*>&    aParameters,
                Matrix<DDRMat>&       aSensitivities)
        {
            MORIS_ERROR(false, "A sensitivity evaluation function was not provided to a user-defined geometry. "
                               "Please make sure that you provide this function, or that sensitivities are not required.");
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
