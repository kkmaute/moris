#include "cl_GEN_User_Defined_Geometry.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        User_Defined_Geometry::User_Defined_Geometry(Matrix<DDRMat>&                aADVs,
                                                     Matrix<DDUMat>                 aGeometryVariableIndices,
                                                     Matrix<DDUMat>                 aADVIndices,
                                                     Matrix<DDRMat>                 aConstantParameters,
                                                     MORIS_GEN_FIELD_FUNCTION       aFieldEvaluationFunction,
                                                     MORIS_GEN_SENSITIVITY_FUNCTION aSensitivityEvaluationFunction,
                                                     sint                           aNumRefinements,
                                                     sint                           aRefinementFunctionIndex,
                                                     sint                           aBSplineMeshIndex,
                                                     real                           aBSplineLowerBound,
                                                     real                           aBSplineUpperBound)
                : Field(aADVs,
                        aGeometryVariableIndices,
                        aADVIndices,
                        aConstantParameters,
                        aNumRefinements,
                        aRefinementFunctionIndex,
                        aBSplineMeshIndex,
                        aBSplineLowerBound,
                        aBSplineUpperBound)
        {
            this->set_user_defined_functions(aFieldEvaluationFunction, aSensitivityEvaluationFunction);
        }

        //--------------------------------------------------------------------------------------------------------------

        User_Defined_Geometry::User_Defined_Geometry(Matrix<DDRMat>                 aConstantParameters,
                                                     MORIS_GEN_FIELD_FUNCTION       aFieldEvaluationFunction,
                                                     sint                           aNumRefinements,
                                                     sint                           aRefinementFunctionIndex,
                                                     sint                           aBSplineMeshIndex,
                                                     real                           aBSplineLowerBound,
                                                     real                           aBSplineUpperBound)
                : Field(aConstantParameters,
                        aNumRefinements,
                        aRefinementFunctionIndex,
                        aBSplineMeshIndex,
                        aBSplineLowerBound,
                        aBSplineUpperBound)
        {
            this->set_user_defined_functions(aFieldEvaluationFunction, nullptr);
        }

        //--------------------------------------------------------------------------------------------------------------

        real User_Defined_Geometry::evaluate_field_value(const Matrix<DDRMat>& aCoordinates)
        {
            return this->evaluate_field_value_user_defined(aCoordinates, mFieldVariables);
        }

        //--------------------------------------------------------------------------------------------------------------

        void User_Defined_Geometry::evaluate_all_sensitivities(const Matrix<DDRMat>& aCoordinates, Matrix<DDRMat>& aSensitivities)
        {
            this->evaluate_sensitivity_user_defined(aCoordinates, mFieldVariables, aSensitivities);
        }

        //--------------------------------------------------------------------------------------------------------------

        void User_Defined_Geometry::set_user_defined_functions(MORIS_GEN_FIELD_FUNCTION aFieldEvaluationFunction,
                                                               MORIS_GEN_SENSITIVITY_FUNCTION aSensitivityEvaluationFunction)
        {
            // Set field evaluation function
            evaluate_field_value_user_defined = aFieldEvaluationFunction;

            // Check field evaluation function
            MORIS_ERROR(std::abs(this->evaluate_field_value_user_defined({{0.0, 0.0, 0.0}}, mFieldVariables))
                    < std::numeric_limits<real>::infinity(),
                    "There is an error in a user-defined geometry field (field evaluates to infinity).");

            // Set sensitivity evaluation function
            if (aSensitivityEvaluationFunction == nullptr)
            {
                evaluate_sensitivity_user_defined = &(User_Defined_Geometry::no_sensitivities);
            }
            else
            {
                evaluate_sensitivity_user_defined = aSensitivityEvaluationFunction;

                // Check sensitivity function
                Matrix<DDRMat> tSensitivities;
                this->evaluate_sensitivity_user_defined({{0.0, 0.0, 0.0}}, mFieldVariables, tSensitivities);
                MORIS_ERROR(tSensitivities.n_rows() == 1,
                        "A user-defined geometry must provide a row vector for sensitivities.");
                MORIS_ERROR(tSensitivities.n_cols() == mFieldVariables.size(),
                        "A user-defined geometry must have a sensitivity vector with a length equal to the total "
                        "number of geometry variables (ADVs + constant parameters).");
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void User_Defined_Geometry::no_sensitivities(const Matrix<DDRMat>&  aCoordinates,
                                                     const Cell<real*>&     aParameters,
                                                     Matrix<DDRMat>&        aSensitivities)
        {
            MORIS_ERROR(false, "A sensitivity evaluation function was not provided to a user-defined geometry. "
                               "Please make sure that you provide this function, or that sensitivities are not required.");
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
