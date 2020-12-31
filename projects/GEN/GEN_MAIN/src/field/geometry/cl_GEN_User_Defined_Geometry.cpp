#include "cl_GEN_User_Defined_Geometry.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        User_Defined_Geometry::User_Defined_Geometry(
                Matrix<DDRMat>           aConstants,
                MORIS_GEN_FIELD_FUNCTION aFieldEvaluationFunction,
                Field_Parameters         aParameters)
                : Field(aConstants, aParameters)
        {
            this->set_user_defined_functions(aFieldEvaluationFunction, nullptr);
        }

        //--------------------------------------------------------------------------------------------------------------

        real User_Defined_Geometry::get_field_value(const Matrix<DDRMat>& aCoordinates)
        {
            return this->get_field_value_user_defined(aCoordinates, mFieldVariables);
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix<DDRMat>& User_Defined_Geometry::get_field_sensitivities(const Matrix<DDRMat>& aCoordinates)
        {
            this->get_field_sensitivities_user_defined(aCoordinates, mFieldVariables, mSensitivities);
            return mSensitivities;
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
                // Sensitivity function was not provided
                get_field_sensitivities_user_defined = &(User_Defined_Geometry::no_sensitivities);
            }
            else
            {
                // Sensitivity function was provided
                get_field_sensitivities_user_defined = aSensitivityEvaluationFunction;

                // Check sensitivity function
                this->get_field_sensitivities_user_defined({{0.0, 0.0, 0.0}}, mFieldVariables, mSensitivities);

                // Check for row vector
                MORIS_ERROR(mSensitivities.n_rows() == 1,
                        "A user-defined geometry must provide a row vector for sensitivities.");

                // Check for size
                MORIS_ERROR(mSensitivities.n_cols() == mFieldVariables.size(),
                        "A user-defined geometry must have a sensitivity vector with a length equal to the total "
                        "number of geometry variables (ADVs + constants).");

                // Check for values not nan/infinity
                for (uint tSensitivityIndex = 0; tSensitivityIndex < mSensitivities.n_cols(); tSensitivityIndex++)
                {
                    MORIS_ERROR(std::isfinite(mSensitivities(tSensitivityIndex)),
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
