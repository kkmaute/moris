#ifndef MORIS_CL_GEN_USER_DEFINED_PROPERTY_HPP
#define MORIS_CL_GEN_USER_DEFINED_PROPERTY_HPP

#include "cl_GEN_Property.hpp"
#include "cl_GEN_Field_Analytic.hpp"
#include "fn_Exec_load_user_library.hpp"

namespace moris
{
    namespace ge
    {
        class User_Defined_Property : public Property, public Field_Analytic
        {
        private:
            MORIS_GEN_FIELD_FUNCTION get_field_value_user_defined;
            MORIS_GEN_SENSITIVITY_FUNCTION get_field_sensitivities_user_defined;

        public:

            /**
             * Constructor
             *
             * @param aADVs Reference to the full advs
             * @param aPropertyVariableIndices Indices of property variables to be filled by the ADVs
             * @param aADVIndices The indices of the ADV vector to fill in the property variables
             * @param aConstants The constant field variables not filled by ADVs
             * @param aFieldDependencies Other created fields that this property depends on
             * @param aFieldEvaluationFunction User-defined function for evaluating the property field
             * @param aSensitivityEvaluationFunction User-defined function for evaluating the field sensitivities
             * @param aParameters Additional parameters
             */
            User_Defined_Property(
                    Matrix<DDRMat>&                aADVs,
                    Matrix<DDUMat>                 aPropertyVariableIndices,
                    Matrix<DDUMat>                 aADVIndices,
                    Matrix<DDRMat>                 aConstants,
                    Cell<std::shared_ptr<Field>>   aFieldDependencies,
                    MORIS_GEN_FIELD_FUNCTION       aFieldEvaluationFunction,
                    MORIS_GEN_SENSITIVITY_FUNCTION aSensitivityEvaluationFunction,
                    Field_Parameters               aParameters = {});

            /**
             * Constructor
             *
             * @param aOwnedADVs Owned distributed ADVs
             * @param aPropertyVariableIndices Indices of property variables to be filled by the ADVs
             * @param aADVIndices The indices of the ADV vector to fill in the property variables
             * @param aConstants The constant field variables not filled by ADVs
             * @param aFieldDependencies Other created fields that this property depends on
             * @param aFieldEvaluationFunction User-defined function for evaluating the property field
             * @param aSensitivityEvaluationFunction User-defined function for evaluating the field sensitivities
             * @param aParameters Additional parameters
             */
            User_Defined_Property(
                    sol::Dist_Vector*              aOwnedADVs,
                    Matrix<DDUMat>                 aPropertyVariableIndices,
                    Matrix<DDUMat>                 aADVIndices,
                    Matrix<DDRMat>                 aConstants,
                    Cell<std::shared_ptr<Field>>   aFieldDependencies,
                    MORIS_GEN_FIELD_FUNCTION       aFieldEvaluationFunction,
                    MORIS_GEN_SENSITIVITY_FUNCTION aSensitivityEvaluationFunction,
                    Field_Parameters               aParameters = {});

            /**
             * Given a node coordinate, returns the field value.
             *
             * @param aCoordinates Coordinate values
             * @return Property value
             */
            real get_field_value(const Matrix<DDRMat>& aCoordinates);

            /**
             * Given a node coordinate, evaluates the sensitivity of the proeprty field with respect to all of the
             * property variables.
             *
             * @param aCoordinates Coordinate values
             * @return Vector of sensitivities
             */
            const Matrix<DDRMat>& get_field_sensitivities(const Matrix<DDRMat>& aCoordinates);

        };
    }
}

#endif //MORIS_CL_GEN_USER_DEFINED_PROPERTY_HPP
