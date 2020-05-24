//
// Created by christopherson on 5/18/20.
//

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
            MORIS_GEN_FIELD_FUNCTION evaluate_field_value_user_defined;
            MORIS_GEN_SENSITIVITY_FUNCTION evaluate_sensitivity_user_defined;

        public:

            /**
             * Constructor
             *
             * @param aADVs Reference to the full advs
             * @param aPropertyVariableIndices Indices of property variables to be filled by the ADVs
             * @param aADVIndices The indices of the ADV vector to fill in the property variables
             * @param aConstantParameters The constant parameters not filled by ADVs
             * @param aPropertyDependencies Other created properties that this property depends on
             */
            User_Defined_Property(Matrix<DDRMat>& aADVs,
                     Matrix<DDUMat> aPropertyVariableIndices,
                     Matrix<DDUMat> aADVIndices,
                     Matrix<DDRMat> aConstantParameters,
                     Cell<std::shared_ptr<Property>> aPropertyDependencies,
                     MORIS_GEN_FIELD_FUNCTION aFieldEvaluationFunction,
                     MORIS_GEN_SENSITIVITY_FUNCTION aSensitivityEvaluationFunction);

            /**
             * Evaluate the property field based on the given coordinates
             *
             * @param aCoordinates Coordinate values
             */
            real evaluate_field_value(const Matrix<DDRMat>& aCoordinates);

            /**
             * Evaluate the sensitivities of the property value with respect to all internal variables
             *
             * @param aCoordinates Coordinate values
             * @param aSensitivities Sensitivity matrix to fill
             */
            void evaluate_all_sensitivities(const Matrix<DDRMat>& aCoordinates, Matrix<DDRMat>& aSensitivities);

        };
    }   // end ge namespace
}   // end moris namespace

#endif //MORIS_CL_GEN_USER_DEFINED_PROPERTY_HPP
