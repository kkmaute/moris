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
             * @param aPropertyDependencies Other created properties that this property depends on@param aNumRefinements The number of refinement steps to use for this geometry
             * @param aRefinementFunctionIndex The index of a user-defined refinement function (-1 = default refinement)
             * @param aBSplineMeshIndex The index of a B-spline mesh for level set discretization (-1 = no B-splines)
             * @param aBSplineLowerBound The lower bound for the B-spline coefficients describing this field
             * @param aBSplineUpperBound The upper bound for the B-spline coefficients describing this field
             */
            User_Defined_Property(Matrix<DDRMat>& aADVs,
                                  Matrix<DDUMat> aPropertyVariableIndices,
                                  Matrix<DDUMat> aADVIndices,
                                  Matrix<DDRMat> aConstantParameters,
                                  Cell<std::shared_ptr<Property>> aPropertyDependencies,
                                  MORIS_GEN_FIELD_FUNCTION aFieldEvaluationFunction,
                                  MORIS_GEN_SENSITIVITY_FUNCTION aSensitivityEvaluationFunction,
                                  sint aNumRefinements = 0,
                                  sint aRefinementFunctionIndex = -1,
                                  sint aBSplineMeshIndex = -1,
                                  real aBSplineLowerBound = -1.0,
                                  real aBSplineUpperBound = 1.0);

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
