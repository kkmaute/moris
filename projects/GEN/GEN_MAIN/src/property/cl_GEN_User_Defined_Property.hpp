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
             * @param aConstantParameters The constant parameters not filled by ADVs
             * @param aName Name of this field for identification
             * @param aFieldDependencies Other created fields that this property depends on
             * @param aNumRefinements The number of refinement steps to use for this property
             * @param aRefinementFunctionIndex The index of a user-defined refinement function (-1 = default refinement)
             * @param aBSplineMeshIndex The index of a B-spline mesh for level set discretization (-1 = no B-splines)
             * @param aBSplineLowerBound The lower bound for the B-spline coefficients describing this field
             * @param aBSplineUpperBound The upper bound for the B-spline coefficients describing this field
             */
            User_Defined_Property(
                    Matrix<DDRMat>&                aADVs,
                    Matrix<DDUMat>                 aPropertyVariableIndices,
                    Matrix<DDUMat>                 aADVIndices,
                    Matrix<DDRMat>                 aConstantParameters,
                    Cell<std::shared_ptr<Field>>   aFieldDependencies,
                    MORIS_GEN_FIELD_FUNCTION       aFieldEvaluationFunction,
                    MORIS_GEN_SENSITIVITY_FUNCTION aSensitivityEvaluationFunction,
                    std::string                    aName = "",
                    Matrix<DDSMat>  aNumRefinements = {{}},
                    Matrix<DDSMat>  aNumPatterns = {{}},
                    sint                           aRefinementFunctionIndex = -1,
                    sint                           aBSplineMeshIndex = -1,
                    real                           aBSplineLowerBound = -1.0,
                    real                           aBSplineUpperBound = 1.0);

            /**
             * Constructor
             *
             * @param aOwnedADVs Owned distributed ADVs
             * @param aPropertyVariableIndices Indices of property variables to be filled by the ADVs
             * @param aADVIndices The indices of the ADV vector to fill in the property variables
             * @param aConstantParameters The constant parameters not filled by ADVs
             * @param aName Name of this field for identification
             * @param aFieldDependencies Other created fields that this property depends on
             * @param aNumRefinements The number of refinement steps to use for this property
             * @param aRefinementFunctionIndex The index of a user-defined refinement function (-1 = default refinement)
             * @param aBSplineMeshIndex The index of a B-spline mesh for level set discretization (-1 = no B-splines)
             * @param aBSplineLowerBound The lower bound for the B-spline coefficients describing this field
             * @param aBSplineUpperBound The upper bound for the B-spline coefficients describing this field
             */
            User_Defined_Property(
                    sol::Dist_Vector*              aOwnedADVs,
                    Matrix<DDUMat>                 aPropertyVariableIndices,
                    Matrix<DDUMat>                 aADVIndices,
                    Matrix<DDRMat>                 aConstantParameters,
                    Cell<std::shared_ptr<Field>>   aFieldDependencies,
                    MORIS_GEN_FIELD_FUNCTION       aFieldEvaluationFunction,
                    MORIS_GEN_SENSITIVITY_FUNCTION aSensitivityEvaluationFunction,
                    std::string                    aName = "",
                    Matrix<DDSMat>  aNumRefinements = {{}},
                    Matrix<DDSMat>  aNumPatterns = {{}},
                    sint                           aRefinementFunctionIndex = -1,
                    sint                           aBSplineMeshIndex = -1,
                    real                           aBSplineLowerBound = -1.0,
                    real                           aBSplineUpperBound = 1.0);

            /**
             * Given a node coordinate, returns the field value.
             *
             * @param aCoordinates Coordinate values
             * @return Property value
             */
            real get_field_value_geometry(uint aNodeIndex,const Matrix<DDRMat>& aCoordinates);

            /**
             * Given a node coordinate, evaluates the sensitivity of the proeprty field with respect to all of the
             * property variables.
             *
             * @param aCoordinates Coordinate values
             * @return Vector of sensitivities
             */
            Matrix<DDRMat> get_field_sensitivities(const Matrix<DDRMat>& aCoordinates);

        };
    }
}

#endif //MORIS_CL_GEN_USER_DEFINED_PROPERTY_HPP
