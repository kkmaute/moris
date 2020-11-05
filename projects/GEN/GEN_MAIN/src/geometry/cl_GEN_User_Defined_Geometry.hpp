#ifndef MORIS_CL_GEN_USER_DEFINED_GEOMETRY_HPP
#define MORIS_CL_GEN_USER_DEFINED_GEOMETRY_HPP

#include "cl_GEN_Geometry.hpp"
#include "cl_GEN_Field_Analytic.hpp"
#include "fn_Exec_load_user_library.hpp"

namespace moris
{
    namespace ge
    {
        class User_Defined_Geometry : public Geometry, public Field_Analytic
        {

        private:
            MORIS_GEN_FIELD_FUNCTION get_field_value_user_defined;
            MORIS_GEN_SENSITIVITY_FUNCTION get_field_sensitivities_user_defined;

        public:

            /**
             * Constructor, sets the pointers to advs and constant parameters for evaluations.
             *
             * @param aADVs Reference to the full advs
             * @param aGeometryVariableIndices Indices of geometry variables to be filled by the ADVs
             * @param aADVIndices The indices of the ADV vector to fill in the geometry variables
             * @param aConstantParameters The constant parameters not filled by ADVs
             * @param aFieldEvaluationFunction User-defined function for evaluating the geometry field
             * @param tSensitivitiesEvaluationFunction User-defined function for evaluating the field sensitivities
             * @param aName Name of this field for identification
             * @param aNumRefinements The number of refinement steps to use for this geometry
             * @param aRefinementFunctionIndex The index of a user-defined refinement function (-1 = default refinement)
             * @param aBSplineMeshIndex Index of a B-spline mesh for discretization (-2 = none, -1 = store nodal values)
             * @param aBSplineLowerBound The lower bound for the B-spline coefficients describing this field
             * @param aBSplineUpperBound The upper bound for the B-spline coefficients describing this field
             */
            User_Defined_Geometry(
                    Matrix<DDRMat>&                aADVs,
                    Matrix<DDUMat>                 aGeometryVariableIndices,
                    Matrix<DDUMat>                 aADVIndices,
                    Matrix<DDRMat>                 aConstantParameters,
                    MORIS_GEN_FIELD_FUNCTION       aFieldEvaluationFunction,
                    MORIS_GEN_SENSITIVITY_FUNCTION aSensitivityEvaluationFunction = nullptr,
                    std::string                    aName = "",
                    Matrix<DDSMat>  aNumRefinements = {{}},
                    Matrix<DDSMat>  aNumPatterns = {{}},
                    sint                           aRefinementFunctionIndex = -1,
                    sint                           aBSplineMeshIndex = -2,
                    real                           aBSplineLowerBound = -1.0,
                    real                           aBSplineUpperBound = 1.0);

            /**
             * Constructor, sets the field variable pointers to ADVs and constant parameters for evaluations.
             *
             * @param aOwnedADVs Pointer to the owned distributed ADVs
             * @param aFieldVariableIndices Indices of geometry variables to be filled by the ADVs
             * @param aADVIndices The indices of the ADV vector to fill in the geometry variables
             * @param aConstantParameters The constant parameters not filled by ADVs
             * @param aFieldEvaluationFunction User-defined function for evaluating the geometry field
             * @param aSensitivityEvaluationFunction User-defined function for evaluating the field sensitivities
             * @param aName Name of this field for identification
             * @param aNumRefinements The number of refinement steps to use for this field
             * @param aRefinementFunctionIndex The index of a user-defined refinement function (-1 = default refinement)
             * @param aBSplineMeshIndex Index of a B-spline mesh for discretization (-2 = none, -1 = store nodal values)
             * @param aBSplineLowerBound The lower bound for the B-spline coefficients describing this field
             * @param aBSplineUpperBound The upper bound for the B-spline coefficients describing this field
             */
            User_Defined_Geometry(
                    sol::Dist_Vector* aOwnedADVs,
                    Matrix<DDUMat>                 aGeometryVariableIndices,
                    Matrix<DDUMat>                 aADVIndices,
                    Matrix<DDRMat>                 aConstantParameters,
                    MORIS_GEN_FIELD_FUNCTION       aFieldEvaluationFunction,
                    MORIS_GEN_SENSITIVITY_FUNCTION aSensitivityEvaluationFunction = nullptr,
                    std::string                    aName = "",
                    Matrix<DDSMat>  aNumRefinements = {{}},
                    Matrix<DDSMat>  aNumPatterns = {{}},
                    sint                           aRefinementFunctionIndex = -1,
                    sint                           aBSplineMeshIndex = -2,
                    real                           aBSplineLowerBound = -1.0,
                    real                           aBSplineUpperBound = 1.0);

            /**
             * Constructor with only constant parameters
             *
             * @param aConstantParameters The constant parameters not filled by ADVs
             * @param aFieldEvaluationFunction User-defined function for evaluating the geometry field
             * @param aName Name of this field for identification
             * @param aNumRefinements The number of refinement steps to use for this geometry
             * @param aRefinementFunctionIndex The index of a user-defined refinement function (-1 = default refinement)
             * @param aBSplineMeshIndex Index of a B-spline mesh for discretization (-2 = none, -1 = store nodal values)
             * @param aBSplineLowerBound The lower bound for the B-spline coefficients describing this field
             * @param aBSplineUpperBound The upper bound for the B-spline coefficients describing this field
             */
            User_Defined_Geometry(
                    Matrix<DDRMat>           aConstantParameters,
                    MORIS_GEN_FIELD_FUNCTION aFieldEvaluationFunction,
                    std::string              aName = "",
                    Matrix<DDSMat>  aNumRefinements = {{}},
                    Matrix<DDSMat>  aNumPatterns = {{}},
                    sint                     aRefinementFunctionIndex = -1,
                    sint                     aBSplineMeshIndex = -2,
                    real                     aBSplineLowerBound = -1.0,
                    real                     aBSplineUpperBound = 1.0);

            /**
             * Given a node coordinate, returns the field value.
             *
             * @param aCoordinates Coordinate values
             * @return Distance to this geometry
             */
            real get_field_value(const Matrix<DDRMat>& aCoordinates);

            /**
             * Given a node coordinate, evaluates the sensitivity of the geometry field with respect to all of the
             * geometry variables.
             *
             * @param aCoordinates Coordinate values
             * @return Vector of sensitivities
             */
            Matrix<DDRMat> get_field_sensitivities(const Matrix<DDRMat>& aCoordinates);

        private:

            /**
             * Sets the user-defined functions. Eliminates redundant code since it's the same logic for all constructors.
             *
             * @param aFieldEvaluationFunction User-defined function for evaluating the geometry field
             * @param tSensitivitiesEvaluationFunction User-defined function for evaluating the field sensitivities
             */
            void set_user_defined_functions(
                    MORIS_GEN_FIELD_FUNCTION       aFieldEvaluationFunction,
                    MORIS_GEN_SENSITIVITY_FUNCTION aSensitivityEvaluationFunction);

            /**
             * Used internally to automatically error out if no sensitivities were provided
             */
            static void no_sensitivities(
                    const Matrix<DDRMat>& aCoordinates,
                    const Cell<real*>&    aParameters,
                    Matrix<DDRMat>&       aSensitivities);

        };
    }
}

#endif //MORIS_CL_GEN_USER_DEFINED_GEOMETRY_HPP
