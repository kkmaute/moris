#ifndef MORIS_CL_GEN_PLANE_HPP
#define MORIS_CL_GEN_PLANE_HPP

#include "cl_GEN_Geometry.hpp"
#include "cl_GEN_Field_Analytic.hpp"

namespace moris
{
    namespace ge
    {
        class Plane : public Geometry, public Field_Analytic
        {
        private:
            real (Plane::*m_eval_field)(const Matrix<DDRMat>&) = nullptr;
            Matrix<DDRMat> (Plane::*m_eval_sensitivity)(const Matrix<DDRMat>&) = nullptr;

        public:
            /**
             * Constructor, sets the pointers to advs and constant parameters for evaluations.
             *
             * @param aADVs Reference to the full advs
             * @param aGeometryVariableIndices Indices of geometry variables to be filled by the ADVs
             * @param aADVIndices The indices of the ADV vector to fill in the geometry variables
             * @param aConstantParameters The constant parameters not filled by ADVs
             * @param aName Name of this field for identification
             * @param aNumRefinements The number of refinement steps to use for this geometry
             * @param aRefinementFunctionIndex The index of a user-defined refinement function (-1 = default refinement)
             * @param aBSplineMeshIndex The index of a B-spline mesh for level set discretization (-1 = no B-splines)
             * @param aBSplineLowerBound The lower bound for the B-spline coefficients describing this field
             * @param aBSplineUpperBound The upper bound for the B-spline coefficients describing this field
             */
            Plane(Matrix<DDRMat>& aADVs,
                  Matrix<DDUMat>  aGeometryVariableIndices,
                  Matrix<DDUMat>  aADVIndices,
                  Matrix<DDRMat>  aConstantParameters,
                  std::string     aName = "",
                  sint            aNumRefinements = 0,
                  sint            aRefinementFunctionIndex = -1,
                  sint            aBSplineMeshIndex = -1,
                  real            aBSplineLowerBound = -1.0,
                  real            aBSplineUpperBound = 1.0);

            /**
             * Constructor, sets the field variable pointers to ADVs and constant parameters for evaluations.
             *
             * @param aOwnedADVs Pointer to the owned distributed ADVs
             * @param aFieldVariableIndices Indices of geometry variables to be filled by the ADVs
             * @param aADVIndices The indices of the ADV vector to fill in the geometry variables
             * @param aConstantParameters The constant parameters not filled by ADVs
             * @param aName Name of this field for identification
             * @param aNumRefinements The number of refinement steps to use for this field
             * @param aRefinementFunctionIndex The index of a user-defined refinement function (-1 = default refinement)
             * @param aBSplineMeshIndex The index of a B-spline mesh for B-spline discretization (-1 = no B-splines)
             * @param aBSplineLowerBound The lower bound for the B-spline coefficients describing this field
             * @param aBSplineUpperBound The upper bound for the B-spline coefficients describing this field
             */
            Plane(sol::Dist_Vector* aOwnedADVs,
                  Matrix<DDUMat>    aGeometryVariableIndices,
                  Matrix<DDUMat>    aADVIndices,
                  Matrix<DDRMat>    aConstantParameters,
                  std::string       aName = "",
                  sint              aNumRefinements = 0,
                  sint              aRefinementFunctionIndex = -1,
                  sint              aBSplineMeshIndex = -1,
                  real              aBSplineLowerBound = -1.0,
                  real              aBSplineUpperBound = 1.0);

            /**
             * Constructor with only constant parameters, 3D
             *
             * @param aXCenter x-coordinate of the center of the plane
             * @param aYCenter y-coordinate of the center of the plane
             * @param aZCenter z-coordinate of the center of the plane
             * @param aXNormal x normal for the plane
             * @param aYNormal y normal for the plane
             * @param aZNormal z normal for the plane
             * @param aName Name of this field for identification
             * @param aNumRefinements The number of refinement steps to use for this geometry
             * @param aRefinementFunctionIndex The index of a user-defined refinement function (-1 = default refinement)
             * @param aBSplineMeshIndex The index of a B-spline mesh for level set discretization (-1 = no B-splines)
             * @param aBSplineLowerBound The lower bound for the B-spline coefficients describing this field
             * @param aBSplineUpperBound The upper bound for the B-spline coefficients describing this field
             */
            Plane(real        aXCenter,
                  real        aYCenter,
                  real        aZCenter,
                  real        aXNormal,
                  real        aYNormal,
                  real        aZNormal,
                  std::string aName = "",
                  sint        aNumRefinements = 0,
                  sint        aRefinementFunctionIndex = -1,
                  sint        aBSplineMeshIndex = -1,
                  real        aBSplineLowerBound = -1.0,
                  real        aBSplineUpperBound = 1.0);

            /**
             * Constructor with only constant parameters, 2D
             *
             * @param aXCenter x-coordinate of the center of the plane
             * @param aYCenter y-coordinate of the center of the plane
             * @param aXNormal x normal for the plane
             * @param aYNormal y normal for the plane
             * @param aName Name of this field for identification
             * @param aNumRefinements The number of refinement steps to use for this geometry
             * @param aRefinementFunctionIndex The index of a user-defined refinement function (-1 = default refinement)
             * @param aBSplineMeshIndex The index of a B-spline mesh for level set discretization (-1 = no B-splines)
             */
            Plane(real        aXCenter,
                  real        aYCenter,
                  real        aXNormal,
                  real        aYNormal,
                  std::string aName = "",
                  sint        aNumRefinements = 0,
                  sint        aRefinementFunctionIndex = -1,
                  sint        aBSplineMeshIndex = -1,
                  real        aBSplineLowerBound = -1.0,
                  real        aBSplineUpperBound = 1.0);

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
             * 2D evaluation for get_field_value
             */
            real eval_field_2d(const Matrix<DDRMat>& aCoordinates);

            /**
             * 3D evaluation for get_field_value
             */
            real eval_field_3d(const Matrix<DDRMat>& aCoordinates);

            /**
             * 2D evaluation for get_field_sensitivities
             */
            Matrix<DDRMat> eval_sensitivity_2d(const Matrix<DDRMat>& aCoordinates);

            /**
             * 3D evaluation for get_field_sensitivities
             */
            Matrix<DDRMat> eval_sensitivity_3d(const Matrix<DDRMat>& aCoordinates);
        };
    }
}


#endif /* MORIS_CL_GEN_PLANE_HPP */
