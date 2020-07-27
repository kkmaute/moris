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
            real ( Plane:: * m_eval_field )(const Matrix<DDRMat>&) = nullptr;
            void ( Plane:: * m_eval_sensitivity )(const Matrix<DDRMat>&, Matrix<DDRMat>&) = nullptr;

        public:
            /**
             * Constructor, sets the pointers to advs and constant parameters for evaluations
             *
             * @param aADVs Reference to the full advs
             * @param aGeometryVariableIndices Indices of geometry variables to be filled by the ADVs
             * @param aADVIndices The indices of the ADV vector to fill in the geometry variables
             * @param aConstantParameters The constant parameters not filled by ADVs
             * @param aNumRefinements The number of refinement steps to use for this geometry
             * @param aRefinementFunctionIndex The index of a user-defined refinement function (-1 = default refinement)
             * @param aBSplineMeshIndex The index of a B-spline mesh for level set discretization (-1 = no B-splines)
             * @param aLevelSetLowerBound The lower bound for a B-spline level set field describing this geometry
             * @param aLevelSetUpperBound The upper bound for a B-spline level set field describing this geometry
             */
            Plane(Matrix<DDRMat>& aADVs,
                  Matrix<DDUMat>  aGeometryVariableIndices,
                  Matrix<DDUMat>  aADVIndices,
                  Matrix<DDRMat>  aConstantParameters,
                  sint            aNumRefinements = 0,
                  sint            aRefinementFunctionIndex = -1,
                  sint            aBSplineMeshIndex = -1,
                  real            aLevelSetLowerBound = -1.0,
                  real            aLevelSetUpperBound = 1.0);

            /**
             * Constructor with only constant parameters, 3D
             *
             * @param aXCenter x-coordinate of the center of the plane
             * @param aYCenter y-coordinate of the center of the plane
             * @param aZCenter z-coordinate of the center of the plane
             * @param aXNormal x normal for the plane
             * @param aYNormal y normal for the plane
             * @param aZNormal z normal for the plane
             * @param aNumRefinements The number of refinement steps to use for this geometry
             * @param aRefinementFunctionIndex The index of a user-defined refinement function (-1 = default refinement)
             * @param aBSplineMeshIndex The index of a B-spline mesh for level set discretization (-1 = no B-splines)
             * @param aLevelSetLowerBound The lower bound for a B-spline level set field describing this geometry
             * @param aLevelSetUpperBound The upper bound for a B-spline level set field describing this geometry
             */
            Plane(real aXCenter,
                  real aYCenter,
                  real aZCenter,
                  real aXNormal,
                  real aYNormal,
                  real aZNormal,
                  sint aNumRefinements = 0,
                  sint aRefinementFunctionIndex = -1,
                  sint aBSplineMeshIndex = -1,
                  real aLevelSetLowerBound = -1.0,
                  real aLevelSetUpperBound = 1.0);

            /**
             * Constructor with only constant parameters, 2D
             *
             * @param aXCenter x-coordinate of the center of the plane
             * @param aYCenter y-coordinate of the center of the plane
             * @param aXNormal x normal for the plane
             * @param aYNormal y normal for the plane
             * @param aNumRefinements The number of refinement steps to use for this geometry
             * @param aRefinementFunctionIndex The index of a user-defined refinement function (-1 = default refinement)
             * @param aBSplineMeshIndex The index of a B-spline mesh for level set discretization (-1 = no B-splines)
             */
            Plane(real aXCenter,
                  real aYCenter,
                  real aXNormal,
                  real aYNormal,
                  sint aNumRefinements = 0,
                  sint aRefinementFunctionIndex = -1,
                  sint aBSplineMeshIndex = -1,
                  real aLevelSetLowerBound = -1.0,
                  real aLevelSetUpperBound = 1.0);

            /**
             * Given a node coordinate, the geometry needs to return the distance to the nearest function.
             *
             * @param aCoordinates vector of coordinate values
             * @return distance to nearest function
             */
            real evaluate_field_value(const Matrix<DDRMat>& aCoordinates);

            /**
             * Given a node coordinate @param aCoordinates, the function returns a matrix of sensitivities of the
             * geometry location with respect to the ADVs
             *
             * @param aCoordinates Vector of coordinate values
             * @param aSensitivities Matrix of sensitivities
             */
            void evaluate_all_sensitivities(const Matrix<DDRMat>& aCoordinates, Matrix<DDRMat>& aSensitivities);


        private:
            /**
             * 2D evaluation for evaluate_field_value
             */
            moris::real eval_field_2d(const Matrix<DDRMat>& aCoordinates);

            /**
             * 3D evaluation for evaluate_field_value
             */
            moris::real eval_field_3d(const Matrix<DDRMat>& aCoordinates);

            /**
             * 2D evaluation for evaluate_all_sensitivities
             */
            void eval_sensitivity_2d(const Matrix<DDRMat>& aCoordinates, Matrix<DDRMat>& aSensitivities);

            /**
             * 3D evaluation for evaluate_all_sensitivities
             */
            void eval_sensitivity_3d(const Matrix<DDRMat>& aCoordinates, Matrix<DDRMat>& aSensitivities);
        };
    }
}


#endif /* MORIS_CL_GEN_PLANE_HPP */
