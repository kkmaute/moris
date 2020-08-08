#ifndef MORIS_CL_GEN_SPHERE_HPP_
#define MORIS_CL_GEN_SPHERE_HPP_

#include <cmath>

#include "cl_GEN_Geometry.hpp"
#include "cl_GEN_Field_Analytic.hpp"

namespace moris
{
    namespace ge
    {
        class Sphere : public Geometry, public Field_Analytic
        {
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
             * @param aBSplineLowerBound The lower bound for the B-spline coefficients describing this field
             * @param aBSplineUpperBound The upper bound for the B-spline coefficients describing this field
             */
            Sphere(Matrix<DDRMat>& aADVs,
                   Matrix<DDUMat>  aGeometryVariableIndices,
                   Matrix<DDUMat>  aADVIndices,
                   Matrix<DDRMat>  aConstantParameters,
                   sint            aNumRefinements = 0,
                   sint            aRefinementFunctionIndex = -1,
                   sint            aBSplineMeshIndex = -1,
                   real            aBSplineLowerBound = -1.0,
                   real            aBSplineUpperBound = 1.0);

            /**
             * Constructor with only constant parameters
             *
             * @param aXCenter x-coordinate of the center of the sphere
             * @param aYCenter y-coordiante of the center of the sphere
             * @param aZCenter z-coordinate of the center of the sphere
             * @param aRadius radius of the sphere
             * @param aNumRefinements The number of refinement steps to use for this geometry
             * @param aRefinementFunctionIndex The index of a user-defined refinement function (-1 = default refinement)
             * @param aBSplineMeshIndex The index of a B-spline mesh for level set discretization (-1 = no B-splines)
             * @param aBSplineLowerBound The lower bound for the B-spline coefficients describing this field
             * @param aBSplineUpperBound The upper bound for the B-spline coefficients describing this field
             */
            Sphere(real aXCenter,
                   real aYCenter,
                   real aZCenter,
                   real aRadius,
                   sint aNumRefinements = 0,
                   sint aRefinementFunctionIndex = -1,
                   sint aBSplineMeshIndex = -1,
                   real aBSplineLowerBound = -1.0,
                   real aBSplineUpperBound = 1.0);

            /**
             * Given a node coordinate, returns the field value.
             *
             * @param aCoordinates Coordinate values
             * @return Distance to this geometry
             */
            real evaluate_field_value(const Matrix<DDRMat>& aCoordinates);

            /**
             * Given a node coordinate, evaluates the sensitivity of the geometry field with respect to all of the
             * geometry variables.
             *
             * @param aCoordinates Coordinate values
             * @param aSensitivities Vector of sensitivities
             */
            void evaluate_all_sensitivities(const Matrix<DDRMat>& aCoordinates, Matrix<DDRMat>& aSensitivities);

        };
    }
}

#endif /* MORIS_CL_GEN_SPHERE_HPP_ */
