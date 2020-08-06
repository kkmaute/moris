#ifndef MORIS_CL_GEN_SUPERELLIPSOID_HPP
#define MORIS_CL_GEN_SUPERELLIPSOID_HPP

#include "cl_GEN_Geometry.hpp"
#include "cl_GEN_Field_Analytic.hpp"

namespace moris
{
    namespace ge
    {
        class Superellipsoid : public Geometry, public Field_Analytic
        {
        private:
            real mEpsilon = 1E-8;

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
            Superellipsoid(Matrix<DDRMat>& aADVs,
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
             * @param aXCenter x-coordinate of the center of the superellipse
             * @param aYCenter y-coordiante of the center of the superellipse
             * @param aXSemidiameter Superellipsoid semi-diameter in the x direction
             * @param aYSemidiameter Superellipsoid semi-diameter in the y direction
             * @param aExponent Superellipsoid exponent
             * @param aNumRefinements The number of refinement steps to use for this geometry
             * @param aRefinementFunctionIndex The index of a user-defined refinement function (-1 = default refinement)
             * @param aBSplineMeshIndex The index of a B-spline mesh for level set discretization (-1 = no B-splines)
             * @param aBSplineLowerBound The lower bound for the B-spline coefficients describing this field
             * @param aBSplineUpperBound The upper bound for the B-spline coefficients describing this field
             */
            Superellipsoid(real aXCenter,
                           real aYCenter,
                           real aZCenter,
                           real aXSemidiameter,
                           real aYSemidiameter,
                           real aZSemidiameter,
                           real aExponent,
                           sint aNumRefinements = 0,
                           sint aRefinementFunctionIndex = -1,
                           sint aBSplineMeshIndex = -1,
                           real aBSplineLowerBound = -1.0,
                           real aBSplineUpperBound = 1.0);

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

        };
    }
}

#endif //MORIS_CL_GEN_SUPERELLIPSOID_HPP
