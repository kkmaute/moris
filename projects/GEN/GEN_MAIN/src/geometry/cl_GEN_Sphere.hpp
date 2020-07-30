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
             * @param aLevelSetLowerBound The lower bound for a B-spline level set field describing this geometry
             * @param aLevelSetUpperBound The upper bound for a B-spline level set field describing this geometry
             */
            Sphere(Matrix<DDRMat>& aADVs,
                   Matrix<DDUMat>  aGeometryVariableIndices,
                   Matrix<DDUMat>  aADVIndices,
                   Matrix<DDRMat>  aConstantParameters,
                   sint            aNumRefinements = 0,
                   sint            aRefinementFunctionIndex = -1,
                   sint            aBSplineMeshIndex = -1,
                   real            aLevelSetLowerBound = -1.0,
                   real            aLevelSetUpperBound = 1.0);

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
             * @param aLevelSetLowerBound The lower bound for a B-spline level set field describing this geometry
             * @param aLevelSetUpperBound The upper bound for a B-spline level set field describing this geometry
             */
            Sphere(real aXCenter,
                   real aYCenter,
                   real aZCenter,
                   real aRadius,
                   sint aNumRefinements = 0,
                   sint aRefinementFunctionIndex = -1,
                   sint aBSplineMeshIndex = -1,
                   real aLevelSetLowerBound = -1.0,
                   real aLevelSetUpperBound = 1.0);

            /**
             * Given a node coordinate, this returns the distance to the nearest portion of the sphere's surface
             *
             * @param aCoordinates vector of coordinate values
             * @return distance to the sphere
             */
            real evaluate_field_value(const moris::Matrix<moris::DDRMat>& aCoordinates);

            /**
             * Given a node coordinate, returns a matrix of all sensitivities with respect to the sphere parameters
             *
             * @param aCoordinates Vector of coordinate values
             * @param aSensitivities Matrix of sensitivities
             */
            void evaluate_all_sensitivities(const moris::Matrix<moris::DDRMat>& aCoordinates, Matrix<DDRMat>& aSensitivities);

        };
    }
}

#endif /* MORIS_CL_GEN_SPHERE_HPP_ */
