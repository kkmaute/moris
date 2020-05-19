#ifndef MORIS_CL_GEN_CIRCLE_HPP
#define MORIS_CL_GEN_CIRCLE_HPP

#include "cl_GEN_Geometry_Analytic.hpp"
#include "cl_Matrix.hpp"

namespace moris
{
    namespace ge
    {
        class Circle : public Geometry_Analytic
        {
        public:

            /**
             * Constructor, sets the pointers to advs and constant parameters for evaluations
             *
             * @param aADVs Reference to the full advs
             * @param aGeometryVariableIndices Indices of geometry variables to be filled by the ADVs
             * @param aADVIndices The indices of the ADV vector to fill in the geometry variables
             * @param aConstantParameters The constant parameters not filled by ADVs
             */
            Circle(Matrix<DDRMat>& aADVs, Matrix<DDUMat> aGeometryVariableIndices, Matrix<DDUMat> aADVIndices, Matrix<DDRMat> aConstantParameters);

            /**
             * Constructor with only constant parameters
             *
             * @param aXCenter x-coordinate of the center of the circle
             * @param aYCenter y-coordiante of the center of the circle
             * @param aRadius radius of the circle
             */
            Circle(real aXCenter, real aYCenter, real aRadius);

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

#endif /* MORIS_CL_GEN_CIRCLE_HPP */


