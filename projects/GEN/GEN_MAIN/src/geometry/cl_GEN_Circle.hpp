#ifndef MORIS_CL_GEN_CIRCLE_HPP
#define MORIS_CL_GEN_CIRCLE_HPP

#include "cl_GEN_Geometry.hpp"
#include "cl_Matrix.hpp"

namespace moris
{
    namespace ge
    {
        class Circle : public Geometry
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
            Circle(Matrix<DDRMat>& aADVs, Matrix<DDUMat> aGeometryVariableIndices, Matrix<DDUMat> aADVIndices, Matrix<DDRMat> aConstantParameters)
            : Geometry(aADVs, aGeometryVariableIndices, aADVIndices, aConstantParameters)
            {
            }

            /**
             * Constructor with only constant parameters
             *
             * @param aXCenter x-coordinate of the center of the circle
             * @param aYCenter y-coordiante of the center of the circle
             * @param aRadius radius of the circle
             */
            Circle(real aXCenter, real aYCenter, real aRadius) : Geometry(Matrix<DDRMat>({{aXCenter, aYCenter, aRadius}}))
            {
            }

            /**
             * Given a node coordinate, the geometry needs to return the distance to the nearest function.
             *
             * @param aCoordinates vector of coordinate values
             * @return distance to nearest function
             */
            real evaluate_field_value(const moris::Matrix<moris::DDRMat>& aCoordinates);

            /**
             * Given a node coordinate @param[in] aCoordinates, the function returns a matrix of relevant node coordinates
             * Where each row represents a design variable and each column is x, y, z sensitivities
             *
             * @param aCoordinates vector of coordinate values
             * @return matrix of sensitivities
             */
            moris::Matrix<moris::DDRMat> evaluate_sensitivity(const moris::Matrix<moris::DDRMat>& aCoordinates);

        };
    }
}

#endif /* MORIS_CL_GEN_CIRCLE_HPP */


