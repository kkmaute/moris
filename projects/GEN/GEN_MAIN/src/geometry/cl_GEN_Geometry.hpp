#ifndef MORIS_CL_GEN_GEOMETRY_ANALYTIC_HPP
#define MORIS_CL_GEN_GEOMETRY_ANALYTIC_HPP

#include "cl_Matrix.hpp"

namespace moris
{
    namespace ge
    {
        class Geometry
        {
        protected:
            Cell<real*> mGeometryVariables;
            Matrix<DDRMat> mConstantParameters;

            /**
             * Constructor, sets the pointers to advs and constant parameters for evaluations
             *
             * @param aADVs Reference to the full advs
             * @param aGeometryVariableIndices Indices of geometry variables to be filled by the ADVs
             * @param aADVIndices The indices of the ADV vector to fill in the geometry variables
             * @param aConstantParameters The constant parameters not filled by ADVs
             */
            Geometry(Matrix<DDRMat>& aADVs, Matrix<DDUMat> aGeometryVariableIndices, Matrix<DDUMat> aADVIndices, Matrix<DDRMat> aConstantParameters);

            /**
             * Constructor for only constant parameters
             *
             * @param aConstantParameters The parameters that define this geometry
             */
            Geometry(Matrix<DDRMat> aConstantParameters);

        public:
            /**
             * Destructor
             */
            ~Geometry()
            {
            }

            /**
             * Given a node coordinate, the geometry needs to return the distance to the nearest function.
             *
             * @param aCoordinates vector of coordinate values
             * @return distance to nearest function
             */
            virtual real evaluate_field_value(const moris::Matrix<moris::DDRMat>& aCoordinates) = 0;

            /**
             * Given a node coordinate @param[in] aCoordinates, the function returns a matrix of relevant node coordinates
             * Where each row represents a design variable and each column is x, y, z sensitivities
             *
             * @param aCoordinates vector of coordinate values
             * @return matrix of sensitivities
             */
            virtual moris::Matrix<moris::DDRMat> evaluate_sensitivity(const moris::Matrix<moris::DDRMat>& aCoordinates) = 0;

        };
    }
}

#endif /* MORIS_CL_GEN_GEOMETRY_ANALYTIC_HPP */
