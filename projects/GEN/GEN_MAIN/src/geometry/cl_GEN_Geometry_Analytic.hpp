#ifndef MORIS_CL_GEN_GEOMETRY_ANALYTIC_HPP
#define MORIS_CL_GEN_GEOMETRY_ANALYTIC_HPP

#include "cl_GEN_Field_Base.hpp"

namespace moris
{
    namespace ge
    {
        class Geometry_Analytic : public Field
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
            Geometry_Analytic(Matrix<DDRMat>& aADVs, Matrix<DDUMat> aGeometryVariableIndices, Matrix<DDUMat> aADVIndices, Matrix<DDRMat> aConstantParameters);

            /**
             * Constructor for only constant parameters
             *
             * @param aConstantParameters The parameters that define this geometry
             */
            Geometry_Analytic(Matrix<DDRMat> aConstantParameters);
        };
    }
}

#endif /* MORIS_CL_GEN_GEOMETRY_ANALYTIC_HPP */
