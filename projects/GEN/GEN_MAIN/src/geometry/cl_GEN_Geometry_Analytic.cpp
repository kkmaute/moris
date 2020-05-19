#include "cl_GEN_Geometry_Analytic.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Geometry_Analytic::Geometry_Analytic(Matrix<DDRMat>& aADVs,
                                             Matrix<DDUMat> aGeometryVariableIndices,
                                             Matrix<DDUMat> aADVIndices,
                                             Matrix<DDRMat> aConstantParameters)
                                             : Field(aADVs,
                                                     aGeometryVariableIndices,
                                                     aADVIndices,
                                                     aConstantParameters)
        {

        }

        //--------------------------------------------------------------------------------------------------------------

        Geometry_Analytic::Geometry_Analytic(Matrix<DDRMat> aConstantParameters)
                : Field(aConstantParameters)
        {

        }

        //--------------------------------------------------------------------------------------------------------------
    }
}
