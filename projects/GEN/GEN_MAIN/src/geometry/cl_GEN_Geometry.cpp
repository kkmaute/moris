#include "cl_GEN_Geometry.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Geometry::Geometry(Matrix<DDRMat>& aADVs, Matrix<DDUMat> aGeometryVariableIndices, Matrix<DDUMat> aADVIndices, Matrix<DDRMat> aConstantParameters)
                : mConstantParameters(aConstantParameters)
        {
            // Check that the number of geometry variables indices equals the number of ADV indices, resize geometry variables
            MORIS_ERROR(aGeometryVariableIndices.length() == aADVIndices.length(),
                    "gen::Geometry: Number of geometry variables indices must equal the number of ADV indices");

            // Resize geometry variables
            uint tNumInputs = aGeometryVariableIndices.length() + mConstantParameters.length();
            mGeometryVariables.resize(tNumInputs);

            // Fill with pointers to ADVs
            for (uint tADVFillIndex = 0; tADVFillIndex < aGeometryVariableIndices.length(); tADVFillIndex++)
            {
                mGeometryVariables(aGeometryVariableIndices(tADVFillIndex)) = &(aADVs(aADVIndices(tADVFillIndex)));
            }

            // Fill with constant parameters
            uint tParameterIndex = 0;
            for (uint tVariableIndex = 0; tVariableIndex < tNumInputs; tVariableIndex++)
            {
                if (mGeometryVariables(tVariableIndex) == nullptr)
                {
                    mGeometryVariables(tVariableIndex) = &(mConstantParameters(tParameterIndex++));
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        Geometry::Geometry(Matrix<DDRMat> aConstantParameters)
                : mConstantParameters(aConstantParameters)
        {
            // Resize geometry variables
            uint tNumInputs = mConstantParameters.length();
            mGeometryVariables.resize(tNumInputs);

            // Fill geometry variables
            for (uint tVariableIndex = 0; tVariableIndex < tNumInputs; tVariableIndex++)
            {
                mGeometryVariables(tVariableIndex) = &(mConstantParameters(tVariableIndex));
            }
        }

        //--------------------------------------------------------------------------------------------------------------
    }
}
