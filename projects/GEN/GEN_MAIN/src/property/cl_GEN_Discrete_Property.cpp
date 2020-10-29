#include "cl_GEN_Discrete_Property.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Discrete_Property::Discrete_Property(Matrix<DDRMat>& aADVs,
                                             Matrix<DDUMat>  aPropertyVariableIndices,
                                             Matrix<DDUMat>  aADVIndices,
                                             Matrix<DDRMat>  aConstantParameters,
                                             Matrix<DDSMat>  aNumRefinements,
                                             Matrix<DDSMat>  aNumPatterns,
                                             sint            aRefinementFunctionIndex,
                                             sint            aBSplineMeshIndex,
                                             real            aBSplineLowerBound,
                                             real            aBSplineUpperBound)
                : Field(aADVs,
                        aPropertyVariableIndices,
                        aADVIndices,
                        aConstantParameters,
                        aNumRefinements,
                        aNumPatterns,
                        aRefinementFunctionIndex,
                        aBSplineMeshIndex,
                        aBSplineLowerBound,
                        aBSplineUpperBound),
                  Field_Discrete(aADVs.length())
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        real Discrete_Property::evaluate_field_value(uint aNodeIndex)
        {
            return *mFieldVariables(aNodeIndex);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Discrete_Property::evaluate_all_sensitivities(uint aNodeIndex, Matrix<DDRMat>& aSensitivities)
        {
            aSensitivities.resize(1, mFieldVariables.size());
            aSensitivities(aNodeIndex) = 1;
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
