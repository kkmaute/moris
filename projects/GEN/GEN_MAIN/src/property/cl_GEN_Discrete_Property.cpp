//
// Created by christopherson on 5/19/20.
//

#include "cl_GEN_Discrete_Property.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Discrete_Property::Discrete_Property(Matrix<DDRMat>& aADVs,
                                             Matrix<DDUMat> aPropertyVariableIndices,
                                             Matrix<DDUMat> aADVIndices,
                                             Matrix<DDRMat> aConstantParameters)
                                             : Field(aADVs,
                                                     aPropertyVariableIndices,
                                                     aADVIndices,
                                                     aConstantParameters)
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
