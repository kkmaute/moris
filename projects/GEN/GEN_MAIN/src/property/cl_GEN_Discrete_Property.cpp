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
                                             sint            aNumRefinements,
                                             sint            aRefinementFunctionIndex,
                                             sint            aBSplineMeshIndex,
                                             real            aBSplineLowerBound,
                                             real            aBSplineUpperBound)
                : Field(aADVs,
                        aPropertyVariableIndices,
                        aADVIndices,
                        aConstantParameters,
                        aNumRefinements,
                        aRefinementFunctionIndex,
                        aBSplineMeshIndex,
                        aBSplineLowerBound,
                        aBSplineUpperBound),
                  Field_Discrete(aADVs.length())
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        Discrete_Property::Discrete_Property(sol::Dist_Vector* aOwnedADVs,
                                             Matrix<DDUMat>    aPropertyVariableIndices,
                                             Matrix<DDUMat>    aADVIndices,
                                             Matrix<DDRMat>    aConstantParameters,
                                             sint              aNumRefinements,
                                             sint              aRefinementFunctionIndex,
                                             sint              aBSplineMeshIndex,
                                             real              aBSplineLowerBound,
                                             real              aBSplineUpperBound)
                : Field(aOwnedADVs,
                        aPropertyVariableIndices,
                        aADVIndices,
                        aConstantParameters,
                        aNumRefinements,
                        aRefinementFunctionIndex,
                        aBSplineMeshIndex,
                        aBSplineLowerBound,
                        aBSplineUpperBound),
                  Field_Discrete(aOwnedADVs->vec_local_length())
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        real Discrete_Property::get_field_value(uint aNodeIndex)
        {
            return *mFieldVariables(aNodeIndex);
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Discrete_Property::get_field_sensitivities(uint aNodeIndex)
        {
            return {{1.0}};
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDSMat> Discrete_Property::get_determining_adv_ids(uint aNodeIndex)
        {
            return {{(sint)aNodeIndex}};
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
