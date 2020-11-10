#include "cl_GEN_Discrete_Property.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Discrete_Property::Discrete_Property(
                Matrix<DDRMat>& aADVs,
                Matrix<DDUMat>  aPropertyVariableIndices,
                Matrix<DDUMat>  aADVIndices,
                Matrix<DDRMat>  aConstantParameters,
                std::string     aName,
                Matrix<DDSMat>  aNumRefinements,
                Matrix<DDSMat>  aRefMeshIndex,
                sint            aRefinementFunctionIndex,
                sint            aBSplineMeshIndex,
                real            aBSplineLowerBound,
                real            aBSplineUpperBound)
                : Field(aADVs,
                        aPropertyVariableIndices,
                        aADVIndices,
                        aConstantParameters,
                        aName,
                        aNumRefinements,
                        aRefMeshIndex,
                        aRefinementFunctionIndex,
                        aBSplineMeshIndex,
                        aBSplineLowerBound,
                        aBSplineUpperBound),
                  Field_Discrete_Integration(aADVs.length())
        {
            mSensitivities = {{1.0}};
        }

        //--------------------------------------------------------------------------------------------------------------

        Discrete_Property::Discrete_Property(
                sol::Dist_Vector* aOwnedADVs,
                Matrix<DDUMat>    aPropertyVariableIndices,
                Matrix<DDUMat>    aADVIndices,
                Matrix<DDRMat>    aConstantParameters,
                std::string       aName,
                Matrix<DDSMat>    aNumRefinements,
                Matrix<DDSMat>    aRefMeshIndex,
                sint              aRefinementFunctionIndex,
                sint              aBSplineMeshIndex,
                real              aBSplineLowerBound,
                real              aBSplineUpperBound)
                : Field(aOwnedADVs,
                        aPropertyVariableIndices,
                        aADVIndices,
                        aConstantParameters,
                        aName,
                        aNumRefinements,
                        aRefMeshIndex,
                        aRefinementFunctionIndex,
                        aBSplineMeshIndex,
                        aBSplineLowerBound,
                        aBSplineUpperBound),
                  Field_Discrete_Integration(aOwnedADVs->vec_local_length())
        {
            mSensitivities = {{1.0}};
        }

        //--------------------------------------------------------------------------------------------------------------

        real Discrete_Property::get_field_value(uint aNodeIndex)
        {
            return *mFieldVariables(aNodeIndex);
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix<DDRMat>& Discrete_Property::get_field_sensitivities(uint aNodeIndex)
        {
            return mSensitivities;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDSMat> Discrete_Property::get_determining_adv_ids(uint aNodeIndex)
        {
            return {{(sint)aNodeIndex}};
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
