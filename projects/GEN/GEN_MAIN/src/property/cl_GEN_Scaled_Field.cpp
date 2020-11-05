#include "cl_GEN_Scaled_Field.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Scaled_Field::Scaled_Field(
                Matrix<DDRMat>&              aADVs,
                Matrix<DDUMat>               aPropertyVariableIndices,
                Matrix<DDUMat>               aADVIndices,
                Matrix<DDRMat>               aConstantParameters,
                Cell<std::shared_ptr<Field>> aFieldDependencies,
                std::string                  aName,
                Matrix<DDSMat>  aNumRefinements,
                Matrix<DDSMat>  aRefinementMeshIndices,
                sint                         aRefinementFunctionIndex,
                sint                         aBSplineMeshIndex,
                real                         aBSplineLowerBound,
                real                         aBSplineUpperBound)
                : Field(aADVs,
                        aPropertyVariableIndices,
                        aADVIndices,
                        aConstantParameters,
                        aName,
                        aNumRefinements,
                        aRefinementMeshIndices,
                        aRefinementFunctionIndex,
                        aBSplineMeshIndex,
                        aBSplineLowerBound,
                        aBSplineUpperBound),
                  Property(aFieldDependencies)
        {
            MORIS_ERROR(mFieldDependencies.size() == 1, "A scaled field property must depend on one field.");
            MORIS_ERROR(mFieldVariables.size() == 1, "A scaled field property must have one scaling factor.");
            MORIS_ERROR(aPropertyVariableIndices.length() == 0,
                    "A scaled field property must have a constant scaling factor for now.");
        }

        //--------------------------------------------------------------------------------------------------------------

        Scaled_Field::Scaled_Field(
                sol::Dist_Vector*            aOwnedADVs,
                Matrix<DDUMat>               aPropertyVariableIndices,
                Matrix<DDUMat>               aADVIndices,
                Matrix<DDRMat>               aConstantParameters,
                Cell<std::shared_ptr<Field>> aFieldDependencies,
                std::string                  aName,
                Matrix<DDSMat>  aNumRefinements,
                Matrix<DDSMat>  aRefinementMeshIndices,
                sint                         aRefinementFunctionIndex,
                sint                         aBSplineMeshIndex,
                real                         aBSplineLowerBound,
                real                         aBSplineUpperBound)
                : Field(aOwnedADVs,
                        aPropertyVariableIndices,
                        aADVIndices,
                        aConstantParameters,
                        aName,
                        aNumRefinements,
                        aRefinementMeshIndices,
                        aRefinementFunctionIndex,
                        aBSplineMeshIndex,
                        aBSplineLowerBound,
                        aBSplineUpperBound),
                  Property(aFieldDependencies)
        {
            MORIS_ERROR(mFieldDependencies.size() == 1, "A scaled field property must depend on one field.");
            MORIS_ERROR(mFieldVariables.size() == 1, "A scaled field property must have one scaling factor.");
            MORIS_ERROR(aPropertyVariableIndices.length() == 0,
                        "A scaled field property must have a constant scaling factor for now.");
        }

        //--------------------------------------------------------------------------------------------------------------

        real Scaled_Field::get_field_value(uint aNodeIndex, const Matrix<DDRMat>& aCoordinates)
        {
            return *mFieldVariables(0) * mFieldDependencies(0)->get_field_value(aNodeIndex, aCoordinates);
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Scaled_Field::get_field_sensitivities(uint aNodeIndex, const Matrix<DDRMat>& aCoordinates)
        {
            return *mFieldVariables(0) * mFieldDependencies(0)->get_field_sensitivities(aNodeIndex, aCoordinates);
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDSMat> Scaled_Field::get_determining_adv_ids(uint aNodeIndex, const Matrix<DDRMat>& aCoordinates)
        {
            return mFieldDependencies(0)->get_determining_adv_ids(aNodeIndex, aCoordinates);
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
