#include "cl_GEN_Level_Set.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Level_Set::Level_Set(Matrix<DDRMat>& aADVs,
                             Matrix<DDUMat>  aGeometryVariableIndices,
                             Matrix<DDUMat>  aADVIndices,
                             Matrix<DDRMat>  aConstantParameters,
                             mtk::Mesh*      aMesh,
                             sint            aNumRefinements = 0,
                             sint            aRefinementFunctionIndex = -1)
                : Field(aADVs, aGeometryVariableIndices, aADVIndices, aConstantParameters),
                  Geometry(aNumRefinements, aRefinementFunctionIndex),
                  mMesh(aMesh),
                  mNumOriginalNodes(aMesh->get_num_nodes())
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        real Level_Set::evaluate_field_value(uint aNodeIndex)
        {
            if (aNodeIndex < mNumOriginalNodes)
            {
                return 0.0;
            }
            else
            {
                MORIS_ASSERT((aNodeIndex - mNumOriginalNodes) < mChildNodes.size(),
                             "A discrete level set field value was requested from a node that this field doesn't know. "
                             "Perhaps a child node was not added to this field?");
                return mChildNodes(aNodeIndex - mNumOriginalNodes)->interpolate_geometry_field_value(this);
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Level_Set::evaluate_all_sensitivities(uint aNodeIndex, Matrix<DDRMat>& aSensitivities)
        {

        }

        //--------------------------------------------------------------------------------------------------------------

        void Level_Set::add_child_node(uint aNodeIndex, std::shared_ptr<Child_Node> aChildNode)
        {
            mChildNodes.push_back(aChildNode);
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}