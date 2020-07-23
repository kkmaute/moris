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
                             sint            aNumRefinements,
                             sint            aRefinementFunctionIndex)
                : Field(aADVs, aGeometryVariableIndices, aADVIndices, aConstantParameters),
                  Geometry(aNumRefinements, aRefinementFunctionIndex),
                  mMesh(aMesh),
                  mNumOriginalNodes(aMesh->get_num_nodes())
        {
            MORIS_ASSERT(mFieldVariables.size() == mMesh->get_num_coeffs(0), "There must be a field variable for each "
                                                                             "B-spline coefficient in a level set geometry.");
        }

        //--------------------------------------------------------------------------------------------------------------

        real Level_Set::evaluate_field_value(uint aNodeIndex)
        {
            // Initialize
            real tResult = 0.0;

            // If node is on original interpolation mesh
            if (aNodeIndex < mNumOriginalNodes)
            {
                Matrix<IndexMat> tIndices = mMesh->get_bspline_inds_of_node_loc_ind(aNodeIndex, EntityRank::BSPLINE);
                Matrix<DDRMat> tMatrix = mMesh->get_t_matrix_of_node_loc_ind(aNodeIndex, EntityRank::BSPLINE);
                for (uint tBspline = 0; tBspline < tIndices.length(); tBspline++)
                {
                    tResult += tMatrix(tBspline) * (*mFieldVariables(tIndices(tBspline)));
                }
            }

            // Otherwise, use child nodes
            else
            {
                MORIS_ASSERT((aNodeIndex - mNumOriginalNodes) < mChildNodes.size(),
                             "A level set field value was requested from a node that this geometry doesn't know. "
                             "Perhaps a child node was not added to this field?");
                tResult = mChildNodes(aNodeIndex - mNumOriginalNodes)->interpolate_geometry_field_value(this);
            }

            // Return result
            return tResult;
        }

        //--------------------------------------------------------------------------------------------------------------

        void Level_Set::evaluate_all_sensitivities(uint aNodeIndex, Matrix<DDRMat>& aSensitivities)
        {
            // Initialize
            aSensitivities.set_size(1, mFieldVariables.size(), 0.0);

            // If node is on original interpolation mesh
            if (aNodeIndex < mNumOriginalNodes)
            {
                Matrix<IndexMat> tIndices = mMesh->get_bspline_inds_of_node_loc_ind(aNodeIndex, EntityRank::BSPLINE);
                Matrix<DDRMat> tMatrix = mMesh->get_t_matrix_of_node_loc_ind(aNodeIndex, EntityRank::BSPLINE);
                for (uint tBspline = 0; tBspline < tIndices.length(); tBspline++)
                {
                    aSensitivities(tIndices(tBspline)) = tMatrix(tBspline);
                }
            }

            // Otherwise, use child nodes
            else
            {
                MORIS_ASSERT((aNodeIndex - mNumOriginalNodes) < mChildNodes.size(),
                             "A level set sensitivity was requested from a node that this geometry doesn't know. "
                             "Perhaps a child node was not added to this field?");
                mChildNodes(aNodeIndex - mNumOriginalNodes)->interpolate_geometry_sensitivity(this, aSensitivities);
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Level_Set::add_child_node(uint aNodeIndex, std::shared_ptr<Child_Node> aChildNode)
        {
            mChildNodes.push_back(aChildNode);
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}