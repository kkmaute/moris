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
                             sint            aRefinementFunctionIndex,
                             uint            aBSplineMeshIndex,
                             real            aLevelSetLowerBound,
                             real            aLevelSetUpperBound)
                : Field(aADVs, aGeometryVariableIndices, aADVIndices, aConstantParameters),
                  Geometry(aNumRefinements,
                           aRefinementFunctionIndex,
                           (sint)aBSplineMeshIndex,
                           aLevelSetLowerBound,
                           aLevelSetUpperBound),
                  mMesh(aMesh),
                  mNumOriginalNodes(aMesh->get_num_nodes())
        {
            // Check that number of variables equals the number of B-spline coefficients
            MORIS_ASSERT(mFieldVariables.size() == mMesh->get_num_coeffs(aBSplineMeshIndex),
                    "There must be a field variable for each B-spline coefficient in a level set geometry.");
        }

        //--------------------------------------------------------------------------------------------------------------

        Level_Set::Level_Set(Matrix<DDRMat>&           aADVs,
                             uint                      aADVIndex,
                             mtk::Mesh*                aMesh,
                             std::shared_ptr<Geometry> aGeometry)
                : Field(aADVs, aADVIndex, aMesh->get_num_coeffs(aGeometry->get_bspline_mesh_index())),
                  Geometry(aGeometry->get_num_refinements(),
                           aGeometry->get_refinement_function_index(),
                           aGeometry->get_bspline_mesh_index(),
                           aGeometry->get_level_set_lower_bound(),
                           aGeometry->get_level_set_upper_bound()),
                  mMesh(aMesh),
                  mNumOriginalNodes(aMesh->get_num_nodes())
        {
            // Check for linear B-splines
            MORIS_ERROR(mNumOriginalNodes == mMesh->get_num_coeffs(this->get_bspline_mesh_index()),
                    "GEN level sets are currently only supported for linear B-splines.");

            // Assign ADVs
            for (uint tNodeIndex = 0; tNodeIndex < mNumOriginalNodes; tNodeIndex++)
            {
                aADVs(aADVIndex++) = aGeometry->evaluate_field_value(tNodeIndex, mMesh->get_node_coordinate(tNodeIndex));
            }
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

        bool Level_Set::conversion_to_level_set()
        {
            return false;
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}