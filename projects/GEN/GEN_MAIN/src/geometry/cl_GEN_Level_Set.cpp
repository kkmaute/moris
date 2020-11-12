#include "cl_GEN_Level_Set.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Mesh_Factory.hpp"
#include "cl_MTK_Mapper.hpp"
#include "cl_SOL_Matrix_Vector_Factory.hpp"
#include "cl_SOL_Dist_Map.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Level_Set::Level_Set(
                Matrix<DDRMat>&          aADVs,
                Matrix<DDUMat>           aGeometryVariableIndices,
                Matrix<DDUMat>           aADVIndices,
                Matrix<DDRMat>           aConstantParameters,
                mtk::Interpolation_Mesh* aMesh,
                std::string              aName,
                Matrix<DDSMat>           aNumRefinements,
                Matrix<DDSMat>           aRefinementMeshIndices,
                sint                     aRefinementFunctionIndex,
                uint                     aBSplineMeshIndex,
                real                     aBSplineLowerBound,
                real                     aBSplineUpperBound)
                : Field(aADVs,
                        aGeometryVariableIndices,
                        aADVIndices,
                        aConstantParameters,
                        aName,
                        aNumRefinements,
                        aRefinementMeshIndices,
                        aRefinementFunctionIndex,
                        aBSplineMeshIndex,
                        aBSplineLowerBound,
                        aBSplineUpperBound),
                  Field_Discrete_Integration(aMesh->get_num_nodes()),
                  mMesh(aMesh)
        {
            // Check that number of variables equals the number of B-spline coefficients
            MORIS_ASSERT(mFieldVariables.size() == mMesh->get_num_coeffs(aBSplineMeshIndex),
                    "There must be a field variable for each B-spline coefficient in a level set geometry.");
        }

        //--------------------------------------------------------------------------------------------------------------

        Level_Set::Level_Set(
                sol::Dist_Vector*         aOwnedADVs,
                const Matrix<DDSMat>&     aOwnedADVIds,
                const Matrix<DDSMat>&     aSharedADVIds,
                uint                      aOwnedADVIdsOffset,
                mtk::Interpolation_Mesh*  aMesh,
                std::shared_ptr<Geometry> aGeometry)
                : Field(aSharedADVIds,
                        aGeometry->get_name(),
                        aGeometry->get_num_refinements(),
                        aGeometry->get_refinement_mesh_indices(),
                        aGeometry->get_refinement_function_index(),
                        aGeometry->get_bspline_mesh_index(),
                        aGeometry->get_bspline_lower_bound(),
                        aGeometry->get_bspline_upper_bound()),
                  Field_Discrete_Integration(aMesh->get_num_nodes()),
                  mMesh(aMesh)
        {
            // Map to B-splines
            Matrix<DDRMat> tTargetField = this->map_to_bsplines(aGeometry);

            // Get B-spline mesh index
            uint tBSplineMeshIndex = this->get_bspline_mesh_index();

            // Assign ADVs
            for (uint tBSplineIndex = 0; tBSplineIndex < mMesh->get_num_coeffs(tBSplineMeshIndex); tBSplineIndex++)
            {
                if (par_rank() == aMesh->get_entity_owner(tBSplineIndex, EntityRank::BSPLINE, tBSplineMeshIndex))
                {
                    // Assign distributed vector element based on B-spline ID and offset
                    (*aOwnedADVs)(aOwnedADVIds(aOwnedADVIdsOffset++)) = tTargetField(tBSplineIndex);
                }
            }

            // Determine number of owned and shared node IDs
            uint tOwnedNodeCount = 0;
            for (uint tNodeIndex = 0; tNodeIndex < mMesh->get_num_nodes(); tNodeIndex++)
            {
                if (par_rank() == aMesh->get_entity_owner(tNodeIndex, EntityRank::NODE, tBSplineMeshIndex))
                {
                    tOwnedNodeCount++;
                }
            }
            Matrix<DDSMat> tOwnedNodeIDs(tOwnedNodeCount, 1);
            Matrix<DDSMat> tSharedNodeIDs(mMesh->get_num_nodes(), 1);

            // Assign owned and shared node IDs
            tOwnedNodeCount = 0;
            for (uint tNodeIndex = 0; tNodeIndex < mMesh->get_num_nodes(); tNodeIndex++)
            {
                sint tNodeID = mMesh->get_glb_entity_id_from_entity_loc_index(
                        tNodeIndex,
                        EntityRank::NODE,
                        tBSplineMeshIndex);
                tSharedNodeIDs(tNodeIndex) = tNodeID;
                if (par_rank() == aMesh->get_entity_owner(tNodeIndex, EntityRank::NODE, tBSplineMeshIndex))
                {
                    tOwnedNodeIDs(tOwnedNodeCount++) = tNodeID;
                }
            }

            // Create owned and shared distributed vectors
            sol::Matrix_Vector_Factory tDistributedFactory;
            sol::Dist_Map* tOwnedNodeMap = tDistributedFactory.create_map(tOwnedNodeIDs);
            sol::Dist_Map* tSharedNodeMap = tDistributedFactory.create_map(tSharedNodeIDs);
            mOwnedNodalValues = tDistributedFactory.create_vector(tOwnedNodeMap, 1, true);
            mSharedNodalValues = tDistributedFactory.create_vector(tSharedNodeMap, 1, true);

            // Import ADVs and assign nodal values
            this->import_advs(aOwnedADVs);

        }

        //--------------------------------------------------------------------------------------------------------------

        Level_Set::~Level_Set()
        {
            delete mOwnedNodalValues;
            delete mSharedNodalValues;
        }

        //--------------------------------------------------------------------------------------------------------------

        real Level_Set::get_field_value(uint aNodeIndex)
        {
            sint tNodeID = mMesh->get_glb_entity_id_from_entity_loc_index(
                    aNodeIndex,
                    EntityRank::NODE,
                    this->get_bspline_mesh_index());

            return (*mSharedNodalValues)(tNodeID);
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix<DDRMat>& Level_Set::get_field_sensitivities(uint aNodeIndex)
        {
            return mMesh->get_t_matrix_of_node_loc_ind(aNodeIndex, this->get_bspline_mesh_index());
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDSMat> Level_Set::get_determining_adv_ids(uint aNodeIndex)
        {
            return mMesh->get_bspline_ids_of_node_loc_ind(aNodeIndex, this->get_bspline_mesh_index());
        }

        //--------------------------------------------------------------------------------------------------------------

        void Level_Set::import_advs(sol::Dist_Vector* aOwnedADVs)
        {
            // Import ADVs as usual
            Field::import_advs(aOwnedADVs);

            // Reset evaluated field
            mOwnedNodalValues->vec_put_scalar(0);

            // Evaluate field at owned nodes
            for (uint tNodeIndex = 0; tNodeIndex < mMesh->get_num_nodes(); tNodeIndex++)
            {
                if (par_rank() == mMesh->get_entity_owner(tNodeIndex, EntityRank::NODE, this->get_bspline_mesh_index()))
                {
                    sint tNodeID = mMesh->get_glb_entity_id_from_entity_loc_index(
                            tNodeIndex,
                            EntityRank::NODE,
                            this->get_bspline_mesh_index());
                    Matrix<IndexMat> tBSplineIndices = mMesh->get_bspline_inds_of_node_loc_ind(tNodeIndex, this->get_bspline_mesh_index());
                    Matrix<DDRMat> tMatrix = mMesh->get_t_matrix_of_node_loc_ind(tNodeIndex, this->get_bspline_mesh_index());
                    for (uint tBSpline = 0; tBSpline < tBSplineIndices.length(); tBSpline++)
                    {
                        (*mOwnedNodalValues)(tNodeID) += tMatrix(tBSpline) * (*mFieldVariables(tBSplineIndices(tBSpline)));
                    }
                }
            }

            // Global assembly
            mOwnedNodalValues->vector_global_assembly();

            // Import nodal values
            mSharedNodalValues->import_local_to_global(*mOwnedNodalValues);
        }

        //--------------------------------------------------------------------------------------------------------------

        bool Level_Set::conversion_to_bsplines()
        {
            return false;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Level_Set::map_to_bsplines(std::shared_ptr<Geometry> aGeometry)
        {
            // Create field
            Matrix<DDRMat> tTargetField(0, 0);

            // Find first node with interpolation
            uint tTestNodeIndex = 0;
            while (not mMesh->node_has_interpolation(tTestNodeIndex, this->get_bspline_mesh_index()))
            {
                tTestNodeIndex++;
            }

            // Check if L2 is needed
            if (mMesh->get_bspline_inds_of_node_loc_ind(tTestNodeIndex, this->get_bspline_mesh_index()).length() > 1)
            {
                // Check for serial
                MORIS_ERROR(par_size() == 1, "L2 not implemented in parallel");

                // Set values ons source field
                Matrix<DDRMat> tSourceField(mNumOriginalNodes, 1);
                for (uint tNodeIndex = 0; tNodeIndex < mNumOriginalNodes; tNodeIndex++)
                {
                    tSourceField(tNodeIndex) = aGeometry->get_field_value(tNodeIndex, mMesh->get_node_coordinate(tNodeIndex));
                }

                // Create integration mesh
                mtk::Integration_Mesh * tIntegrationMesh = create_integration_mesh_from_interpolation_mesh(MeshType::HMR, mMesh);

                // Create mesh manager
                std::shared_ptr<mtk::Mesh_Manager> tMeshManager = std::make_shared<mtk::Mesh_Manager>();

                // Register mesh pair
                uint tMeshIndex = tMeshManager->register_mesh_pair(mMesh, tIntegrationMesh);

                // Use mapper
                mapper::Mapper tMapper(tMeshManager, tMeshIndex, (uint)this->get_bspline_mesh_index());
                tMapper.perform_mapping(tSourceField,
                                        EntityRank::NODE,
                                        tTargetField,
                                        EntityRank::BSPLINE);
            }
            else // Nodal values, no L2
            {
                // Target field needs to be created because map from B-spline ID to node is not known
                tTargetField.resize(mMesh->get_num_coeffs(this->get_bspline_mesh_index()), 1);

                // Assign ADVs directly
                for (uint tNodeIndex = 0; tNodeIndex < mNumOriginalNodes; tNodeIndex++)
                {
                    if (mMesh->node_has_interpolation(tNodeIndex, this->get_bspline_mesh_index()))
                    {
                        // Get B-spline index
                        uint tBSplineIndex = mMesh->get_bspline_inds_of_node_loc_ind(tNodeIndex, this->get_bspline_mesh_index())(0);

                        // Assign target field
                        tTargetField(tBSplineIndex) = aGeometry->get_field_value(tNodeIndex, mMesh->get_node_coordinate(tNodeIndex));
                    }
                }
            }

            // Return mapped field
            return tTargetField;
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
