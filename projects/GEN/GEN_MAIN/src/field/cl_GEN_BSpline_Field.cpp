#include "cl_GEN_BSpline_Field.hpp"
#include "st_MTK_Mesh_Pair.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Mesh_Factory.hpp"
#include "cl_MTK_Field.hpp"
#include "cl_MTK_Mapper.hpp"
#include "cl_SOL_Matrix_Vector_Factory.hpp"
#include "cl_SOL_Dist_Map.hpp"
#include "fn_trans.hpp"

// Logging package
#include "cl_Logger.hpp"
#include "cl_Tracer.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        BSpline_Field::BSpline_Field(
                sol::Dist_Vector*        aOwnedADVs,
                const Matrix<DDSMat>&    aOwnedADVIds,
                const Matrix<DDSMat>&    aSharedADVIds,
                uint                     aOwnedADVIdsOffset,
                mtk::Interpolation_Mesh* aMesh,
                std::shared_ptr<Field>   aField)
                : Field(aSharedADVIds, aField)
                , Field_Discrete_Integration(aMesh->get_num_nodes())
                , mMesh(aMesh)
        {
            // Map to B-splines
            Matrix<DDRMat> tTargetField = this->map_to_bsplines(aField);

            // Get B-spline mesh index
            uint tBSplineMeshIndex = this->get_discretization_mesh_index();

            // Assign ADVs
            for (uint tBSplineIndex = 0; tBSplineIndex < mMesh->get_num_coeffs(tBSplineMeshIndex); tBSplineIndex++)
            {
                if ((uint) par_rank() == aMesh->get_entity_owner(tBSplineIndex, EntityRank::BSPLINE, tBSplineMeshIndex))
                {
                    // Assign distributed vector element based on B-spline ID and offset
                    (*aOwnedADVs)(aOwnedADVIds(aOwnedADVIdsOffset++)) = tTargetField(tBSplineIndex);
                }
            }

            // Determine number of owned and shared node IDs
            uint tOwnedNodeCount = 0;
            for (uint tNodeIndex = 0; tNodeIndex < mMesh->get_num_nodes(); tNodeIndex++)
            {
                if ((uint) par_rank() == aMesh->get_entity_owner(tNodeIndex, EntityRank::NODE, tBSplineMeshIndex))
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
                if ((uint) par_rank() == aMesh->get_entity_owner(tNodeIndex, EntityRank::NODE, tBSplineMeshIndex))
                {
                    tOwnedNodeIDs(tOwnedNodeCount++) = tNodeID;
                }
            }

            // Create owned and shared distributed vectors
            sol::Matrix_Vector_Factory tDistributedFactory;
            sol::Dist_Map* tOwnedNodeMap = tDistributedFactory.create_map(tOwnedNodeIDs);
            sol::Dist_Map* tSharedNodeMap = tDistributedFactory.create_map(tSharedNodeIDs);
            mOwnedNodalValues = tDistributedFactory.create_vector(tOwnedNodeMap, 1, false, true);
            mSharedNodalValues = tDistributedFactory.create_vector(tSharedNodeMap, 1, false, true);

            // Import ADVs and assign nodal values
            this->import_advs(aOwnedADVs);

        }

        //--------------------------------------------------------------------------------------------------------------

        BSpline_Field::~BSpline_Field()
        {
            delete mOwnedNodalValues;
            delete mSharedNodalValues;
        }

        //--------------------------------------------------------------------------------------------------------------

        real BSpline_Field::get_field_value(uint aNodeIndex)
        {
            sint tNodeID = mMesh->get_glb_entity_id_from_entity_loc_index(
                    aNodeIndex,
                    EntityRank::NODE,
                    this->get_discretization_mesh_index());

            return (*mSharedNodalValues)(tNodeID);
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix<DDRMat>& BSpline_Field::get_field_sensitivities(uint aNodeIndex)
        {
            mSensitivities = trans(mMesh->get_t_matrix_of_node_loc_ind(aNodeIndex, this->get_discretization_mesh_index()));
            return mSensitivities;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDSMat> BSpline_Field::get_determining_adv_ids(uint aNodeIndex)
        {
            return mMesh->get_coefficient_IDs_of_node(aNodeIndex, this->get_discretization_mesh_index());
        }

        //--------------------------------------------------------------------------------------------------------------

        void BSpline_Field::import_advs(sol::Dist_Vector* aOwnedADVs)
        {
            // Import ADVs as usual
            Field::import_advs(aOwnedADVs);

            // Reset evaluated field
            mOwnedNodalValues->vec_put_scalar(0);

            // Evaluate field at owned nodes
            for (uint tNodeIndex = 0; tNodeIndex < mMesh->get_num_nodes(); tNodeIndex++)
            {
                if ((uint) par_rank() == mMesh->get_entity_owner(tNodeIndex, EntityRank::NODE, this->get_discretization_mesh_index()))
                {
                    sint tNodeID = mMesh->get_glb_entity_id_from_entity_loc_index(
                            tNodeIndex,
                            EntityRank::NODE,
                            this->get_discretization_mesh_index());
                    Matrix<IndexMat> tBSplineIndices = mMesh->get_coefficient_indices_of_node(tNodeIndex, this->get_discretization_mesh_index());
                    Matrix<DDRMat> tMatrix = mMesh->get_t_matrix_of_node_loc_ind(tNodeIndex, this->get_discretization_mesh_index());
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

        bool BSpline_Field::discretization_intention()
        {
            return false;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> BSpline_Field::map_to_bsplines(std::shared_ptr<Field> aField)
        {
            // Create source field
            Matrix<DDRMat> tNodalValues(mNumOriginalNodes, 1);
            for (uint tNodeIndex = 0; tNodeIndex < mNumOriginalNodes; tNodeIndex++)
            {
                tNodalValues(tNodeIndex) =
                        aField->get_field_value(tNodeIndex, mMesh->get_node_coordinate(tNodeIndex));
            }

            // Create mesh pair
            mtk::Mesh_Pair tMeshPair;
            tMeshPair.mInterpolationMesh = mMesh;
            tMeshPair.mIntegrationMesh = create_integration_mesh_from_interpolation_mesh(MeshType::HMR, mMesh);
            std::shared_ptr<mtk::Mesh_Manager> tMeshManager = std::make_shared<mtk::Mesh_Manager>();
            tMeshManager->register_mesh_pair(tMeshPair);
            mtk::Field* tField = new mtk::Field(tMeshManager, 0, this->get_discretization_mesh_index());

            // Use mapper
            mtk::Mapper tMapper;
            tField->set_nodal_values(tNodalValues);
            tMapper.perform_mapping(tField, EntityRank::NODE, EntityRank::BSPLINE);

            // Get coefficients
            Matrix<DDRMat> tCoefficients = tField->get_coefficients();

            // Clean up
            delete tMeshPair.mIntegrationMesh;
            delete tField;

            // Return mapped field
            return tCoefficients;
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
