#include "cl_GEN_BSpline_Field.hpp"
#include "st_MTK_Mesh_Pair.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Mesh_Factory.hpp"
#include "cl_MTK_Field.hpp"
#include "cl_MTK_Field_Discrete.hpp"
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
                const Matrix<DDUMat>&    aCoefficientIndices,
                const Matrix<DDSMat>&    aSharedADVIds,
                uint                     aADVOffsetID,
                mtk::Interpolation_Mesh* aMesh,
                std::shared_ptr<Field>   aField)
                : Field(aCoefficientIndices, aSharedADVIds, aField)
                , Field_Discrete_Integration(aMesh->get_num_nodes())
                , mSharedADVIds(aSharedADVIds)
                , mMesh(aMesh)
        {
            // Map to B-splines
            Matrix<DDRMat> tTargetField = this->map_to_bsplines(aField);
            MORIS_ERROR(tTargetField.length() == aCoefficientIndices.length(),
                    "MTK mapper is reporting a different number of coefficients than the mesh at the finest level.");

            // Get B-spline mesh index
            uint tDiscretizationMeshIndex = this->get_discretization_mesh_index();

            // Assign ADVs
            for (uint tCoefficient = 0; tCoefficient < aCoefficientIndices.length(); tCoefficient++)
            {
                uint tCoefficientIndex = aCoefficientIndices(tCoefficient);
                if ((uint) par_rank() == mMesh->get_entity_owner(tCoefficientIndex, EntityRank::BSPLINE, tDiscretizationMeshIndex))
                {
                    // Calculate ADV ID using offset
                    sint tADVId = aADVOffsetID + aMesh->get_glb_entity_id_from_entity_loc_index(
                            tCoefficientIndex,
                            EntityRank::BSPLINE,
                            tDiscretizationMeshIndex);

                    // Assign distributed vector element based on ID
                    (*aOwnedADVs)(tADVId) = tTargetField(tCoefficientIndex);
                }
            }

            // Determine number of owned and shared node IDs
            uint tOwnedNodeCount = 0;
            for (uint tNodeIndex = 0; tNodeIndex < mMesh->get_num_nodes(); tNodeIndex++)
            {
                if ((uint) par_rank() == mMesh->get_entity_owner(tNodeIndex, EntityRank::NODE, tDiscretizationMeshIndex))
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
                        tDiscretizationMeshIndex);
                tSharedNodeIDs(tNodeIndex) = tNodeID;
                if ((uint) par_rank() == mMesh->get_entity_owner(tNodeIndex, EntityRank::NODE, tDiscretizationMeshIndex))
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
            return trans( mMesh->get_coefficient_IDs_of_node(aNodeIndex, this->get_discretization_mesh_index()) );
        }

        //--------------------------------------------------------------------------------------------------------------

        void BSpline_Field::import_advs(sol::Dist_Vector* aOwnedADVs)
        {
            // Import ADVs as usual
            Field::import_advs(aOwnedADVs);

            // Reset evaluated field
            mOwnedNodalValues->vec_put_scalar(0);

            moris::Matrix<DDRMat> tCoeff(mFieldVariables.size());

            for( uint Ik = 0; Ik<mFieldVariables.size(); Ik++)
            {
                tCoeff( Ik ) = *mFieldVariables( Ik );
            }

            // Create mesh pair
            mtk::Mesh_Pair tMeshPair(mMesh, create_integration_mesh_from_interpolation_mesh(MeshType::HMR, mMesh));

            mtk::Field_Discrete* tField = new mtk::Field_Discrete( tMeshPair, this->get_discretization_mesh_index());

            // Use mapper
            mtk::Mapper tMapper;
            tField->unlock_field();
            tField->set_coefficients(tCoeff);
            tMapper.perform_mapping(tField, EntityRank::BSPLINE, EntityRank::NODE);

            // Get coefficients
            Matrix<DDRMat> tNodalValues = tField->get_nodal_values();

            // Clean up
            delete tMeshPair.get_integration_mesh();
            delete tField;
            //-----------------------------------------------------

            // Evaluate field at owned nodes
            for (uint tNodeIndex = 0; tNodeIndex < mMesh->get_num_nodes(); tNodeIndex++)
            {
                if ((uint) par_rank() == mMesh->get_entity_owner(tNodeIndex, EntityRank::NODE ))
                {
                    sint tNodeID = mMesh->get_glb_entity_id_from_entity_loc_index(
                            tNodeIndex,
                            EntityRank::NODE );
                    MORIS_ASSERT(tNodalValues( tNodeIndex ) != MORIS_REAL_MAX, "value is MORIS_REAL_MAX, check mapper");
                        (*mOwnedNodalValues)(tNodeID) = tNodalValues( tNodeIndex );
                }
            }

            // Global assembly
            mOwnedNodalValues->vector_global_assembly();

            // Import nodal values
            mSharedNodalValues->import_local_to_global(*mOwnedNodalValues);
        }

        //--------------------------------------------------------------------------------------------------------------

        mtk::Interpolation_Mesh* BSpline_Field::get_mesh()
        {
            return mMesh;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> BSpline_Field::map_to_bsplines(std::shared_ptr<Field> aField)
        {
            // Mapper
            mtk::Mapper tMapper;

            // Output field
            mtk::Mesh_Pair tOutputMeshPair(mMesh, create_integration_mesh_from_interpolation_mesh(MeshType::HMR, mMesh));

            mtk::Field_Discrete* tOutputField = new mtk::Field_Discrete( tOutputMeshPair, this->get_discretization_mesh_index());

            // Input mesh
            mtk::Interpolation_Mesh* tInputMesh = aField->get_mesh();
            if (not tInputMesh)
            {
                tInputMesh = mMesh;
            }

            // Nodal values
            Matrix<DDRMat> tNodalValues(tInputMesh->get_num_nodes(), 1);
            for (uint tNodeIndex = 0; tNodeIndex < tInputMesh->get_num_nodes(); tNodeIndex++)
            {
                tNodalValues(tNodeIndex) =
                        aField->get_field_value(tNodeIndex, mMesh->get_node_coordinate(tNodeIndex));
            }

            // Interpolation not needed
            if (tInputMesh == mMesh)
            {
                tOutputField->unlock_field();
                tOutputField->set_nodal_values(tNodalValues);
                tMapper.map_input_field_to_output_field_2(tOutputField);
            }
            // Interpolation needed
            else
            {
                // Input field
                mtk::Mesh_Pair tInputMeshPair(tInputMesh, create_integration_mesh_from_interpolation_mesh(MeshType::HMR, tInputMesh));

                mtk::Field_Discrete* tInputField = new mtk::Field_Discrete(
                        tInputMeshPair,
                        aField->get_discretization_mesh_index());

                // Do interpolation
                tInputField->unlock_field();
                tInputField->set_nodal_values(tNodalValues);
                tMapper.map_input_field_to_output_field(tInputField, tOutputField);

                delete tInputMeshPair.get_integration_mesh();
            }

            // Get coefficients
            Matrix<DDRMat> tCoefficients = tOutputField->get_coefficients();

            // Clean up
            delete tOutputMeshPair.get_integration_mesh();
            delete tOutputField;

            // Return mapped field
            return tCoefficients;
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
