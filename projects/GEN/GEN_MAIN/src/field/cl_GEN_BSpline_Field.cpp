#include "cl_GEN_BSpline_Field.hpp"
#include "cl_MTK_Mesh_Pair.hpp"
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
                sol::Dist_Vector*      aOwnedADVs,
                const Matrix<DDUMat>&  aCoefficientIndices,
                const Matrix<DDSMat>&  aSharedADVIds,
                uint                   aADVOffsetID,
                mtk::Mesh_Pair         aMeshPair,
                std::shared_ptr<Field> aField)
                : Field(aCoefficientIndices, aSharedADVIds, aMeshPair, aField)
                , Field_Discrete_Integration(aMeshPair.get_interpolation_mesh()->get_num_nodes())
                , mSharedADVIds(aSharedADVIds)
                , mADVOffsetID(aADVOffsetID)
                , mMeshPair(aMeshPair)
        {
            // Map to B-splines
            Matrix<DDRMat> tTargetField = this->map_to_bsplines(aField);

            this->distribute_coeffs(
                    tTargetField,
                    aOwnedADVs,
                    aCoefficientIndices,
                    aSharedADVIds,
                    aADVOffsetID);

            mFieldIsDiscrete = true;
        }

        //--------------------------------------------------------------------------------------------------------------

        BSpline_Field::BSpline_Field(
                sol::Dist_Vector*      aOwnedADVs,
                const Matrix<DDUMat>&  aCoefficientIndices,
                const Matrix<DDSMat>&  aSharedADVIds,
                uint                   aADVOffsetID,
                mtk::Mesh_Pair         aMeshPair,
                std::shared_ptr<Field> aField,
                std::shared_ptr<mtk::Field> aMTKField )
                : Field(aCoefficientIndices, aSharedADVIds, aMeshPair, aField)
                , Field_Discrete_Integration(aMeshPair.get_interpolation_mesh()->get_num_nodes())
                , mSharedADVIds(aSharedADVIds)
                , mADVOffsetID(aADVOffsetID)
                , mMeshPair(aMeshPair)
        {
            // Map to B-splines
            const Matrix<DDRMat> & tTargetField = aMTKField->get_coefficients();

            this->distribute_coeffs(
                    tTargetField,
                    aOwnedADVs,
                    aCoefficientIndices,
                    aSharedADVIds,
                    aADVOffsetID);

            mFieldIsDiscrete = true;
        }

        //--------------------------------------------------------------------------------------------------------------

        BSpline_Field::~BSpline_Field()
        {
            delete mOwnedNodalValues;
            delete mSharedNodalValues;
        }
        //--------------------------------------------------------------------------------------------------------------

        void BSpline_Field::distribute_coeffs(
                const Matrix<DDRMat> & aTargetField,
                sol::Dist_Vector*      aOwnedADVs,
                const Matrix<DDUMat>&  aCoefficientIndices,
                const Matrix<DDSMat>&  aSharedADVIds,
                uint                   aADVOffsetID)
        {
            MORIS_ERROR(aTargetField.length() == aCoefficientIndices.length(),
                    "MTK mapper is reporting a different number of coefficients than the mesh at the finest level. %-5i | %-5i",
                    aTargetField.length(),
                    aCoefficientIndices.length() );

            // Get B-spline mesh index
            mtk::Mesh* tMesh = mMeshPair.get_interpolation_mesh();
            uint tDiscretizationMeshIndex = this->get_discretization_mesh_index();

            // Assign ADVs
            for (uint tCoefficient = 0; tCoefficient < aCoefficientIndices.length(); tCoefficient++)
            {
                uint tCoefficientIndex = aCoefficientIndices(tCoefficient);
                if ((uint) par_rank() == tMesh->get_entity_owner(tCoefficientIndex, EntityRank::BSPLINE, tDiscretizationMeshIndex))
                {
                    // Calculate ADV ID using offset
                    sint tADVId = mADVOffsetID + tMesh->get_glb_entity_id_from_entity_loc_index(
                            tCoefficientIndex,
                            EntityRank::BSPLINE,
                            tDiscretizationMeshIndex);

                    // Assign distributed vector element based on ID
                    (*aOwnedADVs)(tADVId) = aTargetField(tCoefficientIndex);
                }
            }

            // Determine number of owned and shared node IDs
            uint tOwnedNodeCount = 0;
            for (uint tNodeIndex = 0; tNodeIndex < tMesh->get_num_nodes(); tNodeIndex++)
            {
                if ((uint) par_rank() == tMesh->get_entity_owner(tNodeIndex, EntityRank::NODE, tDiscretizationMeshIndex))
                {
                    tOwnedNodeCount++;
                }
            }
            Matrix<DDSMat> tOwnedNodeIDs(tOwnedNodeCount, 1);
            Matrix<DDSMat> tSharedNodeIDs(tMesh->get_num_nodes(), 1);

            // Assign owned and shared node IDs
            tOwnedNodeCount = 0;
            for (uint tNodeIndex = 0; tNodeIndex < tMesh->get_num_nodes(); tNodeIndex++)
            {
                sint tNodeID = tMesh->get_glb_entity_id_from_entity_loc_index(
                        tNodeIndex,
                        EntityRank::NODE,
                        tDiscretizationMeshIndex);
                tSharedNodeIDs(tNodeIndex) = tNodeID;
                if ((uint) par_rank() == tMesh->get_entity_owner(tNodeIndex, EntityRank::NODE, tDiscretizationMeshIndex))
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

        real BSpline_Field::get_field_value(uint aNodeIndex)
        {
            mtk::Mesh* tMesh = mMeshPair.get_interpolation_mesh();
            sint tNodeID = tMesh->get_glb_entity_id_from_entity_loc_index(
                    aNodeIndex,
                    EntityRank::NODE,
                    this->get_discretization_mesh_index());

            return (*mSharedNodalValues)(tNodeID);
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix<DDRMat>& BSpline_Field::get_dfield_dadvs(uint aNodeIndex)
        {
            mtk::Mesh* tMesh = mMeshPair.get_interpolation_mesh();
            mSensitivities = trans(tMesh->get_t_matrix_of_node_loc_ind(aNodeIndex, this->get_discretization_mesh_index()));
            return mSensitivities;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDSMat> BSpline_Field::get_determining_adv_ids(uint aNodeIndex)
        {
            mtk::Mesh* tMesh = mMeshPair.get_interpolation_mesh();
            return trans( tMesh->get_coefficient_IDs_of_node(aNodeIndex, this->get_discretization_mesh_index()) ) + mADVOffsetID;
        }

        void BSpline_Field::get_coefficient_vector()
        {
            uint tNumCoeffs = mFieldVariables.size();

            mCoefficients.set_size( tNumCoeffs, 1 );

            for( uint Ik = 0; Ik<tNumCoeffs; Ik++)
            {
                mCoefficients( Ik ) = *mFieldVariables( Ik );
            }
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

            // Create field
            mtk::Field_Discrete* tField = new mtk::Field_Discrete(mMeshPair, this->get_discretization_mesh_index());

            // Use mapper
            mtk::Mapper tMapper;
            tField->unlock_field();
            tField->set_coefficients(tCoeff);
            tMapper.perform_mapping(tField, EntityRank::BSPLINE, EntityRank::NODE);

            // Get coefficients
            Matrix<DDRMat> tNodalValues = tField->get_values();
            this->unlock_field();
            this->set_values(tNodalValues);

            // Clean up
            delete tField;
            //-----------------------------------------------------

            // Evaluate field at owned nodes
            mtk::Mesh* tMesh = mMeshPair.get_interpolation_mesh();
            for (uint tNodeIndex = 0; tNodeIndex < tMesh->get_num_nodes(); tNodeIndex++)
            {
                if ((uint) par_rank() == tMesh->get_entity_owner(tNodeIndex, EntityRank::NODE ))
                {
                    sint tNodeID = tMesh->get_glb_entity_id_from_entity_loc_index(
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

        mtk::Mesh_Pair BSpline_Field::get_mesh_pair()
        {
            return mMeshPair;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> BSpline_Field::map_to_bsplines(std::shared_ptr<Field> aField)
        {
            // Mapper
            mtk::Mapper tMapper;

            // New mesh
            mtk::Interpolation_Mesh* tMesh = mMeshPair.get_interpolation_mesh();

            // Output field
            mtk::Field_Discrete* tOutputField = new mtk::Field_Discrete(mMeshPair, this->get_discretization_mesh_index());

            // Nodal values
            Matrix<DDRMat> tNodalValues(tMesh->get_num_nodes(), 1);
            for (uint tNodeIndex = 0; tNodeIndex < tMesh->get_num_nodes(); tNodeIndex++)
            {
                tNodalValues(tNodeIndex) =
                        aField->get_field_value(tNodeIndex, tMesh->get_node_coordinate(tNodeIndex));
            }

            tOutputField->unlock_field();
            tOutputField->set_values(tNodalValues);
            tMapper.map_input_field_to_output_field_2(tOutputField);

            // Get coefficients
            Matrix<DDRMat> tCoefficients = tOutputField->get_coefficients();

            // Clean up
            delete tOutputField;

            // Return mapped field
            return tCoefficients;
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
