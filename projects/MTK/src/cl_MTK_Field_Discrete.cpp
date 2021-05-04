
#include "cl_MTK_Field_Discrete.hpp"
#include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src
#include "cl_SOL_Matrix_Vector_Factory.hpp"

namespace moris
{
    namespace mtk
    {
        Field_Discrete::Field_Discrete(
                mtk::Mesh_Pair   aMeshPairs,
                uint     const & aDiscretizationMeshIndex,
                uint     const & mNumberOfFields)
                : Field(aMeshPairs,mNumberOfFields),
                  mDiscretizationMeshIndex(aDiscretizationMeshIndex)
        {
            // update coefficient data
            this->update_coefficent_data();

            mFieldIsDiscrete = true;
        }

        // ----------------------------------------------------------------------------------------------

        Field_Discrete::~Field_Discrete()
        {
            // delete distributed vectors
            delete mOwnedNodalValues;
            delete mSharedNodalValues;
        }

        // ----------------------------------------------------------------------------------------------

        void Field_Discrete::update_coefficent_data()
        {
            // get interpolation mesh
            mtk::Mesh * tIPmesh = mMeshPair.get_interpolation_mesh();

            // get number of nodes
            uint tNumberOfNodes = mNodalValues.n_rows();

            // check for correct number of nodes
            MORIS_ASSERT( tNumberOfNodes == tIPmesh->get_num_nodes(),
                    "Field_Discrete::update_coefficent_data - incorrect number of nodes.\n");

            // get maximum number of coefficients used by mesh
            mMaxNumberOfCoefficients = tIPmesh->get_max_num_coeffs_on_proc(mDiscretizationMeshIndex);

            // allocate auxiliary vectors to collect indices and IDs of used coefficients
            Matrix< DDSMat > tAllCoefIds(mMaxNumberOfCoefficients,1,-1);
            Matrix< DDSMat > tAllCoefIndices(mMaxNumberOfCoefficients,1,-1);

            // set size nodal vectors
            Matrix<DDSMat> tOwnedNodeIDs(tNumberOfNodes, 1);
            Matrix<DDSMat> tSharedNodeIDs(tNumberOfNodes, 1);

            // reset number of coefficients by used field
            mNumberOfCoefficients = 0;

            // counter for number of owned nodes
            uint tOwnedNodeCounter = 0;

            // loop over all nodes to (a) determine which coefficients are used by field by building mesh coefficient index
            // to mesh coefficient id map and (b) determine number of used coefficients
            for (uint tNodeIndex=0;tNodeIndex<mNodalValues.n_rows();++tNodeIndex)
            {
                // get node ID of current node
                tSharedNodeIDs(tNodeIndex) = tIPmesh->get_glb_entity_id_from_entity_loc_index(
                        tNodeIndex,
                        EntityRank::NODE);

                // get ownership of current node
                sint tNodeOwner = tIPmesh->get_entity_owner(tNodeIndex, EntityRank::NODE);

                if ( par_rank() == tNodeOwner )
                {
                    // store node ID in list of owned nodes
                    tOwnedNodeIDs(tOwnedNodeCounter)=tSharedNodeIDs(tNodeIndex);

                    // increment counter of owned nodes
                    tOwnedNodeCounter++;
                }

                // check whether node has an underlying discretization on this processor
                bool tNodeHasDiscretization = tIPmesh->get_mtk_vertex(tNodeIndex).has_interpolation(mDiscretizationMeshIndex);

                // process only nodes that have discretization
                if ( tNodeHasDiscretization )
                {
                    // get indices and IDs from mtk mesh - FIXME: should return const &
                    const Matrix< IndexMat > tCoefIndices = tIPmesh->get_coefficient_indices_of_node(
                            tNodeIndex,
                            mDiscretizationMeshIndex);

                    const Matrix< IdMat > tCoefIds = tIPmesh->get_coefficient_IDs_of_node(
                            tNodeIndex,
                            mDiscretizationMeshIndex);

                    // check that number of indices and ids are the same
                    MORIS_ASSERT( tCoefIds.numel() == tCoefIndices.numel(),
                            "Field_Discrete::update_coefficent_data - numbers of coefficients and ids do not match.\n");

                    // get number of coefficients for current node
                    uint tNumCoefOfNode=tCoefIds.numel();

                    for (uint tCoefIndex=0;tCoefIndex<tNumCoefOfNode;++tCoefIndex)
                    {
                        // get coefficient index
                        moris_index tCurrentIndex = tCoefIndices(tCoefIndex);

                        // check whether mesh coefficient has already been set
                        if ( tAllCoefIds(tCurrentIndex) == -1 )
                        {
                            // increase field coefficient count
                            mNumberOfCoefficients++;

                            // populate mesh index to mesh coefficient id map
                            tAllCoefIds(tCurrentIndex) = tCoefIds(tCoefIndex);
                        }
                        else
                        {
                            // check for consistency
                            MORIS_ASSERT( tAllCoefIds(tCurrentIndex) == tCoefIds(tCoefIndex),
                                    "Field_Discrete::update_coefficent_data - inconsistent index and ids.\n");
                        }
                    }
                }
            }

            // resize list of owned nodes
            tOwnedNodeIDs.resize(tOwnedNodeCounter,1);

            // allocate size for map from mesh indices to field indices
            mMeshToFieldCoefficientIndexMap.set_size(mMaxNumberOfCoefficients,1);

            // allocate size for map from field indices to mesh ID and ownership
            mFieldIndexToMeshCoefficientIdAndOwnerMap.set_size(mNumberOfCoefficients,2);

            // initialize field mesh index counter
            uint tCounter=0;

            // loop over all mesh coefficient indices
            for (uint tCoefIndex=0;tCoefIndex<mMaxNumberOfCoefficients;++tCoefIndex)
            {
                if ( tAllCoefIds(tCoefIndex) > -1 )
                {
                    // assign field coefficient index to mesh index
                    mMeshToFieldCoefficientIndexMap(tCoefIndex)=tCounter;

                    // assign mesh coefficient id to field index
                    mFieldIndexToMeshCoefficientIdAndOwnerMap(tCounter,0)=tAllCoefIds(tCoefIndex);

                    // get ownership of mesh coefficient
                    mFieldIndexToMeshCoefficientIdAndOwnerMap(tCounter,1)=tIPmesh->get_entity_owner(
                            tCoefIndex,
                            EntityRank::BSPLINE,
                            mDiscretizationMeshIndex);

                    // increase counter of field indices
                    tCounter++;
                }
            }

            // check that maps for all used coefficients have been populated
            MORIS_ASSERT ( mNumberOfCoefficients == (sint)tCounter,
                    "Field_Discrete::update_coefficent_data - inconsistent map counter.\n");

            // create distributed vectors for node value computation
            sol::Matrix_Vector_Factory tDistributedFactory;

            sol::Dist_Map* tOwnedNodeMap  = tDistributedFactory.create_map(tOwnedNodeIDs);
            sol::Dist_Map* tSharedNodeMap = tDistributedFactory.create_map(tSharedNodeIDs);

            mOwnedNodalValues  = tDistributedFactory.create_vector(tOwnedNodeMap,  1, false, true);
            mSharedNodalValues = tDistributedFactory.create_vector(tSharedNodeMap, 1, false, true);
        }

        // ----------------------------------------------------------------------------------------------

        uint Field_Discrete::get_discretization_order() const
        {
            // get interpolation mesh
            mtk::Mesh * tIPmesh = mMeshPair.get_interpolation_mesh();

            // return discretization order
            return tIPmesh->get_discretization_order( this->get_discretization_mesh_index() );
        }

        // ----------------------------------------------------------------------------------------------

        moris_index Field_Discrete::get_discretization_mesh_index() const
        {
            MORIS_ASSERT( mDiscretizationMeshIndex > -1,
                    "Field_Discrete::get_discretization_mesh_index - discretization mesh index not set.\n");

            return mDiscretizationMeshIndex;
        }

        //------------------------------------------------------------------------------

        void Field_Discrete::set_coefficient_vector(const Matrix< DDRMat > & aCoefficients)
        {
            // set coefficient initialization flag to true
            mCoefficientsAreInitialized = true;

            // copy coefficients
            mCoefficients = aCoefficients;

            mNumberOfCoefficients = aCoefficients.n_rows();
        }

        // ----------------------------------------------------------------------------------------------

        const Matrix<IdMat> & Field_Discrete::get_coefficient_id_and_owner_vector()
        {
            // check for proper size of mesh coefficient ID vector
            MORIS_ASSERT( (sint)mFieldIndexToMeshCoefficientIdAndOwnerMap.n_rows() == mNumberOfCoefficients,
                    "Field_Discrete::get_coefficient_id_vector - mesh coefficient ID vector has incorrect size.\n");

            // return mesh coefficient ID vector
            return mFieldIndexToMeshCoefficientIdAndOwnerMap;
        }

        // ----------------------------------------------------------------------------------------------

        void Field_Discrete::compute_nodal_values()
        {
            MORIS_ASSERT( mDiscretizationMeshIndex > -1,
                    "Field_Discrete::compute_nodal_values - discretization mesh index not set.\n");

            // get interpolation mesh
            mtk::Mesh * tIPmesh = mMeshPair.get_interpolation_mesh();

            // check that number of coefficient on mesh matches size of coefficient vector
            MORIS_ASSERT( (sint)mCoefficients.n_rows() == mNumberOfCoefficients,
                    "Field_Discrete::compute_nodal_values - number of coefficient on mesh does not match size of coefficient vector.\n");

            // check that coefficient vector has been initialized
            MORIS_ASSERT( mCoefficientsAreInitialized,
                    "Field_Discrete::compute_nodal_values - coefficient vector has not been initialized.\n");

            // make sure that nodal value matrix is properly sized
            mNodalValues.resize( tIPmesh->get_num_nodes(), mNumberOfFields );

            //loop over all nodes
            for (uint tNodeIndex=0;tNodeIndex<mNodalValues.n_rows();++tNodeIndex)
            {
                // get node owner
                sint tNodeOwner = tIPmesh->get_entity_owner(tNodeIndex, EntityRank::NODE );

                // process only owned nodes
                if ( par_rank() == tNodeOwner )
                {
                    MORIS_ASSERT ( tIPmesh->get_mtk_vertex(tNodeIndex).has_interpolation(mDiscretizationMeshIndex),
                            "Field_Discrete::compute_nodal_values - owned node with index %d does not have discretization.\n",
                            tNodeIndex );

                    // get t-matrix and coefficient indices - FIXME: should receive const reference
                    const Matrix<IndexMat> tMeshCoefIndices = tIPmesh->get_coefficient_indices_of_node(
                            tNodeIndex,
                            mDiscretizationMeshIndex);

                    const Matrix<DDRMat> tMatrix = tIPmesh->get_t_matrix_of_node_loc_ind(
                            tNodeIndex,
                            mDiscretizationMeshIndex);

                    // compute nodal value : t-matrix * coefficients
                    real tValue=0.0;

                    for (uint tCoef=0;tCoef<tMeshCoefIndices.numel();++tCoef)
                    {
                        // get field coefficient index
                        moris_index tFieldCoefficientIndex =
                                mMeshToFieldCoefficientIndexMap( tMeshCoefIndices(tCoef) );

                        // compute contribution to nodal value
                        tValue += tMatrix(tCoef) * mCoefficients( tFieldCoefficientIndex );
                    }

                    // get node ID
                    sint tNodeID = tIPmesh->get_glb_entity_id_from_entity_loc_index(
                                                tNodeIndex,
                                                EntityRank::NODE );

                    // copy nodal value on distributed vector
                    (*mOwnedNodalValues)(tNodeID) = tValue;
                }
            }

            // get values for shared and owned nodes
            mSharedNodalValues->import_local_to_global(*mOwnedNodalValues);

            // copy shared and own nodal values onto local nodal field
            for (uint tNodeIndex=0;tNodeIndex<mNodalValues.n_rows();++tNodeIndex)
            {
                // get node ID
                 sint tNodeID = tIPmesh->get_glb_entity_id_from_entity_loc_index(
                                             tNodeIndex,
                                             EntityRank::NODE );

                 // extract nodal value
                 real tValue = (*mSharedNodalValues)(tNodeID);

                // apply nodal value to all fields
                for (uint tFieldIndex=0;tFieldIndex<mNumberOfFields;++tFieldIndex)
                {
                    mNodalValues(tNodeIndex,tFieldIndex)=tValue;
                }
            }

            mUpdateNodalValues=false;
        }

        // ----------------------------------------------------------------------------------------------

        void Field_Discrete::compute_derivatives_of_field_value(
                Matrix< DDRMat >       & aDerivatives,
                Matrix< IndexMat >     & aCoefIndices,
                uint             const & aNodeIndex,
                uint             const & aFieldIndex)
        {
            // check that discretization index is valid
            MORIS_ASSERT( mDiscretizationMeshIndex > -1,
                    "Field_Discrete::compute_nodal_values - discretization mesh index not set.\n");

            // get interpolation mesh
            mtk::Mesh * tIPmesh = mMeshPair.get_interpolation_mesh();

            // check that node has an underlying discretization; return zero size vectors if
            // discretization cannot be accessed
            if ( ! tIPmesh->get_mtk_vertex(aNodeIndex).has_interpolation(mDiscretizationMeshIndex) )
            {
                aDerivatives.set_size(0,0);
                aCoefIndices.set_size(0,0);
                return;
            }

            // get coefficient indices for requested node
            Matrix < IndexMat > tMeshCoefficientIndices = tIPmesh->get_coefficient_indices_of_node(
                    aNodeIndex,
                    mDiscretizationMeshIndex);

            // get number of nodal coefficients
            uint tNumberOfNodalCoefficients = tMeshCoefficientIndices.numel();

            // set size of coefficient index vector to be returned
            aCoefIndices.set_size(tNumberOfNodalCoefficients,1);

            // map mesh coefficient indices onto field coefficient indices
            for (uint tIndex=0;tIndex<tNumberOfNodalCoefficients;++tIndex)
            {
                // extract mesh index
                moris_index tMeshIndex = tMeshCoefficientIndices(tIndex);

                // get field index for given mesh index
                aCoefIndices(tIndex) = mMeshToFieldCoefficientIndexMap( tMeshIndex );
            }

            // get and return t-matrix coefficients for requested node
            aDerivatives = tIPmesh->get_t_matrix_of_node_loc_ind(
                    aNodeIndex,
                    mDiscretizationMeshIndex);
        }

        // ----------------------------------------------------------------------------------------------

    }
}
