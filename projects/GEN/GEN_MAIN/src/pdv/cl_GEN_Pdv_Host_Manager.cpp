#include "cl_GEN_Pdv_Host_Manager.hpp"
#include "cl_SOL_Matrix_Vector_Factory.hpp"
#include "cl_Communication_Tools.hpp"
#include "fn_trans.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Pdv_Host_Manager::Pdv_Host_Manager()
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        Pdv_Host_Manager::~Pdv_Host_Manager()
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::set_owned_adv_ids(Matrix<DDSMat> aOwnedADVIds)
        {
            mOwnedADVIds = aOwnedADVIds;
            mADVIdsSet = true;
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::set_communication_table( const Matrix< IdMat > & aCommTable )
        {
            mCommTable = aCommTable;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix< IdMat > Pdv_Host_Manager::get_communication_table() // FIXME
        {
            return mCommTable;
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::set_num_background_nodes(uint aNumNodes)
        {
            mIntersectionNodes.resize(aNumNodes);
            mNumBackgroundNodesSet = true;
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::reset()
        {
            mIpPdvHosts.clear();
            mIntersectionNodes.clear();
            mOwnedPdvLocalToGlobalMap.resize(0, 0);
            mOwnedAndSharedPdvLocalToGlobalMap.resize(0, 0);
            mNumOwnedPdvs = 0;
            mNumOwnedAndSharedPdvs = 0;
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::get_ip_dv_types_for_set(const moris::moris_index aIPMeshSetIndex, Cell<Cell<PDV_Type>>& aPdvTypes)
        {
            if (mIpPdvTypes.size() > 0) // FIXME
            {
                aPdvTypes = mIpPdvTypes(aIPMeshSetIndex);
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::get_ig_dv_types_for_set(const moris::moris_index aIGMeshSetIndex, Cell<Cell<PDV_Type>>& aPdvTypes)
        {
            if (mIgPdvTypes.size() > 0) // FIXME
            {
                aPdvTypes = mIgPdvTypes(aIGMeshSetIndex);
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::get_ip_unique_dv_types_for_set(const moris_index aIPMeshSetIndex, Cell<PDV_Type>& aPdvTypes)
        {
            if (mUniqueIpPdvTypes.size() > 0) // FIXME
            {
                aPdvTypes =  mUniqueIpPdvTypes(aIPMeshSetIndex);
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::get_ig_unique_dv_types_for_set(const moris::moris_index aIGMeshSetIndex, Cell<PDV_Type>& aPdvTypes)
        {
            if (mUniqueIgPdvTypes.size() > 0) // FIXME
            {
                aPdvTypes =  mUniqueIgPdvTypes(aIGMeshSetIndex);
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::get_ip_pdv_value(
                const Matrix<IndexMat> & aNodeIndices,
                const Cell<PDV_Type>   & aPdvTypes,
                Cell<Matrix<DDRMat>>   & aDvValues )
        {
            // Get the number of node indices requested
            uint tNumIndices = aNodeIndices.numel();

            // Get the number of dv types requested
            uint tNumTypes = aPdvTypes.size();

            // Cell size
            aDvValues.resize(tNumTypes);

            // loop over the node indices
            for (uint tPdvTypeIndex = 0; tPdvTypeIndex < tNumTypes; tPdvTypeIndex++)
            {
                // Matrix size
                aDvValues(tPdvTypeIndex).set_size(tNumIndices, 1);

                // loop over the requested dv types
                for (uint tNodeIndex = 0; tNodeIndex < tNumIndices; tNodeIndex++)
                {
                    aDvValues(tPdvTypeIndex)(tNodeIndex) =
                            mIpPdvHosts(aNodeIndices(tNodeIndex))->get_pdv_value(aPdvTypes(tPdvTypeIndex));
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::get_ip_pdv_value(const Matrix<IndexMat>& aNodeIndices,
                const Cell<PDV_Type>&   aPdvTypes,
                Cell<Matrix<DDRMat>>&   aDvValues,
                Cell<Matrix<DDSMat>>&   aIsActiveDv)
        {
            // Get the number of node indices requested
            uint tNumIndices = aNodeIndices.length();

            // Get the number of dv types requested
            uint tNumTypes = aPdvTypes.size();

            // Cell size
            aDvValues.resize(tNumTypes);
            aIsActiveDv.resize(tNumTypes);

            // loop over the node indices
            for (uint tPdvTypeIndex = 0; tPdvTypeIndex < tNumTypes; tPdvTypeIndex++)
            {
                // Matrix size
                aDvValues(tPdvTypeIndex).set_size(tNumIndices, 1);
                aIsActiveDv(tPdvTypeIndex).set_size(tNumIndices, 1);

                // loop over the requested dv types
                for (uint tNodeIndex = 0; tNodeIndex < tNumIndices; tNodeIndex++)
                {
                    aDvValues(tPdvTypeIndex)(tNodeIndex) =
                            mIpPdvHosts(aNodeIndices(tNodeIndex))->get_pdv_value(aPdvTypes(tPdvTypeIndex));
                    aIsActiveDv(tPdvTypeIndex)(tNodeIndex) =
                            mIpPdvHosts(aNodeIndices(tNodeIndex))->is_active_type(aPdvTypes(tPdvTypeIndex));
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::get_ig_pdv_value(const Matrix<IndexMat>& aNodeIndices,
                const Cell<PDV_Type>&   aPdvTypes,
                Cell<Matrix<DDRMat>>&   aDvValues)
        {
            // Get the number of node indices requested
            uint tNumIndices = aNodeIndices.length();

            // Get the number of dv types requested
            uint tNumTypes = aPdvTypes.size();

            // Cell size
            aDvValues.resize(tNumTypes);

            // loop over the node indices
            for (uint tPdvTypeIndex = 0; tPdvTypeIndex < tNumTypes; tPdvTypeIndex++)
            {
                // Matrix size
                aDvValues(tPdvTypeIndex).set_size(tNumIndices, 1);

                // loop over the requested dv types
                for (uint tNode = 0; tNode < tNumIndices; tNode++)
                {
                    if (mIntersectionNodes(aNodeIndices(tNode)))
                    {
                        aDvValues(tPdvTypeIndex)(tNode) =
                                mIntersectionNodes(aNodeIndices(tNode))->get_coordinate_value(static_cast<uint>(aPdvTypes(tPdvTypeIndex)));
                    }
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::get_ig_pdv_value(const Matrix<IndexMat>& aNodeIndices,
                const Cell<PDV_Type>&   aPdvTypes,
                Cell<Matrix<DDRMat>>&   aDvValues,
                Cell<Matrix<DDSMat>>&   aIsActiveDv )
        {
            // Get the number of node indices requested
            uint tNumIndices = aNodeIndices.length();

            // Get the number of dv types requested
            uint tNumTypes = aPdvTypes.size();

            // Cell size
            aDvValues.resize(tNumTypes);
            aIsActiveDv.resize(tNumTypes);

            // loop over the node indices
            for (uint tPdvTypeIndex = 0; tPdvTypeIndex < tNumTypes; tPdvTypeIndex++)
            {
                // Matrix size
                aDvValues(tPdvTypeIndex).set_size(tNumIndices, 1);
                aIsActiveDv(tPdvTypeIndex).set_size(tNumIndices, 1, 0);

                // loop over the requested dv types
                for (uint tNode = 0; tNode < tNumIndices; tNode++)
                {
                    if (mIntersectionNodes(aNodeIndices(tNode)))
                    {
                        aDvValues(tPdvTypeIndex)(tNode) =
                                mIntersectionNodes(aNodeIndices(tNode))->get_coordinate_value(static_cast<uint>(aPdvTypes(tPdvTypeIndex)));
                        aIsActiveDv(tPdvTypeIndex)(tNode) = 1;
                    }
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix<DDSMat> & Pdv_Host_Manager::get_my_local_global_map()
        {
            return mOwnedPdvLocalToGlobalMap;
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix<DDSMat> & Pdv_Host_Manager::get_my_local_global_overlapping_map()
        {
            return mOwnedAndSharedPdvLocalToGlobalMap;
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::get_ip_dv_ids_for_type_and_ind(const Matrix<IndexMat>& aNodeIndices,
                const Cell<PDV_Type>&   aPdvTypes,
                Cell<Matrix<IdMat>>&    aDvIds)
        {
            // get the number of node indices requested
            uint tNumIndices = aNodeIndices.length();

            // get the number of dv types requested
            uint tNumTypes = aPdvTypes.size();

            // set size for list of dv values
            aDvIds.resize(tNumTypes);

            // loop over the node indices
            for (uint tNode = 0; tNode < tNumIndices; tNode++)
            {
                // loop over the requested dv types
                for (uint tPdvTypeIndex = 0; tPdvTypeIndex < tNumTypes; tPdvTypeIndex++)
                {
                    aDvIds(tPdvTypeIndex).resize(tNumIndices, 1);
                    aDvIds(tPdvTypeIndex)(tNode) = mIpPdvHosts(aNodeIndices(tNode))->get_global_index_for_pdv_type(aPdvTypes(tPdvTypeIndex));
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::get_ig_dv_ids_for_type_and_ind(const Matrix<IndexMat>& aNodeIndices,
                const Cell<PDV_Type>&   aPdvTypes,
                Cell<Matrix<IdMat>>&    aDvIds)
        {
            // get the number of node indices requested
            uint tNumIndices = aNodeIndices.numel();

            // get the number of dv types requested
            uint tNumTypes = aPdvTypes.size();

            // set size for list of dv values
            aDvIds.resize( tNumTypes );

            // loop over the requested dv types
            for ( uint tPdvTypeIndex = 0; tPdvTypeIndex < tNumTypes; tPdvTypeIndex++ )
            {
                aDvIds(tPdvTypeIndex).set_size(tNumIndices, 1, -1);

                // loop over the node indices
                for ( uint tNode = 0; tNode < tNumIndices; tNode++ )
                {
                    if (mIntersectionNodes(aNodeIndices(tNode)))
                    {
                        aDvIds(tPdvTypeIndex)(tNode) =
                                mIntersectionNodes(aNodeIndices(tNode))->get_starting_pdv_id()
                                + static_cast<uint>(aPdvTypes(tPdvTypeIndex));
                    }
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::get_ip_requested_dv_types( Cell< PDV_Type > & aPdvTypes )
        {
            aPdvTypes =  mRequestedIpPdvTypes;
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::get_ig_requested_dv_types( Cell< PDV_Type > & aPdvTypes )
        {
            aPdvTypes =  mRequestedIgPdvTypes;
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::create_interpolation_pdv_hosts(
                const Cell<Matrix<DDSMat>>       & aNodeIndicesPerSet,
                const Cell<Matrix<DDSMat>>       & aNodeIdsPerSet,
                const Cell<Matrix<DDSMat>>       & aNodeOwnersPerSet,
                const Cell<Matrix<DDRMat>>       & aNodeCoordinates,
                const Cell<Cell<Cell<PDV_Type>>> & aPdvTypes)
        {
            // Check that number of sets is consistent
            uint tNumSets = aPdvTypes.size();
            MORIS_ERROR(tNumSets == aNodeIndicesPerSet.size(),
                    "Information passed to Pdv_Host_Manager.create_interpolation_pdv_hosts() does not have a consistent number of sets!");

            // Set PDV types
            mIpPdvTypes = aPdvTypes;
            mUniqueIpPdvTypes.resize(tNumSets);

            // Initialize PDV hosts
            mIpPdvHosts.resize(aNodeCoordinates.size(), nullptr);

            // Create PDV hosts
            for (uint tMeshSetIndex = 0; tMeshSetIndex < tNumSets; tMeshSetIndex++)
            {
                // Get number of unique PDVs
                uint tNumUniquePdvs = 0;
                for (uint tGroupIndex = 0; tGroupIndex < mIpPdvTypes(tMeshSetIndex).size(); tGroupIndex++)
                {
                    tNumUniquePdvs += mIpPdvTypes(tMeshSetIndex)(tGroupIndex).size();
                }
                mUniqueIpPdvTypes(tMeshSetIndex).resize(tNumUniquePdvs);

                // Copy PDV types over These are the pdvs for this set
                uint tUniquePdvIndex = 0;
                for (uint tGroupIndex = 0; tGroupIndex < mIpPdvTypes(tMeshSetIndex).size(); tGroupIndex++)
                {
                    for (uint tPdvIndex = 0; tPdvIndex < mIpPdvTypes(tMeshSetIndex)(tGroupIndex).size(); tPdvIndex++)
                    {
                        mUniqueIpPdvTypes(tMeshSetIndex)(tUniquePdvIndex++) = mIpPdvTypes(tMeshSetIndex)(tGroupIndex)(tPdvIndex);
                    }
                }

                // Create PDV hosts
                for (uint tNodeIndexOnSet = 0; tNodeIndexOnSet < aNodeIndicesPerSet(tMeshSetIndex).numel(); tNodeIndexOnSet++)
                {
                    // Create new host or add unique PDVs

                    moris_index tNodeIndex = aNodeIndicesPerSet(tMeshSetIndex)( tNodeIndexOnSet);
                    moris_id tNodeId       = aNodeIdsPerSet(tMeshSetIndex)( tNodeIndexOnSet);
                    moris_index tNodeOwner = aNodeOwnersPerSet(tMeshSetIndex)( tNodeIndexOnSet);

                    if (mIpPdvHosts(tNodeIndex) == nullptr)
                    {
                        mIpPdvHosts(tNodeIndex)
                                = std::make_shared<Interpolation_Pdv_Host>(
                                        this,
                                        tNodeIndex,
                                        tNodeId,
                                        tNodeOwner,
                                        aNodeCoordinates(aNodeIndicesPerSet(tMeshSetIndex)(tNodeIndexOnSet)),
                                        mUniqueIpPdvTypes(tMeshSetIndex) );
                    }
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::set_integration_pdv_types(Cell<Cell<Cell<PDV_Type>>> aPdvTypes)
        {
            // Check that number of sets is consistent
            uint tNumSets = aPdvTypes.size();

            // Set PDV types
            mIgPdvTypes = aPdvTypes;
            mUniqueIgPdvTypes.resize(tNumSets);

            // Unique PDV types
            for (uint tMeshSetIndex = 0; tMeshSetIndex < tNumSets; tMeshSetIndex++)
            {
                // Get number of unique PDVs
                uint tNumUniquePdvs = 0;
                for (uint tGroupIndex = 0; tGroupIndex < mIgPdvTypes(tMeshSetIndex).size(); tGroupIndex++)
                {
                    tNumUniquePdvs += mIgPdvTypes(tMeshSetIndex)(tGroupIndex).size();
                }
                mUniqueIgPdvTypes(tMeshSetIndex).resize(tNumUniquePdvs);

                // Copy PDV types over
                uint tUniquePdvIndex = 0;
                for (uint tGroupIndex = 0; tGroupIndex < mIgPdvTypes(tMeshSetIndex).size(); tGroupIndex++)
                {
                    for (uint tPdvIndex = 0; tPdvIndex < mIgPdvTypes(tMeshSetIndex)(tGroupIndex).size(); tPdvIndex++)
                    {
                        mUniqueIgPdvTypes(tMeshSetIndex)(tUniquePdvIndex++) = mIgPdvTypes(tMeshSetIndex)(tGroupIndex)(
                                tPdvIndex);
                    }
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::set_intersection_node(uint aNodeIndex, std::shared_ptr<Intersection_Node> aIntersectionNode)
        {
            // Check node index
            MORIS_ASSERT(mNumBackgroundNodesSet,
                         "Number of background nodes must be set in the PDV Host Manager before intersection nodes can be created.");
            MORIS_ASSERT(aNodeIndex == mIntersectionNodes.size(),
                         "Intersection nodes must be added to the PDV Host Manager in order by node index.");

            // Add intersection node
            mIntersectionNodes.push_back(aIntersectionNode);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::update_intersection_node(
                const moris_index & aNodeIndex,
                const moris_index & aNodeId,
                const moris_index & aNodeOwner)
        {
            //MORIS_ASSERT( mIntersectionNodes(aNodeIndex)!= nullptr,
            //        "Pdv_Host_Manager::update_intersection_node(), Intersection node doe not exist.");
            //FIXME the size of this cell should be correct. this is a hack to account for a wrong size
            if (aNodeIndex >= (sint)mIntersectionNodes.size())
            {
                mIntersectionNodes.resize(aNodeIndex + 1, nullptr);
            }
            if( mIntersectionNodes(aNodeIndex)!= nullptr )
            {
                mIntersectionNodes(aNodeIndex)->set_vertex_id(aNodeId);
                mIntersectionNodes(aNodeIndex)->set_vertex_owner(aNodeOwner);
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::set_requested_interpolation_pdv_types(Cell<PDV_Type> aPdvTypes)
        {
            mRequestedIpPdvTypes = aPdvTypes;
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::set_requested_integration_pdv_types(Cell<PDV_Type> aPdvTypes)
        {
            mRequestedIgPdvTypes = aPdvTypes;
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::create_interpolation_pdv(uint aNodeIndex, PDV_Type aPdvType, real aPdvVal)
        {
            mIpPdvHosts(aNodeIndex)->create_pdv(aPdvType, aPdvVal);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::create_interpolation_pdv(uint aNodeIndex, PDV_Type aPdvType, std::shared_ptr<Property> aProperty)
        {
            mIpPdvHosts(aNodeIndex)->create_pdv(aPdvType, aProperty);
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Pdv_Host_Manager::compute_diqi_dadv(const Matrix<DDSMat>& aFullADVIds)
        {
            // Check for ADV IDs
            MORIS_ERROR(mADVIdsSet, "PDV Host Manager must have ADV IDs set before computing sensitivities.");
            
            // Get dIQI/dPDV and dPDV/dADV
            sol::Dist_Vector* tdIQIdPDV = this->get_dQIdp();

            // Create factory for resulting distributed vector
            sol::Matrix_Vector_Factory tDistributedFactory;

            // Create maps
            std::shared_ptr<sol::Dist_Map> tOwnedADVMap = tDistributedFactory.create_map(mOwnedADVIds);
            std::shared_ptr<sol::Dist_Map> tFullADVMap = tDistributedFactory.create_map(aFullADVIds);

            // Create vectors
            sint tNumIQIs = tdIQIdPDV->get_num_vectors();
            sol::Dist_Vector* tdIQIdADV = tDistributedFactory.create_vector(tOwnedADVMap, tNumIQIs);
            sol::Dist_Vector* tFulldIQIdADV = tDistributedFactory.create_vector(tFullADVMap, tNumIQIs);

            // Initialize to zero
            tdIQIdADV->vec_put_scalar(0.0);

            // Loop of interpolation PDV hosts
            for (uint tPDVHostIndex = 0; tPDVHostIndex < mIpPdvHosts.size(); tPDVHostIndex++)
            {
                if (mIpPdvHosts(tPDVHostIndex) and mIpPdvHosts(tPDVHostIndex)->get_pdv_owning_processor() == par_rank())
                {
                    // Get number of PDVs
                    uint tNumPDVsOnHost = mIpPdvHosts(tPDVHostIndex)->get_num_pdvs();

                    // Assemble sensitivities
                    for (uint tPDVIndex = 0; tPDVIndex < tNumPDVsOnHost; tPDVIndex++)
                    {
                        // Get sensitivities
                        Matrix<DDRMat> tHostADVSensitivities = mIpPdvHosts(tPDVHostIndex)->get_sensitivities(tPDVIndex);

                        // Get PDV/ADV IDs
                        moris_id tPDVID = mIpPdvHosts(tPDVHostIndex)->get_pdv_id(tPDVIndex);
                        Matrix<DDSMat> tADVIds = mIpPdvHosts(tPDVHostIndex)->get_determining_adv_ids(tPDVIndex);

                        for (uint tVectorIndex = 0; tVectorIndex < (uint)tNumIQIs; tVectorIndex++)
                        {
                            Matrix<DDRMat> tIndividualSensitivity = (*tdIQIdPDV)(tPDVID, tVectorIndex) * tHostADVSensitivities;

                            // Fill matrix
                            tdIQIdADV->sum_into_global_values(tADVIds, tIndividualSensitivity, tVectorIndex);
                        }
                    }
                }
            }

            // Loop over intersection nodes for inserting
            for (uint tIntersectionIndex = 0; tIntersectionIndex < mIntersectionNodes.size(); tIntersectionIndex++)
            {
                if (mIntersectionNodes(tIntersectionIndex) and mIntersectionNodes(tIntersectionIndex)->get_vertex_owner() == par_rank())
                {
                    // Get starting ID and number of coordinates
                    uint tStartingGlobalIndex = mIntersectionNodes(tIntersectionIndex)->get_starting_pdv_id();
                    uint tNumCoordinates = mIntersectionNodes(tIntersectionIndex)->get_num_pdvs();

                    // Get first parent sensitivities and ADV IDs
                    Matrix<DDRMat> tHostADVSensitivities = mIntersectionNodes(tIntersectionIndex)->get_first_parent_sensitivities();
                    Matrix<DDSMat> tADVIds = mIntersectionNodes(tIntersectionIndex)->get_first_parent_determining_adv_ids();

                    // Assemble second parent
                    for (uint tCoordinateIndex = 0; tCoordinateIndex < tNumCoordinates; tCoordinateIndex++)
                    {
                        moris_id tPDVID = tStartingGlobalIndex + tCoordinateIndex;

                        for (uint tVectorIndex = 0; tVectorIndex < (uint)tNumIQIs; tVectorIndex++)
                        {
                            Matrix<DDRMat> tIndividualSensitivity = (*tdIQIdPDV)(tPDVID, tVectorIndex) * tHostADVSensitivities.get_row(tCoordinateIndex);

                            // Fill matrix
                            tdIQIdADV->sum_into_global_values(tADVIds, tIndividualSensitivity, tVectorIndex);
                        }
                    }

                    // Get second parent sensitivities and ADV IDs
                    tHostADVSensitivities = mIntersectionNodes(tIntersectionIndex)->get_second_parent_sensitivities();
                    tADVIds = mIntersectionNodes(tIntersectionIndex)->get_second_parent_determining_adv_ids();

                    // Assemble first parent
                    for (uint tCoordinateIndex = 0; tCoordinateIndex < tNumCoordinates; tCoordinateIndex++)
                    {
                        moris_id tPDVID = tStartingGlobalIndex + tCoordinateIndex;

                        for (uint tVectorIndex = 0; tVectorIndex < (uint)tNumIQIs; tVectorIndex++)
                        {
                            Matrix<DDRMat> tIndividualSensitivity = (*tdIQIdPDV)(tPDVID, tVectorIndex) * tHostADVSensitivities.get_row(tCoordinateIndex);

                            // Fill matrix
                            tdIQIdADV->sum_into_global_values(tADVIds, tIndividualSensitivity, tVectorIndex);
                        }
                    }

                }
            }

            // Global assembly
            tdIQIdADV->vector_global_assembly();

            // Import
            tFulldIQIdADV->import_local_to_global(*tdIQIdADV);
            
            // Extract values
            Matrix<DDRMat> tFullSensitivity(0, 0);
            tFulldIQIdADV->extract_copy(tFullSensitivity);
            tFullSensitivity = trans(tFullSensitivity);
            
            return tFullSensitivity;
        }

        //--------------------------------------------------------------------------------------------------------------

        sol::Dist_Matrix* Pdv_Host_Manager::compute_dpdv_dadv()
        {
            // Create factory for distributed matrix
            sol::Matrix_Vector_Factory tDistributedFactory;

            // Create row and column maps
            std::shared_ptr<sol::Dist_Map> tOwnedADVMap = tDistributedFactory.create_map(mOwnedADVIds);
            std::shared_ptr<sol::Dist_Map> tOwnedPDVMap = tDistributedFactory.create_map(mOwnedPdvLocalToGlobalMap);

            // Create vector
            sol::Dist_Matrix* tdPDVdADV = tDistributedFactory.create_matrix(tOwnedPDVMap, tOwnedADVMap);

            // Global assembly
            tdPDVdADV->matrix_global_assembly();

            return tdPDVdADV;
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::communicate_dof_types( moris::Cell< enum PDV_Type > & aPdvTypeList )
        {
            // Get processor size
            int tSize = par_size();

            // Get number of local dv types
            moris::sint tNumLocalDvTypes = aPdvTypeList.size();

            // Get number of global dv types
            moris::sint tNumMaxGlobalDvTypes = sum_all( tNumLocalDvTypes );

            if ( par_rank() == 0 )
            {
                // Set size of of pdv type list = number of global types
                mPdvTypeList.resize( tNumMaxGlobalDvTypes );
            }

            // Create list containing the number of local dof types
            moris::Cell < moris::sint > tNumLocalDvTypesList ( tSize );

            // Insert number of local dof types into list containing the number of local dof types
            MPI_Allgather( &tNumLocalDvTypes, 1, MPI_UNSIGNED, (tNumLocalDvTypesList.data()).data(), 1, MPI_UNSIGNED,  MPI_COMM_WORLD );

            // Create list containing the offsets of the local dof types in relation to processor 0
            moris::Cell< moris::sint > tDvTypeOffset( tSize, 0 );

            // Fill the list with the corresponding offsets
            for ( int Ip = 1; Ip < tSize; ++Ip )
            {
                tDvTypeOffset( Ip ) = tDvTypeOffset( Ip-1 ) + tNumLocalDvTypesList( Ip-1 );
            }

            // Assemble list containing all used dof types. Dof types are not unique
            MPI_Gatherv(
                    ((aPdvTypeList.data()).data()),
                    tNumLocalDvTypes,
                    MPI_UNSIGNED,
                    (mPdvTypeList.data()).data(),
                    (tNumLocalDvTypesList.data()).data(),
                    (tDvTypeOffset.data()).data(),
                    MPI_UNSIGNED,
                    0,
                    MPI_COMM_WORLD );

            // Temporary variable for mPdvTypeList size
            moris::uint tPdvTypeListSize;

            if ( par_rank() == 0 )
            {
                // Sort this created list
                std::sort( ( mPdvTypeList.data() ).data(), ( mPdvTypeList.data() ).data() + mPdvTypeList.size() );

                // use std::unique and std::distance to create list containing all used dof types. This list is unique
                auto last = std::unique( ( mPdvTypeList.data() ).data(), ( mPdvTypeList.data() ).data() + mPdvTypeList.size() );
                auto pos  = std::distance( ( mPdvTypeList.data() ).data(), last );

                mPdvTypeList.resize( pos );

                tPdvTypeListSize = mPdvTypeList.size();
            }

            // Bcast size of mPdvTypeList on processor 0
            MPI_Bcast( & tPdvTypeListSize, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD );

            // Resize mPdvTypeList on all processors
            mPdvTypeList.resize( tPdvTypeListSize );

            // Bcast unique mPdvTypeList to all processors
            MPI_Bcast( (mPdvTypeList.data()).data(), mPdvTypeList.size(), MPI_UNSIGNED, 0, MPI_COMM_WORLD );
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::create_dv_type_map()
        {
            //Get number of unique adofs of this equation object
            moris::uint tNumUniquePdvTypes = mPdvTypeList.size();

            // Get maximal dv type enum number
            moris::sint tMaxDvTypeEnumNumber = 0;

            // Loop over all pdv types to get the highest enum index
            for ( moris::uint Ii = 0; Ii < tNumUniquePdvTypes; Ii++ )
            {
                tMaxDvTypeEnumNumber = std::max( tMaxDvTypeEnumNumber, static_cast< int >( mPdvTypeList( Ii ) ) );
            }

            // +1 because c++ is 0 based
            tMaxDvTypeEnumNumber = tMaxDvTypeEnumNumber + 1;

            // Set size of mapping matrix
            mPdvTypeMap.set_size( tMaxDvTypeEnumNumber, 1, -1 );

            // Loop over all pdv types to create the mapping matrix
            for ( moris::uint Ii = 0; Ii < tNumUniquePdvTypes; Ii++ )
            {
                mPdvTypeMap( static_cast< int >( mPdvTypeList( Ii ) ), 0 ) = Ii;
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::communicate_check_if_owned_pdv_exists()
        {
            // Build communication table map to determine the right position for each processor rank. +1 because c++ is 0 based
            Matrix< DDSMat > tCommTableMap ( mCommTable.max() + 1, 1, -1);

            moris::uint tNumCommProcs = mCommTable.numel();

            // Loop over communication table to fill the communication table map
            for ( moris::uint Ik = 0; Ik < tNumCommProcs; Ik++ )
            {
                tCommTableMap( mCommTable( Ik ), 0 ) = Ik;
            }

            // Loop over all different pdv types for IP node pdvs
            for ( moris::uint Ij = 0; Ij < mPdvTypeList.size(); Ij++ )
            {
                enum PDV_Type tPdvType = mPdvTypeList( Ij );

                moris::Cell< Matrix< DDUMat > > tSharedPdvPosGlobal( tNumCommProcs );

                // Set Mat to store number of shared pdvs per processor
                Matrix< DDUMat > tNumSharedPdvsPerProc( tNumCommProcs, 1, 0 );

                // Loop over pdvs per type. Count number of pdvs per proc which have to be communicated
                for ( moris::uint Ib = 0; Ib < mIpPdvHosts.size(); Ib++ )
                {
                    // Check if pdv at this position is not NULL
                    if ( mIpPdvHosts( Ib )->get_pdv_exists( tPdvType ) )
                    {
                        // Check if owning processor is this processor
                        if ( mIpPdvHosts( Ib )->get_pdv_owning_processor() != par_rank() )
                        {
                            // get owning processor
                            moris::moris_id tProcID = mIpPdvHosts( Ib )->get_pdv_owning_processor();

                            moris::sint tProcIdPos = tCommTableMap( tProcID, 0 );

                            MORIS_ASSERT( tProcIdPos != -1, "Dof_Manager::communicate_check_if_owned_pdv_exists: Map returns proc rank -1. Check communication table");

                            // Add +1 to the processor number of shared dofs per processor
                            tNumSharedPdvsPerProc( tProcIdPos, 0) = tNumSharedPdvsPerProc( tProcIdPos, 0) + 1;
                        }
                    }
                }

                // Set size of the moris::Mats in the Cell
                for ( moris::uint Ik = 0; Ik < tNumCommProcs; Ik++ )
                {
                    if ( tNumSharedPdvsPerProc( Ik, 0 ) != 0 )
                    {
                        tSharedPdvPosGlobal( Ik ).set_size( tNumSharedPdvsPerProc( Ik, 0 ), 1);
                    }
                }

                // Temporary Mat to add external adof ids at the next spot in the matrix which will be communicated
                Matrix< DDUMat > tShredPdvPosPerProc( tNumCommProcs, 1, 0 );

                // Loop over pdv per type
                for ( moris::uint Ia = 0; Ia < mIpPdvHosts.size(); Ia++ )
                {
                    // Check if pdv at this position is not NULL
                    if ( mIpPdvHosts( Ia )->get_pdv_exists( tPdvType ) )
                    {
                        // Check if owning processor is this processor
                        if ( mIpPdvHosts( Ia )->get_pdv_owning_processor() != par_rank() )
                        {
                            // Get owning processor
                            moris::uint tProcID = mIpPdvHosts( Ia )->get_pdv_owning_processor();

                            moris::sint tProcIdPos = tCommTableMap( tProcID, 0 );

                            // Add owning processor id to moris::Mat
                            tSharedPdvPosGlobal( tProcIdPos )( tShredPdvPosPerProc( tProcIdPos, 0 ), 0 ) =
                                    mIpPdvHosts( Ia )->get_pdv_vertex_id();

                            tShredPdvPosPerProc( tProcIdPos, 0 ) = tShredPdvPosPerProc( tProcIdPos, 0 ) + 1;
                        }
                    }
                }

                // receiving list
                moris::Cell< Matrix< DDUMat > > tMatsToReceive;

                barrier();

                // Communicate position of shared pdv to the owning processor
                communicate_mats(
                        mCommTable,
                        tSharedPdvPosGlobal,
                        tMatsToReceive );

                // Loop over all Mats set dummy owned pdv
                for ( moris::uint Ik = 0; Ik < tMatsToReceive.size(); Ik++ )
                {
                    for ( moris::uint Ii = 0; Ii < tMatsToReceive( Ik ).numel(); Ii++ )
                    {
                        // Get owned pdv Id
                        auto tIter = mIPVertexIdtoIndMap.find( tMatsToReceive( Ik )( Ii ) );

                        moris::uint tLocalPdvInd = tIter->second;

                        if ( mIpPdvHosts( tLocalPdvInd )->get_pdv_exists( tPdvType )  )
                        {
							MORIS_ERROR(false," add the follwing lines");
                            // FIXME create pdv
                        }
                    }
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::get_num_pdvs()
        {
            // Loop over all different pdv types for IP node pdvs
            for ( moris::uint Ij = 0; Ij < mPdvTypeList.size(); Ij++ )
            {
                enum PDV_Type tPdvType = mPdvTypeList( Ij );

                for ( moris::uint Ib = 0; Ib < mIpPdvHosts.size(); Ib++ )
                {
                    // Check if pdv at this position is not NULL
                    if ( mIpPdvHosts( Ib )->get_pdv_exists( tPdvType ) )
                    {
                        // Check if owning processor is this processor
                        if ( mIpPdvHosts( Ib )->get_pdv_owning_processor() == par_rank() )
                        {
                            mNumOwnedPdvs++;
                        }
                        mNumOwnedAndSharedPdvs++;
                    }
                }
            }

            // Loop over intersection node pdvs
            for ( moris::uint Ij = 0; Ij < mIntersectionNodes.size(); Ij++ )
            {
                if( mIntersectionNodes( Ij ) != nullptr )
                {
                    uint tNumPdvsOnIntersectionNode = mIntersectionNodes( Ij )->get_num_pdvs();

                    if( mIntersectionNodes( Ij )->get_vertex_owner() == par_rank() )
                    {
                        mNumOwnedPdvs = mNumOwnedPdvs + tNumPdvsOnIntersectionNode;
                    }
                    mNumOwnedAndSharedPdvs = mNumOwnedAndSharedPdvs + tNumPdvsOnIntersectionNode;
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        uint Pdv_Host_Manager::communicate_pdv_offsets( const moris::uint & aNumOwnedPdvs )
        {
            MORIS_LOG_INFO( "System has a total of %-5i pdvs.", sum_all(aNumOwnedPdvs) );

            // Get list containing the number of owned pdvs of each processor
            Matrix< DDUMat > tNumOwnedPdvsList;
            comm_gather_and_broadcast( aNumOwnedPdvs, tNumOwnedPdvsList );

            Matrix< DDUMat > tOwnedPdvsOffsetList( tNumOwnedPdvsList.numel(), 1, 0 );

            // Loop over all entries to create the offsets. Starting with 1
            for ( moris::uint Ij = 1; Ij < tOwnedPdvsOffsetList.numel(); Ij++ )
            {
                // Add the number of owned pdvs of the previous processor to the offset of the previous processor
                tOwnedPdvsOffsetList( Ij, 0 ) = tOwnedPdvsOffsetList( Ij-1, 0 ) + tNumOwnedPdvsList( Ij-1, 0 );
            }

            return tOwnedPdvsOffsetList( par_rank(), 0);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::set_owned_pdv_ids( uint aPdvOffset )
        {
            moris::uint tOwnedIdCounter = aPdvOffset;

            tOwnedIdCounter = this->set_owned_interpolation_pdv_ids( tOwnedIdCounter );
            this->set_owned_intersection_node_pdv_ids( tOwnedIdCounter );
        }

        //--------------------------------------------------------------------------------------------------------------

        moris_id Pdv_Host_Manager::set_owned_interpolation_pdv_ids( moris_id aOwnedIdCounter )
        {
            // Loop over all different pdv types for IP node pdvs
            for ( moris::uint Ij = 0; Ij < mPdvTypeList.size(); Ij++ )
            {
                enum PDV_Type tPdvType = mPdvTypeList( Ij );

                // Loop over pdvs per type.
                for ( moris::uint Ib = 0; Ib < mIpPdvHosts.size(); Ib++ )
                {
                    // Check if pdv at this position is not NULL
                    if ( mIpPdvHosts( Ib )->get_pdv_exists( tPdvType ) )
                    {
                        // Check if owning processor is this processor
                        if ( mIpPdvHosts( Ib )->get_pdv_owning_processor() == par_rank() )
                        {
                            mIpPdvHosts( Ib )->set_pdv_id( tPdvType, aOwnedIdCounter );
                            aOwnedIdCounter++;
                        }
                    }
                }
            }
            return aOwnedIdCounter;
        }

        //--------------------------------------------------------------------------------------------------------------

        moris_id Pdv_Host_Manager::set_owned_intersection_node_pdv_ids( moris_id aOwnedIdCounter )
        {
            // Loop over intersection node pdvs
            for ( moris::uint Ij = 0; Ij < mIntersectionNodes.size(); Ij++ )
            {
                if( mIntersectionNodes( Ij ) != nullptr )
                {
                    if( mIntersectionNodes( Ij )->get_vertex_owner() == par_rank() )
                    {
                        mIntersectionNodes( Ij )->set_starting_pdv_id( aOwnedIdCounter );

                        uint tNumPdvsOnIntersectionNode = mIntersectionNodes( Ij )->get_num_pdvs();

                        aOwnedIdCounter = aOwnedIdCounter + tNumPdvsOnIntersectionNode;
                    }
                }
            }

            return aOwnedIdCounter;
        }
        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::communicate_shared_pdv_ids()
        {
            this->communicate_shared_interpolation_pdv_ids();
            this->communicate_shared_intersection_node_pdv_ids();
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::communicate_shared_interpolation_pdv_ids()
        {
            // Build communication table map to determine the right position for each processor rank. +1 because c++ is 0 based
            Matrix< DDSMat > tCommTableMap ( mCommTable.max() + 1, 1, -1);

            moris::uint tNumCommProcs = mCommTable.numel();

            // Loop over communication table to fill the communication table map
            for ( moris::uint Ik = 0; Ik < tNumCommProcs; Ik++ )
            {
                tCommTableMap( mCommTable( Ik ), 0 ) = Ik;
            }

            // Loop over all different pdv types for IP node pdvs
            for ( moris::uint Ij = 0; Ij < mPdvTypeList.size(); Ij++ )
            {
                enum PDV_Type tPdvType = mPdvTypeList( Ij );

                // Set Mat to store number of shared pdv per processor
                Matrix< DDUMat > tNumSharedPdvsPerProc( tNumCommProcs, 1, 0 );

                moris::uint tCounter = 0;
                moris::uint tSharedCounter = 0;

                moris::Cell< Matrix< DDUMat > > tSharedPdvPosGlobal( tNumCommProcs );
                moris::Cell< Matrix< DDUMat > > tSharedPdvPosLocal( tNumCommProcs );

                // Loop over pdvs per type. Count number of pdvs per proc which have to be communicated
                for ( moris::uint Ib = 0; Ib < mIpPdvHosts.size(); Ib++ )
                {
                    // Check if pdv at this position is not NULL
                    if ( mIpPdvHosts( Ib )->get_pdv_exists( tPdvType ) )
                    {
                        // Check if owning processor is this processor
                        if ( mIpPdvHosts( Ib )->get_pdv_owning_processor() != par_rank() )
                        {
                            moris::moris_id tProcID  = mIpPdvHosts( Ib )->get_pdv_owning_processor();

                            moris::sint tProcIdPos = tCommTableMap( tProcID );

                            // Add +1 to the processor number of shared dv per processor
                            tNumSharedPdvsPerProc( tProcIdPos )++;

                            tSharedCounter++;
                        }
                    }
                }

                // Set size of the moris::Mats in the Cell
                for ( moris::uint Ik = 0; Ik < tNumCommProcs; Ik++ )
                {
                    if ( tNumSharedPdvsPerProc( Ik, 0 ) != 0 )
                    {
                        tSharedPdvPosGlobal( Ik ).set_size( tNumSharedPdvsPerProc( Ik, 0 ), 1);
                        tSharedPdvPosLocal( Ik ).set_size( tNumSharedPdvsPerProc( Ik, 0 ), 1);
                    }
                }

                // Temporary Mat to add external pdv ids at the next spot in the matrix which will be communicated
                Matrix< DDUMat > tSharedPdvPosPerProc( tNumCommProcs, 1, 0 );

                // Loop over pdvs
                for ( moris::uint Ib = 0; Ib < mIpPdvHosts.size(); Ib++ )
                {
                    // Check if pdv at this position is not NULL
                    if ( mIpPdvHosts( Ib )->get_pdv_exists( tPdvType ) )
                    {
                        // Check if owning processor is this processor
                        if ( mIpPdvHosts( Ib )->get_pdv_owning_processor() != par_rank() )
                        {
                            moris::moris_id tProcID  = mIpPdvHosts( Ib )->get_pdv_owning_processor();

                            moris::sint tProcIdPos = tCommTableMap( tProcID );

                            // Add owning processor id to moris::Mat
                            tSharedPdvPosGlobal( tProcIdPos )( tSharedPdvPosPerProc( tProcIdPos ) ) =
                                    mIpPdvHosts( Ib )->get_pdv_vertex_id();

                            // Add pdv position to Mat
                            tSharedPdvPosLocal( tProcIdPos ) ( tSharedPdvPosPerProc( tProcIdPos ) ) = tCounter;

                            tSharedPdvPosPerProc( tProcIdPos )++;
                        }
                        tCounter++;
                    }
                }

                // receiving list
                moris::Cell< Matrix< DDUMat > > tMatsToReceive;

                barrier();

                // Communicate position of shared pdvs to the owning processor
                communicate_mats(
                        mCommTable,
                        tSharedPdvPosGlobal,
                        tMatsToReceive );

                // Create List of Mats containing the shared node Ids
                moris::Cell< Matrix< DDUMat > > tSharedPdvIdList( tNumCommProcs );

                // Loop over all Mats setting the size
                for ( moris::uint Ik = 0; Ik < tMatsToReceive.size(); Ik++ )
                {
                    tSharedPdvIdList( Ik ).set_size( tMatsToReceive( Ik ).numel(), 1);
                }

                // Loop over all received positions and get the pdv id of the owning pdv
                for ( moris::uint Ik = 0; Ik < tMatsToReceive.size(); Ik++ )
                {
                    for ( moris::uint Ii = 0; Ii < tMatsToReceive( Ik ).numel(); Ii++ )
                    {
                        // Get owned pdv Id
                        auto tIter = mIPVertexIdtoIndMap.find( tMatsToReceive( Ik )( Ii ) );

                        moris::uint tLocalPdvInd = tIter->second;

                        MORIS_ASSERT( ( mIpPdvHosts( tLocalPdvInd )->get_pdv_owning_processor() ) == par_rank(),
                                "Pdv_Host_Manager::communicate_shared_interpolation_pdv_ids(): Pdv not owned by this processor");

                        tSharedPdvIdList( Ik )( Ii ) = mIpPdvHosts( tLocalPdvInd )->get_pdv_id( tPdvType );
                    }
                }

                moris::Cell< Matrix< DDUMat > > tMatsToReceive2;

                barrier();

                // Communicate owned pdv Id back to the processor with the shared pdv
                communicate_mats(
                        mCommTable,
                        tSharedPdvIdList,
                        tMatsToReceive2 );

                moris::uint tPdvPosCounter = 0;

                Matrix< DDUMat > tListSharedPdvIds( tSharedCounter, 1, MORIS_UINT_MAX );
                Matrix< DDUMat > tListSharedPdvPos( tSharedCounter, 1, MORIS_UINT_MAX );

                // assemble Ids in list of shared pdv ids and assemble the corresponding postions
                for ( moris::uint Ik = 0; Ik < tMatsToReceive2.size(); Ik++ )
                {
                    if( tMatsToReceive2( Ik ).numel() >= 1)
                    {
                        tListSharedPdvIds( {tPdvPosCounter, tPdvPosCounter + tMatsToReceive2( Ik ).numel() -1 }, { 0, 0 } ) =
                                tMatsToReceive2( Ik ).matrix_data();

                        tListSharedPdvPos( {tPdvPosCounter, tPdvPosCounter + tSharedPdvPosLocal( Ik ).numel() -1 }, { 0, 0 } ) =
                                tSharedPdvPosLocal( Ik ).matrix_data();

                        tPdvPosCounter =tPdvPosCounter + tMatsToReceive2( Ik ).numel();
                    }
                }

                if ( tListSharedPdvIds.numel() != 0 )
                {
                    MORIS_ASSERT( tListSharedPdvIds.max() != MORIS_UINT_MAX,
                            "Pdv_Host_Manager::communicate_shared_interpolation_pdv_ids(), communicated Ids not set correctly");

                    MORIS_ASSERT( tListSharedPdvPos.max() != MORIS_UINT_MAX,
                            "Pdv_Host_Manager::communicate_shared_interpolation_pdv_ids(), positions for communicated Ids not set correctly");
                }

                // Set the Id of the shared pdvs
                for ( moris::uint Ij = 0; Ij < tListSharedPdvIds.numel(); Ij++ )
                {
                    mIpPdvHosts( tListSharedPdvPos( Ij ) )->set_pdv_id( tPdvType, tListSharedPdvIds( Ij ) );
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::communicate_shared_intersection_node_pdv_ids()
        {
            // Build communication table map to determine the right position for each processor rank. +1 because c++ is 0 based
            Matrix< DDSMat > tCommTableMap ( mCommTable.max() + 1, 1, -1);

            moris::uint tNumCommProcs = mCommTable.numel();

            // Loop over communication table to fill the communication table map
            for ( moris::uint Ik = 0; Ik < tNumCommProcs; Ik++ )
            {
                tCommTableMap( mCommTable( Ik ), 0 ) = Ik;
            }

            moris::uint tCounter = 0;
            moris::uint tSharedCounter = 0;

            moris::Cell< Matrix< DDUMat > > tSharedPdvPosGlobal( tNumCommProcs );
            moris::Cell< Matrix< DDUMat > > tSharedPdvPosLocal( tNumCommProcs );

            // Set Mat to store number of shared pdv per processor
            Matrix< DDUMat > tNumSharedPdvsPerProc( tNumCommProcs, 1, 0 );

            // Loop over pdvs
            for ( moris::uint Ij = 0; Ij < mIntersectionNodes.size(); Ij++ )
            {
                // Check if pdv at this position is not NULL
                if ( mIntersectionNodes( Ij ) != nullptr )
                {
                    // Check if owning processor is this processor
                    if ( mIntersectionNodes( Ij )->get_vertex_owner() != par_rank() )
                    {
                        // get owning procssor
                        moris::moris_id tProcID = mIntersectionNodes( Ij )->get_vertex_owner();

                        moris::sint tProcIdPos = tCommTableMap( tProcID );

                        // Add +1 to the processor number of shared dv per processor
                        tNumSharedPdvsPerProc( tProcIdPos )++;

                        tSharedCounter++;
                    }
                }
            }

            // Set size of the moris::Mats in the Cell
            for ( moris::uint Ik = 0; Ik < tNumCommProcs; Ik++ )
            {
                if ( tNumSharedPdvsPerProc( Ik, 0 ) != 0 )
                {
                    tSharedPdvPosGlobal( Ik ).set_size( tNumSharedPdvsPerProc( Ik, 0 ), 1);
                    tSharedPdvPosLocal( Ik ).set_size( tNumSharedPdvsPerProc( Ik, 0 ), 1);
                }
            }

            // Temporary Mat to add external pdv ids at the next spot in the matrix which will be communicated
            Matrix< DDUMat > tSharedPdvPosPerProc( tNumCommProcs, 1, 0 );

            // Loop over pdvs
            for ( moris::uint Ij = 0; Ij < mIntersectionNodes.size(); Ij++ )
            {
                // Check if pdv at this position is not NULL
                if ( mIntersectionNodes( Ij ) != nullptr )
                {
                    // Check if owning processor is this processor
                    if ( mIntersectionNodes( Ij )->get_vertex_owner() != par_rank() )
                    {
                        // Get owning processor
                        moris::moris_id tProcID = mIntersectionNodes( Ij )->get_vertex_owner();

                        moris::sint tProcIdPos = tCommTableMap( tProcID );

                        // Add owning processor id to moris::Mat
                        tSharedPdvPosGlobal( tProcIdPos )( tSharedPdvPosPerProc( tProcIdPos ) ) =
                                mIntersectionNodes( Ij )->get_vertex_id();

                        // Add pdv position to Mat
                        tSharedPdvPosLocal( tProcIdPos ) ( tSharedPdvPosPerProc( tProcIdPos ) ) = tCounter;

                        tSharedPdvPosPerProc( tProcIdPos )++;
                    }
                }
                tCounter++;
            }

            // receiving list
            moris::Cell< Matrix< DDUMat > > tMatsToReceive;

            barrier();

            // Communicate position of shared pdvs to the owning processor
            communicate_mats(
                    mCommTable,
                    tSharedPdvPosGlobal,
                    tMatsToReceive );

            // Create List of Mats containing the shared node Ids
            moris::Cell< Matrix< DDUMat > > tSharedPdvIdList( tNumCommProcs );

            // Loop over all Mats setting the size
            for ( moris::uint Ik = 0; Ik < tMatsToReceive.size(); Ik++ )
            {
                tSharedPdvIdList( Ik ).set_size( tMatsToReceive( Ik ).numel(), 1);
            }

            // Loop over all received positions and get the pdv id of the owning pdv
            for ( moris::uint Ik = 0; Ik < tMatsToReceive.size(); Ik++ )
            {
                for ( moris::uint Ii = 0; Ii < tMatsToReceive( Ik ).numel(); Ii++ )
                {
                    // Get owned pdv Id
                    auto tIter = mIGVertexIdtoIndMap.find( tMatsToReceive( Ik )( Ii ) );

                    moris::uint tLocalPdvInd = tIter->second;

                    MORIS_ASSERT( ( mIntersectionNodes( tLocalPdvInd )->get_vertex_owner() ) == par_rank(),
                            "Pdv_Host_Manager::communicate_shared_pdv_ids(): Pdv not owned by this processor");

                    tSharedPdvIdList( Ik )( Ii ) = mIntersectionNodes( tLocalPdvInd )->get_starting_pdv_id();
                }
            }

            moris::Cell< Matrix< DDUMat > > tMatsToReceive2;

            barrier();

            // Communicate owned pdv Id back to the processor with the shared pdv
            communicate_mats(
                    mCommTable,
                    tSharedPdvIdList,
                    tMatsToReceive2 );

            moris::uint tPdvPosCounter = 0;

            Matrix< DDUMat > tListSharedPdvIds( tSharedCounter, 1, MORIS_UINT_MAX );
            Matrix< DDUMat > tListSharedPdvPos( tSharedCounter, 1, MORIS_UINT_MAX );

            // assemble Ids in list of shared pdv ids and assemble the corresponding postions
            for ( moris::uint Ik = 0; Ik < tMatsToReceive2.size(); Ik++ )
            {
                if( tMatsToReceive2( Ik ).numel() >= 1)
                {
                    tListSharedPdvIds( {tPdvPosCounter, tPdvPosCounter + tMatsToReceive2( Ik ).numel() -1 }, { 0, 0 } ) =
                            tMatsToReceive2( Ik ).matrix_data();

                    tListSharedPdvPos( {tPdvPosCounter, tPdvPosCounter + tSharedPdvPosLocal( Ik ).numel() -1 }, { 0, 0 } ) =
                            tSharedPdvPosLocal( Ik ).matrix_data();

                    tPdvPosCounter =tPdvPosCounter + tMatsToReceive2( Ik ).numel();
                }
            }

            if( tListSharedPdvIds.numel() != 0 )
            {
                MORIS_ASSERT( tListSharedPdvIds.max() != MORIS_UINT_MAX,
                        "Pdv_Host_Manager::communicate_shared_pdv_ids(), communicated Ids not set correctly");

                MORIS_ASSERT( tListSharedPdvPos.max() != MORIS_UINT_MAX,
                        "Pdv_Host_Manager::communicate_shared_pdv_ids(), positions for communicated Ids not set correctly");
            }

            // Set the Id of the shared pdvs
            for ( moris::uint Ij = 0; Ij < tListSharedPdvIds.numel(); Ij++ )
            {
                mIntersectionNodes( tListSharedPdvPos( Ij ) )->set_starting_pdv_id( tListSharedPdvIds( Ij ) );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::build_local_to_global_maps()
        {
            mOwnedPdvLocalToGlobalMap.set_size( mNumOwnedPdvs, 1, -1 );
            mOwnedAndSharedPdvLocalToGlobalMap.set_size( mNumOwnedAndSharedPdvs, 1, -1 );

            uint tCounter = 0;
            uint tCounter2 = 0;

            // Loop over all different pdv types for IP node pdvs
             for ( moris::uint Ij = 0; Ij < mPdvTypeList.size(); Ij++ )
             {
                 enum PDV_Type tPdvType = mPdvTypeList( Ij );

                 // Loop over pdvs per type. Count number of pdvs per proc which have to be communicated
                 for ( moris::uint Ib = 0; Ib < mIpPdvHosts.size(); Ib++ )
                 {
                     // Check if pdv at this position is not NULL
                     if ( mIpPdvHosts( Ib )->get_pdv_exists( tPdvType ) )
                     {
                         // Check if owning processor is this processor
                         if ( mIpPdvHosts( Ib )->get_pdv_owning_processor() == par_rank() )
                         {
                             mOwnedPdvLocalToGlobalMap( tCounter ++)  = mIpPdvHosts( Ib )->get_pdv_id( tPdvType );
                         }
                         mOwnedAndSharedPdvLocalToGlobalMap( tCounter2 ++) = mIpPdvHosts( Ib )->get_pdv_id( tPdvType );
                     }
                 }
             }

            // Loop over intersection node pdvs
            for ( moris::uint Ij = 0; Ij < mIntersectionNodes.size(); Ij++ )
            {
                if( mIntersectionNodes( Ij ) != nullptr )
                {
                    uint tNumPdvsOnIntersectionNode = mIntersectionNodes( Ij )->get_num_pdvs();

                    for( moris::uint Ik = 0; Ik < tNumPdvsOnIntersectionNode; Ik++)
                    {
                        if( mIntersectionNodes( Ij )->get_vertex_owner() == par_rank() )
                        {
                            mOwnedPdvLocalToGlobalMap( tCounter ++) = mIntersectionNodes( Ij )->get_starting_pdv_id() + Ik;
                        }
                        mOwnedAndSharedPdvLocalToGlobalMap( tCounter2 ++) = mIntersectionNodes( Ij )->get_starting_pdv_id() + Ik;
                    }
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::create_pdv_ids()
        {
            uint tPdvOffset = 0;

            if ( !(par_size() <= 1) )
            {
                //this->communicate_check_if_owned_pdv_exists();
            }

            this->get_num_pdvs();

            if ( !(par_size() <= 1) )
            {
                tPdvOffset = this->communicate_pdv_offsets( mNumOwnedPdvs );
            }

            this->set_owned_pdv_ids( tPdvOffset );

            if ( !(par_size() <= 1) )
            {
                this->communicate_shared_pdv_ids();
            }

            this->build_local_to_global_maps();
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
