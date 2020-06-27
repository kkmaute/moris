#include "cl_GEN_Pdv_Host_Manager.hpp"
#include "fn_sum.hpp"

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

        void Pdv_Host_Manager::get_ip_dv_types_for_set(const moris::moris_index aIPMeshSetIndex, Cell<Cell<PDV_Type>>& aPdvTypes)
        {
            if (mIpPdvTypes.size() > 0)
            {
                aPdvTypes = mIpPdvTypes(aIPMeshSetIndex);
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::get_ig_dv_types_for_set(const moris::moris_index aIGMeshSetIndex, Cell<Cell<PDV_Type>>& aPdvTypes)
        {
            if (mIgPdvTypes.size() > 0)
            {
                aPdvTypes = mIgPdvTypes(aIGMeshSetIndex);
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::get_ip_unique_dv_types_for_set(const moris_index aIPMeshSetIndex, Cell<PDV_Type>& aPdvTypes)
        {
            if (mUniqueIpPdvTypes.size() > 0)
            {
                aPdvTypes =  mUniqueIpPdvTypes(aIPMeshSetIndex);
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::get_ig_unique_dv_types_for_set(const moris::moris_index aIGMeshSetIndex, Cell<PDV_Type>& aPdvTypes)
        {
            if (mUniqueIgPdvTypes.size() > 0)
            {
                aPdvTypes =  mUniqueIgPdvTypes(aIGMeshSetIndex);
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::get_ip_pdv_value(const Matrix<IndexMat>&   aNodeIndices,
                                                const Cell<PDV_Type>&       aPdvTypes,
                                                Cell<Matrix<DDRMat>>&     aDvValues)
        {
            // get the number of node indices requested
            uint tNumIndices = aNodeIndices.length();

            // get the number of dv types requested
            uint tNumTypes = aPdvTypes.size();

            // set size for list of dv values
            aDvValues.resize(tNumTypes);

            // loop over the node indices
            for (uint iInd = 0; iInd < tNumIndices; iInd++)
            {
                // loop over the requested dv types
                for (uint iType = 0; iType < tNumTypes; iType++)
                {
                    aDvValues(iType).resize(tNumIndices, 1);
                    aDvValues(iType)(iInd) = mIpPdvHosts(aNodeIndices(iInd))->get_pdv_value(aPdvTypes(iType));
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::get_ip_pdv_value(const Matrix<IndexMat>&   aNodeIndices,
                                                const Cell<PDV_Type>&       aPdvTypes,
                                                Cell<Matrix<DDRMat>>&     aDvValues,
                                                Cell<Matrix<DDSMat>>&     aIsActiveDv)
        {
            // get the number of node indices requested
            uint tNumIndices = aNodeIndices.length();

            // get the number of dv types requested
            uint tNumTypes = aPdvTypes.size();

            // set size for list of dv values
            aDvValues.resize( tNumTypes );

            // set size for list of active flags
            aIsActiveDv.resize( tNumTypes );

            // loop over the requested dv types
            for ( uint iType = 0; iType < tNumTypes; iType++ )
            {
                aDvValues( iType ).resize( tNumIndices, 1 );
                aIsActiveDv( iType ).resize( tNumIndices, 1 );

                // loop over the node indices
                for ( uint iInd = 0; iInd < tNumIndices; iInd++ )
                {
                    aDvValues( iType )( iInd )   = mIpPdvHosts( aNodeIndices( iInd ) )->get_pdv_value( aPdvTypes( iType ) );
                    aIsActiveDv( iType )( iInd ) = mIpPdvHosts( aNodeIndices( iInd ) )->is_active_type( aPdvTypes( iType ) );
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::get_ig_pdv_value(const Matrix<IndexMat>&   aNodeIndices,
                                                const Cell<PDV_Type>&       aPdvTypes,
                                                Cell<Matrix<DDRMat>>&     aDvValues)
        {
            // get the number of node indices requested
            uint tNumIndices = aNodeIndices.length();

            // get the number of dv types requested
            uint tNumTypes = aPdvTypes.size();

            // set size for list of dv values
            aDvValues.resize(tNumTypes);

            // loop over the node indices
            for (uint iInd = 0; iInd < tNumIndices; iInd++)
            {
                // loop over the requested dv types
                for (uint iType = 0; iType < tNumTypes; iType++)
                {
                    aDvValues(iType).resize(tNumIndices, 1);
                    aDvValues(iType)(iInd) = mIgPdvHosts(aNodeIndices(iInd))->get_pdv_value(aPdvTypes(iType));
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::get_ig_pdv_value(const Matrix<IndexMat> & aNodeIndices,
                                                const Cell<PDV_Type>   & aPdvTypes,
                                                Cell<Matrix<DDRMat>>   & aDvValues,
                                                Cell<Matrix<DDSMat>>   & aIsActiveDv )
        {
            // get the number of node indices requested
            uint tNumIndices = aNodeIndices.length();

            // get the number of dv types requested
            uint tNumTypes = aPdvTypes.size();

            // set size for list of dv values
            //aDvValues.resize( tNumTypes );

            // set size for list of active flags
            aIsActiveDv.resize( tNumTypes );

            // loop over the requested dv types
            for ( uint iType = 0; iType < tNumTypes; iType++ )
            {
                aIsActiveDv( iType ).resize( tNumTypes, 1 );
                //aDvValues( iType ).resize( tNumTypes, 1 );

                // loop over the node indices
                for ( uint iInd = 0; iInd < tNumIndices; iInd++ )
                {
                    aIsActiveDv( iType )( iInd ) = mIgPdvHosts( aNodeIndices( iInd ) )->is_active_type( aPdvTypes( iType ) );
                    aDvValues( iType )( iInd )   = mIgPdvHosts( aNodeIndices( iInd ) )->get_pdv_value( aPdvTypes( iType ) );
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix<DDSMat> & Pdv_Host_Manager::get_my_local_global_map()
        {
            return mGlobalPdvTypeMap;
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::get_ip_dv_ids_for_type_and_ind(const Matrix<IndexMat>&    aNodeIndices,
                                                              const Cell<PDV_Type>&         aPdvTypes,
                                                              Cell<Matrix<IdMat>>&        aDvIds)
        {
            // get the number of node indices requested
            uint tNumIndices = aNodeIndices.length();

            // get the number of dv types requested
            uint tNumTypes = aPdvTypes.size();

            // set size for list of dv values
            aDvIds.resize(tNumTypes);

            // loop over the node indices
            for (uint iInd = 0; iInd < tNumIndices; iInd++)
            {
                // loop over the requested dv types
                for (uint iType = 0; iType < tNumTypes; iType++)
                {
                    aDvIds(iType).resize(tNumIndices, 1);
                    aDvIds(iType)(iInd) = mIpPdvHosts(aNodeIndices(iInd))->get_global_index_for_pdv_type(aPdvTypes(iType));
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::get_ig_dv_ids_for_type_and_ind(const Matrix<IndexMat>&    aNodeIndices,
                                                              const Cell<PDV_Type>&         aPdvTypes,
                                                              Cell<Matrix<IdMat>>&        aDvIds)
        {
            // get the number of node indices requested
            uint tNumIndices = aNodeIndices.numel();

            // get the number of dv types requested
            uint tNumTypes = aPdvTypes.size();

            // set size for list of dv values
            aDvIds.resize( tNumTypes );

            // loop over the requested dv types
            for ( uint iType = 0; iType < tNumTypes; iType++ )
            {
                aDvIds( iType ).resize( tNumIndices, 1 );

                // loop over the node indices
                for ( uint iInd = 0; iInd < tNumIndices; iInd++ )
                {
                    aDvIds( iType )( iInd ) = mIgPdvHosts( aNodeIndices( iInd ) )->get_global_index_for_pdv_type( aPdvTypes( iType ) );
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

        void Pdv_Host_Manager::create_ip_pdv_hosts(Cell<Matrix<DDSMat>>        aNodeIndicesPerSet,
                                                   Cell<Matrix<DDRMat>>        aNodeCoordinates,
                                                   Cell<Cell<Cell<PDV_Type>>>  aPdvTypes)
        {
            // Check that number of sets is consistent
            uint tNumSets = aPdvTypes.size();
            MORIS_ERROR(tNumSets == aNodeIndicesPerSet.size(),
                    "Information passed to Pdv_Host_Manager.create_ip_pdv_hosts() does not have a consistent number of sets!");

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

                // Copy PDV types over
                uint tUniquePdvIndex = 0;
                for (uint tGroupIndex = 0; tGroupIndex < mIpPdvTypes(tMeshSetIndex).size(); tGroupIndex++)
                {
                    for (uint tPdvIndex = 0; tPdvIndex < mIpPdvTypes(tMeshSetIndex)(tGroupIndex).size(); tPdvIndex++)
                    {
                        mUniqueIpPdvTypes(tMeshSetIndex)(tUniquePdvIndex++) = mIpPdvTypes(tMeshSetIndex)(tGroupIndex)(tPdvIndex);
                    }
                }

                // Create PDV hosts
                for (uint tNodeIndexOnSet = 0; tNodeIndexOnSet < aNodeIndicesPerSet(tMeshSetIndex).length(); tNodeIndexOnSet++)
                {
                    // Create new host or add unique PDVs
                    uint tNumAddedPdvs = tNumUniquePdvs;
                    if (mIpPdvHosts(aNodeIndicesPerSet(tMeshSetIndex)(tNodeIndexOnSet)) == nullptr)
                    {
                        mIpPdvHosts(aNodeIndicesPerSet(tMeshSetIndex)(tNodeIndexOnSet))
                        = std::make_shared<Interpolation_Pdv_Host>(aNodeIndicesPerSet(tMeshSetIndex)(tNodeIndexOnSet),
                                                                   aNodeCoordinates(aNodeIndicesPerSet(tMeshSetIndex)(tNodeIndexOnSet)),
                                                                   mUniqueIpPdvTypes(tMeshSetIndex),
                                                                   mGlobalPdvIndex);
                    }
                    else
                    {
                        tNumAddedPdvs = mIpPdvHosts(aNodeIndicesPerSet(tMeshSetIndex)(tNodeIndexOnSet))->add_pdv_types(mUniqueIpPdvTypes(tMeshSetIndex), mGlobalPdvIndex);
                    }

                    // Resize global map
                    mGlobalPdvTypeMap.resize(mGlobalPdvTypeMap.length() + tNumAddedPdvs, 1);

                    // Update global PDV indices
                    for (uint tPdvIndex = 0; tPdvIndex < tNumAddedPdvs; tPdvIndex++)
                    {
                        mGlobalPdvTypeMap(mGlobalPdvIndex) = mGlobalPdvIndex;
                        mGlobalPdvIndex++;
                    }
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::create_ig_pdv_hosts(Cell<Matrix<DDSMat>>        aNodeIndicesPerSet,
                                                   Cell<Matrix<DDRMat>>        aNodeCoordinates,
                                                   Cell<Cell<Cell<PDV_Type>>>  aPdvTypes,
                                                   Cell<Intersection_Info>     aIntersectionInfo)
        {
            // Check that number of sets is consistent
            uint tNumSets = aPdvTypes.size();
            MORIS_ERROR(tNumSets == aNodeIndicesPerSet.size(),
            "Information passed to Pdv_Host_Manager.create_ig_pdv_hosts() does not have a consistent number of sets!");

            // Set PDV types
            mIgPdvTypes = aPdvTypes;
            mUniqueIgPdvTypes.resize(tNumSets);

            // Initialize PDV hosts
            mIgPdvHosts.resize(aNodeCoordinates.size(), nullptr);

            // Create PDV hosts
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
                        mUniqueIgPdvTypes(tMeshSetIndex)(tUniquePdvIndex++) = mIgPdvTypes(tMeshSetIndex)(tGroupIndex)(tPdvIndex);
                    }
                }

                // Create PDV hosts
                for (uint tNodeIndexOnSet = 0; tNodeIndexOnSet < aNodeIndicesPerSet(tMeshSetIndex).length(); tNodeIndexOnSet++)
                {
                    // Create new host if not created already
                    uint tNumAddedPdvs = tNumUniquePdvs;
                    if (mIgPdvHosts(aNodeIndicesPerSet(tMeshSetIndex)(tNodeIndexOnSet)) == nullptr)
                    {
                        mIgPdvHosts(aNodeIndicesPerSet(tMeshSetIndex)(tNodeIndexOnSet))
                        = std::make_shared<Integration_Pdv_Host>(aNodeIndicesPerSet(tMeshSetIndex)(tNodeIndexOnSet),
                                                                 aNodeCoordinates(aNodeIndicesPerSet(tMeshSetIndex)(tNodeIndexOnSet)),
                                                                 mGlobalPdvIndex);
                        // Resize global map
                        mGlobalPdvTypeMap.resize(mGlobalPdvTypeMap.length() + tNumAddedPdvs, 1);

                        // Update global PDV indices
                        for (uint tPdvIndex = 0; tPdvIndex < tNumAddedPdvs; tPdvIndex++)
                        {
                            mGlobalPdvTypeMap(mGlobalPdvIndex) = mGlobalPdvIndex;
                            mGlobalPdvIndex++;
                        }
                    }
                }
            }

            // Set intersection info
            for (uint tIntersectionIndex = 0; tIntersectionIndex < aIntersectionInfo.size(); tIntersectionIndex++)
            {
                mIgPdvHosts(aIntersectionInfo(tIntersectionIndex).mNodeIndex)->set_as_intersection(
                        this->convert_info_to_intersection(aIntersectionInfo(tIntersectionIndex)));
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::set_ip_requested_dv_types(Cell<PDV_Type>& aPdvTypes)
        {
            mRequestedIpPdvTypes = aPdvTypes;
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::set_ig_requested_dv_types(Cell<PDV_Type>& aPdvTypes)
        {
            mRequestedIgPdvTypes = aPdvTypes;
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::create_ip_pdv(uint aNodeIndex, PDV_Type aPdvType, real aPdvVal)
        {
            mIpPdvHosts(aNodeIndex)->create_pdv(aPdvType, aPdvVal);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::create_ip_pdv(uint aNodeIndex, PDV_Type aPdvType, std::shared_ptr<Property> aPropertyPointer)
        {
            mIpPdvHosts(aNodeIndex)->create_pdv(aPdvType, aPropertyPointer);
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Pdv_Host_Manager::compute_diqi_dadv()
        {
            return this->get_dQidu() * this->compute_dpdv_dadv();
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Pdv_Host_Manager::compute_dpdv_dadv()
        {
            // Information from each host
            Matrix<DDRMat> tTotalAdvSensitivities(0, 0);
            Matrix<DDRMat> tTestSensitivities(0, 0);

            // Need to know the size of the sensitivities
            uint tFirstPdvHostIndex = 0;
            if (mIpPdvHosts.size() > 0)
            {
                while (tTestSensitivities.numel() == 0)
                {
                    mIpPdvHosts(tFirstPdvHostIndex++)->get_all_sensitivities(tTestSensitivities);
                }
            }
            else if (mIgPdvHosts.size() > 0)
            {
                while (tTestSensitivities.numel() == 0)
                {
                    mIgPdvHosts(tFirstPdvHostIndex++)->get_all_sensitivities(tTestSensitivities);
                }
            }
            tTotalAdvSensitivities.resize(mGlobalPdvIndex, tTestSensitivities.n_cols());

            // Loop over PDV hosts
            for (uint tPdvHostIndex = 0; tPdvHostIndex < mIpPdvHosts.size(); tPdvHostIndex++)
            {
                // Get sensitivities
                Matrix<DDRMat> tHostAdvSensitivities(0, 0);
                mIpPdvHosts(tPdvHostIndex)->get_all_sensitivities(tHostAdvSensitivities);
                const Matrix<DDUMat>& tHostPdvIndices = mIpPdvHosts(tPdvHostIndex)->get_all_global_indices();

                // Add to total matrix
                for (uint tRowIndex = 0; tRowIndex < tHostAdvSensitivities.n_rows(); tRowIndex++)
                {
                    tTotalAdvSensitivities.set_row(tHostPdvIndices(tRowIndex), tHostAdvSensitivities.get_row(tRowIndex));
                }
            }
            for (uint tPdvHostIndex = 0; tPdvHostIndex < mIgPdvHosts.size(); tPdvHostIndex++)
            {
                // Get sensitivities
                Matrix<DDRMat> tHostAdvSensitivities(0, 0);
                mIgPdvHosts(tPdvHostIndex)->get_all_sensitivities(tHostAdvSensitivities);
                uint tStartingGlobalIndex = mIgPdvHosts(tPdvHostIndex)->get_starting_global_index();

                // Add to total matrix
                for (uint tRowIndex = 0; tRowIndex < tHostAdvSensitivities.n_rows(); tRowIndex++)
                {
                    tTotalAdvSensitivities.set_row(tStartingGlobalIndex + tRowIndex, tHostAdvSensitivities.get_row(tRowIndex));
                }
            }
            return tTotalAdvSensitivities;
        }

        //--------------------------------------------------------------------------------------------------------------

        std::shared_ptr<Intersection> Pdv_Host_Manager::convert_info_to_intersection(Intersection_Info aIntersectionInfo)
        {
            std::shared_ptr<Intersection> tIntersection = std::make_shared<Intersection>();

            // Copy over geometry
            tIntersection->mGeometry = aIntersectionInfo.mGeometry;

            // Nodes
            tIntersection->mNodes.resize(aIntersectionInfo.mParentNodeIndices.length());
            for (uint tNodeIndex = 0; tNodeIndex < aIntersectionInfo.mParentNodeIndices.length(); tNodeIndex++)
            {
                tIntersection->mNodes(tNodeIndex) = mIgPdvHosts(aIntersectionInfo.mParentNodeIndices(tNodeIndex));
            }

            return tIntersection;
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
