#include "cl_GEN_Pdv_Host_Manager.hpp"
#include "cl_Communication_Tools.hpp"
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

        void Pdv_Host_Manager::set_num_advs(uint aNumADVs)
        {
            mNumADVs = aNumADVs;
            mNumADVsSet = true;
        }

        //--------------------------------------------------------------------------------------------------------------
        
        void Pdv_Host_Manager::set_num_background_nodes(uint aNumNodes)
        {
            mIntersectionNodes.resize(aNumNodes);
            mNumStaticNodesSet = true;
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::reset()
        {
            mIntersectionNodes.resize(0);
            mGlobalPdvTypeMap.resize(0, 0);
            mGlobalPdvIndex = 0;
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

        void Pdv_Host_Manager::get_ip_pdv_value(const Matrix<IndexMat>& aNodeIndices,
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
            // MORIS_ERROR(par_size() == 1, "PDV Host Manager local/global map will not work in parallel.");
            return mGlobalPdvTypeMap;
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
                                mIntersectionNodes(aNodeIndices(tNode))->get_starting_pdv_index()
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

        void Pdv_Host_Manager::create_interpolation_pdv_hosts(Cell<Matrix<DDSMat>>       aNodeIndicesPerSet,
                                                   Cell<Matrix<DDRMat>>       aNodeCoordinates,
                                                   Cell<Cell<Cell<PDV_Type>>> aPdvTypes)
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
            MORIS_ASSERT(mNumStaticNodesSet,
                    "Number of background nodes must be set in the PDV Host Manager before intersection nodes can be created.");
            MORIS_ASSERT(aNodeIndex == mIntersectionNodes.size(),
                    "Intersection nodes must be added to the PDV Host Manager in order by node index.");

            // Add intersection node
            mIntersectionNodes.push_back(aIntersectionNode);

            // Check for if this should be a PDV
            if (mIntersectionNodes(aNodeIndex))
            {
                // Set global index
                mIntersectionNodes(aNodeIndex)->set_starting_pdv_index(mGlobalPdvIndex);

                // Number of PDVs being added
                uint tNumAddedPdvs = aIntersectionNode->get_global_coordinates().length();

                // Resize global map
                mGlobalPdvTypeMap.resize(mGlobalPdvTypeMap.length() + tNumAddedPdvs, 1);

                // Update global index
                for (uint tPdvIndex = 0; tPdvIndex < tNumAddedPdvs; tPdvIndex++)
                {
                    mGlobalPdvTypeMap(mGlobalPdvIndex) = mGlobalPdvIndex;
                    mGlobalPdvIndex++;
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::set_num_integration_nodes(uint aNumNodes)
        {
            mIntersectionNodes.resize(aNumNodes);
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

        void Pdv_Host_Manager::create_interpolation_pdv(uint aNodeIndex, PDV_Type aPdvType, std::shared_ptr<Property> aPropertyPointer)
        {
            mIpPdvHosts(aNodeIndex)->create_pdv(aPdvType, aPropertyPointer);
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Pdv_Host_Manager::compute_diqi_dadv()
        {
            return this->get_dQIdp() * this->compute_dpdv_dadv();
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Pdv_Host_Manager::compute_dpdv_dadv()
        {
            // Check for ADVs
            MORIS_ERROR(mNumADVsSet, "PDV Host Manager must have the number of ADVs set before computing sensitivities.");

            // Total matrix dpdv/dadv
            Matrix<DDRMat> tTotalAdvSensitivities(mGlobalPdvIndex, mNumADVs, 0.0);

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
            for (uint tIntersectionIndex = 0; tIntersectionIndex < mIntersectionNodes.size(); tIntersectionIndex++)
            {
                if (mIntersectionNodes(tIntersectionIndex))
                {
                    // Get sensitivities
                    Matrix<DDRMat> tHostAdvSensitivities(0, 0);
                    mIntersectionNodes(tIntersectionIndex)->get_all_sensitivities(tHostAdvSensitivities);
                    uint tStartingGlobalIndex = mIntersectionNodes(tIntersectionIndex)->get_starting_pdv_index();

                    // Add to total matrix
                    for (uint tRowIndex = 0; tRowIndex < tHostAdvSensitivities.n_rows(); tRowIndex++)
                    {
                        tTotalAdvSensitivities.set_row(tStartingGlobalIndex + tRowIndex, tHostAdvSensitivities.get_row(tRowIndex));
                    }
                }
            }
            return tTotalAdvSensitivities;
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
