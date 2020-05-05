#include "cl_GEN_Pdv_Host_Manager.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Pdv_Host_Manager::Pdv_Host_Manager()
        {
//            // set size for IP dv type enum to index map
//            mIpGlobalPdvTypeMap.set_size(static_cast<size_t>(PDV::UNDEFINED), 1, gNoIndex);
//            // set size for IG dv type enum to index map
//            mIgGlobalPdvTypeMap.set_size(static_cast<size_t>(PDV::Z_COORDINATE) + 1, 1, gNoIndex);
        }

        //--------------------------------------------------------------------------------------------------------------

        Pdv_Host_Manager::~Pdv_Host_Manager()
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::get_ip_dv_types_for_set(const moris::moris_index aIPMeshSetIndex, Cell<Cell<PDV>>& aDvTypes)
        {
            aDvTypes = mIpPdvTypes(aIPMeshSetIndex);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::get_ig_dv_types_for_set(const moris::moris_index aIGMeshSetIndex, Cell<Cell<PDV>>& aDvTypes)
        {
            aDvTypes = mIgPdvTypes(aIGMeshSetIndex);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::get_ip_unique_dv_types_for_set(const moris_index aIPMeshSetIndex, Cell<PDV>& aDvTypes)
        {
            // Get number of unique PDVs
            uint tNumUniquePdvs = 0;
            for (uint tGroupIndex = 0; tGroupIndex < mIpPdvTypes(aIPMeshSetIndex).size(); tGroupIndex++)
            {
                tNumUniquePdvs += mIpPdvTypes(aIPMeshSetIndex)(tGroupIndex).size();
            }

            // Copy PDV types over
            uint tTotalPdvIndex = 0;
            for (uint tGroupIndex = 0; tGroupIndex < mIpPdvTypes(aIPMeshSetIndex).size(); tGroupIndex++)
            {
                for (uint tPdvIndex = 0; tPdvIndex < mIpPdvTypes(aIPMeshSetIndex)(tGroupIndex).size(); tPdvIndex++)
                {
                    aDvTypes(tTotalPdvIndex++) = mIpPdvTypes(aIPMeshSetIndex)(tGroupIndex)(tPdvIndex);
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::get_ig_unique_dv_types_for_set(const moris::moris_index aIGMeshSetIndex, Cell<PDV>& aDvTypes)
        {
            // Get number of unique PDVs
            uint tNumUniquePdvs = 0;
            for (uint tGroupIndex = 0; tGroupIndex < mIgPdvTypes(aIGMeshSetIndex).size(); tGroupIndex++)
            {
                tNumUniquePdvs += mIgPdvTypes(aIGMeshSetIndex)(tGroupIndex).size();
            }

            // Copy PDV types over
            uint tTotalPdvIndex = 0;
            for (uint tGroupIndex = 0; tGroupIndex < mIgPdvTypes(aIGMeshSetIndex).size(); tGroupIndex++)
            {
                for (uint tPdvIndex = 0; tPdvIndex < mIgPdvTypes(aIGMeshSetIndex)(tGroupIndex).size(); tPdvIndex++)
                {
                    aDvTypes(tTotalPdvIndex++) = mIgPdvTypes(aIGMeshSetIndex)(tGroupIndex)(tPdvIndex);
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::get_ip_pdv_value(const Matrix<IndexMat>&   aNodeIndices,
                                                const Cell<PDV>&       aDvTypes,
                                                Cell<Matrix<DDRMat>>&     aDvValues)
        {
            // get the number of node indices requested
            uint tNumIndices = aNodeIndices.length();

            // get the number of dv types requested
            uint tNumTypes = aDvTypes.size();

            // set size for list of dv values
            aDvValues.resize(tNumTypes);

            // loop over the node indices
            for (uint iInd = 0; iInd < tNumIndices; iInd++)
            {
                // loop over the requested dv types
                for (uint iType = 0; iType < tNumTypes; iType++)
                {
                    aDvValues(iType).resize(tNumIndices, 1);
                    aDvValues(iType)(iInd) = this->get_ip_pdv_by_type_and_index(aNodeIndices(iInd), aDvTypes(iType))->get_val()(0, 0);
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::get_ip_pdv_value(const Matrix<IndexMat>&   aNodeIndices,
                                                const Cell<PDV>&       aDvTypes,
                                                Cell<Matrix<DDRMat>>&     aDvValues,
                                                Cell<Matrix<DDSMat>>&     aIsActiveDv)
        {
            get_ip_pdv_value(aNodeIndices, aDvTypes, aDvValues);

            // get the number of node indices requested
            uint tNumIndices = aNodeIndices.length();

            // get the number of dv types requested
            uint tNumTypes = aDvTypes.size();

            // set size for list of active flags
            aIsActiveDv.resize(tNumIndices);

            // loop over the node indices
            for (uint iInd = 0; iInd < tNumIndices; iInd++)
            {
                aIsActiveDv(iInd).resize(tNumTypes, 1);

                // loop over the requested dv types
                for (uint iType = 0; iType < tNumTypes; iType++)
                {
                    aIsActiveDv(iInd)(iType) = this->check_ip_for_active_type(aNodeIndices(iInd), aDvTypes(iType));
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::get_ig_pdv_value(const Matrix<IndexMat>&   aNodeIndices,
                                                const Cell<PDV>&       aDvTypes,
                                                Cell<Matrix<DDRMat>>&     aDvValues)
        {
            // get the number of node indices requested
            uint tNumIndices = aNodeIndices.length();

            // get the number of dv types requested
            uint tNumTypes = aDvTypes.size();

            // set size for list of dv values
            aDvValues.resize(tNumTypes);

            // loop over the node indices
            for (uint iInd = 0; iInd < tNumIndices; iInd++)
            {
                // loop over the requested dv types
                for (uint iType = 0; iType < tNumTypes; iType++)
                {
                    aDvValues(iType).resize(tNumIndices, 1);
                    aDvValues(iType)(iInd) = this->get_ip_pdv_by_type_and_index(aNodeIndices(iInd), aDvTypes(iType))->get_val()(0, 0);
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::get_ig_pdv_value(const Matrix<IndexMat>&   aNodeIndices,
                                                const Cell<PDV>&       aDvTypes,
                                                Cell<Matrix<DDRMat>>&     aDvValues,
                                                Cell<Matrix<DDSMat>>&     aIsActiveDv)
        {
            // get the number of node indices requested
            uint tNumIndices = aNodeIndices.length();

            // get the number of dv types requested
            uint tNumTypes = aDvTypes.size();

            // set size for list of active flags
            aIsActiveDv.resize(tNumIndices);

            // loop over the node indices
            for (uint iInd = 0; iInd < tNumIndices; iInd++)
            {
                aIsActiveDv(iInd).resize(tNumTypes, 1);

                // loop over the requested dv types
                for (uint iType = 0; iType < tNumTypes; iType++)
                {
                    aIsActiveDv(iInd)(iType) = this->check_ig_for_active_type(aNodeIndices(iInd), aDvTypes(iType));
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDSMat> Pdv_Host_Manager::get_ip_local_global_map()
        {
            return mIpGlobalPdvTypeMap;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDSMat> Pdv_Host_Manager::get_ig_local_global_map()
        {
            return mIgGlobalPdvTypeMap;
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::get_ip_dv_ids_for_type_and_ind(const Cell<moris_index>&    aNodeIndices,
                                                              const Cell<PDV>&         aDvTypes,
                                                              Cell<Matrix<IdMat>>&        aDvIds)
        {
            /*
             * - each cell is a row vector of global IDs per each type
             * - return the global ids of the dv type on a specified vertex
             */
        
            uint tNumIndices = aNodeIndices.size();
            uint tNumTypes   = aDvTypes.size();
        
            moris::Cell< uint > tCounter(tNumTypes, 0);
        
            aDvIds.resize(tNumTypes);
            for (uint Ik = 0; Ik <tNumTypes; Ik++)
            {
                aDvIds(Ik).set_size(tNumIndices, 1);
            }
        
            for (uint iType = 0; iType<tNumTypes; iType++)
            {
                for (uint iInd = 0; iInd<tNumIndices; iInd++)
                {
        
                    bool tDvTypeExists = this->check_ip_for_active_type(aNodeIndices(iInd), aDvTypes(iType));   // flag for if the DV type exists on the current host
        
                    if (tDvTypeExists)
                    {
                        aDvIds(iType)(iInd) = this->get_ip_global_index_for_dv_type(aNodeIndices(iInd), aDvTypes(iType));
                        tCounter(iType)++;
                    }
                }
            }
        
            for (uint Ik = 0; Ik < tNumTypes; Ik++)
            {
                aDvIds(Ik).resize(tCounter(Ik), 1);
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::get_ig_dv_ids_for_type_and_ind(const Cell<moris_index>&    aNodeIndices,
                                                              const Cell<PDV>&         aDvTypes,
                                                              Cell<Matrix<IdMat>>&        aDvIds)
        {
            /*
             * - each cell is a row vector of global IDs per each type
             * - return the global ids of the dv type on a specified vertex
             */

            uint tNumIndices = aNodeIndices.size();
            uint tNumTypes   = aDvTypes.size();

            Cell<uint> tCounter(tNumTypes, 0);

            aDvIds.resize(tNumTypes);
            for (uint Ik = 0; Ik < tNumTypes; Ik++)
            {
                aDvIds(Ik).set_size(tNumIndices, 1);
            }

            for (uint iType = 0; iType < tNumTypes; iType++)
            {
                for (uint iInd = 0; iInd < tNumIndices; iInd++)
                {

                    bool tDvTypeExists = this->check_ig_for_active_type(aNodeIndices(iInd), aDvTypes(iType));   // flag for if the DV type exists on the current host

                    if (tDvTypeExists)
                    {
                        aDvIds(iType)(iInd) = this->get_ig_global_index_for_dv_type(aNodeIndices(iInd), aDvTypes(iType));
                        tCounter(iType)++;
                    }
                }
            }

            for (uint Ik = 0; Ik < tNumTypes; Ik++)
            {
                aDvIds(Ik).resize(tCounter(Ik), 1);
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::initialize_ip_hosts(uint aTotalNumVertices)
        {
            // set size for the list of pdv host
            mIpPdvHosts.resize(aTotalNumVertices, nullptr);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::initialize_ig_hosts(uint aTotalNumVertices)
        {
            // set size for the list of pdv host
            mIgPdvHosts.resize(aTotalNumVertices, nullptr);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::create_ip_pdv_host(uint aNumPdvs, moris_index aNodeIndex)
        {
            // if pdv host not assigned yet
            if (mIpPdvHosts(aNodeIndex) == nullptr)
            {
                // create a pdv host
                mIpPdvHosts(aNodeIndex) = std::make_shared<GEN_Pdv_Host>(aNumPdvs, aNodeIndex);
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::create_ig_pdv_host(uint aNumPdvs, moris_index aNodeIndex)
        {
            // if pdv host is not assigned yet
            if (mIgPdvHosts(aNodeIndex) == nullptr)
            {
                // create a pdv host
                mIgPdvHosts(aNodeIndex) = std::make_shared<GEN_Pdv_Host>(aNumPdvs, aNodeIndex);
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::assign_property_to_pdv_type_by_vertex_index(std::shared_ptr<GEN_Property> aPropertyPointer,
                                                                           PDV                        aPdvType,
                                                                           moris_index                   aNodeIndex)
        {
            // create a pdv host for vertex index
            this->create_ip_pdv_host(mIpNumPDVs, aNodeIndex);
        
            // get the pdv host and create the pdv for dv type
            mIpPdvHosts(aNodeIndex)->create_pdv(aPropertyPointer, aPdvType, mIpGlobalPdvTypeMap);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::assign_field_to_pdv_type_by_vertex_index(std::shared_ptr<GEN_Field>   aFieldPointer,
                                                                        PDV                       aPdvType,
                                                                        moris_index                  aNodeIndex)
        {
            // create a pdv host for vertex index
            this->create_ip_pdv_host(mIpNumPDVs, aNodeIndex);
        
            // get the pdv host and create the pdv for dv type
            mIpPdvHosts(aNodeIndex)->create_pdv(aFieldPointer, aPdvType, mIpGlobalPdvTypeMap);
        }

        //--------------------------------------------------------------------------------------------------------------

        std::shared_ptr<GEN_Pdv> Pdv_Host_Manager::get_ip_pdv_by_type_and_index(moris_index aNodeIndex, PDV aPdvType)
        {
            return mIpPdvHosts(aNodeIndex)->get_pdv_by_type(aPdvType, mIpGlobalPdvTypeMap);
        }

        //--------------------------------------------------------------------------------------------------------------

        std::shared_ptr<GEN_Pdv> Pdv_Host_Manager::get_ig_pdv_by_type_and_index(moris_index aNodeIndex, PDV aPdvType)
        {
            return mIpPdvHosts(aNodeIndex)->get_pdv_by_type(aPdvType, mIgGlobalPdvTypeMap);
        }

        //--------------------------------------------------------------------------------------------------------------

        sint Pdv_Host_Manager::check_ip_for_active_type(moris_index aNodeIndex, PDV aPdvType)
        {
            return mIpPdvHosts(aNodeIndex)->is_active_type(aPdvType, mIpGlobalPdvTypeMap);
        }

        //--------------------------------------------------------------------------------------------------------------

        sint Pdv_Host_Manager::check_ig_for_active_type(moris_index aNodeIndex, PDV aPdvType)
        {
            return mIgPdvHosts(aNodeIndex)->is_active_type(aPdvType, mIgGlobalPdvTypeMap);
        }

        //--------------------------------------------------------------------------------------------------------------

        uint Pdv_Host_Manager::get_ip_global_index_for_dv_type(moris_index aNodeIndex, PDV aPdvType)
        {
            return mIpPdvHosts(aNodeIndex)->get_global_index_for_dv_type(aPdvType, mIpGlobalPdvTypeMap);
        }

        //--------------------------------------------------------------------------------------------------------------

        uint Pdv_Host_Manager::get_ig_global_index_for_dv_type(moris_index aNodeIndex, PDV aPdvType)
        {
            return mIgPdvHosts(aNodeIndex)->get_global_index_for_dv_type(aPdvType, mIgGlobalPdvTypeMap);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::set_ip_pdv_types(Cell<Cell<Cell<PDV>>> aPdvTypeList)
        {
            mIpPdvTypes = aPdvTypeList;
//            // communicate dv types
//            this->communicate_ip_dv_types(aPdvTypeList);
//
//            // get number of dv types
//            uint tNumTypes = mIpPdvTypes.size();
//
//            // loop over dv types
//            for(uint i = 0; i < tNumTypes; i++)
//            {
//                // populate the dv type to index map
//                mIpGlobalPdvTypeMap(static_cast<sint>(mIpPdvTypes(i))) = mIpNumPDVs++;
//            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::set_ig_pdv_types(Cell<Cell<Cell<PDV>>> aPdvTypeList)
        {
            mIgPdvTypes = aPdvTypeList;
//            // communicate dv types
//            this->communicate_ig_dv_types(aPdvTypeList);
//
//            // get number of dv types
//            uint tNumTypes = mIgPdvTypes.size();
//
//            // loop over dv types
//            for(uint i = 0; i < tNumTypes; i++)
//            {
//                // populate the dv type to index map
//                mIgGlobalPdvTypeMap(static_cast<sint>(mIgPdvTypes(i))) = mIgNumPDVs++;
//            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::mark_ip_pdv_as_inactive(moris_index aNodeIndex, PDV aPdvType)
        {
            mIpPdvHosts(aNodeIndex)->mark_pdv_as_inactive(aPdvType, mIpGlobalPdvTypeMap);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::mark_ig_pdv_as_inactive(moris_index aNodeIndex, PDV aPdvType)
        {
            mIgPdvHosts(aNodeIndex)->mark_pdv_as_inactive(aPdvType, mIgGlobalPdvTypeMap);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::update_ip_local_to_global_dv_type_map() // FIXME
        {
//            // Get number of unique PDVs
//            uint tNumUniquePdvs = 0;
//            for (uint tGroupIndex = 0; tGroupIndex < mIpPdvTypes(0).size(); tGroupIndex++)
//            {
//                tNumUniquePdvs += mIpPdvTypes(0)(tGroupIndex).size();
//            }
//
//            // get the number of pdv hosts
//            uint tNumHosts = mIpPdvHosts.size();
//
//            // loop over the dv types
//            for(uint iType = 0; iType < tNumUniquePdvs; iType++)
//            {
//                // loop over the pdv hosts
//                for(uint iHost = 0; iHost < tNumHosts; iHost++)
//                {
//                    // if the pdv host was assigned
//                    if(mIpPdvHosts(iHost) != nullptr)
//                    {
//                        // assign a global id to type
//                        mIpPdvHosts(iHost)->assign_id_to_type(mIpGlobalID, mIpPdvTypes(iType), mIpGlobalPdvTypeMap);
//
//                        // update global id counter
//                        mIpGlobalID++;
//                    }
//                }
//            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host_Manager::update_ig_local_to_global_dv_type_map()
        {
//            // Get number of unique PDVs
//            uint tNumUniquePdvs = 0;
//            for (uint tGroupIndex = 0; tGroupIndex < mIpPdvTypes(0).size(); tGroupIndex++)
//            {
//                tNumUniquePdvs += mIpPdvTypes(0)(tGroupIndex).size();
//            }
//
//            // get the number of pdv hosts
//            uint tNumHosts = mIgPdvHosts.size();
//
//            // loop over the dv types
//            for (uint iType = 0; iType < tNumUniquePdvs; iType++)
//            {
//                // loop over the pdv hosts
//                for (uint iHost = 0; iHost < tNumHosts; iHost++)
//                {
//                    // if the pdv host was assigned
//                    if (mIgPdvHosts(iHost) != nullptr)
//                    {
//                        // assign a global id to type
//                        mIgPdvHosts(iHost)->assign_id_to_type(mIgGlobalID, mIgPdvTypes(iType), mIgGlobalPdvTypeMap);
//
//                        // update global id counter
//                        mIgGlobalID++;
//                    }
//                }
//            }
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}