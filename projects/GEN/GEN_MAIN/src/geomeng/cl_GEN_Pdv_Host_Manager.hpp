#ifndef CL_GEN_PDV_HOST_MANAGER_HPP_
#define CL_GEN_PDV_HOST_MANAGER_HPP_

#include "cl_MSI_Design_Variable_Interface.hpp"

// GEN_MAIN
#include "cl_GEN_Pdv_Host.hpp"
#include "cl_GEN_Field.hpp"
// GEN_CORE
#include "cl_GEN_Dv_Enums.hpp"

// CORE
#include "cl_Matrix.hpp"
#include "cl_Communication_Tools.hpp"
#include "cl_Communication_Manager.hpp"
#include "cl_Map.hpp"
#include "fn_sum.hpp"

namespace moris
{
    namespace ge
    {
        class Pdv_Host_Manager : public MSI::Design_Variable_Interface
        {
        private:
            // list of pdv hosts - interpolation nodes
            Cell<std::shared_ptr<GEN_Pdv_Host>> mIpPdvHosts;
            // list of pdv hosts - integration nodes
            Cell<std::shared_ptr<GEN_Pdv_Host>> mIgPdvHosts;
            
            // Cell of vertex indices per mesh set
            Cell<Matrix<DDSMat>> mVertexIndicesOnSet;
            
            // position in map corresponds to the value of the pdv enum
            Matrix<IndexMat> mIpGlobalPdvTypeMap;
            // position in map corresponds to the value of the pdv enum
            Matrix<IndexMat> mIgGlobalPdvTypeMap;
            
            // list containing all the used pdv types
            Cell<Cell<Cell<GEN_DV>>> mIpPdvTypes;
            
            // list containing all the used pdv types
            Cell<Cell<Cell<GEN_DV>>> mIgPdvTypes;
            
            // total number of dv types on interpolation nodes
            uint mIpNumPDVs = 0;
            // total number of dv types on integration nodes
            uint mIgNumPDVs = 0;
            
            // global dv ID for interpolation nodes
            uint mIpGlobalID = 0;
            // global dv ID for integration nodes
            uint mIgGlobalID = 0;
            
        public:
            
            /**
             * constructor
             */
            Pdv_Host_Manager();
            
            /**
             * trivial destructor
             */
            ~Pdv_Host_Manager(){};
            
            /**
             * get dv types for set
             * @param aIPMeshSetIndex integration mesh index
             * @param aDvTypes        list of groups of dv types to fill
             */
            void get_ip_dv_types_for_set(const moris_index aIGMeshSetIndex, Cell<Cell<GEN_DV>>& aDvTypes);
            
            /**
             * get dv types for set
             * @param aIGMeshSetIndex integration mesh index
             * @param aDvTypes        list of groups of dv types to fill
             */
            void get_ig_dv_types_for_set(const moris_index aIGMeshSetIndex, Cell<Cell<GEN_DV>>& aDvTypes);
            
            /**
             * get unique dv types for set
             * @param aIPMeshSetIndex integration mesh index
             * @param aDvTypes        list dv types to fill
             */
            void get_ip_unique_dv_types_for_set(const moris_index aIGMeshSetIndex, Cell<GEN_DV>& aDvTypes);
            
            /**
             * get unique dv types for set
             * @param aIGMeshSetIndex integration mesh index
             * @param aDvTypes        list dv types to fill
             */
            void get_ig_unique_dv_types_for_set(const moris_index aIGMeshSetIndex, Cell<GEN_DV>& aDvTypes);
            
            /**
             * get pdv values for requested vertex indices and dv types
             * @param aNodeIndices list of node indices
             * @param aDvTypes     list of dv types
             * @param aDvValues    list of returned dv values (DvType)(vertexIndex)
             */
            void get_ip_pdv_value(const Matrix<IndexMat>&   aNodeIndices,
                                  const Cell<GEN_DV>&       aDvTypes,
                                  Cell<Matrix<DDRMat>>&     aDvValues);
            
            /**
             * get pdv values for requested vertex indices and dv types
             * @param aNodeIndices list of node indices
             * @param aDvTypes     list of dv types
             * @param aDvValues    list of returned dv values (DvType)(vertexIndex)
             * @param aIsActive    list of if design variable is active (vertexIndex)(DvType)
             */
            void get_ip_pdv_value(const Matrix<IndexMat>&   aNodeIndices,
                                  const Cell<GEN_DV>&       aDvTypes,
                                  Cell<Matrix<DDRMat>>&     aDvValues,
                                  Cell<Matrix<DDSMat>>&     aIsActiveDv);
            
            /**
             * get pdv values for requested vertex indices and dv types
             * @param aNodeIndices list of node indices
             * @param aDvTypes     list of dv types
             * @param aDvValues    list of dv values (DvType)(vertexIndex)
             */
            void get_ig_pdv_value(const Matrix<IndexMat>&   aNodeIndices,
                                  const Cell<GEN_DV>&       aDvTypes,
                                  Cell<Matrix<DDRMat>>&     aDvValues);
            
            /**
             * get pdv values for requested vertex indices and dv types
             * @param aNodeIndices list of node indices
             * @param aDvTypes     list of dv types
             * @param aDvValues    list of dv values (DvType)(vertexIndex)
             * @param aIsActive    list of active design variables (vertexIndex)(DvType)
             */
            void get_ig_pdv_value(const Matrix<IndexMat>&   aNodeIndices,
                                  const Cell<GEN_DV>&       aDvTypes,
                                  Cell<Matrix<DDRMat>>&     aDvValues,
                                  Cell<Matrix<DDSMat>>&     aIsActiveDv);
            
            /**
             * Get the local to global dv type map for interpolation nodes
             *
             * @return Matrix map from dv type to index
             */
            Matrix<DDSMat> get_ip_local_global_map();
            
            /**
             * Get the local to global dv type map for integration nodes
             *
             * @return Matrix map from dv type to index
             */
            Matrix<DDSMat> get_ig_local_global_map();
            
            /**
             * Return local to global DV type map
             *
             * @param aNodeIndex   List of vertex indices
             * @param aDvType        List of Dv types
             * @param aDvIds         List of Dv Ids
             */
            void get_ip_dv_ids_for_type_and_ind(const Cell<moris_index>&    aNodeIndices,
                                                const Cell<GEN_DV>&         aDvTypes,
                                                Cell<Matrix<IdMat>>&        aDvIds);
            
            /**
             * Get local to global DV type map
             *
             * @param aNodeIndex   List of vertex indices
             * @param aDvType        List of Dv types
             * @param aDvIds         List of Dv Ids
             */
            void get_ig_dv_ids_for_type_and_ind(const Cell<moris_index>&    aNodeIndices,
                                                const Cell<GEN_DV>&         aDvTypes,
                                                Cell<Matrix<IdMat>>&        aDvIds);

            /**
             * get requested dv types for sensitivity analysis
             * @param[ in ] aDvTypes list of dv types to fill
             */
            virtual void get_ip_requested_dv_types( Cell< enum GEN_DV > & aDvTypes )
            {
                MORIS_ERROR(false, "Not implemented yet.");
            }

            /**
             * get requested dv types for sensitivity analysis
             * @param[ in ] aDvTypes list of dv types to fill
             */
            virtual void get_ig_requested_dv_types( Cell< enum GEN_DV > & aDvTypes )
            {
                MORIS_ERROR(false, "Not implemented yet.");
            }
            
            /**
             * initialize the list of interpolation pdv hosts
             * 
             * @param aTotalNumVertices number of vertices in IP mesh
             */
            void initialize_ip_hosts(uint aTotalNumVertices);
            
            /**
             * initialize the list of integration pdv hosts
             * 
             * @param aTotalNumVertices number of vertices in IG mesh
             */
            void initialize_ig_hosts(uint aTotalNumVertices);
            
            /**
             * Create PDV host on interpolation vertex
             * 
             * @param aNumPdvs     a number of dv
             * @param aNodeIndex a node index
             */
            void create_ip_pdv_host(uint aNumPdvs, moris_index aNodeIndex);
            
            /**
             * Create PDV host on integration vertex
             * 
             * @param aNumPdvs     a number of dv
             * @param aNodeIndex a node index
             */
            void create_ig_pdv_host(uint aNumPdvs, moris_index aNodeIndex);
            
            /**
             * Assign a GEN property to pdv type by node index
             * 
             * @param aPropertyPointer a GEN property pointer
             * @param aPdvType         a list of dv types
             * @param aNodeIndex     a node index
             */
            void assign_property_to_pdv_type_by_vertex_index(std::shared_ptr<GEN_Property>  aPropertyPointer,
                                                              GEN_DV                        aPdvType,
                                                              moris_index                   aNodeIndex);
            
            /**
             * Assign a GEN Field to pdv type by node index
             * 
             * @param aFieldPointer a GEN Field pointer
             * @param aPdvType      a list of dv types
             * @param aNodeIndex  a node index
             */
            void assign_field_to_pdv_type_by_vertex_index(std::shared_ptr<GEN_Field>   aFieldPointer,
                                                          GEN_DV                       aPdvType,
                                                          moris_index                  aNodeIndex);
            
            /**
             * Get pdv by type and node index
             *
             * @param aNodeIndex     a node index
             * @param aPdvType         a list of dv types
             */
            std::shared_ptr<GEN_Pdv> get_ip_pdv_by_type_and_index(moris_index aNodeIndex, GEN_DV aPdvType);

            /**
             * Get pdv by type and node index
             *
             * @param aNodeIndex     a node index
             * @param aPdvType         a list of dv types
             */
            std::shared_ptr<GEN_Pdv> get_ig_pdv_by_type_and_index(moris_index aNodeIndex, GEN_DV aPdvType);

            /**
             * check for active pdv by type and node index
             * @param aNodeIndex     a node index
             * @param aPdvType         DV type
             */
            sint check_ip_for_active_type(moris_index aNodeIndex, GEN_DV aPdvType);

            /**
             * check for active pdv by type and node index
             * @param aNodeIndex     a node index
             * @param aPdvType         DV type
             */
            sint check_ig_for_active_type(moris_index aNodeIndex, GEN_DV aPdvType);

            /**
             * get global index for dv type
             * @param aNodeIndex a node index
             * @param aPdvType     a list of dv types
             */
            uint get_ip_global_index_for_dv_type(moris_index aNodeIndex, GEN_DV aPdvType);

            /**
             * get global index for dv type
             * @param aNodeIndex a node index
             * @param aPdvType     a list of dv types
             */
            uint get_ig_global_index_for_dv_type(moris_index aNodeIndex, GEN_DV aPdvType);

            /**
             * set dv types
             * @param aPdvTypeList list of dv types
             */
            void set_ip_pdv_types(Cell<Cell<Cell<GEN_DV>>> aPdvTypeList);

            /**
             * set dv types
             * @param aPdvTypeList list of dv types
             */
            void set_ig_pdv_types(Cell<Cell<Cell<GEN_DV>>> aPdvTypeList);

            /**
             * Mark a PDV on an interpolation node as being inactive
             *
             * @param aNodeIndex IP node index
             * @param aPdvType PDV on the node to be marked
             */
            void mark_ip_pdv_as_inactive(moris_index aNodeIndex, GEN_DV aPdvType);

            /**
             * Mark a PDV on an integration node as being inactive
             *
             * @param aNodeIndex IG node index
             * @param aPdvType PDV on the node to be marked
             */
            void mark_ig_pdv_as_inactive(moris_index aNodeIndex, GEN_DV aPdvType);

//            /**
//             * communicate dv types
//             * @param aPdvTypeList a local list of dv types
//             */
//            void communicate_ip_dv_types(moris::Cell<GEN_DV>& aPdvTypeList)
//            {
//                // get processor size
//                int tSize = par_size();
//
//                // get number of local dv types
//                moris::sint tNumLocalDvTypes = aPdvTypeList.size();
//
//                // variable for maximal possible global dv types ???
//                moris::sint tNumMaxGlobalDvTypes;
//
//                // get number of global dv types
//                sum_all(tNumLocalDvTypes, tNumMaxGlobalDvTypes);
//
//                if (par_rank() == 0)
//                {
//                    // set size of dv type list = number of global types
//                    mIpPdvTypes.resize(tNumMaxGlobalDvTypes);
//                }
//
//                // create list containing the number of local dv types
//                moris::Cell <moris::sint> tNumLocalDvTypesList (tSize);
//
//                // insert number of local dv types into list containing the number of local dv types
//                MPI_Allgather(&tNumLocalDvTypes, 1, MPI_UNSIGNED, (tNumLocalDvTypesList.data()).data(), 1, MPI_UNSIGNED,  MPI_COMM_WORLD);
//
//                // create list containing the offsets of the local dv types in relation to processor 0
//                moris::Cell<moris::sint> tDofTypeOffset(tSize, 0);
//
//                // fill the list with the corresponding offsets
//                for (int Ip = 1; Ip <tSize; ++Ip)
//                {
//                    tDofTypeOffset(Ip) = tDofTypeOffset(Ip-1) + tNumLocalDvTypesList(Ip-1);
//                }
//
//                // assemble list containing all used dv types. Dv types are not unique
//                MPI_Gatherv(((aPdvTypeList.data()).data()),
//                               tNumLocalDvTypes,
//                               MPI_UNSIGNED,
//                               (mIpPdvTypes.data()).data(),
//                               (tNumLocalDvTypesList.data()).data(),
//                               (tDofTypeOffset.data()).data(),
//                               MPI_UNSIGNED,
//                               0,
//                               MPI_COMM_WORLD);
//
//                // temporary variable for mIpPdvTypes size
//                moris::uint tPdvTypeListSize;
//
//                if (par_rank() == 0)
//                {
//                    // sort this created list
//                    std::sort((mIpPdvTypes.data()).data(), (mIpPdvTypes.data()).data() + mIpPdvTypes.size());
//
//                    // use std::unique and std::distance to create list containing all used dof types. This list is unique
//                    auto last = std::unique((mIpPdvTypes.data()).data(), (mIpPdvTypes.data()).data() + mIpPdvTypes.size());
//                    auto pos  = std::distance((mIpPdvTypes.data()).data(), last);
//
//                    mIpPdvTypes.resize(pos);
//
//                    tPdvTypeListSize = mIpPdvTypes.size();
//                }
//
//                // bcast size of mIpPdvTypes on processor 0
//                MPI_Bcast(& tPdvTypeListSize, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
//
//                // resize mIpPdvTypes on all processors
//                mIpPdvTypes.resize(tPdvTypeListSize);
//
//                // bcast unique mIpPdvTypes to all processors
//                MPI_Bcast((mIpPdvTypes.data()).data(), mIpPdvTypes.size(), MPI_UNSIGNED, 0, MPI_COMM_WORLD);
//            }
//        //------------------------------------------------------------------------------
//            /**
//             * communicate dv types
//             * @param aPdvTypeList a local list of dv types
//             */
//            void communicate_ig_dv_types(moris::Cell<GEN_DV>& aPdvTypeList)
//            {
//                // get processor size
//                int tSize = par_size();
//
//                // get number of local dv types
//                moris::sint tNumLocalDvTypes = aPdvTypeList.size();
//
//                // variable for maximal possible global dv types ???
//                moris::sint tNumMaxGlobalDvTypes;
//
//                // get number of global dv types
//                sum_all(tNumLocalDvTypes, tNumMaxGlobalDvTypes);
//
//                if (par_rank() == 0)
//                {
//                    // set size of dv type list = number of global types
//                    mIgPdvTypes.resize(tNumMaxGlobalDvTypes);
//                }
//
//                // create list containing the number of local dv types
//                moris::Cell <moris::sint> tNumLocalDvTypesList (tSize);
//
//                // insert number of local dv types into list containing the number of local dv types
//                MPI_Allgather(&tNumLocalDvTypes, 1, MPI_UNSIGNED, (tNumLocalDvTypesList.data()).data(), 1, MPI_UNSIGNED,  MPI_COMM_WORLD);
//
//                // create list containing the offsets of the local dv types in relation to processor 0
//                moris::Cell<moris::sint> tDofTypeOffset(tSize, 0);
//
//                // fill the list with the corresponding offsets
//                for (int Ip = 1; Ip <tSize; ++Ip)
//                {
//                    tDofTypeOffset(Ip) = tDofTypeOffset(Ip-1) + tNumLocalDvTypesList(Ip-1);
//                }
//
//                // assemble list containing all used dv types. Dv types are not unique
//                MPI_Gatherv(((aPdvTypeList.data()).data()),
//                               tNumLocalDvTypes,
//                               MPI_UNSIGNED,
//                               (mIgPdvTypes.data()).data(),
//                               (tNumLocalDvTypesList.data()).data(),
//                               (tDofTypeOffset.data()).data(),
//                               MPI_UNSIGNED,
//                               0,
//                               MPI_COMM_WORLD);
//
//                // temporary variable for mIpPdvTypes size
//                moris::uint tPdvTypeListSize;
//
//                if (par_rank() == 0)
//                {
//                    // sort this created list
//                    std::sort((mIgPdvTypes.data()).data(), (mIgPdvTypes.data()).data() + mIgPdvTypes.size());
//
//                    // use std::unique and std::distance to create list containing all used dof types. This list is unique
//                    auto last = std::unique((mIgPdvTypes.data()).data(), (mIgPdvTypes.data()).data() + mIgPdvTypes.size());
//                    auto pos  = std::distance((mIgPdvTypes.data()).data(), last);
//
//                    mIgPdvTypes.resize(pos);
//
//                    tPdvTypeListSize = mIgPdvTypes.size();
//                }
//
//                // bcast size of mIgPdvTypes on processor 0
//                MPI_Bcast(& tPdvTypeListSize, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
//
//                // resize mIgPdvTypes on all processors
//                mIgPdvTypes.resize(tPdvTypeListSize);
//
//                // bcast unique mIgPdvTypes to all processors
//                MPI_Bcast((mIgPdvTypes.data()).data(), mIgPdvTypes.size(), MPI_UNSIGNED, 0, MPI_COMM_WORLD);
//            }

            /**
             * Update local to global dv type map for interpolation nodes
             */
            void update_ip_local_to_global_dv_type_map();

            /**
             * update local to global dv type map for interpolation nodes
             */
            void update_ig_local_to_global_dv_type_map();
        };
    }   // end ge namespace
}  // end moris namepspace

#endif /* PROJECTS_GEN_SRC_GEOMENG_CL_GEN_PDV_HOST_MANAGER_HPP_ */
