#ifndef CL_GEN_PDV_Type_HOST_MANAGER_HPP_
#define CL_GEN_PDV_Type_HOST_MANAGER_HPP_

#include "cl_MSI_Design_Variable_Interface.hpp"
#include "cl_GEN_Pdv_Host.hpp"
#include "cl_GEN_Field.hpp"
#include "cl_GEN_Pdv_Enums.hpp"
#include "cl_Matrix.hpp"

namespace moris
{
    namespace ge
    {
        class Pdv_Host_Manager : public MSI::Design_Variable_Interface
        {
        private:
            // list of pdv hosts - interpolation nodes
            Cell<std::shared_ptr<Pdv_Host>> mIpPdvHosts;
            Cell<std::shared_ptr<Pdv_Host>> mIgPdvHosts;
            
            // Groups of PDV_Type types used per set
            Cell<Cell<Cell<PDV_Type>>> mIpPdvTypes;
            Cell<Cell<Cell<PDV_Type>>> mIgPdvTypes;

            // Ungrouped PDV_Type types
            Cell<Cell<PDV_Type>> mUniqueIpPdvTypes;
            Cell<Cell<PDV_Type>> mUniqueIgPdvTypes;
            
            // Requested PDV_Type types
            Cell<PDV_Type> mRequestedIpPdvTypes;
            Cell<PDV_Type> mRequestedIgPdvTypes;

            // List of global indices for identifying a given local PDV
            Matrix<IndexMat> mGlobalPdvTypeMap;

            // Requested IQI types
            Cell<std::string> mRequestedIQIs;
            
            // Pdv index
            uint mGlobalPdvIndex = 0;
            
        public:
            
            /**
             * Constructor
             */
            Pdv_Host_Manager();
            
            /**
             * Destructor
             */
            ~Pdv_Host_Manager();
            
            /**
             * Get dv types for set
             *
             * @param aIPMeshSetIndex integration mesh index
             * @param aPdvTypes        list of groups of dv types to fill
             */
            void get_ip_dv_types_for_set(const moris_index aIGMeshSetIndex, Cell<Cell<PDV_Type>>& aPdvTypes);
            
            /**
             * Get dv types for set
             *
             * @param aIGMeshSetIndex integration mesh index
             * @param aPdvTypes        list of groups of dv types to fill
             */
            void get_ig_dv_types_for_set(const moris_index aIGMeshSetIndex, Cell<Cell<PDV_Type>>& aPdvTypes);
            
            /**
             * Get unique dv types for set
             *
             * @param aIPMeshSetIndex integration mesh index
             * @param aPdvTypes        list dv types to fill
             */
            void get_ip_unique_dv_types_for_set(const moris_index aIGMeshSetIndex, Cell<PDV_Type>& aPdvTypes);
            
            /**
             * Get unique dv types for set
             *
             * @param aIGMeshSetIndex integration mesh index
             * @param aPdvTypes        list dv types to fill
             */
            void get_ig_unique_dv_types_for_set(const moris_index aIGMeshSetIndex, Cell<PDV_Type>& aPdvTypes);
            
            /**
             * Get pdv values for requested vertex indices and dv types
             *
             * @param aNodeIndices list of node indices
             * @param aPdvTypes     list of dv types
             * @param aDvValues    list of returned dv values (DvType)(vertexIndex)
             */
            void get_ip_pdv_value(const Matrix<IndexMat>&   aNodeIndices,
                                  const Cell<PDV_Type>&       aPdvTypes,
                                  Cell<Matrix<DDRMat>>&     aDvValues);
            
            /**
             * Get pdv values for requested vertex indices and dv types
             *
             * @param aNodeIndices list of node indices
             * @param aPdvTypes     list of dv types
             * @param aDvValues    list of returned dv values (DvType)(vertexIndex)
             * @param aIsActive    list of if design variable is active (vertexIndex)(DvType)
             */
            void get_ip_pdv_value(const Matrix<IndexMat>&   aNodeIndices,
                                  const Cell<PDV_Type>&       aPdvTypes,
                                  Cell<Matrix<DDRMat>>&     aDvValues,
                                  Cell<Matrix<DDSMat>>&     aIsActiveDv);
            
            /**
             * Get pdv values for requested vertex indices and dv types
             *
             * @param aNodeIndices list of node indices
             * @param aPdvTypes     list of dv types
             * @param aDvValues    list of dv values (DvType)(vertexIndex)
             */
            void get_ig_pdv_value(const Matrix<IndexMat>&   aNodeIndices,
                                  const Cell<PDV_Type>&       aPdvTypes,
                                  Cell<Matrix<DDRMat>>&     aDvValues);
            
            /**
             * Get pdv values for requested vertex indices and dv types
             *
             * @param aNodeIndices list of node indices
             * @param aPdvTypes     list of dv types
             * @param aDvValues    list of dv values (DvType)(vertexIndex)
             * @param aIsActive    list of active design variables (vertexIndex)(DvType)
             */
            void get_ig_pdv_value(const Matrix<IndexMat>&   aNodeIndices,
                                  const Cell<PDV_Type>&       aPdvTypes,
                                  Cell<Matrix<DDRMat>>&     aDvValues,
                                  Cell<Matrix<DDSMat>>&     aIsActiveDv);
            
            /**
             * Get the local to global pdv type map
             *
             * @return Matrix map from pdv type to index
             */
            const Matrix<DDSMat> & get_my_local_global_map();
            
            /**
             * Return local to global DV type map
             *
             * @param aNodeIndex   List of vertex indices
             * @param aPdvType        List of Dv types
             * @param aDvIds         List of Dv Ids
             */
            void get_ip_dv_ids_for_type_and_ind(const Matrix<IndexMat>&    aNodeIndices,
                                                const Cell<PDV_Type>&         aPdvTypes,
                                                Cell<Matrix<IdMat>>&        aDvIds);
            
            /**
             * Get local to global DV type map
             *
             * @param aNodeIndex   List of vertex indices
             * @param aPdvType        List of Dv types
             * @param aDvIds         List of Dv Ids
             */
            void get_ig_dv_ids_for_type_and_ind(const Matrix<IndexMat>&    aNodeIndices,
                                                const Cell<PDV_Type>&         aPdvTypes,
                                                Cell<Matrix<IdMat>>&        aDvIds);

            /**
             * Get requested pdv types on interpolation mesh nodes for sensitivity analysis
             *
             * @param[ in ] aPdvTypes list of dv types to fill
             */
            void get_ip_requested_dv_types( Cell< PDV_Type > & aPdvTypes );

            /**
             * Get requested pdv types on integration mesh nodes for sensitivity analysis
             *
             * @param[ in ] aPdvTypes list of dv types to fill
             */
            void get_ig_requested_dv_types( Cell< PDV_Type > & aPdvTypes );
            
            /**
             * Assign a GEN property to pdv type by node index
             * 
             * @param aPropertyPointer a GEN property pointer
             * @param aPdvType         a list of dv types
             * @param aNodeIndex     a node index
             */
            void assign_property_to_pdv_type_by_vertex_index(std::shared_ptr<GEN_Property>  aPropertyPointer,
                                                              PDV_Type                        aPdvType,
                                                              moris_index                   aNodeIndex);
            
            /**
             * Assign a GEN Field to pdv type by node index
             * 
             * @param aFieldPointer a GEN Field pointer
             * @param aPdvType      a list of dv types
             * @param aNodeIndex  a node index
             */
            void assign_field_to_pdv_type_by_vertex_index(std::shared_ptr<GEN_Field>   aFieldPointer,
                                                          PDV_Type                       aPdvType,
                                                          moris_index                  aNodeIndex);
            
            /**
             * Get pdv by type and node index
             *
             * @param aNodeIndex     a node index
             * @param aPdvType         a list of dv types
             */
            std::shared_ptr<GEN_Pdv> get_ip_pdv_by_type_and_index(moris_index aNodeIndex, PDV_Type aPdvType);

            /**
             * Get pdv by type and node index
             *
             * @param aNodeIndex     a node index
             * @param aPdvType         a list of dv types
             */
            std::shared_ptr<GEN_Pdv> get_ig_pdv_by_type_and_index(moris_index aNodeIndex, PDV_Type aPdvType);

            /**
             * Create the pdv hosts on interpolation nodes based on the pdv types per set
             *
             * @param aTotalNodes The total number of nodes where a PDV_Type host will be created
             * @param aNodeIndicesPerSet The node indices contained on a set
             * @param aPdvTypes The PDV_Type types per set, grouped
             */
            void create_ip_pdv_hosts(uint aTotalNodes, Cell<Matrix<DDSMat>> aNodeIndicesPerSet, Cell<Cell<Cell<PDV_Type>>> aPdvTypes);

            /**
             * Create the pdv hosts on integration nodes based on the pdv types per set
             *
             * @param aTotalNodes The total number of nodes where a PDV_Type host will be created
             * @param aNodeIndicesPerSet The node indices contained on a set
             * @param aPdvTypes The PDV_Type types per set, grouped
             */
            void create_ig_pdv_hosts(uint aTotalNodes, Cell<Matrix<DDSMat>> aNodeIndicesPerSet, Cell<Cell<Cell<PDV_Type>>> aPdvTypes);

            /**
             * Mark a PDV_Type on an interpolation node as being inactive
             *
             * @param aNodeIndex IP node index
             * @param aPdvType PDV_Type on the node to be marked
             */
            void mark_ip_pdv_as_inactive(moris_index aNodeIndex, PDV_Type aPdvType);

            /**
             * Mark a PDV_Type on an integration node as being inactive
             *
             * @param aNodeIndex IG node index
             * @param aPdvType PDV_Type on the node to be marked
             */
            void mark_ig_pdv_as_inactive(moris_index aNodeIndex, PDV_Type aPdvType);
            
            /**
             * Set the requested interpolation node PDV_Type types for sensitivities
             *
             * @param aPdvTypes the pdv types which will be requested by MDL
             */
            void set_ip_requested_dv_types(Cell<PDV_Type>& aPdvTypes);

            /**
             * Set the requested integration node PDV_Type types for sensitivities
             *
             * @param aPdvTypes the pdv types which will be requested by MDL
             */
            void set_ig_requested_dv_types(Cell<PDV_Type>& aPdvTypes);

            /**
             * Create PDV_Type on interpolation mesh node with GEN field
             *
             * @param aNodeIndex Node index for pulling a value from the field
             * @param aPdvType PDV_Type type
             * @param aFieldPointer Pointer to a GEN field
             */
            void create_ip_pdv(uint aNodeIndex, PDV_Type aPdvType, std::shared_ptr<GEN_Field> aFieldPointer);

            /**
             * Create PDV_Type on interpolation mesh node with GEN property
             *
             * @param aNodeIndex Node index
             * @param aPdvType PDV_Type type
             * @param aPropertyPointer Pointer to a GEN property
             */
            void create_ip_pdv(uint aNodeIndex, PDV_Type aPdvType, std::shared_ptr<GEN_Property> aPropertyPointer);

            /**
             * Create PDV_Type on interpolation mesh node with real value
             *
             * @param aNodeIndex Node index
             * @param aPdvType PDV_Type type
             * @param aPdvVal PDV_Type value
             */
            void create_ip_pdv(uint aNodeIndex, PDV_Type aPdvType, moris::real aPdvVal);

            /**
             * Create PDV_Type on integration mesh node with real value
             *
             * @param aNodeIndex Node index
             * @param aPdvType PDV_Type type
             * @param aPdvVal PDV_Type value
             */
            void create_ig_pdv(uint aNodeIndex, PDV_Type aPdvType, moris::real aPdvVal);

//            /**
//             * communicate dv types
//             * @param aPdvTypeList a local list of dv types
//             */
//            void communicate_ip_dv_types(moris::Cell<PDV_Type>& aPdvTypeList)
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
//            void communicate_ig_dv_types(moris::Cell<PDV_Type>& aPdvTypeList)
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

        };
    }   // end ge namespace
}  // end moris namepspace

#endif /* PROJECTS_GEN_SRC_GEOMENG_CL_GEN_PDV_HOST_MANAGER_HPP_ */
