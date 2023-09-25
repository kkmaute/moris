/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Pdv_Host_Manager.hpp
 *
 */

#ifndef MORIS_CL_GEN_PDV_HOST_MANAGER_HPP_
#define MORIS_CL_GEN_PDV_HOST_MANAGER_HPP_

#include "cl_MSI_Design_Variable_Interface.hpp"
#include "cl_GEN_Interpolation_Pdv_Host.hpp"
#include "cl_GEN_Intersection_Node.hpp"
#include "GEN_Data_Types.hpp"

#include "cl_Matrix.hpp"
#include "cl_SOL_Dist_Matrix.hpp"
#include <unordered_map>

namespace moris
{
    namespace ge
    {
        class Pdv_Host_Manager : public MSI::Design_Variable_Interface
        {
            //-------------------------------------------------------------------------------

          private:
            // GEN-MESH map
            moris::Cell< moris_index > mGenMeshMap;
            bool                       mGenMeshMapIsInitialized = false;

            // ADV IDs
            Matrix< DDSMat > mOwnedADVIds;
            bool             mADVIdsSet = false;

            // Static nodes set
            bool mNumBackgroundNodesSet = false;

            // PDV type map
            moris::Cell< enum PDV_Type > mPdvTypeList;    // List containing all used unique dv types.
            Matrix< DDSMat >             mPdvTypeMap;     // Map which maps the unique dv types onto consecutive values.
            Matrix< IdMat >              mCommTable;

            // list of pdv hosts - interpolation nodes
            Cell< std::shared_ptr< Interpolation_Pdv_Host > > mIpPdvHosts;
            Cell< std::shared_ptr< Intersection_Node > >      mIntersectionNodes;

            std::unordered_map< moris_id, moris_index > mIGVertexIdtoIndMap;
            std::unordered_map< moris_id, moris_index > mIPVertexIdtoIndMap;
            std::unordered_map< moris_id, moris_index > mIPBaseVertexIdtoIndMap;

            // Groups of PDV types used per set
            Cell< Cell< Cell< PDV_Type > > > mIpPdvTypes;
            Cell< Cell< Cell< PDV_Type > > > mIgPdvTypes;

            // Ungrouped PDV types
            Cell< Cell< PDV_Type > > mUniqueIpPdvTypes;
            Cell< Cell< PDV_Type > > mUniqueIgPdvTypes;

            // Requested PDV types
            Cell< PDV_Type > mRequestedIpPdvTypes;
            Cell< PDV_Type > mRequestedIgPdvTypes;

            // List of global indices for identifying a given local PDV
            Matrix< IndexMat > mOwnedPdvLocalToGlobalMap;
            Matrix< IndexMat > mOwnedAndSharedPdvLocalToGlobalMap;

            uint mNumOwnedPdvs          = 0;
            uint mNumOwnedAndSharedPdvs = 0;

            // Requested IQI types
            Cell< std::string > mRequestedIQIs;

            //-------------------------------------------------------------------------------

          public:
            //-------------------------------------------------------------------------------

            /**
             * Constructor
             */
            Pdv_Host_Manager();

            //-------------------------------------------------------------------------------

            /**
             * Destructor
             */
            ~Pdv_Host_Manager();

            //-------------------------------------------------------------------------------

            const Matrix< DDSMat >&
            get_pdv_type_map()
            {
                return mPdvTypeMap;
            }

            //-------------------------------------------------------------------------------

            uint
            get_max_num_pdvs()
            {
                return mPdvTypeList.size();
            }

            //-------------------------------------------------------------------------------

            /**
             * Sets the owned ADV IDs.
             *
             * @param aOwnedADVIds Owned ADV IDs
             */
            void set_owned_adv_ids( Matrix< DDSMat > aOwnedADVIds );

            //-------------------------------------------------------------------------------

            /**
             * Sets the number of nodes.
             *
             * @param aNumNodes Number of nodes
             */
            void
            set_num_background_nodes( uint aNumNodes );

            //-------------------------------------------------------------------------------

            void
            set_communication_table( const Matrix< IdMat >& aCommTable );

            //-------------------------------------------------------------------------------

            Matrix< IdMat > get_communication_table();

            //-------------------------------------------------------------------------------

            void
            set_vertex_global_to_local_maps(
                    std::unordered_map< moris_id, moris_index >& aIPVertexGlobaToLocalMap,
                    std::unordered_map< moris_id, moris_index >& aIGVertexGlobaToLocalMap )
            {
                mIPVertexIdtoIndMap = aIPVertexGlobaToLocalMap;
                mIGVertexIdtoIndMap = aIGVertexGlobaToLocalMap;
            };

            //-------------------------------------------------------------------------------

            void
            set_GenMeshMap( moris::Cell< moris_index > aGenMeshMap );

            //-------------------------------------------------------------------------------

            /**
             * Resets the stored information about PDV hosts.
             */
            void
            reset();

            //-------------------------------------------------------------------------------

            /**
             * Get dv types for set
             *
             * @param aIPMeshSetIndex integration mesh index
             * @param aPdvTypes       list of groups of dv types to fill
             */
            void
            get_ip_dv_types_for_set(
                    const moris_index         aIGMeshSetIndex,
                    Cell< Cell< PDV_Type > >& aPdvTypes );

            //-------------------------------------------------------------------------------

            /**
             * Get dv types for set
             *
             * @param aIGMeshSetIndex integration mesh index
             * @param aPdvTypes       list of groups of dv types to fill
             */
            void
            get_ig_dv_types_for_set(
                    const moris_index         aIGMeshSetIndex,
                    Cell< Cell< PDV_Type > >& aPdvTypes );

            //-------------------------------------------------------------------------------

            /**
             * Get unique dv types for set
             *
             * @param aIPMeshSetIndex integration mesh index
             * @param aPdvTypes       list dv types to fill
             */
            void
            get_ip_unique_dv_types_for_set(
                    const moris_index aIGMeshSetIndex,
                    Cell< PDV_Type >& aPdvTypes );

            //-------------------------------------------------------------------------------

            /**
             * Get unique dv types for set
             *
             * @param aIGMeshSetIndex integration mesh index
             * @param aPdvTypes       list dv types to fill
             */
            void
            get_ig_unique_dv_types_for_set(
                    const moris_index aIGMeshSetIndex,
                    Cell< PDV_Type >& aPdvTypes );

            //-------------------------------------------------------------------------------

            /**
             * Get pdv values for requested vertex indices and dv types
             *
             * @param aNodeIndices list of node indices
             * @param aPdvTypes    list of dv types
             * @param aDvValues    list of returned dv values (DvType)(vertexIndex)
             */
            void
            get_ip_pdv_value(
                    const Matrix< IndexMat >& aNodeIndices,
                    const Cell< PDV_Type >&   aPdvTypes,
                    Cell< Matrix< DDRMat > >& aDvValues );

            //-------------------------------------------------------------------------------

            /**
             * Get pdv values for requested vertex indices and dv types
             *
             * @param aNodeIndices list of node indices
             * @param aPdvTypes    list of dv types
             * @param aDvValues    list of returned dv values (DvType)(vertexIndex)
             * @param aIsActive    list of if design variable is active (vertexIndex)(DvType)
             */
            void
            get_ip_pdv_value(
                    const Matrix< IndexMat >& aNodeIndices,
                    const Cell< PDV_Type >&   aPdvTypes,
                    Cell< Matrix< DDRMat > >& aDvValues,
                    Cell< Matrix< DDSMat > >& aIsActiveDv );

            //-------------------------------------------------------------------------------

            /**
             * Get pdv values for requested vertex indices and dv types
             *
             * @param aNodeIndices list of node indices
             * @param aPdvTypes    list of dv types
             * @param aDvValues    list of dv values (DvType)(vertexIndex)
             */
            void
            get_ig_pdv_value(
                    const Matrix< IndexMat >& aNodeIndices,
                    const Cell< PDV_Type >&   aPdvTypes,
                    Cell< Matrix< DDRMat > >& aDvValues );

            //-------------------------------------------------------------------------------

            /**
             * Get pdv values for requested vertex indices and dv types
             *
             * @param aNodeIndices list of node indices
             * @param aPdvTypes    list of dv types
             * @param aDvValues    list of dv values (DvType)(vertexIndex)
             * @param aIsActive    list of active design variables (vertexIndex)(DvType)
             */
            void
            get_ig_pdv_value(
                    const Matrix< IndexMat >& aNodeIndices,
                    const Cell< PDV_Type >&   aPdvTypes,
                    Cell< Matrix< DDRMat > >& aDvValues,
                    Cell< Matrix< DDSMat > >& aIsActiveDv );

            //-------------------------------------------------------------------------------

            /**
             * Get the local to global pdv type map
             *
             * @return Matrix map from pdv type to index
             */
            const Matrix< DDSMat >& get_my_local_global_map();

            const Matrix< DDSMat >& get_my_local_global_overlapping_map();

            //-------------------------------------------------------------------------------

            /**
             * Return local to global DV type map
             *
             * @param aNodeIndex List of vertex indices
             * @param aPdvType   List of Dv types
             * @param aDvIds     List of Dv Ids
             */
            void get_ip_dv_ids_for_type_and_ind(
                    const Matrix< IndexMat >& aNodeIndices,
                    const Cell< PDV_Type >&   aPdvTypes,
                    Cell< Matrix< IdMat > >&  aDvIds );

            //-------------------------------------------------------------------------------

            /**
             * Get local to global DV type map
             *
             * @param aNodeIndex List of vertex indices
             * @param aPdvType   List of Dv types
             * @param aDvIds     List of Dv Ids
             */
            void get_ig_dv_ids_for_type_and_ind(
                    const Matrix< IndexMat >& aNodeIndices,
                    const Cell< PDV_Type >&   aPdvTypes,
                    Cell< Matrix< IdMat > >&  aDvIds );

            //-------------------------------------------------------------------------------

            /**
             * Get requested pdv types on interpolation mesh nodes for sensitivity analysis
             *
             * @param[ in ] aPdvTypes list of dv types to fill
             */
            void get_ip_requested_dv_types( Cell< PDV_Type >& aPdvTypes );

            //-------------------------------------------------------------------------------

            /**
             * Get requested pdv types on integration mesh nodes for sensitivity analysis
             *
             * @param[ in ] aPdvTypes list of dv types to fill
             */
            void get_ig_requested_dv_types( Cell< PDV_Type >& aPdvTypes );

            //-------------------------------------------------------------------------------

            //-------------------------------------------------------------------------------

            /**
             * Create the pdv hosts on interpolation nodes based on the pdv types per set
             *
             * @param aNodeIndicesPerSet The node indices contained on a set
             * @param aNodeCoordinates The node coordinates indexed by node
             * @param aPdvTypes The PDV types per set, grouped
             */
            void set_interpolation_pdv_types(
                    const Cell< Cell< Cell< PDV_Type > > >& aPdvTypes );

            /**
             * Create the pdv hosts on interpolation nodes based on the pdv types per set
             *
             * @param aNodeIndicesPerSet The node indices contained on a set
             * @param aNodeCoordinates The node coordinates indexed by node
             * @param aPdvTypes The PDV types per set, grouped
             */
            void create_interpolation_pdv_hosts(
                    const Cell< Cell< uint > >& aNodeIndicesPerSet,
                    const Cell< Cell< sint > >& aNodeIdsPerSet,
                    const Cell< Cell< uint > >& aNodeOwnersPerSet,
                    const Cell< Matrix< DDRMat > >& aNodeCoordinates );

            //-------------------------------------------------------------------------------

            /**
             * Remove columns of sensitivity values associated with unused (invalid) ADVs
             *
             *@param aADVIds Vector of ADV IDs
             *@param aHostADVSensitivities Matrix of sensitivity values
             */
            void
            remove_sensitivities_of_unused_variables(
                    Matrix< DDSMat >& aADVIds,
                    Matrix< DDRMat >& aHostADVSensitivities );

            //-------------------------------------------------------------------------------

            /**
             * Set the integration PDV types per set.
             *
             * @param aPdvTypes The PDV types per set, grouped
             */
            void set_integration_pdv_types( Cell< Cell< Cell< PDV_Type > > > aPdvTypes );

            //-------------------------------------------------------------------------------

            /**
             * Set an intersection at a node index and assign its starting PDV index for later.
             *
             * @param aNodeIndex Node index
             * @param aIntersectionNode Intersection node admitted by the geometry engine
             */
            void set_intersection_node(
                    uint                                 aNodeIndex,
                    std::shared_ptr< Intersection_Node > aIntersectionNode );

            //-------------------------------------------------------------------------------

            /**
             * Gets an intersection node from the PDV host manager, used so the geometry engine can create new
             * intersection nodes.
             *
             * @param aNodeIndex Node index
             * @return Intersection node already been set
             */
            std::shared_ptr< Intersection_Node > get_intersection_node( uint aNodeIndex );

            //-------------------------------------------------------------------------------

            /**
             * Updates an intersection node with a node ID and node owner.
             *
             * @param aNodeIndex Node index
             * @param aNodeId Node ID
             * @param aNodeOwner Node owner
             */
            void update_intersection_node(
                    const moris_index& aNodeIndex,
                    const moris_index& aNodeId,
                    const moris_index& aNodeOwner );

            //-------------------------------------------------------------------------------

            /**
             * Set the requested interpolation node PDV types for sensitivities
             *
             * @param aPdvTypes the pdv types which will be requested by MDL
             */
            void set_requested_interpolation_pdv_types( Cell< PDV_Type > aPdvTypes );

            //-------------------------------------------------------------------------------

            /**
             * Set the requested integration node PDV types for sensitivities
             *
             * @param aPdvTypes the pdv types which will be requested by MDL
             */
            void set_requested_integration_pdv_types( Cell< PDV_Type > aPdvTypes );

            //-------------------------------------------------------------------------------

            /**
             * Create PDV on interpolation mesh node with real value
             *
             * @param aNodeIndex Node index
             * @param aPdvType PDV type
             * @param aPdvVal PDV value
             */
            void create_interpolation_pdv(
                    uint        aNodeIndex,
                    PDV_Type    aPdvType,
                    moris::real aPdvVal );

            //-------------------------------------------------------------------------------

            /**
             * Create PDV on interpolation mesh node with GEN property
             *
             * @param aNodeIndex Node index
             * @param aPdvType PDV type
             * @param aProperty Pointer to a GEN property
             */
            void create_interpolation_pdv(
                    uint                        aNodeIndex,
                    PDV_Type                    aPdvType,
                    std::shared_ptr< Property > aProperty );

            //-------------------------------------------------------------------------------

            /**
             * Does the necessary chain rule on the IQI derivatives with respect to PDVs which each of the PDV
             * derivatives with respect to the ADVs, to obtain the complete sensitivities.
             *
             * @return Matrix of optimization sensitivities
             */
            Matrix< DDRMat > compute_diqi_dadv( const Matrix< DDSMat >& aFullADVIds );

            //-------------------------------------------------------------------------------

            void communicate_dof_types( moris::Cell< enum PDV_Type >& aPdvTypeList );

            //-------------------------------------------------------------------------------

            void create_dv_type_map();

            //-------------------------------------------------------------------------------

            void create_pdv_ids();

            //-------------------------------------------------------------------------------

            void
            set_dQIdp(
                    Cell< Matrix< DDRMat >* > adQIdp,
                    Matrix< DDSMat >*         aMap );

            //-------------------------------------------------------------------------------

          private:

            /**
             * Counts the number of owned and shared PDVs, including both interpolation PDVs and intersection nodes.
             * These values are stored in mNumOwnedPdvs and mNumOwnedAndSharedPdvs.
             */
            void count_owned_and_shared_pdvs();

            /**
             * Communicate ID offsets, for setting unique PDV IDs across processors
             *
             * @param aNumOwnedIDs Number of owned IDs on this processor
             * @return The offset for this processor
             */
            static uint communicate_offsets( uint aNumOwnedIDs );

            //-------------------------------------------------------------------------------

            void set_owned_pdv_ids( uint aPdvOffset );

            //-------------------------------------------------------------------------------

            moris_id set_owned_interpolation_pdv_ids( moris_id aOwnedIdCounter );

            //-------------------------------------------------------------------------------

            moris_id set_owned_intersection_node_pdv_ids( moris_id aOwnedIdCounter );

            //-------------------------------------------------------------------------------

            void communicate_shared_pdv_ids();

            //-------------------------------------------------------------------------------

            void communicate_shared_interpolation_pdv_ids();

            //-------------------------------------------------------------------------------

            void communicate_shared_intersection_node_pdv_ids();

            //-------------------------------------------------------------------------------

            void build_local_to_global_maps();

            //-------------------------------------------------------------------------------

        };    // class Pdv_Host_Manager

        //-------------------------------------------------------------------------------

    }    // namespace ge
}    // namespace moris

#endif /* MORIS_CL_GEN_PDV_HOST_MANAGER_HPP_ */
