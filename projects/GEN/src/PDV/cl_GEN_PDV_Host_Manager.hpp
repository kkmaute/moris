/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_PDV_Host_Manager.hpp
 *
 */

#pragma once

#include "cl_MSI_Design_Variable_Interface.hpp"
#include "cl_GEN_Interpolation_PDV_Host.hpp"
#include "cl_GEN_Intersection_Node.hpp"
#include "cl_GEN_Floating_Node.hpp"
#include "cl_GEN_Node_Manager.hpp"
#include "GEN_Data_Types.hpp"

#include "cl_Matrix.hpp"
#include "cl_SOL_Dist_Matrix.hpp"
#include <unordered_map>

namespace moris::gen
{
    class PDV_Host_Manager : public MSI::Design_Variable_Interface
            , public std::enable_shared_from_this< PDV_Host_Manager >
    {
      private:
        // Node manager
        Node_Manager& mNodeManager;

        // GEN-MESH map
        Vector< moris_index > mGenMeshMap;
        bool                  mGenMeshMapIsInitialized = false;

        // ADV IDs
        Vector< sint > mOwnedADVIds;
        bool           mADVIdsSet = false;

        // PDV type map
        Vector< PDV_Type > mPDVTypeList;    // List containing all used unique dv types.
        Matrix< DDSMat >   mPDVTypeMap;     // Map which maps the unique dv types onto consecutive values.
        Matrix< IdMat >    mCommTable;

        // list of pdv hosts - interpolation nodes
        Vector< std::shared_ptr< Interpolation_PDV_Host > > mIpPDVHosts;
        Vector< Intersection_Node* >                        mIntersectionNodes;

        std::unordered_map< moris_id, moris_index > mIGVertexIdtoIndMap;
        std::unordered_map< moris_id, moris_index > mIPVertexIdtoIndMap;
        std::unordered_map< moris_id, moris_index > mIPBaseVertexIdtoIndMap;

        // Groups of PDV types used per set
        Vector< Vector< Vector< PDV_Type > > > mIpPDVTypes;
        Vector< Vector< Vector< PDV_Type > > > mIgPDVTypes;

        // Ungrouped PDV types
        Vector< Vector< PDV_Type > > mUniqueIpPDVTypes;
        Vector< Vector< PDV_Type > > mUniqueIgPDVTypes;

        // Requested PDV types
        Vector< PDV_Type > mRequestedIpPDVTypes;
        Vector< PDV_Type > mRequestedIgPDVTypes;

        // List of global indices for identifying a given local PDV
        Matrix< IndexMat > mOwnedPDVLocalToGlobalMap;
        Matrix< IndexMat > mOwnedAndSharedPDVLocalToGlobalMap;

        uint mNumOwnedPDVs          = 0;
        uint mNumOwnedAndSharedPDVs = 0;

      public:
        /**
         * Constructor
         *
         * @param aNodeManager Node manager from the geometry engine
         */
        explicit PDV_Host_Manager( Node_Manager& aNodeManager, Vector< std::string >& aRequestedQIs );

        /**
         * Destructor
         */
        ~PDV_Host_Manager() override = default;

        //-------------------------------------------------------------------------------

        const Matrix< DDSMat >&
        get_pdv_type_map()
        {
            return mPDVTypeMap;
        }

        //-------------------------------------------------------------------------------

        uint
        get_max_num_pdvs()
        {
            return mPDVTypeList.size();
        }

        //-------------------------------------------------------------------------------

        /**
         * Sets the owned ADV IDs.
         *
         * @param aOwnedADVIds Owned ADV IDs
         */
        void set_owned_adv_ids( const Vector< sint >& aOwnedADVIds );

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

        void
        set_GenMeshMap( const Vector< moris_index >& aGenMeshMap ) override;

        /**
         * Resets the stored information about PDV hosts.
         */
        void
        reset();

        /**
         * Get dv types for set
         *
         * @param aIPMeshSetIndex integration mesh index
         * @param aPDVTypes       list of groups of dv types to fill
         */
        void
        get_ip_dv_types_for_set(
                moris_index                   aIGMeshSetIndex,
                Vector< Vector< PDV_Type > >& aPDVTypes ) const override;

        /**
         * Get dv types for set
         *
         * @param aIGMeshSetIndex integration mesh index
         * @param aPDVTypes       list of groups of dv types to fill
         */
        void
        get_ig_dv_types_for_set(
                const moris_index             aIGMeshSetIndex,
                Vector< Vector< PDV_Type > >& aPDVTypes ) const override;

        /**
         * Get unique dv types for set
         *
         * @param aIPMeshSetIndex integration mesh index
         * @param aPDVTypes       list dv types to fill
         */
        void
        get_ip_unique_dv_types_for_set(
                moris_index         aIGMeshSetIndex,
                Vector< PDV_Type >& aPDVTypes ) const override;

        /**
         * Get unique dv types for set
         *
         * @param aIGMeshSetIndex integration mesh index
         * @param aPDVTypes       list dv types to fill
         */
        void
        get_ig_unique_dv_types_for_set(
                const moris_index   aIGMeshSetIndex,
                Vector< PDV_Type >& aPDVTypes ) const override;

        /**
         * Get pdv values for requested vertex indices and dv types
         *
         * @param aNodeIndices list of node indices
         * @param aPDVTypes    list of dv types
         * @param aDvValues    list of returned dv values (DvType)(vertexIndex)
         */
        void
        get_ip_pdv_value(
                const Matrix< IndexMat >&   aNodeIndices,
                const Vector< PDV_Type >&   aPDVTypes,
                Vector< Matrix< DDRMat > >& aDvValues ) const override;

        /**
         * Get pdv values for requested vertex indices and dv types
         *
         * @param aNodeIndices list of node indices
         * @param aPDVTypes    list of dv types
         * @param aDvValues    list of dv values (DvType)(vertexIndex)
         * @param aIsActive    list of active design variables (vertexIndex)(DvType)
         */
        void
        get_ig_pdv_value(
                const Matrix< IndexMat >&   aNodeIndices,
                const Vector< PDV_Type >&   aPDVTypes,
                Vector< Matrix< DDRMat > >& aDvValues,
                Vector< Vector< bool > >&   aIsActiveDv ) const override;

        /**
         * Get the local to global pdv type map
         *
         * @return Matrix map from pdv type to index
         */
        const Matrix< DDSMat >& get_my_local_global_map() override;

        const Matrix< DDSMat >& get_my_local_global_overlapping_map();

        /**
         * Return local to global DV type map
         *
         * @param aNodeIndex List of vertex indices
         * @param aPDVType   List of Dv types
         * @param aDvIds     List of Dv Ids
         */
        void get_ip_dv_ids_for_type_and_ind(
                const Matrix< IndexMat >&  aNodeIndices,
                const Vector< PDV_Type >&  aPDVTypes,
                Vector< Matrix< IdMat > >& aDvIds ) const override;

        /**
         * Get local to global DV type map
         *
         * @param aNodeIndex List of vertex indices
         * @param aPDVType   List of Dv types
         * @param aDvIds     List of Dv Ids
         */
        void get_ig_dv_ids_for_type_and_ind(
                const Matrix< IndexMat >&  aNodeIndices,
                const Vector< PDV_Type >&  aPDVTypes,
                Vector< Matrix< IdMat > >& aDvIds ) const override;

        /**
         * Get requested pdv types on interpolation mesh nodes for sensitivity analysis
         *
         * @param[ in ] aPDVTypes list of dv types to fill
         */
        void get_ip_requested_dv_types( Vector< PDV_Type >& aPDVTypes ) const override;

        /**
         * Get requested pdv types on integration mesh nodes for sensitivity analysis
         *
         * @param[ in ] aPDVTypes list of dv types to fill
         */
        void get_ig_requested_dv_types( Vector< PDV_Type >& aPDVTypes ) override;

        /**
         * Create the pdv hosts on interpolation nodes based on the pdv types per set
         *
         * @param aNodeIndicesPerSet The node indices contained on a set
         * @param aNodeCoordinates The node coordinates indexed by node
         * @param aPDVTypes The PDV types per set, grouped
         */
        void set_interpolation_pdv_types(
                const Vector< Vector< Vector< PDV_Type > > >& aPDVTypes );

        /**
         * Create the pdv hosts on interpolation nodes based on the pdv types per set
         *
         * @param aNodeIndicesPerSet The node indices contained on a set
         * @param aNodeCoordinates The node coordinates indexed by node
         * @param aPDVTypes The PDV types per set, grouped
         */
        void create_interpolation_pdv_hosts(
                const Vector< Vector< uint > >&   aNodeIndicesPerSet,
                const Vector< Vector< sint > >&   aNodeIdsPerSet,
                const Vector< Vector< uint > >&   aNodeOwnersPerSet,
                const Vector< Matrix< DDRMat > >& aNodeCoordinates );

        /**
         * Remove columns of sensitivity values associated with unused (invalid) ADVs
         *
         *@param aADVIds Vector of ADV IDs
         *@param aHostADVSensitivities Matrix of sensitivity values
         */
        void
        remove_sensitivities_of_unused_variables(
                Vector< sint >&   aADVIds,
                Matrix< DDRMat >& aHostADVSensitivities );

        /**
         * Set the integration PDV types per set.
         *
         * @param aPDVTypes The PDV types per set, grouped
         */
        void set_integration_pdv_types( Vector< Vector< Vector< PDV_Type > > > aPDVTypes );

        /**
         * Set the requested interpolation node PDV types for sensitivities
         *
         * @param aPDVTypes the pdv types which will be requested by MDL
         */
        void set_requested_interpolation_pdv_types( const Vector< PDV_Type >& aPDVTypes );

        /**
         * Set the requested integration node PDV types for sensitivities
         *
         * @param aPDVTypes the pdv types which will be requested by MDL
         */
        void set_requested_integration_pdv_types( const Vector< PDV_Type >& aPDVTypes );

        /**
         * Create PDV on interpolation mesh node with real value
         *
         * @param aNodeIndex Node index
         * @param aPDVType PDV type
         * @param aPDVVal PDV value
         */
        void create_interpolation_pdv(
                uint        aNodeIndex,
                PDV_Type    aPDVType,
                moris::real aPDVVal );

        /**
         * Create PDV on interpolation mesh node with GEN property
         *
         * @param aNodeIndex Node index
         * @param aPDVType PDV type
         * @param aProperty Pointer to a GEN property
         */
        void create_interpolation_pdv(
                uint                               aNodeIndex,
                PDV_Type                           aPDVType,
                const std::shared_ptr< Property >& aProperty );

        /**
         * Does the necessary chain rule on the IQI derivatives with respect to PDVs which each of the PDV
         * derivatives with respect to the ADVs, to obtain the complete sensitivities.
         *
         * @return Matrix of optimization sensitivities
         */
        Matrix< DDRMat > compute_dqi_dadv( const Vector< sint >& aFullADVIds, sol::Dist_Vector* aFulldQIdADV );

        void communicate_dof_types( Vector< enum PDV_Type >& aPDVTypeList );

        //-------------------------------------------------------------------------------

        void create_dv_type_map();

        //-------------------------------------------------------------------------------

        void create_pdv_ids();

        //-------------------------------------------------------------------------------

        // void
        // set_dQIdp(
        //         const Vector< Matrix< DDRMat >* >& adQIdp,
        //         Matrix< DDSMat >*                  aMap );`

        //-------------------------------------------------------------------------------

      private:
        /**
         * Counts the number of owned and shared PDVs, including both interpolation PDVs and intersection nodes.
         * These values are stored in mNumOwnedPDVs and mNumOwnedAndSharedPDVs.
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

        void set_owned_pdv_ids( uint aPDVOffset );

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
    };
}    // namespace moris::gen
