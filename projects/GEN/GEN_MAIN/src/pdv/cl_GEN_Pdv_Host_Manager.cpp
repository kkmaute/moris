/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Pdv_Host_Manager.cpp
 *
 */

#include "cl_GEN_Pdv_Host_Manager.hpp"
#include "cl_SOL_Matrix_Vector_Factory.hpp"
#include "cl_Communication_Tools.hpp"
#include "fn_trans.hpp"
#include "fn_eye.hpp"

// detailed logging
#include "cl_Tracer.hpp"

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

        void
        Pdv_Host_Manager::set_owned_adv_ids( Matrix< DDSMat > aOwnedADVIds )
        {
            mOwnedADVIds = aOwnedADVIds;
            mADVIdsSet   = true;
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Pdv_Host_Manager::set_communication_table( const Matrix< IdMat >& aCommTable )
        {
            mCommTable = aCommTable;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix< IdMat >
        Pdv_Host_Manager::get_communication_table()    // FIXME
        {
            return mCommTable;
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Pdv_Host_Manager::set_num_background_nodes( uint aNumNodes )
        {
            mIntersectionNodes.resize( aNumNodes );
            mNumBackgroundNodesSet = true;
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Pdv_Host_Manager::reset()
        {
            mIpPdvHosts.clear();
            mIntersectionNodes.clear();
            mOwnedPdvLocalToGlobalMap.resize( 0, 0 );
            mOwnedAndSharedPdvLocalToGlobalMap.resize( 0, 0 );
            mNumOwnedPdvs          = 0;
            mNumOwnedAndSharedPdvs = 0;
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Pdv_Host_Manager::get_ip_dv_types_for_set(
                const moris::moris_index  aIPMeshSetIndex,
                Cell< Cell< PDV_Type > >& aPdvTypes )
        {
            if ( mIpPdvTypes.size() > 0 )    // FIXME
            {
                aPdvTypes = mIpPdvTypes( aIPMeshSetIndex );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Pdv_Host_Manager::get_ig_dv_types_for_set(
                const moris::moris_index  aIGMeshSetIndex,
                Cell< Cell< PDV_Type > >& aPdvTypes )
        {
            if ( mIgPdvTypes.size() > 0 )    // FIXME
            {
                aPdvTypes = mIgPdvTypes( aIGMeshSetIndex );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Pdv_Host_Manager::get_ip_unique_dv_types_for_set(
                const moris_index aIPMeshSetIndex,
                Cell< PDV_Type >& aPdvTypes )
        {
            if ( mUniqueIpPdvTypes.size() > 0 )    // FIXME
            {
                aPdvTypes = mUniqueIpPdvTypes( aIPMeshSetIndex );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Pdv_Host_Manager::get_ig_unique_dv_types_for_set(
                const moris::moris_index aIGMeshSetIndex,
                Cell< PDV_Type >&        aPdvTypes )
        {
            if ( mUniqueIgPdvTypes.size() > 0 )    // FIXME
            {
                aPdvTypes = mUniqueIgPdvTypes( aIGMeshSetIndex );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Pdv_Host_Manager::get_ip_pdv_value(
                const Matrix< IndexMat >& aNodeIndices,
                const Cell< PDV_Type >&   aPdvTypes,
                Cell< Matrix< DDRMat > >& aDvValues )
        {
            // Get the number of node indices requested
            uint tNumIndices = aNodeIndices.numel();

            // Get the number of dv types requested
            uint tNumTypes = aPdvTypes.size();

            // Cell size
            aDvValues.resize( tNumTypes );

            // loop over the node indices
            for ( uint tPdvTypeIndex = 0; tPdvTypeIndex < tNumTypes; tPdvTypeIndex++ )
            {
                // Matrix size
                aDvValues( tPdvTypeIndex ).set_size( tNumIndices, 1 );

                // loop over the requested dv types
                for ( uint tNodeIndex = 0; tNodeIndex < tNumIndices; tNodeIndex++ )
                {
                    // get node index of PDV host
                    uint tIdnx = aNodeIndices( tNodeIndex );

                    MORIS_ASSERT( mIpPdvHosts.size() > tIdnx,
                            "Pdv_Host_Manager::get_ip_pdv_value - IP PDV host does not exist at node with index. size to small %d\n",
                            tIdnx );

                    // check that PDV host exists
                    MORIS_ASSERT( mIpPdvHosts( tIdnx ),
                            "Pdv_Host_Manager::get_ip_pdv_value - IP PDV host does not exist at node with index %d\n",
                            tIdnx );

                    aDvValues( tPdvTypeIndex )( tNodeIndex ) =
                            mIpPdvHosts( tIdnx )->get_pdv_value( aPdvTypes( tPdvTypeIndex ) );
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Pdv_Host_Manager::set_GenMeshMap( moris::Cell< moris_index > aGenMeshMap )
        {
            mGenMeshMap              = aGenMeshMap;
            mGenMeshMapIsInitialized = true;
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Pdv_Host_Manager::get_ip_pdv_value(
                const Matrix< IndexMat >& aNodeIndices,
                const Cell< PDV_Type >&   aPdvTypes,
                Cell< Matrix< DDRMat > >& aDvValues,
                Cell< Matrix< DDSMat > >& aIsActiveDv )
        {
            // Get the number of node indices requested
            uint tNumIndices = aNodeIndices.length();

            // Get the number of dv types requested
            uint tNumTypes = aPdvTypes.size();

            // Cell size
            aDvValues.resize( tNumTypes );
            aIsActiveDv.resize( tNumTypes );

            // loop over the node indices
            for ( uint tPdvTypeIndex = 0; tPdvTypeIndex < tNumTypes; tPdvTypeIndex++ )
            {
                // Matrix size
                aDvValues( tPdvTypeIndex ).set_size( tNumIndices, 1 );
                aIsActiveDv( tPdvTypeIndex ).set_size( tNumIndices, 1 );

                // loop over the requested dv types
                for ( uint tNodeIndex = 0; tNodeIndex < tNumIndices; tNodeIndex++ )
                {
                    // get node index of PDV host
                    uint tIdnx = aNodeIndices( tNodeIndex );

                    // check that PDV host exists
                    MORIS_ASSERT( mIpPdvHosts( tIdnx ),
                            "Pdv_Host_Manager::get_ip_pdv_value - IP PDV host does not exist at node with index %d\n",
                            tIdnx );

                    aDvValues( tPdvTypeIndex )( tNodeIndex ) =
                            mIpPdvHosts( tIdnx )->get_pdv_value( aPdvTypes( tPdvTypeIndex ) );

                    aIsActiveDv( tPdvTypeIndex )( tNodeIndex ) =
                            mIpPdvHosts( tIdnx )->is_active_type( aPdvTypes( tPdvTypeIndex ) );
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Pdv_Host_Manager::get_ig_pdv_value(
                const Matrix< IndexMat >& aNodeIndices,
                const Cell< PDV_Type >&   aPdvTypes,
                Cell< Matrix< DDRMat > >& aDvValues )
        {
            // Get the number of node indices requested
            uint tNumIndices = aNodeIndices.length();

            // Get the number of dv types requested
            uint tNumTypes = aPdvTypes.size();

            // Cell size
            aDvValues.resize( tNumTypes );

            // loop over the node indices
            for ( uint tPdvTypeIndex = 0; tPdvTypeIndex < tNumTypes; tPdvTypeIndex++ )
            {
                // Matrix size
                aDvValues( tPdvTypeIndex ).set_size( tNumIndices, 1 );

                // loop over the requested dv types
                for ( uint tNode = 0; tNode < tNumIndices; tNode++ )
                {
                    if ( mIntersectionNodes( aNodeIndices( tNode ) ) )
                    {
                        aDvValues( tPdvTypeIndex )( tNode ) =                    //
                                mIntersectionNodes( aNodeIndices( tNode ) )->    //
                                get_coordinate_value( static_cast< uint >( aPdvTypes( tPdvTypeIndex ) ) );
                    }
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Pdv_Host_Manager::get_ig_pdv_value(
                const Matrix< IndexMat >& aNodeIndices,
                const Cell< PDV_Type >&   aPdvTypes,
                Cell< Matrix< DDRMat > >& aDvValues,
                Cell< Matrix< DDSMat > >& aIsActiveDv )
        {
            // Get the number of node indices requested
            uint tNumIndices = aNodeIndices.length();

            // Get the number of dv types requested
            uint tNumTypes = aPdvTypes.size();

            // Cell size
            aDvValues.resize( tNumTypes );
            aIsActiveDv.resize( tNumTypes );

            // loop over the node indices
            for ( uint tPdvTypeIndex = 0; tPdvTypeIndex < tNumTypes; tPdvTypeIndex++ )
            {
                // Matrix size
                aDvValues( tPdvTypeIndex ).set_size( tNumIndices, 1 );
                aIsActiveDv( tPdvTypeIndex ).set_size( tNumIndices, 1, 0 );

                // loop over the requested dv types
                for ( uint tNode = 0; tNode < tNumIndices; tNode++ )
                {
                    // get the node index within GEN, if mtk mesh has been cleaned up, translate vertex index between GEN and MTK
                    uint tGenMeshNodeIndex = aNodeIndices( tNode );
                    if ( mGenMeshMapIsInitialized )
                    {
                        tGenMeshNodeIndex = mGenMeshMap( tGenMeshNodeIndex );
                    }

                    if ( mIntersectionNodes( tGenMeshNodeIndex ) )
                    {
                        aDvValues( tPdvTypeIndex )( tNode ) =                    //
                                mIntersectionNodes( aNodeIndices( tNode ) )->    //
                                get_coordinate_value( static_cast< uint >( aPdvTypes( tPdvTypeIndex ) ) );

                        aIsActiveDv( tPdvTypeIndex )( tNode ) = 1;
                    }
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDSMat >&
        Pdv_Host_Manager::get_my_local_global_map()
        {
            return mOwnedPdvLocalToGlobalMap;
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDSMat >&
        Pdv_Host_Manager::get_my_local_global_overlapping_map()
        {
            return mOwnedAndSharedPdvLocalToGlobalMap;
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Pdv_Host_Manager::get_ip_dv_ids_for_type_and_ind(
                const Matrix< IndexMat >& aNodeIndices,
                const Cell< PDV_Type >&   aPdvTypes,
                Cell< Matrix< IdMat > >&  aDvIds )
        {
            // get the number of node indices requested
            uint tNumIndices = aNodeIndices.length();

            // get the number of dv types requested
            uint tNumTypes = aPdvTypes.size();

            // set size for list of dv values
            aDvIds.resize( tNumTypes );

            // loop over the requested dv types
            for ( uint tPdvTypeIndex = 0; tPdvTypeIndex < tNumTypes; tPdvTypeIndex++ )
            {
                // resize vector for current PDV type
                aDvIds( tPdvTypeIndex ).resize( tNumIndices, 1 );

                // loop over the node indices
                for ( uint tNode = 0; tNode < tNumIndices; tNode++ )
                {
                    // get node index
                    uint tIdnx = aNodeIndices( tNode );

                    // check that PDV host exists
                    MORIS_ASSERT( mIpPdvHosts( tIdnx ),
                            "Pdv_Host_Manager::get_ip_pdv_value - IP PDV host does not exist at node with index %d\n",
                            tIdnx );

                    aDvIds( tPdvTypeIndex )( tNode ) =    //
                            mIpPdvHosts( tIdnx )->        //
                            get_global_index_for_pdv_type( aPdvTypes( tPdvTypeIndex ) );
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Pdv_Host_Manager::get_ig_dv_ids_for_type_and_ind(
                const Matrix< IndexMat >& aNodeIndices,
                const Cell< PDV_Type >&   aPdvTypes,
                Cell< Matrix< IdMat > >&  aDvIds )
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
                aDvIds( tPdvTypeIndex ).set_size( tNumIndices, 1, -1 );

                // loop over the node indices
                for ( uint tNode = 0; tNode < tNumIndices; tNode++ )
                {
                    uint tGenMeshNodeIndex = aNodeIndices( tNode );
                    if ( mGenMeshMapIsInitialized )
                    {
                        tGenMeshNodeIndex = mGenMeshMap( aNodeIndices( tNode ) );
                    }

                    if ( mIntersectionNodes( tGenMeshNodeIndex ) )
                    {
                        aDvIds( tPdvTypeIndex )( tNode ) =
                                mIntersectionNodes( tGenMeshNodeIndex )->get_starting_pdv_id()
                                + static_cast< uint >( aPdvTypes( tPdvTypeIndex ) );
                    }
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Pdv_Host_Manager::get_ip_requested_dv_types( Cell< PDV_Type >& aPdvTypes )
        {
            aPdvTypes = mRequestedIpPdvTypes;
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Pdv_Host_Manager::get_ig_requested_dv_types( Cell< PDV_Type >& aPdvTypes )
        {
            aPdvTypes = mRequestedIgPdvTypes;
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Pdv_Host_Manager::create_interpolation_pdv_hosts(
                const Cell< Matrix< DDSMat > >&         aNodeIndicesPerSet,
                const Cell< Matrix< DDSMat > >&         aNodeIdsPerSet,
                const Cell< Matrix< DDSMat > >&         aNodeOwnersPerSet,
                const Cell< Matrix< DDRMat > >&         aNodeCoordinatesPerSet,
                const Cell< Cell< Cell< PDV_Type > > >& aPdvTypes )
        {
            // Check that number of sets is consistent
            uint tNumSets = aPdvTypes.size();

            MORIS_ERROR( tNumSets == aNodeIndicesPerSet.size(),
                    "Pdv_Host_Manager::create_interpolation_pdv_hosts - inconsistent number of sets!" );

            // Set PDV types
            mIpPdvTypes = aPdvTypes;
            mUniqueIpPdvTypes.resize( tNumSets );

            // determine maximum node index used for sizing the pdv hosts
            moris_index tMax = 0;
            for ( moris::uint iSet = 0; iSet < tNumSets; iSet++ )
            {
                if ( aNodeIndicesPerSet( iSet ).numel() > 0 )
                {
                    tMax = std::max( aNodeIndicesPerSet( iSet ).max(), tMax );
                }
            }

            // Initialize PDV hosts
            mIpPdvHosts.resize( tMax + 1, nullptr );

            // Create PDV hosts
            for ( uint tMeshSetIndex = 0; tMeshSetIndex < tNumSets; tMeshSetIndex++ )
            {
                // Get number of unique PDV types for this set
                uint tNumUniquePdvs = 0;

                for ( uint tGroupIndex = 0; tGroupIndex < mIpPdvTypes( tMeshSetIndex ).size(); tGroupIndex++ )
                {
                    tNumUniquePdvs += mIpPdvTypes( tMeshSetIndex )( tGroupIndex ).size();
                }
                mUniqueIpPdvTypes( tMeshSetIndex ).resize( tNumUniquePdvs );

                // Copy PDV types over These are the pdvs for this set
                uint tUniquePdvIndex = 0;

                for ( uint tGroupIndex = 0; tGroupIndex < mIpPdvTypes( tMeshSetIndex ).size(); tGroupIndex++ )
                {
                    for ( uint tPdvIndex = 0; tPdvIndex < mIpPdvTypes( tMeshSetIndex )( tGroupIndex ).size(); tPdvIndex++ )
                    {
                        mUniqueIpPdvTypes( tMeshSetIndex )( tUniquePdvIndex++ ) = mIpPdvTypes( tMeshSetIndex )( tGroupIndex )( tPdvIndex );
                    }
                }

                // get number of nodes in current set
                uint tNumberOfNodes = aNodeIndicesPerSet( tMeshSetIndex ).numel();

                // Create PDV hosts on interpolation nodes
                for ( uint tNodeIndexOnSet = 0; tNodeIndexOnSet < tNumberOfNodes; tNodeIndexOnSet++ )
                {
                    // Create new host or add unique PDVs
                    moris_index tNodeIndex = aNodeIndicesPerSet( tMeshSetIndex )( tNodeIndexOnSet );
                    moris_id    tNodeId    = aNodeIdsPerSet( tMeshSetIndex )( tNodeIndexOnSet );
                    moris_index tNodeOwner = aNodeOwnersPerSet( tMeshSetIndex )( tNodeIndexOnSet );

                    Matrix< DDRMat > tNodeCoordinates = aNodeCoordinatesPerSet( tMeshSetIndex ).get_row( tNodeIndexOnSet );

                    // Create PDV host unless it already exists
                    // FIXME: why is it that if it exists it already has the same PDVtypes; this needs to be checked
                    if ( mIpPdvHosts( tNodeIndex ) == nullptr )
                    {
                        mIpPdvHosts( tNodeIndex ) =
                                std::make_shared< Interpolation_Pdv_Host >(
                                        this,
                                        tNodeIndex,
                                        tNodeId,
                                        tNodeOwner,
                                        tNodeCoordinates );

                        // create map for global IDs to local Indices for non-unzipped mesh
                        mIPBaseVertexIdtoIndMap[ tNodeId ] = tNodeIndex;
                    }
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Pdv_Host_Manager::set_integration_pdv_types( Cell< Cell< Cell< PDV_Type > > > aPdvTypes )
        {
            // Check that number of sets is consistent
            uint tNumSets = aPdvTypes.size();

            // Set PDV types
            mIgPdvTypes = aPdvTypes;
            mUniqueIgPdvTypes.resize( tNumSets );

            // Unique PDV types
            for ( uint tMeshSetIndex = 0; tMeshSetIndex < tNumSets; tMeshSetIndex++ )
            {
                // Get number of unique PDVs
                uint tNumUniquePdvs = 0;
                for ( uint tGroupIndex = 0; tGroupIndex < mIgPdvTypes( tMeshSetIndex ).size(); tGroupIndex++ )
                {
                    tNumUniquePdvs += mIgPdvTypes( tMeshSetIndex )( tGroupIndex ).size();
                }
                mUniqueIgPdvTypes( tMeshSetIndex ).resize( tNumUniquePdvs );

                // Copy PDV types over
                uint tUniquePdvIndex = 0;
                for ( uint tGroupIndex = 0; tGroupIndex < mIgPdvTypes( tMeshSetIndex ).size(); tGroupIndex++ )
                {
                    for ( uint tPdvIndex = 0; tPdvIndex < mIgPdvTypes( tMeshSetIndex )( tGroupIndex ).size(); tPdvIndex++ )
                    {
                        mUniqueIgPdvTypes( tMeshSetIndex )( tUniquePdvIndex++ ) =
                                mIgPdvTypes( tMeshSetIndex )( tGroupIndex )( tPdvIndex );
                    }
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Pdv_Host_Manager::set_intersection_node(
                uint                                 aNodeIndex,
                std::shared_ptr< Intersection_Node > aIntersectionNode )
        {
            // Check node index
            MORIS_ASSERT( mNumBackgroundNodesSet,
                    "Pdv_Host_Manager::set_intersection_node - %s",
                    "Number of background nodes must be set before intersection nodes can be created." );

            MORIS_ASSERT( aNodeIndex == mIntersectionNodes.size(),
                    "Pdv_Host_Manager::set_intersection_node - %s",
                    "Intersection nodes must be added to the PDV Host Manager in order by node index." );

            // Add intersection node
            mIntersectionNodes.push_back( aIntersectionNode );
        }

        //--------------------------------------------------------------------------------------------------------------

        std::shared_ptr< Intersection_Node >
        Pdv_Host_Manager::get_intersection_node( uint aNodeIndex )
        {
            MORIS_ASSERT( aNodeIndex < mIntersectionNodes.size(), "Node does not exist in PDV host manager." );
            return mIntersectionNodes( aNodeIndex );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Pdv_Host_Manager::update_intersection_node(
                const moris_index& aNodeIndex,
                const moris_index& aNodeId,
                const moris_index& aNodeOwner )
        {
            // MORIS_ASSERT( mIntersectionNodes(aNodeIndex)!= nullptr,
            //         "Pdv_Host_Manager::update_intersection_node(), Intersection node doe not exist.");
            // FIXME the size of this cell should be correct. this is a hack to account for a wrong size
            if ( aNodeIndex >= (sint)mIntersectionNodes.size() )
            {
                mIntersectionNodes.resize( aNodeIndex + 1, nullptr );
            }
            if ( mIntersectionNodes( aNodeIndex ) != nullptr )
            {
                mIntersectionNodes( aNodeIndex )->set_id( aNodeId );
                mIntersectionNodes( aNodeIndex )->set_owner( aNodeOwner );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Pdv_Host_Manager::set_requested_interpolation_pdv_types( Cell< PDV_Type > aPdvTypes )
        {
            mRequestedIpPdvTypes = aPdvTypes;
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Pdv_Host_Manager::set_requested_integration_pdv_types( Cell< PDV_Type > aPdvTypes )
        {
            mRequestedIgPdvTypes = aPdvTypes;
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Pdv_Host_Manager::create_interpolation_pdv(
                uint     aNodeIndex,
                PDV_Type aPdvType,
                real     aPdvVal )
        {
            // Check that PDF host exists
            MORIS_ASSERT( mIpPdvHosts( aNodeIndex ),
                    "Pdv_Host_Manager::create_interpolation_pdv - IP PDV host does not exist at node with index %d\n",
                    aNodeIndex );

            // Create PDV with given type and value
            mIpPdvHosts( aNodeIndex )->create_pdv( aPdvType, aPdvVal );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Pdv_Host_Manager::create_interpolation_pdv(
                uint                        aNodeIndex,
                PDV_Type                    aPdvType,
                std::shared_ptr< Property > aProperty )
        {
            // Check that PDV host exists
            MORIS_ASSERT( mIpPdvHosts( aNodeIndex ),
                    "Pdv_Host_Manager::create_interpolation_pdv - IP PDV host does not exist at node with index %d\n",
                    aNodeIndex );

            // Create PDV with given type and property
            mIpPdvHosts( aNodeIndex )->create_pdv( aPdvType, aProperty );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Pdv_Host_Manager::remove_sensitivities_of_unused_variables(
                Matrix< DDSMat >& aADVIds,
                Matrix< DDRMat >& aHostADVSensitivities )
        {
            // get number of ADV entries and number of sensitivities per IDs
            uint tNumAdvs = aHostADVSensitivities.n_cols();
            uint tNumSens = aHostADVSensitivities.n_rows();

            // check for consistent lengths of ID vector and matrix of sensitivities
            MORIS_ERROR( aADVIds.n_cols() == tNumAdvs,
                    "Pdv_Host_Manager::compute_diqi_dadv - inconsistent number of ADVs (%d) vs. sensitivity values for interpolation PDVs (%zu)",
                    tNumAdvs,
                    aADVIds.n_cols() );

            // skip of ADV vector is empty
            if ( tNumAdvs == 0 ) return;

            // skip if all ADV IDs are valid, i.e. positive
            if ( aADVIds.min() > -1 ) return;

            // initialize counter of used variables
            uint tCounter = 0;

            // loop over all ADV IDs
            for ( uint iv = 0; iv < aADVIds.numel(); ++iv )
            {
                // check for valid ID
                if ( aADVIds( iv ) >= 0 )
                {
                    // copy ADV IDs and sensitivity values
                    aADVIds( tCounter ) = aADVIds( iv );

                    // copy sensitivity values
                    aHostADVSensitivities.get_column( tCounter ) = aHostADVSensitivities.get_column( iv );

                    // increment counter
                    tCounter++;
                }
            }

            // resize ID vector and matrix of sensitivities
            aADVIds.resize( 1, tCounter );
            aHostADVSensitivities.resize( tNumSens, tCounter );
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix< DDRMat >
        Pdv_Host_Manager::compute_diqi_dadv( const Matrix< DDSMat >& aFullADVIds )
        {
            Tracer tTracer( "GEN", "PDV Host Manager", "compute dQi/dadv" );

            // Check for ADV IDs
            MORIS_ERROR( mADVIdsSet,
                    "PDV Host Manager must have ADV IDs set before computing sensitivities." );

            // Get dIQI/dPDV and dPDV/dADV
            sol::Dist_Vector* tdIQIdPDV = this->get_dQIdp();

            // Create factory for resulting distributed vector
            sol::Matrix_Vector_Factory tDistributedFactory;

            // Create maps
            sol::Dist_Map* tOwnedADVMap = tDistributedFactory.create_map( mOwnedADVIds );
            sol::Dist_Map* tFullADVMap  = tDistributedFactory.create_map( aFullADVIds );

            // Create distributed vectors for derivatives of IQIs wrt Advs
            sint tNumIQIs = tdIQIdPDV->get_num_vectors();

            sol::Dist_Vector* tdIQIdADV     = tDistributedFactory.create_vector( tOwnedADVMap, tNumIQIs, false, true );
            sol::Dist_Vector* tFulldIQIdADV = tDistributedFactory.create_vector( tFullADVMap, tNumIQIs, false, true );

            // Initialize derivatives of IQIs wrt Advs to zero
            tdIQIdADV->vec_put_scalar( 0.0 );

            // Loop of interpolation PDV hosts
            for ( uint tPDVHostIndex = 0; tPDVHostIndex < mIpPdvHosts.size(); tPDVHostIndex++ )
            {
                // Check if PDV host exists
                if ( mIpPdvHosts( tPDVHostIndex ) )
                {
                    // Check if processor own PDV host
                    if ( mIpPdvHosts( tPDVHostIndex )->get_pdv_owning_processor() == par_rank() )
                    {
                        // Get number of PDVs
                        uint tNumPDVsOnHost = mIpPdvHosts( tPDVHostIndex )->get_num_pdvs();

                        // Assemble sensitivities
                        for ( uint tPDVIndex = 0; tPDVIndex < tNumPDVsOnHost; tPDVIndex++ )
                        {
                            // Get PDVs
                            moris_id tPDVID = mIpPdvHosts( tPDVHostIndex )->get_pdv_id( tPDVIndex );

                            // FIXME checking if the pdv is defined
                            if ( tPDVID != -1 )
                            {
                                // Get sensitivities
                                Matrix< DDRMat > tHostADVSensitivities =
                                        mIpPdvHosts( tPDVHostIndex )->get_sensitivities( tPDVIndex );

                                // Get ADV IDs
                                Matrix< DDSMat > tADVIds =
                                        mIpPdvHosts( tPDVHostIndex )->get_determining_adv_ids( tPDVIndex );

                                // remove sensitivities wrt unused variables
                                this->remove_sensitivities_of_unused_variables( tADVIds, tHostADVSensitivities );

                                // loop over all IQIs
                                for ( uint tVectorIndex = 0; tVectorIndex < (uint)tNumIQIs; tVectorIndex++ )
                                {
                                    Matrix< DDRMat > tIndividualSensitivity =
                                            ( *tdIQIdPDV )( tPDVID, tVectorIndex ) * tHostADVSensitivities;

                                    // Fill matrix
                                    tdIQIdADV->sum_into_global_values( tADVIds, tIndividualSensitivity, tVectorIndex );
                                }
                            }
                        }
                    }
                }
            }

            // Create ADV Host Sensitivities
            Matrix< DDRMat > tHostADVSensitivities;
            Matrix< DDRMat > tI;
            // Loop over intersection nodes for inserting
            for ( uint tIntersectionIndex = 0; tIntersectionIndex < mIntersectionNodes.size(); tIntersectionIndex++ )
            {
                if ( mIntersectionNodes( tIntersectionIndex ) and mIntersectionNodes( tIntersectionIndex )->get_owner() == par_rank() )
                {
                    // Get starting ID and number of coordinates
                    uint tStartingGlobalIndex = mIntersectionNodes( tIntersectionIndex )->get_starting_pdv_id();
                    uint tNumCoordinates      = mIntersectionNodes( tIntersectionIndex )->get_num_pdvs();

                    // Parent sensitivities and ADV IDs
                    tHostADVSensitivities.set_size( 0.0, 0.0 );
                    eye( tNumCoordinates, tNumCoordinates, tI );
                    mIntersectionNodes( tIntersectionIndex )->get_dcoordinate_dadv( tHostADVSensitivities, tI );

                    Matrix< DDSMat > tADVIds =
                            mIntersectionNodes( tIntersectionIndex )->get_coordinate_determining_adv_ids();

                    // remove sensitivities wrt unused variables
                    this->remove_sensitivities_of_unused_variables( tADVIds, tHostADVSensitivities );

                    // loop overall coordinate directions
                    for ( uint tCoordinateIndex = 0; tCoordinateIndex < tNumCoordinates; tCoordinateIndex++ )
                    {
                        // get PDV ID
                        moris_id tPDVID = tStartingGlobalIndex + tCoordinateIndex;

                        // loop over all IQIs
                        for ( uint tVectorIndex = 0; tVectorIndex < (uint)tNumIQIs; tVectorIndex++ )
                        {
                            Matrix< DDRMat > tIndividualSensitivity =
                                    ( *tdIQIdPDV )( tPDVID, tVectorIndex ) *    //
                                    tHostADVSensitivities.get_row( tCoordinateIndex );

                            // Fill matrix
                            tdIQIdADV->sum_into_global_values( tADVIds, tIndividualSensitivity, tVectorIndex );
                        }
                    }
                }
            }

            // Global assembly
            tdIQIdADV->vector_global_assembly();

            // Import
            tFulldIQIdADV->import_local_to_global( *tdIQIdADV );

            // Extract values
            Matrix< DDRMat > tFullSensitivity( 0, 0 );
            tFulldIQIdADV->extract_copy( tFullSensitivity );
            tFullSensitivity = trans( tFullSensitivity );

            // Clean up
            delete tdIQIdPDV;
            delete tdIQIdADV;
            delete tFulldIQIdADV;

            return tFullSensitivity;
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Pdv_Host_Manager::communicate_dof_types( moris::Cell< enum PDV_Type >& aPdvTypeList )
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
            moris::Cell< moris::sint > tNumLocalDvTypesList( tSize );

            // Insert number of local dof types into list containing the number of local dof types
            MPI_Allgather(
                    &tNumLocalDvTypes,
                    1,
                    MPI_UNSIGNED,
                    ( tNumLocalDvTypesList.data() ).data(),
                    1,
                    MPI_UNSIGNED,
                    MPI_COMM_WORLD );

            // Create list containing the offsets of the local dof types in relation to processor 0
            moris::Cell< moris::sint > tDvTypeOffset( tSize, 0 );

            // Fill the list with the corresponding offsets
            for ( int Ip = 1; Ip < tSize; ++Ip )
            {
                tDvTypeOffset( Ip ) = tDvTypeOffset( Ip - 1 ) + tNumLocalDvTypesList( Ip - 1 );
            }

            // Assemble list containing all used dof types. Dof types are not unique
            MPI_Gatherv(
                    ( ( aPdvTypeList.data() ).data() ),
                    tNumLocalDvTypes,
                    MPI_UNSIGNED,
                    ( mPdvTypeList.data() ).data(),
                    ( tNumLocalDvTypesList.data() ).data(),
                    ( tDvTypeOffset.data() ).data(),
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
            MPI_Bcast( &tPdvTypeListSize, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD );

            // Resize mPdvTypeList on all processors
            mPdvTypeList.resize( tPdvTypeListSize );

            // Bcast unique mPdvTypeList to all processors
            MPI_Bcast( ( mPdvTypeList.data() ).data(), mPdvTypeList.size(), MPI_UNSIGNED, 0, MPI_COMM_WORLD );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Pdv_Host_Manager::create_dv_type_map()
        {
            // Get number of unique adofs of this equation object
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

        void
        Pdv_Host_Manager::communicate_check_if_owned_pdv_exists()
        {
            // Build communication table map to determine the right position for each processor rank.
            Matrix< DDSMat > tCommTableMap( mCommTable.max() + 1, 1, -1 );

            moris::uint tNumCommProcs = mCommTable.numel();

            // Loop over communication table to fill the communication table map
            for ( moris::uint Ik = 0; Ik < tNumCommProcs; Ik++ )
            {
                tCommTableMap( mCommTable( Ik ), 0 ) = Ik;
            }

            // FIXME: cannot have communication within the following loop
            // Loop over all different pdv types for IP node pdvs
            for ( moris::uint Ij = 0; Ij < mPdvTypeList.size(); Ij++ )
            {
                // define current PDV type
                enum PDV_Type tPdvType = mPdvTypeList( Ij );

                // Define vector to store number of shared PDVs
                Matrix< DDUMat > tNumSharedPdvsPerProc( tNumCommProcs, 1, 0 );

                // Loop over pdvs per type. Count number of pdvs per proc which have to be communicated
                for ( moris::uint Ib = 0; Ib < mIpPdvHosts.size(); Ib++ )
                {
                    // Check if PDV host exists
                    if ( mIpPdvHosts( Ib ) )
                    {
                        // Check that PDV of given type exists
                        if ( mIpPdvHosts( Ib )->get_pdv_exists( tPdvType ) )
                        {
                            // get owning processor
                            moris::moris_index tProcIndex = mIpPdvHosts( Ib )->get_pdv_owning_processor();

                            // Check if owning processor is not this processor
                            if ( tProcIndex != par_rank() )
                            {
                                // get position of communicating processor in communication table
                                moris::sint tProcIdPos = tCommTableMap( tProcIndex, 0 );

                                // check that processor exists in communication table
                                MORIS_ASSERT( tProcIdPos != -1,
                                        "Pdv_Host_Manager::communicate_check_if_owned_pdv_exists: Processor does "
                                        "not exist in communication table" );

                                // Add +1 to the processor number of shared dofs per processor
                                tNumSharedPdvsPerProc( tProcIdPos )++;
                            }
                        }
                    }
                }

                // Define for each communication processor a vector for communicating shared PDV Ids
                moris::Cell< Matrix< DDUMat > > tSharedPdvPosGlobal( tNumCommProcs );

                // Set size of vector for communicating shared PDV Ids
                for ( moris::uint Ik = 0; Ik < tNumCommProcs; Ik++ )
                {
                    // Get number of pdvs shared with current processor
                    uint tNumberOfSharedPDVs = tNumSharedPdvsPerProc( Ik, 0 );

                    if ( tNumberOfSharedPDVs != 0 )
                    {
                        tSharedPdvPosGlobal( Ik ).set_size( tNumberOfSharedPDVs, 1 );
                    }
                }

                // Reset vector to store number of shared pdv per processor
                tNumSharedPdvsPerProc.fill( 0 );

                // Loop over pdv per type
                for ( moris::uint Ia = 0; Ia < mIpPdvHosts.size(); Ia++ )
                {
                    // Check if PDV host exists
                    if ( mIpPdvHosts( Ia ) )
                    {
                        // Check if PDV exists
                        if ( mIpPdvHosts( Ia )->get_pdv_exists( tPdvType ) )
                        {
                            // Get owning processor rank
                            moris::moris_index tProcIndex = mIpPdvHosts( Ia )->get_pdv_owning_processor();

                            // Check if owning processor is this processor
                            if ( tProcIndex != par_rank() )
                            {
                                // get position of owning process in communication table
                                moris::sint tProcIdPos = tCommTableMap( tProcIndex, 0 );

                                // get position of next element in node ID list of owning processor
                                uint tProcListPos = tNumSharedPdvsPerProc( tProcIdPos );

                                // Add owning processor id to moris::Mat
                                tSharedPdvPosGlobal( tProcIdPos )( tProcListPos ) = mIpPdvHosts( Ia )->get_id();

                                // Add +1 for position of next element in node ID list of owning processor
                                tNumSharedPdvsPerProc( tProcIdPos )++;
                            }
                        }
                    }
                }

                // receiving list
                moris::Cell< Matrix< DDUMat > > tMatsToReceive;

                // synchronize parallel process
                // FIXME: should not be needed if done correctly
                barrier();

                // Communicate position of shared pdv to the owning processor
                communicate_mats(
                        mCommTable,
                        tSharedPdvPosGlobal,
                        tMatsToReceive );

                // check that number of received vectors is consistent with communication table
                MORIS_ERROR( tMatsToReceive.size() == tNumCommProcs,
                        "Pdv_Host_Manager::communicate_check_if_owned_pdv_exists - incorrect data received\n" );

                // Loop over all communicating processors
                for ( moris::uint Ik = 0; Ik < tNumCommProcs; Ik++ )
                {
                    // Loop over all PDV hosts which are shared
                    for ( moris::uint Ii = 0; Ii < tMatsToReceive( Ik ).numel(); Ii++ )
                    {
                        // requested PDV host id
                        uint tReqPdvHostId = tMatsToReceive( Ik )( Ii );

                        // Get index of PDV host on this processor
                        auto tIter = mIPBaseVertexIdtoIndMap.find( tReqPdvHostId );

                        // Check that host id exists
                        MORIS_ERROR( tIter != mIPBaseVertexIdtoIndMap.end(),
                                "Pdv_Host_Manager::communicate_check_if_owned_pdv_exists - PDV host with ID %d does not exist on Proc %d.\n",
                                tReqPdvHostId,
                                par_rank() );

                        // Get local index of PDV host
                        moris::uint tLocalPdvInd = tIter->second;

                        // Check that PDV host exists
                        MORIS_ERROR( mIpPdvHosts( tLocalPdvInd ),
                                "Pdv_Host_Manager::communicate_check_if_owned_pdv_exists - PDV host with ID %d and local index %d does not exist on Proc %d.\n",
                                tReqPdvHostId,
                                tLocalPdvInd,
                                par_rank() );

                        // Check if PDF of given type exists
                        MORIS_ERROR( mIpPdvHosts( tLocalPdvInd )->get_pdv_exists( tPdvType ),
                                "Pdv_Host_Manager::communicate_check_if_owned_pdv_exists - PDV missing on Node with ID %d on Proc %d.\n",
                                tReqPdvHostId,
                                par_rank() );
                    }
                }
            }

            // synchronize parallel process
            barrier();
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Pdv_Host_Manager::get_num_pdvs()
        {
            // Loop over all different pdv types for IP node pdvs
            for ( moris::uint Ij = 0; Ij < mPdvTypeList.size(); Ij++ )
            {
                enum PDV_Type tPdvType = mPdvTypeList( Ij );

                for ( moris::uint Ib = 0; Ib < mIpPdvHosts.size(); Ib++ )
                {
                    // Check if PDV host exists
                    if ( mIpPdvHosts( Ib ) )
                    {
                        // Check if PDV exists
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
            }

            // Loop over intersection node pdvs
            for ( moris::uint Ij = 0; Ij < mIntersectionNodes.size(); Ij++ )
            // for ( moris::uint Ij = 0; Ij < mGenMeshMap.size(); Ij++ )
            {
                uint tMeshNodeIndex = Ij;
                if ( mGenMeshMapIsInitialized )
                {
                    tMeshNodeIndex = mGenMeshMap( Ij );
                }

                if ( mIntersectionNodes( tMeshNodeIndex ) != nullptr )
                {
                    uint tNumPdvsOnIntersectionNode = mIntersectionNodes( tMeshNodeIndex )->get_num_pdvs();

                    if ( mIntersectionNodes( tMeshNodeIndex )->get_owner() == par_rank() )
                    {
                        mNumOwnedPdvs += tNumPdvsOnIntersectionNode;
                    }
                    mNumOwnedAndSharedPdvs += tNumPdvsOnIntersectionNode;
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        uint
        Pdv_Host_Manager::communicate_pdv_offsets( const moris::uint& aNumOwnedPdvs )
        {
            // Get list containing the number of owned pdvs of each processor
            Matrix< DDUMat > tNumOwnedPdvsList;

            comm_gather_and_broadcast( aNumOwnedPdvs, tNumOwnedPdvsList );

            Matrix< DDUMat > tOwnedPdvsOffsetList( tNumOwnedPdvsList.numel(), 1, 0 );

            // Loop over all entries to create the offsets. Starting with 1
            for ( moris::uint Ij = 1; Ij < tOwnedPdvsOffsetList.numel(); Ij++ )
            {
                // Add the number of owned pdvs of the previous processor to the offset of the previous processor
                tOwnedPdvsOffsetList( Ij, 0 ) = tOwnedPdvsOffsetList( Ij - 1, 0 ) + tNumOwnedPdvsList( Ij - 1, 0 );
            }

            return tOwnedPdvsOffsetList( par_rank(), 0 );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Pdv_Host_Manager::set_owned_pdv_ids( uint aPdvOffset )
        {
            moris::uint tOwnedIdCounter = aPdvOffset;

            tOwnedIdCounter = this->set_owned_interpolation_pdv_ids( tOwnedIdCounter );

            this->set_owned_intersection_node_pdv_ids( tOwnedIdCounter );
        }

        //--------------------------------------------------------------------------------------------------------------

        moris_id
        Pdv_Host_Manager::set_owned_interpolation_pdv_ids( moris_id aOwnedIdCounter )
        {
            moris_id tSaveOffset = aOwnedIdCounter;

            // Loop over all different pdv types for IP node pdvs
            for ( moris::uint Ij = 0; Ij < mPdvTypeList.size(); Ij++ )
            {
                enum PDV_Type tPdvType = mPdvTypeList( Ij );

                // Loop over pdvs per type.
                for ( moris::uint Ib = 0; Ib < mIpPdvHosts.size(); Ib++ )
                {
                    // Check if PDV host exists
                    if ( mIpPdvHosts( Ib ) )
                    {
                        // Check if PDV exists
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
            }

            moris_id tNumOwnedInterpolationIds = aOwnedIdCounter - tSaveOffset;

            MORIS_LOG_INFO( "System has a total of %-5i interpolation pdvs.", sum_all( tNumOwnedInterpolationIds ) );

            return aOwnedIdCounter;
        }

        //--------------------------------------------------------------------------------------------------------------

        moris_id
        Pdv_Host_Manager::set_owned_intersection_node_pdv_ids( moris_id aOwnedIdCounter )
        {
            moris_id tSaveOffset = aOwnedIdCounter;

            // Loop over intersection node pdvs
            for ( moris::uint Ij = 0; Ij < mIntersectionNodes.size(); Ij++ )
            // for ( moris::uint Ij = 0; Ij < mGenMeshMap.size(); Ij++ )
            {
                uint tMeshNodeIndex = Ij;
                if ( mGenMeshMapIsInitialized )
                {
                    tMeshNodeIndex = mGenMeshMap( Ij );
                }

                if ( mIntersectionNodes( tMeshNodeIndex ) != nullptr )
                {
                    if ( mIntersectionNodes( tMeshNodeIndex )->get_owner() == par_rank() )
                    {
                        mIntersectionNodes( tMeshNodeIndex )->set_starting_pdv_id( aOwnedIdCounter );
                        uint tNumPdvsOnIntersectionNode = mIntersectionNodes( tMeshNodeIndex )->get_num_pdvs();
                        aOwnedIdCounter += tNumPdvsOnIntersectionNode;
                    }
                }
            }

            moris_id tNumOwnedInterpolationIds = aOwnedIdCounter - tSaveOffset;

            MORIS_LOG_INFO( "System has a total of %-5i intersection node pdvs.", sum_all( tNumOwnedInterpolationIds ) );

            return aOwnedIdCounter;
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Pdv_Host_Manager::communicate_shared_pdv_ids()
        {
            this->communicate_shared_interpolation_pdv_ids();
            this->communicate_shared_intersection_node_pdv_ids();
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Pdv_Host_Manager::communicate_shared_interpolation_pdv_ids()
        {
            // Build communication table map to determine the right position for each processor rank.
            Matrix< DDSMat > tCommTableMap( mCommTable.max() + 1, 1, -1 );

            moris::uint tNumCommProcs = mCommTable.numel();

            // Loop over communication table to fill the communication table map
            for ( moris::uint Ik = 0; Ik < tNumCommProcs; Ik++ )
            {
                tCommTableMap( mCommTable( Ik ), 0 ) = Ik;
            }

            // FIXME: cannot have communication within following loop
            // Loop over all different pdv types for IP node pdvs
            for ( moris::uint Ij = 0; Ij < mPdvTypeList.size(); Ij++ )
            {
                enum PDV_Type tPdvType = mPdvTypeList( Ij );

                // Define vector to store number of shared pdv per processor
                Matrix< DDUMat > tNumSharedPdvsPerProc( tNumCommProcs, 1, 0 );

                // Loop over pdvs per type. Count number of pdvs per processor which have to be communicated
                for ( moris::uint Ib = 0; Ib < mIpPdvHosts.size(); Ib++ )
                {
                    // Check if PDV host exists
                    if ( mIpPdvHosts( Ib ) )
                    {
                        // Check if PDV exists for given type
                        if ( mIpPdvHosts( Ib )->get_pdv_exists( tPdvType ) )
                        {
                            // Get owning processor rank
                            moris::moris_index tProcIndex = mIpPdvHosts( Ib )->get_pdv_owning_processor();

                            // Check if owning processor is not this processor
                            if ( tProcIndex != par_rank() )
                            {
                                // get position of owning process in communication table
                                moris::sint tProcIdPos = tCommTableMap( tProcIndex );

                                MORIS_ASSERT( tProcIdPos != -1,
                                        "Pdv_Host_Manager::communicate_shared_interpolation_pdv_ids: Processor "
                                        "does not exist in communication table" );

                                // Add +1 to the processor number of shared dv per processor
                                tNumSharedPdvsPerProc( tProcIdPos )++;
                            }
                        }
                    }
                }

                // Define cells of vectors to store PDV host IDs and indices of shared PDvs
                moris::Cell< Matrix< DDUMat > > tSharedPdvPosGlobal( tNumCommProcs );
                moris::Cell< Matrix< DDUMat > > tSharedPdvPosLocal( tNumCommProcs );

                // Set size of vectors to store PDV host IDs and indices of shared PDvs
                for ( moris::uint Ik = 0; Ik < tNumCommProcs; Ik++ )
                {
                    // Get number of pdvs shared with current processor
                    uint tNumberOfSharedPDVs = tNumSharedPdvsPerProc( Ik, 0 );

                    // if there are any PDVs shared with this processor set size of vectors
                    if ( tNumberOfSharedPDVs != 0 )
                    {
                        tSharedPdvPosGlobal( Ik ).set_size( tNumberOfSharedPDVs, 1 );
                        tSharedPdvPosLocal( Ik ).set_size( tNumberOfSharedPDVs, 1 );
                    }
                }

                // Reset vector to store number of shared pdv per processor
                tNumSharedPdvsPerProc.fill( 0 );

                // Loop over pdvs
                for ( moris::uint Ib = 0; Ib < mIpPdvHosts.size(); Ib++ )
                {
                    // Check if PDV host exists
                    if ( mIpPdvHosts( Ib ) )
                    {
                        // Check if PDV exists for given type
                        if ( mIpPdvHosts( Ib )->get_pdv_exists( tPdvType ) )
                        {
                            // Get owning processor rank
                            moris::moris_index tProcIndex = mIpPdvHosts( Ib )->get_pdv_owning_processor();

                            // Check if owning processor is not this processor
                            if ( tProcIndex != par_rank() )
                            {
                                // get position of owning process in communication table
                                moris::sint tProcIdPos = tCommTableMap( tProcIndex );

                                // get position of next element in node ID list of owning processor
                                uint tProcListPos = tNumSharedPdvsPerProc( tProcIdPos );

                                // Add Id of PDV host to global list of owning processor
                                tSharedPdvPosGlobal( tProcIdPos )( tProcListPos ) = mIpPdvHosts( Ib )->get_id();

                                // Add local position of existing pdv hosts to local list of owning processor
                                tSharedPdvPosLocal( tProcIdPos )( tProcListPos ) = Ib;

                                // Add +1 for position of next element in node ID list of owning processor
                                tNumSharedPdvsPerProc( tProcIdPos )++;
                            }
                        }
                    }
                }

                // receiving list
                moris::Cell< Matrix< DDUMat > > tMatsToReceive;

                // FIXME: should not be needed if done correctly
                barrier();

                // Communicate position of shared pdvs to the owning processor
                communicate_mats(
                        mCommTable,
                        tSharedPdvPosGlobal,
                        tMatsToReceive );

                // check that number of received vectors is consistent with communication table
                MORIS_ERROR( tMatsToReceive.size() == tNumCommProcs,
                        "Pdv_Host_Manager::communicate_shared_interpolation_pdv_ids - incorrect data received\n" );

                // Create List of Mats containing the shared PDV host Ids
                moris::Cell< Matrix< DDUMat > > tSharedPdvIdList( tNumCommProcs );

                // Loop over all communicating processors, set size of vectors and initialize with MORIS_UINT_MAX
                for ( moris::uint Ik = 0; Ik < tNumCommProcs; Ik++ )
                {
                    tSharedPdvIdList( Ik ).set_size( tMatsToReceive( Ik ).numel(), 1, MORIS_UINT_MAX );
                }

                // Loop over all processors with PDV hosts owned by this processor
                for ( moris::uint Ik = 0; Ik < tNumCommProcs; Ik++ )
                {
                    // Loop over all PDV hosts for which IDs are requested
                    for ( moris::uint Ii = 0; Ii < tMatsToReceive( Ik ).numel(); Ii++ )
                    {
                        // requested PDV host id
                        uint tReqPdvHostId = tMatsToReceive( Ik )( Ii );

                        // Get index of PDV host on this processor
                        auto tIter = mIPBaseVertexIdtoIndMap.find( tReqPdvHostId );

                        moris::uint tLocalPdvInd = tIter->second;

                        // Check that PDV host exists
                        MORIS_ASSERT( mIpPdvHosts( tLocalPdvInd ),
                                "Pdv_Host_Manager::communicate_shared_interpolation_pdv_ids() - %s",
                                "Pdv host does not exist on this processor" );

                        // Check that PDV Id is valid
                        MORIS_ASSERT( mIpPdvHosts( tLocalPdvInd )->get_pdv_id( tPdvType ) != gNoID,
                                "Pdv_Host_Manager::communicate_shared_interpolation_pdv_ids() - %s",
                                "Local Pdv Id is invalid" );

                        // Check that PDV host is indeed owned by this processor
                        MORIS_ASSERT( ( mIpPdvHosts( tLocalPdvInd )->get_pdv_owning_processor() ) == par_rank(),
                                "Pdv_Host_Manager::communicate_shared_interpolation_pdv_ids() - %s",
                                "Pdv not owned by this processor" );

                        // store Id of PDV for given type on list of requesting processor
                        tSharedPdvIdList( Ik )( Ii ) = mIpPdvHosts( tLocalPdvInd )->get_pdv_id( tPdvType );
                    }
                }

                // receiving list
                moris::Cell< Matrix< DDUMat > > tMatsToReceive2;

                // FIXME: should not be needed if done correctly
                barrier();

                // Communicate owned pdv Id back to the processor with the shared pdv
                communicate_mats(
                        mCommTable,
                        tSharedPdvIdList,
                        tMatsToReceive2 );

                // check that number of received vectors is consistent with communication table
                MORIS_ERROR( tMatsToReceive2.size() == tNumCommProcs,
                        "Pdv_Host_Manager::communicate_shared_interpolation_pdv_ids - incorrect data received\n" );

                // Loop over all communication processors and assign received PDV Ids
                for ( moris::uint Ik = 0; Ik < tNumCommProcs; Ik++ )
                {
                    // number of PDV Ids sent by communicating processor
                    uint tNumberOfReceivedIds = tMatsToReceive2( Ik ).numel();

                    // Check that number of received PDV Ids is consistent with original request
                    MORIS_ERROR( tSharedPdvPosLocal( Ik ).numel() == tNumberOfReceivedIds,
                            "Pdv_Host_Manager::communicate_shared_interpolation_pdv_ids - %s",
                            "mismatch between requested and received PDV Ids.\n" );

                    // loop over all received PDV Ids
                    for ( uint Ii = 0; Ii < tNumberOfReceivedIds; ++Ii )
                    {
                        // get PDV host index
                        uint tPdvHostIndex = tSharedPdvPosLocal( Ik )( Ii );

                        // Check that PDV host exists
                        MORIS_ASSERT( mIpPdvHosts( tPdvHostIndex ),
                                "Pdv_Host_Manager::communicate_shared_interpolation_pdv_ids - %s",
                                "requesting PDV host does not exist any longer.\n" );

                        // Get received PDV Id
                        uint tReceivedPdvId = tMatsToReceive2( Ik )( Ii );

                        // Check that received Id is valid
                        MORIS_ASSERT( tReceivedPdvId != MORIS_UINT_MAX,
                                "Pdv_Host_Manager::communicate_shared_interpolation_pdv_ids() - %s",
                                "received Ids is invalid.\n" );

                        // Set received PDV Id for current type
                        mIpPdvHosts( tPdvHostIndex )->set_pdv_id( tPdvType, tReceivedPdvId );
                    }
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Pdv_Host_Manager::communicate_shared_intersection_node_pdv_ids()
        {
            // Build communication table map to determine the right position for each processor rank.
            Matrix< DDSMat > tCommTableMap( mCommTable.max() + 1, 1, -1 );

            moris::uint tNumCommProcs = mCommTable.numel();

            // Loop over communication table to fill the communication table map
            for ( moris::uint Ik = 0; Ik < tNumCommProcs; Ik++ )
            {
                tCommTableMap( mCommTable( Ik ), 0 ) = Ik;
            }

            moris::uint tCounter       = 0;
            moris::uint tSharedCounter = 0;

            moris::Cell< Matrix< DDUMat > > tSharedPdvPosGlobal( tNumCommProcs );
            moris::Cell< Matrix< DDUMat > > tSharedPdvPosLocal( tNumCommProcs );

            // Set Mat to store number of shared pdv per processor
            Matrix< DDUMat > tNumSharedPdvsPerProc( tNumCommProcs, 1, 0 );

            // Loop over pdvs
            for ( moris::uint Ij = 0; Ij < mIntersectionNodes.size(); Ij++ )
            // for ( moris::uint Ij = 0; Ij < mGenMeshMap.size(); Ij++ )
            {

                uint tMeshNodeIndex = Ij;
                if ( mGenMeshMapIsInitialized )
                {
                    tMeshNodeIndex = mGenMeshMap( Ij );
                }

                // Check if pdv at this position is not NULL
                if ( mIntersectionNodes( tMeshNodeIndex ) != nullptr )
                {
                    // Check if owning processor is this processor
                    if ( mIntersectionNodes( tMeshNodeIndex )->get_owner() != par_rank() )
                    {
                        // get owning processor
                        moris::moris_id tProcID = mIntersectionNodes( tMeshNodeIndex )->get_owner();

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
                    tSharedPdvPosGlobal( Ik ).set_size( tNumSharedPdvsPerProc( Ik, 0 ), 1 );
                    tSharedPdvPosLocal( Ik ).set_size( tNumSharedPdvsPerProc( Ik, 0 ), 1 );
                }
            }

            // Temporary Mat to add external pdv ids at the next spot in the matrix which will be communicated
            Matrix< DDUMat > tSharedPdvPosPerProc( tNumCommProcs, 1, 0 );

            // Loop over pdvs
            for ( moris::uint Ij = 0; Ij < mIntersectionNodes.size(); Ij++ )
            // for ( moris::uint Ij = 0; Ij < mGenMeshMap.size(); Ij++ )
            {
                uint tMeshNodeIndex = Ij;
                if ( mGenMeshMapIsInitialized )
                {
                    tMeshNodeIndex = mGenMeshMap( Ij );
                }

                // Check if pdv at this position is not NULL
                if ( mIntersectionNodes( tMeshNodeIndex ) != nullptr )
                {
                    // Check if owning processor is this processor
                    if ( mIntersectionNodes( tMeshNodeIndex )->get_owner() != par_rank() )
                    {
                        // Get owning processor
                        moris::moris_id tProcID = mIntersectionNodes( tMeshNodeIndex )->get_owner();

                        moris::sint tProcIdPos = tCommTableMap( tProcID );

                        // Add owning processor id to moris::Mat
                        tSharedPdvPosGlobal( tProcIdPos )( tSharedPdvPosPerProc( tProcIdPos ) ) =
                                mIntersectionNodes( tMeshNodeIndex )->get_id();

                        // Add pdv position to Mat
                        tSharedPdvPosLocal( tProcIdPos )( tSharedPdvPosPerProc( tProcIdPos ) ) = tCounter;

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
                tSharedPdvIdList( Ik ).set_size( tMatsToReceive( Ik ).numel(), 1 );
            }

            // Loop over all received positions and get the pdv id of the owning pdv
            for ( moris::uint Ik = 0; Ik < tMatsToReceive.size(); Ik++ )
            {
                for ( moris::uint Ii = 0; Ii < tMatsToReceive( Ik ).numel(); Ii++ )
                {
                    // Get owned pdv Id
                    auto tIter = mIGVertexIdtoIndMap.find( tMatsToReceive( Ik )( Ii ) );

                    moris::uint tLocalPdvInd = tIter->second;

                    uint tMeshNodeIndex = tLocalPdvInd;
                    if ( mGenMeshMapIsInitialized )
                    {
                        tMeshNodeIndex = mGenMeshMap( tLocalPdvInd );
                    }

                    MORIS_ASSERT( ( mIntersectionNodes( tMeshNodeIndex )->get_owner() ) == par_rank(),
                            "Pdv_Host_Manager::communicate_shared_pdv_ids(): Pdv not owned by this processor" );

                    tSharedPdvIdList( Ik )( Ii ) = mIntersectionNodes( tMeshNodeIndex )->get_starting_pdv_id();
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

            // assemble Ids in list of shared pdv ids and assemble the corresponding positions
            for ( moris::uint Ik = 0; Ik < tMatsToReceive2.size(); Ik++ )
            {
                if ( tMatsToReceive2( Ik ).numel() >= 1 )
                {
                    tListSharedPdvIds( { tPdvPosCounter, tPdvPosCounter + tMatsToReceive2( Ik ).numel() - 1 }, { 0, 0 } ) =
                            tMatsToReceive2( Ik ).matrix_data();

                    tListSharedPdvPos( { tPdvPosCounter, tPdvPosCounter + tSharedPdvPosLocal( Ik ).numel() - 1 }, { 0, 0 } ) =
                            tSharedPdvPosLocal( Ik ).matrix_data();

                    tPdvPosCounter = tPdvPosCounter + tMatsToReceive2( Ik ).numel();
                }
            }

            if ( tListSharedPdvIds.numel() != 0 )
            {
                MORIS_ASSERT( tListSharedPdvIds.max() != MORIS_UINT_MAX,
                        "Pdv_Host_Manager::communicate_shared_pdv_ids(), communicated Ids not set correctly" );

                MORIS_ASSERT( tListSharedPdvPos.max() != MORIS_UINT_MAX,
                        "Pdv_Host_Manager::communicate_shared_pdv_ids(), positions for communicated Ids not set correctly" );
            }

            // Set the Id of the shared pdvs
            for ( moris::uint Ij = 0; Ij < tListSharedPdvIds.numel(); Ij++ )
            {
                uint tMeshNodeIndex = tListSharedPdvPos( Ij );
                if ( mGenMeshMapIsInitialized )
                {
                    tMeshNodeIndex = mGenMeshMap( tListSharedPdvPos( Ij ) );
                }

                mIntersectionNodes( tMeshNodeIndex )->set_starting_pdv_id( tListSharedPdvIds( Ij ) );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Pdv_Host_Manager::build_local_to_global_maps()
        {
            mOwnedPdvLocalToGlobalMap.set_size( mNumOwnedPdvs, 1, -1 );
            mOwnedAndSharedPdvLocalToGlobalMap.set_size( mNumOwnedAndSharedPdvs, 1, -1 );

            uint tCounter  = 0;
            uint tCounter2 = 0;

            // Loop over all different pdv types for IP node pdvs
            for ( moris::uint Ij = 0; Ij < mPdvTypeList.size(); Ij++ )
            {
                enum PDV_Type tPdvType = mPdvTypeList( Ij );

                // Loop over pdvs per type. Count number of pdvs per proc which have to be communicated
                for ( moris::uint Ib = 0; Ib < mIpPdvHosts.size(); Ib++ )
                {
                    // Check if PDV host exists
                    if ( mIpPdvHosts( Ib ) )
                    {
                        // Check if PDV exists for given type
                        if ( mIpPdvHosts( Ib )->get_pdv_exists( tPdvType ) )
                        {
                            // Check if owning processor is this processor
                            if ( mIpPdvHosts( Ib )->get_pdv_owning_processor() == par_rank() )
                            {
                                mOwnedPdvLocalToGlobalMap( tCounter++ ) = mIpPdvHosts( Ib )->get_pdv_id( tPdvType );
                            }
                            mOwnedAndSharedPdvLocalToGlobalMap( tCounter2++ ) = mIpPdvHosts( Ib )->get_pdv_id( tPdvType );
                        }
                    }
                }
            }

            // Loop over intersection node pdvs
            for ( moris::uint Ij = 0; Ij < mIntersectionNodes.size(); Ij++ )
            // for ( moris::uint Ij = 0; Ij < mGenMeshMap.size(); Ij++ )
            {
                uint tMeshNodeIndex = Ij;
                if ( mGenMeshMapIsInitialized )
                {
                    tMeshNodeIndex = mGenMeshMap( Ij );
                }

                if ( mIntersectionNodes( tMeshNodeIndex ) != nullptr )
                {
                    uint tNumPdvsOnIntersectionNode = mIntersectionNodes( tMeshNodeIndex )->get_num_pdvs();

                    for ( moris::uint Ik = 0; Ik < tNumPdvsOnIntersectionNode; Ik++ )
                    {
                        if ( mIntersectionNodes( tMeshNodeIndex )->get_owner() == par_rank() )
                        {
                            mOwnedPdvLocalToGlobalMap( tCounter++ ) = mIntersectionNodes( tMeshNodeIndex )->get_starting_pdv_id() + Ik;
                        }

                        mOwnedAndSharedPdvLocalToGlobalMap( tCounter2++ ) = mIntersectionNodes( tMeshNodeIndex )->get_starting_pdv_id() + Ik;
                    }
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Pdv_Host_Manager::create_pdv_ids()
        {
            // FIXME comments are missing ; no explanation at all what is going on here
            uint tPdvOffset = 0;

            if ( !( par_size() <= 1 ) )
            {
                this->communicate_check_if_owned_pdv_exists();
            }

            this->get_num_pdvs();

            MORIS_LOG_INFO( "System has a total of %-5i pdvs.", sum_all( mNumOwnedPdvs ) );

            if ( !( par_size() <= 1 ) )
            {
                tPdvOffset = this->communicate_pdv_offsets( mNumOwnedPdvs );
            }

            this->set_owned_pdv_ids( tPdvOffset );

            if ( !( par_size() <= 1 ) )
            {
                this->communicate_shared_pdv_ids();
            }

            this->build_local_to_global_maps();
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Pdv_Host_Manager::set_dQIdp( Cell< Matrix< DDRMat >* > adQIdp,
                Matrix< DDSMat >*                              aMap )
        {
            // Number
            sint tNumIQIs = adQIdp.size();

            // Create factory for resulting distributed vector
            sol::Matrix_Vector_Factory tDistributedFactory;

            // node map from
            sol::Dist_Map* tMap = tDistributedFactory.create_map( this->get_my_local_global_map() );

            // allocate dist vector
            sol::Dist_Vector* tdQIDp = tDistributedFactory.create_vector( tMap, tNumIQIs, false, true );

            for ( uint iIQI = 0; iIQI < (uint)tNumIQIs; iIQI++ )
            {
                MORIS_ASSERT( mIntersectionNodes.size() == adQIdp( iIQI )->n_rows(), "Dimension mismatch, rows nodes" );

                // iterate through intersection vertices
                for ( uint iNodes = 0; iNodes < mIntersectionNodes.size(); iNodes++ )
                {
                    if ( mIntersectionNodes( iNodes ) != nullptr )
                    {
                        moris_id tStartingPDVId = mIntersectionNodes( iNodes )->get_starting_pdv_id();

                        moris::Matrix< DDSMat > tPDVIds( 1, mIntersectionNodes( iNodes )->get_num_pdvs() );

                        for ( moris::uint iPdv = 0; iPdv < mIntersectionNodes( iNodes )->get_num_pdvs(); iPdv++ )
                        {
                            tPDVIds( iPdv ) = tStartingPDVId + iPdv;
                        }

                        Matrix< DDRMat > tIndividualSensitivity = adQIdp( iIQI )->get_row( iNodes );

                        tdQIDp->sum_into_global_values( tPDVIds, tIndividualSensitivity, iIQI );
                    }
                }
            }


            this->set_dQIdp_dist_vect( tdQIDp );
        }

        //--------------------------------------------------------------------------------------------------------------

    }    // namespace ge
}    // namespace moris
