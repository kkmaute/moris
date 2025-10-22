/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_PDV_Host_Manager.cpp
 *
 */

#include "cl_GEN_PDV_Host_Manager.hpp"

#include <utility>
#include "cl_SOL_Matrix_Vector_Factory.hpp"
#include "cl_Communication_Tools.hpp"
#include "fn_trans.hpp"
#include "fn_eye.hpp"

// detailed logging
#include "cl_Tracer.hpp"

namespace moris::gen
{
    //--------------------------------------------------------------------------------------------------------------

    PDV_Host_Manager::PDV_Host_Manager(
            Node_Manager&          aNodeManager,
            Vector< std::string >& aRequestedQIs )
            : MSI::Design_Variable_Interface( aRequestedQIs )
            , mNodeManager( aNodeManager )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    PDV_Host_Manager::set_owned_adv_ids( const Vector< sint >& aOwnedADVIds )
    {
        mOwnedADVIds = aOwnedADVIds;
        mADVIdsSet   = true;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    PDV_Host_Manager::set_communication_table( const Matrix< IdMat >& aCommTable )
    {
        mCommTable = aCommTable;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< IdMat >
    PDV_Host_Manager::get_communication_table()    // FIXME
    {
        return mCommTable;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    PDV_Host_Manager::reset()
    {
        mIpPDVHosts.clear();
        mOwnedPDVLocalToGlobalMap.resize( 0, 0 );
        mOwnedAndSharedPDVLocalToGlobalMap.resize( 0, 0 );
        mNumOwnedPDVs          = 0;
        mNumOwnedAndSharedPDVs = 0;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    PDV_Host_Manager::get_ip_dv_types_for_set(
            const moris::moris_index      aIPMeshSetIndex,
            Vector< Vector< PDV_Type > >& aPDVTypes ) const
    {
        if ( mIpPDVTypes.size() > 0 )    // FIXME
        {
            aPDVTypes = mIpPDVTypes( aIPMeshSetIndex );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    PDV_Host_Manager::get_ig_dv_types_for_set(
            const moris_index             aIGMeshSetIndex,
            Vector< Vector< PDV_Type > >& aPDVTypes ) const
    {
        if ( mIgPDVTypes.size() > 0 )    // FIXME
        {
            aPDVTypes = mIgPDVTypes( aIGMeshSetIndex );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    PDV_Host_Manager::get_ip_unique_dv_types_for_set(
            const moris_index   aIPMeshSetIndex,
            Vector< PDV_Type >& aPDVTypes ) const
    {
        if ( mUniqueIpPDVTypes.size() > 0 )    // FIXME
        {
            aPDVTypes = mUniqueIpPDVTypes( aIPMeshSetIndex );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    PDV_Host_Manager::get_ig_unique_dv_types_for_set(
            const moris::moris_index aIGMeshSetIndex,
            Vector< PDV_Type >&      aPDVTypes ) const
    {
        if ( mUniqueIgPDVTypes.size() > 0 )    // FIXME
        {
            aPDVTypes = mUniqueIgPDVTypes( aIGMeshSetIndex );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    PDV_Host_Manager::get_ip_pdv_value(
            const Matrix< IndexMat >&   aNodeIndices,
            const Vector< PDV_Type >&   aPDVTypes,
            Vector< Matrix< DDRMat > >& aDvValues ) const
    {
        // Get the number of node indices requested
        uint tNumIndices = aNodeIndices.numel();

        // Get the number of dv types requested
        uint tNumTypes = aPDVTypes.size();

        // Cell size
        aDvValues.resize( tNumTypes );

        // loop over the node indices
        for ( uint tPDVTypeIndex = 0; tPDVTypeIndex < tNumTypes; tPDVTypeIndex++ )
        {
            // Matrix size
            aDvValues( tPDVTypeIndex ).set_size( tNumIndices, 1 );

            // loop over the requested dv types
            for ( uint tNodeIndex = 0; tNodeIndex < tNumIndices; tNodeIndex++ )
            {
                // get node index of PDV host
                uint tIdnx = aNodeIndices( tNodeIndex );

                MORIS_ASSERT( mIpPDVHosts.size() > tIdnx,
                        "PDV_Host_Manager::get_ip_pdv_value - IP PDV host does not exist at node with index. size to small %d\n",
                        tIdnx );

                // check that PDV host exists
                MORIS_ASSERT( mIpPDVHosts( tIdnx ),
                        "PDV_Host_Manager::get_ip_pdv_value - IP PDV host does not exist at node with index %d\n",
                        tIdnx );

                aDvValues( tPDVTypeIndex )( tNodeIndex ) =
                        mIpPDVHosts( tIdnx )->get_pdv_value( aPDVTypes( tPDVTypeIndex ) );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    PDV_Host_Manager::set_GenMeshMap( const Vector< moris_index >& aGenMeshMap )
    {
        mGenMeshMap              = aGenMeshMap;
        mGenMeshMapIsInitialized = true;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    PDV_Host_Manager::get_ig_pdv_value(
            const Vector< moris_index >& aNodeIndices,
            const Vector< PDV_Type >&    aPDVTypes,
            Vector< Matrix< DDRMat > >&  aDvValues,
            Vector< Vector< bool > >&    aIsActiveDv ) const
    {
        // Get the number of node indices requested
        uint tNumIndices = aNodeIndices.size();

        // Get the number of dv types requested
        uint tNumTypes = aPDVTypes.size();

        // Cell size
        aDvValues.resize( tNumTypes );
        aIsActiveDv.resize( tNumTypes );

        // loop over the node indices
        for ( uint tPDVTypeIndex = 0; tPDVTypeIndex < tNumTypes; tPDVTypeIndex++ )
        {
            // Matrix size
            aDvValues( tPDVTypeIndex ).set_size( tNumIndices, 1 );
            aIsActiveDv( tPDVTypeIndex ).resize( tNumIndices, false );

            // loop over the requested dv types
            for ( uint tNode = 0; tNode < tNumIndices; tNode++ )
            {
                // get the node index within GEN, if mtk mesh has been cleaned up, translate vertex index between GEN and MTK
                uint tGenMeshNodeIndex = aNodeIndices( tNode );
                if ( mGenMeshMapIsInitialized )
                {
                    tGenMeshNodeIndex = mGenMeshMap( tGenMeshNodeIndex );
                }

                if ( mNodeManager.node_depends_on_advs( tGenMeshNodeIndex ) )
                {
                    aDvValues( tPDVTypeIndex )( tNode )   = mNodeManager.get_node_coordinate_value( aNodeIndices( tNode ), static_cast< uint >( aPDVTypes( tPDVTypeIndex ) ) );
                    aIsActiveDv( tPDVTypeIndex )( tNode ) = true;
                }
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDSMat >&
    PDV_Host_Manager::get_my_local_global_map()
    {
        return mOwnedPDVLocalToGlobalMap;
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDSMat >&
    PDV_Host_Manager::get_my_local_global_overlapping_map()
    {
        return mOwnedAndSharedPDVLocalToGlobalMap;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    PDV_Host_Manager::get_ip_dv_ids_for_type_and_ind(
            const Matrix< IndexMat >&  aNodeIndices,
            const Vector< PDV_Type >&  aPDVTypes,
            Vector< Matrix< IdMat > >& aDvIds ) const
    {
        // get the number of node indices requested
        uint tNumIndices = aNodeIndices.length();

        // get the number of dv types requested
        uint tNumTypes = aPDVTypes.size();

        // set size for list of dv values
        aDvIds.resize( tNumTypes );

        // loop over the requested dv types
        for ( uint tPDVTypeIndex = 0; tPDVTypeIndex < tNumTypes; tPDVTypeIndex++ )
        {
            // resize vector for current PDV type
            aDvIds( tPDVTypeIndex ).resize( tNumIndices, 1 );

            // loop over the node indices
            for ( uint tNode = 0; tNode < tNumIndices; tNode++ )
            {
                // get node index
                uint tIdnx = aNodeIndices( tNode );

                // check that PDV host exists
                MORIS_ASSERT( mIpPDVHosts( tIdnx ),
                        "PDV_Host_Manager::get_ip_pdv_value - IP PDV host does not exist at node with index %d\n",
                        tIdnx );

                aDvIds( tPDVTypeIndex )( tNode ) =    //
                        mIpPDVHosts( tIdnx )->        //
                        get_global_index_for_pdv_type( aPDVTypes( tPDVTypeIndex ) );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    PDV_Host_Manager::get_ig_dv_ids_for_type_and_ind(
            const Vector< moris_index >&     aNodeIndices,
            const Vector< PDV_Type >&        aPDVTypes,
            Vector< Vector< moris_index > >& aDvIds ) const
    {
        // get the number of node indices requested
        uint tNumIndices = aNodeIndices.size();

        // get the number of dv types requested
        uint tNumTypes = aPDVTypes.size();

        // set size for list of dv values
        aDvIds.resize( tNumTypes );

        // loop over the requested dv types
        for ( uint tPDVTypeIndex = 0; tPDVTypeIndex < tNumTypes; tPDVTypeIndex++ )
        {
            aDvIds( tPDVTypeIndex ).resize( tNumIndices, -1 );

            // loop over the node indices
            for ( uint tNode = 0; tNode < tNumIndices; tNode++ )
            {
                uint tGenMeshNodeIndex = aNodeIndices( tNode );
                if ( mGenMeshMapIsInitialized )
                {
                    tGenMeshNodeIndex = mGenMeshMap( aNodeIndices( tNode ) );
                }

                if ( mNodeManager.node_depends_on_advs( tGenMeshNodeIndex ) )
                {
                    aDvIds( tPDVTypeIndex )( tNode ) = mNodeManager.get_derived_node_starting_pdv_id( tGenMeshNodeIndex )
                                                     + static_cast< uint >( aPDVTypes( tPDVTypeIndex ) );
                }
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    PDV_Host_Manager::get_ip_requested_dv_types( Vector< PDV_Type >& aPDVTypes ) const
    {
        aPDVTypes = mRequestedIpPDVTypes;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    PDV_Host_Manager::get_ig_requested_dv_types( Vector< PDV_Type >& aPDVTypes )
    {
        aPDVTypes = mRequestedIgPDVTypes;
    }

    //--------------------------------------------------------------------------------------------------------------

    void PDV_Host_Manager::set_interpolation_pdv_types(
            const Vector< Vector< Vector< PDV_Type > > >& aPDVTypes )
    {
        // Get number of sets
        uint tNumSets = aPDVTypes.size();

        // Set PDV types
        mIpPDVTypes = aPDVTypes;
        mUniqueIpPDVTypes.resize( tNumSets );

        // Loop over each mesh set
        for ( uint tMeshSetIndex = 0; tMeshSetIndex < tNumSets; tMeshSetIndex++ )
        {
            // Get number of unique PDV types for this set
            uint tNumUniquePDVs = 0;
            for ( uint tGroupIndex = 0; tGroupIndex < mIpPDVTypes( tMeshSetIndex ).size(); tGroupIndex++ )
            {
                tNumUniquePDVs += mIpPDVTypes( tMeshSetIndex )( tGroupIndex ).size();
            }
            mUniqueIpPDVTypes( tMeshSetIndex ).resize( tNumUniquePDVs );

            // Copy PDV types over to unique list that doesn't consider grouping
            uint tUniquePDVIndex = 0;
            for ( uint tGroupIndex = 0; tGroupIndex < mIpPDVTypes( tMeshSetIndex ).size(); tGroupIndex++ )
            {
                for ( uint tPDVIndex = 0; tPDVIndex < mIpPDVTypes( tMeshSetIndex )( tGroupIndex ).size(); tPDVIndex++ )
                {
                    mUniqueIpPDVTypes( tMeshSetIndex )( tUniquePDVIndex++ ) = mIpPDVTypes( tMeshSetIndex )( tGroupIndex )( tPDVIndex );
                }
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    PDV_Host_Manager::create_interpolation_pdv_hosts(
            const Vector< Vector< uint > >&   aNodeIndicesPerSet,
            const Vector< Vector< sint > >&   aNodeIdsPerSet,
            const Vector< Vector< uint > >&   aNodeOwnersPerSet,
            const Vector< Matrix< DDRMat > >& aNodeCoordinatesPerSet )
    {
        // Check that number of sets is consistent
        uint tNumSets = mIpPDVTypes.size();
        MORIS_ERROR( tNumSets == aNodeIndicesPerSet.size(),
                "PDV_Host_Manager::create_interpolation_pdv_hosts - inconsistent number of sets!" );

        // determine maximum node index used for sizing the pdv hosts
        uint tMax = 0;
        for ( uint iSetIndex = 0; iSetIndex < tNumSets; iSetIndex++ )
        {
            if ( aNodeIndicesPerSet( iSetIndex ).size() > 0 )
            {
                tMax = std::max( *std::max_element( aNodeIndicesPerSet( iSetIndex ).begin(), aNodeIndicesPerSet( iSetIndex ).end() ), tMax );
            }
        }

        // Initialize PDV hosts
        mIpPDVHosts.resize( tMax + 1, nullptr );

        // Create PDV hosts
        for ( uint tMeshSetIndex = 0; tMeshSetIndex < tNumSets; tMeshSetIndex++ )
        {
            // Get number of unique PDV types for this set
            uint tNumUniquePDVs = 0;

            for ( uint tGroupIndex = 0; tGroupIndex < mIpPDVTypes( tMeshSetIndex ).size(); tGroupIndex++ )
            {
                tNumUniquePDVs += mIpPDVTypes( tMeshSetIndex )( tGroupIndex ).size();
            }
            mUniqueIpPDVTypes( tMeshSetIndex ).resize( tNumUniquePDVs );

            // Copy PDV types over to unique list that doesn't consider grouping
            uint tUniquePDVIndex = 0;
            for ( uint tGroupIndex = 0; tGroupIndex < mIpPDVTypes( tMeshSetIndex ).size(); tGroupIndex++ )
            {
                for ( uint tPDVIndex = 0; tPDVIndex < mIpPDVTypes( tMeshSetIndex )( tGroupIndex ).size(); tPDVIndex++ )
                {
                    mUniqueIpPDVTypes( tMeshSetIndex )( tUniquePDVIndex++ ) = mIpPDVTypes( tMeshSetIndex )( tGroupIndex )( tPDVIndex );
                }
            }

            // get number of nodes in current set
            uint tNumberOfNodes = aNodeIndicesPerSet( tMeshSetIndex ).size();

            // Create PDV hosts on interpolation nodes
            for ( uint tNodeIndexOnSet = 0; tNodeIndexOnSet < tNumberOfNodes; tNodeIndexOnSet++ )
            {
                // Create new host or add unique PDVs
                moris_index      tNodeIndex       = aNodeIndicesPerSet( tMeshSetIndex )( tNodeIndexOnSet );
                moris_id         tNodeId          = aNodeIdsPerSet( tMeshSetIndex )( tNodeIndexOnSet );
                moris_index      tNodeOwner       = aNodeOwnersPerSet( tMeshSetIndex )( tNodeIndexOnSet );
                Matrix< DDRMat > tNodeCoordinates = aNodeCoordinatesPerSet( tMeshSetIndex ).get_row( tNodeIndexOnSet );

                // Create PDV host unless it already exists
                if ( mIpPDVHosts( tNodeIndex ) == nullptr and tNodeId >= 0 )
                {
                    mIpPDVHosts( tNodeIndex ) =
                            std::make_shared< Interpolation_PDV_Host >(
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
    PDV_Host_Manager::set_integration_pdv_types( Vector< Vector< Vector< PDV_Type > > > aPDVTypes )
    {
        // Check that number of sets is consistent
        uint tNumSets = aPDVTypes.size();

        // Set PDV types
        mIgPDVTypes = aPDVTypes;
        mUniqueIgPDVTypes.resize( tNumSets );

        // Unique PDV types
        for ( uint tMeshSetIndex = 0; tMeshSetIndex < tNumSets; tMeshSetIndex++ )
        {
            // Get number of unique PDVs
            uint tNumUniquePDVs = 0;
            for ( uint tGroupIndex = 0; tGroupIndex < mIgPDVTypes( tMeshSetIndex ).size(); tGroupIndex++ )
            {
                tNumUniquePDVs += mIgPDVTypes( tMeshSetIndex )( tGroupIndex ).size();
            }
            mUniqueIgPDVTypes( tMeshSetIndex ).resize( tNumUniquePDVs );

            // Copy PDV types over
            uint tUniquePDVIndex = 0;
            for ( uint tGroupIndex = 0; tGroupIndex < mIgPDVTypes( tMeshSetIndex ).size(); tGroupIndex++ )
            {
                for ( uint tPDVIndex = 0; tPDVIndex < mIgPDVTypes( tMeshSetIndex )( tGroupIndex ).size(); tPDVIndex++ )
                {
                    mUniqueIgPDVTypes( tMeshSetIndex )( tUniquePDVIndex++ ) =
                            mIgPDVTypes( tMeshSetIndex )( tGroupIndex )( tPDVIndex );
                }
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    PDV_Host_Manager::set_requested_interpolation_pdv_types( const Vector< PDV_Type >& aPDVTypes )
    {
        mRequestedIpPDVTypes = aPDVTypes;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    PDV_Host_Manager::set_requested_integration_pdv_types( const Vector< PDV_Type >& aPDVTypes )
    {
        mRequestedIgPDVTypes = aPDVTypes;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    PDV_Host_Manager::create_interpolation_pdv(
            uint     aNodeIndex,
            PDV_Type aPDVType,
            real     aPDVVal )
    {
        // Check that PDF host exists
        MORIS_ASSERT( mIpPDVHosts( aNodeIndex ),
                "PDV_Host_Manager::create_interpolation_pdv - IP PDV host does not exist at node with index %d\n",
                aNodeIndex );

        // Create PDV with given type and value
        mIpPDVHosts( aNodeIndex )->create_pdv( aPDVType, aPDVVal );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    PDV_Host_Manager::create_interpolation_pdv(
            uint                               aNodeIndex,
            PDV_Type                           aPDVType,
            const std::shared_ptr< Property >& aProperty )
    {
        // Check that PDV host exists
        MORIS_ASSERT( mIpPDVHosts( aNodeIndex ),
                "PDV_Host_Manager::create_interpolation_pdv - IP PDV host does not exist at node with index %d\n",
                aNodeIndex );

        // Create PDV with given type and property
        mIpPDVHosts( aNodeIndex )->create_pdv( aPDVType, aProperty );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    PDV_Host_Manager::remove_sensitivities_of_unused_variables(
            Vector< sint >&   aADVIds,
            Matrix< DDRMat >& aHostADVSensitivities )
    {
        // get number of ADV entries and number of sensitivities per IDs
        uint tNumADVs = aHostADVSensitivities.n_cols();
        uint tNumSens = aHostADVSensitivities.n_rows();

        // check for consistent lengths of ID vector and matrix of sensitivities
        MORIS_ERROR( aADVIds.size() == tNumADVs,
                "PDV_Host_Manager::remove_sensitivities_of_unused_variables - inconsistent number of ADVs (%d) vs. sensitivity values for interpolation PDVs (%zu)",
                tNumADVs,
                aADVIds.size() );

        // skip of ADV vector is empty
        if ( tNumADVs == 0 ) return;

        // skip if all ADV IDs are valid, i.e. positive
        if ( aADVIds.min() > -1 ) return;

        // initialize counter of used variables
        uint tCounter = 0;

        // loop over all ADV IDs
        for ( uint iv = 0; iv < aADVIds.size(); ++iv )
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
        aADVIds.resize( tCounter );
        aHostADVSensitivities.resize( tNumSens, tCounter );
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    PDV_Host_Manager::compute_dqi_dadv( const Vector< sint >& aFullADVIds, sol::Dist_Vector* aFulldQIdADV )
    {
        Tracer tTracer( "GEN", "PDV Host Manager", "compute dQi/dadv" );

        // Check for ADV IDs
        MORIS_ERROR( mADVIdsSet,
                "PDV Host Manager must have ADV IDs set before computing sensitivities." );

        // Modules that have QIs
        Vector< Module_Type > tQIModules = { Module_Type::XTK, Module_Type::FEM };

        // Create factory for resulting distributed vector
        sol::Matrix_Vector_Factory tDistributedFactory;

        // Create maps
        sol::Dist_Map* tOwnedADVMap = tDistributedFactory.create_map( mOwnedADVIds );
        sol::Dist_Map* tFullADVMap  = tDistributedFactory.create_map( aFullADVIds );

        // Initialize vectors for the total number of QIs. Initialize counter for global QI index
        sol::Dist_Vector* tdIQIdADV     = tDistributedFactory.create_vector( tOwnedADVMap, mRequestedQIs.size(), false, true );
        sol::Dist_Vector* tdFulldQIdADV = tDistributedFactory.create_vector( tFullADVMap, mRequestedQIs.size(), false, true );

        // Initialize derivatives of IQIs wrt Advs to zero
        tdIQIdADV->vec_put_scalar( 0.0 );

        // Indices of the requested QIs
        Vector< uint > tRequestedQIIndices = this->get_requested_QIs< uint >();

        // Loop over all modules that have QIs
        for ( const auto& tModule : tQIModules )
        {
            // Get the sensitivities for this module
            sol::Dist_Vector* tdQI = this->get_dQIdp( tModule );

            if ( tdQI == nullptr )
            {
                continue;    // no QIs have sensitivites for this module
            }

            // Get the requested QI indices for this module
            Vector< uint > tModuleRequestedQIIndices = this->get_requested_QIs< uint >( tModule );

            // Loop of interpolation PDV hosts
            for ( uint tPDVHostIndex = 0; tPDVHostIndex < mIpPDVHosts.size(); tPDVHostIndex++ )
            {
                // Check if PDV host exists and is owned by this processor
                if ( mIpPDVHosts( tPDVHostIndex ) and mIpPDVHosts( tPDVHostIndex )->get_pdv_owning_processor() == par_rank() )
                {
                    // Get number of PDVs
                    uint tNumPDVsOnHost = mIpPDVHosts( tPDVHostIndex )->get_num_pdvs();

                    // Assemble sensitivities
                    for ( uint tPDVIndex = 0; tPDVIndex < tNumPDVsOnHost; tPDVIndex++ )
                    {
                        // Get PDVs
                        moris_id tPDVID = mIpPDVHosts( tPDVHostIndex )->get_pdv_id( tPDVIndex );

                        // FIXME checking if the pdv is defined
                        if ( tPDVID != -1 )
                        {
                            // Get sensitivities
                            Matrix< DDRMat > tHostADVSensitivities =
                                    mIpPDVHosts( tPDVHostIndex )->get_sensitivities( tPDVIndex );

                            // Get ADV IDs
                            Vector< sint > tADVIds = mIpPDVHosts( tPDVHostIndex )->get_determining_adv_ids( tPDVIndex );

                            // remove sensitivities wrt unused variables
                            this->remove_sensitivities_of_unused_variables( tADVIds, tHostADVSensitivities );

                            // Loop over requested QIs in this module
                            for ( uint iModuleQIIndex = 0; iModuleQIIndex < tModuleRequestedQIIndices.size(); iModuleQIIndex++ )
                            {
                                // Get the index in the full requested list
                                uint tQIIndex = std::distance( tRequestedQIIndices.cbegin(), std::find( tRequestedQIIndices.cbegin(), tRequestedQIIndices.cend(), tModuleRequestedQIIndices( iModuleQIIndex ) ) );

                                Matrix< DDRMat > tIndividualSensitivity =
                                        ( *tdQI )( tPDVID, iModuleQIIndex ) * tHostADVSensitivities;

                                // Fill matrix
                                tdIQIdADV->sum_into_global_values( tADVIds, tIndividualSensitivity, tQIIndex );
                            }
                        }
                    }
                }
            }

            // Create ADV Host Sensitivities
            Matrix< DDRMat > tHostADVSensitivities;
            Matrix< DDRMat > tI;

            // Loop over intersection nodes for inserting
            for ( uint iNodeIndex = mNodeManager.get_number_of_background_nodes(); iNodeIndex < mNodeManager.get_total_number_of_nodes(); iNodeIndex++ )
            {
                if ( mNodeManager.node_depends_on_advs( iNodeIndex ) and mNodeManager.get_derived_node_owner( iNodeIndex ) == par_rank() )
                {
                    // Get starting ID and number of coordinates
                    uint tStartingGlobalIndex = mNodeManager.get_derived_node_starting_pdv_id( iNodeIndex );
                    uint tNumCoordinates      = mNodeManager.get_number_of_derived_node_pdvs( iNodeIndex );

                    // Parent sensitivities and ADV IDs
                    tHostADVSensitivities.set_size( 0, 0 );
                    eye( tNumCoordinates, tNumCoordinates, tI );
                    mNodeManager.append_dcoordinate_dadv_from_derived_node( iNodeIndex, tHostADVSensitivities, tI );
                    Vector< sint > tADVIds = mNodeManager.get_coordinate_determining_adv_ids_from_derived_node( iNodeIndex );

                    // remove sensitivities wrt unused variables
                    this->remove_sensitivities_of_unused_variables( tADVIds, tHostADVSensitivities );

                    // loop overall coordinate directions
                    for ( uint tCoordinateIndex = 0; tCoordinateIndex < tNumCoordinates; tCoordinateIndex++ )
                    {
                        // get PDV ID
                        moris_id tPDVID = tStartingGlobalIndex + tCoordinateIndex;

                        // Loop over requested QIs for this module
                        for ( uint iModuleQIIndex = 0; iModuleQIIndex < tModuleRequestedQIIndices.size(); iModuleQIIndex++ )
                        {
                            // Get the index in the full requested list
                            uint tQIIndex = std::distance( tRequestedQIIndices.cbegin(), std::find( tRequestedQIIndices.cbegin(), tRequestedQIIndices.cend(), tModuleRequestedQIIndices( iModuleQIIndex ) ) );

                            Matrix< DDRMat > tIndividualSensitivity =
                                    ( *tdQI )( tPDVID, iModuleQIIndex ) * tHostADVSensitivities.get_row( tCoordinateIndex );

                            // Fill matrix
                            tdIQIdADV->sum_into_global_values( tADVIds, tIndividualSensitivity, tQIIndex );
                        }
                    }
                }
            }
        }    // end loop over modules

        // Global assembly
        tdIQIdADV->vector_global_assembly();

        // Import
        tdFulldQIdADV->import_local_to_global( *tdIQIdADV );

        // Add to any GQI sensitivities
        aFulldQIdADV->vec_plus_vec( 1.0, *tdFulldQIdADV, 1.0 );

        // Extract values
        Matrix< DDRMat > tFullSensitivity( 0, 0 );
        aFulldQIdADV->extract_copy( tFullSensitivity );
        tFullSensitivity = trans( tFullSensitivity );

        // Clean up
        delete tdIQIdADV;
        delete tdFulldQIdADV;

        return tFullSensitivity;
    }

    //--------------------------------------------------------------------------------------------------------------

    void PDV_Host_Manager::communicate_dof_types( Vector< enum PDV_Type >& aPDVTypeList )
    {
        // Get processor size
        int tSize = par_size();

        // Get number of local dv types
        moris::sint tNumLocalDvTypes = aPDVTypeList.size();

        // Get number of global dv types
        moris::sint tNumMaxGlobalDvTypes = sum_all( tNumLocalDvTypes );

        if ( par_rank() == 0 )
        {
            // Set size of of pdv type list = number of global types
            mPDVTypeList.resize( tNumMaxGlobalDvTypes );
        }

        // Create list containing the number of local dof types
        Vector< moris::sint > tNumLocalDvTypesList( tSize );

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
        Vector< moris::sint > tDvTypeOffset( tSize, 0 );

        // Fill the list with the corresponding offsets
        for ( int Ip = 1; Ip < tSize; ++Ip )
        {
            tDvTypeOffset( Ip ) = tDvTypeOffset( Ip - 1 ) + tNumLocalDvTypesList( Ip - 1 );
        }

        // Assemble list containing all used dof types. Dof types are not unique
        MPI_Gatherv(
                ( ( aPDVTypeList.data() ).data() ),
                tNumLocalDvTypes,
                MPI_UNSIGNED,
                ( mPDVTypeList.data() ).data(),
                ( tNumLocalDvTypesList.data() ).data(),
                ( tDvTypeOffset.data() ).data(),
                MPI_UNSIGNED,
                0,
                MPI_COMM_WORLD );

        // Temporary variable for mPDVTypeList size
        moris::uint tPDVTypeListSize;

        if ( par_rank() == 0 )
        {
            // Sort this created list
            std::sort( ( mPDVTypeList.data() ).data(), ( mPDVTypeList.data() ).data() + mPDVTypeList.size() );

            // use std::unique and std::distance to create list containing all used dof types. This list is unique
            auto last = std::unique( ( mPDVTypeList.data() ).data(), ( mPDVTypeList.data() ).data() + mPDVTypeList.size() );
            auto pos  = std::distance( ( mPDVTypeList.data() ).data(), last );

            mPDVTypeList.resize( pos );

            tPDVTypeListSize = mPDVTypeList.size();
        }

        // Bcast size of mPDVTypeList on processor 0
        MPI_Bcast( &tPDVTypeListSize, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD );

        // Resize mPDVTypeList on all processors
        mPDVTypeList.resize( tPDVTypeListSize );

        // Bcast unique mPDVTypeList to all processors
        MPI_Bcast( ( mPDVTypeList.data() ).data(), mPDVTypeList.size(), MPI_UNSIGNED, 0, MPI_COMM_WORLD );
    }

    //--------------------------------------------------------------------------------------------------------------

    void PDV_Host_Manager::create_dv_type_map()
    {
        // Get number of unique adofs of this equation object
        moris::uint tNumUniquePDVTypes = mPDVTypeList.size();

        // Get maximal dv type enum number
        moris::sint tMaxDvTypeEnumNumber = 0;

        // Loop over all pdv types to get the highest enum index
        for ( moris::uint Ii = 0; Ii < tNumUniquePDVTypes; Ii++ )
        {
            tMaxDvTypeEnumNumber = std::max( tMaxDvTypeEnumNumber, static_cast< int >( mPDVTypeList( Ii ) ) );
        }

        // +1 because c++ is 0 based
        tMaxDvTypeEnumNumber = tMaxDvTypeEnumNumber + 1;

        // Set size of mapping matrix
        mPDVTypeMap.set_size( tMaxDvTypeEnumNumber, 1, -1 );

        // Loop over all pdv types to create the mapping matrix
        for ( moris::uint Ii = 0; Ii < tNumUniquePDVTypes; Ii++ )
        {
            mPDVTypeMap( static_cast< int >( mPDVTypeList( Ii ) ), 0 ) = Ii;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void PDV_Host_Manager::count_owned_and_shared_pdvs()
    {
        // Loop over all different pdv types for IP node pdvs
        for ( moris::uint Ij = 0; Ij < mPDVTypeList.size(); Ij++ )
        {
            enum PDV_Type tPDVType = mPDVTypeList( Ij );

            for ( moris::uint Ib = 0; Ib < mIpPDVHosts.size(); Ib++ )
            {
                // Check if PDV host exists
                if ( mIpPDVHosts( Ib ) )
                {
                    // Check if PDV exists
                    if ( mIpPDVHosts( Ib )->get_pdv_exists( tPDVType ) )
                    {
                        // Check if owning processor is this processor
                        if ( mIpPDVHosts( Ib )->get_pdv_owning_processor() == par_rank() )
                        {
                            mNumOwnedPDVs++;
                        }
                        mNumOwnedAndSharedPDVs++;
                    }
                }
            }
        }

        // Loop over intersection node pdvs
        for ( uint iNodeIndex = mNodeManager.get_number_of_background_nodes(); iNodeIndex < mNodeManager.get_total_number_of_nodes(); iNodeIndex++ )
        {
            uint tMeshNodeIndex = iNodeIndex;
            if ( mGenMeshMapIsInitialized )
            {
                tMeshNodeIndex = mGenMeshMap( iNodeIndex );
            }

            if ( mNodeManager.node_depends_on_advs( tMeshNodeIndex ) )
            {
                uint tNumPDVsOnIntersectionNode = mNodeManager.get_number_of_derived_node_pdvs( tMeshNodeIndex );

                if ( mNodeManager.get_derived_node_owner( tMeshNodeIndex ) == par_rank() )
                {
                    mNumOwnedPDVs += tNumPDVsOnIntersectionNode;
                }
                mNumOwnedAndSharedPDVs += tNumPDVsOnIntersectionNode;
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    uint PDV_Host_Manager::communicate_offsets( uint aNumOwnedIDs )
    {
        // Create list containing the number of owned IDs on each processor
        Matrix< DDUMat > tNumOwnedIDsPerProcessor;

        // Broadcast number of owned IDs
        allgather_scalar( aNumOwnedIDs, tNumOwnedIDsPerProcessor );

        // Create ID offset list
        Matrix< DDUMat > tOffsetList( tNumOwnedIDsPerProcessor.numel(), 1, 0 );

        // Loop over all entries to create the offsets. Starting with 1
        for ( uint iProcessorIndex = 1; iProcessorIndex < tOffsetList.numel(); iProcessorIndex++ )
        {
            // Add the number of owned IDs of the previous processor to the offset of the previous processor
            tOffsetList( iProcessorIndex ) = tOffsetList( iProcessorIndex - 1 ) + tNumOwnedIDsPerProcessor( iProcessorIndex - 1 );
        }

        return tOffsetList( par_rank() );
    }

    //--------------------------------------------------------------------------------------------------------------

    void PDV_Host_Manager::set_owned_pdv_ids( uint aPDVOffset )
    {
        moris::uint tOwnedIdCounter = aPDVOffset;

        tOwnedIdCounter = this->set_owned_interpolation_pdv_ids( tOwnedIdCounter );

        this->set_owned_intersection_node_pdv_ids( tOwnedIdCounter );
    }

    //--------------------------------------------------------------------------------------------------------------

    moris_id
    PDV_Host_Manager::set_owned_interpolation_pdv_ids( moris_id aOwnedIdCounter )
    {
        moris_id tSaveOffset = aOwnedIdCounter;

        // Loop over all different pdv types for IP node pdvs
        for ( moris::uint Ij = 0; Ij < mPDVTypeList.size(); Ij++ )
        {
            enum PDV_Type tPDVType = mPDVTypeList( Ij );

            // Loop over pdvs per type.
            for ( moris::uint Ib = 0; Ib < mIpPDVHosts.size(); Ib++ )
            {
                // Check if PDV host exists
                if ( mIpPDVHosts( Ib ) )
                {
                    // Check if PDV exists
                    if ( mIpPDVHosts( Ib )->get_pdv_exists( tPDVType ) )
                    {
                        // Check if owning processor is this processor
                        if ( mIpPDVHosts( Ib )->get_pdv_owning_processor() == par_rank() )
                        {
                            mIpPDVHosts( Ib )->set_pdv_id( tPDVType, aOwnedIdCounter );
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
    PDV_Host_Manager::set_owned_intersection_node_pdv_ids( moris_id aOwnedIdCounter )
    {
        moris_id tSaveOffset = aOwnedIdCounter;

        // Loop over intersection node pdvs
        for ( uint iNodeIndex = mNodeManager.get_number_of_background_nodes(); iNodeIndex < mNodeManager.get_total_number_of_nodes(); iNodeIndex++ )
        {
            uint tMeshNodeIndex = iNodeIndex;
            if ( mGenMeshMapIsInitialized )
            {
                tMeshNodeIndex = mGenMeshMap( iNodeIndex );
            }

            if ( mNodeManager.node_depends_on_advs( tMeshNodeIndex ) and mNodeManager.get_derived_node_owner( tMeshNodeIndex ) == par_rank() )
            {
                mNodeManager.set_derived_node_starting_pdv_id( tMeshNodeIndex, aOwnedIdCounter );
                aOwnedIdCounter += mNodeManager.get_number_of_derived_node_pdvs( tMeshNodeIndex );
            }
        }

        moris_id tNumOwnedInterpolationIds = aOwnedIdCounter - tSaveOffset;

        MORIS_LOG_INFO( "System has a total of %-5i intersection node pdvs.", sum_all( tNumOwnedInterpolationIds ) );

        return aOwnedIdCounter;
    }

    //--------------------------------------------------------------------------------------------------------------

    void PDV_Host_Manager::communicate_shared_pdv_ids()
    {
        this->communicate_shared_interpolation_pdv_ids();
        this->communicate_shared_intersection_node_pdv_ids();
    }

    //--------------------------------------------------------------------------------------------------------------

    void PDV_Host_Manager::communicate_shared_interpolation_pdv_ids()
    {
        // Build communication table map to determine the right position for each processor rank.
        Vector< moris_id > tCommTableMap = build_communication_table_map( mCommTable );
        moris::uint        tNumCommProcs = mCommTable.numel();

        // FIXME: cannot have communication within following loop
        // Loop over all different pdv types for IP node pdvs
        for ( uint iPDVTypeIndex = 0; iPDVTypeIndex < mPDVTypeList.size(); iPDVTypeIndex++ )
        {
            enum PDV_Type tPDVType = mPDVTypeList( iPDVTypeIndex );

            // Define vector to store number of shared pdv per processor
            Matrix< DDUMat > tNumSharedPDVsPerProc( tNumCommProcs, 1, 0 );

            // Loop over pdvs per type. Count number of pdvs per processor which have to be communicated
            for ( uint iPDVHostIndex = 0; iPDVHostIndex < mIpPDVHosts.size(); iPDVHostIndex++ )
            {
                // Check if PDV host exists
                if ( mIpPDVHosts( iPDVHostIndex ) )
                {
                    // Check if PDV exists for given type
                    if ( mIpPDVHosts( iPDVHostIndex )->get_pdv_exists( tPDVType ) )
                    {
                        // Get owning processor rank
                        moris::moris_index tProcIndex = mIpPDVHosts( iPDVHostIndex )->get_pdv_owning_processor();

                        // Check if owning processor is not this processor
                        if ( tProcIndex != par_rank() )
                        {
                            // get position of owning process in communication table
                            moris::sint tProcIdPos = tCommTableMap( tProcIndex );

                            MORIS_ASSERT( tProcIdPos != -1,
                                    "PDV_Host_Manager::communicate_shared_interpolation_pdv_ids: Processor "
                                    "does not exist in communication table" );

                            // Add +1 to the processor number of shared dv per processor
                            tNumSharedPDVsPerProc( tProcIdPos )++;
                        }
                    }
                }
            }

            // Define cells of vectors to store PDV host IDs and indices of shared PDvs
            Vector< Matrix< DDUMat > > tSharedPDVHostIds( tNumCommProcs );
            Vector< Matrix< DDUMat > > tSharedPDVHostIndices( tNumCommProcs );

            // Set size of vectors to store PDV host IDs and indices of shared PDvs
            for ( uint iCommunicationProcIndex = 0; iCommunicationProcIndex < tNumCommProcs; iCommunicationProcIndex++ )
            {
                // Get number of pdvs shared with current processor
                uint tNumberOfSharedPDVs = tNumSharedPDVsPerProc( iCommunicationProcIndex );

                // if there are any PDVs shared with this processor set size of vectors
                if ( tNumberOfSharedPDVs != 0 )
                {
                    tSharedPDVHostIds( iCommunicationProcIndex ).set_size( tNumberOfSharedPDVs, 1 );
                    tSharedPDVHostIndices( iCommunicationProcIndex ).set_size( tNumberOfSharedPDVs, 1 );
                }
            }

            // Reset vector to store number of shared pdv per processor
            tNumSharedPDVsPerProc.fill( 0 );

            // Loop over pdvs
            for ( uint iPDVHostIndex = 0; iPDVHostIndex < mIpPDVHosts.size(); iPDVHostIndex++ )
            {
                // Check if PDV host exists
                if ( mIpPDVHosts( iPDVHostIndex ) )
                {
                    // Check if PDV exists for given type
                    if ( mIpPDVHosts( iPDVHostIndex )->get_pdv_exists( tPDVType ) )
                    {
                        // Get owning processor rank
                        moris::moris_index tProcIndex = mIpPDVHosts( iPDVHostIndex )->get_pdv_owning_processor();

                        // Check if owning processor is not this processor
                        if ( tProcIndex != par_rank() )
                        {
                            // get position of owning process in communication table
                            moris::sint tProcIdPos = tCommTableMap( tProcIndex );

                            // get position of next element in node ID list of owning processor
                            uint tProcListPos = tNumSharedPDVsPerProc( tProcIdPos );

                            // Add Id of PDV host to global list of owning processor
                            tSharedPDVHostIds( tProcIdPos )( tProcListPos ) = mIpPDVHosts( iPDVHostIndex )->get_id();

                            // Add local position of existing pdv hosts to local list of owning processor
                            tSharedPDVHostIndices( tProcIdPos )( tProcListPos ) = iPDVHostIndex;

                            // Add +1 for position of next element in node ID list of owning processor
                            tNumSharedPDVsPerProc( tProcIdPos )++;
                        }
                    }
                }
            }

            // Receive IDs of PDV hosts that this processor owns
            Vector< Matrix< DDUMat > > tOwnedIds;

            // FIXME: should not be needed if done correctly
            barrier();

            // Communicate position of shared PDVs to the owning processor
            communicate_mats(
                    mCommTable,
                    tSharedPDVHostIds,
                    tOwnedIds );

            // check that number of received vectors is consistent with communication table
            MORIS_ERROR( tOwnedIds.size() == tNumCommProcs,
                    "PDV_Host_Manager::communicate_shared_interpolation_pdv_ids - incorrect data received\n" );

            // Loop over all processors with PDV hosts owned by this processor
            for ( uint iCommunicationProcIndex = 0; iCommunicationProcIndex < tNumCommProcs; iCommunicationProcIndex++ )
            {
                // Loop over all PDV hosts for which IDs are requested
                for ( uint iOwnedPDVIndex = 0; iOwnedPDVIndex < tOwnedIds( iCommunicationProcIndex ).numel(); iOwnedPDVIndex++ )
                {
                    // requested PDV host id
                    uint tReqPDVHostId = tOwnedIds( iCommunicationProcIndex )( iOwnedPDVIndex );

                    // Get index of PDV host on this processor
                    auto        tIter         = mIPBaseVertexIdtoIndMap.find( tReqPDVHostId );
                    moris::uint tPDVHostIndex = tIter->second;

                    // Check that host id exists
                    MORIS_ASSERT( tIter != mIPBaseVertexIdtoIndMap.end(),
                            "PDV_Host_Manager::communicate_shared_interpolation_pdv_ids() - PDV host with ID %d does not exist on Proc %d.\n",
                            tReqPDVHostId,
                            par_rank() );

                    // Check that PDV host exists
                    MORIS_ASSERT( mIpPDVHosts( tPDVHostIndex ),
                            "PDV_Host_Manager::communicate_shared_interpolation_pdv_ids() - %s",
                            "PDV host does not exist on this processor" );

                    // Check that PDV host is indeed owned by this processor
                    MORIS_ASSERT( ( mIpPDVHosts( tPDVHostIndex )->get_pdv_owning_processor() ) == par_rank(),
                            "PDV_Host_Manager::communicate_shared_interpolation_pdv_ids() - %s",
                            "PDV not owned by this processor" );

                    // Check if PDF of given type exists
                    MORIS_ASSERT( mIpPDVHosts( tPDVHostIndex )->get_pdv_exists( tPDVType ),
                            "PDV_Host_Manager::communicate_check_if_owned_pdv_exists - PDV missing on Node with ID %d on Proc %d.\n",
                            tReqPDVHostId,
                            par_rank() );

                    // Send back PDV ID
                    tOwnedIds( iCommunicationProcIndex )( iOwnedPDVIndex ) = mIpPDVHosts( tPDVHostIndex )->get_pdv_id( tPDVType );
                }
            }

            // receiving list
            Vector< Matrix< DDUMat > > tSharedPDVIds;

            // FIXME: should not be needed if done correctly
            barrier();

            // Communicate owned PDV ID back to the processor with the shared pdv
            communicate_mats(
                    mCommTable,
                    tOwnedIds,
                    tSharedPDVIds );

            // check that number of received vectors is consistent with communication table
            MORIS_ERROR( tSharedPDVIds.size() == tNumCommProcs,
                    "PDV_Host_Manager::communicate_shared_interpolation_pdv_ids - incorrect data received\n" );

            // Loop over all communication processors and assign received PDV Ids
            for ( moris::uint iCommunicationProcIndex = 0; iCommunicationProcIndex < tNumCommProcs; iCommunicationProcIndex++ )
            {
                // number of PDV Ids sent by communicating processor
                uint tNumberOfReceivedIds = tSharedPDVIds( iCommunicationProcIndex ).numel();

                // Check that number of received PDV Ids is consistent with original request
                MORIS_ERROR( tSharedPDVHostIndices( iCommunicationProcIndex ).numel() == tNumberOfReceivedIds,
                        "PDV_Host_Manager::communicate_shared_interpolation_pdv_ids - %s",
                        "mismatch between requested and received PDV Ids.\n" );

                // loop over all received PDV Ids
                for ( uint iSharedPDVIdIndex = 0; iSharedPDVIdIndex < tNumberOfReceivedIds; ++iSharedPDVIdIndex )
                {
                    // get PDV host index
                    uint tPDVHostIndex = tSharedPDVHostIndices( iCommunicationProcIndex )( iSharedPDVIdIndex );

                    // Check that PDV host exists
                    MORIS_ASSERT( mIpPDVHosts( tPDVHostIndex ),
                            "PDV_Host_Manager::communicate_shared_interpolation_pdv_ids - %s",
                            "requesting PDV host does not exist any longer.\n" );

                    // Get received PDV Id
                    uint tReceivedPDVId = tSharedPDVIds( iCommunicationProcIndex )( iSharedPDVIdIndex );

                    // Check that received Id is valid
                    MORIS_ASSERT( tReceivedPDVId != MORIS_UINT_MAX,
                            "PDV_Host_Manager::communicate_shared_interpolation_pdv_ids() - %s",
                            "received Ids is invalid.\n" );

                    // Set received PDV Id for current type
                    mIpPDVHosts( tPDVHostIndex )->set_pdv_id( tPDVType, tReceivedPDVId );
                }
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void PDV_Host_Manager::communicate_shared_intersection_node_pdv_ids()
    {
        // Build communication table map to determine the right position for each processor rank.
        Vector< moris_id > tCommTableMap  = build_communication_table_map( mCommTable );
        moris::uint        tNumCommProcs  = mCommTable.numel();
        moris::uint        tSharedCounter = 0;

        Vector< Matrix< DDUMat > > tSharedPDVIds( tNumCommProcs );
        Vector< Matrix< DDUMat > > tSharedPDVPosLocal( tNumCommProcs );

        // Set Mat to store number of shared pdv per processor
        Matrix< DDUMat > tNumSharedPDVsPerProc( tNumCommProcs, 1, 0 );

        // Loop over pdvs
        for ( uint iNodeIndex = mNodeManager.get_number_of_background_nodes(); iNodeIndex < mNodeManager.get_total_number_of_nodes(); iNodeIndex++ )
        {
            uint tMeshNodeIndex = iNodeIndex;
            if ( mGenMeshMapIsInitialized )
            {
                tMeshNodeIndex = mGenMeshMap( iNodeIndex );
            }

            // Get derived node owner
            moris_index tProcIndex = mNodeManager.get_derived_node_owner( tMeshNodeIndex );

            // Check if node depends on ADVs
            if ( mNodeManager.node_depends_on_advs( tMeshNodeIndex ) and tProcIndex != par_rank() )
            {
                // Get proc position
                sint tProcIdPos = tCommTableMap( tProcIndex );

                // Add +1 to the processor number of shared dv per processor
                tNumSharedPDVsPerProc( tProcIdPos )++;

                tSharedCounter++;
            }
        }

        // Set size of the moris::Mats in the Cell
        for ( moris::uint Ik = 0; Ik < tNumCommProcs; Ik++ )
        {
            if ( tNumSharedPDVsPerProc( Ik, 0 ) != 0 )
            {
                tSharedPDVIds( Ik ).set_size( tNumSharedPDVsPerProc( Ik, 0 ), 1 );
                tSharedPDVPosLocal( Ik ).set_size( tNumSharedPDVsPerProc( Ik, 0 ), 1 );
            }
        }

        // Temporary Mat to add external pdv ids at the next spot in the matrix which will be communicated
        Matrix< DDUMat > tSharedPDVPosPerProc( tNumCommProcs, 1, 0 );

        // Loop over all nodes
        for ( uint iNodeIndex = mNodeManager.get_number_of_background_nodes(); iNodeIndex < mNodeManager.get_total_number_of_nodes(); iNodeIndex++ )
        {
            uint tMeshNodeIndex = iNodeIndex;
            if ( mGenMeshMapIsInitialized )
            {
                tMeshNodeIndex = mGenMeshMap( iNodeIndex );
            }

            // Check that node depends on ADVs
            if ( mNodeManager.node_depends_on_advs( tMeshNodeIndex ) )
            {
                // Get node owner
                moris_id tNodeOwner = mNodeManager.get_derived_node_owner( iNodeIndex );

                // Check that owner is not this proc
                if ( tNodeOwner != par_rank() )
                {
                    // Get owning processor position
                    moris::sint tProcIdPos = tCommTableMap( tNodeOwner );

                    // Add owning processor id to moris::Mat
                    tSharedPDVIds( tProcIdPos )( tSharedPDVPosPerProc( tProcIdPos ) ) = mNodeManager.get_derived_node_id( tMeshNodeIndex );

                    // Add pdv position to Mat
                    tSharedPDVPosLocal( tProcIdPos )( tSharedPDVPosPerProc( tProcIdPos ) ) = tMeshNodeIndex;
                    tSharedPDVPosPerProc( tProcIdPos )++;
                }
            }
        }

        // receiving list
        Vector< Matrix< DDUMat > > tMatsToReceive;
        barrier();

        // Communicate position of shared pdvs to the owning processor
        communicate_mats(
                mCommTable,
                tSharedPDVIds,
                tMatsToReceive );

        // Create List of Mats containing the shared node Ids
        Vector< Matrix< DDUMat > > tSharedPDVIdList( tNumCommProcs );

        // Loop over all Mats setting the size
        for ( moris::uint Ik = 0; Ik < tMatsToReceive.size(); Ik++ )
        {
            tSharedPDVIdList( Ik ).set_size( tMatsToReceive( Ik ).numel(), 1 );
        }

        // Loop over all received positions and get the pdv id of the owning pdv
        for ( moris::uint Ik = 0; Ik < tMatsToReceive.size(); Ik++ )
        {
            for ( moris::uint Ii = 0; Ii < tMatsToReceive( Ik ).numel(); Ii++ )
            {
                // Get owned pdv Id
                auto tIter = mIGVertexIdtoIndMap.find( tMatsToReceive( Ik )( Ii ) );

                moris::uint tLocalPDVInd = tIter->second;

                uint tMeshNodeIndex = tLocalPDVInd;
                if ( mGenMeshMapIsInitialized )
                {
                    tMeshNodeIndex = mGenMeshMap( tLocalPDVInd );
                }

                MORIS_ASSERT( mNodeManager.get_derived_node_owner( tMeshNodeIndex ) == par_rank(),
                        "PDV_Host_Manager::communicate_shared_pdv_ids(): PDV not owned by this processor" );

                tSharedPDVIdList( Ik )( Ii ) = mNodeManager.get_derived_node_starting_pdv_id( tMeshNodeIndex );
            }
        }

        Vector< Matrix< DDUMat > > tMatsToReceive2;

        barrier();

        // Communicate owned pdv Id back to the processor with the shared pdv
        communicate_mats(
                mCommTable,
                tSharedPDVIdList,
                tMatsToReceive2 );

        moris::uint tPDVPosCounter = 0;

        Matrix< DDUMat > tListSharedPDVIds( tSharedCounter, 1, MORIS_UINT_MAX );
        Matrix< DDUMat > tListSharedPDVPos( tSharedCounter, 1, MORIS_UINT_MAX );

        // assemble Ids in list of shared pdv ids and assemble the corresponding positions
        for ( moris::uint Ik = 0; Ik < tMatsToReceive2.size(); Ik++ )
        {
            if ( tMatsToReceive2( Ik ).numel() >= 1 )
            {
                tListSharedPDVIds( { tPDVPosCounter, tPDVPosCounter + tMatsToReceive2( Ik ).numel() - 1 }, { 0, 0 } ) =
                        tMatsToReceive2( Ik ).matrix_data();

                tListSharedPDVPos( { tPDVPosCounter, tPDVPosCounter + tSharedPDVPosLocal( Ik ).numel() - 1 }, { 0, 0 } ) =
                        tSharedPDVPosLocal( Ik ).matrix_data();

                tPDVPosCounter = tPDVPosCounter + tMatsToReceive2( Ik ).numel();
            }
        }

        if ( tListSharedPDVIds.numel() != 0 )
        {
            MORIS_ASSERT( tListSharedPDVIds.max() != MORIS_UINT_MAX,
                    "PDV_Host_Manager::communicate_shared_pdv_ids(), communicated Ids not set correctly" );

            MORIS_ASSERT( tListSharedPDVPos.max() != MORIS_UINT_MAX,
                    "PDV_Host_Manager::communicate_shared_pdv_ids(), positions for communicated Ids not set correctly" );
        }

        // Set the Id of the shared pdvs
        for ( moris::uint Ij = 0; Ij < tListSharedPDVIds.numel(); Ij++ )
        {
            uint tMeshNodeIndex = tListSharedPDVPos( Ij );
            if ( mGenMeshMapIsInitialized )
            {
                tMeshNodeIndex = mGenMeshMap( tListSharedPDVPos( Ij ) );
            }

            // Set starting PDV ID
            if ( not mNodeManager.is_background_node( tMeshNodeIndex ) )
            {
                mNodeManager.set_derived_node_starting_pdv_id( tMeshNodeIndex, tListSharedPDVIds( Ij ) );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void PDV_Host_Manager::build_local_to_global_maps()
    {
        mOwnedPDVLocalToGlobalMap.set_size( mNumOwnedPDVs, 1, -1 );
        mOwnedAndSharedPDVLocalToGlobalMap.set_size( mNumOwnedAndSharedPDVs, 1, -1 );

        uint tOwnedNodeCounter          = 0;
        uint tOwnedAndSharedNodeCounter = 0;

        // Loop over all different pdv types for IP node pdvs
        for ( moris::uint Ij = 0; Ij < mPDVTypeList.size(); Ij++ )
        {
            enum PDV_Type tPDVType = mPDVTypeList( Ij );

            // Loop over pdvs per type. Count number of pdvs per proc which have to be communicated
            for ( moris::uint Ib = 0; Ib < mIpPDVHosts.size(); Ib++ )
            {
                // Check if PDV host exists
                if ( mIpPDVHosts( Ib ) )
                {
                    // Check if PDV exists for given type
                    if ( mIpPDVHosts( Ib )->get_pdv_exists( tPDVType ) )
                    {
                        // Check if owning processor is this processor
                        if ( mIpPDVHosts( Ib )->get_pdv_owning_processor() == par_rank() )
                        {
                            mOwnedPDVLocalToGlobalMap( tOwnedNodeCounter++ ) = mIpPDVHosts( Ib )->get_pdv_id( tPDVType );
                        }
                        mOwnedAndSharedPDVLocalToGlobalMap( tOwnedAndSharedNodeCounter++ ) = mIpPDVHosts( Ib )->get_pdv_id( tPDVType );
                    }
                }
            }
        }

        // Loop over intersection node pdvs
        for ( uint iNodeIndex = mNodeManager.get_number_of_background_nodes(); iNodeIndex < mNodeManager.get_total_number_of_nodes(); iNodeIndex++ )
        {
            uint tMeshNodeIndex = iNodeIndex;
            if ( mGenMeshMapIsInitialized )
            {
                tMeshNodeIndex = mGenMeshMap( iNodeIndex );
            }

            if ( mNodeManager.node_depends_on_advs( tMeshNodeIndex ) )
            {
                uint tNumPDVsOnIntersectionNode = mNodeManager.get_number_of_derived_node_pdvs( tMeshNodeIndex );

                for ( uint iCoordinateIndex = 0; iCoordinateIndex < tNumPDVsOnIntersectionNode; iCoordinateIndex++ )
                {
                    if ( mNodeManager.get_derived_node_owner( tMeshNodeIndex ) == par_rank() )
                    {
                        mOwnedPDVLocalToGlobalMap( tOwnedNodeCounter++ ) = mNodeManager.get_derived_node_starting_pdv_id( tMeshNodeIndex ) + iCoordinateIndex;
                    }

                    mOwnedAndSharedPDVLocalToGlobalMap( tOwnedAndSharedNodeCounter++ ) = mNodeManager.get_derived_node_starting_pdv_id( tMeshNodeIndex ) + iCoordinateIndex;
                }
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void PDV_Host_Manager::create_pdv_ids()
    {
        // Start with no PDV offset
        uint tPDVOffset = 0;

        // Count and store the number of owned and shared PDVs on this processor
        this->count_owned_and_shared_pdvs();
        MORIS_LOG_INFO( "System has a total of %-5i pdvs.", sum_all( mNumOwnedPDVs ) );

        // If parallel, need to communicate PDV offsets to each processor
        if ( par_size() > 1 )
        {
            tPDVOffset = this->communicate_offsets( mNumOwnedPDVs );
        }

        // Set owned PDV IDs based on the PDV offset
        this->set_owned_pdv_ids( tPDVOffset );

        // Communicate owned PDV IDs to processors that share this ID
        if ( par_size() > 1 )
        {
            this->communicate_shared_pdv_ids();
        }

        // Build local to global maps
        this->build_local_to_global_maps();
    }

    //--------------------------------------------------------------------------------------------------------------

    // void
    // PDV_Host_Manager::set_dQIdp(
    //         const Vector< Matrix< DDRMat >* >& adQIdp,
    //         Matrix< DDSMat >*                  aMap )
    // {
    //     // Number
    //     sint tNumIQIs = adQIdp.size();

    //     // Create factory for resulting distributed vector
    //     sol::Matrix_Vector_Factory tDistributedFactory;

    //     // node map from
    //     sol::Dist_Map* tMap = tDistributedFactory.create_map( this->get_my_local_global_map() );

    //     // allocate dist vector
    //     sol::Dist_Vector* tdQIDp = tDistributedFactory.create_vector( tMap, tNumIQIs, false, true );

    //     for ( uint iIQI = 0; iIQI < (uint)tNumIQIs; iIQI++ )
    //     {
    //         // iterate through intersection vertices
    //         for ( uint iNodeIndex = mNodeManager.get_number_of_background_nodes(); iNodeIndex < mNodeManager.get_total_number_of_nodes(); iNodeIndex++ )
    //         {
    //             if ( mNodeManager.node_depends_on_advs( iNodeIndex ) )
    //             {
    //                 // Get number of PDVs and starting ID
    //                 moris_id tStartingPDVId = mNodeManager.get_derived_node_starting_pdv_id( iNodeIndex );
    //                 uint     tNumberOfPDVs  = mNodeManager.get_number_of_derived_node_pdvs( iNodeIndex );

    //                 moris::Matrix< DDSMat > tPDVIds( 1, tNumberOfPDVs );
    //                 for ( moris::uint iPDV = 0; iPDV < tNumberOfPDVs; iPDV++ )
    //                 {
    //                     tPDVIds( iPDV ) = tStartingPDVId + iPDV;
    //                 }

    //                 Matrix< DDRMat > tIndividualSensitivity = adQIdp( iIQI )->get_row( iNodeIndex );

    //                 tdQIDp->sum_into_global_values( tPDVIds, tIndividualSensitivity, iIQI );
    //             }
    //         }
    //     }

    //     this->set_dQIdp_dist_vect( tdQIDp );
    // }
}    // namespace moris::gen