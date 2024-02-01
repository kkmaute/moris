/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Interpolation_PDV_Host.cpp
 *
 */

#include "cl_GEN_Interpolation_PDV_Host.hpp"
#include "cl_GEN_PDV_Value.hpp"
#include "cl_GEN_PDV_Property.hpp"
#include "cl_GEN_PDV_Host_Manager.hpp"

#include "fn_stringify_matrix.hpp"

namespace moris::gen
{

    //--------------------------------------------------------------------------------------------------------------

    Interpolation_PDV_Host::Interpolation_PDV_Host(
            PDV_Host_Manager*       aPDVHostManager,
            const moris_index&      aNodeIndex,
            const moris_id&         aNodeId,
            const moris_index&      aNodeOwner,
            const Matrix< DDRMat >& aCoordinates )
    {
        // check for aPDVHostManager is valid
        MORIS_ERROR( aPDVHostManager != nullptr,
                "Interpolation_PDV_Host::Interpolation_PDV_Host - PDVHostManager does not exist.\n" );

        // set PDV host manager
        mPDVHostManager = aPDVHostManager;

        // store nodal data
        mNodeIndex = aNodeIndex;
        mNodeId    = aNodeId;
        mNodeOwner = aNodeOwner;

        // store coordinates
        mCoordinates = aCoordinates;

        // allocate size for storing PDVs
        mPDVs.resize( mPDVHostManager->get_max_num_pdvs(), nullptr );
    }

    //--------------------------------------------------------------------------------------------------------------

    Interpolation_PDV_Host::~Interpolation_PDV_Host()
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    uint
    Interpolation_PDV_Host::get_num_pdvs()
    {
        return mPDVs.size();
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Interpolation_PDV_Host::create_pdv( PDV_Type aPDVType, real aPDVVal )
    {
        const Matrix< DDSMat >& tPDVTypeMap = mPDVHostManager->get_pdv_type_map();

        sint tPDVIndex = tPDVTypeMap( static_cast< sint >( aPDVType ) );

        // Check PDV type
        MORIS_ASSERT( tPDVIndex != -1,
                "Interpolation_PDV_Host::create_pdv - PDV type does not exist at node with index %d.\n",
                mNodeIndex );

        // FIXME: not clear whether finding an existing PDV indicates an error or not;
        //        see also creating a PDV with a pointer
        if ( mPDVs( tPDVIndex ) == nullptr )
        {
            // Create a pdv with pdv value
            mPDVs( tPDVIndex ) = std::make_shared< PDV_Value >( aPDVVal );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Interpolation_PDV_Host::set_pdv_id(
            PDV_Type       aPDVType,
            const moris_id aCounterId )
    {
        const Matrix< DDSMat >& tPDVTypeMap = mPDVHostManager->get_pdv_type_map();

        sint tPDVIndex = tPDVTypeMap( static_cast< uint >( aPDVType ) );

        // Check PDV type
        MORIS_ASSERT( tPDVIndex != -1,
                "Interpolation_PDV_Host::set_pdv_id - PDV type does not exist at node with index %d.\n",
                mNodeIndex );

        mPDVs( tPDVIndex )->set_id( aCounterId );
    }

    //--------------------------------------------------------------------------------------------------------------

    moris_id
    Interpolation_PDV_Host::get_pdv_id( PDV_Type aPDVType )
    {
        const Matrix< DDSMat >& tPDVTypeMap = mPDVHostManager->get_pdv_type_map();

        sint tPDVIndex = tPDVTypeMap( static_cast< uint >( aPDVType ) );

        // Check PDV type
        MORIS_ASSERT( tPDVIndex != -1,
                "Interpolation_PDV_Host::get_pdv_id - PDV type does not exist at node with index %d.\n",
                mNodeIndex );

        return this->get_pdv_id( tPDVIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    moris_id
    Interpolation_PDV_Host::get_pdv_id( uint aPDVIndex )
    {
        // return PDV Id if PDV for given index exists; otherwise return -1
        if ( mPDVs( aPDVIndex ) )
        {
            return mPDVs( aPDVIndex )->get_id();
        }

        return -1;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Interpolation_PDV_Host::create_pdv(
            PDV_Type                    aPDVType,
            std::shared_ptr< Property > aPropertyPointer )
    {
        // Check that PDV has not already been created
        MORIS_ASSERT( aPropertyPointer != nullptr,
                "Interpolation_PDV_Host::create_pdv - property pointer is nullptr.\n" );

        const Matrix< DDSMat >& tPDVTypeMap = mPDVHostManager->get_pdv_type_map();

        sint tPDVIndex = tPDVTypeMap( static_cast< sint >( aPDVType ) );

        // Check PDV type
        MORIS_ASSERT( tPDVIndex != -1,
                "Interpolation_PDV_Host::create_pdv - PDV type does not exist at node with index %d.\n",
                mNodeIndex );

        // FIXME: not clear whether finding an existing PDV indicates an error or not;
        //        see also creating a PDV with a value
        // Check that PDV has not already been created
        // MORIS_ASSERT( mPDVs(tPDVIndex) == nullptr,
        //        "Interpolation_PDV_Host::create_pdv - PDV has already been created.\n");

        // Create a pdv with property pointer
        mPDVs( tPDVIndex ) = std::make_shared< PDV_Property >( aPropertyPointer );
    }

    //--------------------------------------------------------------------------------------------------------------

    bool
    Interpolation_PDV_Host::is_active_type( PDV_Type aPDVType )
    {
        const Matrix< DDSMat >& tPDVTypeMap = mPDVHostManager->get_pdv_type_map();

        sint tPDVIndex = tPDVTypeMap( static_cast< sint >( aPDVType ) );

        // Check PDV type
        MORIS_ASSERT( tPDVIndex != -1,
                "Interpolation_PDV_Host::is_active_type - PDV type does not exist at node with index %d.\n",
                mNodeIndex );

        // check that PDV pointer is valid
        MORIS_ASSERT( mPDVs( tPDVIndex ) != nullptr,
                "Interpolation_PDV_Host::is_active_type - PDV does not exist at node with index %d.\n",
                mNodeIndex );

        // return if active PDV
        return mPDVs( tPDVIndex )->is_active();
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Interpolation_PDV_Host::set_global_index_for_pdv_type(
            PDV_Type aPDVType,
            moris_id aId )
    {
        const Matrix< DDSMat >& tPDVTypeMap = mPDVHostManager->get_pdv_type_map();

        sint tPDVIndex = tPDVTypeMap( static_cast< sint >( aPDVType ) );

        // Check PDV type
        MORIS_ASSERT( tPDVIndex != -1,
                "Interpolation_PDV_Host::set_global_index_for_pdv_type - PDV type does not exist at node with index %d.\n",
                mNodeIndex );

        // check that PDV pointer is valid
        MORIS_ASSERT( mPDVs( tPDVIndex ) != nullptr,
                "Interpolation_PDV_Host::set_global_index_for_pdv_type - PDV does not exist at node with index %d.\n",
                mNodeIndex );

        mPDVs( tPDVIndex )->set_id( aId );
    }

    //--------------------------------------------------------------------------------------------------------------

    uint
    Interpolation_PDV_Host::get_global_index_for_pdv_type( PDV_Type aPDVType )
    {
        const Matrix< DDSMat >& tPDVTypeMap = mPDVHostManager->get_pdv_type_map();

        sint tPDVIndex = tPDVTypeMap( static_cast< sint >( aPDVType ) );

        // Check PDV type
        MORIS_ASSERT( tPDVIndex != -1,
                "Interpolation_PDV_Host::get_global_index_for_pdv_type - PDV type does not exist at node with index %d.\n",
                mNodeIndex );

        // check that PDV pointer is valid
        MORIS_ASSERT( mPDVs( tPDVIndex ) != nullptr,
                "Interpolation_PDV_Host::get_global_index_for_pdv_type - PDV does not exist at node with index %d.\n",
                mNodeIndex );

        // Return id from map
        return mPDVs( tPDVIndex )->get_id();
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDSMat >
    Interpolation_PDV_Host::get_all_global_indices()
    {
        // Initialize global IDs
        Matrix< DDSMat > tGlobPDVId( mPDVs.size(), 1 );

        // Loop over PDVs and get IDs
        uint tCounter = 0;
        for ( auto tPDV : mPDVs )
        {
            if ( tPDV != nullptr )
            {
                tGlobPDVId( tCounter++ ) = tPDV->get_id();
            }
        }

        // Resize IDs
        tGlobPDVId.resize( tCounter, 1 );

        // Return IDs
        return tGlobPDVId;
    }

    //--------------------------------------------------------------------------------------------------------------

    real
    Interpolation_PDV_Host::get_pdv_value( PDV_Type aPDVType )
    {
        const Matrix< DDSMat >& tPDVTypeMap = mPDVHostManager->get_pdv_type_map();

        sint tPDVIndex = tPDVTypeMap( static_cast< sint >( aPDVType ) );

        // Check PDV type
        MORIS_ASSERT( tPDVIndex != -1,
                "Interpolation_PDV_Host::get_global_index_for_pdv_type - PDV type doesn't exist." );

        // Check whether PDV exists
#ifdef MORIS_HAVE_DEBUG
        if ( mPDVs( tPDVIndex ) == nullptr )
        {
            std::cout << "PDV not found at current node index." << std::endl;
            std::cout << "PDV node index #" << mNodeIndex << " (ID: " << mNodeId << ")" << std::endl;
            std::cout << "PDV index: " << tPDVIndex << std::endl;
            std::cout << "Node Coordinates: " << ios::stringify_log( mCoordinates ) << std::endl;
            MORIS_ERROR( false, "Interpolation_PDV_Host::get_pdv_value() - PDV does not exist at node with index %d.\n", mNodeIndex );
        }
#endif    // MORIS_HAVE_DEBUG

        // Return value
        return mPDVs( tPDVIndex )->get_value( mNodeIndex, mCoordinates );
    }

    //--------------------------------------------------------------------------------------------------------------

    bool
    Interpolation_PDV_Host::get_pdv_exists( PDV_Type aPDVType )
    {
        const Matrix< DDSMat >& tPDVTypeMap = mPDVHostManager->get_pdv_type_map();

        sint tPDVIndex = tPDVTypeMap( static_cast< sint >( aPDVType ) );

        // Check PDV type
        MORIS_ASSERT( tPDVIndex != -1,
                "Interpolation_PDV_Host::get_pdv_exists - PDV type %d that doesn't exist on this host.",
                static_cast< sint >( aPDVType ) );

        // if PDV exists return true
        if ( mPDVs( tPDVIndex ) )
        {
            return true;
        }

        return false;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    Interpolation_PDV_Host::get_sensitivities( uint aPDVIndex )
    {
        // If PDV exists and is active, ask it for sensitivities; otherwise, return zero matrix
        if ( mPDVs( aPDVIndex ) and mPDVs( aPDVIndex )->is_active() )
        {
            return mPDVs( aPDVIndex )->get_sensitivities( mNodeIndex, mCoordinates );
        }

        return Matrix< DDRMat >( 0, 0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDSMat >
    Interpolation_PDV_Host::get_determining_adv_ids( uint aPDVIndex )
    {
        // If PDV exists and is active, ask it for depending ADV IDs; Otherwise, return zero matrix
        if ( mPDVs( aPDVIndex ) and mPDVs( aPDVIndex )->is_active() )
        {
            return mPDVs( aPDVIndex )->get_determining_adv_ids( mNodeIndex, mCoordinates );
        }

        return Matrix< DDSMat >( 0, 0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Interpolation_PDV_Host::print_state()
    {
        const Matrix< DDSMat >& tPDVTypeMap = mPDVHostManager->get_pdv_type_map();

        std::cout << "--------------------------------------------------------------\n";
        std::cout << " Interpolation_PDV_Host: \n";
        std::cout << " Current processor rank: " << par_rank() << std::endl;
        std::cout << " mNodeId:                " << mNodeId << std::endl;
        std::cout << " mNodeIndex:             " << mNodeIndex << std::endl;
        std::cout << " mNodeOwner:             " << mNodeOwner << std::endl;
        std::cout << " Number of PDV types:    " << tPDVTypeMap.n_rows() << std::endl;
        std::cout << " Number of PDVs:         " << mPDVs.size() << std::endl;

        print( mCoordinates, "mCoordinates" );

        for ( uint i = 0; i < tPDVTypeMap.n_rows(); ++i )
        {
            sint tIndex = tPDVTypeMap( i );

            if ( tIndex >= 0 )
            {
                if ( mPDVs( tIndex ) )
                {
                    std::cout << "PDV of type " << i << " - PDV - ID :" << mPDVs( tIndex )->get_id() << std::endl;
                }
                else
                {
                    std::cout << "PDV of type " << i << " does not exist\n";
                }
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

}
