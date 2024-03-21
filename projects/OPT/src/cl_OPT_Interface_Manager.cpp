/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_OPT_Interface_Manager.cpp
 *
 */

#include "cl_OPT_Interface_Manager.hpp"
#include "fn_sum.hpp"
#include "fn_Parsing_Tools.hpp"
#include "cl_Communication_Tools.hpp"

namespace moris
{
    namespace opt
    {

        // -------------------------------------------------------------------------------------------------------------

        Interface_Manager::Interface_Manager(
                Parameter_List                                  aParameterList,
                Vector< std::shared_ptr< Criteria_Interface > > aInterfaces )
                : mInterfaces( aInterfaces )
        {
            // Set number of interfaces
            mNumInterfaces = aInterfaces.size();

            // Set parameters
            mSharedADVs = aParameterList.get< bool >( "shared_advs" );
            mParallel   = aParameterList.get< bool >( "parallel" );
            string_to_mat( aParameterList.get< std::string >( "num_processors_per_interface" ), mProcessorBoundaries );

            // Check processors
            if ( mParallel )
            {
                MORIS_ERROR( sum( mProcessorBoundaries ) <= par_size(),
                        "Sum of number of processors per interface must not exceed number of available processors" );
            }

            // Change from number of processors to processor boundaries
            for ( uint tInterfaceIndex = 1; tInterfaceIndex < mNumInterfaces; tInterfaceIndex++ )
            {
                mProcessorBoundaries( tInterfaceIndex ) += mProcessorBoundaries( tInterfaceIndex - 1 );
            }
        }

        // -------------------------------------------------------------------------------------------------------------

        void
        Interface_Manager::initialize(
                Matrix< DDRMat >& aGlobalADVs,
                Matrix< DDRMat >& aGlobalLowerBounds,
                Matrix< DDRMat >& aGlobalUpperBounds,
                Matrix< IdMat >& )
        {
            // Set up ADVs
            Matrix< IdMat > tDummy;
            mInterfaces( 0 )->initialize( aGlobalADVs, aGlobalLowerBounds, aGlobalUpperBounds, tDummy );

            uint tCurrentGlobalADV = aGlobalADVs.length();

            Matrix< DDRMat > tLocalADVs;
            Matrix< DDRMat > tLocalLowerBounds;
            Matrix< DDRMat > tLocalUpperBounds;

            // ADVs per interface
            mNumADVsPerInterface.set_size( mNumInterfaces, 1, 0 );
            mNumADVsPerInterface( 0 ) = tCurrentGlobalADV;

            // Set upper limit of lower bounds and lower limit of upper bounds if shared
            if ( mSharedADVs )
            {
                for ( uint tInterfaceIndex = 1; tInterfaceIndex < mNumInterfaces; tInterfaceIndex++ )
                {
                    // Get the local bounds
                    Matrix< IdMat > tDummy;
                    mInterfaces( tInterfaceIndex )->initialize( tLocalADVs, tLocalLowerBounds, tLocalUpperBounds, tDummy );

                    // Compare with current global bounds
                    for ( uint tADVIndex = 0; tADVIndex < tCurrentGlobalADV; tADVIndex++ )
                    {
                        aGlobalLowerBounds( tADVIndex ) = std::max( aGlobalLowerBounds( tADVIndex ), tLocalLowerBounds( tADVIndex ) );
                        aGlobalUpperBounds( tADVIndex ) = std::min( aGlobalUpperBounds( tADVIndex ), tLocalUpperBounds( tADVIndex ) );
                    }
                }
            }

            // ADVs are not shared
            else
            {
                // Loop through local ADVs
                for ( uint tInterfaceIndex = 1; tInterfaceIndex < mNumInterfaces; tInterfaceIndex++ )
                {
                    // Get the local ADVs and bounds
                    Matrix< IdMat > tDummy;
                    mInterfaces( tInterfaceIndex )->initialize( tLocalADVs, tLocalLowerBounds, tLocalUpperBounds, tDummy );
                    mNumADVsPerInterface( tInterfaceIndex ) = tLocalADVs.length();

                    // Put into the global ADVs
                    aGlobalADVs.resize( tCurrentGlobalADV + mNumADVsPerInterface( tInterfaceIndex ), 1 );
                    aGlobalLowerBounds.resize( tCurrentGlobalADV + mNumADVsPerInterface( tInterfaceIndex ), 1 );
                    aGlobalUpperBounds.resize( tCurrentGlobalADV + mNumADVsPerInterface( tInterfaceIndex ), 1 );

                    for ( uint tADVIndex = 0; tADVIndex < mNumADVsPerInterface( tInterfaceIndex ); tADVIndex++ )
                    {
                        aGlobalADVs( tCurrentGlobalADV + tADVIndex )        = tLocalADVs( tADVIndex );
                        aGlobalLowerBounds( tCurrentGlobalADV + tADVIndex ) = tLocalLowerBounds( tADVIndex );
                        aGlobalUpperBounds( tCurrentGlobalADV + tADVIndex ) = tLocalUpperBounds( tADVIndex );
                    }

                    // Update number of global ADVs
                    tCurrentGlobalADV += mNumADVsPerInterface( tInterfaceIndex );
                }
            }
        }

        // -------------------------------------------------------------------------------------------------------------

        Matrix< DDRMat >
        Interface_Manager::perform( Matrix< DDRMat >& aNewADVs )
        {
            // Set up global criteria
            Matrix< DDRMat > tGlobalCriteria( 1, 1 );
            Matrix< DDRMat > tLocalCriteria;

            uint tCurrentGlobalCriteria = 0;

            mNumCriteriaPerInterface.set_size( mNumInterfaces, 1 );

            // Local ADVs
            Matrix< DDRMat > tLocalADVs( 0, 0 );

            // Get criteria in parallel if requested
            Vector< Matrix< DDRMat > > tReceiveMats;
            if ( mParallel )
            {
                // Assign new comm color and rank
                int tColor = 0;
                while ( par_rank() >= mProcessorBoundaries( tColor ) and tColor < int( mNumInterfaces ) )
                {
                    tColor++;
                }
                int tKey = mProcessorBoundaries( tColor ) - par_rank() - 1;

                // Split and get local criteria
                comm_split( tColor, tKey, "criteria_communicator" );

                tLocalADVs = get_local_advs( aNewADVs, tColor );

                tLocalCriteria = mInterfaces( tColor )->get_criteria( tLocalADVs );

                barrier( "criteria_evaluation" );
                comm_join();

                // Set up communication list and cells of mats
                Matrix< IdMat >          tCommunicationList( 0, 0 );
                Vector< Matrix< DDRMat > > tSendMats( 0 );

                // Send to other procs if I was rank 0 after split
                if ( !tKey )
                {
                    tCommunicationList.resize( par_size(), 1 );
                    tSendMats.resize( par_size() );
                    for ( int tParRank = 0; tParRank < par_size(); tParRank++ )
                    {
                        tSendMats( tParRank )          = tLocalCriteria;
                        tCommunicationList( tParRank ) = tParRank;
                    }
                }

                // Otherwise just receive from rank 0 procs
                else
                {
                    tCommunicationList.resize( mNumInterfaces, 1 );
                    tSendMats.resize( mNumInterfaces );
                    for ( uint tInterfaceIndex = 0; tInterfaceIndex < mNumInterfaces; tInterfaceIndex++ )
                    {
                        tCommunicationList( tInterfaceIndex ) = mProcessorBoundaries( tInterfaceIndex ) - 1;
                    }
                }

                // Communicate local criteria
                communicate_mats( tCommunicationList, tSendMats, tReceiveMats );

                // Only use necessary info if I was rank 0
                if ( !tKey )
                {
                    for ( uint tInterfaceIndex = 0; tInterfaceIndex < mNumInterfaces; tInterfaceIndex++ )
                    {
                        tReceiveMats( tInterfaceIndex ) = tReceiveMats( mProcessorBoundaries( tInterfaceIndex ) - 1 );
                    }
                    tReceiveMats( tColor ) = tLocalCriteria;
                }
            }

            // Get criteria in serial
            else
            {
                tReceiveMats.resize( mNumInterfaces );
                for ( uint tInterfaceIndex = 0; tInterfaceIndex < mNumInterfaces; tInterfaceIndex++ )
                {
                    tLocalADVs                      = get_local_advs( aNewADVs, tInterfaceIndex );
                    tReceiveMats( tInterfaceIndex ) = mInterfaces( tInterfaceIndex )->get_criteria( tLocalADVs );
                }
            }

            // Assemble on procs
            for ( uint tInterfaceIndex = 0; tInterfaceIndex < mNumInterfaces; tInterfaceIndex++ )
            {
                // Get the local criteria
                tLocalCriteria                              = tReceiveMats( tInterfaceIndex );
                mNumCriteriaPerInterface( tInterfaceIndex ) = tLocalCriteria.length();

                // Put into the global criteria
                tGlobalCriteria.resize( tCurrentGlobalCriteria + mNumCriteriaPerInterface( tInterfaceIndex ), 1 );
                for ( uint tCriteriaIndex = 0; tCriteriaIndex < mNumCriteriaPerInterface( tInterfaceIndex ); tCriteriaIndex++ )
                {
                    tGlobalCriteria( tCurrentGlobalCriteria + tCriteriaIndex ) = tLocalCriteria( tCriteriaIndex );
                }

                // Update number of global criteria
                tCurrentGlobalCriteria += mNumCriteriaPerInterface( tInterfaceIndex );
            }
            return tGlobalCriteria;
        }

        // -------------------------------------------------------------------------------------------------------------

        Matrix< DDRMat >
        Interface_Manager::compute_dcriteria_dadv()
        {
            // Set up global criteria gradients
            Matrix< DDRMat > tGlobalCriteriaGradients( sum( mNumCriteriaPerInterface ), sum( mNumADVsPerInterface ), 0.0 );
            Matrix< DDRMat > tLocalCriteriaGradients;
            uint             tCurrentGlobalCriteria = 0;
            uint             tCurrentGlobalADV      = 0;
            uint             tCurrentLocalADVs      = 0;

            // Get criteria gradients in parallel
            Vector< Matrix< DDRMat > > tReceiveMats;
            if ( mParallel )
            {
                // Assign new comm color and rank
                int tColor = 0;
                while ( par_rank() >= mProcessorBoundaries( tColor ) and tColor < int( mNumInterfaces ) )
                {
                    tColor++;
                }
                int tKey = mProcessorBoundaries( tColor ) - par_rank() - 1;

                // Split and get local criteria
                comm_split( tColor, tKey, "criteria_communicator" );
                tLocalCriteriaGradients = mInterfaces( tColor )->get_dcriteria_dadv();
                barrier( "d_criteria_dadv" );
                comm_join();

                // Set up communication list and cells of mats
                Matrix< IdMat >          tCommunicationList( 0, 0 );
                Vector< Matrix< DDRMat > > tSendMats( 0 );

                // Send to other procs if I was rank 0 after split
                if ( !tKey )
                {
                    tCommunicationList.resize( par_size(), 1 );
                    tSendMats.resize( par_size() );
                    for ( int tParRank = 0; tParRank < par_size(); tParRank++ )
                    {
                        tSendMats( tParRank )          = tLocalCriteriaGradients;
                        tCommunicationList( tParRank ) = tParRank;
                    }
                }

                // Otherwise just receive from rank 0 procs
                else
                {
                    tCommunicationList.resize( mNumInterfaces, 1 );
                    tSendMats.resize( mNumInterfaces );
                    for ( uint tInterfaceIndex = 0; tInterfaceIndex < mNumInterfaces; tInterfaceIndex++ )
                    {
                        tCommunicationList( tInterfaceIndex ) = mProcessorBoundaries( tInterfaceIndex ) - 1;
                    }
                }

                // Communicate local criteria
                communicate_mats( tCommunicationList, tSendMats, tReceiveMats );

                // Only use necessary info if I was rank 0
                if ( !tKey )
                {
                    for ( uint tInterfaceIndex = 0; tInterfaceIndex < mNumInterfaces; tInterfaceIndex++ )
                    {
                        tReceiveMats( tInterfaceIndex ) = tReceiveMats( mProcessorBoundaries( tInterfaceIndex ) - 1 );
                    }
                    tReceiveMats( tColor ) = tLocalCriteriaGradients;
                }
            }

            // Get local criteria gradients in serial
            else
            {
                tReceiveMats.resize( mNumInterfaces );
                for ( uint tInterfaceIndex = 0; tInterfaceIndex < mNumInterfaces; tInterfaceIndex++ )
                {
                    tReceiveMats( tInterfaceIndex ) = mInterfaces( tInterfaceIndex )->get_dcriteria_dadv();
                }
            }

            // Assemble into global list
            for ( uint tInterfaceIndex = 0; tInterfaceIndex < mNumInterfaces; tInterfaceIndex++ )
            {
                // Get the local criteria
                if ( mSharedADVs )
                {
                    tCurrentLocalADVs = tReceiveMats( tInterfaceIndex ).n_cols();
                }
                else
                {
                    tCurrentLocalADVs = mNumADVsPerInterface( tInterfaceIndex );
                }

                // Put into the global criteria
                tGlobalCriteriaGradients( { tCurrentGlobalCriteria, tCurrentGlobalCriteria + mNumCriteriaPerInterface( tInterfaceIndex ) - 1 },
                        { tCurrentGlobalADV, tCurrentGlobalADV + tCurrentLocalADVs - 1 } ) = tReceiveMats( tInterfaceIndex )( { 0, mNumCriteriaPerInterface( tInterfaceIndex ) - 1 },
                        { 0, tCurrentLocalADVs - 1 } );

                // Update number of global critera/advs
                tCurrentGlobalCriteria += mNumCriteriaPerInterface( tInterfaceIndex );
                if ( !mSharedADVs )
                {
                    tCurrentGlobalADV += mNumADVsPerInterface( tInterfaceIndex );
                }
            }
            return tGlobalCriteriaGradients;
        }

        // -------------------------------------------------------------------------------------------------------------

        Matrix< DDRMat >
        Interface_Manager::get_local_advs(
                Matrix< DDRMat > aGlobalADVs,
                uint             tThisInterface )
        {
            Matrix< DDRMat > tLocalADVs( 0, 0 );

            // ADVs are shared
            if ( mSharedADVs )
            {
                tLocalADVs = aGlobalADVs;
            }

            // ADVs are not shared
            else
            {
                tLocalADVs.resize( mNumADVsPerInterface( tThisInterface ), 1 );
                uint tGlobalADVIndex = 0;
                for ( uint tInterfaceIndex = 0; tInterfaceIndex < tThisInterface; tInterfaceIndex++ )
                {
                    tGlobalADVIndex += mNumADVsPerInterface( tInterfaceIndex );
                }
                for ( uint tADVIndex = 0; tADVIndex < mNumADVsPerInterface( tThisInterface ); tADVIndex++ )
                {
                    tLocalADVs( tADVIndex ) = aGlobalADVs( tGlobalADVIndex++ );
                }
            }

            return tLocalADVs;
        }

        // -------------------------------------------------------------------------------------------------------------

    }    // namespace opt
}    // namespace moris
