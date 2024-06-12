/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Design.cpp
 *
 */

#include "cl_GEN_Design.hpp"

namespace moris::gen
{

    //--------------------------------------------------------------------------------------------------------------

    Design_Parameters::Design_Parameters( const Parameter_List& aParameterList )
            : mNumberOfRefinements( aParameterList.get< Vector< uint > >( "number_of_refinements" ) )
            , mRefinementMeshIndices( aParameterList.get< Vector< uint > >( "refinement_mesh_index" ) )
            , mRefinementFunctionIndex( aParameterList.get< sint >( "refinement_function_index" ) )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    Design_Parameters::Design_Parameters()
            : mNumberOfRefinements( {} )
            , mRefinementMeshIndices( {} )
            , mRefinementFunctionIndex( -1 )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    Design::Design( Design_Parameters aParameters )
            : mParameters( std::move( aParameters ) )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    const Vector< uint >& Design::get_num_refinements()
    {
        return mParameters.mNumberOfRefinements;
    }

    //--------------------------------------------------------------------------------------------------------------

    const Vector< uint >&
    Design::get_refinement_mesh_indices()
    {
        return mParameters.mRefinementMeshIndices;
    }

    //--------------------------------------------------------------------------------------------------------------

    sint
    Design::get_refinement_function_index()
    {
        return mParameters.mRefinementFunctionIndex;
    }

    //--------------------------------------------------------------------------------------------------------------

    sint
    Design::append_adv_info(
            mtk::Interpolation_Mesh* aMesh,
            Vector< sint >&          aOwnedADVIds,
            Matrix< IdMat >&         aOwnedijklIDs,
            sint                     aOffsetID,
            Vector< real >&          aLowerBounds,
            Vector< real >&          aUpperBounds )
    {
        // Store the ADV offset ID for this design
        mOffsetID = aOffsetID;

        // initialize counter for number coefficients
        uint tNumberOfCoefficients = 0;

        // Determine if level set will be created
        if ( this->intended_discretization() )
        {
            // Get discretization mesh index
            uint tDiscretizationMeshIndex = this->get_discretization_mesh_index();

            uint tMaxNumberOfCoefficients = aMesh->get_max_num_coeffs_on_proc( tDiscretizationMeshIndex );

            Matrix< IdMat >    tAllCoefIds( tMaxNumberOfCoefficients, 1, gNoID );
            Matrix< IndexMat > tAllCoefIndices( tMaxNumberOfCoefficients, 1, gNoIndex );
            Matrix< IdMat >    tAllCoefOwners( tMaxNumberOfCoefficients, 1, gNoID );
            Matrix< IdMat >    tAllCoefijklIDs( tMaxNumberOfCoefficients, 1, gNoID );

            for ( uint iNodeIndex = 0; iNodeIndex < aMesh->get_num_nodes(); iNodeIndex++ )
            {
                // check whether node has an underlying discretization on this processor
                bool tNodeHasDiscretization =
                        aMesh->get_mtk_vertex( iNodeIndex ).has_interpolation( tDiscretizationMeshIndex );

                // process only nodes that have discretization
                if ( tNodeHasDiscretization )
                {
                    // get indices and IDs from mtk mesh - FIXME: should return const &
                    const Matrix< IndexMat > tCoefIndices = aMesh->get_coefficient_indices_of_node(
                            iNodeIndex,
                            tDiscretizationMeshIndex );

                    const Matrix< IdMat > tCoefIds = aMesh->get_coefficient_IDs_of_node(
                            iNodeIndex,
                            tDiscretizationMeshIndex );

                    const Matrix< IdMat > tCoefOwners = aMesh->get_coefficient_owners_of_node(
                            iNodeIndex,
                            tDiscretizationMeshIndex );

                    Matrix< IdMat > tCoeffijklIDs;

                    if ( mtk::MeshType::HMR == aMesh->get_mesh_type() )
                    {
                        tCoeffijklIDs = aMesh->get_coefficient_ijkl_IDs_of_node(
                                iNodeIndex,
                                tDiscretizationMeshIndex );
                    }

                    // check that number of indices and ids are the same
                    MORIS_ASSERT( tCoefIds.numel() == tCoefIndices.numel(),
                            "distribute_advs - numbers of coefficients and ids do not match.\n" );

                    // get number of coefficients for current node
                    uint tNumCoefOfNode = tCoefIds.numel();

                    for ( uint tCoefIndex = 0; tCoefIndex < tNumCoefOfNode; ++tCoefIndex )
                    {
                        // get coefficient index
                        moris_index tCurrentIndex = tCoefIndices( tCoefIndex );

                        // check whether mesh coefficient has already been set
                        if ( tAllCoefIds( tCurrentIndex ) == -1 )
                        {
                            // increase field coefficient count
                            tNumberOfCoefficients++;

                            // populate mesh index to mesh coefficient id map
                            tAllCoefIds( tCurrentIndex ) = tCoefIds( tCoefIndex );

                            tAllCoefOwners( tCurrentIndex ) = tCoefOwners( tCoefIndex );

                            if ( aMesh->get_mesh_type() == mtk::MeshType::HMR )
                            {
                                tAllCoefijklIDs( tCurrentIndex ) = tCoeffijklIDs( tCoefIndex );
                            }
                        }
                        else
                        {
                            // check for consistency
                            MORIS_ASSERT( tAllCoefIds( tCurrentIndex ) == tCoefIds( tCoefIndex ),
                                    "distribute_advs - inconsistent index and ids.\n" );
                        }
                    }
                }
            }

            if ( par_size() > 1 )
            {
                communicate_missing_owned_coefficients(
                        aMesh,
                        tAllCoefIds,
                        tAllCoefOwners,
                        tAllCoefijklIDs,
                        tNumberOfCoefficients,
                        tDiscretizationMeshIndex,
                        aMesh->get_mesh_type() );
            }

            // Count number of owned coefficients
            uint tOwnedCounter = 0;
            for ( uint iCoefficientIndex = 0; iCoefficientIndex < tAllCoefIds.numel(); iCoefficientIndex++ )
            {
                if ( tAllCoefIds( iCoefficientIndex ) != gNoID && tAllCoefOwners( iCoefficientIndex ) == par_rank() )
                {
                    tOwnedCounter++;
                }
            }

            // Create vectors of owned coefficients
            Matrix< DDUMat > tOwnedCoefficients( tOwnedCounter, 1 );

            // Set owned coefficients
            tOwnedCounter = 0;
            for ( uint Ik = 0; Ik < tAllCoefIds.numel(); Ik++ )
            {
                if ( tAllCoefIds( Ik ) != gNoID && tAllCoefOwners( Ik ) == par_rank() )
                {
                    tOwnedCoefficients( tOwnedCounter++ ) = Ik;
                }
            }

            // Sizes of ID vectors
            uint tNumOwnedADVs         = aOwnedADVIds.size();
            uint tNumOwnedCoefficients = tOwnedCoefficients.numel();

            // Resize ID lists and bounds
            aOwnedADVIds.resize( tNumOwnedADVs + tNumOwnedCoefficients );
            aLowerBounds.resize( tNumOwnedADVs + tNumOwnedCoefficients );
            aUpperBounds.resize( tNumOwnedADVs + tNumOwnedCoefficients );
            aOwnedijklIDs.resize( tNumOwnedADVs + tNumOwnedCoefficients, 1 );
            Vector< sint > tSharedADVIds( tAllCoefIds.length() );

            // Add owned coefficients to lists
            for ( uint iOwnedCoefficient = 0; iOwnedCoefficient < tNumOwnedCoefficients; iOwnedCoefficient++ )
            {
                // Set the ADV ID as the offset plus the entity ID
                sint tADVId = aOffsetID
                            + aMesh->get_glb_entity_id_from_entity_loc_index(
                                    tOwnedCoefficients( iOwnedCoefficient ),
                                    mtk::EntityRank::BSPLINE,
                                    tDiscretizationMeshIndex );

                MORIS_ASSERT( tADVId - aOffsetID == tAllCoefIds( tOwnedCoefficients( iOwnedCoefficient ) ), "check if this is a problem" );

                aOwnedADVIds( tNumOwnedADVs + iOwnedCoefficient ) = tADVId;
                aLowerBounds( tNumOwnedADVs + iOwnedCoefficient ) = this->get_discretization_lower_bound();
                aUpperBounds( tNumOwnedADVs + iOwnedCoefficient ) = this->get_discretization_upper_bound();

                if ( aMesh->get_mesh_type() == mtk::MeshType::HMR )
                {
                    aOwnedijklIDs( tNumOwnedADVs + iOwnedCoefficient ) = tAllCoefijklIDs( tOwnedCoefficients( iOwnedCoefficient ) );
                }
            }

            // Add shared coefficients to field-specific list
            for ( uint iSharedCoefficientIndex = 0; iSharedCoefficientIndex < tAllCoefIds.length(); iSharedCoefficientIndex++ )
            {
                // Set the ADV ID as the offset plus the entity ID
                tSharedADVIds( iSharedCoefficientIndex ) = aOffsetID + tAllCoefIds( iSharedCoefficientIndex );
            }

            // Append the shared ADV IDs
            mSharedADVIDs.push_back( tSharedADVIds );

            // Update offset based on maximum ID
            return aOffsetID + aMesh->get_max_entity_id( mtk::EntityRank::BSPLINE, tDiscretizationMeshIndex );
        }
        return aOffsetID;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Design::communicate_missing_owned_coefficients(
            mtk::Interpolation_Mesh* aMesh,
            Matrix< IdMat >&         aAllCoefIds,
            Matrix< IdMat >&         aAllCoefOwners,
            Matrix< IdMat >&         aAllCoefijklIds,
            uint                     aNumCoeff,
            uint                     aDiscretizationMeshIndex,
            mtk::MeshType            aMeshType )
    {

        Matrix< IdMat > tCommTable = aMesh->get_communication_table();

        // Build communication table map to determine the right position for each processor rank. +1 because c++ is 0 based
        Matrix< DDSMat > tCommTableMap( tCommTable.max() + 1, 1, -1 );

        moris::uint tNumCommProcs = tCommTable.numel();

        // Loop over communication table to fill the communication table map
        for ( moris::uint Ik = 0; Ik < tNumCommProcs; Ik++ )
        {
            tCommTableMap( tCommTable( Ik ) ) = Ik;
        }

        Vector< Matrix< IdMat > > tSharedCoeffsPosGlobal( tNumCommProcs );
        Vector< Matrix< IdMat > > tSharedCoeffsijklIdGlobal( tNumCommProcs );

        // Set Mat to store number of shared coeffs per processor
        Matrix< DDUMat > tNumSharedCoeffsPerProc( tNumCommProcs, 1, 0 );

        // Count number of coeffs per proc which have to be communicated
        for ( moris::uint Ib = 0; Ib < aAllCoefIds.numel(); Ib++ )
        {
            // Check if coeffs at this position is not NULL
            if ( aAllCoefIds( Ib ) != gNoID && aAllCoefOwners( Ib ) != par_rank() )
            {

                // get owning processor
                moris::moris_id tProcID = aAllCoefOwners( Ib );

                moris::sint tProcIdPos = tCommTableMap( tProcID );

                MORIS_ASSERT( tProcIdPos != gNoID,
                        "communicate_missing_owned_coefficients: Map returns proc rank -1. Check communication table" );

                // Add +1 to the processor number of shared coeffs per processor
                tNumSharedCoeffsPerProc( tProcIdPos )++;
            }
        }

        // Set size of the moris::Mats in the Cell
        for ( moris::uint Ik = 0; Ik < tNumCommProcs; Ik++ )
        {
            if ( tNumSharedCoeffsPerProc( Ik ) != 0 )
            {
                tSharedCoeffsPosGlobal( Ik ).set_size( tNumSharedCoeffsPerProc( Ik ), 1 );
                tSharedCoeffsijklIdGlobal( Ik ).set_size( tNumSharedCoeffsPerProc( Ik ), 1 );
            }
        }

        // Temporary Mat to add external coeffs ids at the next spot in the matrix which will be communicated
        Matrix< DDUMat > tShredCoeffPosPerProc( tNumCommProcs, 1, 0 );

        // Loop over coeffs per type
        for ( moris::uint Ia = 0; Ia < aAllCoefIds.numel(); Ia++ )
        {
            // Check if coeffs at this position is not NULL
            if ( aAllCoefIds( Ia ) != gNoID && aAllCoefOwners( Ia ) != par_rank() )
            {
                // Get owning processor
                moris::uint tProcID = aAllCoefOwners( Ia );

                moris::sint tProcIdPos = tCommTableMap( tProcID );

                // Add owning processor id to moris::Mat
                tSharedCoeffsPosGlobal( tProcIdPos )( tShredCoeffPosPerProc( tProcIdPos ) ) =
                        aAllCoefIds( Ia );

                if ( aMeshType == mtk::MeshType::HMR )
                {
                    tSharedCoeffsijklIdGlobal( tProcIdPos )( tShredCoeffPosPerProc( tProcIdPos ) ) =
                            aAllCoefijklIds( Ia );
                }

                tShredCoeffPosPerProc( tProcIdPos )++;
            }
        }

        // receiving list
        Vector< Matrix< IdMat > > tMatsToReceive;
        Vector< Matrix< IdMat > > tMatsToReceiveijklID;

        barrier();

        // Communicate position of shared adofs to the owning processor
        communicate_mats(
                tCommTable,
                tSharedCoeffsPosGlobal,
                tMatsToReceive );

        barrier();

        if ( aMeshType == mtk::MeshType::HMR )
        {
            communicate_mats(
                    tCommTable,
                    tSharedCoeffsijklIdGlobal,
                    tMatsToReceiveijklID );

            MORIS_ASSERT( tMatsToReceiveijklID.size() == tMatsToReceive.size(), "size must be the same" );
        }

        map< moris_id, moris_index > tCoeffGlobalToLocalMap;
        aMesh->get_adof_map(
                aDiscretizationMeshIndex,
                tCoeffGlobalToLocalMap );

        // Loop over all Mats set dummy owned coeffs
        for ( moris::uint Ik = 0; Ik < tMatsToReceive.size(); Ik++ )
        {
            for ( moris::uint Ii = 0; Ii < tMatsToReceive( Ik ).numel(); Ii++ )
            {
                // Get owned coeff Index
                moris_id    tID            = tMatsToReceive( Ik )( Ii );
                moris_index tLocalCoeffInd = tCoeffGlobalToLocalMap.find( tID );

                if ( aAllCoefIds( tLocalCoeffInd ) == gNoID )
                {
                    aAllCoefIds( tLocalCoeffInd )    = tID;
                    aAllCoefOwners( tLocalCoeffInd ) = par_rank();

                    if ( aMeshType == mtk::MeshType::HMR )
                    {
                        aAllCoefijklIds( tLocalCoeffInd ) = tMatsToReceiveijklID( Ik )( Ii );
                    }

                    aNumCoeff++;
                }

                MORIS_ASSERT( aAllCoefIds( tLocalCoeffInd ) == tID,
                        "communicate_missing_owned_coefficients( ), coefficient IDs are not parallel consistent" );

                if ( aMeshType == mtk::MeshType::HMR )
                {
                    MORIS_ASSERT( aAllCoefijklIds( tLocalCoeffInd ) == tMatsToReceiveijklID( Ik )( Ii ),
                            "communicate_missing_owned_coefficients( ), coefficient ijkl IDs are not parallel consistent" );
                }
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

}    // namespace moris::gen
