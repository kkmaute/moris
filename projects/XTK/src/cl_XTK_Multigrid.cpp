/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Multigrid.cpp
 *
 */

#include "cl_XTK_Multigrid.hpp"

#include <iomanip>
#include <iostream>
#include <string>
#include <unordered_map>
#include <fstream>

#include "cl_Communication_Tools.hpp"
#include "linalg_typedefs.hpp"
#include "fn_assert.hpp"
#include "fn_sort.hpp"
#include "moris_typedefs.hpp"
#include "cl_MTK_Cell_Info_Factory.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh_Core.hpp"
#include "xtk_typedefs.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
#include "cl_XTK_Model.hpp"
#include "cl_MTK_Writer_Exodus.hpp"

namespace xtk
{

    Multigrid::Multigrid( xtk::Model* aXTKModelPtr )
            : mXTKModelPtr( aXTKModelPtr )
            , mMeshIndex( 0 )
    {
    }

    //------------------------------------------------------------------------------

    void
    Multigrid::create_fine_to_coarse_relationship()
    {
        // set size of fine to coarse list
        mFineBasisToCoarseBasis.resize( mNumBasis );

        moris::mtk::Interpolation_Mesh& tInterpolationMesh = mXTKModelPtr->get_background_mesh();

        // get num bg basis
        uint tNumBGBasis = tInterpolationMesh.get_num_basis( 0 );

        // loop over bg basis
        for ( uint Ik = 0; Ik < tNumBGBasis; Ik++ )
        {
            // Basis indices are consecutive and correspond to Ik

            // get enriched basis for this bg basis
            const Matrix< IndexMat >& tEnrichedCoeffsForBackgroundCoeffs = mXTKModelPtr->mEnrichedInterpMesh( 0 )->get_enriched_coefficients_at_background_coefficient( mMeshIndex, Ik );

            // loop over enriched basis for background basis Ik
            for ( uint Ia = 0; Ia < tEnrichedCoeffsForBackgroundCoeffs.numel(); Ia++ )
            {
                // get num coarse basis for this basis
                uint tNumCoarseBasis = tInterpolationMesh.get_num_coarse_basis_of_basis( 0, Ik );

                // get basis Ia
                moris_index tEnrichedBasisInd = tEnrichedCoeffsForBackgroundCoeffs( Ia );

                // resize for coarse basis for this fine basis
                mFineBasisToCoarseBasis( tEnrichedBasisInd ).set_size( tNumCoarseBasis, 1, -1 );

                // get subphases for enriched basis
                const Vector< moris::Matrix< moris::IndexMat > >& tSubphaseIndForEnrichedBasis = mXTKModelPtr->mEnrichment->get_subphases_loc_inds_in_enriched_basis();

                // get FIRST sub-phase index of basis. First because we assume the fine basis support is complete within the coarse one
                moris_index tFirstSubphaseInSupportIndex = tSubphaseIndForEnrichedBasis( tEnrichedBasisInd )( 0 );

                // get bg basis interpolating into tFirstSubphaseInSupportIndex ( Interpolation cell index corresponds to subphase index)
                Vector< moris_index > tBasisForSubphaseIndex = mXTKModelPtr->mEnrichment->mEnrichmentData( mMeshIndex ).mSubphaseBGBasisIndices( tFirstSubphaseInSupportIndex );

                // get enrichment level for bg basis interpolating into tFirstSubphaseInSupportIndex ( Interpolation cell index corresponds to subphase index)
                Vector< moris_index > tBasisEnrLevForSubphaseIndex = mXTKModelPtr->mEnrichment->mEnrichmentData( mMeshIndex ).mSubphaseBGBasisEnrLev( tFirstSubphaseInSupportIndex );

                // build map which maps bg basis index to entry in tBasisForSubphaseIndex/tBasisEnrLevForSubphaseIndex
                Mini_Map< moris_id, moris_id > tSubPhaseBasisMap = mXTKModelPtr->mEnrichment->construct_subphase_basis_to_basis_map( tBasisForSubphaseIndex );

                // loop over coarse basis
                for ( uint Ii = 0; Ii < tNumCoarseBasis; Ii++ )
                {
                    // get coarse bg basis
                    moris_index tCoarseBasisIndex = tInterpolationMesh.get_coarse_basis_index_of_basis( 0, Ik, Ii );

                    // find enrichment level for this coarse bg basis and subphase
                    moris_index tEnrichmentLev = tBasisEnrLevForSubphaseIndex( tSubPhaseBasisMap.find( tCoarseBasisIndex )->second );

                    // get enriched basis for this coarse bg basis and enrichment level
                    const Matrix< IndexMat >& tCoarseEnrichedCoeffsForCoarseBackgroundCoeffs = mXTKModelPtr->mEnrichedInterpMesh( 0 )
                                                                                                      ->get_enriched_coefficients_at_background_coefficient( mMeshIndex, tCoarseBasisIndex );

                    moris_index tEnrichedCoarseBasisIndex = tCoarseEnrichedCoeffsForCoarseBackgroundCoeffs( tEnrichmentLev );

                    // add enriched coarse basis to list
                    mFineBasisToCoarseBasis( tEnrichedBasisInd )( Ii ) = tEnrichedCoarseBasisIndex;
                }
            }
        }

        //        print( mFineBasisToCoarseBasis,"mFineBasisToCoarseBasis");
    }

    //------------------------------------------------------------------------------

    void
    Multigrid::create_coarse_to_fine_relationship()
    {
        // set size of fine to coarse list
        mCoarseBasisToFineBasis.resize( mNumBasis );

        moris::Vector< uint > tCounter( mNumBasis, 0 );

        for ( uint Ik = 0; Ik < mFineBasisToCoarseBasis.size(); Ik++ )
        {
            for ( uint Ii = 0; Ii < mFineBasisToCoarseBasis( Ik ).numel(); Ii++ )
            {
                tCounter( mFineBasisToCoarseBasis( Ik )( Ii ) )++;
            }
        }

        for ( uint Ik = 0; Ik < mCoarseBasisToFineBasis.size(); Ik++ )
        {
            mCoarseBasisToFineBasis( Ik ).set_size( tCounter( Ik ), 1, 0 );

            tCounter( Ik ) = 0;
        }

        for ( uint Ik = 0; Ik < mFineBasisToCoarseBasis.size(); Ik++ )
        {
            for ( uint Ii = 0; Ii < mFineBasisToCoarseBasis( Ik ).numel(); Ii++ )
            {
                mCoarseBasisToFineBasis( mFineBasisToCoarseBasis( Ik )( Ii ) )( tCounter( mFineBasisToCoarseBasis( Ik )( Ii ) ) ) = Ik;

                tCounter( mFineBasisToCoarseBasis( Ik )( Ii ) )++;
            }
        }

        //        print( mCoarseBasisToFineBasis,"mCoarseBasisToFineBasis");
    }

    //------------------------------------------------------------------------------

    void
    Multigrid::create_coarse_to_fine_weights()
    {
        moris::mtk::Interpolation_Mesh& tInterpolationMesh = mXTKModelPtr->get_background_mesh();

        mCoarseBasisToFineBasisWeights.resize( mCoarseBasisToFineBasis.size() );

        for ( uint Ik = 0; Ik < mCoarseBasisToFineBasis.size(); Ik++ )
        {
            mCoarseBasisToFineBasisWeights( Ik ).set_size( mCoarseBasisToFineBasis( Ik ).numel(), 1, MORIS_REAL_MAX );

            moris_index tBackgroundIndex = mEnrichedBasisToBackgroundBasis( Ik );

            moris::Matrix< DDRMat > tWeights = tInterpolationMesh.get_fine_basis_weights_of_basis( 0, tBackgroundIndex );

            moris::Matrix< DDSMat > tBGBasis = tInterpolationMesh.get_fine_basis_inds_of_basis( 0, tBackgroundIndex );

            moris::map< moris_index, sint > tBasisToPositionMap;
            for ( uint Ii = 0; Ii < tBGBasis.numel(); Ii++ )
            {
                tBasisToPositionMap[ tBGBasis( Ii ) ] = Ii;
            }

            for ( uint Ii = 0; Ii < mCoarseBasisToFineBasis( Ik ).numel(); Ii++ )
            {
                moris_index tFineBackgroundIndex = mEnrichedBasisToBackgroundBasis( mCoarseBasisToFineBasis( Ik )( Ii ) );

                sint tPos = tBasisToPositionMap.find( tFineBackgroundIndex );

                mCoarseBasisToFineBasisWeights( Ik )( Ii ) = tWeights( tPos );
            }
        }
        //        print( mCoarseBasisToFineBasisWeights,"mCoarseBasisToFineBasisWeights");
    }

    //------------------------------------------------------------------------------

    void
    Multigrid::build_enriched_coeff_to_background_coeff_map()
    {
        // get num enriched basis
        mNumBasis = mXTKModelPtr->mEnrichedInterpMesh( 0 )->get_max_num_coeffs_on_proc( 0 );

        // set size
        mEnrichedBasisToBackgroundBasis.resize( mNumBasis, -1 );

        // get background basis to enriched basis list. ( name of get function is misleading )
        const Vector< Matrix< IndexMat > >& tBackgroundCoeffsToEnrichedCoeffs = mXTKModelPtr->mEnrichedInterpMesh( 0 )
                                                                                      ->get_enriched_coefficients_to_background_coefficients( mMeshIndex );

        for ( uint Ik = 0; Ik < tBackgroundCoeffsToEnrichedCoeffs.size(); Ik++ )
        {
            for ( uint Ii = 0; Ii < tBackgroundCoeffsToEnrichedCoeffs( Ik ).numel(); Ii++ )
            {
                MORIS_ASSERT( mEnrichedBasisToBackgroundBasis( tBackgroundCoeffsToEnrichedCoeffs( Ik )( Ii ) ) == -1,
                        " Multigrid::build_enriched_coeff_to_background_coeff_map(), Enriched Basis index for two background basis indices" );

                mEnrichedBasisToBackgroundBasis( tBackgroundCoeffsToEnrichedCoeffs( Ik )( Ii ) ) = Ik;
            }
        }
    }

    //------------------------------------------------------------------------------
    // #ifdef MORIS_HAVE_DEBUG
    //    void Multigrid::save_to_vtk( const std::string & aFilePath )
    //    {
    //        // start timer
    //        tic tTimer;
    //
    //        // modify filename
    //        std::string tFilePath =  parallelize_path( aFilePath );
    //
    //        // open the file
    //        std::ofstream tFile( tFilePath, std::ios::binary );
    //
    //        // containers
    //        float tFChar = 0;
    //        int   tIChar = 0;
    //
    //        tFile << "# vtk DataFile Version 3.0" << std::endl;
    //        tFile << "GO BUFFS!" << std::endl;
    //        tFile << "BINARY" << std::endl;
    //
    //       // get my rank
    ////       moris_id tMyRank = par_rank();
    //
    //       uint tNumberOfBasis = mEnrichedBasisCoords.n_rows();
    //
    //       // write node data
    //       tFile << "DATASET UNSTRUCTURED_GRID" << std::endl;
    //       tFile << "POINTS " << tNumberOfBasis << " float"  << std::endl;
    //
    //       // ask settings for number of dimensions
    //       auto tNumberOfDimensions =   mXTKModelPtr->get_spatial_dim();
    //
    //       if ( tNumberOfDimensions == 2 )
    //       {
    //           for( auto tCoords : mEnrichedBasisCoords )
    //           {
    //               // write coordinates to mesh
    //               tFChar = swap_byte_endian( (float) tCoords.data()[ 0 ] );
    //               tFile.write( (char*) &tFChar, sizeof(float));
    //               tFChar = swap_byte_endian( (float) tCoords.data()[ 1 ] );
    //               tFile.write( (char*) &tFChar, sizeof(float));
    //               tFChar = swap_byte_endian( (float) 0 );
    //               tFile.write( (char*) &tFChar, sizeof(float));
    //           }
    //       }
    //       else if ( tNumberOfDimensions == 3 )
    //       {
    //           for( auto tCoords : mEnrichedBasisCoords )
    //           {
    //               // write coordinates to mesh
    //               tFChar = swap_byte_endian( (float) tCoords.data()[ 0 ] );
    //               tFile.write( (char*) &tFChar, sizeof(float));
    //               tFChar = swap_byte_endian( (float) tCoords.data()[ 1 ] );
    //               tFile.write( (char*) &tFChar, sizeof(float));
    //               tFChar = swap_byte_endian( (float) tCoords.data()[ 2 ] );
    //               tFile.write( (char*) &tFChar, sizeof(float));
    //           }
    //       }
    //
    //       tFile << std::endl;
    //
    //       // write each basis as its own element
    //       tFile << "CELLS " << tNumberOfBasis << " " << 2*tNumberOfBasis << std::endl;
    //
    //       int tOne = swap_byte_endian( (int) 1 );
    //
    //       // reset counter
    //       int tCount = 0;
    //
    //       for( auto tCoords : mEnrichedBasisCoords )
    //       {
    //               tIChar = swap_byte_endian( tCount );
    //               tFile.write( ( char* ) &tOne, sizeof(int));
    //               tFile.write( ( char *) &tIChar, sizeof(int));
    //
    //               ++tCount;
    //       }
    //
    //       // write cell types
    //       tFile << "CELL_TYPES " << tNumberOfBasis << std::endl;
    //       tIChar = swap_byte_endian( (int) 2 );
    //       for ( luint k = 0; k < tNumberOfBasis; ++k)
    //       {
    //           tFile.write( (char*) &tIChar, sizeof( int ) );
    //       }
    //
    //       // write node data
    //       tFile << "POINT_DATA " << tNumberOfBasis << std::endl;
    //
    //       // write state
    //       tFile << "SCALARS STATE int" << std::endl;
    //       tFile << "LOOKUP_TABLE default" << std::endl;
    //       for( auto tState : mEnrichedBasisStatus )
    //       {
    //               tIChar = swap_byte_endian( (int)tState );
    //
    //               tFile.write( ( char *) &tIChar, sizeof(int));
    //       }
    //       tFile << std::endl;
    //
    //       // write basis ID
    //       tFile << "SCALARS Ind int" << std::endl;
    //       tFile << "LOOKUP_TABLE default" << std::endl;
    //       for ( uint Ik = 0; Ik < tNumberOfBasis; ++Ik)
    //       {
    //               tIChar = swap_byte_endian( (int) Ik );
    //
    //               tFile.write( ( char *) &tIChar, sizeof(int));
    //       }
    //       tFile << std::endl;
    //
    //       // write basis level
    //       tFile << "SCALARS LEVEL int" << std::endl;
    //       tFile << "LOOKUP_TABLE default" << std::endl;
    //       for( auto tLevel : mEnrichedBasisLevel )
    //       {
    //               tIChar = swap_byte_endian( (int) tLevel );
    //
    //               tFile.write( ( char *) &tIChar, sizeof(int));
    //       }
    //       tFile << std::endl;
    //
    ////       // write basis owner
    ////       tFile << "SCALARS OWNER int" << std::endl;
    ////       tFile << "LOOKUP_TABLE default" << std::endl;
    ////       for( auto tOwner : mEnrichedBasisOwner )
    ////       {
    ////           tIChar = swap_byte_endian( (int) tBasis->get_owner() );
    ////
    ////           tFile.write( ( char *) &tIChar, sizeof(int));
    ////       }
    ////       tFile << std::endl;
    //
    //       // close the output file
    //       tFile.close();
    //
    //       // stop timer
    //       real tElapsedTime = tTimer.toc<moris::chronos::milliseconds>().wall;
    //
    //       // print output
    //       MORIS_LOG_INFO( "%s Created VTK debug file.\n               Mesh has %lu basis.\n               Creation took %5.3f seconds.\n\n",
    //               proc_string().c_str(),
    //               ( long unsigned int ) tNumberOfBasis,
    //               ( double ) tElapsedTime / 1000 );
    //    }
    // #endif

    void
    Multigrid::build_basis_exodus_information( std::string aName )
    {
        moris::mtk::Interpolation_Mesh& tInterpolationMesh = mXTKModelPtr->get_background_mesh();

        // get num enriched basis
        uint tNumEnrichedBasis = mXTKModelPtr->mEnrichedInterpMesh( 0 )->get_max_num_coeffs_on_proc( 0 );

#ifdef MORIS_HAVE_DEBUG
        mEnrichedBasisCoords.set_size( tNumEnrichedBasis, mXTKModelPtr->get_spatial_dim() );
        mEnrichedBasisStatus.set_size( tNumEnrichedBasis, 1 );
#endif
        mEnrichedBasisLevel.set_size( tNumEnrichedBasis, 1 );

        for ( uint Ik = 0; Ik < mNumBasis; Ik++ )
        {
            moris_index tBackgroundIndex = mEnrichedBasisToBackgroundBasis( Ik );

            mEnrichedBasisLevel( Ik ) = tInterpolationMesh.get_basis_level( 0, tBackgroundIndex );
#ifdef MORIS_HAVE_DEBUG
            mEnrichedBasisStatus( Ik )                                                     = tInterpolationMesh.get_basis_status( 0, tBackgroundIndex );
            mEnrichedBasisCoords( { Ik, Ik }, { 0, mXTKModelPtr->get_spatial_dim() - 1 } ) = tInterpolationMesh.get_basis_coords( 0, tBackgroundIndex ).matrix_data();
#endif
        }

#ifdef MORIS_HAVE_DEBUG
        // Create writer/file
        moris::mtk::Writer_Exodus tWriter;
        std::string               tMorisRoot = std::getenv( "MORISOUTPUT" );

        MORIS_ERROR( tMorisRoot.size() > 0,
                "Environment variable MORISOUTPUT not set." );

        std::string tTempName = "multigrid_basis_temp.exo";
        tWriter.write_points( tMorisRoot, aName, tMorisRoot, tTempName, mEnrichedBasisCoords );

        // Create fields
        moris::Vector< std::string > tPointFieldNames( 2 );
        tPointFieldNames( 0 ) = "Enriched Basis Level";
        tPointFieldNames( 1 ) = "Enriched Basis Status";
        tWriter.set_point_fields( tPointFieldNames );

        // Write the fields
        tWriter.set_time( 0.0 );
        tWriter.write_point_field( "Enriched Basis Level", mEnrichedBasisLevel );
        tWriter.write_point_field( "Enriched Basis Status", mEnrichedBasisStatus );

        // Close file
        tWriter.close_file();
#endif
        //        print( mEnrichedBasisCoords,"mEnrichedBasisCoords");
    }
}    // namespace xtk
