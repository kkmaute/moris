/*
 * cl_XTK_Multigrid.cpp
 *
 *  Created on: Feb 18, 2019
 *      Author: Schmidt
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
#include "typedefs.hpp"
#include "cl_MTK_Cell_Info_Factory.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh_Core.hpp"
#include "xtk_typedefs.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
#include "cl_XTK_Model.hpp"

namespace xtk
{

    Multigrid::Multigrid( xtk::Model * aXTKModelPtr ) : mXTKModelPtr( aXTKModelPtr )
    {

    }

//------------------------------------------------------------------------------

    void Multigrid::create_fine_to_coarse_relationship()
    {
        // set size of fine to coarse list
        mFineBasisToCoarseBasis.resize( mNumBasis );

        moris::mtk::Interpolation_Mesh & tInterpolationMesh = mXTKModelPtr->get_background_mesh().get_mesh_data();

        // get num bg basis
        uint tNumBGBasis = tInterpolationMesh.get_num_basis( 0 );

        // loop over bg basis
        for ( uint Ik = 0; Ik < tNumBGBasis; Ik ++ )
        {
            // Basis indices are consecutive and correpond to Ik

            // get enriched basis for this bg basis
            const Matrix<IndexMat> & tEnrichedCoeffsForBackroundCoeffs = mXTKModelPtr->mEnrichedInterpMesh( 0 )
                                                       ->get_enriched_coefficients_at_background_coefficient( Ik );

            // loop over enriched basis for background basis Ik
            for ( uint Ia = 0; Ia < tEnrichedCoeffsForBackroundCoeffs.numel(); Ia ++ )
            {
                // get num coarse basis for this basis
                uint tNumCoarseBasis = tInterpolationMesh.get_num_coarse_basis_of_basis( 0, Ik );

                // get basis Ia
                moris_index tEnrichedBasisInd = tEnrichedCoeffsForBackroundCoeffs( Ia );

                // resize for coarse basis for this fine basis
                mFineBasisToCoarseBasis( tEnrichedBasisInd ).set_size( tNumCoarseBasis, 1, -1 );

                // get subphases for enriched basis
                const Cell< moris::Matrix< moris::IndexMat > > & tSubphaseIndForEnrichedBasis = mXTKModelPtr->mEnrichment
                                                      ->get_subphases_loc_inds_in_enriched_basis();

                // get FIRST sub-phase index of basis. First because we assume the fine basis support is complete within the coarse one
                moris_index tFirstSubphaseInSupportIndex = tSubphaseIndForEnrichedBasis( tEnrichedBasisInd )( 0 );

                // check if subphase is in child mesh. If it is continue with child mesh data. Otherwise with enriched mesh data
                if( mXTKModelPtr->subphase_is_in_child_mesh( tFirstSubphaseInSupportIndex ) )
                {
                    // get list with child meshes in cut mesh. cut mesh is collection of child meshes. subphase corresponds to child mesh entry
                    const Matrix<IndexMat> & tSubPhaseIndexToChildMesh = mXTKModelPtr->get_cut_mesh().get_subphase_to_child_mesh_connectivity();

                    // get child mesh index
                    moris_index tChildMeshIndex = tSubPhaseIndexToChildMesh( tFirstSubphaseInSupportIndex );

                    // get child mesh using child mesh index
                    Child_Mesh & tCM = mXTKModelPtr->get_cut_mesh().get_child_mesh( tChildMeshIndex );

                    // get local index of this subphase in the child mesh
                    moris_index tSubphaseLocIndex = tCM.get_subphase_loc_index( tFirstSubphaseInSupportIndex );

                    // Subphase basis of the current subphase in child mesh
                    Cell< moris_index > const & tSubPhaseBasis = tCM.get_subphase_basis_indices( tSubphaseLocIndex );

                    Cell< moris_index > const & tSubphaseBasisEnrichLev = tCM.get_subphase_basis_enrichment_levels( tSubphaseLocIndex );

                    // construct a map between bg basis index and index relative to the subphase cluster
                    std::unordered_map< moris_id, moris_id > tSubPhaseBasisMap = mXTKModelPtr->mEnrichment->construct_subphase_basis_to_basis_map( tSubPhaseBasis );

                    // loop over coarse background basis
                    for ( uint Ii = 0; Ii < tNumCoarseBasis; Ii ++ )
                    {
                        // get coarse basis iondex
                        moris_id tCoarseBasisIndex = tInterpolationMesh.get_coarse_basis_index_of_basis( 0, Ik, Ii );

                        // find enriched coarse basis index for this coarse bg basis index
                        moris_id tEnrichedBasisLocation = tSubPhaseBasisMap.find( tCoarseBasisIndex )->second;

                        moris_index tEnrichmentLevel = tSubphaseBasisEnrichLev( tEnrichedBasisLocation );

                        // get enriched basis for this coarse bg basis and enrichment level
                        const Matrix<IndexMat> & tCoarseEnrichedCoeffsForCoarseBackroundCoeffs = mXTKModelPtr->mEnrichedInterpMesh( 0 )
                                                                   ->get_enriched_coefficients_at_background_coefficient( tCoarseBasisIndex );

                        moris_index tEnrichedCoarseBasisIndex = tCoarseEnrichedCoeffsForCoarseBackroundCoeffs( tEnrichmentLevel );

                        // add enriched coarse basis index to list
                        mFineBasisToCoarseBasis( tEnrichedBasisInd )( Ii ) = tEnrichedCoarseBasisIndex;

//                        std::cout<<tCoarseBasisIndex<<" tCoarseBasisIndex"<<std::endl;
//                        std::cout<<tEnrichedBasisIndex<<" tEnrichedBasisIndex"<<std::endl;
                    }
                }
                else
                {
                    // get bg basis interpolating into tFirstSubphaseInSupportIndex ( Interpolation cell index correponds to subphase index)
                    Cell< moris_index > tBasisForSubphaseIndex = mXTKModelPtr->mEnrichment->mInterpCellBasis( tFirstSubphaseInSupportIndex );

                    // get enrichment level for bg basis interpolating into tFirstSubphaseInSupportIndex ( Interpolation cell index correponds to subphase index)
                    Cell< moris_index > tBasisEnrLevForSubphaseIndex = mXTKModelPtr->mEnrichment->mInterpCellBasisEnrLev( tFirstSubphaseInSupportIndex );

                    // build map which maps bg basis index to entry in tBasisForSubphaseIndex/tBasisEnrLevForSubphaseIndex
                    std::unordered_map< moris_id, moris_id > tSubPhaseBasisMap = mXTKModelPtr->mEnrichment->construct_subphase_basis_to_basis_map( tBasisForSubphaseIndex );

                    // loop over coarse basis
                    for ( uint Ii = 0; Ii < tNumCoarseBasis; Ii ++ )
                    {
                        // get coarse bg basis
                        moris_index tCoarseBasisIndex = tInterpolationMesh.get_coarse_basis_index_of_basis( 0, Ik, Ii );

                        // find enrichment level for this coarse bg basis and subphase
                        moris_index tEnrichmentLev = tBasisEnrLevForSubphaseIndex( tSubPhaseBasisMap.find( tCoarseBasisIndex )->second );

                        // get enriched basis for this coarse bg basis and enrichment level
                        const Matrix<IndexMat> & tCoarseEnrichedCoeffsForCoarseBackroundCoeffs = mXTKModelPtr->mEnrichedInterpMesh( 0 )
                                                                   ->get_enriched_coefficients_at_background_coefficient( tCoarseBasisIndex );

                        moris_index tEnrichedCoarseBasisIndex = tCoarseEnrichedCoeffsForCoarseBackroundCoeffs( tEnrichmentLev );

                        // add enriched coarse basis to list
                        mFineBasisToCoarseBasis( tEnrichedBasisInd )( Ii ) = tEnrichedCoarseBasisIndex;
                    }
                }
            }
        }

//        print( mFineBasisToCoarseBasis,"mFineBasisToCoarseBasis");
    }

//------------------------------------------------------------------------------

    void Multigrid::create_coarse_to_fine_relationship()
    {
        // set size of fine to coarse list
        mCoarseBasisToFineBasis.resize( mNumBasis );

        moris::Cell< uint > tCounter( mNumBasis, 0 );

        for( uint Ik = 0; Ik < mFineBasisToCoarseBasis.size(); Ik++ )
        {
            for( uint Ii = 0; Ii < mFineBasisToCoarseBasis( Ik ).numel(); Ii++ )
            {
                tCounter( mFineBasisToCoarseBasis( Ik )( Ii ) )++;
            }
        }

        for( uint Ik = 0; Ik < mCoarseBasisToFineBasis.size(); Ik++ )
        {
            mCoarseBasisToFineBasis( Ik ).set_size( tCounter( Ik ), 1, 0 );

            tCounter( Ik ) = 0;
        }

        for( uint Ik = 0; Ik < mFineBasisToCoarseBasis.size(); Ik++ )
        {
            for( uint Ii = 0; Ii < mFineBasisToCoarseBasis( Ik ).numel(); Ii++ )
            {
                mCoarseBasisToFineBasis( mFineBasisToCoarseBasis( Ik )( Ii ) )( tCounter( mFineBasisToCoarseBasis( Ik )( Ii ) ) ) = Ik;

                tCounter( mFineBasisToCoarseBasis( Ik )( Ii ) )++;
            }
        }

//        print( mCoarseBasisToFineBasis,"mCoarseBasisToFineBasis");
    }

//------------------------------------------------------------------------------

    void Multigrid::create_coarse_to_fine_weights()
    {
        moris::mtk::Interpolation_Mesh & tInterpolationMesh = mXTKModelPtr->get_background_mesh().get_mesh_data();

        mCoarseBasisToFineBasisWeights.resize( mCoarseBasisToFineBasis.size() );

        for( uint Ik = 0; Ik < mCoarseBasisToFineBasis.size(); Ik++ )
        {
            mCoarseBasisToFineBasisWeights( Ik ).set_size( mCoarseBasisToFineBasis( Ik ).numel(), 1, MORIS_REAL_MAX );

            moris_index tBackgroundIndex = mEnrichedBasisToBackgroundBasis( Ik );

            moris::Matrix< DDRMat > tWeights = tInterpolationMesh.get_fine_basis_weights_of_basis ( 0, tBackgroundIndex );

            moris::Matrix< DDSMat > tBGBasis = tInterpolationMesh.get_fine_basis_inds_of_basis ( 0, tBackgroundIndex );

            moris::map< moris_index, sint > tBasisToPositionMap;
            for( uint Ii = 0; Ii < tBGBasis.numel(); Ii++ )
            {
                tBasisToPositionMap[ tBGBasis( Ii ) ] = Ii;
            }

            for( uint Ii = 0; Ii < mCoarseBasisToFineBasis( Ik ).numel(); Ii++ )
            {
                moris_index tFineBackgroundIndex = mEnrichedBasisToBackgroundBasis( mCoarseBasisToFineBasis( Ik )( Ii ) );

                sint tPos = tBasisToPositionMap.find( tFineBackgroundIndex );

                mCoarseBasisToFineBasisWeights( Ik )( Ii ) = tWeights( tPos );
            }
        }
//        print( mCoarseBasisToFineBasisWeights,"mCoarseBasisToFineBasisWeights");
    }

//------------------------------------------------------------------------------

    void Multigrid::build_enriched_coeff_to_background_coeff_map()
    {
        // get num enriched basis
    	mNumBasis = mXTKModelPtr->mEnrichedInterpMesh( 0 )->get_num_coeffs( 0 );

        // set size
        mEnrichedBasisToBackgroundBasis.resize( mNumBasis, -1 );

        // get background basis to enriched basis list. ( name of get function is misleading )
        const Cell<Matrix<IndexMat>> & tBackgroundCoeffsToEnrichedCoeffs = mXTKModelPtr->mEnrichedInterpMesh( 0 )
                                                   ->get_enriched_coefficients_to_background_coefficients();

        for( uint Ik = 0; Ik < tBackgroundCoeffsToEnrichedCoeffs.size(); Ik++ )
        {
            for( uint Ii = 0; Ii < tBackgroundCoeffsToEnrichedCoeffs( Ik ).numel(); Ii++ )
            {
                MORIS_ASSERT( mEnrichedBasisToBackgroundBasis( tBackgroundCoeffsToEnrichedCoeffs( Ik )( Ii ) ) == -1,
                        " Multigrid::build_enriched_coeff_to_background_coeff_map(), Enriched Basis index for two background basis indices");

                mEnrichedBasisToBackgroundBasis( tBackgroundCoeffsToEnrichedCoeffs( Ik )( Ii ) ) = Ik;
            }
        }
    }

//------------------------------------------------------------------------------
#ifdef DEBUG
    void Multigrid::save_to_vtk( const std::string & aFilePath )
    {
        // start timer
        tic tTimer;

        // modify filename
        std::string tFilePath =  parallelize_path( aFilePath );

        // open the file
        std::ofstream tFile( tFilePath, std::ios::binary );

        // containers
        float tFChar = 0;
        int   tIChar = 0;

        tFile << "# vtk DataFile Version 3.0" << std::endl;
        tFile << "GO BUFFS!" << std::endl;
        tFile << "BINARY" << std::endl;

       // get my rank
//       moris_id tMyRank = par_rank();

       uint tNumberOfBasis = mEnrichedBasisCoords.size();

       // write node data
       tFile << "DATASET UNSTRUCTURED_GRID" << std::endl;
       tFile << "POINTS " << tNumberOfBasis << " float"  << std::endl;

       // ask settings for numner of dimensions
       auto tNumberOfDimensions =   mXTKModelPtr->get_spatial_dim();

       if ( tNumberOfDimensions == 2 )
       {
           for( auto tCoords : mEnrichedBasisCoords )
           {
               // write coordinates to mesh
               tFChar = swap_byte_endian( (float) tCoords.data()[ 0 ] );
               tFile.write( (char*) &tFChar, sizeof(float));
               tFChar = swap_byte_endian( (float) tCoords.data()[ 1 ] );
               tFile.write( (char*) &tFChar, sizeof(float));
               tFChar = swap_byte_endian( (float) 0 );
               tFile.write( (char*) &tFChar, sizeof(float));
           }
       }
       else if ( tNumberOfDimensions == 3 )
       {
           for( auto tCoords : mEnrichedBasisCoords )
           {
               // write coordinates to mesh
               tFChar = swap_byte_endian( (float) tCoords.data()[ 0 ] );
               tFile.write( (char*) &tFChar, sizeof(float));
               tFChar = swap_byte_endian( (float) tCoords.data()[ 1 ] );
               tFile.write( (char*) &tFChar, sizeof(float));
               tFChar = swap_byte_endian( (float) tCoords.data()[ 2 ] );
               tFile.write( (char*) &tFChar, sizeof(float));
           }
       }

       tFile << std::endl;

       // write each basis as its own element
       tFile << "CELLS " << tNumberOfBasis << " " << 2*tNumberOfBasis << std::endl;

       int tOne = swap_byte_endian( (int) 1 );

       // reset counter
       int tCount = 0;

       for( auto tCoords : mEnrichedBasisCoords )
       {
               tIChar = swap_byte_endian( tCount );
               tFile.write( ( char* ) &tOne, sizeof(int));
               tFile.write( ( char *) &tIChar, sizeof(int));

               ++tCount;
       }

       // write cell types
       tFile << "CELL_TYPES " << tNumberOfBasis << std::endl;
       tIChar = swap_byte_endian( (int) 2 );
       for ( luint k = 0; k < tNumberOfBasis; ++k)
       {
           tFile.write( (char*) &tIChar, sizeof( int ) );
       }

       // write node data
       tFile << "POINT_DATA " << tNumberOfBasis << std::endl;

       // write state
       tFile << "SCALARS STATE int" << std::endl;
       tFile << "LOOKUP_TABLE default" << std::endl;
       for( auto tState : mEnrichedBasisStatus )
       {
               tIChar = swap_byte_endian( (int)tState );

               tFile.write( ( char *) &tIChar, sizeof(int));
       }
       tFile << std::endl;

       // write basis ID
       tFile << "SCALARS Ind int" << std::endl;
       tFile << "LOOKUP_TABLE default" << std::endl;
       for ( uint Ik = 0; Ik < tNumberOfBasis; ++Ik)
       {
               tIChar = swap_byte_endian( (int) Ik );

               tFile.write( ( char *) &tIChar, sizeof(int));
       }
       tFile << std::endl;

       // write basis level
       tFile << "SCALARS LEVEL int" << std::endl;
       tFile << "LOOKUP_TABLE default" << std::endl;
       for( auto tLevel : mEnrichedBasisLevel )
       {
               tIChar = swap_byte_endian( (int) tLevel );

               tFile.write( ( char *) &tIChar, sizeof(int));
       }
       tFile << std::endl;

//       // write basis owner
//       tFile << "SCALARS OWNER int" << std::endl;
//       tFile << "LOOKUP_TABLE default" << std::endl;
//       for( auto tOwner : mEnrichedBasisOwner )
//       {
//           tIChar = swap_byte_endian( (int) tBasis->get_owner() );
//
//           tFile.write( ( char *) &tIChar, sizeof(int));
//       }
//       tFile << std::endl;

       // close the output file
       tFile.close();

       // stop timer
       real tElapsedTime = tTimer.toc<moris::chronos::milliseconds>().wall;

       // print output
       MORIS_LOG_INFO( "%s Created VTK debug file.\n               Mesh has %lu basis.\n               Creation took %5.3f seconds.\n\n",
               proc_string().c_str(),
               ( long unsigned int ) tNumberOfBasis,
               ( double ) tElapsedTime / 1000 );
    }

//------------------------------------------------------------------------------
#endif
    void Multigrid::build_basis_exodus_information()
    {
        moris::mtk::Interpolation_Mesh & tInterpolationMesh = mXTKModelPtr->get_background_mesh().get_mesh_data();
#ifdef DEBUG
        mEnrichedBasisCoords.resize( mNumBasis );
        mEnrichedBasisStatus.set_size( mNumBasis, 1 );
#endif
        mEnrichedBasisLevel .set_size( mNumBasis, 1 );


        for( uint Ik = 0; Ik < mNumBasis; Ik++ )
        {
            moris_index tBackgroundIndex = mEnrichedBasisToBackgroundBasis( Ik );

            mEnrichedBasisLevel ( Ik ) = tInterpolationMesh.get_basis_level ( 0, tBackgroundIndex);
#ifdef DEBUG
            mEnrichedBasisStatus( Ik ) = tInterpolationMesh.get_basis_status( 0, tBackgroundIndex);
            mEnrichedBasisCoords( Ik ) = tInterpolationMesh.get_basis_coords( 0, tBackgroundIndex);
#endif
        }

//        print( mEnrichedBasisCoords,"mEnrichedBasisCoords");
    }

//#endif
}


