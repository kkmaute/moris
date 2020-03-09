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

namespace xtk {

    Multigrid::Multigrid(xtk::Model *aXTKModelPtr) : mXTKModelPtr(aXTKModelPtr) {
//        mXTKModelPtr
    }

//------------------------------------------------------------------------------

    void Multigrid::create_fine_to_coarse_relationship() {
        // get num enriched basis
        uint tNumEnrichedBasis = mXTKModelPtr->mEnrichedInterpMesh(0)->get_num_coeffs(0);

        // set size of fine to coarse list
        mFineBasisToCoarseBasis.resize(tNumEnrichedBasis);

        moris::mtk::Interpolation_Mesh &tInterpolationMesh = mXTKModelPtr->get_background_mesh().get_mesh_data();

        // get num bg basis
        uint tNumBGBasis = tInterpolationMesh.get_num_basis(0);

        std::cout << tNumBGBasis << " tNumBGBasis" << std::endl;

        // loop over bg basis
        for (uint Ik = 0; Ik < tNumBGBasis; Ik++) {
            // Basis indices are consecutive and correpond to Ik

            // get enriched basis for this bg basis
            const Matrix <IndexMat> &tEnrichedCoeffsForBackroundCoeffs = mXTKModelPtr->mEnrichedInterpMesh(0)
                    ->get_enriched_coefficients_at_background_coefficient(Ik);

            // loop over enriched basis for background basis Ik
            for (uint Ia = 0; Ia < tEnrichedCoeffsForBackroundCoeffs.numel(); Ia++) {
                // get num coarse basis for this basis
                uint tNumCoarseBasis = tInterpolationMesh.get_num_coarse_basis_of_basis(0, Ik);

                // get basis Ia
                moris_index tEnrichedBasisInd = tEnrichedCoeffsForBackroundCoeffs(Ia);

                // resize for coarse basis for this fine basis
                mFineBasisToCoarseBasis(tEnrichedBasisInd).set_size(tNumCoarseBasis, 1, -1);

                // get subphases for enriched basis
                const Cell <moris::Matrix<moris::IndexMat>> &tSubphaseIndForEnrichedBasis = mXTKModelPtr->mEnrichment
                        ->get_subphases_loc_inds_in_enriched_basis();

                // get FIRST sub-phase index of basis. First because we assume the fine basis support is complete within the coarse one
                moris_index tFirstSubphaseInSupportIndex = tSubphaseIndForEnrichedBasis(tEnrichedBasisInd)(0);

                // check if subphase is in child mesh. If it is continue with child mesh data. Otherwise with enriched mesh data
                if (mXTKModelPtr->subphase_is_in_child_mesh(tFirstSubphaseInSupportIndex)) {
                    // get list with child meshes in cut mesh. cut mesh is collection of child meshes. subphase corresponds to child mesh entry
                    const Matrix <IndexMat> &tSubPhaseIndexToChildMesh = mXTKModelPtr->get_cut_mesh().get_subphase_to_child_mesh_connectivity();

                    // get child mesh index
                    moris_index tChildMeshIndex = tSubPhaseIndexToChildMesh(tFirstSubphaseInSupportIndex);

                    // get child mesh using child mesh index
                    Child_Mesh &tCM = mXTKModelPtr->get_cut_mesh().get_child_mesh(tChildMeshIndex);

                    // get local index of this subphase in the child mesh
                    moris_index tSubphaseLocIndex = tCM.get_subphase_loc_index(tFirstSubphaseInSupportIndex);

                    // Subphase basis of the current subphase in child mesh
                    Cell <moris_index> const &tSubPhaseBasis = tCM.get_subphase_basis_indices(tSubphaseLocIndex);

                    Cell <moris_index> const &tSubphaseBasisEnrichLev = tCM.get_subphase_basis_enrichment_levels(
                            tSubphaseLocIndex);

                    // construct a map between bg basis index and index relative to the subphase cluster
                    std::unordered_map <moris_id, moris_id> tSubPhaseBasisMap = mXTKModelPtr->mEnrichment->construct_subphase_basis_to_basis_map(
                            tSubPhaseBasis);

                    // loop over coarse background basis
                    for (uint Ii = 0; Ii < tNumCoarseBasis; Ii++) {
                        // get coarse basis iondex
                        moris_id tCoarseBasisIndex = tInterpolationMesh.get_coarse_basis_index_of_basis(0, Ik, Ii);

                        // find enriched coarse basis index for this coarse bg basis index
                        moris_id tEnrichedBasisLocation = tSubPhaseBasisMap.find(tCoarseBasisIndex)->second;

                        moris_index tEnrichmentLevel = tSubphaseBasisEnrichLev(tEnrichedBasisLocation);

                        // get enriched basis for this coarse bg basis and enrichment level
                        const Matrix <IndexMat> &tCoarseEnrichedCoeffsForCoarseBackroundCoeffs = mXTKModelPtr->mEnrichedInterpMesh(
                                        0)
                                ->get_enriched_coefficients_at_background_coefficient(tCoarseBasisIndex);

                        moris_index tEnrichedCoarseBasisIndex = tCoarseEnrichedCoeffsForCoarseBackroundCoeffs(
                                tEnrichmentLevel);

                        // add enriched coarse basis index to list
                        mFineBasisToCoarseBasis(tEnrichedBasisInd)(Ii) = tEnrichedCoarseBasisIndex;

//                        std::cout<<tCoarseBasisIndex<<" tCoarseBasisIndex"<<std::endl;
//                        std::cout<<tEnrichedBasisIndex<<" tEnrichedBasisIndex"<<std::endl;
                    }
                } else {
                    // get bg basis interpolating into tFirstSubphaseInSupportIndex ( Interpolation cell index correponds to subphase index)
                    Cell <moris_index> tBasisForSubphaseIndex = mXTKModelPtr->mEnrichment->mInterpCellBasis(
                            tFirstSubphaseInSupportIndex);

                    // get enrichment level for bg basis interpolating into tFirstSubphaseInSupportIndex ( Interpolation cell index correponds to subphase index)
                    Cell <moris_index> tBasisEnrLevForSubphaseIndex = mXTKModelPtr->mEnrichment->mInterpCellBasisEnrLev(
                            tFirstSubphaseInSupportIndex);

                    // build map which maps bg basis index to entry in tBasisForSubphaseIndex/tBasisEnrLevForSubphaseIndex
                    std::unordered_map <moris_id, moris_id> tSubPhaseBasisMap = mXTKModelPtr->mEnrichment->construct_subphase_basis_to_basis_map(
                            tBasisForSubphaseIndex);

                    // loop over coarse basis
                    for (uint Ii = 0; Ii < tNumCoarseBasis; Ii++) {
                        // get coarse bg basis
                        moris_index tCoarseBasisIndex = tInterpolationMesh.get_coarse_basis_index_of_basis(0, Ik, Ii);

                        // find enrichment level for this coarse bg basis and subphase
                        moris_index tEnrichmentLev = tBasisEnrLevForSubphaseIndex(
                                tSubPhaseBasisMap.find(tCoarseBasisIndex)->second);

                        // get enriched basis for this coarse bg basis and enrichment level
                        const Matrix <IndexMat> &tCoarseEnrichedCoeffsForCoarseBackroundCoeffs = mXTKModelPtr->mEnrichedInterpMesh(
                                        0)
                                ->get_enriched_coefficients_at_background_coefficient(tCoarseBasisIndex);

                        moris_index tEnrichedCoarseBasisIndex = tCoarseEnrichedCoeffsForCoarseBackroundCoeffs(
                                tEnrichmentLev);

                        // add enriched coarse basis to list
                        mFineBasisToCoarseBasis(tEnrichedBasisInd)(Ii) = tEnrichedCoarseBasisIndex;
                    }
                }
            }
        }

        print(mFineBasisToCoarseBasis, "mFineBasisToCoarseBasis");
    }

//------------------------------------------------------------------------------

    void Multigrid::create_coarse_to_fine_relationship() {
//        // get num enriched basis
//        uint tNumEnrichedBasis = mXTKModelPtr->mEnrichedInterpMesh( 0 )->get_num_coeffs( 0 );
//
//        // set size of coarse to fine list
//        mCoarseBasisToFineBasis.resize( tNumEnrichedBasis );
//
//        moris::mtk::Interpolation_Mesh & tInterpolationMesh = mXTKModelPtr->get_background_mesh().get_mesh_data();
//
//        // get num bg basis
//        uint tNumBGBasis = tInterpolationMesh.get_num_basis( 0 );

    }

//------------------------------------------------------------------------------

    void Multigrid::build_enriched_coeff_to_background_coeff_map() {
        // get num enriched basis
        uint tNumEnrichedBasis = mXTKModelPtr->mEnrichedInterpMesh(0)->get_num_coeffs(0);

        std::cout << tNumEnrichedBasis << " tNumEnrichedBasis " << std::endl;

        // set size
        mEnrichedBasisToBackgroundBasis.resize(tNumEnrichedBasis, -1);

        // get background basis to enriched basis list. ( name of get function is misleading )
        const Cell <Matrix<IndexMat>> &tBackgroundCoeffsToEnrichedCoeffs = mXTKModelPtr->mEnrichedInterpMesh(0)
                ->get_enriched_coefficients_to_background_coefficients();

        for (uint Ik = 0; Ik < tBackgroundCoeffsToEnrichedCoeffs.size(); Ik++) {
            for (uint Ii = 0; Ii < tBackgroundCoeffsToEnrichedCoeffs(Ik).numel(); Ii++) {
                MORIS_ASSERT(mEnrichedBasisToBackgroundBasis(tBackgroundCoeffsToEnrichedCoeffs(Ik)(Ii)) == -1,
                             " Multigrid::build_enriched_coeff_to_background_coeff_map(), Enriched Basis index for two background basis indices");

                mEnrichedBasisToBackgroundBasis(tBackgroundCoeffsToEnrichedCoeffs(Ik)(Ii)) = Ik;
            }
        }
    }


//------------------------------------------------------------------------------
}
#ifdef DEBUG
#include "cl_MTK_Writer_Exodus.hpp"

    void xtk::Multigrid::build_basis_exodus_information()
    {
        moris::mtk::Interpolation_Mesh & tInterpolationMesh = mXTKModelPtr->get_background_mesh().get_mesh_data();

        // get num enriched basis
        uint tNumEnrichedBasis = mXTKModelPtr->mEnrichedInterpMesh( 0 )->get_num_coeffs( 0 );

        mEnrichedBasisCoords.set_size( tNumEnrichedBasis, mXTKModelPtr->get_spatial_dim() );
        mEnrichedBasisLevel .set_size( tNumEnrichedBasis, 1 );
        mEnrichedBasisStatus.set_size( tNumEnrichedBasis, 1 );

        for( uint Ik = 0; Ik < tNumEnrichedBasis; Ik++ )
        {
            moris_index tBackgroundIndex = mEnrichedBasisToBackgroundBasis( Ik );

            mEnrichedBasisLevel ( Ik ) = tInterpolationMesh.get_basis_level ( 0, tBackgroundIndex);
            mEnrichedBasisStatus( Ik ) = tInterpolationMesh.get_basis_status( 0, tBackgroundIndex);
            mEnrichedBasisCoords( {Ik, Ik}, {0, mXTKModelPtr->get_spatial_dim() - 1} ) = tInterpolationMesh.get_basis_coords( 0, tBackgroundIndex).matrix_data();
        }

        // Create writer/file
        moris::mtk::Writer_Exodus tWriter;
        std::string tMorisRoot = std::getenv("MORISOUTPUT");
        tWriter.write_points(tMorisRoot, "test.exo", mEnrichedBasisCoords);

        // Create fields
        moris::Cell<std::string> tPointFieldNames(2);
        tPointFieldNames(0) = "Enriched Basis Level";
        tPointFieldNames(1) = "Enriched Basis Status";
        tWriter.set_point_fields(tPointFieldNames);

        // Write the fields
        tWriter.set_time(0.0);
        tWriter.write_point_field("Enriched Basis Level", mEnrichedBasisLevel);
        tWriter.write_point_field("Enriched Basis Status", mEnrichedBasisStatus);

        // Close file
        tWriter.close_file();

//        print( mEnrichedBasisCoords,"mEnrichedBasisCoords");
    }

#endif


