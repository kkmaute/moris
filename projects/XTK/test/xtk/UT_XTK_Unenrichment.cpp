/*
 * Copyright (c) 2022 University of Colorado
 *Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_XTK_Unenrichment.cpp
 *
 */


#include "catch.hpp"
#include "paths.hpp"
#include "cl_Matrix.hpp"

// implementations to test
#include "cl_MTK_Mesh_Factory.hpp"

#define protected public
#define private public
#include "cl_XTK_Model.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
#undef protected
#undef private

#include "cl_XTK_Enriched_Integration_Mesh.hpp"

#include "cl_HMR_Mesh_Interpolation.hpp"
#include "cl_HMR.hpp"
#include "cl_GEN_Circle.hpp"
#include "cl_GEN_Sphere.hpp"
#include "cl_GEN_Plane.hpp"
#include "fn_PRM_HMR_Parameters.hpp"
#include "fn_PRM_XTK_Parameters.hpp"
#include "cl_MTK_Intersection_Detect.hpp"
#include "cl_MTK_Writer_Exodus.hpp"

#include "cl_MTK_Intersection_Mesh.hpp"

namespace xtk
{
    TEST_CASE( "XTK UnEnrichment", "[XTK],[XTK_UnEnrichment]" )
    {
        // the setup of the problem is : a 2*1 mesh with 2*1 elements where the right element is cut, so enriched and unzipped.
        // the total number of basis should be 10 after enrichment but after unenrichment it should be 12 again.

        if ( par_size() == 1 or par_size() == 2 )
        {
            uint tLagrangeMeshIndex = 0;

            // HMR parameter list
            moris::ParameterList tParameters = moris::prm::create_hmr_parameter_list();

            tParameters.set( "number_of_elements_per_dimension", "2, 1" );
            tParameters.set( "domain_dimensions", "2, 1" );
            tParameters.set( "domain_offset", "0.0, 0.0" );
            tParameters.set( "domain_sidesets", "1,2,3,4" );
            tParameters.set( "lagrange_output_meshes", "0" );

            tParameters.set( "lagrange_orders", "2" );
            tParameters.set( "lagrange_pattern", "0" );
            tParameters.set( "bspline_orders", "2" );
            tParameters.set( "bspline_pattern", "0" );

            tParameters.set( "lagrange_to_bspline", "0" );

            tParameters.set( "truncate_bsplines", 1 );
            tParameters.set( "refinement_buffer", 1 );
            tParameters.set( "staircase_buffer", 1 );
            tParameters.set( "initial_refinement", "0" );
            tParameters.set( "initial_refinement_pattern", "0" );

            tParameters.set( "use_multigrid", 0 );
            tParameters.set( "severity_level", 0 );

            // create the hmr mesh
            hmr::HMR tHMR( tParameters );

            tHMR.perform_initial_refinement();

            tHMR.finalize();

            moris::hmr::Interpolation_Mesh_HMR* tInterpolationMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );

            // xtk parameter list
            moris::ParameterList tXTKParameters = moris::prm::create_xtk_parameter_list();
            tXTKParameters.set( "decompose", true );
            tXTKParameters.set( "decomposition_type", "conformal" );
            tXTKParameters.set( "enrich", true );
            tXTKParameters.set( "basis_rank", "bspline" );
            tXTKParameters.set( "enrich_mesh_indices", "0" );
            tXTKParameters.set( "ghost_stab", true );
            tXTKParameters.set( "multigrid", false );
            tXTKParameters.set( "verbose", true );
            tXTKParameters.set( "print_enriched_ig_mesh", false );
            tXTKParameters.set( "exodus_output_XTK_ig_mesh", false );
            tXTKParameters.set( "high_to_low_dbl_side_sets", true );
            tXTKParameters.set( "write_cell_enrichments_levels", false );

            // define the sphere such that it is non interacting
            Vector< std::shared_ptr< moris::ge::Geometry > > tGeometry( 1 );
            tGeometry( 0 ) = std::make_shared< moris::ge::Plane >( 1.5, 0.5, 1.0, 0.0 );

            // define ge engine
            moris::ge::Geometry_Engine_Parameters tGeometryEngineParameters;
            tGeometryEngineParameters.mGeometries = tGeometry;
            moris::ge::Geometry_Engine tGeometryEngine( tInterpolationMesh, tGeometryEngineParameters );

            // dimension of the
            size_t tModelDimension = 2;

            Model tXTKModel( tModelDimension, tInterpolationMesh, &tGeometryEngine );
            tXTKModel.mVerbose = false;

            // Specify decomposition Method and Cut Mesh ---------------------------------------
            Vector< enum Subdivision_Method > tDecompositionMethods = { Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3 };
            tXTKModel.decompose( tDecompositionMethods );

            // perform basis enrichment
            tXTKModel.perform_basis_enrichment( mtk::EntityRank::BSPLINE, 0 );

            // perform unenrichment
            tXTKModel.perform_unenrichment( { { 0 } } );

            // construct ghost
            tXTKModel.construct_face_oriented_ghost_penalization_cells();

            // get the enriched interpolation mesh
            xtk::Enriched_Interpolation_Mesh& tInterpMesh = tXTKModel.get_enriched_interp_mesh();

            // mesh indices
            moris::Matrix< moris::IndexMat > tMeshIndices( 1, 1, 0 );

            // check if the unenrichment mesh indices are set
            REQUIRE( tInterpMesh.mUnenrichedMeshIndices( 0 ) == tMeshIndices( 0 ) );

            // get the adof map of the non-enriched
            moris::map< moris::moris_id, moris::moris_index > tUnenrichedMap;
            tInterpMesh.get_adof_map( tMeshIndices( 0 ), tUnenrichedMap );

            // get the hmr map which is basically
            moris::map< moris::moris_id, moris::moris_index > tHMRMap;
            tXTKModel.mBackgroundMesh->get_adof_map( tMeshIndices( 0 ), tHMRMap );

            // loop over the vertex enrichments to see if the t-matrics are actually overwritten
            for ( Vertex_Enrichment* iVertexEnrichment : tInterpMesh.mInterpVertEnrichment( tInterpMesh.mUnenrichedMeshIndices( 0 ) ) )
            {
                // if it is not nullptr
                if ( iVertexEnrichment->has_interpolation() )
                {
                    mtk::Vertex_Interpolation const * tBaseVertex = iVertexEnrichment->get_base_vertex_interpolation();
                    if ( tBaseVertex )
                    {
                        moris::Matrix< moris::IndexMat >         tBaseIndices = tBaseVertex->get_indices();
                        moris::Matrix< moris::IndexMat > const & tIndices     = iVertexEnrichment->get_basis_indices();

                        moris::Matrix< moris::IdMat >         tBaseIds = tBaseVertex->get_ids();
                        moris::Matrix< moris::IdMat > const & tIds     = iVertexEnrichment->get_basis_ids();


                        moris::Matrix< moris::IdMat >         tBaseOwners = tBaseVertex->get_owners();
                        moris::Matrix< moris::IdMat > const & tOwners     = iVertexEnrichment->get_owners();


                        bool tSameIndex = std::equal( tBaseIndices.begin(), tBaseIndices.end(), tIndices.begin(),    //
                                []( moris_index aBaseIndex, moris_index aIndex ) -> bool { return aBaseIndex == aIndex; } );

                        bool tSameId = std::equal( tBaseIds.begin(), tBaseIds.end(), tIds.begin(),    //
                                []( moris_id aBaseId, moris_id aId ) -> bool { return aBaseId == aId; } );


                        bool tSameOwner = std::equal( tBaseOwners.begin(), tBaseOwners.end(), tOwners.begin(),    //
                                []( moris_id aBaseOwner, moris_id aOwner ) -> bool { return aBaseOwner == aOwner; } );


                        bool tItisOverwritten = tSameIndex and tSameId and tSameOwner;

                        CHECK( tItisOverwritten );
                    }

                    // if it has interpolation data but does not have a base vertex them it is created during the ghost and it is in auro
                    // we will check that ids created during the ghost
                    else
                    {
                        moris::Matrix< moris::IndexMat > const & tIndices = iVertexEnrichment->get_basis_indices();
                        moris::Matrix< moris::IdMat > const &    tIds     = iVertexEnrichment->get_basis_ids();
                        moris::Matrix< moris::IdMat > const &    tOwners  = iVertexEnrichment->get_owners();

                        // if there is a vertex that belongs to the other processor, it should be added to the map and
                        // the indices and id's should not exist in the map beforehand
                        for ( uint iCounter = 0; iCounter < tOwners.numel(); iCounter++ )
                        {
                            // if it is not owned by the processor and it does not exist in hmr map
                            // the second check is to ensure we are not overwriting thr map
                            if ( tOwners( iCounter ) != par_rank() and !tHMRMap.key_exists( tIds( iCounter ) ) )
                            {
                                // add the pair to the map
                                tHMRMap[ tIds( iCounter ) ] = tIndices( iCounter );
                            }
                        }
                    }
                }
            }

            // initialize the max index that exist
            moris::moris_index tMaxIndex = -1;

            // loop over the maps and check fi the map created by xtk and modified hmr map are the same
            for ( const auto& iHMRGlobalToLocal : tHMRMap )
            {
                // get the indices from two maps
                moris_index tHMRIndex         = iHMRGlobalToLocal.second;
                moris_index tNonEnirhcedIndex = tUnenrichedMap.find( iHMRGlobalToLocal.first );

                // compare the indices to be the same
                CHECK( tHMRIndex == tNonEnirhcedIndex );

                // update the max index
                tMaxIndex = std::max( tMaxIndex, tHMRIndex );
            }

            // this ensures that the indices assigned during the ghost are continuation of the HMR index
            CHECK( tInterpMesh.get_max_num_coeffs_on_proc( tMeshIndices( 0 ) ) == (uint)( tMaxIndex + 1 ) );
        }
    }
}    // namespace xtk