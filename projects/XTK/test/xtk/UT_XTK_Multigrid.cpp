/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_XTK_Multigrid.cpp
 *
 */

#include <memory>
#include <mpi.h>
#include "catch.hpp"

// XTKL: Mesh Includes
#include "cl_MTK_Mesh.hpp"
#include "fn_verify_tet_topology.hpp"
#include "fn_write_element_ownership_as_field.hpp"

// XTKL: Geometry  Include
#include "cl_Logger.hpp"

// XTKL: Container includes
#include "cl_Vector.hpp"

// XTKL: Linear Algebra Includes

#include "cl_Matrix.hpp"
#include "cl_XTK_Matrix_Base_Utilities.hpp"
#include "linalg_typedefs.hpp"

#include "geometry/cl_Mesh_Field_Geometry.hpp"
#include "geometry/cl_Plane.hpp"
#include "cl_MGE_Geometry_Engine.hpp"

#include "cl_XTK_Model.hpp"
#include "cl_XTK_Enums.hpp"
#include "cl_XTK_Cut_Mesh.hpp"
#include "cl_XTK_Enrichment.hpp"
#include "cl_MTK_Mesh_XTK_Impl.hpp"

#include "xtk_typedefs.hpp"
#include "geometry/cl_Geom_Field.hpp"

#include "cl_HMR_Parameters.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Field.hpp"

moris::real
LevelSetFunction( const moris::Matrix< DDRMat >& aPoint )
{
    // return norm( aPoint ) - 0.9;
    return aPoint( 0 ) - 10;
}

namespace moris::xtk
{
    TEST_CASE( "XTK Multigrid", "[xtk_multigrid]" )
    {
        if ( par_size() == 1 )
        {
            // order for this example
            moris::uint tOrder = 1;

            // create parameter object
            moris::hmr::Parameters tParameters;
            tParameters.set_number_of_elements_per_dimension( 1, 1, 1 );
            tParameters.set_verbose( false );
            tParameters.set_multigrid( true );
            tParameters.set_bspline_truncation( true );
            tParameters.set_refinement_buffer( 1 );

            // create HMR object
            moris::hmr::HMR tHMR( tParameters );

            uint tNumberOfElements = tHMR.get_database()->get_background_mesh()->get_number_of_active_elements_on_proc();

            // flag all elements
            for ( uint e = 0; e < tNumberOfElements; ++e )
            {
                tHMR.flag_element( e );
            }

            // flag first element for refinement
            // tHMR.flag_element( 0 );
            tHMR.perform_refinement_based_on_working_pattern( 0 );

            //            tNumberOfElements = tHMR.get_database()->get_background_mesh()->get_number_of_active_elements_on_proc();
            //
            //            // flag all elements
            //            for( uint e=0; e<tNumberOfElements; ++e )
            //            {
            //                tHMR.flag_element( e );
            //            }
            //
            //            //tHMR.flag_element( 0 );
            //            tHMR.perform_refinement_based_on_working_pattern( moris::hmr::RefinementMode::SIMPLE );
            //            tHMR.update_refinement_pattern();
            //
            //            tNumberOfElements = tHMR.get_database()->get_background_mesh()->get_number_of_active_elements_on_proc();
            //
            //            // flag all elements
            //            for( uint e=0; e<tNumberOfElements; ++e )
            //            {
            //                tHMR.flag_element( e );
            //            }
            //
            //            //tHMR.flag_element( 0 );
            //            tHMR.perform_refinement_based_on_working_pattern( moris::hmr::RefinementMode::SIMPLE );
            //            tHMR.update_refinement_pattern();

            tHMR.finalize();

            std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tOrder );

            // create field
            std::shared_ptr< moris::hmr::Field > tField = tMesh->create_field( "Circle", tOrder );

            // evaluate node values
            tField->evaluate_scalar_function( LevelSetFunction );

            tHMR.save_bsplines_to_vtk( "BSplines1.vtk" );

            //            tHMR.save_to_exodus( "xtk_test_sphere1.exo" );

            //-----------------------------------------------------------------------------------------------------

            for ( moris::uint i = 0; i < tMesh->get_num_entities( EntityRank::NODE ); i++ )
            {
                std::cout << i << " iterator" << std::endl;
                moris::print( tMesh->get_coefficient_indices_of_node( i, 0 ), "bspline inds" );
                // moris::print(tMesh->get_t_matrix_of_node_loc_ind(i,EntityRank::BSPLINE),"bspline t matrix");
            }

            //    moris::moris_id tElementInd = 8;
            //    Matrix<IndexMat> tFaces = tMesh->get_entity_connected_to_entity_loc_inds(tElementInd,EntityRank::ELEMENT,EntityRank::FACE);
            //    Matrix<IndexMat> tNeighborElements = tMesh->get_elements_connected_to_element_and_face_ind_loc_inds(tElementInd);
            //    moris::print(tNeighborElements,"tNeighborElements");
            //    moris::print(tFaces,"my faces");
            //    for(moris::uint i =0 ; i <tNeighborElements.n_cols(); i++)
            //    {
            //        moris::print(tMesh->get_entity_connected_to_entity_loc_inds(tNeighborElements(0,i),EntityRank::ELEMENT,EntityRank::FACE),"neighbors faces");
            //    }

            //    xtk::Geom_Field tFieldAsGeom(tField);

            moris::mtk::Mesh*        tMeshData = moris::mtk::create_interpolation_mesh( mtk::MeshType::STK, "xtk_test_sphere1.exo" );
            std::string              tLSFName  = "Circle";
            xtk::Mesh_Field_Geometry tLevelSetMesh( tMeshData, { tLSFName } );

            // Tell the geometry engine about the discrete field mesh and how to interpret phases
            Geometry_Engine tGeometryEngine( tLevelSetMesh );

            // Tell the XTK model that it should decompose with a C_HIERARCHY_TET4, on the same mesh that the level set field is defined on.
            size_t                            tModelDimension       = 3;
            Vector< enum Subdivision_Method > tDecompositionMethods = { Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4 };
            Model                             tXTKModel( tModelDimension, tMeshData, tGeometryEngine );
            tXTKModel.mVerbose = false;

            tXTKModel.set_HMR_mesh_ptr( tMesh );

            // Do the cutting
            tXTKModel.decompose( tDecompositionMethods );

            tXTKModel.perform_basis_enrichment();

            tXTKModel.perform_multilevel_enrichment_internal();

            //            moris::mtk::Mesh* tXTKMTK = tXTKModel.get_xtk_as_mtk();
            //
            //            // to implement
            //            xtk::Model* tXTKModelFromMTK = tXTKMTK.get_xtk_model();
            //
            //            Enrichment const & tEnrichment = tXTKModelFromMTK->get_basis_enrichment();

            Output_Options tOutputOptions;
            tOutputOptions.mAddNodeSets = false;
            tOutputOptions.mAddSideSets = false;

            moris::mtk::Mesh* tCutMeshData = tXTKModel.get_output_mesh( tOutputOptions );

            // std::string tPrefix = std::getenv("MORISOUTPUT");
            std::string tMeshOutputFile = "xtk_hmr_output.e";
            tCutMeshData->create_output_mesh( tMeshOutputFile );
            delete tCutMeshData;
            delete tMeshData;
        }
    }

}    // namespace moris::xtk
