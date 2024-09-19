/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_XTK_HMR_Multigrid.cpp
 *
 */

#include "catch.hpp"

#include "cl_XTK_Model.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"

#include "moris_typedefs.hpp"
#include "HDF5_Tools.hpp"
#include "paths.hpp"

#include "cl_MTK_Mesh_Manager.hpp"

#include "cl_MTK_Vertex.hpp"    //MTK
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh.hpp"

#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Integration_Mesh_STK.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Writer_Exodus.hpp"

#include "cl_Matrix.hpp"    //LINALG
#include "linalg_typedefs.hpp"
#include "fn_equal_to.hpp"    // ALG/src

#include "cl_HMR_Mesh_Interpolation.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Background_Mesh.hpp"      //HMR/src
#include "cl_HMR_BSpline_Mesh_Base.hpp"    //HMR/src
#include "cl_HMR_Element.hpp"              //HMR/src
#include "cl_HMR_Factory.hpp"              //HMR/src
#include "cl_HMR_Field.hpp"
#include "cl_HMR_Lagrange_Mesh_Base.hpp"    //HMR/src
#include "cl_HMR_Parameters.hpp"            //HMR/src

#include "cl_GEN_Line.hpp"

#include "fn_norm.hpp"

namespace moris::xtk
{
    moris::real
    CircleFuncXTKHMR2D( const moris::Matrix< DDRMat > &aPoint )
    {

        moris::real mXCenter = 0;
        moris::real mYCenter = 0;
        moris::real mRadius  = 1.1;

        return ( aPoint( 0 ) - mXCenter ) * ( aPoint( 0 ) - mXCenter )
             + ( aPoint( 1 ) - mYCenter ) * ( aPoint( 1 ) - mYCenter )
             - ( mRadius * mRadius );
    }

    moris::real
    PlaneFuncXTKHMR2D( const moris::Matrix< DDRMat > &aPoint )
    {
        return aPoint( 0 ) - 0.511;
    }

    TEST_CASE( "2D XTK WITH HMR MULLTIGRID 11", "[XTK_HMR_Multigrid]" )
    {

        if ( par_size() <= 1 )
        {
            std::string tFieldName = "Cylinder";

            moris::uint tLagrangeMeshIndex = 0;

            moris::hmr::Parameters tParameters;

            tParameters.set_number_of_elements_per_dimension( 10, 5 );
            tParameters.set_domain_dimensions( 2, 1 );
            tParameters.set_domain_offset( { { -1.0 }, { -0.5 } } );
            tParameters.set_bspline_truncation( true );

            tParameters.set_output_meshes( { { 0 } } );

            tParameters.set_lagrange_orders( { 1 } );
            tParameters.set_lagrange_patterns( { 0 } );

            tParameters.set_bspline_orders( { 1 } );
            tParameters.set_bspline_patterns( { 0 } );

            tParameters.set_side_sets( { { 1 }, { 2 }, { 3 }, { 4 } } );

            tParameters.set_union_pattern( 2 );
            tParameters.set_working_pattern( 3 );

            tParameters.set_refinement_buffer( 2 );
            tParameters.set_staircase_buffer( 2 );

            tParameters.set_multigrid( true );

            hmr::HMR tHMR( tParameters );

            std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

            // create field
            std::shared_ptr< moris::hmr::Field > tField = tMesh->create_field( tFieldName, tLagrangeMeshIndex );

            tField->evaluate_scalar_function( PlaneFuncXTKHMR2D );

            for ( uint k = 0; k < 2; ++k )
            {
                tHMR.flag_surface_elements_on_working_pattern( tField );
                tHMR.perform_refinement_based_on_working_pattern( 0 );

                tField->evaluate_scalar_function( PlaneFuncXTKHMR2D );
            }

            tHMR.finalize();

            tHMR.calculate_bspline_coordinates( tLagrangeMeshIndex, 0 );

            tHMR.save_to_exodus( 0, "./xtk_exo/xtk_hmr_2d_ip.e" );

            //         tHMR.save_bsplines_to_vtk( "./xtk_exo/Bspline.vtk", 0, 0 );

            hmr::Interpolation_Mesh_HMR *tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );

            Vector< std::shared_ptr< moris::gen::Level_Set_Geometry > > tGeometryVector( 1 );
            tGeometryVector( 0 ) = std::make_shared< moris::gen::Line >( 0.511, 0.0, 1.0, 0.0 );

            size_t                                 tModelDimension = 2;
            moris::gen::Geometry_Engine_Parameters tGeometryEngineParameters;
            tGeometryEngineParameters.mGeometries = tGeometryVector;
            moris::gen::Geometry_Engine tGeometryEngine( tInterpMesh, tGeometryEngineParameters );
            Model                       tXTKModel( tModelDimension, tInterpMesh, &tGeometryEngine );
            tXTKModel.mVerbose = false;

            // Specify decomposition Method and Cut Mesh ---------------------------------------
            Vector< enum Subdivision_Method > tDecompositionMethods = { Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3 };
            tXTKModel.decompose( tDecompositionMethods );

            tXTKModel.perform_basis_enrichment( mtk::EntityRank::BSPLINE, 0 );

            tXTKModel.construct_multigrid();

            // get meshes
            xtk::Enriched_Interpolation_Mesh &tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();

            uint tNumBasis = tEnrInterpMesh.get_num_basis( 0 );

            std::string tMorisRoot    = moris::get_base_moris_dir();
            std::string tHdf5FilePath = tMorisRoot + "/projects/XTK/test/xtk/data/Reference_Multigrid.hdf5";

            //------------------------------------------------------------------------------
            //    write solution ( uncomment this if you want to recreate solution files )
            //------------------------------------------------------------------------------

            //        // create file
            //        hid_t tFileID = create_hdf5_file( tHdf5FilePath );
            //
            //        // error handler
            //        herr_t tStatus = 0;
            //
            //        for( uint Ik = 0; Ik<tNumBasis; Ik++ )
            //        {
            //            std::string ChildBasisForCoarseBasis   = "ChildBasisForCoarseBasis_"   + std::to_string( Ik );
            //            std::string BasisWeightsForCoarseBasis = "BasisWeightsForCoarseBasis_" + std::to_string( Ik );
            //            std::string ParentBasisForFineBasis    = "ParentBasisForFineBasis_"    + std::to_string( Ik );
            //
            //            moris::Matrix< DDSMat > tMatInd     = tEnrInterpMesh.get_fine_basis_inds_of_basis( 0, Ik );
            //            moris::Matrix< DDRMat > tMatWeights = tEnrInterpMesh.get_fine_basis_weights_of_basis( 0, Ik );
            //            moris::Matrix< DDSMat > tMatFineToCoarseInd( tEnrInterpMesh.get_num_coarse_basis_of_basis( 0, Ik ),1);
            //
            //            for( uint Ii = 0; Ii<tEnrInterpMesh.get_num_coarse_basis_of_basis( 0, Ik ); Ii++ )
            //            {
            //            	tMatFineToCoarseInd( Ii )=tEnrInterpMesh.get_coarse_basis_index_of_basis( 0, Ik, Ii );
            //            }
            //
            //            // save data
            //            save_matrix_to_hdf5_file( tFileID, ChildBasisForCoarseBasis  , tMatInd            , tStatus );
            //            save_matrix_to_hdf5_file( tFileID, BasisWeightsForCoarseBasis, tMatWeights        , tStatus );
            //            save_matrix_to_hdf5_file( tFileID, ParentBasisForFineBasis   , tMatFineToCoarseInd, tStatus );
            //        }

            //------------------------------------------------------------------------------
            //    check solution
            //------------------------------------------------------------------------------

            // open file
            hid_t tFileID = open_hdf5_file( tHdf5FilePath );

            // error handler
            herr_t tStatus = 0;

            moris::Matrix< DDSMat > tMatIndRef;
            moris::Matrix< DDRMat > tMatWeightsRef;
            moris::Matrix< DDSMat > tMatFineToCoarseIndRef;

            for ( uint Ik = 0; Ik < tNumBasis; Ik++ )
            {
                std::string ChildBasisForCoarseBasis   = "ChildBasisForCoarseBasis_" + std::to_string( Ik );
                std::string BasisWeightsForCoarseBasis = "BasisWeightsForCoarseBasis_" + std::to_string( Ik );
                std::string ParentBasisForFineBasis    = "ParentBasisForFineBasis_" + std::to_string( Ik );

                moris::Matrix< DDSMat > tMatInd     = tEnrInterpMesh.get_fine_basis_inds_of_basis( 0, Ik );
                moris::Matrix< DDRMat > tMatWeights = tEnrInterpMesh.get_fine_basis_weights_of_basis( 0, Ik );
                moris::Matrix< DDSMat > tMatFineToCoarseInd( tEnrInterpMesh.get_num_coarse_basis_of_basis( 0, Ik ), 1 );

                for ( uint Ii = 0; Ii < tEnrInterpMesh.get_num_coarse_basis_of_basis( 0, Ik ); Ii++ )
                {
                    tMatFineToCoarseInd( Ii ) = tEnrInterpMesh.get_coarse_basis_index_of_basis( 0, Ik, Ii );
                }

                // read solution from file
                load_matrix_from_hdf5_file( tFileID, ChildBasisForCoarseBasis, tMatIndRef, tStatus );
                load_matrix_from_hdf5_file( tFileID, BasisWeightsForCoarseBasis, tMatWeightsRef, tStatus );
                load_matrix_from_hdf5_file( tFileID, ParentBasisForFineBasis, tMatFineToCoarseIndRef, tStatus );

                bool tCheck = true;
                for ( uint Ik = 0; Ik < tMatInd.numel(); Ik++ )
                {
                    if ( tMatInd( Ik ) != tMatIndRef( Ik ) ) { tCheck = false; }
                }
                for ( uint Ik = 0; Ik < tMatWeights.numel(); Ik++ )
                {
                    if ( tMatWeights( Ik ) != tMatWeightsRef( Ik ) ) { tCheck = false; }
                }
                for ( uint Ik = 0; Ik < tMatFineToCoarseInd.numel(); Ik++ )
                {
                    if ( tMatFineToCoarseInd( Ik ) != tMatFineToCoarseIndRef( Ik ) ) { tCheck = false; }
                }
                CHECK( tCheck );
            }

            // close file
            close_hdf5_file( tFileID );
            delete tInterpMesh;

            /*        // output to exodus file ----------------------------------------------------------
                    xtk::Enrichment const & tEnrichment = tXTKModel.get_basis_enrichment();

                    moris_index tSSIndex = tEnrIgMesh.create_side_set_from_dbl_side_set(1,"ghost_ss_p0");
                    tEnrIgMesh.create_block_set_from_cells_of_side_set(tSSIndex,"ghost_bs_p0", mtk::CellTopology::QUAD4);

                     // Declare the fields related to enrichment strategy in output options
                     Vector<std::string> tEnrichmentFieldNames = tEnrichment.get_cell_enrichment_field_names();

                    // output solution and meshes
                    xtk::Output_Options tOutputOptions;
                    tOutputOptions.mAddNodeSets = false;
                    tOutputOptions.mAddSideSets = true;
                    tOutputOptions.mAddClusters = false;

                    // add solution field to integration mesh
                    std::string tIntegSolFieldName = "solution";
                    tOutputOptions.mRealNodeExternalFieldNames = {tIntegSolFieldName};
                    tOutputOptions.mRealElementExternalFieldNames = tEnrichmentFieldNames;

                    moris::mtk::Integration_Mesh* tIntegMesh1 = tXTKModel.get_output_mesh(tOutputOptions);

                    tEnrichment.write_cell_enrichment_to_fields(tEnrichmentFieldNames,tIntegMesh1);

                    std::string tMeshOutputFile ="./xtk_exo/xtk_hmr_2d_ig_stk.e";
                    tIntegMesh1->create_output_mesh(tMeshOutputFile);

                    // Write mesh
                    moris::mtk::Writer_Exodus writer(&tEnrIgMesh);
                    writer.write_mesh("", "xtk_hmr_2d_ig_multigrid.exo", "", "xtk_temp.exo");

                    // Write the fields
                    writer.set_time(0.0);
                    writer.close_file();

                    delete tIntegMesh1;
            */
        }
    }

}    // namespace moris::xtk
