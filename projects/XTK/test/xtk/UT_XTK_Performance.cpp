///*
// * UT_XTK_HMR.cpp
// *
// *  Created on: Aug 29, 2019
// *      Author: doble
// */
//
//#include "catch.hpp"
//
//#include "cl_XTK_Model.hpp"
//#include "cl_XTK_Enriched_Integration_Mesh.hpp"
//#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
//#include "cl_XTK_Vertex_Enrichment.hpp"
//#include "cl_XTK_Hole_Seeder.hpp"
//#include "cl_GEN_Sphere_Box.hpp"
////#include "cl_Geom_Field.hpp"
//#include "typedefs.hpp"
//
//#include "cl_MTK_Mesh_Manager.hpp"
//
//#include "cl_MTK_Vertex.hpp"    //MTK
//#include "cl_MTK_Cell.hpp"
//#include "cl_MTK_Enums.hpp"
//#include "cl_MTK_Mesh.hpp"
//#include "cl_MTK_Mesh_Checker.hpp"
//
//#include "cl_MTK_Mesh_Manager.hpp"
//#include "cl_MTK_Integration_Mesh_STK.hpp"
//#include "cl_MTK_Interpolation_Mesh.hpp"
//#include "cl_MTK_Integration_Mesh.hpp"
//#include "cl_MTK_Writer_Exodus.hpp"
//
//#include "cl_Matrix.hpp"        //LINALG
//#include "linalg_typedefs.hpp"
//#include "fn_equal_to.hpp" // ALG/src
//
//#include "cl_HMR_Mesh_Interpolation.hpp"
//#include "cl_HMR.hpp"
//#include "cl_HMR_Background_Mesh.hpp" //HMR/src
//#include "cl_HMR_BSpline_Mesh_Base.hpp" //HMR/src
//#include "cl_HMR_Element.hpp" //HMR/src
//#include "cl_HMR_Factory.hpp" //HMR/src
//#include "cl_HMR_Field.hpp"
//#include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src
//#include "cl_HMR_Parameters.hpp" //HMR/src
//
//#include "cl_GEN_User_Defined_Geometry.hpp"
//
//#include "fn_norm.hpp"
//#include "fn_compute_interface_surface_area.hpp"
//
//#include <map>
//namespace moris
//{
//
//    // User defined gyroid
//    moris::real
//    LevelSetSphereGyroid(const moris::Matrix< moris::DDRMat > & aPoint )
//    {
//        moris::real tFuncVal = std::sin(aPoint(0))*std::cos(aPoint(1))+
//                        std::sin(aPoint(1))*std::cos(aPoint(2))+
//                        std::sin(aPoint(2))*std::cos(aPoint(0));
//        return tFuncVal;
//    }
//
//    moris::real LevelSetSphereGyroidGeometry(const moris::Matrix<moris::DDRMat>& aCoordinates, const moris::Cell<moris::real*>& aParameters)
//    {
//        return LevelSetSphereGyroid(aCoordinates);
//    }
//
//TEST_CASE("XTK Performance Test","[XTK_PERF]")
//{
//        // sweep parameters
//        moris::real tInitRad    = 3.15;
//        moris::real tInitNexp   = 2;
//        moris::uint tInitNumX   = 2;
//        moris::uint tInitNumY   = 2;
//        moris::uint tInitNumZ   = 2;
//        moris::uint tInitNumRef = 0;
//
//        // num steps
//        moris::uint tNumRadSteps  = 1;
//        moris::uint tNumNexpSteps = 1;
//        moris::uint tNumNumXSteps = 1;
//        moris::uint tNumNumYSteps = 1;
//        moris::uint tNumNumZSteps = 1;
//        moris::uint tNumRefSteps  = 1;
//
//        // increment
//        moris::real tRadInc  = 0.1;
//        moris::real tNexpInc = 0.0;
//        moris::uint tNumXInc = 1;
//        moris::uint tNumYInc = 0;
//        moris::uint tNumZInc = 0;
//        moris::uint tRefinementInc = 1;
//
//        // Keep track of how many runs for hdf file naming
//        moris::uint tCount = 0;
//
//
//        moris::uint iOrder = 2;
//        for(moris::uint iR = 0; iR < tNumRadSteps; iR++)
//        {
//            // geometric info
//            moris::real tCurrentRad    = tInitRad;
//            moris::real tCurrentNexp   = tInitNexp;
//            moris::uint tCurrentNumX   = tInitNumX + iR*tNumXInc;
//            moris::uint tCurrentNumY   = tInitNumY + iR*tNumYInc;
//            moris::uint tCurrentNumZ   = tInitNumZ + iR*tNumZInc;
//
//            // number of geometries
//            moris::uint tCurrentNumGeom = tCurrentNumX*tCurrentNumY*tCurrentNumZ;
//
//            std::string tFieldName = "Geom";
//
//            moris::uint tLagrangeMeshIndex = 0;
//            moris::uint tBSplineMeshIndex = 0;
//
//            moris::hmr::Parameters tParameters;
//
//            tParameters.set_number_of_elements_per_dimension( { {50}, {50}, {50} } );
//            tParameters.set_domain_dimensions({ {10}, {10}, {10} });
//            tParameters.set_domain_offset({ {0.0}, {0.0}, {0.0} });
//            tParameters.set_side_sets({ {5}, {6} });
//            tParameters.set_bspline_truncation( true );
//            tParameters.set_lagrange_orders  ( { {iOrder} });
//            tParameters.set_lagrange_patterns( { {0} });
//            tParameters.set_bspline_orders   ( { {iOrder} } );
//            tParameters.set_bspline_patterns ( { {0} } );
//            tParameters.set_output_meshes( { {0} } );
//            tParameters.set_staircase_buffer( 1 );
//            tParameters.set_initial_refinement( 0 );
//            tParameters.set_number_aura( true );
//            Cell< Matrix< DDSMat > > tLagrangeToBSplineMesh( 1 );
//            tLagrangeToBSplineMesh( 0 ) = { {0} };
//            tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );
//
//            // create the HMR object by passing the settings to the constructor
//            moris::hmr::HMR tHMR( tParameters );
//
//            std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );
//
//            // seed the field
//            xtk::Hole_Seeder tHoleSeeder(
//                    tMesh.get(),
//                    tCurrentRad,
//                    tCurrentRad,
//                    tCurrentRad,
//                    tCurrentNexp,
//                    tCurrentNumX,
//                    tCurrentNumY,
//                    tCurrentNumZ);
//
//            // perform hole seeding
//            tHoleSeeder.seed_field();
//
//            // get the spheres
//            moris::Cell<std::shared_ptr<moris::ge::Sphere_Box>> & tGeometries = tHoleSeeder.get_seeded_geometies();
//
//            // create fie
//            moris::Cell<std::shared_ptr< moris::hmr::Field >> tFields ( tGeometries.size()+1 );
//
//            for(moris::uint iGeom = 0; iGeom  < tGeometries.size(); iGeom++)
//            {
//                tFields(iGeom) = tMesh->create_field( tFieldName + "_"+ std::to_string(iGeom), tLagrangeMeshIndex );
//                tFields(iGeom)->evaluate_scalar_function( *tGeometries(iGeom) );
//            }
//
//            std::shared_ptr< moris::hmr::Field > tGyroidField = tMesh->create_field( "Gyroid", tLagrangeMeshIndex );
//            tGyroidField->evaluate_scalar_function( LevelSetSphereGyroid );
//
//            for( uint k=0; k<tInitNumRef; ++k )
//            {
//                for(moris::uint iGeom = 0; iGeom  < tGeometries.size(); iGeom++)
//                {
//                    tHMR.flag_surface_elements_on_working_pattern( tFields(iGeom) );
//                }
//
//                tHMR.flag_surface_elements_on_working_pattern( tGyroidField );
//
//
//                tHMR.perform_refinement_based_on_working_pattern( 0 );
//
//                for(moris::uint iGeom = 0; iGeom  < tGeometries.size(); iGeom++)
//                {
//                    tFields(iGeom)->evaluate_scalar_function( *tGeometries(iGeom) );
//                }
//                tGyroidField->evaluate_scalar_function( LevelSetSphereGyroid );
//            }
//
//            tHMR.finalize();
//
//            hmr::Interpolation_Mesh_HMR * tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex  );
//
//            // convert geoms to cell of base class
//            moris::Cell<std::shared_ptr<moris::ge::Geometry>> tGeometryVector(tGeometries.size());
//            for(moris::uint iGeom = 0; iGeom  < tGeometries.size(); iGeom++)
//            {
//                tGeometryVector(iGeom) = tGeometries(iGeom);
//            }
//
//            tGeometryVector.push_back( std::make_shared<moris::ge::User_Defined_Geometry>(Matrix<DDRMat>(0, 0), &(LevelSetSphereGyroidGeometry)));
//
//
//            size_t tModelDimension = 3;
//            moris::ge::Phase_Table tPhaseTable (tGeometryVector.size());
//            moris::ge::Geometry_Engine tGeometryEngine(tGeometryVector,tPhaseTable,tModelDimension);
//            xtk::Model tXTKModel(tModelDimension,tInterpMesh,&tGeometryEngine);
//            tXTKModel.mVerbose  =  true;
//
//            // Tell the XTK model that it should decompose with a C_HIERARCHY_TET4, on the same mesh that the level set field is defined on.
//            Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8,Subdivision_Method::C_HIERARCHY_TET4};
//
//            // Do the cutting
//            tXTKModel.decompose(tDecompositionMethods);
//
//            // Perform the enrichment
//            tXTKModel.perform_basis_enrichment(EntityRank::BSPLINE,0);
//
//            // perform ghost stabilization
//            tXTKModel.construct_face_oriented_ghost_penalization_cells();
//
//            // print the timing data to console
//            tXTKModel.print_timing_data();
//
//            // save the timing data to a file
////            std::string tFileName = "./xtk_timing_"+std::to_string(tCount) + "_"+ std::to_string(par_size()) + ".hdf";
//            std::string tFileName = "./xtk_timing_"+ std::to_string(par_size()) + ".hdf";
//            tXTKModel.save_timing_to_hdf5(tFileName);
//
//            // save the geometric data to a file
////            std::string tStatsFileName = "./xtk_model_stats_"+std::to_string(tCount) + "_"+ std::to_string(par_size()) + ".hdf";
//            std::string tStatsFileName = "./xtk_model_stats_"+ std::to_string(par_size()) + ".hdf";
//            tXTKModel.save_model_statistics_to_file(tStatsFileName);
//
//
//            xtk::Enriched_Integration_Mesh   & tEnrIntegMesh = tXTKModel.get_enriched_integ_mesh();
//            std::string tEnrIgMeshFileName = "./xtk_exo/performance_test.e";
//            tEnrIntegMesh.deactivate_empty_sets();
//            // Write mesh
//            moris::mtk::Writer_Exodus writer(&tEnrIntegMesh);
//            writer.write_mesh("", tEnrIgMeshFileName, "", "xtk_temp.exo");
//
//            // Write the fields
//            writer.set_time(0.0);
//            writer.close_file();
//
//
//
//            // increment sweep count
//            tCount++;
//            delete tInterpMesh;
//        }
//
//
//}
//}
//
//
