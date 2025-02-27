/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_XTK_HMR_2D.cpp
 *
 */

// #include "catch.hpp"

// #include "cl_XTK_Model.hpp"
// #include "cl_XTK_Enriched_Integration_Mesh.hpp"

// #include "moris_typedefs.hpp"

// #include "cl_MTK_Mesh_Manager.hpp"

// #include "cl_MTK_Vertex.hpp"    //MTK
// #include "cl_MTK_Cell.hpp"
// #include "cl_MTK_Enums.hpp"
// #include "cl_MTK_Mesh.hpp"

// #include "cl_MTK_Mesh_Manager.hpp"
// #include "cl_MTK_Integration_Mesh_STK.hpp"
// #include "cl_MTK_Interpolation_Mesh.hpp"
// #include "cl_MTK_Integration_Mesh.hpp"
// #include "cl_MTK_Writer_Exodus.hpp"

// #include "cl_Matrix.hpp"        //LINALG
// #include "linalg_typedefs.hpp"
// #include "fn_equal_to.hpp" // ALG/src

// #include "cl_HMR_Mesh_Interpolation.hpp"
// #include "cl_HMR.hpp"
// #include "cl_HMR_Background_Mesh.hpp" //HMR/src
// #include "cl_HMR_BSpline_Mesh_Base.hpp" //HMR/src
// #include "cl_HMR_Element.hpp" //HMR/src
// #include "cl_HMR_Factory.hpp" //HMR/src
// #include "cl_HMR_Field.hpp"
// #include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src
// #include "cl_HMR_Parameters.hpp" //HMR/src

// #include "cl_GEN_Plane.hpp"
// #include "cl_GEN_Circle.hpp"

// #include "fn_norm.hpp"

// namespace moris::xtk
// {
// moris::real
// CircleFuncXTKHMR2D(const moris::Matrix< DDRMat > & aPoint )
// {

//     moris::real mXCenter = 0;
//     moris::real mYCenter = 0;
//     moris::real mRadius = 1.1;

//     return    (aPoint(0) - mXCenter) * (aPoint(0) - mXCenter)
//             + (aPoint(1) - mYCenter) * (aPoint(1) - mYCenter)
//             - (mRadius * mRadius);
// }

// moris::real
// PlaneFuncXTKHMR2D(const moris::Matrix< DDRMat > & aPoint )
// {
//     // Get variables
//     real tXCenter = 0.0;
//     real tYCenter = 0.0;
//     real tXNormal = 1.0;
//     real tYNormal = 1.0;

//     // Evaluate field value
//     return tXNormal * (aPoint(0) - tXCenter) + tYNormal * (aPoint(1) - tYCenter);

// }

// TEST_CASE("2D XTK WITH HMR","[XTK_HMR_2D]")
// {
//     if(par_size()<=2)
//     {
//         for( moris::uint iOrder = 1; iOrder < 4; iOrder ++)
//         {

//             std::string tFieldName = "Cylinder";

//             moris::uint tLagrangeMeshIndex = 0;

//             moris::hmr::Parameters tParameters;

//             tParameters.set_number_of_elements_per_dimension( { {24}, {24}} );
//             tParameters.set_domain_dimensions({ {2}, {2} });
//             tParameters.set_domain_offset({ {-1.0}, {-1.0} });
//             tParameters.set_bspline_truncation( true );

//             tParameters.set_output_meshes({{ {0} }} );

//             tParameters.set_lagrange_orders  ( { {iOrder} });
//             tParameters.set_lagrange_patterns({ {0} });

//             tParameters.set_bspline_orders   ( { {iOrder} } );
//             tParameters.set_bspline_patterns ( { {0} } );

//             tParameters.set_create_side_sets( true );

//             tParameters.set_union_pattern( 2 );
//             tParameters.set_working_pattern( 3 );

//             tParameters.set_refinement_buffer( 2 );
//             tParameters.set_staircase_buffer( 2);
//             tParameters.set_number_aura(true);

//             Vector< Matrix< DDSMat > > tLagrangeToBSplineMesh( 1 );
//             tLagrangeToBSplineMesh( 0 ) = { {0} };

//             tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

//             hmr::HMR tHMR( tParameters );

//             std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

//             // create field
//             std::shared_ptr< moris::hmr::Field > tField = tMesh->create_field( tFieldName, tLagrangeMeshIndex );

//             tField->evaluate_scalar_function( CircleFuncXTKHMR2D );

//             for( uint k=0; k<3; ++k )
//             {
//                 tHMR.flag_surface_elements_on_working_pattern( tField );
//                 tHMR.perform_refinement_based_on_working_pattern( 0 );

//                 tField->evaluate_scalar_function( CircleFuncXTKHMR2D );
//             }

//             tHMR.finalize();

//             hmr::Interpolation_Mesh_HMR * tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex  );

//             Vector< std::shared_ptr<moris::gen::Level_Set_Geometry> > tGeometryVector(1);
//             tGeometryVector(0) = std::make_shared<moris::gen::Circle>(0.0, 0.0, 1.1);

//             size_t tModelDimension = 2;
//             moris::gen::Geometry_Engine_Parameters tGeometryEngineParameters;
//             tGeometryEngineParameters.mGeometries = tGeometryVector;
//             moris::gen::Geometry_Engine tGeometryEngine(tInterpMesh, tGeometryEngineParameters);
//             Model tXTKModel(tModelDimension, tInterpMesh, &tGeometryEngine);
//             tXTKModel.mVerbose  =  false;

//             //Specify decomposition Method and Cut Mesh ---------------------------------------
//             Vector<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3};
//             tXTKModel.decompose(tDecompositionMethods);

//             tXTKModel.perform_basis_enrichment( mtk::EntityRank::BSPLINE,0);

//             tXTKModel.construct_face_oriented_ghost_penalization_cells();

//             delete tInterpMesh;
//         }
//     }
// }

// TEST_CASE("2D XTK WITH HMR WEIRD INTERSECTION","[XTK_HMR_2D_WI]")
// {
//     if(par_size()<=1)
//     {
//         std::string tFieldName = "Cylinder";

//         moris::uint tLagrangeMeshIndex = 0;

//         moris::hmr::Parameters tParameters;

//         tParameters.set_number_of_elements_per_dimension( { {2}, {2}} );
//         tParameters.set_domain_dimensions({ {2}, {2} });
//         tParameters.set_domain_offset({ {-1.0}, {-1.0} });
//         tParameters.set_bspline_truncation( true );

//         tParameters.set_output_meshes( {{ {0} }} );

//         tParameters.set_lagrange_orders  ( { {1} });
//         tParameters.set_lagrange_patterns({ {0} });

//         tParameters.set_bspline_orders   ( { {1} } );
//         tParameters.set_bspline_patterns ( { {0} } );

//         tParameters.set_create_side_sets( true );

//         tParameters.set_union_pattern( 2 );
//         tParameters.set_working_pattern( 3 );

//         tParameters.set_refinement_buffer( 2 );
//         tParameters.set_staircase_buffer( 2 );

//         Vector< Matrix< DDSMat > > tLagrangeToBSplineMesh( 1 );
//         tLagrangeToBSplineMesh( 0 ) = { {0} };

//         tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

//         hmr::HMR tHMR( tParameters );

//         std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

//         // create field
//         std::shared_ptr< moris::hmr::Field > tField = tMesh->create_field( tFieldName, tLagrangeMeshIndex );

//         tField->evaluate_scalar_function( CircleFuncXTKHMR2D );

//         for( uint k=0; k<2; ++k )
//         {
//             tHMR.flag_surface_elements_on_working_pattern( tField );
//             tHMR.perform_refinement_based_on_working_pattern( 0 );

//             tField->evaluate_scalar_function( CircleFuncXTKHMR2D );
//         }

//         tHMR.finalize();

//         tHMR.save_to_exodus( 0, "./xtk_exo/xtk_hmr_wi_2d_ip.e" );

//         hmr::Interpolation_Mesh_HMR * tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex  );

//         // create a plane which intentionally intersects from fine to coarse
//         moris::Matrix<DDRMat> tCenters = {{ 0.1,0.1 }};
//         moris::Matrix<DDRMat> tNormals = {{ 1.0,0.0 }};
//         Vector<std::shared_ptr<moris::gen::Level_Set_Geometry>> tGeometry(1);
//         tGeometry(0) = std::make_shared<moris::gen::Plane>(tCenters(0), tCenters(1), tNormals(0), tNormals(1));

//         size_t tModelDimension = 2;
//         moris::gen::Geometry_Engine_Parameters tGeometryEngineParameters;
//         tGeometryEngineParameters.mGeometries = tGeometry;
//         moris::gen::Geometry_Engine tGeometryEngine(tInterpMesh, tGeometryEngineParameters);
//         Model tXTKModel(tModelDimension, tInterpMesh, &tGeometryEngine);
//         tXTKModel.mVerbose  =  false;

//         //Specify decomposition Method and Cut Mesh ---------------------------------------
//         Vector<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3};
//         bool tSuccess = tXTKModel.decompose(tDecompositionMethods);

//         CHECK(!tSuccess);

//     }
// }

// TEST_CASE("2D Conformal Coincident Subdivision","[CM_2D_LIN_COIN]")
// {

//     if(par_size()<=1)
//     {
//         std::string tFieldName = "Cylinder";

//         moris::uint tLagrangeMeshIndex = 0;

//         moris::hmr::Parameters tParameters;

//         tParameters.set_number_of_elements_per_dimension( { {1}, {1}} );
//         tParameters.set_domain_dimensions({ {2}, {2} });
//         tParameters.set_domain_offset({ {-1.0}, {-1.0} });
//         tParameters.set_bspline_truncation( true );

//         tParameters.set_output_meshes( {{ {0} }} );

//         tParameters.set_lagrange_orders  ( { {1} });
//         tParameters.set_lagrange_patterns({ {0} });

//         tParameters.set_bspline_orders   ( { {1} } );
//         tParameters.set_bspline_patterns ( { {0} } );

//         tParameters.set_create_side_sets( true );

//         tParameters.set_union_pattern( 2 );
//         tParameters.set_working_pattern( 3 );

//         tParameters.set_refinement_buffer( 2 );
//         tParameters.set_staircase_buffer( 2 );

//         Vector< Matrix< DDSMat > > tLagrangeToBSplineMesh( 1 );
//         tLagrangeToBSplineMesh( 0 ) = { {0} };

//         tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

//         hmr::HMR tHMR( tParameters );

//         std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

//         // create field
//         std::shared_ptr< moris::hmr::Field > tField = tMesh->create_field( tFieldName, tLagrangeMeshIndex );

//         tField->evaluate_scalar_function( PlaneFuncXTKHMR2D );

//         tHMR.finalize();

//         tHMR.save_to_exodus( 0, "./xtk_exo/xtk_hmr_2d_coincidence.e" );

//         hmr::Interpolation_Mesh_HMR * tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex  );

//         // create a plane which intentionally intersects from fine to coarse
//         moris::Matrix<DDRMat> tCenters = {{ 0.0,0.0 }};
//         moris::Matrix<DDRMat> tNormals = {{ 1.0,0.0 }};
//         Vector<std::shared_ptr<moris::gen::Level_Set_Geometry>> tGeometry(2);
//         tGeometry(0) = std::make_shared<moris::gen::Plane>(tCenters(0), tCenters(1), tNormals(0), tNormals(1)); // center vertical
//         tGeometry(1) = std::make_shared<moris::gen::Plane>(tCenters(0), tCenters(1), tNormals(1), tNormals(0)); // center horizontal

//         size_t tModelDimension = 2;
//         moris::gen::Geometry_Engine_Parameters tGeometryEngineParameters;
//         tGeometryEngineParameters.mGeometries = tGeometry;
//         tGeometryEngineParameters.mIsocontourTolerance = 1E-8;
//         moris::gen::Geometry_Engine tGeometryEngine(tInterpMesh, tGeometryEngineParameters);
//         Model tXTKModel(tModelDimension, tInterpMesh, &tGeometryEngine);
//         tXTKModel.mVerbose  =  true;

//         //Specify decomposition Method and Cut Mesh ---------------------------------------
//         Vector<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3};
//         tXTKModel.decompose(tDecompositionMethods);

//         tXTKModel.perform_basis_enrichment( mtk::EntityRank::BSPLINE,0);

//         // Write mesh
//         xtk::Enriched_Integration_Mesh & tEnrIgMesh = tXTKModel.get_enriched_integ_mesh(0);

//         tEnrIgMesh.deactivate_empty_sets();

//         moris::mtk::Writer_Exodus writer(&tEnrIgMesh);
//         writer.write_mesh("", "./xtk_2d_coincident.exo", "", "temp.exo");

//         delete tInterpMesh;
//     }
// }

// }
