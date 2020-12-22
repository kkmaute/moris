// /*
//  * UT_XTK_HMR_2D.cpp
//  *
//  *  Created on: Sep 10, 2019
//  *      Author: doble
//  */

// #include "catch.hpp"

// #include "cl_XTK_Model.hpp"
// #include "cl_XTK_Enriched_Integration_Mesh.hpp"

// #include "typedefs.hpp"

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

// namespace xtk
// {
// moris::real
// CircleFunc_PHASE0(const moris::Matrix< moris::DDRMat > & aPoint )
// {

//     moris::real mXCenter = -1;
//     moris::real mYCenter = 0;
//     moris::real mRadius = 1.1111;

//     return    (aPoint(0) - mXCenter) * (aPoint(0) - mXCenter)
//             + (aPoint(1) - mYCenter) * (aPoint(1) - mYCenter)
//             - (mRadius * mRadius);
// }

// moris::real
// CircleFunc_PHASE1(const moris::Matrix< moris::DDRMat > & aPoint )
// {

//     moris::real mXCenter = 0.0;
//     moris::real mYCenter = 0;
//     moris::real mRadius = 1.1111;

//     return    (aPoint(0) - mXCenter) * (aPoint(0) - mXCenter)
//             + (aPoint(1) - mYCenter) * (aPoint(1) - mYCenter)
//             - (mRadius * mRadius);
// }

// moris::real
// CircleFunc_PHASE2(const moris::Matrix< moris::DDRMat > & aPoint )
// {

//     moris::real mXCenter = 1.0;
//     moris::real mYCenter = 0;
//     moris::real mRadius = 1.1111;

//     return    (aPoint(0) - mXCenter) * (aPoint(0) - mXCenter)
//             + (aPoint(1) - mYCenter) * (aPoint(1) - mYCenter)
//             - (mRadius * mRadius);
// }

// TEST_CASE("Phase Table Scratch","[XTK_PHASE_TABLE]")
// {
//     moris::uint iOrder =  1;

//             std::string tFieldName0 = "Sphere0";
//             std::string tFieldName1 = "Sphere1";
//             std::string tFieldName2 = "Sphere2";

//             moris::uint tLagrangeMeshIndex = 0;

//             moris::hmr::Parameters tParameters;

//             tParameters.set_number_of_elements_per_dimension( { {60}, {20}} );
//             tParameters.set_domain_dimensions({ {6}, {4} });
//             tParameters.set_domain_offset({ {-3.0}, {-2.0} });
//             tParameters.set_bspline_truncation( true );

//             tParameters.set_output_meshes( { {0} } );

//             tParameters.set_lagrange_orders  ( { {iOrder} });
//             tParameters.set_lagrange_patterns({ {0} });

//             tParameters.set_bspline_orders   ( { {iOrder} } );
//             tParameters.set_bspline_patterns ( { {0} } );

//             tParameters.set_side_sets({{1},{2},{3},{4} });

//             tParameters.set_union_pattern( 2 );
//             tParameters.set_working_pattern( 3 );

//             tParameters.set_refinement_buffer( 2 );
//             tParameters.set_staircase_buffer( 2);
//             tParameters.set_number_aura(true);

//             Cell< Matrix< DDSMat > > tLagrangeToBSplineMesh( 1 );
//             tLagrangeToBSplineMesh( 0 ) = { {0} };

//             tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

//             hmr::HMR tHMR( tParameters );

//             std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

//             // create field
//             std::shared_ptr< moris::hmr::Field > tField0 = tMesh->create_field( tFieldName0, tLagrangeMeshIndex );
//             std::shared_ptr< moris::hmr::Field > tField1 = tMesh->create_field( tFieldName1, tLagrangeMeshIndex );
//             std::shared_ptr< moris::hmr::Field > tField2 = tMesh->create_field( tFieldName2, tLagrangeMeshIndex );

//             tField0->evaluate_scalar_function( CircleFunc_PHASE0 );
//             tField1->evaluate_scalar_function( CircleFunc_PHASE1 );
//             tField2->evaluate_scalar_function( CircleFunc_PHASE2 );

//             for( uint k=0; k<3; ++k )
//             {
                
//             tField0->evaluate_scalar_function( CircleFunc_PHASE0 );
//             tField1->evaluate_scalar_function( CircleFunc_PHASE1 );
//             tField2->evaluate_scalar_function( CircleFunc_PHASE2 );

//                 tHMR.flag_surface_elements_on_working_pattern( tField0 );
//                 tHMR.flag_surface_elements_on_working_pattern( tField1 );
//                 tHMR.flag_surface_elements_on_working_pattern( tField2 );

//                 tHMR.perform_refinement_based_on_working_pattern( 0 );
//             }

//             tField0->evaluate_scalar_function( CircleFunc_PHASE0 );
//             tField1->evaluate_scalar_function( CircleFunc_PHASE1 );
//             tField2->evaluate_scalar_function( CircleFunc_PHASE2 );

//             tHMR.finalize();

//             hmr::Interpolation_Mesh_HMR * tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex  );


//             moris::Cell< std::shared_ptr<moris::ge::Geometry> > tGeometryVector(3);
//             tGeometryVector(0) = std::make_shared<moris::ge::Circle>(-1.0, 0.0, 1.111);
//             tGeometryVector(1) = std::make_shared<moris::ge::Circle>( 0.0, 0.0, 1.111);
//             tGeometryVector(2) = std::make_shared<moris::ge::Circle>( 1.0, 0.0, 1.111);

//             size_t tModelDimension = 2;

//             // setup the phase table
//             moris_index tNumGeom      = tGeometryVector.size();
//             // moris::ge::Phase_Table tPhaseTable( tNumGeom ); 
//             //moris_index tNumBulkPhase = 3;
//             Matrix<DDUMat> tGeomIndexToBulkPhase = {{0,0,2,2,0,0,2,1}};
//             moris::ge::Phase_Table tPhaseTable( tNumGeom, tGeomIndexToBulkPhase );

//             tPhaseTable.print();

//             moris::ge::Geometry_Engine tGeometryEngine(tGeometryVector, tPhaseTable, tInterpMesh);
//             Model tXTKModel(tModelDimension, tInterpMesh, &tGeometryEngine);
//             tXTKModel.mVerbose  =  true;

//             //Specify decomposition Method and Cut Mesh ---------------------------------------
//             Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3};
//             tXTKModel.decompose(tDecompositionMethods);

//             tXTKModel.perform_basis_enrichment(EntityRank::BSPLINE,0);


//             tXTKModel.construct_face_oriented_ghost_penalization_cells();


//             // output to exodus file ----------------------------------------------------------
//            // Write mesh
//            moris::mtk::Writer_Exodus writer(&tXTKModel.get_enriched_integ_mesh(0));
//            writer.write_mesh("", "./xtk_exo/phase_table_2n.exo", "", "temp.exo");
//            // Write the fields
//            writer.set_time(0.0);
//            writer.close_file();


//             delete tInterpMesh;

// }

// }
