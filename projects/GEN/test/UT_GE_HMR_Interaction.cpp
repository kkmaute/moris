#include <catch.hpp>

#include "cl_HMR.hpp"
#include "cl_HMR_Background_Mesh.hpp" //HMR/src
#include "cl_HMR_Background_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_BSpline_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Element.hpp" //HMR/src
#include "cl_HMR_Factory.hpp" //HMR/src
#include "cl_HMR_Field.hpp"
#include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Parameters.hpp" //HMR/src

#include "cl_Communication_Manager.hpp" //COM/src
#include "cl_Communication_Tools.hpp" //COM/src
#include "typedefs.hpp" //COR/src
#include "cl_Matrix.hpp" //LINALG/src

#include "HDF5_Tools.hpp"

// GE includes
#include "cl_GE_Core.hpp"
#include "cl_GE_Factory.hpp"
#include "cl_GE_Intersection_Object_Line.hpp"
#include "cl_GE_Node.hpp"

// LINALG includes
#include "cl_Matrix.hpp"
#include "fn_all_true.hpp"
#include "fn_equal_to.hpp"
#include "linalg_typedefs.hpp"
#include "op_equal_equal.hpp"

// MTK includes
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_Mesh_Factory.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Mesh_Tools.hpp"
#include "cl_MTK_Scalar_Field_Info.hpp"
#include "cl_MTK_Vertex.hpp"

#include "cl_HMR_Factory.hpp" //HMR/src
#include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src

using namespace moris;
using namespace hmr;
using namespace ge;

TEST_CASE("GE_HMR_Interaction","[moris],[GE],[GE_HMR_Interaction]")
{
    if(par_size() == 1)
    {
        for( moris::uint tOrder=1; tOrder<=1; tOrder++ )
        {
            uint tLagrangeMeshIndex = 0;
            // empty container for B-Spline meshes
            moris::Cell< moris::hmr::BSpline_Mesh_Base* > tBSplineMeshes;

            // create settings object
            moris::hmr::Parameters tParameters;

            tParameters.set_number_of_elements_per_dimension( { {4}, {4} } );

            // B-Spline truncation is turned on by default.
            // It is recommended to leave this setting as is.
            tParameters.set_bspline_truncation( true );

            tParameters.set_lagrange_orders  ( { {tOrder} });
            tParameters.set_lagrange_patterns({ {2} });

            tParameters.set_bspline_orders   ( { {tOrder}, {tOrder}} );
            tParameters.set_bspline_patterns ( { {0}, {1}} );

            tParameters.set_staircase_buffer( tOrder );

            tParameters.set_initial_refinement( 0 );

            Cell< Matrix< DDUMat > > tLagrangeToBSplineMesh( 1 );
            tLagrangeToBSplineMesh( 0 ) = { {0}, {1} };

            tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

            //------------------------------------------------------------------------------
            //  HMR Initialization
            //------------------------------------------------------------------------------

            // create the HMR object by passing the settings to the constructor
            moris::hmr::HMR tHMR( tParameters );

            // std::shared_ptr< Database >
            auto tDatabase = tHMR.get_database();

            // manually select output pattern
            tDatabase->set_activation_pattern( 0 );

            tHMR.perform_initial_refinement( 0 );

            tDatabase->get_background_mesh()->save_to_vtk("Bachgroundmesh_0_initial.vtk");

            // manually select output pattern
            tDatabase->set_activation_pattern( 1 );

            // refine the last element three times
            // fixme: change this to 2
            for( uint tLevel = 0; tLevel < 2; ++tLevel )
            {
                tDatabase->get_background_mesh()->get_element( 0 )->put_on_refinement_queue();

                // manually refine, do not reset pattern
                tDatabase->get_background_mesh()->perform_refinement( 1 );
            }

            tDatabase->get_background_mesh()->save_to_vtk("Bachgroundmesh_1_initial.vtk");



            tDatabase->unite_patterns( 0, 1, 2 );

            tDatabase->get_background_mesh()->save_to_vtk("Bachgroundmesh_2_initial.vtk");

            tDatabase->update_bspline_meshes();
            tDatabase->update_lagrange_meshes();
            // calculate T-Matrices etc
            tDatabase->finalize();

            std::shared_ptr< hmr::Interpolation_Mesh_HMR > tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );
            std::shared_ptr< moris::hmr::Integration_Mesh_HMR > tIntegrationMesh = tHMR.create_integration_mesh( 1, 2,*tInterpMesh );

            mtk::Mesh_Manager tMesh;
            uint tMeshIndex = tMesh.register_mesh_pair( tInterpMesh.get(), tIntegrationMesh.get() );

            //--------------------------------------------------------------------------------------------------

            // input parameters for the circle LS
            moris::Cell< real > tCircleInputs = {{0},{0},{0.9}};
            //------------------------------------------------------------------------------

            Ge_Factory tFactory;
            std::shared_ptr< Geometry > tGeom = tFactory.set_geometry_type( GeomType::ANALYTIC );

            tGeom->set_my_mesh( &tMesh );
            tGeom->set_my_constants(tCircleInputs);

            tGeom->set_analytical_function(AnalyticType::CIRCLE);
            tGeom->set_analytical_function_dphi_dx(AnalyticType::CIRCLE);

            GE_Core tGeometryEngine;
            moris_index tMyGeomIndex = tGeometryEngine.set_geometry( tGeom );

            uint tNumOfIPNodes = tIntegrationMesh->get_num_nodes();

            Matrix< DDRMat > tFieldData( tNumOfIPNodes,1, 0.0);
            for( uint n=0; n<tNumOfIPNodes; n++)
            {
                tFieldData(n)        = tGeometryEngine.get_field_vals( tMyGeomIndex, n )( 0 );     //FIXME
            }

            print(tFieldData, "tFieldData");

            tHMR.flag_surface_elements( tFieldData );

            tHMR.perform_refinement( moris::hmr::RefinementMode::SIMPLE, 1 );
//            tHMR.update_refinement_pattern( 1 );

            moris::Cell< moris::hmr::BSpline_Mesh_Base* > tBSplineMeshes1;

            // create factory
            moris::hmr::Factory tFactory_HMR;

            moris::hmr::Lagrange_Mesh_Base* tLagrangeMesh =  tFactory_HMR.create_lagrange_mesh( &tParameters,
                                                                                                tDatabase->get_background_mesh(),
                                                                                                 tBSplineMeshes1,
                                                                                                 1,
                                                                                                 1 );

             // output to exodus
             STK * tSTK = tLagrangeMesh->create_stk_object(0);
             tSTK->save_to_file( "GE_HMR_Mesh.g");
             delete tSTK;


//        // create first order Lagrange mesh
//        moris::hmr::Lagrange_Mesh_Base* tLagrangeMesh_1 =  tFactory.create_lagrange_mesh( tParameters,
//                                                                                          tBackgroundMesh,
//                                                                                          tBSplineMeshes,
//                                                                                          0,
//                                                                                          1 );
//        // create first order Lagrange mesh
//        moris::hmr::Lagrange_Mesh_Base* tLagrangeMesh_2 =  tFactory.create_lagrange_mesh( tParameters,
//                                                                                          tBackgroundMesh,
//                                                                                          tBSplineMeshes,
//                                                                                          1,
//                                                                                          1 );
//
//        REQUIRE( tLagrangeMesh_1->get_number_of_nodes_on_proc()  == 43 );
//        REQUIRE( tLagrangeMesh_2->get_number_of_nodes_on_proc()  == 30 );
//
//
//        // output to exodus
//        STK * tSTK = tLagrangeMesh_2->create_stk_object(0);
//        tSTK->save_to_file( "cccccc.g");
//        delete tSTK;
//
//        // Check some basis coordinates of Lagrange mesh 1
//        const moris::real* tXYZ_1 = tLagrangeMesh_1->get_node_by_index( 2 )->get_xyz( );
//        REQUIRE( tXYZ_1[0]  == 0.0625 );    REQUIRE( tXYZ_1[1]  == 0.0625 );
//        const moris::real* tXYZ_2 = tLagrangeMesh_1->get_node_by_index( 13 )->get_xyz( );
//        REQUIRE( tXYZ_2[0]  == 0.25 );    REQUIRE( tXYZ_2[1]  == 0.25 );
//        const moris::real* tXYZ_3 = tLagrangeMesh_1->get_node_by_index( 15 )->get_xyz( );
//        REQUIRE( tXYZ_3[0]  == 0.375 );    REQUIRE( tXYZ_3[1]  == 0.125 );
//        const moris::real* tXYZ_4 = tLagrangeMesh_1->get_node_by_index( 36 )->get_xyz( );
//        REQUIRE( tXYZ_4[0]  == 0.75 );    REQUIRE( tXYZ_4[1]  == 0.75 );
//        const moris::real* tXYZ_14 = tLagrangeMesh_1->get_node_by_index( 41 )->get_xyz( );
//        REQUIRE( tXYZ_14[0]  == 0.75 );    REQUIRE( tXYZ_14[1]  == 1.0 );
//
////        // Check some basis coordinates of Lagrange mesh 2
////        const moris::real* tXYZ_5 = tLagrangeMesh_2->get_node_by_index( 2 )->get_xyz( );
////        REQUIRE( tXYZ_5[0]  == 0.25 );    REQUIRE( tXYZ_5[1]  == 0.25 );
////        const moris::real* tXYZ_6 = tLagrangeMesh_2->get_node_by_index( 13 )->get_xyz( );
////        REQUIRE( tXYZ_6[0]  == 0.75 );    REQUIRE( tXYZ_6[1]  == 0.5 );
////        const moris::real* tXYZ_7 = tLagrangeMesh_2->get_node_by_index( 15 )->get_xyz( );
////        REQUIRE( tXYZ_7[0]  == 0.25 );    REQUIRE( tXYZ_7[1]  == 0.75 );
////        const moris::real* tXYZ_8 = tLagrangeMesh_2->get_node_by_index( 36 )->get_xyz( );
////        REQUIRE( tXYZ_8[0]  == 0.875 );    REQUIRE( tXYZ_8[1]  == 1.0 );
//
//
//        // delete mesh
//        delete tLagrangeMesh_1;
//        delete tLagrangeMesh_2;
    }
    }
}

