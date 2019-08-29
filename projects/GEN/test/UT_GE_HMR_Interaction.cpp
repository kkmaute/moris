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

using namespace moris;
using namespace hmr;

TEST_CASE("GE_HMR_Interaction","[moris],[GE],[GE_HMR_Interaction]")
{
    if(par_size() == 1)
    {
        for( moris::uint tOrder=1; tOrder<=1; tOrder++ )
        {
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

            tParameters.set_initial_refinement( 2 );

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
            for( uint tLevel = 0; tLevel < 1; ++tLevel )
            {
                tDatabase->get_background_mesh()->get_element( 0 )->put_on_refinement_queue();

                // manually refine, do not reset pattern
                tDatabase->get_background_mesh()->perform_refinement( 1 );
            }

            tDatabase->get_background_mesh()->save_to_vtk("Bachgroundmesh_1_initial.vtk");

            tDatabase->unite_patterns( 0, 1, 2);

            tDatabase->update_bspline_meshes();
            tDatabase->update_lagrange_meshes();
            // calculate T-Matrices etc
            tDatabase->finalize();


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

