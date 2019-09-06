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

TEST_CASE("HMR_IO","[moris],[hmr],[HMR_IO]")
{
    if(par_size() == 1)
    {
        for( moris::uint tOrder=1; tOrder<=3; tOrder++ )
        {
            uint tLagrangeMeshIndex = 0;

            // create settings object
            moris::hmr::Parameters tParameters;

            tParameters.set_number_of_elements_per_dimension( { {4}, {4} } );

            tParameters.set_domain_dimensions({ {1}, {1} });
            tParameters.set_domain_offset({ {-0.5}, {-0.5} });

            // B-Spline truncation is turned on by default.
            // It is recommended to leave this setting as is.
            tParameters.set_bspline_truncation( true );

            tParameters.set_lagrange_orders  ( { {tOrder} });
            tParameters.set_lagrange_patterns({ {2} });

            tParameters.set_bspline_orders   ( { {tOrder}, {tOrder}} );
            tParameters.set_bspline_patterns ( { {0}, {1}} );

            tParameters.set_staircase_buffer( 3 );
            tParameters.set_refinement_buffer( 3 );

            tParameters.set_initial_refinement( 1 );

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

            // manually select output pattern
            tDatabase->set_activation_pattern( 1 );

            // refine the last element three times
            // fixme: change this to 2
            for( uint tLevel = 0; tLevel < 3; ++tLevel )
            {
                tDatabase->get_background_mesh()->get_element( 0 )->put_on_refinement_queue();

                // manually refine, do not reset pattern
                tDatabase->get_background_mesh()->perform_refinement( 1 );
            }

            tDatabase->unite_patterns( 0, 1, 2 );

            tDatabase->update_bspline_meshes();
            tDatabase->update_lagrange_meshes();
            // calculate T-Matrices etc
            tDatabase->finalize();

            tHMR.save_to_hdf5( "Mesh_Data_test.hdf5", 0 );

            moris::hmr::HMR tHMR_Input( "Mesh_Data_test.hdf5" );

            // empty container for B-Spline meshes
            moris::Cell< moris::hmr::BSpline_Mesh_Base* > tBSplineMeshes;

            // create factory
            moris::hmr::Factory tFactory;

            // create first order Lagrange mesh
            moris::hmr::Lagrange_Mesh_Base* tLagrangeMesh =  tFactory.create_lagrange_mesh( tHMR_Input.get_database()->get_parameters(),
                                                                                              tHMR_Input.get_database()->get_background_mesh(),
                                                                                              tBSplineMeshes,
                                                                                              2,
                                                                                              1 );

            REQUIRE( tLagrangeMesh->get_number_of_nodes_on_proc()  == 116 );

            // Check some basis coordinates of Lagrange mesh 1
            const moris::real* tXYZ_1 = tLagrangeMesh->get_node_by_index( 2 )->get_xyz( );
            REQUIRE( tXYZ_1[0]  == -0.46875 );    REQUIRE( tXYZ_1[1]  == -0.46875 );
            const moris::real* tXYZ_2 = tLagrangeMesh->get_node_by_index( 71 )->get_xyz( );
            REQUIRE( tXYZ_2[0]  == 0.0 );    REQUIRE( tXYZ_2[1]  == 0.0 );
            const moris::real* tXYZ_3 = tLagrangeMesh->get_node_by_index( 109 )->get_xyz( );
            REQUIRE( tXYZ_3[0]  == 0.25 );    REQUIRE( tXYZ_3[1]  == 0.375 );
    }
    }
}

