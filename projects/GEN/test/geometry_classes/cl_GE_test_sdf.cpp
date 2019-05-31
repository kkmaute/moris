/*
 * cl_GE_test_sdf.cpp
 *
 *  Created on: May 19, 2019
 *      Author: sonne
 */

#include "catch.hpp"

// GE includes
#include "cl_GE_Core.hpp"
#include "cl_GE_Factory.hpp"

// HMR includes
#include "cl_HMR.hpp"
#include "cl_HMR_Field.hpp"
#include "cl_HMR_Mesh.hpp"
#include "cl_HMR_Parameters.hpp"
#include "fn_HMR_Exec_perform_mapping.hpp"

// LINALG includes
#include "cl_Matrix.hpp"
#include "fn_all_true.hpp"
#include "fn_equal_to.hpp"
#include "linalg_typedefs.hpp"
#include "op_equal_equal.hpp"

//------------------------------------------------------------------------------

using namespace moris;
using namespace ge;

TEST_CASE("sdf_functionalities_test","[GE],[sdf_functionalities]")
{
    if(par_size()<=1)
    {
        uint tLagrangeOrder = 2;
        hmr::ParameterList tParameters = hmr::create_hmr_parameter_list();
        tParameters.set( "number_of_elements_per_dimension", "10, 10, 10" );
        tParameters.set( "domain_dimensions",                "5.6, 2.6, 3.4" );
        tParameters.set( "domain_offset",                    "-4.9, 3.25, -1.7" );
        tParameters.set( "lagrange_orders",                  "2" );
        tParameters.set( "verbose", 1 ); // verbosity flag for debugging and other outputs

        hmr::HMR tHMR( tParameters );
        // create MTK mesh object and SDF field
        std::shared_ptr< hmr::Mesh > tMesh = tHMR.create_mesh(tLagrangeOrder);
        std::shared_ptr< hmr::Field > tField = tMesh->create_field( "SDF", 2);
        tHMR.finalize();

        std::shared_ptr< moris::hmr::Interpolation_Mesh_HMR > tInterpolationMesh = tHMR.create_interpolation_mesh(tLagrangeOrder, tHMR.mParameters->get_lagrange_output_pattern());
        std::shared_ptr< moris::hmr::Integration_Mesh_HMR >   tIntegrationMesh   = tHMR.create_integration_mesh(tLagrangeOrder, tHMR.mParameters->get_lagrange_output_pattern(),*tInterpolationMesh);

        // place the pair in mesh manager
        mtk::Mesh_Manager tMeshManager;
        uint tMeshIndex = tMeshManager.register_mesh_pair(tInterpolationMesh.get(),tIntegrationMesh.get());

        // get path for STL file to load
        std::string tObjectPath = "/projects/HMR/tutorials/bracket.obj";
        tObjectPath = std::getenv("MORISROOT") + tObjectPath;
    //------------------------------------------------------------------------------

        Ge_Factory tFactory;
        std::shared_ptr< Geometry > tSDF = tFactory.set_geometry_type(GeomType::SDF);

        tSDF->set_my_mesh( &tMeshManager );
        tSDF->add_hmr_field( tField );
        tSDF->initialize_sdf( tObjectPath, tInterpolationMesh );

        GE_Core tGeomEng;
        tGeomEng.set_geometry(tSDF);
    //------------------------------------------------------------------------------
    // CHECK if the node values returned are the same as the ones in the SDF field
    //------------------------------------------------------------------------------
        uint tNumNodes = tSDF->get_my_mesh()->get_interpolation_mesh( 0 )->get_num_nodes();
        Matrix< DDRMat > tScalarMat(tNumNodes,1, 0.0);
        for(uint n=0; n < tNumNodes; ++n)
        {
            Matrix< DDRMat > tTempMat = tGeomEng.get_field_vals(0,n);
            tScalarMat(n) = tTempMat(0,0);
        }

        bool tMatrixMatch = all_true( tField->get_node_values() == tScalarMat );
        CHECK( tMatrixMatch );

//        tHMR.save_to_exodus( "SDF_test.exo" );
    }
}
