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
        uint tLagrangeMeshIndex = 0;
        hmr::ParameterList tParameters = hmr::create_hmr_parameter_list();
        tParameters.set( "number_of_elements_per_dimension", "10, 10, 10" );
        tParameters.set( "domain_dimensions",                "5.6, 2.6, 3.4" );
        tParameters.set( "domain_offset",                    "-4.9, 3.25, -1.7" );

        tParameters.set( "lagrange_orders", "2" );
        tParameters.set( "lagrange_pattern", "0" );
        tParameters.set( "bspline_orders", "2" );
        tParameters.set( "bspline_pattern", "0" );

        tParameters.set( "lagrange_output_meshes", "0" );

        tParameters.set( "lagrange_to_bspline", "0" );


        hmr::HMR tHMR( tParameters );
        // create MTK mesh object and SDF field
        std::shared_ptr< hmr::Mesh > tMesh = tHMR.create_mesh(tLagrangeMeshIndex);
        std::shared_ptr< hmr::Field > tField = tMesh->create_field( "SDF", tLagrangeMeshIndex);
        tHMR.finalize();

        std::shared_ptr< moris::hmr::Interpolation_Mesh_HMR > tInterpolationMesh = tHMR.create_interpolation_mesh(tLagrangeMeshIndex);
        std::shared_ptr< moris::hmr::Integration_Mesh_HMR >   tIntegrationMesh   = tHMR.create_integration_mesh(2, 0,*tInterpolationMesh);

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
        moris_index tSubIndex = tSDF->add_hmr_field( tField );
        tSDF->initialize_sdf( tObjectPath, tInterpolationMesh );

        GE_Core tGeomEng;
        moris_index tMyGeomIndex = tGeomEng.set_geometry(tSDF);
    //------------------------------------------------------------------------------
    // CHECK if the node values returned are the same as the ones in the SDF field
    //------------------------------------------------------------------------------
        //fixme needs to be updated to use the PDVInfoObject implementation

        PDV_Info* tPDVInfo = tGeomEng.get_pdv_info_pointer( tMyGeomIndex );

        Matrix< DDRMat > tLSVals = tPDVInfo->get_field_vals( tSubIndex );           // phi


        bool tMatrixMatch = all_true( tField->get_node_values() == tLSVals );
        CHECK( tMatrixMatch );

//        tHMR.save_to_exodus( "SDF_test.exo" );
    }
}
