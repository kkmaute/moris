/*
 * cl_GE_test_discrete.cpp
 *
 *  Created on: Apr 22, 2019
 *      Author: sonne
 */

#include "catch.hpp"

//------------------------------------------------------------------------------
// GE includes
#include "cl_GE_Core.hpp"
#include "cl_GE_Element.hpp"
#include "cl_GE_Factory.hpp"
#include "cl_GE_Node.hpp"

// HMR includes
#include "cl_HMR.hpp"
#include "cl_HMR_Field.hpp"
#include "cl_HMR_Mesh.hpp"
#include "cl_HMR_Parameters.hpp"
#include "fn_HMR_Exec_perform_mapping.hpp"

// linalg includes
#include "cl_Matrix.hpp"
#include "fn_all_true.hpp"
#include "fn_equal_to.hpp"
#include "linalg_typedefs.hpp"
#include "op_equal_equal.hpp"

// MTK includes
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Vertex.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_Tools.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Scalar_Field_Info.hpp"
//------------------------------------------------------------------------------

using namespace moris;
using namespace ge;

TEST_CASE("discrete_functionalities_test_01","[GE],[discrete_functionalities_hmr_mesh]")
{
    if(par_size()<=1)
    {
//        // discrete scalar level set values
//        Matrix< DDRMat > tLevelSetVals( 25,1 );
//        tLevelSetVals(0) =  1.64; tLevelSetVals(1) =   0.64;
//        tLevelSetVals(2) = -0.36; tLevelSetVals(3) =   0.64;
//        tLevelSetVals(4) =  0.89; tLevelSetVals(5) =  -0.11;
//        tLevelSetVals(6) = -0.11; tLevelSetVals(7) =   0.89;
//        tLevelSetVals(8) =  0.14; tLevelSetVals(9) =   1.64;
//        tLevelSetVals(10) = 0.64; tLevelSetVals(11) =  0.89;
//        tLevelSetVals(12) = 0.89; tLevelSetVals(13) = -0.11;
//        tLevelSetVals(14) = 0.14; tLevelSetVals(15) =  0.64;
//        tLevelSetVals(16) = 1.64; tLevelSetVals(17) = -0.11;
//        tLevelSetVals(18) = 0.89; tLevelSetVals(19) =  0.89;
//        tLevelSetVals(20) = 0.14; tLevelSetVals(21) =  1.64;
//        tLevelSetVals(22) = 0.89; tLevelSetVals(23) =  0.89;
//        tLevelSetVals(24) = 0.14;
////------------------------------------------------------------
//        uint tMeshOrder = 1;    // Legrange order = B-spline order = 2
//
//        hmr::ParameterList tParameters = hmr::create_hmr_parameter_list();
//        tParameters.set( "number_of_elements_per_dimension", "2, 2" );
//        tParameters.set( "bspline_orders", "2" );
//        tParameters.set( "lagrange_orders", "2" );
////        tParameters.set( "verbose", 0 );                    // displays debug info when turned on
//        tParameters.set( "domain_dimensions", "2, 2" );     // if no dimensions set, hmr assumes 1x1
//        tParameters.set( "domain_offset", "-1, -1" );       // offset so that the center node is at ( 0, 0 )
//
//        hmr::HMR tHMR( tParameters );
//        // create MTK mesh object and SDF field
//        std::shared_ptr< hmr::Mesh > tMesh = tHMR.create_mesh(tMeshOrder);
//        std::shared_ptr< hmr::Field > tTargetField = tMesh->create_field( "targetField", tMeshOrder );
//        std::shared_ptr< hmr::Field > tOutputField = tMesh->create_field( "outputField", tMeshOrder );
//
//        std::shared_ptr< moris::hmr::Interpolation_Mesh_HMR > tInterpolationMesh = tHMR.create_interpolation_mesh( tMeshOrder );
//        std::shared_ptr< moris::hmr::Integration_Mesh_HMR >   tIntegrationMesh   = tHMR.create_integration_mesh( tMeshOrder, tHMR.mParameters->get_lagrange_output_pattern(), *tInterpolationMesh.get() );
//
//        // place the pair in mesh manager
//        mtk::Mesh_Manager tMeshManager;
//        uint tMeshIndex = tMeshManager.register_mesh_pair(tInterpolationMesh.get(),tIntegrationMesh.get());
//
//        tTargetField->put_scalar_values_on_field(tLevelSetVals);
//
//        tHMR.finalize();
//
////------------------------------------------------------------
//
//        Ge_Factory tFactory;
//        std::shared_ptr< Geometry > tCircle = tFactory.set_geometry_type(GeomType::DISCRETE);
//        tCircle->set_my_mesh(&tMeshManager);
//        tCircle->set_my_output_field(tOutputField);
//        tCircle->set_my_target_field(tTargetField);
//
//        GE_Core tGeomEng;
//        moris_index tGeomIndex = tGeomEng.set_geometry( tCircle );
//
////print(tGeomEng.get_field_vals(0,15), "B-spline coefficient: ");
//
////        tHMR.save_to_exodus( "circle01.exo" );
    }
}
