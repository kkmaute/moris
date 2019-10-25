/*
 * cl_Fiber_Problem.cpp
 *
 *  Created on: Oct 22, 2019
 *      Author: sonne
 */


#include "catch.hpp"
#include "HDF5_Tools.hpp"

// HMR includes
#include "cl_HMR.hpp"
#include "cl_HMR_Database.hpp"
#include "cl_HMR_Field.hpp"
#include "HMR_Globals.hpp"

#include "fn_cylinder_with_end_caps.hpp"

using namespace moris;
namespace ge
{

TEST_CASE("fiber_problem_test_01","[GE],[fiber_test]")
{
/*
 * ------------------------------------------------------------------------------
 *                              STK Mesh from string
 * ------------------------------------------------------------------------------
 */
//    const std::string fileName = "generated:100x100x1";
//
//    // Declare scalar node field
//    moris::mtk::Scalar_Field_Info<DDRMat> tNodeField1;
//    std::string tFieldName1 = "fiber01";
//    tNodeField1.set_field_name(tFieldName1);
//    tNodeField1.set_field_entity_rank(EntityRank::NODE);
//
//    // Initialize field information container
//    moris::mtk::MtkFieldsInfo tFieldsInfo;
//    // Place the node field into the field info container
//    add_field_for_mesh_input(&tNodeField1,tFieldsInfo);
//    // Declare some supplementary fields
//    mtk::MtkMeshData tMeshData;
//    tMeshData.FieldsInfo = &tFieldsInfo;
//    // create mesh pair
//    mtk::Interpolation_Mesh* tInterpMesh1 = create_interpolation_mesh( MeshType::STK, fileName, &tMeshData );
//    mtk::Integration_Mesh*   tIntegMesh1  = create_integration_mesh_from_interpolation_mesh(MeshType::STK,tInterpMesh1);
//
//    // place the pair in mesh manager
//    mtk::Mesh_Manager tMeshManager;
//    uint tMeshIndex = tMeshManager.register_mesh_pair(tInterpMesh1,tIntegMesh1);

/*
 * ------------------------------------------------------------------------------
 *                      HMR Mesh Using Parameter Options
 * ------------------------------------------------------------------------------
 */
        uint tLagrangeMeshIndex = 0;
        //  HMR Parameters setup
        hmr::ParameterList tParameters = hmr::create_hmr_parameter_list();

//        tParameters.set( "number_of_elements_per_dimension", "100, 100, 2" );
//        tParameters.set( "domain_dimensions",                "160, 80, 8" );
//        tParameters.set( "domain_offset",                    "-80, -40, -4" );
        tParameters.set( "number_of_elements_per_dimension", "100, 100, 2" );
        tParameters.set( "domain_dimensions",                "80, 20, 1" );
        tParameters.set( "domain_offset",                    "-40, 0, -0.5" );

        tParameters.set( "truncate_bsplines", 1 );
        tParameters.set( "lagrange_orders", "1" );
        tParameters.set( "lagrange_pattern", "0" );
        tParameters.set( "bspline_orders", "1" );
        tParameters.set( "bspline_pattern", "0" );

        tParameters.set( "lagrange_output_meshes", "0" );


        tParameters.set( "lagrange_to_bspline", "0" );

        tParameters.set( "use_multigrid", 0 );

        tParameters.set( "refinement_buffer", 1 );
        tParameters.set( "staircase_buffer", 1 );

        //  HMR Initialization
        moris::hmr::HMR tHMR( tParameters );

        // std::shared_ptr< Database >
        auto tDatabase = tHMR.get_database();

        tHMR.perform_initial_refinement( 0 );

        std::shared_ptr< moris::hmr::Mesh > tMesh01 = tHMR.create_mesh( tLagrangeMeshIndex );   // HMR Lagrange mesh

//        tDatabase->get_background_mesh()->save_to_vtk("Bachgroundmesh_2_initial.vtk");
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
moris::tic tTimer0;
        std::shared_ptr< hmr::Field > tField = tMesh01->create_field( "fibers", tLagrangeMeshIndex);
moris::tic tTimer1;
        tField->evaluate_scalar_function( moris::ge::cylinderWithEndCaps );
real tElapsedTime1 = tTimer1.toc<moris::chronos::milliseconds>().wall;
tElapsedTime1 /= 1000;
std::cout<<"==============================================="<< std::endl;
std::cout<<"Total to evaluate scalar function once: "<< tElapsedTime1 << std::endl << std::endl;
std::cout<<"==============================================="<< std::endl;
        for( uint k=0; k<2; ++k )
        {
            tHMR.flag_surface_elements_on_working_pattern( tField );
            tHMR.perform_refinement_based_on_working_pattern( 0 );

            tField->evaluate_scalar_function( moris::ge::cylinderWithEndCaps );
        }
        tDatabase->get_background_mesh()->save_to_vtk("Backgroundmesh_initial.vtk");
//------------------------------------------------------------------------------

    tDatabase->update_bspline_meshes();
    tDatabase->update_lagrange_meshes();

    // calculate T-Matrices etc
    tDatabase->finalize();
//------------------------------------------------------------------------------
    tField->evaluate_scalar_function( moris::ge::cylinderWithEndCaps );
    tHMR.save_to_exodus( 0, "fibers.g" );

real tElapsedTime0 = tTimer0.toc<moris::chronos::milliseconds>().wall;
tElapsedTime0 /= 1000;
std::cout<<"==============================================="<< std::endl;
std::cout<<"Total time for evaluation: "<< tElapsedTime0 << std::endl << std::endl;
std::cout<<"==============================================="<< std::endl;
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//    for(uint k=0; k<2; ++k)
//    {
//        std::cout<<"refinement "<<k<<std::endl;
//
//        std::shared_ptr< hmr::Interpolation_Mesh_HMR > tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );
//        std::shared_ptr< hmr::Integration_Mesh_HMR > tIntegrationMesh = tHMR.create_integration_mesh( 1, 0 ,*tInterpMesh );
//
//        mtk::Mesh_Manager tMesh;
//        uint tMeshIndex = tMesh.register_mesh_pair( tInterpMesh.get(), tIntegrationMesh.get() );
//
//        Ge_Factory tFactory;
//        std::shared_ptr< Geometry > tGeom = tFactory.set_geometry_type( GeomType::ANALYTIC );
//
//        tGeom->set_my_mesh( &tMesh );
//
//        moris_index tSubIndex = tGeom->set_analytical_function( cylinderWithEndCaps );
//
//        GE_Core tGeometryEngine;
//        moris_index tMyGeomIndex = tGeometryEngine.set_geometry( tGeom );
//
//        Matrix< DDRMat > tFieldData = tGeometryEngine.get_field_vals( tMyGeomIndex, tSubIndex );
//
//        tHMR.based_on_field_put_elements_on_queue( tFieldData, tLagrangeMeshIndex );
//
//        tDatabase->get_background_mesh()->perform_refinement( 0 );
//        tHMR.update_refinement_pattern( 0 );
//        tDatabase->finalize();
//    }

//    moris::Cell< moris::hmr::BSpline_Mesh_Base* > tBSplineMeshes1;
//
//    // create factory
//    moris::hmr::Factory tFactory_HMR;
//
//    hmr::Parameters* tParams = tHMR.get_parameters();
//    moris::hmr::Lagrange_Mesh_Base* tLagrangeMesh =  tFactory_HMR.create_lagrange_mesh( tParams,
//                                                                                        tDatabase->get_background_mesh(),
//                                                                                         tBSplineMeshes1,
//                                                                                         tLagrangeMeshIndex,
//                                                                                         1 );
//
//     // output to exodus
//     hmr::STK * tSTK = tLagrangeMesh->create_stk_object(0);
//     tSTK->save_to_file( "fiber_cylinder.g");
//     delete tSTK;
//     delete tLagrangeMesh;

    //------------------------------------------------------------------------------

//    Ge_Factory tFactory;
//    std::shared_ptr< Geometry > tGeom = tFactory.set_geometry_type( GeomType::ANALYTIC );
//
//    tGeom->set_my_mesh( &tMeshManager );
//
//    moris_index tSubIndex = tGeom->set_analytical_function( fiber00 );
//    GE_Core tGeometryEngine;
//    moris_index tMyGeomIndex = tGeometryEngine.set_geometry( tGeom );
//
//    uint tNumOfIPNodes = tIntegMesh1->get_num_nodes();
//
//    Matrix< DDRMat > tFieldData = tGeometryEngine.get_field_vals( tMyGeomIndex, tSubIndex );

//    tInterpMesh1->add_mesh_field_real_scalar_data_loc_inds( tFieldName1, EntityRank::NODE, tFieldData );

//    std::string tOutputFile = "./ge_fiber00.exo";
//    tInterpMesh1->create_output_mesh(tOutputFile);
//    // clean up
//    delete tInterpMesh1;
//    delete tIntegMesh1;
}

}   // ge namespace
