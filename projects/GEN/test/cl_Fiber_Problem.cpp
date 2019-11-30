/*
 * cl_Fiber_Problem.cpp
 *
 *  Created on: Oct 22, 2019
 *      Author: sonne
 */

#include "catch.hpp"
#include "HDF5_Tools.hpp"

#include "typedefs.hpp"

// LINALG
#include "fn_norm.hpp"

// HMR includes
#include "cl_HMR.hpp"
#include "cl_HMR_Database.hpp"
#include "cl_HMR_Field.hpp"
#include "HMR_Globals.hpp"

#include "cl_XTK_Model.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"

#include "cl_GE_Geometry_Library.hpp"
#include "cl_GE_Core.hpp"
#include "cl_GE_Factory.hpp"

#include "../src/geometry/cl_GEN_Circle.hpp"
#include "../src/geometry/cl_GEN_Cylinder_With_End_Caps.hpp"
#include "../src/geometry/cl_GEN_Geom_Data.hpp"
#include "../src/geometry/cl_GEN_Geom_Field.hpp"
#include "../src/geometry/cl_GEN_Geometry.hpp"
#include "../src/geometry/cl_GEN_Multi_Geometry.hpp"
#include "../src/geomeng/cl_GEN_Geometry_Engine.hpp"

#include "cl_MTK_Vertex.hpp"    //MTK
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh.hpp"

#include "cl_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_Tools.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Scalar_Field_Info.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"

#include "cl_FEM_CM_Factory.hpp"              //FEM/INT/src
#include "cl_FEM_Constitutive_User_Defined_Info.hpp"      //FEM/INT/src
#include "cl_FEM_Element_Factory.hpp"          //FEM/INT/src
#include "cl_FEM_ElementProxy.hpp"             //FEM/INT/src
#include "cl_FEM_IWG_Factory.hpp"              //FEM/INT/src
#include "cl_FEM_IWG_User_Defined_Info.hpp"              //FEM/INT/src
#include "cl_FEM_Node_Base.hpp"                //FEM/INT/src
#include "cl_FEM_NodeProxy.hpp"                //FEM/INT/src
#include "cl_FEM_Property_User_Defined_Info.hpp"              //FEM/INT/src
#include "cl_FEM_Set_User_Info.hpp"              //FEM/INT/src

#include "cl_MDL_Model.hpp"

#include "cl_DLA_Solver_Factory.hpp"
#include "cl_DLA_Solver_Interface.hpp"

#include "cl_NLA_Nonlinear_Solver_Factory.hpp"
#include "cl_NLA_Nonlinear_Solver.hpp"
#include "cl_NLA_Nonlinear_Problem.hpp"
#include "cl_MSI_Solver_Interface.hpp"
#include "cl_MSI_Equation_Object.hpp"
#include "cl_MSI_Model_Solver_Interface.hpp"
#include "cl_DLA_Linear_Solver_Aztec.hpp"
#include "cl_DLA_Linear_Solver.hpp"

#include "cl_TSA_Time_Solver_Factory.hpp"
#include "cl_TSA_Monolithic_Time_Solver.hpp"
#include "cl_TSA_Time_Solver.hpp"

#include "linalg_typedefs.hpp"

using namespace moris;
namespace ge
{

Matrix< DDRMat > tConstValFunction2MatMDL( moris::Cell< Matrix< DDRMat > >         & aParameters,
                                           moris::Cell< fem::Field_Interpolator* > & aDofFI,
                                           moris::Cell< fem::Field_Interpolator* > & aDvFI,
                                           fem::Geometry_Interpolator              * aGeometryInterpolator )
        {
            return aParameters( 0 );
        }

Matrix< DDRMat > tConstValFunction( moris::Cell< Matrix< DDRMat > >         & aCoeff,
                                    moris::Cell< fem::Field_Interpolator* > & aDofFieldInterpolator,
                                    moris::Cell< fem::Field_Interpolator* > & aDvFieldInterpolator,
                                    fem::Geometry_Interpolator              * aGeometryInterpolator )
        {
            return aCoeff( 0 );
        }


real testShellFunc( const moris::Matrix< moris::DDRMat > & aPoint, const moris::Cell< moris::real > aConst )
{
    real thickness = 0.1;

    auto circleSrf = norm( aPoint ) - 1;

    auto shell = moris::ge::shell_tor( circleSrf, -thickness );

    return shell;
}

int user_defined_refinement(       hmr::Element             * aElement,
                             const Cell< Matrix< DDRMat > > & aElementLocalValues,
                                   hmr::ParameterList       & aParameters )
{
    int aDoRefine = -1;

    // abs variable field threshold
    real lsth = 0.0;

    // abs variable field bandwidth (absolute)
    real lsbwabs = 0.3;

    // maximum refinement level
    uint maxlevel = 4;

    // min refinement level
    uint minlevel = 0;

    // max refinement level along interface
    uint maxifcref = 4;

    // max refinement level in volume
    uint maxvolref = 0;

    // current refinement level of element
    uint curlevel = aElement->get_level();

    // refinement strategy
    if ( aElementLocalValues( 0 ).min() <= lsth + lsbwabs )
    {
        // for volume refinement
        if ( aElementLocalValues( 0 ).min() <= lsth - lsbwabs )
        {
            if( curlevel < maxvolref && curlevel < maxlevel )
            {
                aDoRefine = 1; // refine
            }
            else if ( curlevel ==  maxvolref || curlevel == minlevel )
            {
                aDoRefine = 0; // keep
            }
            else
            {
                aDoRefine = -1; // coarsen
            }
        }
        // for interface refinement
        else
        {
            if( curlevel < maxifcref && curlevel < maxlevel )
            {
                aDoRefine = 1; // refine
            }
            else if ( curlevel ==  maxifcref || curlevel == minlevel )
            {
                aDoRefine = 0; // keep
            }
            else
            {
                aDoRefine = -1; // coarsen
            }
        }
    }
    else
    {
        if( curlevel <  minlevel )
        {
            aDoRefine = 1; // refine
        }
        else if ( curlevel == minlevel )
        {
            aDoRefine = 0; // keep
        }
    }

    return aDoRefine;
}

TEST_CASE("fiber_problem_test_04", "[GE],[fiber_test_BCs]")
{
//    uint tLagrangeMeshIndex = 0;
//    //  HMR Parameters setup
//    hmr::ParameterList tParameters = hmr::create_hmr_parameter_list();
//
//    tParameters.set( "number_of_elements_per_dimension", "80, 40, 16" );
//    tParameters.set( "domain_dimensions",                "40, 20, 8" );
//    tParameters.set( "domain_offset",                    "-0, -0, -0" );
//
//    tParameters.set( "domain_sidesets",                    "1, 2, 3, 4, 5, 6" );
//
//    tParameters.set( "truncate_bsplines", 1 );
//    tParameters.set( "lagrange_orders", "1" );
//    tParameters.set( "lagrange_pattern", "0" );
//    tParameters.set( "bspline_orders", "1" );
//    tParameters.set( "bspline_pattern", "0" );
//
//    tParameters.set( "lagrange_output_meshes", "0" );
//    tParameters.set( "lagrange_input_meshes", "0" );
//
//
//    tParameters.set( "lagrange_to_bspline", "0" );
//
//    tParameters.set( "use_multigrid", 0 );
//
//    tParameters.set( "refinement_buffer", 1 );
//    tParameters.set( "staircase_buffer", 1 );
//
//    tParameters.insert( "initial_refinement", 0 );
//
//    //  HMR Initialization
//    moris::hmr::HMR tHMR( tParameters );
//
//    auto tDatabase = tHMR.get_database(); // std::shared_ptr< Database >
//
//    tHMR.perform_initial_refinement( 0 );
//    //------------------------------------------------------------------------------
//
//    tDatabase->update_bspline_meshes();
//    tDatabase->update_lagrange_meshes();
//
//    uint tNumberOfFibers      = 10;
//    uint tNumberOfRefinements = 1;
//    std::cout<<"-------------------------------------"<<std::endl;
//    std::cout<<"number of fibers being analyzed:  "<<tNumberOfFibers<<std::endl;
//    std::cout<<"number of refinements performed:  "<<tNumberOfRefinements<<std::endl;
//    std::cout<<"-------------------------------------"<<std::endl;
//    Cell< Matrix< DDRMat > > tFieldData( 1 );
//
//    std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );
//
//    moris::tic tTimer0;
//    for( uint k=0; k<tNumberOfRefinements; ++k )
//    {
//        moris::ge::GEN_CylinderWithEndCaps    tFibers( tNumberOfFibers );
//        moris::Cell<moris::ge::GEN_Geometry*> tGeometryVector = {&tFibers};
//
//        moris::ge::GEN_Phase_Table         tPhaseTable( tGeometryVector.size(),  Phase_Table_Structure::EXP_BASE_2 );
//        moris::ge::GEN_Geometry_Engine     tGENGeometryEngine( tGeometryVector,tPhaseTable,3 );
//
//        moris_index tMeshIndex = tGENGeometryEngine.set_mesh( tMesh );
//
//        tFieldData( 0 ) = tGENGeometryEngine.get_cylinder_vals( tMeshIndex, &tFibers, tNumberOfFibers );
//
//        tHMR.user_defined_flagging( user_defined_refinement, tFieldData, tParameters, 0 );
//
//        tHMR.perform_refinement_based_on_working_pattern( 0, false );
//    }
//
//    tHMR.finalize();
//
//    moris::ge::GEN_CylinderWithEndCaps    tFibers( tNumberOfFibers );
//    moris::Cell<moris::ge::GEN_Geometry*> tGeometryVector_temp = {&tFibers};
//
//    moris::ge::GEN_Phase_Table         tPhaseTable_temp( tGeometryVector_temp.size(),  Phase_Table_Structure::EXP_BASE_2 );
//    moris::ge::GEN_Geometry_Engine     tGENGeometryEngine_temp( tGeometryVector_temp, tPhaseTable_temp, 3 );
//
//    moris_index tMeshIndex = tGENGeometryEngine_temp.set_mesh( tMesh );
//
//    tFieldData( 0 ) = tGENGeometryEngine_temp.get_cylinder_vals( tMeshIndex, &tFibers, tNumberOfFibers );
//
//    std::shared_ptr< hmr::Interpolation_Mesh_HMR >      tInterpMesh      = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );
//    std::shared_ptr< moris::hmr::Integration_Mesh_HMR > tIntegrationMesh = tHMR.create_integration_mesh( 1, 0,*tInterpMesh );
//
//    mtk::Mesh_Manager tMesh1;
//
//    real tElapsedTime0 = tTimer0.toc<moris::chronos::milliseconds>().wall;
//    tElapsedTime0 /= 1000;
//    std::cout<<"==============================================="<< std::endl;
//    std::cout<<"Total time for evaluation: "<< tElapsedTime0 << std::endl << std::endl;
//    std::cout<<"==============================================="<< std::endl;
//    //------------------------------------------------------------------------------
//    moris::ge::GEN_Geom_Data tGeomData( tFieldData( 0 ) );  // fiber LS values as field data
//    //==============================
//    real tRadius  = 10;
//    real tXcenter = 0;
//    real tYcenter = 0;
//    moris::ge::Circle tHoleGeom( tRadius, tXcenter, tYcenter );// hole geometry
//
//    moris::Cell<moris::ge::GEN_Geometry*> tGeometryVector = { &tGeomData, &tHoleGeom };
////    moris::ge::Multi_Geometry tAssembly(tGeometryVector);    // contains all the geometries on the problem
//
//    size_t tModelDimension = 3;
//    moris::ge::GEN_Phase_Table      tPhaseTable( tGeometryVector.size(),  Phase_Table_Structure::EXP_BASE_2 );
//    moris::ge::GEN_Geometry_Engine  tGENGeometryEngine( tGeometryVector, tPhaseTable, tModelDimension );
//    xtk::Model                      tXTKModel( tModelDimension, tInterpMesh.get(), tGENGeometryEngine );
//    tXTKModel.mVerbose = false;
//
//    Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4};
//    tXTKModel.decompose(tDecompositionMethods);
//
////    //=============================== temporary ============================================
////    // output solution and meshes
////    xtk::Output_Options tOutputOptions1;
////    tOutputOptions1.mAddNodeSets = false;
////    tOutputOptions1.mAddSideSets = false;
////    tOutputOptions1.mAddClusters = false;
////
////    std::string tIntegSolFieldNameUX = "UX";
////    std::string tIntegSolFieldNameUY = "UY";
////    std::string tIntegSolFieldNameUZ = "UZ";
////    tOutputOptions1.mRealNodeExternalFieldNames = { tIntegSolFieldNameUX, tIntegSolFieldNameUY, tIntegSolFieldNameUZ };
////
////    moris::mtk::Integration_Mesh* tIntegMesh11 = tXTKModel.get_output_mesh(tOutputOptions1);
////
////    std::string tMeshOutputFile1 = "./fibersGeomCheck.e";
////    tIntegMesh11->create_output_mesh(tMeshOutputFile1);
////    delete tIntegMesh11;
////    //============================= end temporary ===========================================
////    std::string tMeshOutputFile1 = "./output_fibers_29Oct.e";
////    tIntegMesh11->create_output_mesh(tMeshOutputFile1);
//
//    tXTKModel.perform_basis_enrichment(EntityRank::BSPLINE_1,0);
//
//    xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
//    xtk::Enriched_Integration_Mesh   & tEnrIntegMesh = tXTKModel.get_enriched_integ_mesh();
//
//    // place the pair in mesh manager
//    mtk::Mesh_Manager tMeshManager;
//    tMeshManager.register_mesh_pair( &tEnrInterpMesh, &tEnrIntegMesh);
//    //------------------------------------------------------------------------------
//    // create the material properties
//    std::shared_ptr< fem::Property > tPropEModPlate = std::make_shared< fem::Property >();
////    tPropEModPlate->set_parameters( { {{ 1220000 }} } );   // 1.22*10^6 Pa
//    tPropEModPlate->set_parameters( { {{ 100 }} } );
//    tPropEModPlate->set_val_function( tConstValFunction );
//
//    std::shared_ptr< fem::Property > tPropNuPlate = std::make_shared< fem::Property >();
////    tPropNuPlate->set_parameters( { {{ 0.4 }} } );   // Poinson's ratio
//    tPropNuPlate->set_parameters( { {{ 0.0 }} } );
//    tPropNuPlate->set_val_function( tConstValFunction );
//
//
//    std::shared_ptr< fem::Property > tPropEModFibers = std::make_shared< fem::Property >();
////    tPropEModFibers->set_parameters( { {{ 1030230000 }} } );  // 1030.23*10^6 Pa
//    tPropEModFibers->set_parameters( { {{ 100 }} } );
//    tPropEModFibers->set_val_function( tConstValFunction );
//
//    std::shared_ptr< fem::Property > tPropNuFibers = std::make_shared< fem::Property >();
////    tPropNuFibers->set_parameters( { {{ 0.4 }} } );          // Poinson's ratio
//    tPropNuFibers->set_parameters( { {{ 0.0 }} } );          // Poinson's ratio
//    tPropNuFibers->set_val_function( tConstValFunction );
//
//    //==============================
//    // create the boundary conditions and loading properties
//    //==============================
//    // symmetry boundary conditions
//    std::shared_ptr< fem::Property > tPropDirichletUX = std::make_shared< fem::Property >();    // fix displacement to be zero in x
//    tPropDirichletUX->set_parameters( { {{ 0.0 }} } );
//    tPropDirichletUX->set_val_function( tConstValFunction );
//    tPropDirichletUX->set_dof_type( MSI::Dof_Type::UX );
//
//    std::shared_ptr< fem::Property > tPropDirichletUY = std::make_shared< fem::Property >();    // fix displacement to be zero in y
//    tPropDirichletUY->set_parameters( { {{ 0.0 }} } );
//    tPropDirichletUY->set_val_function( tConstValFunction );
//    tPropDirichletUY->set_dof_type( MSI::Dof_Type::UY );
//
//    std::shared_ptr< fem::Property > tPropDirichletUZ = std::make_shared< fem::Property >();    // fix displacement to be zero in z
//    tPropDirichletUZ->set_parameters( { {{ 0.0 }} } );
//    tPropDirichletUZ->set_val_function( tConstValFunction );
//    tPropDirichletUZ->set_dof_type( MSI::Dof_Type::UZ );
//
//    //==============================
//    // loading on free end
//    std::shared_ptr< fem::Property > tPropNeumann = std::make_shared< fem::Property >();
//    tPropNeumann->set_parameters( {{{ 1000000 } , { 0.0 }, {0.0}}} );  // 1.0 MPa traction applied in the x-direction
//    tPropNeumann->set_val_function( tConstValFunction );
//
//    //==============================
//    // define constitutive models
//    fem::CM_Factory tCMFactory;
//
//    std::shared_ptr< fem::Constitutive_Model > tCMPlate = tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
//    tCMPlate->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY , MSI::Dof_Type::UZ }} );
//    tCMPlate->set_properties( { tPropEModPlate, tPropNuPlate } );
//    tCMPlate->set_space_dim( 3 );
//
//    std::shared_ptr< fem::Constitutive_Model > tCMFibers = tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
//    tCMFibers->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY , MSI::Dof_Type::UZ }} );
//    tCMFibers->set_properties( { tPropEModFibers, tPropNuFibers } );
//    tCMFibers->set_space_dim( 3 );
//
//    //==============================
//    // define the IWGs
//    fem::IWG_Factory tIWGFactory;
//
//    std::shared_ptr< fem::IWG > tIWGPlate = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_BULK );
//    tIWGPlate->set_residual_dof_type( { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ } );
//    tIWGPlate->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }} );
//    tIWGPlate->set_constitutive_models( { tCMPlate }, mtk::Master_Slave::MASTER );
//
//    std::shared_ptr< fem::IWG > tIWGFibers = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_BULK );
//    tIWGFibers->set_residual_dof_type( { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ } );
//    tIWGFibers->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }} );
//    tIWGFibers->set_constitutive_models( { tCMFibers }, mtk::Master_Slave::MASTER );
//
//    std::shared_ptr< fem::IWG > tIWGDirichletFixedUx = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_DIRICHLET );
//    tIWGDirichletFixedUx->set_residual_dof_type( { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ } );
//    tIWGDirichletFixedUx->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }} );
//    tIWGDirichletFixedUx->set_constitutive_models( { tCMPlate }, mtk::Master_Slave::MASTER );
//    tIWGDirichletFixedUx->set_properties( { tPropDirichletUX });
//
//    std::shared_ptr< fem::IWG > tIWGDirichletFixedUy = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_DIRICHLET );
//    tIWGDirichletFixedUy->set_residual_dof_type( { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ } );
//    tIWGDirichletFixedUy->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }} );
//    tIWGDirichletFixedUy->set_constitutive_models( { tCMPlate }, mtk::Master_Slave::MASTER );
//    tIWGDirichletFixedUy->set_properties( { tPropDirichletUY });
//
//    std::shared_ptr< fem::IWG > tIWGDirichletFixedUz = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_DIRICHLET );
//    tIWGDirichletFixedUz->set_residual_dof_type( { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ } );
//    tIWGDirichletFixedUz->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }} );
//    tIWGDirichletFixedUz->set_constitutive_models( { tCMPlate }, mtk::Master_Slave::MASTER );
//    tIWGDirichletFixedUz->set_properties( { tPropDirichletUZ });
//
//    std::shared_ptr< fem::IWG > tIWGNeumann = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_NEUMANN );
//    tIWGNeumann->set_residual_dof_type( { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ } );
//    tIWGNeumann->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }} );
//    tIWGNeumann->set_properties( { tPropNeumann }, mtk::Master_Slave::MASTER );
//
//    //==============================
//    // create a list of active block-sets
//    std::string tInterfaceSideSetName    = tEnrIntegMesh.get_interface_side_set_name( 0, 0, 1 );
//    std::string tDblInterfaceSideSetName = tEnrIntegMesh.get_dbl_interface_side_set_name(0,1);
//
//    //==============================
//    // define set info
//    // fibers - phase 1
//    // hole   - phase 2
//    // plate  - phase 3
//
//    fem::Set_User_Info tBulkPlate;
//    tBulkPlate.set_mesh_index( tEnrIntegMesh.get_block_set_index("HMR_dummy_n_p0") );
//    tBulkPlate.set_set_type( fem::Element_Type::BULK );
//    tBulkPlate.set_IWGs( { tIWGPlate } );
//
//    fem::Set_User_Info tBulkFibers;
//    tBulkFibers.set_mesh_index( tEnrIntegMesh.get_block_set_index("HMR_dummy_c_p0") );
//    tBulkFibers.set_set_type( fem::Element_Type::BULK );
//    tBulkFibers.set_IWGs( { tIWGPlate } );
//
//    //==============================
//    // symmetry boundary conditions on side-set 4 ( fix Ux = 0 )
//    fem::Set_User_Info tSetDirichletFixed01;
//    tSetDirichletFixed01.set_mesh_index( tEnrIntegMesh.get_side_set_index("SideSet_4_n_p1") );
//    tSetDirichletFixed01.set_set_type( fem::Element_Type::SIDESET );
////    tSetDirichletFixed01.set_IWGs( { tIWGDirichletFixedUx } );
//    tSetDirichletFixed01.set_IWGs( { tIWGDirichletFixedUx, tIWGDirichletFixedUz } );
//
//    fem::Set_User_Info tSetDirichletFixed02;
//    tSetDirichletFixed02.set_mesh_index( tEnrIntegMesh.get_side_set_index("SideSet_4_c_p1") );
//    tSetDirichletFixed02.set_set_type( fem::Element_Type::SIDESET );
////    tSetDirichletFixed02.set_IWGs( { tIWGDirichletFixedUx } );
//    tSetDirichletFixed02.set_IWGs( { tIWGDirichletFixedUx, tIWGDirichletFixedUz } );
//
//    fem::Set_User_Info tSetDirichletFixed03;
//    tSetDirichletFixed03.set_mesh_index( tEnrIntegMesh.get_side_set_index("SideSet_4_c_p2") );
//    tSetDirichletFixed03.set_set_type( fem::Element_Type::SIDESET );
////    tSetDirichletFixed03.set_IWGs( { tIWGDirichletFixedUx } );
//    tSetDirichletFixed03.set_IWGs( { tIWGDirichletFixedUx, tIWGDirichletFixedUz } );
//
//    fem::Set_User_Info tSetDirichletFixed04;
//    tSetDirichletFixed04.set_mesh_index( tEnrIntegMesh.get_side_set_index("SideSet_4_c_p2") );
//    tSetDirichletFixed04.set_set_type( fem::Element_Type::SIDESET );
////    tSetDirichletFixed04.set_IWGs( { tIWGDirichletFixedUx } );
//    tSetDirichletFixed04.set_IWGs( { tIWGDirichletFixedUx, tIWGDirichletFixedUz } );
//
//    fem::Set_User_Info tSetDirichletFixed05;
//    tSetDirichletFixed05.set_mesh_index( tEnrIntegMesh.get_side_set_index("SideSet_4_c_p3") );
//    tSetDirichletFixed05.set_set_type( fem::Element_Type::SIDESET );
////    tSetDirichletFixed05.set_IWGs( { tIWGDirichletFixedUx } );
//    tSetDirichletFixed05.set_IWGs( { tIWGDirichletFixedUx, tIWGDirichletFixedUz } );
//
//    fem::Set_User_Info tSetDirichletFixed06;
//    tSetDirichletFixed06.set_mesh_index( tEnrIntegMesh.get_side_set_index("SideSet_4_c_p3") );
//    tSetDirichletFixed06.set_set_type( fem::Element_Type::SIDESET );
////    tSetDirichletFixed06.set_IWGs( { tIWGDirichletFixedUx } );
//    tSetDirichletFixed06.set_IWGs( { tIWGDirichletFixedUx, tIWGDirichletFixedUz } );
//
//    //==============================
//    // symmetry boundary conditions on side-set 1 ( fix Uy = 0 )
//    fem::Set_User_Info tSetDirichletFixed07;
//    tSetDirichletFixed07.set_mesh_index( tEnrIntegMesh.get_side_set_index("SideSet_1_n_p1") );
//    tSetDirichletFixed07.set_set_type( fem::Element_Type::SIDESET );
//    tSetDirichletFixed07.set_IWGs( { tIWGDirichletFixedUy } );
//
//    fem::Set_User_Info tSetDirichletFixed08;
//    tSetDirichletFixed08.set_mesh_index( tEnrIntegMesh.get_side_set_index("SideSet_1_c_p1") );
//    tSetDirichletFixed08.set_set_type( fem::Element_Type::SIDESET );
//    tSetDirichletFixed08.set_IWGs( { tIWGDirichletFixedUy } );
//
//    fem::Set_User_Info tSetDirichletFixed09;
//    tSetDirichletFixed09.set_mesh_index( tEnrIntegMesh.get_side_set_index("SideSet_1_n_p2") );
//    tSetDirichletFixed09.set_set_type( fem::Element_Type::SIDESET );
//    tSetDirichletFixed09.set_IWGs( { tIWGDirichletFixedUy } );
//
//    fem::Set_User_Info tSetDirichletFixed10;
//    tSetDirichletFixed10.set_mesh_index( tEnrIntegMesh.get_side_set_index("SideSet_1_c_p2") );
//    tSetDirichletFixed10.set_set_type( fem::Element_Type::SIDESET );
//    tSetDirichletFixed10.set_IWGs( { tIWGDirichletFixedUy } );
//
//    fem::Set_User_Info tSetDirichletFixed11;
//    tSetDirichletFixed11.set_mesh_index( tEnrIntegMesh.get_side_set_index("SideSet_1_n_p3") );
//    tSetDirichletFixed11.set_set_type( fem::Element_Type::SIDESET );
//    tSetDirichletFixed11.set_IWGs( { tIWGDirichletFixedUy } );
//
//    fem::Set_User_Info tSetDirichletFixed12;
//    tSetDirichletFixed12.set_mesh_index( tEnrIntegMesh.get_side_set_index("SideSet_1_c_p3") );
//    tSetDirichletFixed12.set_set_type( fem::Element_Type::SIDESET );
//    tSetDirichletFixed12.set_IWGs( { tIWGDirichletFixedUy } );
//
//    //==============================
//    // symmetry boundary conditions on side-set 5 ( fix Uz = 0 )
////    fem::Set_User_Info tSetDirichletFixed13;
////    tSetDirichletFixed13.set_mesh_index( tEnrIntegMesh.get_side_set_index("SideSet_5_n_p1") );
////    tSetDirichletFixed13.set_set_type( fem::Element_Type::SIDESET );
////    tSetDirichletFixed13.set_IWGs( { tIWGDirichletFixedUz } );
////
////    fem::Set_User_Info tSetDirichletFixed14;
////    tSetDirichletFixed14.set_mesh_index( tEnrIntegMesh.get_side_set_index("SideSet_5_c_p1") );
////    tSetDirichletFixed14.set_set_type( fem::Element_Type::SIDESET );
////    tSetDirichletFixed14.set_IWGs( { tIWGDirichletFixedUz } );
////
////    fem::Set_User_Info tSetDirichletFixed15;
////    tSetDirichletFixed15.set_mesh_index( tEnrIntegMesh.get_side_set_index("SideSet_5_n_p2") );
////    tSetDirichletFixed15.set_set_type( fem::Element_Type::SIDESET );
////    tSetDirichletFixed15.set_IWGs( { tIWGDirichletFixedUz } );
////
////    fem::Set_User_Info tSetDirichletFixed16;
////    tSetDirichletFixed16.set_mesh_index( tEnrIntegMesh.get_side_set_index("SideSet_5_c_p2") );
////    tSetDirichletFixed16.set_set_type( fem::Element_Type::SIDESET );
////    tSetDirichletFixed16.set_IWGs( { tIWGDirichletFixedUz } );
////
////    fem::Set_User_Info tSetDirichletFixed17;
////    tSetDirichletFixed17.set_mesh_index( tEnrIntegMesh.get_side_set_index("SideSet_5_n_p3") );
////    tSetDirichletFixed17.set_set_type( fem::Element_Type::SIDESET );
////    tSetDirichletFixed17.set_IWGs( { tIWGDirichletFixedUz } );
////
////    fem::Set_User_Info tSetDirichletFixed18;
////    tSetDirichletFixed18.set_mesh_index( tEnrIntegMesh.get_side_set_index("SideSet_5_c_p3") );
////    tSetDirichletFixed18.set_set_type( fem::Element_Type::SIDESET );
////    tSetDirichletFixed18.set_IWGs( { tIWGDirichletFixedUz } );
//
//    //==============================
//    // Neumann load on side-set 2
//    fem::Set_User_Info tSetNeumann01;
//    tSetNeumann01.set_mesh_index( tEnrIntegMesh.get_side_set_index("SideSet_2_n_p1") );
//    tSetNeumann01.set_set_type( fem::Element_Type::SIDESET );
//    tSetNeumann01.set_IWGs( { tIWGNeumann } );
//
//    fem::Set_User_Info tSetNeumann02;
//    tSetNeumann02.set_mesh_index( tEnrIntegMesh.get_side_set_index("SideSet_2_c_p1") );
//    tSetNeumann02.set_set_type( fem::Element_Type::SIDESET );
//    tSetNeumann02.set_IWGs( { tIWGNeumann } );
//
//    fem::Set_User_Info tSetNeumann03;
//    tSetNeumann03.set_mesh_index( tEnrIntegMesh.get_side_set_index("SideSet_2_n_p3") );
//    tSetNeumann03.set_set_type( fem::Element_Type::SIDESET );
//    tSetNeumann03.set_IWGs( { tIWGNeumann } );
//
//    fem::Set_User_Info tSetNeumann04;
//    tSetNeumann04.set_mesh_index( tEnrIntegMesh.get_side_set_index("SideSet_2_c_p3") );
//    tSetNeumann04.set_set_type( fem::Element_Type::SIDESET );
//    tSetNeumann04.set_IWGs( { tIWGNeumann } );
//
//    //==============================
//    // create a cell of set info
////    moris::Cell< fem::Set_User_Info > tSetInfo( 24 );
//    moris::Cell< fem::Set_User_Info > tSetInfo( 18 );
//    tSetInfo( 0 )  = tBulkPlate;
//    tSetInfo( 1 )  = tBulkFibers;
//    tSetInfo( 2 )  = tSetDirichletFixed01;
//    tSetInfo( 3 )  = tSetDirichletFixed02;
//    tSetInfo( 4 )  = tSetDirichletFixed03;
//    tSetInfo( 5 )  = tSetDirichletFixed04;
//    tSetInfo( 6 )  = tSetDirichletFixed05;
//    tSetInfo( 7 )  = tSetDirichletFixed06;
//    tSetInfo( 8 )  = tSetDirichletFixed07;
//    tSetInfo( 9 )  = tSetDirichletFixed08;
//    tSetInfo( 10 ) = tSetDirichletFixed09;
//    tSetInfo( 11 ) = tSetDirichletFixed10;
//    tSetInfo( 12 ) = tSetDirichletFixed11;
//    tSetInfo( 13 ) = tSetDirichletFixed12;
//    tSetInfo( 14 ) = tSetNeumann01;
//    tSetInfo( 15 ) = tSetNeumann02;
//    tSetInfo( 16 ) = tSetNeumann03;
//    tSetInfo( 17 ) = tSetNeumann04;
////    tSetInfo( 14 ) = tSetDirichletFixed13;
////    tSetInfo( 15 ) = tSetDirichletFixed14;
////    tSetInfo( 16 ) = tSetDirichletFixed15;
////    tSetInfo( 17 ) = tSetDirichletFixed16;
////    tSetInfo( 18 ) = tSetDirichletFixed17;
////    tSetInfo( 19 ) = tSetDirichletFixed18;
////    tSetInfo( 20 ) = tSetNeumann01;
////    tSetInfo( 21 ) = tSetNeumann02;
////    tSetInfo( 22 ) = tSetNeumann03;
////    tSetInfo( 23 ) = tSetNeumann04;
//
//    //------------------------------------------------------------------------------
//    // create model
//    mdl::Model * tModel = new mdl::Model( &tMeshManager,
//                                          0,
//                                          tSetInfo,
//                                          0,
//                                          false );
//
//    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//    // STEP 1: create linear solver and algorithm
//    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//
//    dla::Solver_Factory  tSolFactory;
//    std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( SolverType::AZTEC_IMPL );
//    //            std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( SolverType::PETSC );
//
//    tLinearSolverAlgorithm->set_param("AZ_diagnostics") = AZ_none;
//    tLinearSolverAlgorithm->set_param("AZ_output") = AZ_none;
//    tLinearSolverAlgorithm->set_param("AZ_max_iter") = 1000;
//    tLinearSolverAlgorithm->set_param("AZ_solver") = AZ_gmres;
//    tLinearSolverAlgorithm->set_param("AZ_subdomain_solve") = AZ_ilu;
//    tLinearSolverAlgorithm->set_param("AZ_graph_fill") = 5;
//
//    //    tLinearSolverAlgorithm->set_param("Use_ML_Prec") = true;
//
//    dla::Linear_Solver tLinSolver;
//
//    tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );
//
//    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//    // STEP 2: create nonlinear solver and algorithm
//    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//    NLA::Nonlinear_Problem * tNonlinearProblem =  new NLA::Nonlinear_Problem( tModel->get_solver_interface(), 0, true, MapType::Epetra );
////        NLA::Nonlinear_Problem * tNonlinearProblem =  new NLA::Nonlinear_Problem( tModel->get_solver_interface(), 0, true, MapType::Petsc );
//
//    NLA::Nonlinear_Solver_Factory tNonlinFactory;
//    std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );
//
//    tNonlinearSolverAlgorithm->set_param("NLA_max_iter")   = 100;
////        tNonlinearSolverAlgorithm->set_param("NLA_hard_break") = false;
////        tNonlinearSolverAlgorithm->set_param("NLA_max_lin_solver_restarts") = 2;
////        tNonlinearSolverAlgorithm->set_param("NLA_rebuild_jacobian") = true;
//
//    tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );
//
//    NLA::Nonlinear_Solver tNonlinearSolver;
//
//    tNonlinearSolver.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );
//
//    tNonlinearSolver.solve( tNonlinearProblem );
//
//    std::cout<<" Solution Vector "<<std::endl;
//    tNonlinearProblem->get_full_vector()->print();
//
//    // output solution and meshes
//    xtk::Output_Options tOutputOptions;
//    tOutputOptions.mAddNodeSets = false;
//    tOutputOptions.mAddSideSets = true;
//    tOutputOptions.mAddClusters = false;
//
//    // add solution field to integration mesh
//    std::string tIntegSolFieldNameUX = "UX";
//    std::string tIntegSolFieldNameUY = "UY";
//    std::string tIntegSolFieldNameUZ = "UZ";
//    tOutputOptions.mRealNodeExternalFieldNames = {tIntegSolFieldNameUX, tIntegSolFieldNameUY, tIntegSolFieldNameUZ};
//
//    moris::mtk::Integration_Mesh* tIntegMesh1 = tXTKModel.get_output_mesh(tOutputOptions);
//
//    // Write to Integration mesh for visualization
//    Matrix<DDRMat> tIntegSolUX = tModel->get_solution_for_integration_mesh_output( MSI::Dof_Type::UX );
//    Matrix<DDRMat> tIntegSolUY = tModel->get_solution_for_integration_mesh_output( MSI::Dof_Type::UY );
//    Matrix<DDRMat> tIntegSolUZ = tModel->get_solution_for_integration_mesh_output( MSI::Dof_Type::UZ );
//
//
//    Matrix<DDRMat> tSTKIntegSolUX(tIntegMesh1->get_num_entities(EntityRank::NODE),1);
//    Matrix<DDRMat> tSTKIntegSolUY(tIntegMesh1->get_num_entities(EntityRank::NODE),1);
//    Matrix<DDRMat> tSTKIntegSolUZ(tIntegMesh1->get_num_entities(EntityRank::NODE),1);
//
//    for(moris::uint i = 0; i < tIntegMesh1->get_num_entities(EntityRank::NODE); i++)
//    {
//        moris::moris_id tID = tIntegMesh1->get_glb_entity_id_from_entity_loc_index(i,EntityRank::NODE);
//        tSTKIntegSolUX(i) = tIntegSolUX(tEnrIntegMesh.get_loc_entity_ind_from_entity_glb_id(tID,EntityRank::NODE));
//        tSTKIntegSolUY(i) = tIntegSolUY(tEnrIntegMesh.get_loc_entity_ind_from_entity_glb_id(tID,EntityRank::NODE));
//        tSTKIntegSolUZ(i) = tIntegSolUZ(tEnrIntegMesh.get_loc_entity_ind_from_entity_glb_id(tID,EntityRank::NODE));
//    }
//
//    // add solution field to integration mesh
//    tIntegMesh1->add_mesh_field_real_scalar_data_loc_inds(tIntegSolFieldNameUX,EntityRank::NODE,tSTKIntegSolUX);
//    tIntegMesh1->add_mesh_field_real_scalar_data_loc_inds(tIntegSolFieldNameUY,EntityRank::NODE,tSTKIntegSolUY);
//    tIntegMesh1->add_mesh_field_real_scalar_data_loc_inds(tIntegSolFieldNameUZ,EntityRank::NODE,tSTKIntegSolUZ);
//
//    //------------------------------------------------------------------------------
//    // output solution to mesh
//    std::string tMeshOutputFile = "fiber_symm_BCs.e";
//
//    tIntegMesh1->create_output_mesh(tMeshOutputFile);
//    delete tIntegMesh1;
}

}   // ge namespace
