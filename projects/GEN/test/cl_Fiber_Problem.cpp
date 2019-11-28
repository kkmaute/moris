///*
// * cl_Fiber_Problem.cpp
// *
// *  Created on: Oct 22, 2019
// *      Author: sonne
// */
//
//#include "catch.hpp"
//#include "HDF5_Tools.hpp"
//
//#include "typedefs.hpp"
//
//// LINALG
//#include "fn_norm.hpp"
//
//// HMR includes
//#include "cl_HMR.hpp"
//#include "cl_HMR_Database.hpp"
//#include "cl_HMR_Field.hpp"
//#include "HMR_Globals.hpp"
//
//#include "cl_XTK_Model.hpp"
//#include "cl_XTK_Enriched_Integration_Mesh.hpp"
//#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
//
//#include "cl_GE_Geometry_Library.hpp"
//#include "cl_GE_Core.hpp"
//#include "cl_GE_Factory.hpp"
//
//#include "../src/new/geometry/cl_GEN_Cylinder_With_End_Caps.hpp"
//#include "../src/new/geometry/cl_GEN_Geom_Data.hpp"
//#include "../src/new/geometry/cl_GEN_Geom_Field.hpp"
//#include "../src/new/geometry/cl_GEN_Geometry.hpp"
//#include "../src/new/geomeng/cl_GEN_Geometry_Engine.hpp"
//
//#include "cl_MTK_Vertex.hpp"    //MTK
//#include "cl_MTK_Cell.hpp"
//#include "cl_MTK_Enums.hpp"
//#include "cl_MTK_Mesh.hpp"
//
//#include "cl_Mesh_Factory.hpp"
//#include "cl_MTK_Mesh_Tools.hpp"
//#include "cl_MTK_Mesh_Data_Input.hpp"
//#include "cl_MTK_Scalar_Field_Info.hpp"
//#include "cl_MTK_Mesh.hpp"
//#include "cl_MTK_Mesh_Manager.hpp"
//#include "cl_MTK_Interpolation_Mesh.hpp"
//#include "cl_MTK_Integration_Mesh.hpp"
//
//#include "cl_FEM_NodeProxy.hpp"                //FEM/INT/src
//#include "cl_FEM_ElementProxy.hpp"             //FEM/INT/src
//#include "cl_FEM_Node_Base.hpp"                //FEM/INT/src
//#include "cl_FEM_Element_Factory.hpp"          //FEM/INT/src
//#include "cl_FEM_IWG_Factory.hpp"              //FEM/INT/src
//#include "cl_FEM_Property_User_Defined_Info.hpp"              //FEM/INT/src
//#include "cl_FEM_IWG_User_Defined_Info.hpp"              //FEM/INT/src
//#include "cl_FEM_Constitutive_User_Defined_Info.hpp"      //FEM/INT/src
//
//#include "cl_MDL_Model.hpp"
//
//#include "cl_DLA_Solver_Factory.hpp"
//#include "cl_DLA_Solver_Interface.hpp"
//
//#include "cl_NLA_Nonlinear_Solver_Factory.hpp"
//#include "cl_NLA_Nonlinear_Solver.hpp"
//#include "cl_NLA_Nonlinear_Problem.hpp"
//#include "cl_MSI_Solver_Interface.hpp"
//#include "cl_MSI_Equation_Object.hpp"
//#include "cl_MSI_Model_Solver_Interface.hpp"
//#include "cl_DLA_Linear_Solver_Aztec.hpp"
//#include "cl_DLA_Linear_Solver.hpp"
//
//#include "cl_TSA_Time_Solver_Factory.hpp"
//#include "cl_TSA_Monolithic_Time_Solver.hpp"
//#include "cl_TSA_Time_Solver.hpp"
//
//using namespace moris;
//namespace ge
//{
//
//Matrix< DDRMat > tConstValFunction2MatMDL( moris::Cell< Matrix< DDRMat > >         & aParameters,
//                                           moris::Cell< fem::Field_Interpolator* > & aDofFI,
//                                           moris::Cell< fem::Field_Interpolator* > & aDvFI,
//                                           fem::Geometry_Interpolator              * aGeometryInterpolator )
//        {
//            return aParameters( 0 );
//        }
//
//Matrix< DDRMat > tConstValFunction( moris::Cell< Matrix< DDRMat > >         & aCoeff,
//                                    moris::Cell< fem::Field_Interpolator* > & aDofFieldInterpolator,
//                                    moris::Cell< fem::Field_Interpolator* > & aDvFieldInterpolator,
//                                    fem::Geometry_Interpolator              * aGeometryInterpolator )
//        {
//            return aCoeff( 0 );
//        }
//
//
//real testShellFunc( const moris::Matrix< moris::DDRMat > & aPoint, const moris::Cell< moris::real > aConst )
//{
//    real thickness = 0.1;
//
//    auto circleSrf = norm( aPoint ) - 1;
//
//    auto shell = moris::ge::shell_tor( circleSrf, -thickness );
//
//    return shell;
//}
//
//int user_defined_refinement(       hmr::Element             * aElement,
//                             const Cell< Matrix< DDRMat > > & aElementLocalValues,
//                                   hmr::ParameterList       & aParameters )
//{
//    int aDoRefine = -1;
//
//    // abs variable field threshold
//    real lsth = 0.0;
//
//    // abs variable field bandwidth (absolute)
//    real lsbwabs = 0.3;
//
//    // maximum refinement level
//    uint maxlevel = 4;
//
//    // min refinement level
//    uint minlevel = 0;
//
//    // max refinement level along interface
//    uint maxifcref = 4;
//
//    // max refinement level in volume
//    uint maxvolref = 0;
//
//    // current refinement level of element
//    uint curlevel = aElement->get_level();
//
//    // refinement strategy
//    if ( aElementLocalValues( 0 ).min() <= lsth + lsbwabs )
//    {
//        // for volume refinement
//        if ( aElementLocalValues( 0 ).min() <= lsth - lsbwabs )
//        {
//            if( curlevel < maxvolref && curlevel < maxlevel )
//            {
//                aDoRefine = 1; // refine
//            }
//            else if ( curlevel ==  maxvolref || curlevel == minlevel )
//            {
//                aDoRefine = 0; // keep
//            }
//            else
//            {
//                aDoRefine = -1; // coarsen
//            }
//        }
//        // for interface refinement
//        else
//        {
//            if( curlevel < maxifcref && curlevel < maxlevel )
//            {
//                aDoRefine = 1; // refine
//            }
//            else if ( curlevel ==  maxifcref || curlevel == minlevel )
//            {
//                aDoRefine = 0; // keep
//            }
//            else
//            {
//                aDoRefine = -1; // coarsen
//            }
//        }
//    }
//    else
//    {
//        if( curlevel <  minlevel )
//        {
//            aDoRefine = 1; // refine
//        }
//        else if ( curlevel == minlevel )
//        {
//            aDoRefine = 0; // keep
//        }
//    }
//
//    return aDoRefine;
//}
//
//TEST_CASE("fiber_problem_test_01","[GE],[fiber_test_heat]")
//{
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
//    uint tNumberOfFibers = 1;
//    std::cout<<"-------------------------------------"<<std::endl;
//    std::cout<<"number of fibers being analyzed:  "<<tNumberOfFibers<<std::endl;
//    std::cout<<"-------------------------------------"<<std::endl;
//    Cell< Matrix< DDRMat > > tFieldData( 1 );
//
//    std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );
//
//    moris::tic tTimer0;
//    for( uint k=0; k<2; ++k )
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
//    std::shared_ptr< hmr::Interpolation_Mesh_HMR > tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );
//    std::shared_ptr< moris::hmr::Integration_Mesh_HMR > tIntegrationMesh = tHMR.create_integration_mesh( 1, 0,*tInterpMesh );
//
//    mtk::Mesh_Manager tMesh1;
//
//    real tElapsedTime0 = tTimer0.toc<moris::chronos::milliseconds>().wall;
//    tElapsedTime0 /= 1000;
//    std::cout<<"==============================================="<< std::endl;
//    std::cout<<"Total time for evaluation: "<< tElapsedTime0 << std::endl << std::endl;
//    std::cout<<"==============================================="<< std::endl;
//
//    moris::ge::GEN_Geom_Data tGeomData( tFieldData( 0 ) );
//    moris::Cell<moris::ge::GEN_Geometry*> tGeometryVector = {&tGeomData};
//
//    size_t tModelDimension = 3;
//    moris::ge::GEN_Phase_Table      tPhaseTable( tGeometryVector.size(),  Phase_Table_Structure::EXP_BASE_2 );
//    moris::ge::GEN_Geometry_Engine  tGENGeometryEngine( tGeometryVector,tPhaseTable,tModelDimension );
//    xtk::Model                      tXTKModel( tModelDimension,tInterpMesh.get(),tGENGeometryEngine );
//    tXTKModel.mVerbose = false;
//
//    Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4};
//    tXTKModel.decompose(tDecompositionMethods);
//
//    // output solution and meshes
//
//    xtk::Output_Options tOutputOptions1;
//    tOutputOptions1.mAddNodeSets = false;
//    tOutputOptions1.mAddSideSets = false;
//    tOutputOptions1.mAddClusters = false;
//
//    std::string tIntegSolFieldName1 = "solution";
//    tOutputOptions1.mRealNodeExternalFieldNames = {tIntegSolFieldName1};
//
//    moris::mtk::Integration_Mesh* tIntegMesh11 = tXTKModel.get_output_mesh(tOutputOptions1);
//
////    std::string tMeshOutputFile1 = "./output_fibers_29Oct.e";
////    tIntegMesh11->create_output_mesh(tMeshOutputFile1);
//
//    tXTKModel.perform_basis_enrichment(EntityRank::BSPLINE_1,0);
//
//    xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
//    xtk::Enriched_Integration_Mesh   & tEnrIntegMesh = tXTKModel.get_enriched_integ_mesh();
//
//    //==============================
//
//    // place the pair in mesh manager
//    mtk::Mesh_Manager tMeshManager;
//    tMeshManager.register_mesh_pair( &tEnrInterpMesh, &tEnrIntegMesh);
//
//    // create IWG user defined info
//    Cell< Cell< fem::IWG_User_Defined_Info > > tIWGUserDefinedInfo( 10 );
//    tIWGUserDefinedInfo( 0 ).resize( 1 );
//    tIWGUserDefinedInfo( 0 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_BULK,
//                                                                { MSI::Dof_Type::TEMP },
//                                                                {{ MSI::Dof_Type::TEMP }},
//                                                                {fem::Property_Type::TEMP_LOAD },
//                                                                { fem::Constitutive_Type::DIFF_LIN_ISO } );
//    tIWGUserDefinedInfo( 1 ).resize( 1 );
//    tIWGUserDefinedInfo( 1 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_BULK,
//                                                                { MSI::Dof_Type::TEMP },
//                                                                {{ MSI::Dof_Type::TEMP }},
//                                                                {fem::Property_Type::TEMP_LOAD },
//                                                                { fem::Constitutive_Type::DIFF_LIN_ISO } );
//    tIWGUserDefinedInfo( 2 ).resize( 1 );
//    tIWGUserDefinedInfo( 2 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_BULK,
//                                                                { MSI::Dof_Type::TEMP },
//                                                                {{ MSI::Dof_Type::TEMP }},
//                                                                {fem::Property_Type::TEMP_LOAD },
//                                                                { fem::Constitutive_Type::DIFF_LIN_ISO } );
//    tIWGUserDefinedInfo( 3 ).resize( 1 );
//    tIWGUserDefinedInfo( 3 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_BULK,
//                                                                { MSI::Dof_Type::TEMP },
//                                                                {{ MSI::Dof_Type::TEMP }},
//                                                                {fem::Property_Type::TEMP_LOAD },
//                                                                { fem::Constitutive_Type::DIFF_LIN_ISO } );
//    tIWGUserDefinedInfo( 4 ).resize( 1 );
//    tIWGUserDefinedInfo( 4 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_DIRICHLET,
//                                                                { MSI::Dof_Type::TEMP },
//                                                                {{ MSI::Dof_Type::TEMP }},
//                                                                { fem::Property_Type::TEMP_DIRICHLET },
//                                                                { fem::Constitutive_Type::DIFF_LIN_ISO } );
//    tIWGUserDefinedInfo( 5 ).resize( 1 );
//    tIWGUserDefinedInfo( 5 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_NEUMANN,
//                                                                { MSI::Dof_Type::TEMP },
//                                                                {{ MSI::Dof_Type::TEMP }},
//                                                                { fem::Property_Type::TEMP_NEUMANN },
//                                                                moris::Cell< fem::Constitutive_Type >( 0 ) );
//
//    tIWGUserDefinedInfo( 6 ).resize( 1 );
//    tIWGUserDefinedInfo( 6 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_INTERFACE,
//                                                                { MSI::Dof_Type::TEMP },
//                                                                {{ MSI::Dof_Type::TEMP }},
//                                                                Cell< fem::Property_Type >( 0 ),
//                                                                {fem::Constitutive_Type::DIFF_LIN_ISO },
//                                                                {{ MSI::Dof_Type::TEMP }},
//                                                                Cell< fem::Property_Type >( 0 ),
//                                                                {fem::Constitutive_Type::DIFF_LIN_ISO } );
//    tIWGUserDefinedInfo( 7 ).resize( 1 );
//    tIWGUserDefinedInfo( 7 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_DIRICHLET,
//                                                                { MSI::Dof_Type::TEMP },
//                                                                {{ MSI::Dof_Type::TEMP }},
//                                                                { fem::Property_Type::TEMP_DIRICHLET },
//                                                                { fem::Constitutive_Type::DIFF_LIN_ISO } );
//
//    tIWGUserDefinedInfo( 8 ).resize( 1 );
//    tIWGUserDefinedInfo( 8 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_DIRICHLET,
//                                                                { MSI::Dof_Type::TEMP },
//                                                                {{ MSI::Dof_Type::TEMP }},
//                                                                { fem::Property_Type::TEMP_DIRICHLET },
//                                                                { fem::Constitutive_Type::DIFF_LIN_ISO } );
//
//    tIWGUserDefinedInfo( 9 ).resize( 1 );
//    tIWGUserDefinedInfo( 9 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_DIRICHLET,
//                                                                { MSI::Dof_Type::TEMP },
//                                                                {{ MSI::Dof_Type::TEMP }},
//                                                                { fem::Property_Type::TEMP_DIRICHLET },
//                                                                { fem::Constitutive_Type::DIFF_LIN_ISO } );
//
//    // create the property user defined infos
//    fem::Property_User_Defined_Info tConductivity( fem::Property_Type::CONDUCTIVITY,
//                                                   Cell< Cell< MSI::Dof_Type > >( 0 ),
//                                                   {{{ 1.0 }}},
//                                                   tConstValFunction2MatMDL,
//                                                   Cell< fem::PropertyFunc >( 0 ) );
//
//    fem::Property_User_Defined_Info tConductivity2( fem::Property_Type::CONDUCTIVITY,
//                                                    Cell< Cell< MSI::Dof_Type > >( 0 ),
//                                                    {{{ 100.0 }}},
//                                                    tConstValFunction2MatMDL,
//                                                    Cell< fem::PropertyFunc >( 0 ) );
//
//    fem::Property_User_Defined_Info tTempDirichlet( fem::Property_Type::TEMP_DIRICHLET,
//                                                    Cell< Cell< MSI::Dof_Type > >( 0 ),
//                                                    {{{ 100.0 }}},
//                                                    tConstValFunction2MatMDL,
//                                                    Cell< fem::PropertyFunc >( 0 ) );
//    fem::Property_User_Defined_Info tTempNeumann( fem::Property_Type::TEMP_NEUMANN,
//                                                  Cell< Cell< MSI::Dof_Type > >( 0 ),
//                                                  {{{ 20.0 }}},
//                                                  tConstValFunction2MatMDL,
//                                                  Cell< fem::PropertyFunc >( 0 ) );
//
//    fem::Property_User_Defined_Info tTempLoad1( fem::Property_Type::TEMP_LOAD,
//                                                Cell< Cell< MSI::Dof_Type > >( 0 ),
//                                                {{{ 0.0 }}},
//                                                tConstValFunction2MatMDL,
//                                                Cell< fem::PropertyFunc >( 0 ) );
//
//    fem::Property_User_Defined_Info tTempLoad2( fem::Property_Type::TEMP_LOAD,
//                                                Cell< Cell< MSI::Dof_Type > >( 0 ),
//                                                {{{ 0.0 }}},
//                                                tConstValFunction2MatMDL,
//                                                Cell< fem::PropertyFunc >( 0 ) );
//
//    // create property user defined info
//    Cell< Cell< Cell< fem::Property_User_Defined_Info > > > tPropertyUserDefinedInfo( 10 );
//    tPropertyUserDefinedInfo( 0 ).resize( 1 );
//    tPropertyUserDefinedInfo( 0 )( 0 ).resize( 2 );
//    tPropertyUserDefinedInfo( 0 )( 0 )( 0 ) = tConductivity;
//    tPropertyUserDefinedInfo( 0 )( 0 )( 1 ) = tTempLoad1;
//
//    tPropertyUserDefinedInfo( 1 ).resize( 1 );
//    tPropertyUserDefinedInfo( 1 )( 0 ).resize( 2 );
//    tPropertyUserDefinedInfo( 1 )( 0 )( 0 ) = tConductivity;
//    tPropertyUserDefinedInfo( 1 )( 0 )( 1 ) = tTempLoad1;
//
//    tPropertyUserDefinedInfo( 2 ).resize( 1 );
//    tPropertyUserDefinedInfo( 2 )( 0 ).resize( 2 );
//    tPropertyUserDefinedInfo( 2 )( 0 )( 0 ) = tConductivity2;
//    tPropertyUserDefinedInfo( 2 )( 0 )( 1 ) = tTempLoad2;
//
//    tPropertyUserDefinedInfo( 3 ).resize( 1 );
//    tPropertyUserDefinedInfo( 3 )( 0 ).resize( 2 );
//    tPropertyUserDefinedInfo( 3 )( 0 )( 0 ) = tConductivity2;
//    tPropertyUserDefinedInfo( 3 )( 0 )( 1 ) = tTempLoad2;
//
//    tPropertyUserDefinedInfo( 4 ).resize( 1 );
//    tPropertyUserDefinedInfo( 4 )( 0 ).resize( 2 );
//    tPropertyUserDefinedInfo( 4 )( 0 )( 0 ) = tConductivity2;
//    tPropertyUserDefinedInfo( 4 )( 0 )( 1 ) = tTempDirichlet;
//
//    tPropertyUserDefinedInfo( 5 ).resize( 1 );
//    tPropertyUserDefinedInfo( 5 )( 0 ).resize( 1 );
//    tPropertyUserDefinedInfo( 5 )( 0 )( 0 ) = tTempNeumann;
//
//    tPropertyUserDefinedInfo( 6 ).resize( 2 );
//    tPropertyUserDefinedInfo( 6 )( 0 ).resize( 1 );
//    tPropertyUserDefinedInfo( 6 )( 0 )( 0 ) = tConductivity;
//    tPropertyUserDefinedInfo( 6 )( 1 ).resize( 1 );
//    tPropertyUserDefinedInfo( 6 )( 1 )( 0 ) = tConductivity2;
//
//    tPropertyUserDefinedInfo( 7 ).resize( 1 );
//    tPropertyUserDefinedInfo( 7 )( 0 ).resize( 2 );
//    tPropertyUserDefinedInfo( 7 )( 0 )( 0 ) = tConductivity2;
//    tPropertyUserDefinedInfo( 7 )( 0 )( 1 ) = tTempDirichlet;
//
//    tPropertyUserDefinedInfo( 8 ).resize( 1 );
//    tPropertyUserDefinedInfo( 8 )( 0 ).resize( 2 );
//    tPropertyUserDefinedInfo( 8 )( 0 )( 0 ) = tConductivity;
//    tPropertyUserDefinedInfo( 8 )( 0 )( 1 ) = tTempDirichlet;
//
//    tPropertyUserDefinedInfo( 9 ).resize( 1 );
//    tPropertyUserDefinedInfo( 9 )( 0 ).resize( 2 );
//    tPropertyUserDefinedInfo( 9 )( 0 )( 0 ) = tConductivity;
//    tPropertyUserDefinedInfo( 9 )( 0 )( 1 ) = tTempDirichlet;
//
//
//    // create constitutive user defined info
//    fem::Constitutive_User_Defined_Info tDiffLinIso( fem::Constitutive_Type::DIFF_LIN_ISO,
//            {{ MSI::Dof_Type::TEMP }},
//            { fem::Property_Type::CONDUCTIVITY } );
//    // create constitutive user defined info
//    Cell< Cell< Cell< fem::Constitutive_User_Defined_Info > > > tConstitutiveUserDefinedInfo(10);
//    tConstitutiveUserDefinedInfo( 0 ).resize( 1 );
//    tConstitutiveUserDefinedInfo( 0 )( 0 ).resize( 1 );
//    tConstitutiveUserDefinedInfo( 0 )( 0 )( 0 ) = tDiffLinIso;
//
//    tConstitutiveUserDefinedInfo( 1 ).resize( 1 );
//    tConstitutiveUserDefinedInfo( 1 )( 0 ).resize( 1 );
//    tConstitutiveUserDefinedInfo( 1 )( 0 )( 0 ) = tDiffLinIso;
//
//    tConstitutiveUserDefinedInfo( 2 ).resize( 1 );
//    tConstitutiveUserDefinedInfo( 2 )( 0 ).resize( 1 );
//    tConstitutiveUserDefinedInfo( 2 )( 0 )( 0 ) = tDiffLinIso;
//
//    tConstitutiveUserDefinedInfo( 3 ).resize( 1 );
//    tConstitutiveUserDefinedInfo( 3 )( 0 ).resize( 1 );
//    tConstitutiveUserDefinedInfo( 3 )( 0 )( 0 ) = tDiffLinIso;
//
//
//    tConstitutiveUserDefinedInfo( 4 ).resize( 1 );
//    tConstitutiveUserDefinedInfo( 4 )( 0 ).resize( 1 );
//    tConstitutiveUserDefinedInfo( 4 )( 0 )( 0 ) = tDiffLinIso;
//
//    tConstitutiveUserDefinedInfo( 5 ).resize( 1 );
//
//    tConstitutiveUserDefinedInfo( 6 ).resize( 2 );
//    tConstitutiveUserDefinedInfo( 6 )( 0 ).resize( 1 );
//    tConstitutiveUserDefinedInfo( 6 )( 0 )( 0 ) = tDiffLinIso;
//    tConstitutiveUserDefinedInfo( 6 )( 1 ).resize( 1 );
//    tConstitutiveUserDefinedInfo( 6 )( 1 )( 0 ) = tDiffLinIso;
//
//
//    tConstitutiveUserDefinedInfo( 7 ).resize( 1 );
//    tConstitutiveUserDefinedInfo( 7 )( 0 ).resize( 1 );
//    tConstitutiveUserDefinedInfo( 7 )( 0 )( 0 ) = tDiffLinIso;
//
//    tConstitutiveUserDefinedInfo( 8 ).resize( 1 );
//    tConstitutiveUserDefinedInfo( 8 )( 0 ).resize( 1 );
//    tConstitutiveUserDefinedInfo( 8 )( 0 )( 0 ) = tDiffLinIso;
//
//    tConstitutiveUserDefinedInfo( 9 ).resize( 1 );
//    tConstitutiveUserDefinedInfo( 9 )( 0 ).resize( 1 );
//    tConstitutiveUserDefinedInfo( 9 )( 0 )( 0 ) = tDiffLinIso;
//
//
//
//    // create a list of active block-sets
//    std::string tInterfaceSideSetName = tEnrIntegMesh.get_interface_side_set_name(0,0,1);
//    std::string tDblInterfaceSideSetName = tEnrIntegMesh.get_dbl_interface_side_set_name(0,1);
//    moris::Cell< moris_index >  tSetList = {  tEnrIntegMesh.get_block_set_index("HMR_dummy_c_p0"),
//            tEnrIntegMesh.get_block_set_index("HMR_dummy_n_p0"),
//            tEnrIntegMesh.get_block_set_index("HMR_dummy_c_p1"),
//            tEnrIntegMesh.get_block_set_index("HMR_dummy_n_p1"),
//            tEnrIntegMesh.get_side_set_index("SideSet_2_n_p1"),
//            tEnrIntegMesh.get_side_set_index("SideSet_4_n_p1"),
//            tEnrIntegMesh.get_double_sided_set_index(tDblInterfaceSideSetName),
//            tEnrIntegMesh.get_side_set_index("SideSet_2_c_p1"),
//            tEnrIntegMesh.get_side_set_index("SideSet_2_c_p0"),
//            tEnrIntegMesh.get_side_set_index("SideSet_2_n_p0")};
//
//    moris::Cell< fem::Element_Type > tSetTypeList = { fem::Element_Type::BULK,
//            fem::Element_Type::BULK,
//            fem::Element_Type::BULK,
//            fem::Element_Type::BULK,
//            fem::Element_Type::SIDESET,
//            fem::Element_Type::SIDESET,
//            fem::Element_Type::DOUBLE_SIDESET,
//            fem::Element_Type::SIDESET,
//            fem::Element_Type::SIDESET,
//            fem::Element_Type::SIDESET};
//
//
//    // create model
//    mdl::Model * tModel = new mdl::Model( &tMeshManager,
//            0,
//            tSetList,
//            tSetTypeList,
//            tIWGUserDefinedInfo,
//            tPropertyUserDefinedInfo,
//            tConstitutiveUserDefinedInfo,
//            0,
//            false);
//
//    moris::Cell< enum MSI::Dof_Type > tDofTypes1( 1, MSI::Dof_Type::TEMP );
//
//    dla::Solver_Factory  tSolFactory;
//    std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( SolverType::AZTEC_IMPL );
//
//    tLinearSolverAlgorithm->set_param("AZ_diagnostics") = AZ_none;
//    tLinearSolverAlgorithm->set_param("AZ_output") = AZ_none;
//
//    dla::Linear_Solver tLinSolver;
//
//    tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );
//
//    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//    // STEP 2: create nonlinear solver and algorithm
//    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//
//    NLA::Nonlinear_Solver_Factory tNonlinFactory;
//    std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );
//
//    tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );
//
//    NLA::Nonlinear_Solver tNonlinearSolver;
//    tNonlinearSolver.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );
//
//    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//    // STEP 3: create time Solver and algorithm
//    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//    tsa::Time_Solver_Factory tTimeSolverFactory;
//    std::shared_ptr< tsa::Time_Solver_Algorithm > tTimeSolverAlgorithm = tTimeSolverFactory.create_time_solver( tsa::TimeSolverType::MONOLITHIC );
//
//    tTimeSolverAlgorithm->set_nonlinear_solver( &tNonlinearSolver );
//
//    tsa::Time_Solver tTimeSolver;
//
//    tTimeSolver.set_time_solver_algorithm( tTimeSolverAlgorithm );
//
//    NLA::SOL_Warehouse tSolverWarehouse;
//
//    tSolverWarehouse.set_solver_interface(tModel->get_solver_interface());
//
//    tNonlinearSolver.set_solver_warehouse( &tSolverWarehouse );
//    tTimeSolver.set_solver_warehouse( &tSolverWarehouse );
//
//    //               tNonlinearSolver.set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }} );
//    //               tTimeSolver.set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }} );
//    tNonlinearSolver.set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
//    tTimeSolver.set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
//
//    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//    // STEP 4: Solve and check
//    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//
//    tTimeSolver.solve();
//
//
//    Matrix<DDRMat> tFullSol;
//    tTimeSolver.get_full_solution(tFullSol);
//
//
//    Matrix<DDRMat> tIntegSol = tModel->get_solution_for_integration_mesh_output( MSI::Dof_Type::TEMP );
//
//    Matrix<DDRMat> tSTKIntegSol(tIntegMesh11->get_num_entities(EntityRank::NODE),1);
//
//    for(moris::uint i = 0; i < tIntegMesh11->get_num_entities(EntityRank::NODE); i++)
//    {
//        moris::moris_id tID = tIntegMesh11->get_glb_entity_id_from_entity_loc_index(i,EntityRank::NODE);
//        moris::moris_index tMyIndex = tEnrIntegMesh.get_loc_entity_ind_from_entity_glb_id(tID,EntityRank::NODE);
//
//        tSTKIntegSol(i) = tIntegSol(tMyIndex);
//    }
//
//
//    // add solution field to integration mesh
//    tIntegMesh11->add_mesh_field_real_scalar_data_loc_inds(tIntegSolFieldName1,EntityRank::NODE,tSTKIntegSol);
//
//
//    std::string tMeshOutputFile = "./mdl_exo/fiber_problem_with_XTK" + std::to_string(1) + "_b"+std::to_string(0)+".e";
//    tIntegMesh11->create_output_mesh(tMeshOutputFile);
//
//    delete tIntegMesh11;
//
//}
//
//TEST_CASE("fiber_problem_test_02","[GE],[fiber_test_displ]")
//{
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
//    uint tNumberOfFibers = 1;
//    std::cout<<"-------------------------------------"<<std::endl;
//    std::cout<<"number of fibers being analyzed:  "<<tNumberOfFibers<<std::endl;
//    std::cout<<"-------------------------------------"<<std::endl;
//    Cell< Matrix< DDRMat > > tFieldData( 1 );
//
//    std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );
//
//    moris::tic tTimer0;
//    for( uint k=0; k<1; ++k )
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
//    std::shared_ptr< hmr::Interpolation_Mesh_HMR > tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );
//    std::shared_ptr< moris::hmr::Integration_Mesh_HMR > tIntegrationMesh = tHMR.create_integration_mesh( 1, 0,*tInterpMesh );
//
//    mtk::Mesh_Manager tMesh1;
//
//    real tElapsedTime0 = tTimer0.toc<moris::chronos::milliseconds>().wall;
//    tElapsedTime0 /= 1000;
//    std::cout<<"==============================================="<< std::endl;
//    std::cout<<"Total time for evaluation: "<< tElapsedTime0 << std::endl << std::endl;
//    std::cout<<"==============================================="<< std::endl;
//
//    moris::ge::GEN_Geom_Data tGeomData( tFieldData( 0 ) );
//    moris::Cell<moris::ge::GEN_Geometry*> tGeometryVector = {&tGeomData};
//
//    size_t tModelDimension = 3;
//    moris::ge::GEN_Phase_Table      tPhaseTable( tGeometryVector.size(),  Phase_Table_Structure::EXP_BASE_2 );
//    moris::ge::GEN_Geometry_Engine  tGENGeometryEngine( tGeometryVector,tPhaseTable,tModelDimension );
//    xtk::Model                      tXTKModel( tModelDimension,tInterpMesh.get(),tGENGeometryEngine );
//    tXTKModel.mVerbose = false;
//
//    Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4};
//    tXTKModel.decompose(tDecompositionMethods);
//
//    // output solution and meshes
//
//    xtk::Output_Options tOutputOptions1;
//    tOutputOptions1.mAddNodeSets = false;
//    tOutputOptions1.mAddSideSets = false;
//    tOutputOptions1.mAddClusters = false;
//
//    std::string tIntegSolFieldNameUX = "UX";
//    std::string tIntegSolFieldNameUY = "UY";
//    std::string tIntegSolFieldNameUZ = "UZ";
//    tOutputOptions1.mRealNodeExternalFieldNames = { tIntegSolFieldNameUX, tIntegSolFieldNameUY, tIntegSolFieldNameUZ };
//
//    moris::mtk::Integration_Mesh* tIntegMesh11 = tXTKModel.get_output_mesh(tOutputOptions1);
//
////    std::string tMeshOutputFile1 = "./output_fibers_29Oct.e";
////    tIntegMesh11->create_output_mesh(tMeshOutputFile1);
//
//    tXTKModel.perform_basis_enrichment(EntityRank::BSPLINE_1,0);
//
//    xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
//    xtk::Enriched_Integration_Mesh   & tEnrIntegMesh = tXTKModel.get_enriched_integ_mesh();
//
//    //==============================
//
//    // place the pair in mesh manager
//    mtk::Mesh_Manager tMeshManager;
//    tMeshManager.register_mesh_pair( &tEnrInterpMesh, &tEnrIntegMesh);
//
//    // create IWG user defined info
//    Cell< Cell< fem::IWG_User_Defined_Info > > tIWGUserDefinedInfo( 10 );
//    tIWGUserDefinedInfo( 0 ).resize( 1 );
//    tIWGUserDefinedInfo( 0 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::STRUC_LINEAR_BULK,
//                                                                { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ },
//                                                                {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }},
//                                                                moris::Cell< fem::Property_Type >( 0 ),
//                                                                { fem::Constitutive_Type::STRUC_LIN_ISO } );
//    tIWGUserDefinedInfo( 1 ).resize( 1 );
//    tIWGUserDefinedInfo( 1 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::STRUC_LINEAR_BULK,
//                                                                { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ },
//                                                                {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }},
//                                                                moris::Cell< fem::Property_Type >( 0 ),
//                                                                { fem::Constitutive_Type::STRUC_LIN_ISO } );
//    tIWGUserDefinedInfo( 2 ).resize( 1 );
//    tIWGUserDefinedInfo( 2 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::STRUC_LINEAR_BULK,
//                                                                { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ },
//                                                                {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }},
//                                                                moris::Cell< fem::Property_Type >( 0 ),
//                                                                { fem::Constitutive_Type::STRUC_LIN_ISO } );
//    tIWGUserDefinedInfo( 3 ).resize( 1 );
//    tIWGUserDefinedInfo( 3 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_BULK,
//                                                                { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ },
//                                                                {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }},
//                                                                moris::Cell< fem::Property_Type >( 0 ),
//                                                                { fem::Constitutive_Type::STRUC_LIN_ISO } );
//    tIWGUserDefinedInfo( 4 ).resize( 1 );
//    tIWGUserDefinedInfo( 4 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::STRUC_LINEAR_DIRICHLET,
//                                                                { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ },
//                                                                {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }},
//                                                                { fem::Property_Type::STRUC_DIRICHLET },
//                                                                { fem::Constitutive_Type::STRUC_LIN_ISO } );
//    tIWGUserDefinedInfo( 5 ).resize( 1 );
//    tIWGUserDefinedInfo( 5 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::STRUC_LINEAR_NEUMANN,
//                                                                { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ },
//                                                                {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }},
//                                                                { fem::Property_Type::STRUC_NEUMANN },
//                                                                moris::Cell< fem::Constitutive_Type >( 0 ) );
//    tIWGUserDefinedInfo( 6 ).resize( 1 );
//    tIWGUserDefinedInfo( 6 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::STRUC_LINEAR_INTERFACE,
//                                                                { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ },
//                                                                {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }},
//                                                                Cell< fem::Property_Type >( 0 ),
//                                                                {fem::Constitutive_Type::STRUC_LIN_ISO },
//                                                                {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }},
//                                                                Cell< fem::Property_Type >( 0 ),
//                                                                {fem::Constitutive_Type::STRUC_LIN_ISO } );
//    tIWGUserDefinedInfo( 7 ).resize( 1 );
//    tIWGUserDefinedInfo( 7 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::STRUC_LINEAR_DIRICHLET,
//                                                                { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ },
//                                                                {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }},
//                                                                { fem::Property_Type::STRUC_DIRICHLET },
//                                                                { fem::Constitutive_Type::STRUC_LIN_ISO } );
//
//    tIWGUserDefinedInfo( 8 ).resize( 1 );
//    tIWGUserDefinedInfo( 8 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::STRUC_LINEAR_DIRICHLET,
//                                                                { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ },
//                                                                {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }},
//                                                                { fem::Property_Type::STRUC_DIRICHLET },
//                                                                { fem::Constitutive_Type::STRUC_LIN_ISO } );
//
//    tIWGUserDefinedInfo( 9 ).resize( 1 );
//    tIWGUserDefinedInfo( 9 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::STRUC_LINEAR_DIRICHLET,
//                                                                { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ },
//                                                                {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }},
//                                                                { fem::Property_Type::STRUC_DIRICHLET },
//                                                                { fem::Constitutive_Type::STRUC_LIN_ISO } );
//
//    // create the property user defined infos
//    fem::Property_User_Defined_Info tYoungs_Modulus( fem::Property_Type::YOUNGS_MODULUS,
//                                                     Cell< Cell< MSI::Dof_Type > >( 0 ),
//                                                     {{{ 1.0 }}},
//                                                     tConstValFunction,
//                                                     Cell< fem::PropertyFunc >( 0 ) );
//
//    fem::Property_User_Defined_Info tYoungs_Modulus2( fem::Property_Type::YOUNGS_MODULUS,
//                                                      Cell< Cell< MSI::Dof_Type > >( 0 ),
//                                                      {{{ 1.0 }}},
//                                                      tConstValFunction,
//                                                      Cell< fem::PropertyFunc >( 0 ) );
//
//    fem::Property_User_Defined_Info tPoissons_Ratio( fem::Property_Type::POISSONS_RATIO,
//                                                     Cell< Cell< MSI::Dof_Type > >( 0 ),
//                                                     {{{ 0.0 }}},
//                                                     tConstValFunction,
//                                                     Cell< fem::PropertyFunc >( 0 ) );
//
//    fem::Property_User_Defined_Info tStrucDirichlet( fem::Property_Type::STRUC_DIRICHLET,
//                                                     Cell< Cell< MSI::Dof_Type > >(0),
//                                                     {{{ 0.0 }, { 0.0 }, { 0.0 }}},
//                                                     tConstValFunction,
//                                                     Cell< fem::PropertyFunc >( 0 ));
//
//    fem::Property_User_Defined_Info tStrucNeumann( fem::Property_Type::STRUC_NEUMANN,
//                                                   Cell< Cell< MSI::Dof_Type > >( 0 ),
//                                                   {{{ 1.0 }, {0.0}, {0.0}}},
//                                                   tConstValFunction,
//                                                   Cell< fem::PropertyFunc >( 0 ) );
//
////    fem::Property_User_Defined_Info tTempLoad1( fem::Property_Type::TEMP_LOAD,
////            Cell< Cell< MSI::Dof_Type > >( 0 ),
////            {{{ 0.0 }}},
////            tConstValFunction2MatMDL,
////            Cell< fem::PropertyFunc >( 0 ) );
////
////    fem::Property_User_Defined_Info tTempLoad2( fem::Property_Type::TEMP_LOAD,
////            Cell< Cell< MSI::Dof_Type > >( 0 ),
////            {{{ 0.0 }}},
////            tConstValFunction2MatMDL,
////            Cell< fem::PropertyFunc >( 0 ) );
//
//    // create property user defined info
//    Cell< Cell< Cell< fem::Property_User_Defined_Info > > > tPropertyUserDefinedInfo( 10 );
//    tPropertyUserDefinedInfo( 0 ).resize( 1 );
//    tPropertyUserDefinedInfo( 0 )( 0 ).resize( 2 );
//    tPropertyUserDefinedInfo( 0 )( 0 )( 0 ) = tYoungs_Modulus;
//    tPropertyUserDefinedInfo( 0 )( 0 )( 1 ) = tPoissons_Ratio;
//
//    tPropertyUserDefinedInfo( 1 ).resize( 1 );
//    tPropertyUserDefinedInfo( 1 )( 0 ).resize( 2 );
//    tPropertyUserDefinedInfo( 1 )( 0 )( 0 ) = tYoungs_Modulus;
//    tPropertyUserDefinedInfo( 1 )( 0 )( 1 ) = tPoissons_Ratio;
//
//    tPropertyUserDefinedInfo( 2 ).resize( 1 );
//    tPropertyUserDefinedInfo( 2 )( 0 ).resize( 2 );
//    tPropertyUserDefinedInfo( 2 )( 0 )( 0 ) = tYoungs_Modulus2;
//    tPropertyUserDefinedInfo( 2 )( 0 )( 1 ) = tPoissons_Ratio;
//
//    tPropertyUserDefinedInfo( 3 ).resize( 1 );
//    tPropertyUserDefinedInfo( 3 )( 0 ).resize( 2 );
//    tPropertyUserDefinedInfo( 3 )( 0 )( 0 ) = tYoungs_Modulus2;
//    tPropertyUserDefinedInfo( 3 )( 0 )( 1 ) = tPoissons_Ratio;
//
//    tPropertyUserDefinedInfo( 4 ).resize( 1 );
//    tPropertyUserDefinedInfo( 4 )( 0 ).resize( 3 );
//    tPropertyUserDefinedInfo( 4 )( 0 )( 0 ) = tYoungs_Modulus2;
//    tPropertyUserDefinedInfo( 4 )( 0 )( 1 ) = tPoissons_Ratio;
//    tPropertyUserDefinedInfo( 4 )( 0 )( 2 ) = tStrucDirichlet;
//
//
//    tPropertyUserDefinedInfo( 5 ).resize( 1 );
//    tPropertyUserDefinedInfo( 5 )( 0 ).resize( 1 );
//    tPropertyUserDefinedInfo( 5 )( 0 )( 0 ) = tStrucNeumann;
//
//    tPropertyUserDefinedInfo( 6 ).resize( 2 );
//    tPropertyUserDefinedInfo( 6 )( 0 ).resize( 2 );
//    tPropertyUserDefinedInfo( 6 )( 0 )( 0 ) = tYoungs_Modulus;
//    tPropertyUserDefinedInfo( 6 )( 0 )( 1 ) = tPoissons_Ratio;
//    tPropertyUserDefinedInfo( 6 )( 1 ).resize( 2 );
//    tPropertyUserDefinedInfo( 6 )( 1 )( 0 ) = tYoungs_Modulus2;
//    tPropertyUserDefinedInfo( 6 )( 1 )( 1 ) = tPoissons_Ratio;
//
//    tPropertyUserDefinedInfo( 7 ).resize( 1 );
//    tPropertyUserDefinedInfo( 7 )( 0 ).resize( 3 );
//    tPropertyUserDefinedInfo( 7 )( 0 )( 0 ) = tYoungs_Modulus2;
//    tPropertyUserDefinedInfo( 7 )( 0 )( 1 ) = tPoissons_Ratio;
//    tPropertyUserDefinedInfo( 7 )( 0 )( 2 ) = tStrucDirichlet;
//
//    tPropertyUserDefinedInfo( 8 ).resize( 1 );
//    tPropertyUserDefinedInfo( 8 )( 0 ).resize( 3 );
//    tPropertyUserDefinedInfo( 8 )( 0 )( 0 ) = tYoungs_Modulus;
//    tPropertyUserDefinedInfo( 8 )( 0 )( 1 ) = tPoissons_Ratio;
//    tPropertyUserDefinedInfo( 8 )( 0 )( 2 ) = tStrucDirichlet;
//
//    tPropertyUserDefinedInfo( 9 ).resize( 1 );
//    tPropertyUserDefinedInfo( 9 )( 0 ).resize( 3 );
//    tPropertyUserDefinedInfo( 9 )( 0 )( 0 ) = tYoungs_Modulus;
//    tPropertyUserDefinedInfo( 9 )( 0 )( 1 ) = tPoissons_Ratio;
//    tPropertyUserDefinedInfo( 9 )( 0 )( 2 ) = tStrucDirichlet;
//
//    // create constitutive user defined info
//    fem::Constitutive_User_Defined_Info tDiffLinIso( fem::Constitutive_Type::STRUC_LIN_ISO,
//                                                     {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }},
//                                                     { fem::Property_Type::YOUNGS_MODULUS, fem::Property_Type::POISSONS_RATIO } );
//    // create constitutive user defined info
//    Cell< Cell< Cell< fem::Constitutive_User_Defined_Info > > > tConstitutiveUserDefinedInfo(10);
//    tConstitutiveUserDefinedInfo( 0 ).resize( 1 );
//    tConstitutiveUserDefinedInfo( 0 )( 0 ).resize( 1 );
//    tConstitutiveUserDefinedInfo( 0 )( 0 )( 0 ) = tDiffLinIso;
//
//    tConstitutiveUserDefinedInfo( 1 ).resize( 1 );
//    tConstitutiveUserDefinedInfo( 1 )( 0 ).resize( 1 );
//    tConstitutiveUserDefinedInfo( 1 )( 0 )( 0 ) = tDiffLinIso;
//
//    tConstitutiveUserDefinedInfo( 2 ).resize( 1 );
//    tConstitutiveUserDefinedInfo( 2 )( 0 ).resize( 1 );
//    tConstitutiveUserDefinedInfo( 2 )( 0 )( 0 ) = tDiffLinIso;
//
//    tConstitutiveUserDefinedInfo( 3 ).resize( 1 );
//    tConstitutiveUserDefinedInfo( 3 )( 0 ).resize( 1 );
//    tConstitutiveUserDefinedInfo( 3 )( 0 )( 0 ) = tDiffLinIso;
//
//
//    tConstitutiveUserDefinedInfo( 4 ).resize( 1 );
//    tConstitutiveUserDefinedInfo( 4 )( 0 ).resize( 1 );
//    tConstitutiveUserDefinedInfo( 4 )( 0 )( 0 ) = tDiffLinIso;
//
//    tConstitutiveUserDefinedInfo( 5 ).resize( 1 );
//
//    tConstitutiveUserDefinedInfo( 6 ).resize( 2 );
//    tConstitutiveUserDefinedInfo( 6 )( 0 ).resize( 1 );
//    tConstitutiveUserDefinedInfo( 6 )( 0 )( 0 ) = tDiffLinIso;
//    tConstitutiveUserDefinedInfo( 6 )( 1 ).resize( 1 );
//    tConstitutiveUserDefinedInfo( 6 )( 1 )( 0 ) = tDiffLinIso;
//
//
//    tConstitutiveUserDefinedInfo( 7 ).resize( 1 );
//    tConstitutiveUserDefinedInfo( 7 )( 0 ).resize( 1 );
//    tConstitutiveUserDefinedInfo( 7 )( 0 )( 0 ) = tDiffLinIso;
//
//    tConstitutiveUserDefinedInfo( 8 ).resize( 1 );
//    tConstitutiveUserDefinedInfo( 8 )( 0 ).resize( 1 );
//    tConstitutiveUserDefinedInfo( 8 )( 0 )( 0 ) = tDiffLinIso;
//
//    tConstitutiveUserDefinedInfo( 9 ).resize( 1 );
//    tConstitutiveUserDefinedInfo( 9 )( 0 ).resize( 1 );
//    tConstitutiveUserDefinedInfo( 9 )( 0 )( 0 ) = tDiffLinIso;
//
//    // create a list of active block-sets
//    std::string tInterfaceSideSetName = tEnrIntegMesh.get_interface_side_set_name(0,0,1);
//    std::string tDblInterfaceSideSetName = tEnrIntegMesh.get_dbl_interface_side_set_name(0,1);
//    moris::Cell< moris_index >  tSetList = {  tEnrIntegMesh.get_block_set_index("HMR_dummy_c_p0"),
//                                              tEnrIntegMesh.get_block_set_index("HMR_dummy_n_p0"),
//                                              tEnrIntegMesh.get_block_set_index("HMR_dummy_c_p1"),
//                                              tEnrIntegMesh.get_block_set_index("HMR_dummy_n_p1"),
//                                              tEnrIntegMesh.get_side_set_index("SideSet_2_n_p1"),
//                                              tEnrIntegMesh.get_side_set_index("SideSet_4_n_p1"),
//                                              tEnrIntegMesh.get_double_sided_set_index(tDblInterfaceSideSetName),
//                                              tEnrIntegMesh.get_side_set_index("SideSet_2_c_p1"),
//                                              tEnrIntegMesh.get_side_set_index("SideSet_2_c_p0"),
//                                              tEnrIntegMesh.get_side_set_index("SideSet_2_n_p0")};
//
//    moris::Cell< fem::Element_Type > tSetTypeList = { fem::Element_Type::BULK,
//                                                      fem::Element_Type::BULK,
//                                                      fem::Element_Type::BULK,
//                                                      fem::Element_Type::BULK,
//                                                      fem::Element_Type::SIDESET,
//                                                      fem::Element_Type::SIDESET,
//                                                      fem::Element_Type::DOUBLE_SIDESET,
//                                                      fem::Element_Type::SIDESET,
//                                                      fem::Element_Type::SIDESET,
//                                                      fem::Element_Type::SIDESET};
//
//
//    // create model
//    mdl::Model * tModel = new mdl::Model( &tMeshManager,
//                                          0,
//                                          tSetList,
//                                          tSetTypeList,
//                                          tIWGUserDefinedInfo,
//                                          tPropertyUserDefinedInfo,
//                                          tConstitutiveUserDefinedInfo,
//                                          0,
//                                          false);
//
//    moris::Cell< enum MSI::Dof_Type > tDofTypes1( 3 );
//    tDofTypes1( 0 ) = MSI::Dof_Type::UX;
//    tDofTypes1( 1 ) = MSI::Dof_Type::UY;
//    tDofTypes1( 2 ) = MSI::Dof_Type::UZ;
//
//    dla::Solver_Factory  tSolFactory;
//    std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( SolverType::AZTEC_IMPL );
//
//    tLinearSolverAlgorithm->set_param("AZ_diagnostics") = AZ_none;
//    tLinearSolverAlgorithm->set_param("AZ_output") = AZ_none;
//    tLinearSolverAlgorithm->set_param("AZ_ilut_fill") = 10.0;         // fill for preconditioner
//    tLinearSolverAlgorithm->set_param("AZ_orthog") = AZ_modified;   // only better if in serial
//    tLinearSolverAlgorithm->set_param("AZ_solver") = AZ_gmres;
//
//    dla::Linear_Solver tLinSolver;
//
//    tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );
//
//    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//    // STEP 2: create nonlinear solver and algorithm
//    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//
//    NLA::Nonlinear_Solver_Factory tNonlinFactory;
//    std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );
//
//    tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );
//
//    NLA::Nonlinear_Solver tNonlinearSolver;
//    tNonlinearSolver.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );
//
//    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//    // STEP 3: create time Solver and algorithm
//    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//    tsa::Time_Solver_Factory tTimeSolverFactory;
//    std::shared_ptr< tsa::Time_Solver_Algorithm > tTimeSolverAlgorithm = tTimeSolverFactory.create_time_solver( tsa::TimeSolverType::MONOLITHIC );
//
//    tTimeSolverAlgorithm->set_nonlinear_solver( &tNonlinearSolver );
//
//    tsa::Time_Solver tTimeSolver;
//
//    tTimeSolver.set_time_solver_algorithm( tTimeSolverAlgorithm );
//
//    NLA::SOL_Warehouse tSolverWarehouse;
//
//    tSolverWarehouse.set_solver_interface(tModel->get_solver_interface());
//
//    tNonlinearSolver.set_solver_warehouse( &tSolverWarehouse );
//    tTimeSolver.set_solver_warehouse( &tSolverWarehouse );
//
//    tNonlinearSolver.set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }} );
//    tTimeSolver.set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }} );
//
//    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//    // STEP 4: Solve and check
//    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//
//    tTimeSolver.solve();
//
//    Matrix<DDRMat> tFullSol;
//    tTimeSolver.get_full_solution(tFullSol);
//
//    Matrix<DDRMat> tIntegSolUX = tModel->get_solution_for_integration_mesh_output( MSI::Dof_Type::UX );
//    Matrix<DDRMat> tIntegSolUY = tModel->get_solution_for_integration_mesh_output( MSI::Dof_Type::UY );
//    Matrix<DDRMat> tIntegSolUZ = tModel->get_solution_for_integration_mesh_output( MSI::Dof_Type::UZ );
//
//    Matrix<DDRMat> tSTKIntegSolUX(tIntegMesh11->get_num_entities(EntityRank::NODE),1);
//    Matrix<DDRMat> tSTKIntegSolUY(tIntegMesh11->get_num_entities(EntityRank::NODE),1);
//    Matrix<DDRMat> tSTKIntegSolUZ(tIntegMesh11->get_num_entities(EntityRank::NODE),1);
//
//    for(moris::uint i = 0; i < tIntegMesh11->get_num_entities(EntityRank::NODE); i++)
//    {
//        moris::moris_id tID = tIntegMesh11->get_glb_entity_id_from_entity_loc_index(i,EntityRank::NODE);
//        moris::moris_index tMyIndex = tEnrIntegMesh.get_loc_entity_ind_from_entity_glb_id(tID,EntityRank::NODE);
//
//        tSTKIntegSolUX(i) = tIntegSolUX(tMyIndex);
//        tSTKIntegSolUY(i) = tIntegSolUY(tMyIndex);
//        tSTKIntegSolUZ(i) = tIntegSolUZ(tMyIndex);
//    }
//
//
//    // add solution field to integration mesh
//    tIntegMesh11->add_mesh_field_real_scalar_data_loc_inds(tIntegSolFieldNameUX,EntityRank::NODE,tSTKIntegSolUX);
//    tIntegMesh11->add_mesh_field_real_scalar_data_loc_inds(tIntegSolFieldNameUY,EntityRank::NODE,tSTKIntegSolUY);
//    tIntegMesh11->add_mesh_field_real_scalar_data_loc_inds(tIntegSolFieldNameUZ,EntityRank::NODE,tSTKIntegSolUZ);
//
//    std::string tMeshOutputFile = "./mdl_exo/fiber_problem_displ_with_XTK" + std::to_string(1) + "_b"+std::to_string(0)+".e";
//    tIntegMesh11->create_output_mesh(tMeshOutputFile);
//
//    delete tIntegMesh11;
//
//}
//
//TEST_CASE("fiber_problem_test_03","[GE],[fiber_test_dummy]")
//{
//    uint tLagrangeMeshIndex = 0;
//    //  HMR Parameters setup
//    hmr::ParameterList tParameters = hmr::create_hmr_parameter_list();
//
//    tParameters.set( "number_of_elements_per_dimension", "80, 40, 16" );
//    tParameters.set( "domain_dimensions",                "40, 20, 8" );
//    tParameters.set( "domain_offset",                    "-0, -0, -0" );
//
//    tParameters.set( "domain_sidesets",                  "1, 2, 3, 4, 5, 6" );
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
//    uint tNumberOfFibers = 1;
//    std::cout<<"-------------------------------------"<<std::endl;
//    std::cout<<"number of fibers being analyzed:  "<<tNumberOfFibers<<std::endl;
//    std::cout<<"-------------------------------------"<<std::endl;
//    Cell< Matrix< DDRMat > > tFieldData( 1 );
//
//    std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );
//    //------------------------------------------------------------------------------
//
//    moris::tic tTimer0;
//    for( uint k=0; k<1; ++k )
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
//    std::shared_ptr< hmr::Interpolation_Mesh_HMR > tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );
//    std::shared_ptr< moris::hmr::Integration_Mesh_HMR > tIntegrationMesh = tHMR.create_integration_mesh( 1, 0,*tInterpMesh );
//
//    mtk::Mesh_Manager tMesh1;
//
//    real tElapsedTime0 = tTimer0.toc<moris::chronos::milliseconds>().wall;
//    tElapsedTime0 /= 1000;
//    std::cout<<"==============================================="<< std::endl;
//    std::cout<<"Total time for evaluation: "<< tElapsedTime0 << std::endl << std::endl;
//    std::cout<<"==============================================="<< std::endl;
//
//    moris::ge::GEN_Geom_Data tGeomData( tFieldData( 0 ) );
//    moris::Cell<moris::ge::GEN_Geometry*> tGeometryVector = {&tGeomData};
//
//    size_t tModelDimension = 3;
//    moris::ge::GEN_Phase_Table      tPhaseTable( tGeometryVector.size(),  Phase_Table_Structure::EXP_BASE_2 );
//    moris::ge::GEN_Geometry_Engine  tGENGeometryEngine( tGeometryVector,tPhaseTable,tModelDimension );
//    xtk::Model                      tXTKModel( tModelDimension,tInterpMesh.get(),tGENGeometryEngine );
//    tXTKModel.mVerbose = false;
//
//    Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4};
//    tXTKModel.decompose(tDecompositionMethods);
//
//    // output solution and meshes
//
//    xtk::Output_Options tOutputOptions1;
//    tOutputOptions1.mAddNodeSets = false;
//    tOutputOptions1.mAddSideSets = false;
//    tOutputOptions1.mAddClusters = false;
//
//    std::string tIntegSolFieldNameUX = "UX";
//    std::string tIntegSolFieldNameUY = "UY";
//    std::string tIntegSolFieldNameUZ = "UZ";
//    tOutputOptions1.mRealNodeExternalFieldNames = { tIntegSolFieldNameUX, tIntegSolFieldNameUY, tIntegSolFieldNameUZ };
//
//    moris::mtk::Integration_Mesh* tIntegMesh11 = tXTKModel.get_output_mesh(tOutputOptions1);
//
//    std::string tMeshOutputFile1 = "./output_fibers_dummy.e";
//    tIntegMesh11->create_output_mesh(tMeshOutputFile1);
//
//
//
//
//    delete tIntegMesh11;
//}
//
//}   // ge namespace
