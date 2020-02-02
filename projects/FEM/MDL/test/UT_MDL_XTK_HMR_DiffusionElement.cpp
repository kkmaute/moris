/*
 * UT_MDL_XTK_HMR_DiffusionElement.cpp
 *
 *  Created on: Jun 27, 2019
 *      Author: doble
 */

#include "catch.hpp"

#include "cl_XTK_Model.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
#include "cl_Geom_Field.hpp"
#include "typedefs.hpp"

#include "cl_MTK_Mesh_Manager.hpp"

#include "cl_MTK_Vertex.hpp"    //MTK
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh.hpp"

#include "cl_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_Tools.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Scalar_Field_Info.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Mesh_Data_STK.hpp"
#include "cl_MTK_Mesh_Core_STK.hpp"
#include "cl_MTK_Interpolation_Mesh_STK.hpp"
#include "cl_MTK_Integration_Mesh_STK.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Double_Side_Cluster.hpp"
#include "cl_MTK_Double_Side_Cluster_Input.hpp"
#include "cl_MTK_Side_Cluster.hpp"
#include "cl_MTK_Side_Cluster_Input.hpp"
#include "cl_MTK_Writer_Exodus.hpp"

#include "cl_Matrix.hpp"        //LINALG
#include "linalg_typedefs.hpp"
#include "fn_equal_to.hpp" // ALG/src

#include "cl_FEM_NodeProxy.hpp"                //FEM/INT/src
#include "cl_FEM_ElementProxy.hpp"             //FEM/INT/src
#include "cl_FEM_Node_Base.hpp"                //FEM/INT/src
#include "cl_FEM_Element_Factory.hpp"          //FEM/INT/src
#include "cl_FEM_IWG_Factory.hpp"              //FEM/INT/src
#include "cl_FEM_IQI_Factory.hpp"              //FEM/INT/src
#include "cl_FEM_CM_Factory.hpp"              //FEM/INT/src
#include "cl_FEM_SP_Factory.hpp"              //FEM/INT/src
#include "cl_FEM_Set_User_Info.hpp"              //FEM/INT/src

#include "cl_MDL_Model.hpp"

#include "cl_VIS_Factory.hpp"
#include "cl_VIS_Visualization_Mesh.hpp"
#include "cl_VIS_Output_Manager.hpp"

#include "cl_HMR_Mesh_Interpolation.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Background_Mesh.hpp" //HMR/src
#include "cl_HMR_BSpline_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Element.hpp" //HMR/src
#include "cl_HMR_Factory.hpp" //HMR/src
#include "cl_HMR_Field.hpp"
#include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Parameters.hpp" //HMR/src

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

#include "fn_norm.hpp"

#include "../projects/GEN/src/geometry/cl_GEN_Geom_Field.hpp"
#include "../projects/GEN/src/geometry/cl_GEN_Geometry.hpp"

namespace moris
{
moris::real
LevelSetPlaneFunction( const moris::Matrix< moris::DDRMat > & aPoint )
{

    real mXn = 0;
    real mYn = 0;
    real mZn = 1.0;
    real mXc = 1.011;
    real mYc = 1.011;
    real mZc = 1.411;
    return mXn*(aPoint(0)-mXc) + mYn*(aPoint(1)-mYc) + mZn*(aPoint(2)-mZc);
}

moris::real
LevelSetSphereFunction( const moris::Matrix< moris::DDRMat > & aPoint )
{

    real mR  = 0.77;
    real mXc = 0.0;
    real mYc = 0.0;
    real mZc = 0.0;
    return   (aPoint( 0) - mXc) * (aPoint( 0) - mXc)
            + (aPoint( 1) - mYc) * (aPoint( 1) - mYc)
            + (aPoint( 2) - mZc) * (aPoint( 2) - mZc)
            - (mR * mR);

}

moris::real
LevelSetSphereCylinder(const moris::Matrix< moris::DDRMat > & aPoint )
{
    moris::Matrix<moris::DDRMat> aCenter = {{0.0},{0.0},{0.0}};
    moris::Matrix<moris::DDRMat> aAxis   = {{0.0},{1.0},{0.0}};
    moris::real aRad = 0.77;
    moris::real aLength = 5;

    MORIS_ASSERT(aCenter.numel() == 3,"Centers need to have length 3");
    MORIS_ASSERT(aAxis.numel() == 3, "axis need to have length 3");

    Cell<moris::real> relativePosition = {(aPoint(0) - aCenter(0)),(aPoint(1) - aCenter(1)),(aPoint(2) - aCenter(2))};
    moris::real lsFromLeft = (relativePosition(0)*(-aAxis(0)) + relativePosition(1)*(-aAxis(1))+ relativePosition(2)*(-aAxis(2))) - aLength/2.0;
    moris::real lsFromRight = (relativePosition(0)*(aAxis(0)) + relativePosition(1)*(aAxis(1))+ relativePosition(2)*(aAxis(2))) - aLength/2.0;

    moris::real axialCrd = (relativePosition(0)*(aAxis(0)) + relativePosition(1)*(aAxis(1))+ relativePosition(2)*(aAxis(2)));
    Cell<moris::real> radDir = {(relativePosition(0) - aAxis(0)*axialCrd), (relativePosition(1) - aAxis(1)*axialCrd),(relativePosition(2) - aAxis(2)*axialCrd)};
    moris::real radDist = std::pow(radDir(0)*radDir(0)+radDir(1)*radDir(1)+radDir(2)*radDir(2), 0.5);
    moris::real lsFromRad = radDist - aRad;

    return -std::max(std::max(lsFromLeft, lsFromRight), lsFromRad);
}

moris::real
LevelSetFunction_star( const moris::Matrix< moris::DDRMat > & aPoint )
{
    moris::real tPhi = std::atan2( aPoint( 0 ), aPoint( 2 ) );

    moris::real tLevelSetVaue = 0.5 + 0.1 * std::sin( 5 * tPhi ) - std::sqrt( std::pow( aPoint( 0 ), 2 ) + std::pow( aPoint( 2 ), 2 ) );

    return tLevelSetVaue;
}

Matrix< DDRMat > tConstValFunction_MDL_XTK_HMR( moris::Cell< Matrix< DDRMat > >         & aCoeff,
                                                moris::Cell< fem::Field_Interpolator* > & aDofFI,
                                                moris::Cell< fem::Field_Interpolator* > & aDvFI,
                                                fem::Geometry_Interpolator             * aGeometryInterpolator )
{
    return aCoeff( 0 );
}

bool tSolverOutputCriteria( moris::tsa::Time_Solver * )
{
    return true;
}

TEST_CASE("HMR Interpolation STK Cut Diffusion Model Lag Order 2","[XTK_HMR_STK_DIFF]")
{
    if(par_size() == 1)
    {
        std::string tFieldName = "Circle";

        moris::uint tLagrangeMeshIndex = 0;
        moris::uint tBSplineMeshIndex = 0;

        moris::hmr::Parameters tParameters;

        tParameters.set_number_of_elements_per_dimension( { {4}, {4}, {4} } );
        tParameters.set_domain_dimensions({ {1}, {1}, {2} });
        tParameters.set_domain_offset({ {0.0}, {0.0}, {0.0} });
        tParameters.set_bspline_truncation( true );
        tParameters.set_side_sets({ {5}, {6} });

        tParameters.set_output_meshes( { {0} } );

        tParameters.set_lagrange_orders  ( { {1} });
        tParameters.set_lagrange_patterns({ {0} });

        tParameters.set_bspline_orders   ( { {1} } );
        tParameters.set_bspline_patterns ( { {0} } );

        tParameters.set_union_pattern( 2 );
        tParameters.set_working_pattern( 3 );

        tParameters.set_refinement_buffer( 2 );
        tParameters.set_staircase_buffer( 2);

        Cell< Matrix< DDUMat > > tLagrangeToBSplineMesh( 1 );
        tLagrangeToBSplineMesh( 0 ) = { {0} };

        tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

        hmr::HMR tHMR( tParameters );

        std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

        // create field
        std::shared_ptr< moris::hmr::Field > tField = tMesh->create_field( tFieldName, tLagrangeMeshIndex );

        tField->evaluate_scalar_function( LevelSetPlaneFunction );

        for( uint k=0; k<0; ++k )
        {
            tHMR.flag_surface_elements_on_working_pattern( tField );
            tHMR.perform_refinement_based_on_working_pattern( 0 );

            tField->evaluate_scalar_function( LevelSetPlaneFunction );
        }

        tHMR.finalize();

        tHMR.save_to_exodus( 0, "./mdl_exo/xtk_hmr_bar_hole_interp_l1_b1.e" );

        std::shared_ptr< hmr::Interpolation_Mesh_HMR > tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex  );

        moris::ge::GEN_Geom_Field tFieldAsGeom(tField);

        moris::Cell<moris::ge::GEN_Geometry*> tGeometryVector = {&tFieldAsGeom};

        // Tell the geometry engine about the discrete field mesh and how to interpret phases
        moris::ge::GEN_Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
        moris::ge::GEN_Geometry_Engine tGeometryEngine(tGeometryVector,tPhaseTable);

        // Tell the XTK model that it should decompose with a C_HIERARCHY_TET4, on the same mesh that the level set field is defined on.
        size_t tModelDimension = 3;
        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8,Subdivision_Method::C_HIERARCHY_TET4};
        xtk::Model tXTKModel(tModelDimension,tInterpMesh.get(),tGeometryEngine);
        tXTKModel.mSameMesh = true;
        tXTKModel.mVerbose = false;

        // Do the cutting
        tXTKModel.decompose(tDecompositionMethods);

        xtk::Output_Options tOutputOptions;
        tOutputOptions.mAddNodeSets = false;
        tOutputOptions.mAddSideSets = true;
        tOutputOptions.mAddClusters = true;

        // add solution field to integration mesh
        std::string tIntegSolFieldName = "solution";
        tOutputOptions.mRealNodeExternalFieldNames = {tIntegSolFieldName};

        moris::mtk::Integration_Mesh* tIntegMesh1 = tXTKModel.get_output_mesh(tOutputOptions);

        // place the pair in mesh manager
        mtk::Mesh_Manager tMeshManager;
        tMeshManager.register_mesh_pair(tInterpMesh.get(), tIntegMesh1);

        //------------------------------------------------------------------------------
        // create the properties
        std::shared_ptr< fem::Property > tPropConductivity = std::make_shared< fem::Property >();
        tPropConductivity->set_parameters( { {{ 1.0 }} } );
        tPropConductivity->set_val_function( tConstValFunction_MDL_XTK_HMR );

        std::shared_ptr< fem::Property > tPropDirichlet = std::make_shared< fem::Property >();
        tPropDirichlet->set_parameters( { {{ 5.0 }} } );
        tPropDirichlet->set_val_function( tConstValFunction_MDL_XTK_HMR );

        std::shared_ptr< fem::Property > tPropNeumann = std::make_shared< fem::Property >();
        tPropNeumann->set_parameters( { {{ 20.0 }} } );
        tPropNeumann->set_val_function( tConstValFunction_MDL_XTK_HMR );

        std::shared_ptr< fem::Property > tPropTempLoad = std::make_shared< fem::Property >();
        tPropTempLoad->set_parameters( { {{ 0.0 }} } );
        tPropTempLoad->set_val_function( tConstValFunction_MDL_XTK_HMR );

        // define constitutive models
        fem::CM_Factory tCMFactory;

        std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIso = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
        tCMDiffLinIso->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} ); // FIXME through the factory?
        tCMDiffLinIso->set_property( tPropConductivity, "Conductivity" );
        tCMDiffLinIso->set_space_dim( 3 );

        // define stabilization parameters
        fem::SP_Factory tSPFactory;
        std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche = tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
        tSPDirichletNitsche->set_parameters( { {{ 1.0 }} } );
        tSPDirichletNitsche->set_property( tPropConductivity, "Material", mtk::Master_Slave::MASTER );

        // define the IWGs
        fem::IWG_Factory tIWGFactory;

        std::shared_ptr< fem::IWG > tIWGBulk = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
        tIWGBulk->set_residual_dof_type( { MSI::Dof_Type::TEMP } );
        tIWGBulk->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tIWGBulk->set_constitutive_model( tCMDiffLinIso, "DiffLinIso", mtk::Master_Slave::MASTER );
        tIWGBulk->set_property( tPropTempLoad, "Load", mtk::Master_Slave::MASTER );

        std::shared_ptr< fem::IWG > tIWGDirichlet = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_DIRICHLET );
        tIWGDirichlet->set_residual_dof_type( { MSI::Dof_Type::TEMP } );
        tIWGDirichlet->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tIWGDirichlet->set_stabilization_parameter( tSPDirichletNitsche, "DirichletNitsche" );
        tIWGDirichlet->set_constitutive_model( tCMDiffLinIso, "DiffLinIso", mtk::Master_Slave::MASTER );
        tIWGDirichlet->set_property( tPropDirichlet, "Dirichlet", mtk::Master_Slave::MASTER );

        std::shared_ptr< fem::IWG > tIWGNeumann = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_NEUMANN );
        tIWGNeumann->set_residual_dof_type( { MSI::Dof_Type::TEMP } );
        tIWGNeumann->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tIWGNeumann->set_property( tPropNeumann, "Neumann", mtk::Master_Slave::MASTER );

        // define set info
        fem::Set_User_Info tSetBulk1;
        tSetBulk1.set_mesh_index( tIntegMesh1->get_set_index_by_name("child_0") );
        tSetBulk1.set_IWGs( { tIWGBulk } );

        fem::Set_User_Info tSetBulk2;
        tSetBulk2.set_mesh_index(  tIntegMesh1->get_set_index_by_name("parent_0") );
        tSetBulk2.set_IWGs( { tIWGBulk } );

        fem::Set_User_Info tSetDirichlet;
        tSetDirichlet.set_mesh_index( tIntegMesh1->get_set_index_by_name("iside_g_0_p0_0_p1_1") );
        tSetDirichlet.set_IWGs( { tIWGDirichlet } );

        fem::Set_User_Info tSetNeumann;
        tSetNeumann.set_mesh_index( tIntegMesh1->get_set_index_by_name("SideSet_1") );
        tSetNeumann.set_IWGs( { tIWGNeumann } );

        // create a cell of set info
        moris::Cell< fem::Set_User_Info > tSetInfo( 4 );
        tSetInfo( 0 ) = tSetBulk1;
        tSetInfo( 1 ) = tSetBulk2;
        tSetInfo( 2 ) = tSetDirichlet;
        tSetInfo( 3 ) = tSetNeumann;

        // create model
        mdl::Model * tModel = new mdl::Model( &tMeshManager,
                                               tBSplineMeshIndex,
                                               tSetInfo );

        moris::Cell< enum MSI::Dof_Type > tDofTypes1( 1, MSI::Dof_Type::TEMP );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 1: create linear solver and algorithm
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        dla::Solver_Factory  tSolFactory;
        std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( SolverType::AZTEC_IMPL );

        tLinearSolverAlgorithm->set_param("AZ_diagnostics") = AZ_none;
        tLinearSolverAlgorithm->set_param("AZ_output") = AZ_none;

        dla::Linear_Solver tLinSolver;

        tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 2: create nonlinear solver and algorithm
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        NLA::Nonlinear_Solver_Factory tNonlinFactory;
        std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );

//        tNonlinearSolverAlgorithm->set_param("NLA_max_iter")   = 10;
//        tNonlinearSolverAlgorithm->set_param("NLA_hard_break") = false;
//        tNonlinearSolverAlgorithm->set_param("NLA_max_lin_solver_restarts") = 2;
//        tNonlinearSolverAlgorithm->set_param("NLA_rebuild_jacobian") = true;

        tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );

        NLA::Nonlinear_Solver tNonlinearSolver;

        tNonlinearSolver.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 3: create time Solver and algorithm
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        tsa::Time_Solver_Factory tTimeSolverFactory;
        std::shared_ptr< tsa::Time_Solver_Algorithm > tTimeSolverAlgorithm = tTimeSolverFactory.create_time_solver( tsa::TimeSolverType::MONOLITHIC );

        tTimeSolverAlgorithm->set_nonlinear_solver( &tNonlinearSolver );

        tsa::Time_Solver tTimeSolver;

        tTimeSolver.set_time_solver_algorithm( tTimeSolverAlgorithm );

        NLA::SOL_Warehouse tSolverWarehouse;

        tSolverWarehouse.set_solver_interface(tModel->get_solver_interface());

        tNonlinearSolver.set_solver_warehouse( &tSolverWarehouse );
        tTimeSolver.set_solver_warehouse( &tSolverWarehouse );

        tNonlinearSolver.set_dof_type_list( tDofTypes1 );
        tTimeSolver.set_dof_type_list( tDofTypes1 );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 4: Solve and check
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        tTimeSolver.solve();


        // TODO: add gold solution data for this problem

        // Write to Integration mesh for visualization
        Matrix<DDRMat> tIntegSol = tModel->get_solution_for_integration_mesh_output( MSI::Dof_Type::TEMP );


        // add solution field to integration mesh
        tIntegMesh1->add_mesh_field_real_scalar_data_loc_inds(tIntegSolFieldName,EntityRank::NODE,tIntegSol);


        Matrix<DDRMat> tFullSol;
        tTimeSolver.get_full_solution(tFullSol);
//
        // verify solution
//        CHECK(norm(tSolution11 - tGoldSolution)<1e-08);
        tModel->output_solution( "Circle" );
        tField->put_scalar_values_on_field( tModel->get_mSolHMR() );
        tHMR.save_to_exodus( 0, "./mdl_exo/xtk_hmr_stk_bar_plane_interp_l2_b2.e" );


        // output solution and meshes
        std::string tMeshOutputFile = "./mdl_exo/xtk_hmr_stk_bar_hole_integ.e";
        tIntegMesh1->create_output_mesh(tMeshOutputFile);

        //    delete tInterpMesh1;
        delete tModel;
        delete tIntegMesh1;
    }
}

TEST_CASE("HMR Interpolation XTK Cut Diffusion Model Lag Order 2","[XTK_HMR_DIFF]")
{
    if(par_size() == 1)
    {
        std::string tFieldName = "Cylinder";


        moris::uint tLagrangeMeshIndex = 0;
        moris::uint tBSplineMeshIndex = 0;

        moris::hmr::Parameters tParameters;

        tParameters.set_number_of_elements_per_dimension( { {4}, {4}, {4} } );
        tParameters.set_domain_dimensions({ {1}, {1}, {2} });
        tParameters.set_domain_offset({ {0.0}, {0.0}, {0.0} });
        tParameters.set_bspline_truncation( true );
        tParameters.set_side_sets({ {5}, {6} });

        tParameters.set_output_meshes( { {0} } );

        tParameters.set_lagrange_orders  ( { {1} });
        tParameters.set_lagrange_patterns({ {0} });

        tParameters.set_bspline_orders   ( { {1} } );
        tParameters.set_bspline_patterns ( { {0} } );

        tParameters.set_union_pattern( 2 );
        tParameters.set_working_pattern( 3 );

        tParameters.set_refinement_buffer( 2 );
        tParameters.set_staircase_buffer( 2);

        tParameters.set_number_aura(  true );

        Cell< Matrix< DDUMat > > tLagrangeToBSplineMesh( 1 );
        tLagrangeToBSplineMesh( 0 ) = { {0} };

        tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

        hmr::HMR tHMR( tParameters );

        std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

        // create field
        std::shared_ptr< moris::hmr::Field > tField = tMesh->create_field( tFieldName, tLagrangeMeshIndex );

        tField->evaluate_scalar_function( LevelSetPlaneFunction );

        for( uint k=0; k<1; ++k )
        {
//            tHMR.finalize();
            tHMR.flag_surface_elements_on_working_pattern( tField );
            tHMR.perform_refinement_based_on_working_pattern( 0 );

            tField->evaluate_scalar_function( LevelSetPlaneFunction );
        }

        tHMR.finalize();


        tHMR.save_to_exodus( 0, "./mdl_exo/xtk_hmr_bar_plane_interp_l2_b2.e" );

        std::shared_ptr< hmr::Interpolation_Mesh_HMR > tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex  );

        moris::ge::GEN_Geom_Field tFieldAsGeom(tField);

        moris::Cell<moris::ge::GEN_Geometry*> tGeometryVector = {&tFieldAsGeom};

        // Tell the geometry engine about the discrete field mesh and how to interpret phases
        moris::ge::GEN_Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
        moris::ge::GEN_Geometry_Engine tGeometryEngine(tGeometryVector,tPhaseTable);

        // Tell the XTK model that it should decompose with a C_HIERARCHY_TET4, on the same mesh that the level set field is defined on.
        size_t tModelDimension = 3;
        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8,Subdivision_Method::C_HIERARCHY_TET4};
        xtk::Model tXTKModel(tModelDimension,tInterpMesh.get(),tGeometryEngine);
        tXTKModel.mSameMesh = true;
        tXTKModel.mVerbose = false;

        // Do the cutting
        tXTKModel.decompose(tDecompositionMethods);

        // Perform the enrichment
        tXTKModel.perform_basis_enrichment(EntityRank::BSPLINE,0);
//        tXTKModel.construct_face_oriented_ghost_penalization_cells();

        xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
        xtk::Enriched_Integration_Mesh   & tEnrIntegMesh = tXTKModel.get_enriched_integ_mesh();

//        tEnrInterpMesh.print_vertex_interpolation();

        // place the pair in mesh manager
        mtk::Mesh_Manager tMeshManager;
        tMeshManager.register_mesh_pair(&tEnrInterpMesh, &tEnrIntegMesh);

        //------------------------------------------------------------------------------
         // create the properties
         std::shared_ptr< fem::Property > tPropConductivity = std::make_shared< fem::Property >();
         tPropConductivity->set_parameters( { {{ 1.0 }} } );
         tPropConductivity->set_val_function( tConstValFunction_MDL_XTK_HMR );

         std::shared_ptr< fem::Property > tPropDirichlet = std::make_shared< fem::Property >();
         tPropDirichlet->set_parameters( { {{ 5.0 }} } );
         tPropDirichlet->set_val_function( tConstValFunction_MDL_XTK_HMR );

         std::shared_ptr< fem::Property > tPropNeumann = std::make_shared< fem::Property >();
         tPropNeumann->set_parameters( { {{ 20.0 }} } );
         tPropNeumann->set_val_function( tConstValFunction_MDL_XTK_HMR );

         std::shared_ptr< fem::Property > tPropTempLoad = std::make_shared< fem::Property >();
         tPropTempLoad->set_parameters( { {{ 0.0 }} } );
         tPropTempLoad->set_val_function( tConstValFunction_MDL_XTK_HMR );

         // define constitutive models
         fem::CM_Factory tCMFactory;

         std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIso = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
         tCMDiffLinIso->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} ); // FIXME through the factory?
         tCMDiffLinIso->set_property( tPropConductivity, "Conductivity" );
         tCMDiffLinIso->set_space_dim( 3 );

         // define stabilization parameters
         fem::SP_Factory tSPFactory;
         std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche = tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
         tSPDirichletNitsche->set_parameters( { {{ 1.0 }} } );
         tSPDirichletNitsche->set_property( tPropConductivity, "Material", mtk::Master_Slave::MASTER );

         // define the IWGs
         fem::IWG_Factory tIWGFactory;

         std::shared_ptr< fem::IWG > tIWGBulk = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
         tIWGBulk->set_residual_dof_type( { MSI::Dof_Type::TEMP } );
         tIWGBulk->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
         tIWGBulk->set_constitutive_model( tCMDiffLinIso, "DiffLinIso", mtk::Master_Slave::MASTER );
         tIWGBulk->set_property( tPropTempLoad, "Load", mtk::Master_Slave::MASTER );

         std::shared_ptr< fem::IWG > tIWGDirichlet = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_DIRICHLET );
         tIWGDirichlet->set_residual_dof_type( { MSI::Dof_Type::TEMP } );
         tIWGDirichlet->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
         tIWGDirichlet->set_stabilization_parameter( tSPDirichletNitsche, "DirichletNitsche" );
         tIWGDirichlet->set_constitutive_model( tCMDiffLinIso, "DiffLinIso", mtk::Master_Slave::MASTER );
         tIWGDirichlet->set_property( tPropDirichlet, "Dirichlet", mtk::Master_Slave::MASTER );

         std::shared_ptr< fem::IWG > tIWGNeumann = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_NEUMANN );
         tIWGNeumann->set_residual_dof_type( { MSI::Dof_Type::TEMP } );
         tIWGNeumann->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
         tIWGNeumann->set_property( tPropNeumann, "Neumann", mtk::Master_Slave::MASTER );

         // define the IQIs
         fem::IQI_Factory tIQIFactory;

         std::shared_ptr< fem::IQI > tIQITEMP = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
         tIQITEMP->set_output_type( vis::Output_Type::TEMP );
         tIQITEMP->set_dof_type_list( { {MSI::Dof_Type::TEMP} }, mtk::Master_Slave::MASTER );
         tIQITEMP->set_output_type_index( 0 );

         // define set info
         fem::Set_User_Info tSetBulk1;
         tSetBulk1.set_mesh_index( tEnrIntegMesh.get_set_index_by_name("HMR_dummy_c_p0") );
         tSetBulk1.set_IWGs( { tIWGBulk } );
         tSetBulk1.set_IQIs( { tIQITEMP } );

         fem::Set_User_Info tSetBulk2;
         tSetBulk2.set_mesh_index( tEnrIntegMesh.get_set_index_by_name("HMR_dummy_n_p0") );
         tSetBulk2.set_IWGs( { tIWGBulk } );
         tSetBulk2.set_IQIs( { tIQITEMP } );

         fem::Set_User_Info tSetDirichlet;
         std::string tInterfaceSideSetName = tEnrIntegMesh.get_interface_side_set_name( 0, 0, 1 );
         tSetDirichlet.set_mesh_index( tEnrIntegMesh.get_set_index_by_name(tInterfaceSideSetName) );
         tSetDirichlet.set_IWGs( { tIWGDirichlet } );

         fem::Set_User_Info tSetNeumann;
         tSetNeumann.set_mesh_index( tEnrIntegMesh.get_set_index_by_name("SideSet_1_n_p0") );
         tSetNeumann.set_IWGs( { tIWGNeumann } );

         // create a cell of set info
         moris::Cell< fem::Set_User_Info > tSetInfo( 4 );
         tSetInfo( 0 ) = tSetBulk1;
         tSetInfo( 1 ) = tSetBulk2;
         tSetInfo( 2 ) = tSetDirichlet;
         tSetInfo( 3 ) = tSetNeumann;

         // create model
         mdl::Model * tModel = new mdl::Model( &tMeshManager,
                                                tBSplineMeshIndex,
                                                tSetInfo );

         // --------------------------------------------------------------------------------------
         // Define outputs

         vis::Output_Manager tOutputData;

         tOutputData.set_outputs( 0,
//                                         VIS_Mesh_Type::STANDARD,
        		 vis::VIS_Mesh_Type::OVERLAPPING_INTERFACE,
                                  "XTK_HMR_DIFF.exo",
                                  { "HMR_dummy_c_p0", "HMR_dummy_n_p0"},
                                  { "Temp" },
                                  {   vis::Field_Type::NODAL },
                                  { vis::Output_Type::TEMP } );

         tModel->set_output_manager( &tOutputData );

         // --------------------------------------------------------------------------------------

        moris::Cell< enum MSI::Dof_Type > tDofTypes1( 1, MSI::Dof_Type::TEMP );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 1: create linear solver and algorithm
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        dla::Solver_Factory  tSolFactory;
        std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( SolverType::AZTEC_IMPL );

        tLinearSolverAlgorithm->set_param("AZ_diagnostics") = AZ_none;
        tLinearSolverAlgorithm->set_param("AZ_output") = AZ_none;

        dla::Linear_Solver tLinSolver;

        tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 2: create nonlinear solver and algorithm
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        NLA::Nonlinear_Solver_Factory tNonlinFactory;
        std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );

//        tNonlinearSolverAlgorithm->set_param("NLA_max_iter")   = 10;
//        tNonlinearSolverAlgorithm->set_param("NLA_hard_break") = false;
//        tNonlinearSolverAlgorithm->set_param("NLA_max_lin_solver_restarts") = 2;
//        tNonlinearSolverAlgorithm->set_param("NLA_rebuild_jacobian") = true;

        tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );

        NLA::Nonlinear_Solver tNonlinearSolver;

        tNonlinearSolver.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 3: create time Solver and algorithm
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        tsa::Time_Solver_Factory tTimeSolverFactory;
        std::shared_ptr< tsa::Time_Solver_Algorithm > tTimeSolverAlgorithm = tTimeSolverFactory.create_time_solver( tsa::TimeSolverType::MONOLITHIC );

        tTimeSolverAlgorithm->set_nonlinear_solver( &tNonlinearSolver );

        tsa::Time_Solver tTimeSolver;

        tTimeSolver.set_time_solver_algorithm( tTimeSolverAlgorithm );

        NLA::SOL_Warehouse tSolverWarehouse;

        tSolverWarehouse.set_solver_interface(tModel->get_solver_interface());

        tNonlinearSolver.set_solver_warehouse( &tSolverWarehouse );
        tTimeSolver.set_solver_warehouse( &tSolverWarehouse );

        tNonlinearSolver.set_dof_type_list( tDofTypes1 );
        tTimeSolver.set_dof_type_list( tDofTypes1 );

        tTimeSolver.set_output( 0, tSolverOutputCriteria );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 4: Solve and check
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        tTimeSolver.solve();


        // TODO: add gold solution data for this problem
//        Matrix<DDRMat> tFullSol;
//        tTimeSolver.get_full_solution(tFullSol);
//        print(tFullSol,"Full Solution");

        // verify solution
//        CHECK(norm(tSolution11 - tGoldSolution)<1e-08);

        // output solution and meshes
        xtk::Output_Options tOutputOptions;
        tOutputOptions.mAddNodeSets = false;
        tOutputOptions.mAddSideSets = true;
        tOutputOptions.mAddClusters = false;

        // add solution field to integration mesh
        std::string tIntegSolFieldName = "solution";
        tOutputOptions.mRealNodeExternalFieldNames = {tIntegSolFieldName};

        moris::mtk::Integration_Mesh* tIntegMesh1 = tXTKModel.get_output_mesh(tOutputOptions);

        // Write to Integration mesh for visualization
        Matrix<DDRMat> tIntegSol = tModel->get_solution_for_integration_mesh_output( MSI::Dof_Type::TEMP );

        Matrix<DDRMat> tSTKIntegSol(tIntegMesh1->get_num_entities(EntityRank::NODE),1);

        for(moris::uint i = 0; i < tIntegMesh1->get_num_entities(EntityRank::NODE); i++)
        {
            moris::moris_id tID = tIntegMesh1->get_glb_entity_id_from_entity_loc_index(i,EntityRank::NODE);
            tSTKIntegSol(i) = tIntegSol(tEnrIntegMesh.get_loc_entity_ind_from_entity_glb_id(tID,EntityRank::NODE));
        }

        // crate field in integration mesh
        moris::moris_index tFieldIndex = tEnrIntegMesh.create_field("Solution",EntityRank::NODE);
        tEnrIntegMesh.add_field_data(tFieldIndex,EntityRank::NODE,tSTKIntegSol);

        // add solution field to integration mesh
        tIntegMesh1->add_mesh_field_real_scalar_data_loc_inds(tIntegSolFieldName,EntityRank::NODE,tSTKIntegSol);


//        Matrix<DDRMat> tFullSol;
//        tTimeSolver.get_full_solution(tFullSol);

        std::string tMeshOutputFile = "./mdl_exo/xtk_hmr_bar_hole_integ.e";
        tIntegMesh1->create_output_mesh(tMeshOutputFile);

        delete tModel;
        delete tIntegMesh1;
    }
}



TEST_CASE("HMR Interpolation XTK Cut Diffusion Model Multigrid","[XTK_HMR_DIFF_MULTIGRID]")
{
//    //FIXME Timesolver for petsc
//
//    if( par_size() == 1 )
//    {
//        gLogger.set_severity_level( 0 );
//
//        std::string tFieldName = "Cylinder";
//
//        // start timer
//        tic tTimer_HMR;
//
//        moris::uint tLagrangeMeshIndex = 0;
//        moris::uint tBSplineMeshIndex = 0;
//
//        moris::hmr::Parameters tParameters;
//
//        tParameters.set_number_of_elements_per_dimension( { {2}, {2}, {4} } );
//        tParameters.set_domain_dimensions({ {2}, {2}, {4} });
//        tParameters.set_domain_offset({ {-1.0}, {-1.0}, {-2.0} });
//        tParameters.set_bspline_truncation( true );
//        tParameters.set_side_sets({ {5}, {6} });
//
//        tParameters.set_multigrid( true );
//
//        tParameters.set_output_meshes( { {0} } );
//
//        tParameters.set_lagrange_orders  ( { {1} });
//        tParameters.set_lagrange_patterns({ {0} });
//
//        tParameters.set_bspline_orders   ( { {1} } );
//        tParameters.set_bspline_patterns ( { {0} } );
//
//        tParameters.set_union_pattern( 2 );
//        tParameters.set_working_pattern( 3 );
//
//        tParameters.set_refinement_buffer( 2 );
//        tParameters.set_staircase_buffer( 2 );
//
//        Cell< Matrix< DDUMat > > tLagrangeToBSplineMesh( 1 );
//        tLagrangeToBSplineMesh( 0 ) = { {0} };
//
//        tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );
//
//        hmr::HMR tHMR( tParameters );
//
//        std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );
//
//        // create field
//        std::shared_ptr< moris::hmr::Field > tField = tMesh->create_field( tFieldName, tLagrangeMeshIndex );
//
//        tField->evaluate_scalar_function( LevelSetSphereCylinder );
//
//        for( uint k=0; k<2; ++k )
//        {
//            tHMR.flag_surface_elements_on_working_pattern( tField );
//            tHMR.perform_refinement_based_on_working_pattern( 0 );
//
//            tField->evaluate_scalar_function( LevelSetSphereCylinder );
//        }
//
//        tHMR.finalize();
//
//        // stop timer
//        real tElapsedTime = tTimer_HMR.toc<moris::chronos::milliseconds>().wall;
//
//        MORIS_LOG_INFO( " HMR took %5.3f seconds.\n", ( double ) tElapsedTime / 1000);
//
////        tHMR.save_to_exodus( "./mdl_exo/xtk_hmr_bar_hole_interp_l1_b1.e" );
//
//        std::shared_ptr< hmr::Interpolation_Mesh_HMR > tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex  );
//
//        // start timer
//        tic tTimer_XTK;
//
//        moris::ge::GEN_Geom_Field tFieldAsGeom(tField);
//
//        moris::Cell<moris::ge::GEN_Geometry*> tGeometryVector = {&tFieldAsGeom};
//
//        // Tell the geometry engine about the discrete field mesh and how to interpret phases
//        moris::ge::GEN_Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
//        moris::ge::GEN_Geometry_Engine tGeometryEngine(tGeometryVector,tPhaseTable);
//
//        // Tell the XTK model that it should decompose with a C_HIERARCHY_TET4, on the same mesh that the level set field is defined on.
//        size_t tModelDimension = 3;
//        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8,Subdivision_Method::C_HIERARCHY_TET4};
//        xtk::Model tXTKModel(tModelDimension,tInterpMesh.get(),tGeometryEngine);
//        tXTKModel.mSameMesh = true;
//        tXTKModel.mVerbose = false;
//
//        // Do the cutting
//        tXTKModel.decompose(tDecompositionMethods);
//
//        xtk::Output_Options tOutputOptions;
//        tOutputOptions.mAddNodeSets = false;
//        tOutputOptions.mAddSideSets = true;
//        tOutputOptions.mAddClusters = true;
//
//        // add solution field to integration mesh
//        std::string tIntegSolFieldName = "solution";
//        tOutputOptions.mRealNodeExternalFieldNames = {tIntegSolFieldName};
//
//        moris::mtk::Integration_Mesh* tIntegMesh1 = tXTKModel.get_output_mesh(tOutputOptions);
//
//        // place the pair in mesh manager
//        mtk::Mesh_Manager tMeshManager;
//        tMeshManager.register_mesh_pair(tInterpMesh.get(), tIntegMesh1);
//
//        //------------------------------------------------------------------------------
//         // create the properties
//         std::shared_ptr< fem::Property > tPropConductivity = std::make_shared< fem::Property >();
//         tPropConductivity->set_parameters( { {{ 1.0 }} } );
//         tPropConductivity->set_val_function( tConstValFunction_MDL_XTK_HMR );
//
//         std::shared_ptr< fem::Property > tPropDirichlet = std::make_shared< fem::Property >();
//         tPropDirichlet->set_parameters( { {{ 5.0 }} } );
//         tPropDirichlet->set_val_function( tConstValFunction_MDL_XTK_HMR );
//
//         std::shared_ptr< fem::Property > tPropNeumann = std::make_shared< fem::Property >();
//         tPropNeumann->set_parameters( { {{ 20.0 }} } );
//         tPropNeumann->set_val_function( tConstValFunction_MDL_XTK_HMR );
//
//         std::shared_ptr< fem::Property > tPropTempLoad = std::make_shared< fem::Property >();
//         tPropTempLoad->set_parameters( { {{ 0.0 }} } );
//         tPropTempLoad->set_val_function( tConstValFunction_MDL_XTK_HMR );
//
//         // define constitutive models
//         fem::CM_Factory tCMFactory;
//
//         std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIso = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
//         tCMDiffLinIso->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} ); // FIXME through the factory?
//         tCMDiffLinIso->set_property( tPropConductivity, "Conductivity" );
//         tCMDiffLinIso->set_space_dim( 3 );
//
//         // define stabilization parameters
//         fem::SP_Factory tSPFactory;
//         std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche = tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
//         tSPDirichletNitsche->set_parameters( { {{ 1.0 }} } );
//         tSPDirichletNitsche->set_property( tPropConductivity, "Material", mtk::Master_Slave::MASTER );
//
//         // define the IWGs
//         fem::IWG_Factory tIWGFactory;
//
//         std::shared_ptr< fem::IWG > tIWGBulk = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
//         tIWGBulk->set_residual_dof_type( { MSI::Dof_Type::TEMP } );
//         tIWGBulk->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
//         tIWGBulk->set_constitutive_model( tCMDiffLinIso, "DiffLinIso", mtk::Master_Slave::MASTER );
//         tIWGBulk->set_property( tPropTempLoad, "Load", mtk::Master_Slave::MASTER );
//
//         std::shared_ptr< fem::IWG > tIWGDirichlet = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_DIRICHLET );
//         tIWGDirichlet->set_residual_dof_type( { MSI::Dof_Type::TEMP } );
//         tIWGDirichlet->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
//         tIWGDirichlet->set_stabilization_parameter( tSPDirichletNitsche, "DirichletNitsche" );
//         tIWGDirichlet->set_constitutive_model( tCMDiffLinIso, "DiffLinIso", mtk::Master_Slave::MASTER );
//         tIWGDirichlet->set_property( tPropDirichlet, "Dirichlet", mtk::Master_Slave::MASTER );
//
//         std::shared_ptr< fem::IWG > tIWGNeumann = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_NEUMANN );
//         tIWGNeumann->set_residual_dof_type( { MSI::Dof_Type::TEMP } );
//         tIWGNeumann->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
//         tIWGNeumann->set_property( tPropNeumann, "Neumann", mtk::Master_Slave::MASTER );
//
//         // define set info
//         fem::Set_User_Info tSetBulk1;
//         tSetBulk1.set_mesh_index( 4 );
//         tSetBulk1.set_IWGs( { tIWGBulk } );
//
//         fem::Set_User_Info tSetBulk2;
//         tSetBulk2.set_mesh_index( 5 );
//         tSetBulk2.set_IWGs( { tIWGBulk } );
//
//         fem::Set_User_Info tSetDirichlet1;
//         tSetDirichlet1.set_mesh_index( 1 );
//         tSetDirichlet1.set_IWGs( { tIWGDirichlet } );
//
//         fem::Set_User_Info tSetDirichlet2;
//         tSetDirichlet2.set_mesh_index( 3 );
//         tSetDirichlet2.set_IWGs( { tIWGDirichlet } );
//
//         fem::Set_User_Info tSetNeumann;
//         tSetNeumann.set_mesh_index( 0 );
//         tSetNeumann.set_IWGs( { tIWGNeumann } );
//
//         // create a cell of set info
//         moris::Cell< fem::Set_User_Info > tSetInfo( 5 );
//         tSetInfo( 0 ) = tSetBulk1;
//         tSetInfo( 1 ) = tSetBulk2;
//         tSetInfo( 2 ) = tSetDirichlet1;
//         tSetInfo( 3 ) = tSetDirichlet2;
//         tSetInfo( 4 ) = tSetNeumann;
//
//         // create model
//         mdl::Model * tModel = new mdl::Model( &tMeshManager,
//                                                tBSplineMeshIndex,
//                                                tSetInfo,
//                                                0, true );
//
//        // stop timer
//        real tElapsedTime1 = tTimer_XTK.toc<moris::chronos::milliseconds>().wall;
//
//        MORIS_LOG_INFO( " XTK took %5.3f seconds.\n", ( double ) tElapsedTime1 / 1000);
//
//        moris::Cell< enum MSI::Dof_Type > tDofTypes1( 1, MSI::Dof_Type::TEMP );
//
//        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//        // STEP 1: create linear solver and algorithm
//        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//
//        dla::Solver_Factory  tSolFactory;
//
//        // create linear solver
//        std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( SolverType::PETSC );
//
//        tLinearSolverAlgorithm->set_param("KSPType") = std::string( KSPFGMRES );
//        tLinearSolverAlgorithm->set_param("PCType")  = std::string( PCMG );
////        tLinearSolverAlgorithm->set_param("PCType")  = std::string( PCILU );
//        tLinearSolverAlgorithm->set_param("ILUFill")  = 3;
//
//        dla::Linear_Solver tLinSolver;
//
//        tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );
//
//        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//        // STEP 2: create nonlinear solver and algorithm
//        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//        NLA::Nonlinear_Problem * tNonlinearProblem =  new NLA::Nonlinear_Problem( tModel->get_solver_interface(), 0, true, MapType::Petsc );
//
//        // create factory for nonlinear solver
//        NLA::Nonlinear_Solver_Factory tNonlinFactory;
//
//        // create nonlinear solver
//        std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );
//
//        NLA::Nonlinear_Solver  tNonlinearSolver;
//
//        tNonlinearSolverAlgorithm->set_param("NLA_max_iter")                = 10;
//        tNonlinearSolverAlgorithm->set_param("NLA_hard_break")              = false;
//        tNonlinearSolverAlgorithm->set_param("NLA_max_lin_solver_restarts") = 2;
//        tNonlinearSolverAlgorithm->set_param("NLA_rebuild_jacobian")        = true;
//
//        // set manager and settings
//        tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );
//
//        tNonlinearSolver.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );
//
//        tNonlinearSolver.solve( tNonlinearProblem );
//
//        // temporary array for solver
//        Matrix< DDRMat > tSolution;
//        tNonlinearSolverAlgorithm->get_full_solution( tSolution );
//
//        CHECK( equal_to( tSolution( 0, 0 ), 4.991691079008517, 1.0e+08 ) );
//        CHECK( equal_to( tSolution( 1, 0 ), 5.000966490616513, 1.0e+08 ) );
//        CHECK( equal_to( tSolution( 2, 0 ), 4.991691079008363, 1.0e+08 ) );
//        CHECK( equal_to( tSolution( 3, 0 ), 5.000966490616352, 1.0e+08 ) );
//        CHECK( equal_to( tSolution( 4, 0 ), 17.06241419809754, 1.0e+08 ) );
//        CHECK( equal_to( tSolution( 5, 0 ), 17.02038634782764, 1.0e+08 ) );
//        CHECK( equal_to( tSolution( 6, 0 ), 17.06241419809900, 1.0e+08 ) );
//        CHECK( equal_to( tSolution( 7, 0 ), 17.02038634782907, 1.0e+08 ) );
//        CHECK( equal_to( tSolution( 8, 0 ), 5.006375939765809, 1.0e+08 ) );
//        CHECK( equal_to( tSolution( 382, 0 ), 51.10753242793825, 1.0e+08 ) );
//        CHECK( equal_to( tSolution( 461, 0 ), 17.06241419810062, 1.0e+08 ) );
//        CHECK( equal_to( tSolution( 505, 0 ), 29.33523440842293, 1.0e+08 ) );

////        // start timer
////        tic tTimer_XTK1;
////
////        // Write to Integration mesh for visualization
////        Matrix<DDRMat> tIntegSol = tModel->get_solution_for_integration_mesh_output( MSI::Dof_Type::TEMP );
////
////        // add solution field to integration mesh
////        tIntegMesh1->add_mesh_field_real_scalar_data_loc_inds(tIntegSolFieldName,EntityRank::NODE,tIntegSol);
////
////        // output solution and meshes
////        std::string tMeshOutputFile = "./mdl_exo/xtk_hmr_bar_hole_integ.e";
////        tIntegMesh1->create_output_mesh(tMeshOutputFile);
////
////        // stop timer
////        real tElapsedTime2 = tTimer_XTK1.toc<moris::chronos::milliseconds>().wall;
////
////        MORIS_LOG_INFO( " output took %5.3f seconds.\n", ( double ) tElapsedTime2 / 1000);
//
//        //    delete tInterpMesh1;
//        delete tNonlinearProblem;
//        delete tModel;
//        delete tIntegMesh1;
//    }
}

//TEST_CASE("HMR Interpolation XTK Cut Diffusion Model Multigrid Star","[XTK_HMR_DIFF_STAR]")
//{
//    if( par_size() == 1 )
//    {
//        gLogger.set_severity_level( 0 );
//
//        std::string tFieldName = "Cylinder";
//
//        // start timer
//        tic tTimer_HMR;
//
//        moris::uint tLagrangeMeshIndex = 0;
//        moris::uint tBSplineMeshIndex = 0;
//
//        moris::hmr::Parameters tParameters;
//
//        tParameters.set_number_of_elements_per_dimension( { {2}, {1}, {2} } );
//        tParameters.set_domain_dimensions({ {2}, {2}, {2} });
//        tParameters.set_domain_offset({ {-1.0}, {-1.0}, {-1.0} });
//        tParameters.set_bspline_truncation( true );
//        tParameters.set_side_sets({ {5}, {6} });
//
//        tParameters.set_multigrid( true );
//
//        tParameters.set_output_meshes( { {0} } );
//
//        tParameters.set_lagrange_orders  ( { {1} });
//        tParameters.set_lagrange_patterns({ {0} });
//
//        tParameters.set_bspline_orders   ( { {1} } );
//        tParameters.set_bspline_patterns ( { {0} } );
//
//        tParameters.set_union_pattern( 2 );
//        tParameters.set_working_pattern( 3 );
//
//        tParameters.set_refinement_buffer( 1 );
//        tParameters.set_staircase_buffer( 2 );
//
//        Cell< Matrix< DDUMat > > tLagrangeToBSplineMesh( 1 );
//        tLagrangeToBSplineMesh( 0 ) = { {0} };
//
//        tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );
//
//        hmr::HMR tHMR( tParameters );
//
//        std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );
//
//        // create field
//        std::shared_ptr< moris::hmr::Field > tField = tMesh->create_field( tFieldName, tLagrangeMeshIndex );
//
//        tField->evaluate_scalar_function( LevelSetFunction_star );
//
//        for( uint k=0; k<4; ++k )
//        {
//            tHMR.flag_surface_elements_on_working_pattern( tField );
//            tHMR.perform_refinement_based_on_working_pattern( 0 );
//
//            tField->evaluate_scalar_function( LevelSetFunction_star );
//        }
//
//        tHMR.finalize();
//
//        tHMR.save_to_exodus( 0, "Mesh_1111.exo" );
//
//        // stop timer
//        real tElapsedTime = tTimer_HMR.toc<moris::chronos::milliseconds>().wall;
//
//        MORIS_LOG_INFO( " HMR took %5.3f seconds.\n", ( double ) tElapsedTime / 1000);
//
////        tHMR.save_to_exodus( "./mdl_exo/xtk_hmr_bar_hole_interp_l1_b1.e" );
//
//        std::shared_ptr< hmr::Interpolation_Mesh_HMR > tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex  );
//
//        // start timer
//        tic tTimer_XTK;
//
//        xtk::Geom_Field tFieldAsGeom(tField);
//
//        moris::Cell<xtk::Geometry*> tGeometryVector = {&tFieldAsGeom};
//
//        // Tell the geometry engine about the discrete field mesh and how to interpret phases
//        xtk::Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
//        xtk::Geometry_Engine tGeometryEngine(tGeometryVector,tPhaseTable);
//
//        // Tell the XTK model that it should decompose with a C_HIERARCHY_TET4, on the same mesh that the level set field is defined on.
//        size_t tModelDimension = 3;
//        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8,Subdivision_Method::C_HIERARCHY_TET4};
//        xtk::Model tXTKModel(tModelDimension,tInterpMesh.get(),tGeometryEngine);
//        tXTKModel.mSameMesh = true;
//        tXTKModel.mVerbose = false;
//
//        // Do the cutting
//        tXTKModel.decompose(tDecompositionMethods);
//
//        xtk::Output_Options tOutputOptions;
//        tOutputOptions.mAddNodeSets = false;
//        tOutputOptions.mAddSideSets = true;
//        tOutputOptions.mAddClusters = true;
//
//        // add solution field to integration mesh
//        std::string tIntegSolFieldName = "solution";
//        tOutputOptions.mRealNodeExternalFieldNames = {tIntegSolFieldName};
//
//        moris::mtk::Integration_Mesh* tIntegMesh1 = tXTKModel.get_output_mesh(tOutputOptions);
//
//        // place the pair in mesh manager
//        mtk::Mesh_Manager tMeshManager;
//        tMeshManager.register_mesh_pair(tInterpMesh.get(), tIntegMesh1);
//
//        // stop timer
//        real tElapsedTime1 = tTimer_XTK.toc<moris::chronos::milliseconds>().wall;
//
//        MORIS_LOG_INFO( " XTK took %5.3f seconds.\n", ( double ) tElapsedTime1 / 1000);
//
//        // create a list of IWG type
//        Cell< Cell< fem::IWG_Type > >tIWGTypeList( 5 );
//        tIWGTypeList( 0 ).resize( 1, fem::IWG_Type::SPATIALDIFF_BULK );
//        tIWGTypeList( 1 ).resize( 1, fem::IWG_Type::SPATIALDIFF_BULK );
//        tIWGTypeList( 2 ).resize( 1, fem::IWG_Type::SPATIALDIFF_DIRICHLET );
//        tIWGTypeList( 3 ).resize( 1, fem::IWG_Type::SPATIALDIFF_DIRICHLET );
//        tIWGTypeList( 4 ).resize( 1, fem::IWG_Type::SPATIALDIFF_NEUMANN );
//
//        // number of groups of IWgs
//        uint tNumSets = tIWGTypeList.size();
//
//        // list of residual dof type
//        moris::Cell< moris::Cell< moris::Cell< MSI::Dof_Type > > > tResidualDofType( tNumSets );
//        tResidualDofType( 0 ).resize( tIWGTypeList( 0 ).size(), { MSI::Dof_Type::TEMP } );
//        tResidualDofType( 1 ).resize( tIWGTypeList( 1 ).size(), { MSI::Dof_Type::TEMP } );
//        tResidualDofType( 2 ).resize( tIWGTypeList( 2 ).size(), { MSI::Dof_Type::TEMP } );
//        tResidualDofType( 3 ).resize( tIWGTypeList( 3 ).size(), { MSI::Dof_Type::TEMP } );
//        tResidualDofType( 4 ).resize( tIWGTypeList( 4 ).size(), { MSI::Dof_Type::TEMP } );
//
//        // list of IWG master dof dependencies
//        moris::Cell< moris::Cell< moris::Cell< moris::Cell< MSI::Dof_Type > > > > tMasterDofTypes( tNumSets );
//        tMasterDofTypes( 0 ).resize( tIWGTypeList( 0 ).size(), {{ MSI::Dof_Type::TEMP }} );
//        tMasterDofTypes( 1 ).resize( tIWGTypeList( 1 ).size(), {{ MSI::Dof_Type::TEMP }} );
//        tMasterDofTypes( 2 ).resize( tIWGTypeList( 2 ).size(), {{ MSI::Dof_Type::TEMP }} );
//        tMasterDofTypes( 3 ).resize( tIWGTypeList( 3 ).size(), {{ MSI::Dof_Type::TEMP }} );
//        tMasterDofTypes( 4 ).resize( tIWGTypeList( 4 ).size(), {{ MSI::Dof_Type::TEMP }} );
//
//        // list of IWG master property dependencies
//        moris::Cell< moris::Cell< moris::Cell< fem::Property_Type > > > tMasterPropTypes( tNumSets );
//        tMasterPropTypes( 0 ).resize( tIWGTypeList( 0 ).size(), { fem::Property_Type::CONDUCTIVITY } );
//        tMasterPropTypes( 1 ).resize( tIWGTypeList( 1 ).size(), { fem::Property_Type::CONDUCTIVITY } );
//        tMasterPropTypes( 2 ).resize( tIWGTypeList( 2 ).size(), { fem::Property_Type::CONDUCTIVITY, fem::Property_Type::TEMP_DIRICHLET } );
//        tMasterPropTypes( 3 ).resize( tIWGTypeList( 3 ).size(), { fem::Property_Type::CONDUCTIVITY, fem::Property_Type::TEMP_DIRICHLET } );
//        tMasterPropTypes( 4 ).resize( tIWGTypeList( 4 ).size(), { fem::Property_Type::TEMP_NEUMANN } );
//
//        // build an IWG user defined info
//        fem::IWG_User_Defined_Info tIWGUserDefinedInfo( tIWGTypeList,
//                                                        tResidualDofType,
//                                                        tMasterDofTypes,
//                                                        tMasterPropTypes );
//
//        // create a list of active block-sets
//        moris::Cell< moris_index >  tSetList = { 4, 5, 1, 3, 0 };
//
//        moris::Cell< fem::Element_Type > tSetTypeList = { fem::Element_Type::BULK,
//                                                          fem::Element_Type::BULK,
//                                                          fem::Element_Type::SIDESET,
//                                                          fem::Element_Type::SIDESET,
//                                                          fem::Element_Type::SIDESET };
//
//        // list of property type
//        Cell< fem::Property_Type > tPropertyTypeList = {{ fem::Property_Type::CONDUCTIVITY   },
//                                                        { fem::Property_Type::TEMP_DIRICHLET },
//                                                        { fem::Property_Type::TEMP_NEUMANN   }};
//
//        // list of property dependencies
//        Cell< Cell< Cell< MSI::Dof_Type > > > tPropertyDofList( 3 );
//
//        // list of the property coefficients
//        Cell< Cell< Matrix< DDRMat > > > tCoeffList( 3 );
//        tCoeffList( 0 ).resize( 1 );
//        tCoeffList( 0 )( 0 )= {{ 1.0 }};
//        tCoeffList( 1 ).resize( 1 );
//        tCoeffList( 1 )( 0 )= {{ 5.0 }};
//        tCoeffList( 2 ).resize( 1 );
//        tCoeffList( 2 )( 0 )= {{ 20.0 }};
//
//        // cast free function into std::function
//        fem::PropertyFunc tValFunction0 = tConstValFunction_MDL_XTK_HMR;
//
//        // create the list with function pointers for the value
//        Cell< fem::PropertyFunc > tValFuncList( 3, tValFunction0 );
//
//        // create the list with cell of function pointers for the derivatives
//        Cell< Cell< fem::PropertyFunc > > tDerFuncList( 3 );
//
//        // collect properties info
//        fem::Property_User_Defined_Info tPropertyUserDefinedInfo( tPropertyTypeList,
//                                                                  tPropertyDofList,
//                                                                  tCoeffList,
//                                                                  tValFuncList,
//                                                                  tDerFuncList );
//
//        // create model
//        mdl::Model * tModel = new mdl::Model( &tMeshManager,
//                                              tBSplineMeshIndex,
//                                              &tIWGUserDefinedInfo,
//                                              tSetList, tSetTypeList,
//                                              &tPropertyUserDefinedInfo,
//                                              0,
//                                              true);
//
//        moris::Cell< enum MSI::Dof_Type > tDofTypes1( 1, MSI::Dof_Type::TEMP );
//
//        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//        // STEP 1: create linear solver and algorithm
//        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//
//        dla::Solver_Factory  tSolFactory;
//
//        // create linear solver
//        std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( SolverType::PETSC );
//
//        tLinearSolverAlgorithm->set_param("KSPType") = std::string( KSPFGMRES );
//        tLinearSolverAlgorithm->set_param("PCType")  = std::string( PCMG );
////        tLinearSolverAlgorithm->set_param("PCType")  = std::string( PCILU );
//        tLinearSolverAlgorithm->set_param("ILUFill")  = 3;
//
//        dla::Linear_Solver tLinSolver;
//
//        tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );
//
//        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//        // STEP 2: create nonlinear solver and algorithm
//        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//        NLA::Nonlinear_Problem * tNonlinearProblem =  new NLA::Nonlinear_Problem( tModel->get_solver_interface(), 0, true, MapType::Petsc );
//
//        // create factory for nonlinear solver
//        NLA::Nonlinear_Solver_Factory tNonlinFactory;
//
//        // create nonlinear solver
//        std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );
//
//        NLA::Nonlinear_Solver  tNonlinearSolver;
//
//        tNonlinearSolverAlgorithm->set_param("NLA_max_iter")                = 10;
//        tNonlinearSolverAlgorithm->set_param("NLA_hard_break")              = false;
//        tNonlinearSolverAlgorithm->set_param("NLA_max_lin_solver_restarts") = 2;
//        tNonlinearSolverAlgorithm->set_param("NLA_rebuild_jacobian")        = true;
//
//        // set manager and settings
//        tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );
//
//        tNonlinearSolver.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );
//
//        tNonlinearSolver.solve( tNonlinearProblem );
//
//        // temporary array for solver
//        Matrix< DDRMat > tSolution;
//        tNonlinearSolverAlgorithm->get_full_solution( tSolution );
//
////        CHECK( equal_to( tSolution( 0, 0 ), 4.991691079008517, 1.0e+08 ) );
//
//
//        // start timer
//        tic tTimer_XTK1;
//
//        // Write to Integration mesh for visualization
//        Matrix<DDRMat> tIntegSol = tModel->get_solution_for_integration_mesh_output( MSI::Dof_Type::TEMP );
//
//        // add solution field to integration mesh
//        tIntegMesh1->add_mesh_field_real_scalar_data_loc_inds(tIntegSolFieldName,EntityRank::NODE,tIntegSol);
//
//        // output solution and meshes
//        std::string tMeshOutputFile = "./mdl_exo/xtk_hmr_bar_hole_integ.e";
//        tIntegMesh1->create_output_mesh(tMeshOutputFile);
//
//        // stop timer
//        real tElapsedTime2 = tTimer_XTK1.toc<moris::chronos::milliseconds>().wall;
//
//        MORIS_LOG_INFO( " output took %5.3f seconds.\n", ( double ) tElapsedTime2 / 1000);
//
//        //    delete tInterpMesh1;
//        delete tNonlinearProblem;
//        delete tModel;
//        delete tIntegMesh1;
//    }
//}

}


