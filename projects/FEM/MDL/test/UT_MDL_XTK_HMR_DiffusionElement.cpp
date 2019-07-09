/*
 * UT_MDL_XTK_HMR_DiffusionElement.cpp
 *
 *  Created on: Jun 27, 2019
 *      Author: doble
 */

#include "catch.hpp"

#include "cl_XTK_Model.hpp"

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

#include "cl_Matrix.hpp"        //LINALG
#include "linalg_typedefs.hpp"
#include "fn_equal_to.hpp" // ALG/src

#include "cl_FEM_NodeProxy.hpp"                //FEM/INT/src
#include "cl_FEM_ElementProxy.hpp"             //FEM/INT/src
#include "cl_FEM_Node_Base.hpp"                //FEM/INT/src
#include "cl_FEM_Element_Factory.hpp"          //FEM/INT/src
#include "cl_FEM_IWG_Factory.hpp"              //FEM/INT/src

#include "cl_MDL_Model.hpp"

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


namespace moris
{
moris::real
LevelSetPlaneFunction( const moris::Matrix< moris::DDRMat > & aPoint )
{

    real mXn = 0;
    real mYn = 0;
    real mZn = 1.0;
    real mXc = 1.0;
    real mYc = 1.0;
    real mZc = 3.51;
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

TEST_CASE("HMR Interpolation XTK Cut Diffusion Model Lag Order 1","[XTK_HMR_DIFF]")
{
    if(par_size() == 1)
    {
        moris::uint tBplineOrder = 1;
        moris::uint tLagrangeOrder = 1;
        moris::uint tMyCoeff = 1;

        hmr::ParameterList tParameters = hmr::create_hmr_parameter_list();

        tParameters.set( "number_of_elements_per_dimension", "2, 2, 4" );
        tParameters.set( "domain_dimensions", "2, 2, 4" );
        tParameters.set( "domain_offset", "-1.0, -1.0, -2.0" );
        tParameters.set( "domain_sidesets", "5,6");
        tParameters.set( "verbose", 0 );
        tParameters.set( "truncate_bsplines", 1 );
        tParameters.set( "bspline_orders", "1" );
        tParameters.set( "lagrange_orders", "1" );
        tParameters.set( "use_multigrid", 0 );
        tParameters.set( "refinement_buffer", 2 );
        tParameters.set( "staircase_buffer", 2 );

        hmr::HMR tHMR( tParameters );

        std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeOrder );

        // create field
        std::shared_ptr< moris::hmr::Field > tField = tMesh->create_field( "Sphere", tBplineOrder );

        tField->evaluate_scalar_function( LevelSetSphereCylinder );

        for( uint k=0; k<1; ++k )
        {
            tHMR.flag_surface_elements( tField );
            tHMR.perform_refinement( moris::hmr::RefinementMode::SIMPLE );
            tHMR.update_refinement_pattern();

            tField->evaluate_scalar_function( LevelSetSphereCylinder );
        }

        tHMR.finalize();

        tHMR.save_to_exodus( "./mdl_exo/xtk_hmr_bar_mesh_bm_bo2.e" );

        std::shared_ptr< hmr::Interpolation_Mesh_HMR > tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeOrder, tHMR.mParameters->get_lagrange_output_pattern()  );

        xtk::Geom_Field tFieldAsGeom(tField);

        moris::Cell<xtk::Geometry*> tGeometryVector = {&tFieldAsGeom};

        // Tell the geometry engine about the discrete field mesh and how to interpret phases
        xtk::Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
        xtk::Geometry_Engine tGeometryEngine(tGeometryVector,tPhaseTable);

        // Tell the XTK model that it should decompose with a C_HIERARCHY_TET4, on the same mesh that the level set field is defined on.
        size_t tModelDimension = 3;
        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8,Subdivision_Method::C_HIERARCHY_TET4};
        xtk::Model tXTKModel(tModelDimension,tInterpMesh.get(),tGeometryEngine);
        tXTKModel.mSameMesh = true;
        tXTKModel.mVerbose = true;

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

        std::string tMeshOutputFile = "./mdl_exo/xtk_hmr_bar_mesh_no_sol.e";
        tIntegMesh1->create_output_mesh(tMeshOutputFile);

        // place the pair in mesh manager
        mtk::Mesh_Manager tMeshManager;
        tMeshManager.register_mesh_pair(tInterpMesh.get(), tIntegMesh1);

        // create a list of IWG type
        Cell< Cell< fem::IWG_Type > >tIWGTypeList( 4 );
        tIWGTypeList( 0 ).resize( 1, fem::IWG_Type::SPATIALDIFF_BULK );
        tIWGTypeList( 1 ).resize( 1, fem::IWG_Type::SPATIALDIFF_DIRICHLET );
        tIWGTypeList( 2 ).resize( 1, fem::IWG_Type::SPATIALDIFF_NEUMANN );

        // create a list of active block-sets
        moris::Cell< moris_index >  tBlocksetList = { 4, 5 };

        // create a list of active side-sets
        moris::Cell< moris_index >  tSidesetList = { 1, 3 };

        std::cout<<"Set name 1 = "<<tIntegMesh1->get_side_set_label(1)<<std::endl;
        std::cout<<"Set name 3 = "<<tIntegMesh1->get_side_set_label(3)<<std::endl;

        // create a list of BC type for the side-sets
        moris::Cell< fem::BC_Type > tSidesetBCTypeList = { fem::BC_Type::DIRICHLET,
                                                           fem::BC_Type::NEUMANN };

        // create a list of active double side-sets
        moris::Cell< moris_index >  tDoubleSidesetList = {  };

        // create model
        mdl::Model * tModel = new mdl::Model( &tMeshManager, tBplineOrder, tIWGTypeList,
                                              tBlocksetList, tSidesetList,
                                              tSidesetBCTypeList,
                                              tDoubleSidesetList );

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

        tNonlinearSolverAlgorithm->set_param("NLA_max_iter")   = 10;
        tNonlinearSolverAlgorithm->set_param("NLA_hard_break") = false;
        tNonlinearSolverAlgorithm->set_param("NLA_max_lin_solver_restarts") = 2;
        tNonlinearSolverAlgorithm->set_param("NLA_rebuild_jacobian") = true;

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

        moris::Matrix< DDRMat > tSolution11;
        tTimeSolver.get_full_solution( tSolution11 );


        Matrix<DDRMat> tGoldSolution =  {{+2.50e+01},
                                         {+2.50e+01},
                                         {+2.50e+01},
                                         {+2.50e+01},
                                         {+4.50e+01},
                                         {+4.50e+01},
                                         {+4.50e+01},
                                         {+4.50e+01},
                                         {+6.50e+01},
                                         {+6.50e+01},
                                         {+6.50e+01},
                                         {+6.50e+01},
                                         {+8.50e+01},
                                         {+8.50e+01},
                                         {+8.50e+01},
                                         {+8.50e+01},
                                         {+5.00e+00},
                                         {+5.00e+00},
                                         {+5.00e+00},
                                         {+5.00e+00}};

        moris::print_fancy(tSolution11,"tSolution11");
        std::cout<<"Min = " <<tSolution11.min()<<std::endl;
        std::cout<<"Max = " <<tSolution11.max()<<std::endl;

        Matrix<DDRMat> tIntegSol = tModel->get_solution_for_integration_mesh_output( MSI::Dof_Type::TEMP );


        // add solution field to integration mesh
        tIntegMesh1->add_mesh_field_real_scalar_data_loc_inds(tIntegSolFieldName,EntityRank::NODE,tIntegSol);


        // verify solution
        //    CHECK(norm(tSolution11 - tGoldSolution)<1e-08);

        // output solution and meshes

        //    tInterpMesh1->add_mesh_field_real_scalar_data_loc_inds(tFieldName1,EntityRank::NODE,tSolution11);
        //
        //    std::string tOutputInterp = "./mdl_exo/xtk_mdl_interp.exo";
        //    tInterpMesh1->create_output_mesh(tOutputInterp);

        tMeshOutputFile = "./mdl_exo/xtk_hmr_bar_mesh.e";
        tIntegMesh1->create_output_mesh(tMeshOutputFile);

        //    delete tInterpMesh1;
        delete tModel;
        delete tIntegMesh1;
    }
}


}


