/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MDL_Diffusion_Phase_Change.cpp
 *
 */

#include "catch.hpp"

#include "cl_MTK_Vertex.hpp"    //MTK
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh.hpp"

#include "cl_MTK_Mesh_Factory.hpp"
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
#include "fn_inv.hpp" // to allow for matrix inverse call

#include <stdlib.h> // to include random number generator

#include "cl_FEM_NodeProxy.hpp"                //FEM/INT/src
#include "cl_FEM_ElementProxy.hpp"             //FEM/INT/src
#include "cl_FEM_Node_Base.hpp"                //FEM/INT/src
#include "cl_FEM_Element_Factory.hpp"          //FEM/INT/src
#include "cl_FEM_IWG_Factory.hpp"              //FEM/INT/src
#include "cl_FEM_CM_Factory.hpp"               //FEM/INT/src
#include "cl_FEM_SP_Factory.hpp"               //FEM/INT/src
#include "cl_FEM_Set_User_Info.hpp"            //FEM/INT/src

#include "cl_MDL_Model.hpp"

#include "cl_HMR.hpp"
#include "cl_HMR_Background_Mesh.hpp" //HMR/src
#include "cl_HMR_BSpline_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Element.hpp" //HMR/src
#include "cl_HMR_Factory.hpp" //HMR/src
#include "cl_HMR_Field.hpp"
#include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Parameters.hpp" //HMR/src
#include "fn_PRM_HMR_Parameters.hpp" // paramter list for hmr mesh

#include "cl_DLA_Solver_Factory.hpp"
#include "cl_DLA_Solver_Interface.hpp"

#include "cl_NLA_Nonlinear_Solver_Factory.hpp"
#include "cl_NLA_Nonlinear_Solver.hpp"
#include "cl_NLA_Nonlinear_Problem.hpp"
#include "cl_NLA_Nonlinear_Algorithm.hpp"
#include "cl_MSI_Solver_Interface.hpp"
#include "cl_MSI_Equation_Object.hpp"
#include "cl_MSI_Model_Solver_Interface.hpp"
#include "cl_DLA_Linear_Solver_Aztec.hpp"
#include "cl_DLA_Linear_Solver.hpp"
#include "cl_SOL_Warehouse.hpp"

#include "cl_TSA_Time_Solver_Factory.hpp"
#include "cl_TSA_Monolithic_Time_Solver.hpp"
#include "cl_TSA_Time_Solver.hpp"

#include "fn_norm.hpp"

namespace moris
{
namespace mdl
{

void tConstantValueFunction( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
                           moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
                           moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 );
}

void tSelectValueFunction( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
                           moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
                           moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = {{ aParameters( 0 )( 0 ),                      0.0,                     0.0 },
                  {                    0.0,    aParameters( 0 )( 1 ),                     0.0 },
                  {                    0.0,                      0.0,    aParameters( 0 )( 2 )}};
}

void function_TEST_Diffusion_Phase_Change()
{
    /* %%%%%%%%%%%%%%%%%%%%%%%%%% CREATE MESH FOR PROBLEM %%%%%%%%%%%%%%%%%%%%%%%%%% */

    /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SET UP MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

    /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SET UP SOLVER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

    moris::Cell< enum MSI::Dof_Type > tDofTypes( 1 );
    tDofTypes( 0 ) = MSI::Dof_Type::TEMP;

     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     // STEP 1: create linear solver and algorithm
     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

     dla::Solver_Factory  tSolverFactory;
     std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolverFactory.create_solver( sol::SolverType::AMESOS_IMPL );

     //if needed, specify parameters for algorithm

     dla::Linear_Solver tLinearSolver;
     tLinearSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     // STEP 2: create nonlinear solver and algorithm
     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

     NLA::Nonlinear_Solver_Factory tNonlinearFactory;
     std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm =
             tNonlinearFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );

     // specify parameters
     tNonlinearSolverAlgorithm->set_param("NLA_max_iter")   = 10;
     tNonlinearSolverAlgorithm->set_param("NLA_hard_break") = false;
     tNonlinearSolverAlgorithm->set_param("NLA_max_lin_solver_restarts") = 2;
     tNonlinearSolverAlgorithm->set_param("NLA_rebuild_jacobian") = true;

     tNonlinearSolverAlgorithm->set_linear_solver( &tLinearSolver );

     NLA::Nonlinear_Solver tNonlinearSolver;
     tNonlinearSolver.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );

     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     // STEP 3: create time Solver and algorithm
     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     tsa::Time_Solver_Factory tTimeSolverFactory;
     std::shared_ptr< tsa::Time_Solver_Algorithm > tTimeSolverAlgorithm =
             tTimeSolverFactory.create_time_solver( tsa::TimeSolverType::MONOLITHIC );

     tTimeSolverAlgorithm->set_nonlinear_solver( &tNonlinearSolver );

     // create time solver
     tsa::Time_Solver tTimeSolver;

     tTimeSolver.set_time_solver_algorithm( tTimeSolverAlgorithm );

     sol::SOL_Warehouse tSolverWarehouse( tModel->get_solver_interface() );
     tNonlinearSolver.set_solver_warehouse( &tSolverWarehouse );
     tTimeSolver.set_solver_warehouse( &tSolverWarehouse );

     tNonlinearSolver.set_dof_type_list( tDofTypes );
     tTimeSolver.set_dof_type_list( tDofTypes );

     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     // STEP 4: Solve and check
     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

     tTimeSolver.solve();

     moris::Matrix< DDRMat > tFinalSolution;

     tTimeSolver.get_full_solution( tFinalSolution );

} //End main_function

//--------------------------------------------------------------------------------------------------------------------//
TEST_CASE( "DiffusionPhaseChange", "[moris],[mdl],[DiffusionPhaseChange]" )
{
    uint p_size = moris::par_size();

    if( p_size == 1 ) // specify it is a serial test only
    {
        function_TEST_Diffusion_Phase_Change();
    }
}

//--------------------------------------------------------------------------------------------------------------------//
}/* namespace mdl */
}/* namespace moris */

