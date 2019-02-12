/*
 * cl_TSA_Time_Solver_Test.cpp
 *
 *  Created on: Jan 21, 2018
 *      Author: schmidt
 */
#include "catch.hpp"
#include "fn_equal_to.hpp"
#include "typedefs.hpp"
#include "cl_Matrix.hpp"

#include "linalg_typedefs.hpp"
#include "cl_Communication_Tools.hpp"

#include "cl_MSI_Dof_Type_Enums.hpp"

#define protected public
#define private   public
#include "cl_TSA_Staggered_Time_Solver.hpp"
#include "cl_TSA_Monolithic_Time_Solver.hpp"
#include "cl_TSA_Solver_Interface_Proxy2.hpp"
#undef protected
#undef private

namespace moris
{
namespace tsa
{
    TEST_CASE("TimeSolverStaggered","[TSA],[TimeSolverStaggered]")
    {
//        dla::Solver_Factory  tSolFactory;
//        std::shared_ptr< dla::Linear_Solver > tLinSolver1 = tSolFactory.create_solver( SolverType::AZTEC_IMPL );
//        std::shared_ptr< dla::Linear_Solver > tLinSolver2 = tSolFactory.create_solver( SolverType::AZTEC_IMPL );
//
//        tLinSolver1->set_param("AZ_diagnostics") = AZ_none;
//        tLinSolver1->set_param("AZ_output") = AZ_none;
//        tLinSolver2->set_param("AZ_solver") = AZ_gmres;
//        tLinSolver2->set_param("AZ_precond") = AZ_dom_decomp;
//
//        dla::Linear_Solver_Manager * tLinSolManager = new dla::Linear_Solver_Manager();
//
//        tLinSolManager->set_linear_solver( 0, tLinSolver1 );
//        tLinSolManager->set_linear_solver( 1, tLinSolver2 );
//
//        Nonlinear_Solver_Factory tNonlinFactory;
//        std::shared_ptr< Nonlinear_Algorithm > tNonLinSolver = tNonlinFactory.create_nonlinear_solver( NonlinearSolverType::NEWTON_SOLVER );

        Time_Solver * tTimesolver1 = new Staggered_Time_Solver();
        Time_Solver * tTimesolver2 = new Monolithic_Time_Solver();
        Time_Solver * tTimesolver3 = new Monolithic_Time_Solver();

        // Create solver interface
        Solver_Interface * tSolverInput = new TSA_Solver_Interface_Proxy_II();

        // Create solver database
        NLA::SOL_Warehouse tSolverWarehouse( tSolverInput );

        NLA::Nonlinear_Solver tNonlinearSolver1( NLA::NonlinearSolverType::NEWTON_SOLVER );
        NLA::Nonlinear_Solver tNonlinearSolver2( NLA::NonlinearSolverType::NEWTON_SOLVER );

        moris::Cell< enum MSI::Dof_Type > tDofTypes1( 1 );
        tDofTypes1( 0 ) = MSI::Dof_Type::TEMP;

        moris::Cell< enum MSI::Dof_Type > tDofTypes2( 1 );
        tDofTypes2( 0 ) = MSI::Dof_Type::UX;

        tNonlinearSolver1.set_dof_type_list( tDofTypes1 );
        tNonlinearSolver2.set_dof_type_list( tDofTypes2 );

        tNonlinearSolver1.set_solver_warehouse( &tSolverWarehouse );
        tNonlinearSolver2.set_solver_warehouse( &tSolverWarehouse );

        tTimesolver1->set_dof_type_list( tDofTypes1 );
        tTimesolver1->set_dof_type_list( tDofTypes2 );
        tTimesolver2->set_dof_type_list( tDofTypes1 );
        tTimesolver3->set_dof_type_list( tDofTypes2 );

        tTimesolver2->set_nonlinear_algorithm( & tNonlinearSolver1 );
        tTimesolver3->set_nonlinear_algorithm( & tNonlinearSolver2 );

        tTimesolver1->set_time_solver( tTimesolver2 );
        tTimesolver1->set_time_solver( tTimesolver3 );

        tTimesolver1->set_solver_warehouse( &tSolverWarehouse );
        tTimesolver2->set_solver_warehouse( &tSolverWarehouse );
        tTimesolver3->set_solver_warehouse( &tSolverWarehouse );

        tTimesolver1 -> solve();

        Matrix< DDRMat > tSol;
        tTimesolver1 ->get_full_solution( tSol );

        print( tSol, "tSol");
        CHECK( equal_to( tSol( 0, 0 ), -8.869937049794211e-01, 1.0e+08 ) );
        CHECK( equal_to( tSol( 1, 0 ), -8.869937049794211e-01, 1.0e+08 ) );
    }
}
}

