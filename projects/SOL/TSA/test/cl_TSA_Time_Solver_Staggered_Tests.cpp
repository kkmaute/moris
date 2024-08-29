/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_TSA_Time_Solver_Staggered_Tests.cpp
 *
 */

#include "catch.hpp"
#include "fn_equal_to.hpp"
#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"

#include "linalg_typedefs.hpp"
#include "cl_Communication_Tools.hpp"

#include "cl_DLA_Solver_Factory.hpp"
#include "cl_NLA_Nonlinear_Solver_Factory.hpp"
#include "cl_NLA_Nonlinear_Algorithm.hpp"
#include "cl_TSA_Time_Solver_Factory.hpp"

#include "cl_MSI_Dof_Type_Enums.hpp"

#define protected public
#define private   public
#include "cl_TSA_Staggered_Time_Solver.hpp"
#include "cl_TSA_Monolithic_Time_Solver.hpp"
#include "cl_TSA_Solver_Interface_Proxy2.hpp"
#include "cl_TSA_Time_Solver.hpp"
#undef protected
#undef private

#include "cl_DLA_Linear_Solver.hpp"
#include "cl_SOL_Warehouse.hpp"
#include "cl_NLA_Nonlinear_Solver.hpp"

namespace moris::tsa
{
    TEST_CASE("TimeSolverStaggered","[TSA],[TimeSolverStaggered]")
    {

//        if ( par_size() == 1 )
//        {
//        /*!
//         * Create solver warehouse
//         *
//         * \code{.cpp}
//         * Vector< enum MSI::Dof_Type > tDofTypes1( 1 );            tDofTypes1( 0 ) = MSI::Dof_Type::TEMP;
//         * Vector< enum MSI::Dof_Type > tDofTypes2( 1 );            tDofTypes2( 0 ) = MSI::Dof_Type::UX;
//         * \endcode
//         */
//        Vector< enum MSI::Dof_Type > tDofTypes1( 1 );            tDofTypes1( 0 ) = MSI::Dof_Type::TEMP;
//        Vector< enum MSI::Dof_Type > tDofTypes2( 1 );            tDofTypes2( 0 ) = MSI::Dof_Type::UX;
//
//        /*!
//         * Create solver factory and linear solver algorithms
//         *
//         * \code{.cpp}
//         * dla::Solver_Factory  tSolFactory;
//         * std::shared_ptr< dla::Linear_Solver_Algorithm > tLinSolverAlgorithm_1 = tSolFactory.create_solver( sol::SolverType::AZTEC_IMPL );
//         * std::shared_ptr< dla::Linear_Solver_Algorithm > tLinSolverAlgorithm_2 = tSolFactory.create_solver( sol::SolverType::AZTEC_IMPL );
//         * \endcode
//         */
//        dla::Solver_Factory  tSolFactory;
//        std::shared_ptr< dla::Linear_Solver_Algorithm > tLinSolverAlgorithm_1 = tSolFactory.create_solver( sol::SolverType::AZTEC_IMPL );
//        std::shared_ptr< dla::Linear_Solver_Algorithm > tLinSolverAlgorithm_2 = tSolFactory.create_solver( sol::SolverType::AZTEC_IMPL );
//
//        /*!
//         * Set linear solver algorithm parameter settings
//         *
//         * \code{.cpp}
//         * tLinSolverAlgorithm_1->set_param("AZ_diagnostics") = AZ_none;
//         * tLinSolverAlgorithm_1->set_param("AZ_output") = AZ_none;
//         * tLinSolverAlgorithm_1->set_param("AZ_solver") = AZ_gmres;
//         * tLinSolverAlgorithm_1->set_param("AZ_precond") = AZ_dom_decomp;
//         * \endcode
//         */
//        tLinSolverAlgorithm_1->set_param("AZ_diagnostics") = AZ_none;
//        tLinSolverAlgorithm_1->set_param("AZ_output") = AZ_none;
//        tLinSolverAlgorithm_1->set_param("AZ_solver") = AZ_gmres;
//        tLinSolverAlgorithm_1->set_param("AZ_precond") = AZ_dom_decomp;
//
//        /*!
//         * Create linear solver
//         *
//         * \code{.cpp}
//         * dla::Linear_Solver * tLinearSolver = new dla::Linear_Solver();
//         * \endcode
//         */
//        dla::Linear_Solver * tLinearSolver = new dla::Linear_Solver();
//
//        /*!
//         * Set linear solver algorithms to linear solver
//         *
//         * \code{.cpp}
//         * tLinearSolver->set_linear_algorithm( 0, tLinSolverAlgorithm_1 );
//         * tLinearSolver->set_linear_algorithm( 1, tLinSolverAlgorithm_2 );
//         * \endcode
//         */
//        tLinearSolver->set_linear_algorithm( 0, tLinSolverAlgorithm_1 );
//        tLinearSolver->set_linear_algorithm( 1, tLinSolverAlgorithm_2 );
//
//        /*!
//         * Create nonlinear solver factory and nonlinear solver algorithms
//         *
//         * \code{.cpp}
//         * NLA::Nonlinear_Solver_Factory tNonlinFactory;
//         * std::shared_ptr< NLA::Nonlinear_Algorithm > tNonLinSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );
//         * \endcode
//         */
//        NLA::Nonlinear_Solver_Factory tNonlinFactory;
//        std::shared_ptr< NLA::Nonlinear_Algorithm > tNonLinSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );
//
//        /*!
//         * Set nonlinear solver algorithm parameters
//         *
//         * \code{.cpp}
//         * tNonLinSolverAlgorithm->set_param("NLA_max_iter")   = 10;
//         * tNonLinSolverAlgorithm->set_param("NLA_hard_break") = false;
//         * tNonLinSolverAlgorithm->set_param("NLA_max_lin_solver_restarts") = 2;
//         * tNonLinSolverAlgorithm->set_param("NLA_rebuild_jacobian") = true;
//         * \endcode
//         */
//        tNonLinSolverAlgorithm->set_param("NLA_max_iter")   = 10;
//        tNonLinSolverAlgorithm->set_param("NLA_hard_break") = false;
//        tNonLinSolverAlgorithm->set_param("NLA_max_lin_solver_restarts") = 2;
//        tNonLinSolverAlgorithm->set_param("NLA_rebuild_jacobian") = true;
//
//        /*!
//         * Set linear solver to the nonlinear solver algorithm
//         *
//         * \code{.cpp}
//         * tNonLinSolverAlgorithm->set_linear_solver( tLinearSolver );
//         * \endcode
//         */
//        tNonLinSolverAlgorithm->set_linear_solver( tLinearSolver );
//
//        /*!
//         * Create nonlinear solver
//         *
//         * \code{.cpp}
//         * NLA::Nonlinear_Solver tNonlinearSolver1( NLA::NonlinearSolverType::NEWTON_SOLVER );
//         * NLA::Nonlinear_Solver tNonlinearSolver2( NLA::NonlinearSolverType::NEWTON_SOLVER );
//         * \endcode
//         */
//        NLA::Nonlinear_Solver tNonlinearSolver1;
//        NLA::Nonlinear_Solver tNonlinearSolver2;
//
//        /*!
//         * Set nonlinear algorithm to nonlinear solver
//         *
//         * \code{.cpp}
//         * tNonlinearSolver1.set_nonlinear_algorithm( tNonLinSolverAlgorithm, 0 );
//         * tNonlinearSolver2.set_nonlinear_algorithm( tNonLinSolverAlgorithm, 0 );
//         * \endcode
//         */
//        tNonlinearSolver1.set_nonlinear_algorithm( tNonLinSolverAlgorithm, 0 );
//        tNonlinearSolver2.set_nonlinear_algorithm( tNonLinSolverAlgorithm, 0 );
//
//        /*!
//         * Create nonlinear solver factory and nonlinear solver algorithms
//         *
//         * \code{.cpp}
//         * Time_Solver_Factory tTimeSolverFactory;
//         * std::shared_ptr< dla::Time_Solver_Algorithm > tTimesolverAlgorithm_1 = tTimeSolverFactory.create_time_solver( TimeSolverType::STAGGERED );
//         * std::shared_ptr< dla::Time_Solver_Algorithm > tTimesolverAlgorithm_2 = tTimeSolverFactory.create_time_solver( TimeSolverType::MONOLITHIC );
//         * std::shared_ptr< dla::Time_Solver_Algorithm > tTimesolverAlgorithm_3 = tTimeSolverFactory.create_time_solver( TimeSolverType::MONOLITHIC );
//         * \endcode
//         */
//        Time_Solver_Factory tTimeSolverFactory;
//        std::shared_ptr< Time_Solver_Algorithm > tTimeSolverAlgorithm_1 = tTimeSolverFactory.create_time_solver( TimeSolverType::STAGGERED );
//        std::shared_ptr< Time_Solver_Algorithm > tTimeSolverAlgorithm_2 = tTimeSolverFactory.create_time_solver( TimeSolverType::MONOLITHIC );
//        std::shared_ptr< Time_Solver_Algorithm > tTimeSolverAlgorithm_3 = tTimeSolverFactory.create_time_solver( TimeSolverType::MONOLITHIC );
//
//        /*!
//         * Set nonlinear solver to time solver algorithm
//         *
//         * \code{.cpp}
//         * tTimesolver2->set_nonlinear_solver( & tNonlinearSolver1 );  ;
//         * tTimesolver3->set_nonlinear_solver( & tNonlinearSolver2 );
//         * \endcode
//         */
//        tTimeSolverAlgorithm_2->set_nonlinear_solver( & tNonlinearSolver1 );
//        tTimeSolverAlgorithm_3->set_nonlinear_solver( & tNonlinearSolver2 );
//
//        tTimeSolverAlgorithm_2->set_param("TSA_Num_Time_Steps")   = 1000;
//        tTimeSolverAlgorithm_2->set_param("TSA_Time_Frame")       = 10.0;
//        tTimeSolverAlgorithm_3->set_param("TSA_Num_Time_Steps")   = 1000;
//        tTimeSolverAlgorithm_3->set_param("TSA_Time_Frame")       = 10.0;
//
//        /*!
//         * Create time solver
//         *
//         * \code{.cpp}
//         * Time_Solver tTimeSolver_Main;
//         * Time_Solver tTimeSolver_2;
//         * Time_Solver tTimeSolver_3;
//         * \endcode
//         */
//        Time_Solver tTimeSolver_Main;
//        Time_Solver tTimeSolver_2;
//        Time_Solver tTimeSolver_3;
//
//        /*!
//         * Set time solver algorithms to time solvers
//         *
//         * \code{.cpp}
//         * tTimeSolver_Main.set_time_solver_algorithm( tTimeSolverAlgorithm_1 );
//         * tTimeSolver_2.set_time_solver_algorithm( tTimeSolverAlgorithm_2 );
//         * tTimeSolver_3.set_time_solver_algorithm( tTimeSolverAlgorithm_3 );
//         * \endcode
//         */
//        tTimeSolver_Main.set_time_solver_algorithm( tTimeSolverAlgorithm_1 );
//        tTimeSolver_2.set_time_solver_algorithm( tTimeSolverAlgorithm_2 );
//        tTimeSolver_3.set_time_solver_algorithm( tTimeSolverAlgorithm_3 );
//
//        /*!
//         * Set time solvers on which a staggered time solver is operating on.
//         *
//         * \code{.cpp}
//         * tTimeSolver_Main.set_sub_time_solver( &tTimeSolver_2 );
//         * tTimeSolver_Main.set_sub_time_solver( &tTimeSolver_3 );
//         * \endcode
//         */
//        tTimeSolver_Main.set_sub_time_solver( &tTimeSolver_2 );
//        tTimeSolver_Main.set_sub_time_solver( &tTimeSolver_3 );
//
//        /*!
//         * Create solver interface
//         *
//         * \code{.cpp}
//         * Solver_Interface * tSolverInput = new TSA_Solver_Interface_Proxy_II();
//         * \endcode
//         */
//        Solver_Interface * tSolverInput = new TSA_Solver_Interface_Proxy_II();
//
//        /*!
//         * Create solver warehouse
//         *
//         * \code{.cpp}
//         * sol::SOL_Warehouse tSolverWarehouse( tSolverInput );
//         * \endcode
//         */
//        sol::SOL_Warehouse tSolverWarehouse( tSolverInput );
//
//        /*!
//         * Create solver warehouse
//         *
//         * \code{.cpp}
//         * tNonlinearSolver1.set_dof_type_list( tDofTypes1 );
//         * tNonlinearSolver2.set_dof_type_list( tDofTypes2 );
//         * tTimeSolver_Main.set_dof_type_list( tDofTypes1 );
//         * tTimeSolver_Main.set_dof_type_list( tDofTypes2 );
//         * tTimeSolver_2.set_dof_type_list( tDofTypes1 );
//         * tTimeSolver_3.set_dof_type_list( tDofTypes2 );
//         * \endcode
//         */
//        tNonlinearSolver1.set_dof_type_list( tDofTypes1 );
//        tNonlinearSolver2.set_dof_type_list( tDofTypes2 );
//        tTimeSolver_Main.set_dof_type_list( tDofTypes1 );
//        tTimeSolver_Main.set_dof_type_list( tDofTypes2 );
//        tTimeSolver_2.set_dof_type_list( tDofTypes1 );
//        tTimeSolver_3.set_dof_type_list( tDofTypes2 );
//
//        /*!
//         * Set warehouse to all time and nonlinear solvers
//         *
//         * \code{.cpp}
//         * tNonlinearSolver1.set_solver_warehouse( &tSolverWarehouse );
//         * tNonlinearSolver2.set_solver_warehouse( &tSolverWarehouse );
//         * tTimeSolver_Main.set_solver_warehouse( &tSolverWarehouse );
//         * tTimeSolver_2.set_solver_warehouse( &tSolverWarehouse );
//         * tTimeSolver_3.set_solver_warehouse( &tSolverWarehouse );
//         * \endcode
//         */
//        tNonlinearSolver1.set_solver_warehouse( &tSolverWarehouse );
//        tNonlinearSolver2.set_solver_warehouse( &tSolverWarehouse );
//        tTimeSolver_Main.set_solver_warehouse( &tSolverWarehouse );
//        tTimeSolver_2.set_solver_warehouse( &tSolverWarehouse );
//        tTimeSolver_3.set_solver_warehouse( &tSolverWarehouse );
//
//        /*!
//         * Solve time system
//         *
//         * \code{.cpp}
//         * tTimesolver1 -> solve();
//         * \endcode
//         */
//        tTimeSolver_Main.solve();
//
//        /*!
//         * Get_solution
//         *
//         * \code{.cpp}
//         * Matrix< DDRMat > tSol;
//         * tTimeSolver_Main.get_full_solution( tSol );
//         * \endcode
//         */
//        Matrix< DDRMat > tSol;
//        tTimeSolver_Main.get_full_solution( tSol );
//
//        print( tSol, "tSol");
//        CHECK( equal_to( tSol( 0, 0 ), -8.869937049794211e-01, 1.0e+08 ) );
//        CHECK( equal_to( tSol( 1, 0 ), -8.869937049794211e-01, 1.0e+08 ) );
//        }
    }
}
