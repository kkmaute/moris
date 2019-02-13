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

#include "cl_DLA_Solver_Factory.hpp"
#include "cl_NLA_Nonlinear_Solver_Factory.hpp"
#include "cl_TSA_Time_Solver_Factory.hpp"

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
        /*!
         * Create solver factory and linear solver algorithms
         *
         * \code{.cpp}
         * dla::Solver_Factory  tSolFactory;
         * std::shared_ptr< dla::Linear_Solver_Algorithm > tLinSolverAlgorithm_1 = tSolFactory.create_solver( SolverType::AZTEC_IMPL );
         * std::shared_ptr< dla::Linear_Solver_Algorithm > tLinSolverAlgorithm_2 = tSolFactory.create_solver( SolverType::AZTEC_IMPL );
         * \endcode
         */
        dla::Solver_Factory  tSolFactory;
        std::shared_ptr< dla::Linear_Solver_Algorithm > tLinSolverAlgorithm_1 = tSolFactory.create_solver( SolverType::AZTEC_IMPL );
        std::shared_ptr< dla::Linear_Solver_Algorithm > tLinSolverAlgorithm_2 = tSolFactory.create_solver( SolverType::AZTEC_IMPL );

        /*!
         * Set linear solver algorithm parameter settings
         *
         * \code{.cpp}
         * tLinSolverAlgorithm_1->set_param("AZ_diagnostics") = AZ_none;
         * tLinSolverAlgorithm_1->set_param("AZ_output") = AZ_none;
         * tLinSolverAlgorithm_1->set_param("AZ_solver") = AZ_gmres;
         * tLinSolverAlgorithm_1->set_param("AZ_precond") = AZ_dom_decomp;
         * \endcode
         */
        tLinSolverAlgorithm_1->set_param("AZ_diagnostics") = AZ_none;
        tLinSolverAlgorithm_1->set_param("AZ_output") = AZ_none;
        tLinSolverAlgorithm_1->set_param("AZ_solver") = AZ_gmres;
        tLinSolverAlgorithm_1->set_param("AZ_precond") = AZ_dom_decomp;

        /*!
         * Create linear solver
         *
         * \code{.cpp}
         * dla::Linear_Solver * tLinearSolver = new dla::Linear_Solver();
         * \endcode
         */
        dla::Linear_Solver * tLinearSolver = new dla::Linear_Solver();

        /*!
         * Set linear solver algorithms to linear solver
         *
         * \code{.cpp}
         * tLinearSolver->set_linear_algorithm( 0, tLinSolverAlgorithm_1 );
         * tLinearSolver->set_linear_algorithm( 1, tLinSolverAlgorithm_2 );
         * \endcode
         */
        tLinearSolver->set_linear_algorithm( 0, tLinSolverAlgorithm_1 );
        tLinearSolver->set_linear_algorithm( 1, tLinSolverAlgorithm_2 );

        /*!
         * Create nonlinear solver factory and nonlinear solver algorithms
         *
         * \code{.cpp}
         * NLA::Nonlinear_Solver_Factory tNonlinFactory;
         * std::shared_ptr< NLA::Nonlinear_Algorithm > tNonLinSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );
         * \endcode
         */
        NLA::Nonlinear_Solver_Factory tNonlinFactory;
        std::shared_ptr< NLA::Nonlinear_Algorithm > tNonLinSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );

        /*!
         * Set nonlinear solver algorithm parameters
         *
         * \code{.cpp}
         * tNonLinSolverAlgorithm->set_param("NLA_max_iter")   = 10;
         * tNonLinSolverAlgorithm->set_param("NLA_hard_break") = false;
         * tNonLinSolverAlgorithm->set_param("NLA_max_lin_solver_restarts") = 2;
         * tNonLinSolverAlgorithm->set_param("NLA_rebuild_jacobian") = true;
         * \endcode
         */
        tNonLinSolverAlgorithm->set_param("NLA_max_iter")   = 10;
        tNonLinSolverAlgorithm->set_param("NLA_hard_break") = false;
        tNonLinSolverAlgorithm->set_param("NLA_max_lin_solver_restarts") = 2;
        tNonLinSolverAlgorithm->set_param("NLA_rebuild_jacobian") = true;

        /*!
         * Set linear solver to the nonlinear solver algorithm
         *
         * \code{.cpp}
         * tNonLinSolverAlgorithm->set_linear_solver( tLinearSolver );
         * \endcode
         */
        tNonLinSolverAlgorithm->set_linear_solver( tLinearSolver );

        /*!
         * Create nonlinear solver
         *
         * \code{.cpp}
         * NLA::Nonlinear_Solver tNonlinearSolver1( NLA::NonlinearSolverType::NEWTON_SOLVER );
         * NLA::Nonlinear_Solver tNonlinearSolver2( NLA::NonlinearSolverType::NEWTON_SOLVER );
         * \endcode
         */
        NLA::Nonlinear_Solver tNonlinearSolver1;
        NLA::Nonlinear_Solver tNonlinearSolver2;

        /*!
         * Set nonlinear algorithm to nonlinear solver
         *
         * \code{.cpp}
         * tNonlinearSolver1.set_nonlinear_algorithm( tNonLinSolverAlgorithm, 0 );
         * tNonlinearSolver2.set_nonlinear_algorithm( tNonLinSolverAlgorithm, 0 );
         * \endcode
         */
        tNonlinearSolver1.set_nonlinear_algorithm( tNonLinSolverAlgorithm, 0 );
        tNonlinearSolver2.set_nonlinear_algorithm( tNonLinSolverAlgorithm, 0 );

        /*!
         * Create nonlinear solver factory and nonlinear solver algorithms
         *
         * \code{.cpp}
         * Time_Solver_Factory tTimeSolverFactory;
         * std::shared_ptr< dla::Time_Solver > tTimesolver1 = tTimeSolverFactory.create_time_solver( TimeSolverType::STAGGERED );
         * std::shared_ptr< dla::Time_Solver > tTimesolver2 = tTimeSolverFactory.create_time_solver( TimeSolverType::MONOLITHIC );
         * std::shared_ptr< dla::Time_Solver > tTimesolver3 = tTimeSolverFactory.create_time_solver( TimeSolverType::MONOLITHIC );
         * \endcode
         */
        Time_Solver_Factory tTimeSolverFactory;
         Time_Solver * tTimesolver1 = tTimeSolverFactory.create_time_solver( TimeSolverType::STAGGERED );
         Time_Solver * tTimesolver2 = tTimeSolverFactory.create_time_solver( TimeSolverType::MONOLITHIC );
         Time_Solver * tTimesolver3 = tTimeSolverFactory.create_time_solver( TimeSolverType::MONOLITHIC );
//        Time_Solver * tTimesolver1 = new Staggered_Time_Solver();
//        Time_Solver * tTimesolver2 = new Monolithic_Time_Solver();
//        Time_Solver * tTimesolver3 = new Monolithic_Time_Solver();

        /*!
         * Create solver interface
         *
         * \code{.cpp}
         * Solver_Interface * tSolverInput = new TSA_Solver_Interface_Proxy_II();
         * \endcode
         */
        Solver_Interface * tSolverInput = new TSA_Solver_Interface_Proxy_II();

        /*!
         * Create solver warehouse
         *
         * \code{.cpp}
         * NLA::SOL_Warehouse tSolverWarehouse( tSolverInput );
         * \endcode
         */
        NLA::SOL_Warehouse tSolverWarehouse( tSolverInput );

        /*!
         * Create solver warehouse
         *
         * \code{.cpp}
         * moris::Cell< enum MSI::Dof_Type > tDofTypes1( 1 );            tDofTypes1( 0 ) = MSI::Dof_Type::TEMP;
         * moris::Cell< enum MSI::Dof_Type > tDofTypes2( 1 );            tDofTypes2( 0 ) = MSI::Dof_Type::UX;
         * \endcode
         */
        moris::Cell< enum MSI::Dof_Type > tDofTypes1( 1 );            tDofTypes1( 0 ) = MSI::Dof_Type::TEMP;
        moris::Cell< enum MSI::Dof_Type > tDofTypes2( 1 );            tDofTypes2( 0 ) = MSI::Dof_Type::UX;

        /*!
         * Create solver warehouse
         *
         * \code{.cpp}
         * tNonlinearSolver1.set_dof_type_list( tDofTypes1 );
         * tNonlinearSolver2.set_dof_type_list( tDofTypes2 );
         *tTimesolver1->set_dof_type_list( tDofTypes1 );
         *tTimesolver1->set_dof_type_list( tDofTypes2 );
         *tTimesolver2->set_dof_type_list( tDofTypes1 );
         *tTimesolver3->set_dof_type_list( tDofTypes2 );
         * \endcode
         */
        tNonlinearSolver1.set_dof_type_list( tDofTypes1 );
        tNonlinearSolver2.set_dof_type_list( tDofTypes2 );
        tTimesolver1->set_dof_type_list( tDofTypes1 );
        tTimesolver1->set_dof_type_list( tDofTypes2 );
        tTimesolver2->set_dof_type_list( tDofTypes1 );
        tTimesolver3->set_dof_type_list( tDofTypes2 );

        /*!
         * Set time solvers on which a staggered time solver is operating on.
         *
         * \code{.cpp}
         * tTimesolver1->set_time_solver( tTimesolver2 );    ;
         * tTimesolver1->set_time_solver( tTimesolver3 );
         * \endcode
         */
        tTimesolver1->set_time_solver( tTimesolver2 );
        tTimesolver1->set_time_solver( tTimesolver3 );

        /*!
         * Set nonlinear solver to time solver
         *
         * \code{.cpp}
         * tTimesolver2->set_nonlinear_solver( & tNonlinearSolver1 );  ;
         * tTimesolver3->set_nonlinear_solver( & tNonlinearSolver2 );
         * \endcode
         */
        tTimesolver2->set_nonlinear_solver( & tNonlinearSolver1 );
        tTimesolver3->set_nonlinear_solver( & tNonlinearSolver2 );

        /*!
         * Set warehouse to all time and nonliner solvers
         *
         * \code{.cpp}
         * tNonlinearSolver1.set_solver_warehouse( &tSolverWarehouse )
         * tNonlinearSolver2.set_solver_warehouse( &tSolverWarehouse )
         * tTimesolver1->set_solver_warehouse( &tSolverWarehouse );
         * tTimesolver2->set_solver_warehouse( &tSolverWarehouse );
         * tTimesolver3->set_solver_warehouse( &tSolverWarehouse );
         * \endcode
         */
        tNonlinearSolver1.set_solver_warehouse( &tSolverWarehouse );
        tNonlinearSolver2.set_solver_warehouse( &tSolverWarehouse );
        tTimesolver1->set_solver_warehouse( &tSolverWarehouse );
        tTimesolver2->set_solver_warehouse( &tSolverWarehouse );
        tTimesolver3->set_solver_warehouse( &tSolverWarehouse );

        /*!
         * Solve time system
         *
         * \code{.cpp}
         * tTimesolver1 -> solve();
         * \endcode
         */
        tTimesolver1 -> solve();

        Matrix< DDRMat > tSol;
        tTimesolver1 ->get_full_solution( tSol );

        print( tSol, "tSol");
        CHECK( equal_to( tSol( 0, 0 ), -8.869937049794211e-01, 1.0e+08 ) );
        CHECK( equal_to( tSol( 1, 0 ), -8.869937049794211e-01, 1.0e+08 ) );
    }
}
}

