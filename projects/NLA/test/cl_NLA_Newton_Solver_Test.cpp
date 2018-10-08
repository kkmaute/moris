/*
 * cl_NLA_Newton_Solver.cpp
 *
 *  Created on: Sep 21, 2018
 *      Author: schmidt
 */
#ifdef MORIS_HAVE_PARALLEL
 #include <mpi.h>
#endif

#include "catch.hpp"
#include "fn_equal_to.hpp"
#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Communication_Tools.hpp"

#include "cl_Solver_Factory.hpp"
#include "cl_Linear_Solver_Aztec.hpp"

#include "cl_NLA_Nonlinear_Solver_Factory.hpp"

#define protected public
#define private   public
#include "cl_NLA_Newton_Solver.hpp"
#undef protected
#undef private

#include "cl_NLA_Solver_Input_Test.cpp"

namespace moris
{
    namespace NLA
    {
    TEST_CASE("Newton Solver Test 1","[NLA],[NLA_Test1]")
    {
        /*!
         * <b> Step 1: Preperation: Create nonlinear and linear solver </b>
         */

        /*!
         * Create nonlinear solver factory. Build nonlinear solver.
         *
         * \code{.cpp}
         * Nonlinear_Solver_Factory tNonlinFactory;
         * std::shared_ptr< Nonlinear_Solver > tNonLinSolver = tNonlinFactory.create_nonlinear_solver( NonlinearSolverType::NEWTON_SOLVER );
         * \endcode
         */
        Nonlinear_Solver_Factory tNonlinFactory;
        std::shared_ptr< Nonlinear_Solver > tNonLinSolver = tNonlinFactory.create_nonlinear_solver( NonlinearSolverType::NEWTON_SOLVER );

        /*!
         * Build linear solver interface. Test problem is:
         *
         * |10        1.2*x1^2-6*X1| |x1| = |0.4 -10x1   -0.4x2^3+5x2^2|
         * |1.2*x2^2-10*x2       10| |x2|   |0.15-0.4x1^3+3x1^2  -10x2 |
         *
         * \code{.cpp}
         * Solver_Interface* tSolverInput = new NLA_Solver_Input_Test( tNonLinSolver );
         * \endcode
         */
        Solver_Interface * tSolverInput = new NLA_Solver_Input_Test( tNonLinSolver );

        /*!
         * Create linear solver factory. Build linear solver.
         *
         * \code{.cpp}
         *         Solver_Factory  tSolFactory;
         * std::shared_ptr< Linear_Solver > tLin = tSolFactory.create_solver( tSolverInput, SolverType::AZTEC_IMPL );
         * \endcode
         */
        Solver_Factory  tSolFactory;
        std::shared_ptr< Linear_Solver > tLin = tSolFactory.create_solver( tSolverInput, SolverType::AZTEC_IMPL );

        /*!
         * Assign linear solver to nonlinear solver
         *
         * \code{.cpp}
         * tNonLinSolver->set_linear_solver( tLin );
         * \endcode
         */
        tNonLinSolver->set_linear_solver( tLin );

        /*!
         * Set linear and nonlinear solver parameters.
         *
         * \code{.cpp}
         * tLin->set_param("AZ_diagnostics") = AZ_none;
         * tLin->set_param("AZ_output")      = AZ_none;
         *
         * tNonLinSolver->set_param("NLA_max_iter")   = 10;
         * tNonLinSolver->set_param("NLA_hard_break") = false;
         * \endcode
         */
        tLin->set_param("AZ_diagnostics") = AZ_none;
        tLin->set_param("AZ_output")      = AZ_none;

        tNonLinSolver->set_param("NLA_max_iter")   = 10;
        tNonLinSolver->set_param("NLA_hard_break") = false;

        /*!
         * Solve nonlinear system
         *
         * \code{.cpp}
         * tNonLinSolver->solver_nonlinear_system();
         * \endcode
         */
        tNonLinSolver->solver_nonlinear_system();

        /*!
         * Get Solution
         *
         * \code{.cpp}
         * tNonLinSolver->solver_nonlinear_system();
         * \endcode
         */
        Matrix< DDSMat > tGlobalIndExtract( 2, 1, 0);
        tGlobalIndExtract( 1, 0 ) = 1;
        Matrix< DDRMat > tMyValues;

        tNonLinSolver->extract_my_values( 2, tGlobalIndExtract, 0, tMyValues);

        CHECK( equal_to( tMyValues( 0, 0 ), 0.04011965, 1.0e+08 ) );
        CHECK( equal_to( tMyValues( 1, 0 ), 0.0154803, 1.0e+08 ) );
    }
    }
}

