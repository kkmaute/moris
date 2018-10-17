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

#include "cl_DLA_Solver_Factory.hpp"
#include "cl_DLA_Linear_Solver_Aztec.hpp"

#include "cl_NLA_Nonlinear_Solver_Factory.hpp"

#define protected public
#define private   public
#include "cl_NLA_Newton_Solver.hpp"
#undef protected
#undef private

#include "cl_NLA_Solver_Interface_Proxy.hpp"

namespace moris
{
    Matrix< DDRMat > test_residual(       Matrix< DDRMat > tMyValues,
                                    const moris::uint      aEquationObjectInd )
    {
    moris::real lambda = 2;
    moris::sint tN = 3;

    Matrix< DDRMat > tF( tN*tN*5, 1, 0.0);

    moris::real hx     = 1.0/(tN-1);
    moris::real hy     = 1.0/(tN-1);
    moris::real sc     = hx*hy*lambda;
    moris::real hxdhy  = hx/hy;
    moris::real hydhx  = hy/hx;

    for ( moris::sint j = 0; j < tN; j++ )
    {
        for ( moris::sint i = 0; i < tN; i++ )
        {
            if (i == 0 || j == 0 || i == tN-1 || j == tN-1)
            {
                 tF( (j*tN) + i, 0 ) = tMyValues((j*tN) + i, 0 );
            }
            else
            {
                moris::real u    = tMyValues((j*tN) + i, 0 );
                moris::real uxx  = (2.0*u - tMyValues((j*tN) + i-1, 0 ) - tMyValues((j*tN)+i+1, 0 ) )*hydhx;
                moris::real uyy  = (2.0*u - tMyValues(((j-1)*tN)+i, 0 ) - tMyValues(((j+1)*tN)+i, 0 ))*hxdhy;
                tF((j*tN) + i, 0 ) = uxx + uyy - sc*std::exp(u);
            }
        }
    }
    Matrix< DDRMat > tResidual( 1, 1, 0.0);

    tResidual( 0, 0 ) = tF( aEquationObjectInd, 0);

    return tResidual;
    }


    namespace NLA
    {
    TEST_CASE("Newton Solver Test 1","[NLA],[NLA_Test1]")
    {
        /*!
         * <b> Step 1: Create proxy interface and nonlinear solver </b>
         */

        /*!
         * Build linear solver interface. For testing we use the well known, non-convex "Rosenbrock" benchmark function. <br>
         * Rosenbrock, H.H. (1960). "An automatic method for finding the greatest or least value of a function". The Computer Journal. 3 (3): 175–184
         *
         * The nonlinear problem is stated as followed.
         * \f[ R = \begin{bmatrix} 0.4 -10x_1   -0.4x_2^3+5x_2^2 \\ 0.15-0.4x_1^3+3x_1^2  -10x_2 \\ \end{bmatrix} = 0 \f]
         *
         * Therfore the corresponding linearized system can be stated as followed.
         *
         *    \f[ \begin{bmatrix} 10 & 1.2x_1^2-6X_1 \\ 1.2x_2^2-10x_2 & 10 \\ \end{bmatrix}
         *    \begin{bmatrix} x_1 \\ x_2 \\ \end{bmatrix} =
         *    \begin{bmatrix} 0.4 -10x_1   -0.4x_2^3+5x_2^2 \\ 0.15-0.4x_1^3+3x_1^2  -10x_2 \\ \end{bmatrix} \f]
         *
         * \code{.cpp}
         * Solver_Interface* tSolverInput = new NLA_Solver_Interface_Proxy();
         * \endcode
         */
        Solver_Interface * tSolverInput = new NLA_Solver_Interface_Proxy();

        /*!
         * Create nonlinear problem class
         *
         * \code{.cpp}
         * Nonlinear_Problem * tNonlinearProblem = new Nonlinear_Problem( tSolverInput );
         * \endcode
         */
        Nonlinear_Problem * tNonlinearProblem = new Nonlinear_Problem( tSolverInput );

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
        //std::shared_ptr< Nonlinear_Solver > tNonLinSolver = tNonlinFactory.create_nonlinear_solver( tSolverInput, NonlinearSolverType::NEWTON_SOLVER );
        //tNonLinSolver->set_nonlinear_problem( tNonlinearProblem );

        /*!
         * Set nonlinear solver parameters
         *
         * \code{.cpp}
         * tNonLinSolver->set_param("NLA_max_iter")   = 10;
         * tNonLinSolver->set_param("NLA_hard_break") = false;
         * tNonLinSolver->set_param("NLA_max_lin_solver_restarts") = 2;
         * \endcode
         */
        tNonLinSolver->set_param("NLA_max_iter")   = 10;
        tNonLinSolver->set_param("NLA_hard_break") = false;
        tNonLinSolver->set_param("NLA_max_lin_solver_restarts") = 2;

        /*!
         * Build linear solver factory and linear solvers.
         *
         * \code{.cpp}
         * dla::Solver_Factory  tSolFactory;
         * std::shared_ptr< dla::Linear_Solver > tLinSolver1 = tSolFactory.create_solver( SolverType::AZTEC_IMPL );
         * std::shared_ptr< dla::Linear_Solver > tLinSolver2 = tSolFactory.create_solver( SolverType::AZTEC_IMPL );
         * \endcode
         */
        dla::Solver_Factory  tSolFactory;
        std::shared_ptr< dla::Linear_Solver > tLinSolver1 = tSolFactory.create_solver( SolverType::AZTEC_IMPL );
        std::shared_ptr< dla::Linear_Solver > tLinSolver2 = tSolFactory.create_solver( SolverType::AZTEC_IMPL );

        /*!
         * Set linear solver options
         *
         * \code{.cpp}
         * tLinSolver1->set_param("AZ_solver") = AZ_cg;
         * tLinSolver1->set_param("AZ_precond") = AZ_none;
         * tLinSolver2->set_param("AZ_solver") = AZ_gmres;
         * tLinSolver2->set_param("AZ_precond") = AZ_dom_decomp;
         * \endcode
         */
        tLinSolver1->set_param("AZ_solver") = AZ_cg;
        tLinSolver1->set_param("AZ_precond") = AZ_none;
        tLinSolver1->set_param("AZ_diagnostics") = AZ_none;
        tLinSolver1->set_param("AZ_output") = AZ_none;

        tLinSolver2->set_param("AZ_solver") = AZ_gmres;
        tLinSolver2->set_param("AZ_precond") = AZ_dom_decomp;

        /*!
         * Set linear solver to linear solver manager
         *
         * \code{.cpp}
         * tNonLinSolver->set_linear_solver( 0, tLinSolver1 );
         * tNonLinSolver->set_linear_solver( 1, tLinSolver2 );
         * \endcode
         */
        tNonLinSolver->set_linear_solver( 0, tLinSolver1 );
        tNonLinSolver->set_linear_solver( 1, tLinSolver2 );

        /*!
         * Solve nonlinear system, passing in the nonlinear problem
         *
         * \code{.cpp}
         * tNonLinSolver->solver_nonlinear_system( tNonlinearProblem );
         * \endcode
         */
        tNonLinSolver->solver_nonlinear_system( tNonlinearProblem );

        /*!
         * Get Solution
         *
         * \code{.cpp}
         * tNonLinSolver->extract_my_values( 2, tGlobalIndExtract, 0, tMyValues);
         * \endcode
         */
        Matrix< DDSMat > tGlobalIndExtract( 2, 1, 0);
        tGlobalIndExtract( 1, 0 ) = 1;
        Matrix< DDRMat > tMyValues;

        tNonLinSolver->extract_my_values( 2, tGlobalIndExtract, 0, tMyValues);

        std::cout<<tMyValues(0,0)<<std::endl;
        std::cout<<tMyValues(1,0)<<std::endl;
        CHECK( equal_to( tMyValues( 0, 0 ), 0.04011965, 1.0e+08 ) );
        CHECK( equal_to( tMyValues( 1, 0 ), 0.0154803, 1.0e+08 ) );
    }
}
}

