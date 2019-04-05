/*
 * cl_NLA_Arc_Length_Solver_Test.cpp
 *
 *  Created on: Apr 5, 2019
 *      Author: sonne
 */

#include "catch.hpp"
#include "fn_equal_to.hpp"
#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Communication_Tools.hpp"

#include "cl_DLA_Solver_Factory.hpp"
#include "cl_DLA_Linear_Solver_Aztec.hpp"
#include "cl_DLA_Linear_Solver.hpp"
#include "cl_Vector.hpp"

#define protected public
#define private   public
#include "cl_NLA_Nonlinear_Solver.hpp"
#include "cl_NLA_Nonlinear_Solver_Factory.hpp"
#include "cl_NLA_Newton_Solver.hpp"
#include "cl_NLA_Nonlinear_Problem.hpp"
#include "cl_NLA_Solver_Interface_Proxy.hpp"
#undef protected
#undef private

namespace moris
{
//------------------------------------------------------------------------------
Matrix< DDRMat >
cubic_residual( const moris::sint        aNX,
                const moris::sint        aNY,
                const moris::real        aLambda,
                const Matrix< DDRMat > & tMyValues,
                const moris::uint        aEquationObjectInd )
{
    Matrix< DDRMat > tRes(1,1);
    tRes(0,0) = aLambda*8 - (0.2*std::pow(tMyValues(0,0),3) - 2.1*std::pow(tMyValues(0,0),2) + 6*tMyValues(0,0));
    return tRes;
}

Matrix< DDRMat >
cubic_jacobian( const moris::sint        aNX,
                const moris::sint        aNY,
                const Matrix< DDRMat > & tMyValues,
                const moris::uint        aEquationObjectInd )
{
    Matrix< DDRMat > tJac(1,1);
    tJac(0,0) = 0.6*std::pow(tMyValues(0,0),2) - 4.2*tMyValues(0,0) + 6;
    return tJac;
}

Matrix< DDSMat > test_topo2( const moris::sint aNX,
                             const moris::sint aNY,
                             const moris::uint aEquationObjectInd )
{
moris::Matrix< moris::DDSMat > tTopo( 1, 1, -1 );

for ( moris::sint Ik = 0; Ik < 1; Ik++ )
{
    tTopo( Ik, 0 ) = Ik;
}

return tTopo;
}
//------------------------------------------------------------------------------

namespace NLA
{
    TEST_CASE("Arc_Length Solver","[NLA],[NLA_Arc_Length]")
    {
        if ( par_size() == 1 )
        {
            // Inputs are NLA_Solver_Interface_Proxy( Number of Dofs,
            //                                        Number of elements,
            //                                        Dummy,
            //                                        Dummy,
            //                                        Residual function pointer,
            //                                        Jacobian function pointer,
            //                                        Topology function pointer );

            Solver_Interface * tSolverInput = new NLA_Solver_Interface_Proxy( 1, 1, 1, 1, cubic_residual, cubic_jacobian, test_topo2 );

            // specify the linear solver and create the linear solver manager
            dla::Linear_Solver * tLinSolManager = new dla::Linear_Solver();
            Nonlinear_Solver  tNonLinSolManager;


            // create nonlinear problem class
            Nonlinear_Problem * tNonlinearProblem = new Nonlinear_Problem( tSolverInput );

            // create nonlinear solver factory. Build nonlinear solver.
            Nonlinear_Solver_Factory tNonlinFactory;
            std::shared_ptr< Nonlinear_Algorithm > tNonLinSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( NonlinearSolverType::ARC_LENGTH_SOLVER );

            // set nonlinear solver parameters
            tNonLinSolverAlgorithm->set_linear_solver( tLinSolManager );

            tNonLinSolverAlgorithm->set_param("NLA_max_iter")   = 5;
            tNonLinSolverAlgorithm->set_param("NLA_hard_break") = false;
            tNonLinSolverAlgorithm->set_param("NLA_max_lin_solver_restarts") = 2;

            tNonLinSolManager.set_nonlinear_algorithm( tNonLinSolverAlgorithm, 0 );

            // build linear solver factory and linear solvers
            dla::Solver_Factory  tSolFactory;
            std::shared_ptr< dla::Linear_Solver_Algorithm > tLinSolver1 = tSolFactory.create_solver( SolverType::AZTEC_IMPL );
            std::shared_ptr< dla::Linear_Solver_Algorithm > tLinSolver2 = tSolFactory.create_solver( SolverType::AZTEC_IMPL );

            // set linear solver options
            tLinSolver1->set_param("AZ_diagnostics") = AZ_none;
            tLinSolver1->set_param("AZ_output") = AZ_none;
            tLinSolver2->set_param("AZ_solver") = AZ_gmres;
            tLinSolver2->set_param("AZ_precond") = AZ_dom_decomp;

            // set linear solver to linear solver manager
            tLinSolManager->set_linear_algorithm( 0, tLinSolver1 );
            tLinSolManager->set_linear_algorithm( 1, tLinSolver2 );

            // solve nonlinear system, passing in the nonlinear problem
            tNonLinSolManager.solve( tNonlinearProblem );

            // get solution
            Matrix< DDSMat > tGlobalIndExtract( 2, 1, 0);
            tGlobalIndExtract( 0, 0 ) = 1;
            Matrix< DDRMat > tMyValues;

            tNonLinSolverAlgorithm->get_full_solution( tMyValues );
            print(tMyValues, "tMyValues");
        }
    } // end arc-length solver test case

} // end NLA namespace
} // end moris namespace
