/*
 * cl_NLA_Arc_Length_Solver_Test.cpp
 *
 *  Created on: Apr 5, 2019
 *      Author: sonne
 */

#include "catch.hpp"
#include "typedefs.hpp"
#include "cl_Communication_Tools.hpp"

//LINALG includes
#include "cl_Matrix.hpp"
#include "fn_equal_to.hpp"
#include "fn_trans.hpp"
#include "linalg_typedefs.hpp"
#include "op_times.hpp"

//DLA includes
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
// sDOF functions
Matrix< DDRMat >
cubic_residual( const moris::sint        aNX,
                const moris::sint        aNY,
                const moris::real        aLambda,
                const Matrix< DDRMat > & aMyValues,
                const moris::uint        aEquationObjectInd )
{
    Matrix< DDRMat > tRes(1,1);
    tRes(0,0) = aLambda*8 - (0.2*std::pow(aMyValues(0,0),3) - 2.1*std::pow(aMyValues(0,0),2) + 6*aMyValues(0,0));
    return tRes;
}

Matrix< DDRMat >
cubic_jacobian( const moris::sint        aNX,
                const moris::sint        aNY,
                const Matrix< DDRMat > & aMyValues,
                const moris::uint        aEquationObjectInd )
{
    Matrix< DDRMat > tJac(1,1);
    tJac(0,0) = 0.6*std::pow(aMyValues(0,0),2) - 4.2*aMyValues(0,0) + 6;
    return tJac;
}

Matrix< DDSMat > test_topo_sDOF( const moris::sint aNX,
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
// mDOF functions
// these are being created/used to verify the solver against a previously completed problems using MATLAB

/*
 * @brief residual function for the mDOF case, this will be used independently on each DOF
 *
 * @param[in] aNX                - dummy
 * @param[in] aNY                - dummy
 * @param[in] aLambda            - current lambda value for iteration
 * @param[in] aMyValues          - current functional values
 * @param[in] aEquationObjectInd - index of current element
 *
 * @param[out] tRes - determined residual value
 */
Matrix< DDRMat >
residual_mDOF( const moris::sint        aNX,
               const moris::sint        aNY,
               const moris::real        aLambda,
               const Matrix< DDRMat > & aMyValues,
               const moris::uint        aEquationObjectInd )
{
    Matrix< DDRMat > tRes;
    // note: aMyValues are the displacement solutions
    //       vector is passed in as a (2x1) so no need to transpose
    /*----------------------------------
     * material and elemental parameters
     * ---------------------------------
     */
    real tSigSat   = 1000000;       // material param
    real tBParam   = 100;           // material param
    real tElLength = 0.5;           // element length
    real tArea     = 0.0003;        // cross-sectional area
    real tJac      = tElLength/2;   // elemental jacobian
    real tWeight   = 2;             // weighting value (assuming single integration point in element)

    Matrix< DDRMat > tBMat(1,2);                                    // define B matrix
    tBMat(0,0) = -1/tElLength;
    tBMat(0,1) =  1/tElLength;

    Matrix< DDRMat > tExternalForce(2,1, 0.0);                           // define external force vector
//    tExternalForce(0,0) = 1.50;
    tExternalForce(0,0) = 0.75;
    tExternalForce(1,0) = 300.75;

    if (aEquationObjectInd==0)      // return the (1x1) residual for the 1st element (only 1-DOF)
    {
        Matrix< DDRMat > tNewVals(2,1, 0.0);
        tNewVals(1,0) = aMyValues(1,0);

        Matrix< DDRMat > tEpsilon = (0.5)*tBMat*tNewVals;                               // determine strain value (this is a scalar), no need to transpose aMyValues here
        real tSigma = tSigSat*(1-std::exp(-tBParam*tEpsilon(0,0)));                     // determine stess value from nonlinear relation
        Matrix< DDRMat > tInternalForce = trans(tBMat)*tSigma*tArea*tJac*tWeight;       // define internal force vector

        tRes.set_size(1,1);
        tRes(0,0) = aLambda*tExternalForce(0,0) - tInternalForce(1,0);
    }
    else if (aEquationObjectInd==1)      // return the (2x1) residual for the 2nd element (2-DOF)
    {
        Matrix< DDRMat > tEpsilon = tBMat*aMyValues;                                    // determine strain value (this is a scalar), no need to transpose aMyValues here
        real tSigma = tSigSat*(1-std::exp(-tBParam*tEpsilon(0,0)));                     // determine stess value from nonlinear relation
        Matrix< DDRMat > tInternalForce = trans(tBMat)*tSigma*tArea*tJac*tWeight;       // define internal force vector

        tRes.set_size(2,1);
        tRes =  aLambda*tExternalForce - tInternalForce;
    }
    else
    {
        MORIS_ASSERT(false,"residual_mDOF(): this test is only setup to have two element");
    }

    return tRes;
}
/*
 * @brief jacobian function for the mDOF case, this will be used independently on each DOF
 *
 * @param[in] aNX             - dummy
 * @param[in] aNY             - dummy
 * @param[in] aMyValues       - current functional values
 * @param[in] aEquationObject - index of current element
 *
 * @param[out] tJac - determined jacobian value
 */
Matrix< DDRMat >
jacobian_mDOF( const moris::sint        aNX,
               const moris::sint        aNY,
               const Matrix< DDRMat > & aMyValues,
               const moris::uint        aEquationObjectInd )
{
    Matrix< DDRMat > tJacobian;
    // note: aMyValues are the displacement solutions
    //       matrix is passed in as a (2x1) so no need to transpose
    /*----------------------------------
     * material and elemental parameters
     * ---------------------------------
     */
    real tSigSat   = 1000000;       // material param
    real tBParam   = 100;           // material param
    real tElLength = 0.5;           // element length
    real tArea     = 0.0003;        // cross-sectional area
    real tJac      = tElLength/2;   // elemental jacobian
    real tWeight   = 2;             // weighting value (using single integration point in element)

    Matrix< DDRMat > tBMat(1,2);                                        // define B matrix
    tBMat(0,0) = -1/tElLength;
    tBMat(0,1) =  1/tElLength;

    if (aEquationObjectInd==0)      // return the (1x1) jacobian for the 1st element (1-DOF)
    {
        Matrix< DDRMat > tNewVals(2,1, 0.0);
        tNewVals(1,0) = aMyValues(1,0);

        Matrix< DDRMat > tEpsilon = tBMat*tNewVals;                        // determine strain value (this is a scalar), no need to transpose aMyValues here
        real tDSigDEps = tBParam*tSigSat*std::exp(-tBParam*tEpsilon(0,0));  // nonlinear material tangent

        Matrix< DDRMat > tTempJac = tDSigDEps*tArea*trans(tBMat)*(tBMat)*tJac*tWeight;
        tJacobian.set_size(1,1);
        tJacobian(0,0) = tTempJac(0,0);
    }
    else if (aEquationObjectInd==1)      // return the (2x2) jacobian for the 2nd element (2-DOF)
    {
        Matrix< DDRMat > tEpsilon = tBMat*aMyValues;                        // determine strain value (this is a scalar), no need to transpose aMyValues here
        real tDSigDEps = tBParam*tSigSat*std::exp(-tBParam*tEpsilon(0,0));  // nonlinear material tangent

        tJacobian.set_size(2,2);
        tJacobian = tDSigDEps*tArea*trans(tBMat)*(tBMat)*tJac*tWeight;
    }
    else
    {
        MORIS_ASSERT(false,"jacobian_mDOF(): this test is only setup to have two element");
    }

    return tJacobian;
}
/*
 * @brief topology function to specify the element to DOF conectivity
 *
 * @param[in] aNX             - dummy
 * @param[in] aNY             - dummy
 * @param[in] aEquationObject - index of current element
 *
 * @param[out] tJac - determined jacobian value
 */
Matrix< DDSMat > test_topo_mDOF( const moris::sint aNX,
                                 const moris::sint aNY,
                                 const moris::uint aEquationObjectInd )
{
    moris::Matrix< moris::DDSMat > tTopo;

    if (aEquationObjectInd==0)
    {
        tTopo.set_size(1,1);
        tTopo(0,0) = 0;
    }
    else if (aEquationObjectInd==1)
    {
        tTopo.set_size(2,1);
        tTopo(0,0) = 0;
        tTopo(1,0) = 1;
    }
    else
    {
        MORIS_ASSERT(false,"test_topo_mdof(): this test is only setup to have two element");
    }

    return tTopo;
}

//------------------------------------------------------------------------------

namespace NLA
{
    TEST_CASE("Arc_Length_Solver_sDOF","[NLA],[NLA_Arc_Length_sDOF]")
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

            Solver_Interface * tSolverInput = new NLA_Solver_Interface_Proxy( 1, 1, 1, 1, cubic_residual, cubic_jacobian, test_topo_sDOF );

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
    } // end arc-length solver test case sDOF

//------------------------------------------------------------------------------
    TEST_CASE("Arc_Length_Solver_mDOF","[NLA],[NLA_Arc_Length_mDOF]")
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

            Solver_Interface * tSolverInput = new NLA_Solver_Interface_Proxy( 2, 2, 1, 1, residual_mDOF, jacobian_mDOF, test_topo_mDOF );

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

            tNonLinSolverAlgorithm->set_param("NLA_max_iter")   = 10;
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
    } // end arc-length solver test case mDOF

//------------------------------------------------------------------------------
} // end NLA namespace
} // end moris namespace
