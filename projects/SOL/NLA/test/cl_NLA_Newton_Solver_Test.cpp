/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_NLA_Newton_Solver_Test.cpp
 *
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
#include "cl_SOL_Dist_Vector.hpp"

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
    Matrix< DDRMat > test_residual1(
            const moris::sint        aNX,
            const moris::sint        aNY,
            const moris::real        aLambda,
            const Matrix< DDRMat > & tMyValues,
            const moris::uint        aEquationObjectInd )
    {
        Matrix< DDRMat > tResidual( 2, 1, 0.0);
        tResidual( 0, 0 ) = (0.4 - 10*tMyValues( 0, 0 ) - 0.4*std::pow(tMyValues( 1, 0 ),3) + 5*std::pow(tMyValues( 1, 0 ),2));
        tResidual( 1, 0 ) = (0.15 - 0.4*std::pow(tMyValues( 0, 0 ),3) + 3*std::pow(tMyValues( 0, 0 ),2) - 10*tMyValues( 1, 0 ));

        return tResidual;
    }

    Matrix< DDRMat > test_jacobian1(
            const moris::sint        aNX,
            const moris::sint        aNY,
            const Matrix< DDRMat > & tMyValues,
            const moris::uint        aEquationObjectInd )
    {
        Matrix< DDRMat > tJacobian( 2, 2, 0.0);

        tJacobian( 0, 0 ) = -10;
        tJacobian( 0, 1 ) = -1.2*std::pow(tMyValues( 0, 0 ),2)+6*tMyValues( 0, 0 );
        tJacobian( 1, 0 ) = -1.2*std::pow(tMyValues( 1, 0 ),2)+10*tMyValues( 1, 0 );
        tJacobian( 1, 1 ) = -10;

        return tJacobian;
    }

    Matrix< DDSMat > test_topo1(
            const moris::sint aNX,
            const moris::sint aNY,
            const moris::uint aEquationObjectInd )
    {
        moris::Matrix< moris::DDSMat > tTopo( 2, 1, -1 );

        for ( moris::sint Ik = 0; Ik < 2; Ik++ )
        {
            tTopo( Ik, 0 ) = Ik;
        }

        return tTopo;
    }

    Matrix< DDRMat > test_residual_bratu(
            const moris::sint        aNX,
            const moris::sint        aNY,
            const moris::real        aLambda,
            const Matrix< DDRMat > & tMyValues,
            const moris::uint        aEquationObjectInd )
    {
        moris::sint tEquationObjectOffset = par_rank() * ( aNX * aNX / par_size() );

        moris::sint j = std::floor( ( aEquationObjectInd + tEquationObjectOffset ) / aNX );
        moris::sint i = ( aEquationObjectInd + tEquationObjectOffset ) - ( j * aNX );
        moris::real lambda = 2;
        moris::sint tN = aNX;

        Matrix< DDRMat > tF( 1, 1, 0.0);

        moris::real hx     = 1.0/(tN-1);
        moris::real hy     = 1.0/(tN-1);
        moris::real sc     = hx*hy*lambda;
        moris::real hxdhy  = hx/hy;
        moris::real hydhx  = hy/hx;

        if (i == 0 || j == 0 || i == tN-1 || j == tN-1)
        {
            tF( 0, 0 ) = tMyValues((j*tN) + i, 0 );
        }
        else
        {
            moris::real u    = tMyValues((j*tN) + i, 0 );
            moris::real uxx  = (2.0*u - tMyValues((j*tN) + i-1, 0 ) - tMyValues((j*tN)+i+1, 0 ) )*hydhx;
            moris::real uyy  = (2.0*u - tMyValues(((j-1)*tN)+i, 0 ) - tMyValues(((j+1)*tN)+i, 0 ))*hxdhy;
            tF(0, 0 ) = uxx + uyy - sc*std::exp(u);
        }

        Matrix< DDRMat > tResidual( 5, 1, 0.0);

        tResidual( 2, 0 ) = tF( 0, 0 );

        return tResidual;
    }

    Matrix< DDRMat > test_jacobian_bratu(
            const moris::sint        aNX,
            const moris::sint        aNY,
            const Matrix< DDRMat > & tMyValues,
            const moris::uint        aEquationObjectInd )
    {
        moris::sint tEquationObjectOffset = par_rank() * ( aNX * aNX / par_size() );

        moris::sint j = std::floor( ( aEquationObjectInd + tEquationObjectOffset ) / aNX );
        moris::sint i = ( aEquationObjectInd + tEquationObjectOffset ) - ( j * aNX );
        moris::real lambda = 2;
        moris::sint tN = aNX;

        Matrix< DDRMat > tVal( 1, 5, 0.0);

        moris::real hx     = 1.0/(tN-1);
        moris::real hy     = 1.0/(tN-1);
        moris::real sc     = hx*hy*lambda;
        moris::real hxdhy  = hx/hy;
        moris::real hydhx  = hy/hx;

        if (i == 0 || j == 0 || i == tN-1 || j == tN-1)
        {
            tVal( 0, 0 ) = 0.0;
            tVal( 0, 1 ) = 0.0;
            tVal( 0, 2 ) = 1.0;
            tVal( 0, 3 ) = 0.0;
            tVal( 0, 4 ) = 0.0;
        }
        else
        {
            moris::real u    = tMyValues(( j*tN ) + i, 0 );

            tVal( 0, 0) = - hxdhy;
            tVal( 0, 1 ) = - hydhx;
            tVal( 0, 2 ) = 2.0 * ( hxdhy + hydhx ) - sc*std::exp( u ) ;
            tVal( 0, 3 ) = - hydhx;
            tVal( 0, 4 ) = -hxdhy;
        }

        moris::Matrix< moris::DDRMat > tJacobian( 5, 5, 0.0 );

        for ( moris::sint Ik = 0; Ik < 5; Ik++ )
        {
            tJacobian( 2, Ik ) = tVal( 0, Ik );
        }

        return tJacobian;
    }

    Matrix< DDSMat > test_topo_bratu(
            const moris::sint aNX,
            const moris::sint aNY,
            const moris::uint aEquationObjectInd )
    {

        moris::sint tEquationObjectOffset = par_rank() * ( aNX * aNX / par_size() );

        moris::sint j = std::floor( ( aEquationObjectInd + tEquationObjectOffset ) / aNX );
        moris::sint i = ( aEquationObjectInd + tEquationObjectOffset ) - ( j * aNX );
        moris::sint tN = aNX;

        Matrix< DDSMat > tCol( 1, 5, 0);

        tCol( 0, 0 ) = (((j-1)*tN) + i );
        tCol( 0, 1 ) = ((j*tN) + i -1 );
        tCol( 0, 2 ) = ((j*tN) + i );
        tCol( 0, 3 ) = ((j*tN) + i +1 );
        tCol( 0, 4 ) = (((j+1)*tN) + i );

        moris::Matrix< moris::DDSMat > tTopo( 5, 1, -1 );

        for ( moris::sint Ik = 0; Ik < 5; Ik++ )
        {
            if ( tCol( 0, Ik ) < 0 )
            {
                tCol( 0, Ik ) = (tN*tN);
            }
            else if ( tCol( 0, Ik ) > (tN*tN)-1 )
            {
                tCol( 0, Ik ) = (tN*tN);
            }
            tTopo( Ik, 0 ) = tCol( 0, Ik );
        }

        return tTopo;
    }

    //------------------------------------------------------------------------------

    namespace NLA
    {
        TEST_CASE("Newton Solver Test 1","[NLA],[NLA_Test1]")
        {
            if ( par_size() == 1 )
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
                 * // Inputs are NLA_Solver_Interface_Proxy( Number of Dofs,
                 * //                                        Number of elements,
                 * //                                        Dummy,
                 * //                                        Dummy,
                 * //                                        Residual function pointer,
                 * //                                        Jacobian function pointer,
                 * //                                        Topology function pointer );
                 * Solver_Interface * tSolverInput = new NLA_Solver_Interface_Proxy( 2, 1, 1, 1, test_residual1, test_jacobian1, test_topo1 );
                 * \endcode
                 */
                Solver_Interface * tSolverInput = new NLA_Solver_Interface_Proxy( 2, 1, 1, 1, test_residual1, test_jacobian1, test_topo1 );

                /*!
                 * specify linear solver and create linear solver manager
                 *
                 * \code{.cpp}
                 * dla::Linear_Solver * tLinSolManager = new dla::Linear_Solver();
                 * Nonlinear_Solver  tNonLinSolManager;
                 * \endcode
                 */

                dla::Linear_Solver * tLinSolManager = new dla::Linear_Solver();
                Nonlinear_Solver  tNonLinSolManager;
                tNonLinSolManager.set_solver_interface( tSolverInput );

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
                 * std::shared_ptr< Nonlinear_Algorithm > tNonlLinSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( NonlinearSolverType::NEWTON_SOLVER );
                 *
                 * tNonlLinSolverAlgorithm->set_linear_solver( tLinSolManager );
                 * \endcode
                 */
                Nonlinear_Solver_Factory tNonlinFactory;
                std::shared_ptr< Nonlinear_Algorithm > tNonlLinSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( NonlinearSolverType::NEWTON_SOLVER );

                tNonlLinSolverAlgorithm->set_linear_solver( tLinSolManager );
                /*!
                 * Set nonlinear solver parameters
                 *
                 * \code{.cpp}
                 * tNonlLinSolverAlgorithm->set_param("NLA_max_iter")   = 10;
                 * tNonlLinSolverAlgorithm->set_param("NLA_hard_break") = false;
                 * tNonlLinSolverAlgorithm->set_param("NLA_max_lin_solver_restarts") = 2;
                 *
                 * tNonLinSolManager.set_nonlinear_algorithm( tNonlLinSolverAlgorithm, 0 );
                 * \endcode
                 */
                tNonlLinSolverAlgorithm->set_param("NLA_max_iter")   = 10;
                tNonlLinSolverAlgorithm->set_param("NLA_hard_break") = false;
                tNonlLinSolverAlgorithm->set_param("NLA_max_lin_solver_restarts") = 2;
                tNonlLinSolverAlgorithm->set_param("NLA_rebuild_jacobian") = true;

                tNonLinSolManager.set_nonlinear_algorithm( tNonlLinSolverAlgorithm, 0 );

                /*!
                 * Build linear solver factory and linear solvers.
                 *
                 * \code{.cpp}
                 * dla::Solver_Factory  tSolFactory;
                 * std::shared_ptr< dla::Linear_Solver_Algorithm > tLinSolver1 = tSolFactory.create_solver( sol::SolverType::AZTEC_IMPL );
                 * std::shared_ptr< dla::Linear_Solver_Algorithm > tLinSolver2 = tSolFactory.create_solver( sol::SolverType::AZTEC_IMPL );
                 * \endcode
                 */
                dla::Solver_Factory  tSolFactory;
                std::shared_ptr< dla::Linear_Solver_Algorithm > tLinSolver1 = tSolFactory.create_solver( sol::SolverType::AZTEC_IMPL );
                std::shared_ptr< dla::Linear_Solver_Algorithm > tLinSolver2 = tSolFactory.create_solver( sol::SolverType::AZTEC_IMPL );

                /*!
                 * Set linear solver options
                 *
                 * \code{.cpp}
                 * tLinSolver1->set_param("AZ_diagnostics") = AZ_none;
                 * tLinSolver1->set_param("AZ_output") = AZ_none;
                 * tLinSolver2->set_param("AZ_solver") = AZ_gmres;
                 * tLinSolver2->set_param("AZ_precond") = AZ_dom_decomp;
                 * \endcode
                 */
                tLinSolver1->set_param("AZ_diagnostics") = AZ_none;
                tLinSolver1->set_param("AZ_output") = AZ_none;
                tLinSolver2->set_param("AZ_solver") = AZ_gmres;
                tLinSolver2->set_param("AZ_precond") = AZ_dom_decomp;

                /*!
                 * Set linear solver to linear solver manager
                 *
                 * \code{.cpp}
                 * tLinSolManager->set_linear_algorithm( 0, tLinSolver1 );
                 * tLinSolManager->set_linear_algorithm( 1, tLinSolver2 );
                 * \endcode
                 */
                tLinSolManager->set_linear_algorithm( 0, tLinSolver1 );
                tLinSolManager->set_linear_algorithm( 1, tLinSolver2 );

                /*!
                 * <b> Step 2: Solve nonlinear system </b>
                 */

                /*!
                 * Solve nonlinear system, passing in the nonlinear problem
                 *
                 * \code{.cpp}
                 * tNonLinSolManager.solve( tNonlinearProblem );
                 * \endcode
                 */
                tNonLinSolManager.solve( tNonlinearProblem );

                /*!
                 * Get Solution
                 *
                 * \code{.cpp}
                 * Matrix< DDSMat > tGlobalIndExtract( 2, 1, 0);
                 * tGlobalIndExtract( 1, 0 ) = 1;
                 * Matrix< DDRMat > tMyValues;
                 *
                 * tNonlLinSolverAlgorithm->extract_my_values( 2, tGlobalIndExtract, 0, tMyValues);
                 * \endcode
                 */
                Matrix< DDSMat > tGlobalIndExtract( 2, 1, 0);
                tGlobalIndExtract( 1, 0 ) = 1;
                moris::Cell< Matrix< DDRMat > > tMyValues;

                tNonlLinSolverAlgorithm->extract_my_values( 2, tGlobalIndExtract, 0, tMyValues);

                CHECK( equal_to( tMyValues( 0 )( 0, 0 ), 0.04011965, 1.0e+08 ) );
                CHECK( equal_to( tMyValues( 0 )( 1, 0 ), 0.0154803, 1.0e+08 ) );

                delete( tNonlinearProblem );
                delete( tLinSolManager );
                delete( tSolverInput );
            }
        }

        //------------------------------------------------------------------------------

        TEST_CASE("Newton Solver Test Amesos","[NLA],[NLA_Test_Amesos]")
        {
            if ( par_size() == 1 )
            {
                Solver_Interface * tSolverInput = new NLA_Solver_Interface_Proxy( 2, 1, 1, 1, test_residual1, test_jacobian1, test_topo1 );

                dla::Linear_Solver * tLinSolManager = new dla::Linear_Solver();
                Nonlinear_Solver  tNonLinSolManager;
                tNonLinSolManager.set_solver_interface( tSolverInput );

                Nonlinear_Problem * tNonlinearProblem = new Nonlinear_Problem( tSolverInput );

                Nonlinear_Solver_Factory tNonlinFactory;
                std::shared_ptr< Nonlinear_Algorithm > tNonlLinSolverAlgorithm =
                        tNonlinFactory.create_nonlinear_solver( NonlinearSolverType::NEWTON_SOLVER );

                tNonlLinSolverAlgorithm->set_linear_solver( tLinSolManager );

                tNonlLinSolverAlgorithm->set_param("NLA_max_iter")   = 10;
                tNonlLinSolverAlgorithm->set_param("NLA_hard_break") = false;
                tNonlLinSolverAlgorithm->set_param("NLA_max_lin_solver_restarts") = 2;
                tNonlLinSolverAlgorithm->set_param("NLA_rebuild_jacobian") = true;

                tNonLinSolManager.set_nonlinear_algorithm( tNonlLinSolverAlgorithm, 0 );

                dla::Solver_Factory  tSolFactory;
                std::shared_ptr< dla::Linear_Solver_Algorithm > tLinSolver1 = tSolFactory.create_solver( sol::SolverType::AMESOS_IMPL );
                std::shared_ptr< dla::Linear_Solver_Algorithm > tLinSolver2 = tSolFactory.create_solver( sol::SolverType::AMESOS_IMPL );

                tLinSolManager->set_linear_algorithm( 0, tLinSolver1 );
                tLinSolManager->set_linear_algorithm( 1, tLinSolver2 );

                tNonLinSolManager.solve( tNonlinearProblem );

                Matrix< DDSMat > tGlobalIndExtract( 2, 1, 0);
                tGlobalIndExtract( 1, 0 ) = 1;
                moris::Cell< Matrix< DDRMat > > tMyValues;

                tNonlLinSolverAlgorithm->extract_my_values( 2, tGlobalIndExtract, 0, tMyValues);

                CHECK( equal_to( tMyValues( 0 )( 0, 0 ), 0.04011965, 1.0e+08 ) );
                CHECK( equal_to( tMyValues( 0 )( 1, 0 ), 0.0154803, 1.0e+08 ) );

                delete( tNonlinearProblem );
                delete( tLinSolManager );
                delete( tSolverInput );
            }
        }

        //------------------------------------------------------------------------------

        TEST_CASE("Newton Solver Test Petsc","[NLA],[NLA_Test_Petsc]")
        {
            if ( par_size() == 1 )
            {

                Solver_Interface * tSolverInput = new NLA_Solver_Interface_Proxy( 2, 1, 1, 1, test_residual1, test_jacobian1, test_topo1 );

                dla::Linear_Solver  tLinSolManager;
                Nonlinear_Solver  tNonLinSolManager;
                tNonLinSolManager.set_solver_interface( tSolverInput );

                Nonlinear_Problem tNonlinearProblem( tSolverInput, 0,true, sol::MapType::Petsc );

                Nonlinear_Solver_Factory tNonlinFactory;
                std::shared_ptr< Nonlinear_Algorithm > tNonlLinSolverAlgorithm =
                        tNonlinFactory.create_nonlinear_solver( NonlinearSolverType::NEWTON_SOLVER );

                tNonlLinSolverAlgorithm->set_linear_solver( &tLinSolManager );

                tNonlLinSolverAlgorithm->set_param("NLA_max_iter")   = 10;
                tNonlLinSolverAlgorithm->set_param("NLA_hard_break") = false;
                tNonlLinSolverAlgorithm->set_param("NLA_max_lin_solver_restarts") = 2;
                tNonlLinSolverAlgorithm->set_param("NLA_rebuild_jacobian") = true;

                tNonLinSolManager.set_nonlinear_algorithm( tNonlLinSolverAlgorithm, 0 );

                dla::Solver_Factory  tSolFactory;
                std::shared_ptr< dla::Linear_Solver_Algorithm > tLinSolver1 = tSolFactory.create_solver( sol::SolverType::PETSC );
                std::shared_ptr< dla::Linear_Solver_Algorithm > tLinSolver2 = tSolFactory.create_solver( sol::SolverType::PETSC );

                tLinSolManager.set_linear_algorithm( 0, tLinSolver1 );
                tLinSolManager.set_linear_algorithm( 1, tLinSolver2 );

                tNonLinSolManager.solve( &tNonlinearProblem );

                Matrix< DDSMat > tGlobalIndExtract( 2, 1, 0);
                tGlobalIndExtract( 1, 0 ) = 1;
                moris::Cell< Matrix< DDRMat > > tMyValues;

                tNonlLinSolverAlgorithm->extract_my_values( 2, tGlobalIndExtract, 0, tMyValues);

                CHECK( equal_to( tMyValues( 0 )( 0, 0 ), 0.04011965, 1.0e+08 ) );
                CHECK( equal_to( tMyValues( 0 )( 1, 0 ), 0.0154803, 1.0e+08 ) );

                //        delete( tNonlinearProblem );
                delete( tSolverInput );
            }
        }

        //------------------------------------------------------------------------------

        TEST_CASE("Newton Solver Test InvResidual Relaxation","[NLA],[NLA_Test_InvResidual_Relaxation]")
        {
            if ( par_size() == 1 )
            {
                Solver_Interface * tSolverInput = new NLA_Solver_Interface_Proxy( 2, 1, 1, 1, test_residual1, test_jacobian1, test_topo1 );

                dla::Linear_Solver * tLinSolManager = new dla::Linear_Solver();
                Nonlinear_Solver  tNonLinSolManager;
                tNonLinSolManager.set_solver_interface( tSolverInput );

                Nonlinear_Problem * tNonlinearProblem = new Nonlinear_Problem( tSolverInput );

                Nonlinear_Solver_Factory tNonlinFactory;
                std::shared_ptr< Nonlinear_Algorithm > tNonlLinSolverAlgorithm =
                        tNonlinFactory.create_nonlinear_solver( NonlinearSolverType::NEWTON_SOLVER );

                tNonlLinSolverAlgorithm->set_linear_solver( tLinSolManager );

                tNonlLinSolverAlgorithm->set_param("NLA_max_iter")   = 20;
                tNonlLinSolverAlgorithm->set_param("NLA_hard_break") = false;
                tNonlLinSolverAlgorithm->set_param("NLA_max_lin_solver_restarts") = 2;
                tNonlLinSolverAlgorithm->set_param("NLA_rebuild_jacobian") = true;
                tNonlLinSolverAlgorithm->set_param("NLA_relaxation_strategy") = static_cast< uint >( sol::SolverRelaxationType::InvResNorm );
                tNonlLinSolverAlgorithm->set_param("NLA_relaxation_parameter") = 0.95;

                tNonLinSolManager.set_nonlinear_algorithm( tNonlLinSolverAlgorithm, 0 );

                dla::Solver_Factory  tSolFactory;
                std::shared_ptr< dla::Linear_Solver_Algorithm > tLinSolver1 = tSolFactory.create_solver( sol::SolverType::AMESOS_IMPL );
                std::shared_ptr< dla::Linear_Solver_Algorithm > tLinSolver2 = tSolFactory.create_solver( sol::SolverType::AMESOS_IMPL );

                tLinSolManager->set_linear_algorithm( 0, tLinSolver1 );
                tLinSolManager->set_linear_algorithm( 1, tLinSolver2 );

                tNonLinSolManager.solve( tNonlinearProblem );

                Matrix< DDSMat > tGlobalIndExtract( 2, 1, 0);
                tGlobalIndExtract( 1, 0 ) = 1;
                moris::Cell< Matrix< DDRMat > > tMyValues;

                tNonlLinSolverAlgorithm->extract_my_values( 2, tGlobalIndExtract, 0, tMyValues);

                CHECK( equal_to( tMyValues( 0 )( 0, 0 ), 0.04011965, 1.0e+08 ) );
                CHECK( equal_to( tMyValues( 0 )( 1, 0 ), 0.0154803, 1.0e+08 ) );

                delete( tNonlinearProblem );
                delete( tLinSolManager );
                delete( tSolverInput );
            }
        }

        //------------------------------------------------------------------------------

        TEST_CASE("Newton Solver Test InvResidual Adaptive Relaxation","[NLA],[NLA_Test_InvResidual_Adaptive_Relaxation]")
        {
            if ( par_size() == 1 )
            {
                Solver_Interface * tSolverInput = new NLA_Solver_Interface_Proxy( 2, 1, 1, 1, test_residual1, test_jacobian1, test_topo1 );

                dla::Linear_Solver * tLinSolManager = new dla::Linear_Solver();
                Nonlinear_Solver  tNonLinSolManager;
                tNonLinSolManager.set_solver_interface( tSolverInput );

                Nonlinear_Problem * tNonlinearProblem = new Nonlinear_Problem( tSolverInput );

                Nonlinear_Solver_Factory tNonlinFactory;
                std::shared_ptr< Nonlinear_Algorithm > tNonlLinSolverAlgorithm =
                        tNonlinFactory.create_nonlinear_solver( NonlinearSolverType::NEWTON_SOLVER );

                tNonlLinSolverAlgorithm->set_linear_solver( tLinSolManager );

                tNonlLinSolverAlgorithm->set_param("NLA_max_iter")   = 20;
                tNonlLinSolverAlgorithm->set_param("NLA_hard_break") = false;
                tNonlLinSolverAlgorithm->set_param("NLA_max_lin_solver_restarts") = 2;
                tNonlLinSolverAlgorithm->set_param("NLA_rebuild_jacobian") = true;
                tNonlLinSolverAlgorithm->set_param("NLA_relaxation_strategy") = static_cast< uint >( sol::SolverRelaxationType::InvResNormAdaptive );
                tNonlLinSolverAlgorithm->set_param("NLA_relaxation_parameter") = 1.0;
                tNonlLinSolverAlgorithm->set_param("NLA_relaxation_damping") = 0.95;

                tNonLinSolManager.set_nonlinear_algorithm( tNonlLinSolverAlgorithm, 0 );

                dla::Solver_Factory  tSolFactory;
                std::shared_ptr< dla::Linear_Solver_Algorithm > tLinSolver1 = tSolFactory.create_solver( sol::SolverType::AMESOS_IMPL );
                std::shared_ptr< dla::Linear_Solver_Algorithm > tLinSolver2 = tSolFactory.create_solver( sol::SolverType::AMESOS_IMPL );

                tLinSolManager->set_linear_algorithm( 0, tLinSolver1 );
                tLinSolManager->set_linear_algorithm( 1, tLinSolver2 );

                tNonLinSolManager.solve( tNonlinearProblem );

                Matrix< DDSMat > tGlobalIndExtract( 2, 1, 0);
                tGlobalIndExtract( 1, 0 ) = 1;
                moris::Cell< Matrix< DDRMat > > tMyValues;

                tNonlLinSolverAlgorithm->extract_my_values( 2, tGlobalIndExtract, 0, tMyValues);

                CHECK( equal_to( tMyValues( 0 )( 0, 0 ), 0.04011965, 1.0e+08 ) );
                CHECK( equal_to( tMyValues( 0 )( 1, 0 ), 0.0154803, 1.0e+08 ) );

                delete( tNonlinearProblem );
                delete( tLinSolManager );
                delete( tSolverInput );
            }
        }

        //------------------------------------------------------------------------------

        TEST_CASE("Newton Solver Test 2","[NLA],[NLA_Bratu]")
        {
            if ( par_size() == 4 )
            {
                moris::sint tNumDofsInXandY = 400;
                moris::uint tNumDofs        = (moris::uint)(tNumDofsInXandY*tNumDofsInXandY);
                moris::uint tNumElements    = tNumDofs/par_size();

                /*!
                 * <b> Step 1: Create proxy interface and nonlinear solver </b>
                 */

                /*!
                 * Build linear solver interface. For testing we use the well known, Bratu problem. <br>
                 * Bratu, G. (1914). "Sur les équations intégrales non linéaires". Bulletin de la S. M. F. tome 42: 113–142 <br>
                 * Mohsen, A. (2014). "A simple solution of the Bratu problem". Computer and Mathematics with Applications. 67: 26–33
                 *
                 *  The Bratu problem is used as a solid fuel ignition (SFI) problem. This problem is modeled by the partial differential equation
                 *
                 *  \f[ \Delta u - \lambda e^u = 0 \qquad on \quad \Omega: 0< x,y < 1 \f]
                 *
                 *  with boundary conditions:
                 *
                 *  \f[ u = 0 \quad on \quad \delta \Omega\f]
                 *
                 * \code{.cpp}
                 * // Inputs are NLA_Solver_Interface_Proxy( Number of Dofs,
                 * //                                        Number of elements,
                 * //                                        Number nodes in X,
                 * //                                        Number nodes in Y,
                 * //                                        Residual function pointer,
                 * //                                        Jacobian function pointer,
                 * //                                        Topology function pointer );
                 * Solver_Interface * tSolverInput = new NLA_Solver_Interface_Proxy( tNumDofs, tNumElements, tNumDofsInXandY, tNumDofsInXandY, test_residual_bratu, test_jacobian_bratu, test_topo_bratu );
                 * \endcode
                 */
                Solver_Interface * tSolverInput = new NLA_Solver_Interface_Proxy(
                        tNumDofs,
                        tNumElements,
                        tNumDofsInXandY,
                        tNumDofsInXandY,
                        test_residual_bratu,
                        test_jacobian_bratu,
                        test_topo_bratu );

                /*!
                 * Build linear solver
                 *
                 * \code{.cpp}
                 * dla::Linear_Solver * tLinSolManager = new dla::Linear_Solver();
                 * \endcode
                 */
                dla::Linear_Solver * tLinSolManager = new dla::Linear_Solver();
                Nonlinear_Solver  tNonLinSolManager;
                tNonLinSolManager.set_solver_interface( tSolverInput );

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
                 * std::shared_ptr< Nonlinear_Algorithm > tNonlLinSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( NonlinearSolverType::NEWTON_SOLVER );
                 * \endcode
                 */
                Nonlinear_Solver_Factory tNonlinFactory;
                std::shared_ptr< Nonlinear_Algorithm > tNonlLinSolverAlgorithm =
                        tNonlinFactory.create_nonlinear_solver( NonlinearSolverType::NEWTON_SOLVER );

                /*!
                 * Assign linear solver manager to nonlinear solver
                 *
                 * \code{.cpp}
                 * NonLinSolver->set_linear_solver( tLinSolManager );
                 * \endcode
                 */
                tNonlLinSolverAlgorithm->set_linear_solver( tLinSolManager );
                /*!
                 * Set nonlinear solver parameters
                 *
                 * \code{.cpp}
                 * tNonlLinSolverAlgorithm->set_param("NLA_max_iter")   = 10;
                 * tNonlLinSolverAlgorithm->set_param("NLA_hard_break") = false;
                 * tNonlLinSolverAlgorithm->set_param("NLA_max_lin_solver_restarts") = 2;
                 * tNonlLinSolverAlgorithm->set_param("NLA_rebuild_jacobian") = false;
                 * tNonlLinSolverAlgorithm->set_param("NLA_restart")    = 2;
                 * \endcode
                 */
                tNonlLinSolverAlgorithm->set_param("NLA_max_iter")   = 20;
                tNonlLinSolverAlgorithm->set_param("NLA_hard_break") = false;
                tNonlLinSolverAlgorithm->set_param("NLA_max_lin_solver_restarts") = 2;

                tNonLinSolManager.set_nonlinear_algorithm( tNonlLinSolverAlgorithm, 0 );

                /*!
                 * Build linear solver factory and linear solvers.
                 *
                 * \code{.cpp}
                 * dla::Solver_Factory  tSolFactory;
                 * std::shared_ptr< dla::Linear_Solver_Algorithm > tLinSolver1 = tSolFactory.create_solver( sol::SolverType::AZTEC_IMPL );
                 * std::shared_ptr< dla::Linear_Solver_Algorithm > tLinSolver2 = tSolFactory.create_solver( sol::SolverType::AZTEC_IMPL );
                 * \endcode
                 */
                dla::Solver_Factory  tSolFactory;
                std::shared_ptr< dla::Linear_Solver_Algorithm > tLinSolver1 = tSolFactory.create_solver( sol::SolverType::AMESOS_IMPL );
                std::shared_ptr< dla::Linear_Solver_Algorithm > tLinSolver2 = tSolFactory.create_solver( sol::SolverType::AZTEC_IMPL );

                /*!
                 * Set linear solver options
                 *
                 * \code{.cpp}
                 * tLinSolver2->set_param("AZ_solver") = AZ_gmres;
                 * tLinSolver2->set_param("AZ_precond") = AZ_dom_decomp;
                 * tLinSolver1->set_param("AZ_pre_calc") = AZ_reuse;
                 * tLinSolver1->set_param("AZ_keep_info") = 1;
                 * \endcode
                 */
                //        tLinSolver1->set_param("AZ_diagnostics") = AZ_none;
                //        tLinSolver1->set_param("AZ_output") = AZ_none;
                //        tLinSolver1->set_param("AZ_keep_info") = 1;
                //        //tLinSolver1->set_param("AZ_pre_calc") = AZ_reuse;
                //        tLinSolver1->set_param("AZ_graph_fill") = 5;
                //
                //        tLinSolver1->set_param("ml_prec_type") = "SA";
                //tLinSolver1->set_param("prec_reuse") = true;

                /*!
                 * Set linear solver to linear solver manager
                 *
                 * \code{.cpp}
                 * tLinSolManager->set_linear_algorithm( 0, tLinSolver1 );
                 * tLinSolManager->set_linear_algorithm( 1, tLinSolver2 );
                 * \endcode
                 */
                tLinSolManager->set_linear_algorithm( 0, tLinSolver1 );
                tLinSolManager->set_linear_algorithm( 1, tLinSolver2 );

                /*!
                 * Solve nonlinear system, passing in the nonlinear problem
                 *
                 * \code{.cpp}
                 * tNonLinSolManager.solve( tNonlinearProblem );
                 * \endcode
                 */
                tNonLinSolManager.solve( tNonlinearProblem );

                Matrix< DDSMat > tGlobalIndExtract( 6, 1, 1046);
                tGlobalIndExtract( 1, 0 ) = 13077;
                tGlobalIndExtract( 2, 0 ) = 15089;
                tGlobalIndExtract( 3, 0 ) = 1032;
                tGlobalIndExtract( 4, 0 ) = 777;
                tGlobalIndExtract( 5, 0 ) = 9999;
                moris::Cell< Matrix< DDRMat > > tMyValues;

                tNonlinearProblem->extract_my_values( 6, tGlobalIndExtract, 0, tMyValues);

                if( par_rank() == 0 )
                {
                    CHECK( equal_to( tMyValues( 0 )( 0, 0 ), 0.00354517, 1.0e+08 ) );
                    CHECK( equal_to( tMyValues( 0 )( 1, 0 ), 0.0470188, 1.0e+08 ) );
                    CHECK( equal_to( tMyValues( 0 )( 2, 0 ), 0.05089045, 1.0e+08 ) );
                    CHECK( equal_to( tMyValues( 0 )( 3, 0 ), 0.00361458, 1.0e+08 ) );
                    CHECK( equal_to( tMyValues( 0 )( 4, 0 ), 0.000600491, 1.0e+08 ) );
                    CHECK( equal_to( tMyValues( 0 )( 5, 0 ), 0.0, 1.0e+08 ) );
                }
                if( par_rank() == 2 )
                {
                    CHECK( equal_to( tMyValues( 0 )( 0, 0 ), 0.00354517, 1.0e+08 ) );
                    CHECK( equal_to( tMyValues( 0 )( 1, 0 ), 0.0470188, 1.0e+08 ) );
                    CHECK( equal_to( tMyValues( 0 )( 2, 0 ), 0.05089045, 1.0e+08 ) );
                    CHECK( equal_to( tMyValues( 0 )( 3, 0 ), 0.00361458, 1.0e+08 ) );
                    CHECK( equal_to( tMyValues( 0 )( 4, 0 ), 0.000600491, 1.0e+08 ) );
                    CHECK( equal_to( tMyValues( 0 )( 5, 0 ), 0.0, 1.0e+08 ) );
                }
            }
        }
    }
}

