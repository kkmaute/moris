/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_NLA_Arc_Length.cpp
 *
 */

#include "cl_Communication_Tools.hpp"
//------------------------------------------------------------------------------
// nla includes
#include "cl_NLA_Arc_Length.hpp"
#include "cl_NLA_Convergence.hpp"
#include "cl_NLA_Nonlinear_Solver.hpp"
//------------------------------------------------------------------------------
// dla includes
#include "cl_SOL_Matrix_Vector_Factory.hpp"
#include "cl_DLA_Linear_Solver_Algorithm.hpp"
#include "cl_DLA_Solver_Interface.hpp"
#include "cl_SOL_Enums.hpp"
#include "cl_SOL_Dist_Matrix.hpp"
#include "cl_SOL_Dist_Vector.hpp"
//------------------------------------------------------------------------------
// linalg includes
#include "cl_Matrix.hpp"
//------------------------------------------------------------------------------
// std includes
#include <ctime>
#include <cmath>
//------------------------------------------------------------------------------

using namespace moris;
using namespace NLA;
using namespace dla;

//--------------------------------------------------------------------------------------------------------------------------

Arc_Length_Solver::Arc_Length_Solver()
        : Nonlinear_Algorithm()
{
    mLinSolverManager = new dla::Linear_Solver();

    this->set_arc_params();
}

//--------------------------------------------------------------------------------------------------------------------------

Arc_Length_Solver::Arc_Length_Solver( dla::Linear_Solver *aLinSolver )
        : Nonlinear_Algorithm()
{
    mLinSolverManager = aLinSolver;

    this->set_arc_params();
}

//--------------------------------------------------------------------------------------------------------------------------

Arc_Length_Solver::~Arc_Length_Solver()
{
}

//--------------------------------------------------------------------------------------------------------------------------
void Arc_Length_Solver::solver_nonlinear_system( Nonlinear_Problem *aNonlinearProblem )
{
    MORIS_ASSERT( mMyTimeSolverAlgorithm != nullptr, "Arc_Length_Solver::solver_nonlinear_system(): please set the time solver for the arc length algorithm (should be monolithic)" );
    //------------------------------------------------------------------------------

    moris::sint tMaxIts     = mParameterListNonlinearSolver.get< moris::sint >( "NLA_max_iter" );
    moris::real tRelaxation = mParameterListNonlinearSolver.get< moris::real >( "NLA_relaxation_parameter" );

    bool        tIsConverged     = false;
    bool        tHardBreak       = false;
    moris::real tMaxNewTime      = 0.0;
    moris::real tMaxAssemblyTime = 0.0;

    //------------------------------------------------------------------------------
    moris::sint tTimeIter = mMyNonLinSolverManager->get_time_step_iter();
    //------------------------------------------------------------------------------
    /*
     * 1) build jacobian
     * 2) use jacobian to solve for d_tilde
     * 3) store jacobian_0 vals, d_tilde0 vals, and arc function denominator
     * 4) run initialization procedures
     */
    //------------------------------------------------------------------------------

    mNonlinearProblem->get_linearized_problem()->assemble_jacobian();    // reassemble jacobian
    mNonlinearProblem->get_linearized_problem()->assemble_residual();    // reassemble residual

    mJac->get_diagonal( *mJacVal );    // fill vector mJacVal with diagonal values of the Jacobian matrix
    if ( tTimeIter < 2 )
    {    // since mD_tilde is only used in the first initialization step, only compute it while timeStep < 2, note: step count starts at 0
        mGlobalRHS->vec_plus_vec( 1, *mFext, 0 );
        this->solve_linear_system( tTimeIter, tHardBreak );    // inv(Ktilde)*Fext

        mD_tilde->vec_plus_vec( 1, *mNonlinearProblem->get_linearized_problem()->get_full_solver_LHS(), 0 );
        mNonlinearProblem->get_linearized_problem()->assemble_residual();    // rebuild residual since it was changed above
    }

    sint tSize = mDeltaD->vec_local_length();
    if ( tTimeIter == 0 )
    {    // store K_tilde0 diagonal values, d_tilde0 vector, and f_arc denominator ( all constant throughout )
        mJacVal0->vec_plus_vec( 1, *mJacVal, 0 );
        mD_tilde0->vec_plus_vec( 1, *mD_tilde, 0 );
        for ( sint i = 0; i < tSize; i++ )
        {
            mArcDenom += mD_tilde0->get_values_pointer()[ i ] * mJacVal0->get_values_pointer()[ i ] * mD_tilde0->get_values_pointer()[ i ];
        }
    }

    if ( tTimeIter < 2 )
    {
        //------------------------------------------------------------------------------
        // procedure 1
        //------------------------------------------------------------------------------
        mArcNumer = 0.0;
        for ( sint i = 0; i < tSize; i++ )
        {
            mArcNumer += mD_tilde->get_values_pointer()[ i ] * mJacVal0->get_values_pointer()[ i ] * mD_tilde->get_values_pointer()[ i ];
        }
        mF_tilde = std::sqrt( ( 1 - mB ) * ( mArcNumer / mArcDenom ) + mB );

        mDeltaLambda = mDeltaA / mF_tilde;
        mDeltaD->vec_plus_vec( 1, *mD_tilde, 0 );
        mDeltaD->scale_vector( mDeltaLambda );    // DeltaD = DeltaLambda*Dtilde

        mLambdaK = mMyTimeSolverAlgorithm->get_new_lambda() + mDeltaLambda;

        mDK->vec_plus_vec( 1, *mDSolve, 0 );
        mDK->vec_plus_vec( 1, *mDeltaD, 1 );    // D_k = D_n + DeltaD
    }
    else
    {
        //------------------------------------------------------------------------------
        // procedure 2
        //------------------------------------------------------------------------------
        mLambdaK = mLambdaSolveNMinus2 - 3 * mLambdaSolveNMinus1 + 3 * mMyTimeSolverAlgorithm->get_new_lambda();

        mDK->vec_plus_vec( 1, *mDSolveNMinus2, 0 );
        mDK->vec_plus_vec( -3, *mDSolveNMinus1, 1 );
        mDK->vec_plus_vec( 3, *mDSolve, 1 );

        mDeltaLambda = mLambdaK - mMyTimeSolverAlgorithm->get_new_lambda();

        mDeltaD->vec_plus_vec( 1, *mDK, 0 );
        mDeltaD->vec_plus_vec( -1, *mDSolve, 1 );    // DeltaD=D_d-D_n
    }
    // update arc length function, residual, and save initial residual for convergence check
    mArcNumer = 0;
    for ( sint i = 0; i < tSize; i++ )
    {
        mArcNumer += mDeltaD->get_values_pointer()[ i ] * mJacVal0->get_values_pointer()[ i ] * mDeltaD->get_values_pointer()[ i ];
    }
    mFArc = std::sqrt( ( 1 - mB ) * ( mArcNumer / mArcDenom ) + mB * std::pow( mDeltaLambda, 2 ) );

    mR0 = mGlobalRHS->vec_norm2()( 0 );
    mNonlinearProblem->get_linearized_problem()->get_full_solver_LHS()->vec_put_scalar( 0.0 );
    //------------------------------------------------------------------------------
    // Arc Length loop
    //------------------------------------------------------------------------------
    sint tIter = 1;
    while ( ( ( std::abs( mGlobalRHS->vec_norm2()( 0 ) / mR0 ) > mResTol ) || ( ( mFArc - mDeltaA ) / mDeltaA > mForTol ) ) && ( tIter <= tMaxIts ) )
    {
        // pause if needed
        if ( mMyNonLinSolverManager->get_solver_interface()->is_forward_analysis() )
        {
            mForwardPauseFunction();
        }
        else
        {
            mSensitivityPauseFunction();
        }

        clock_t tArcLengthLoopStart = clock();
        clock_t tStartAssemblyTime  = clock();

        // assemble RHS and Jacobian
        mMyTimeSolverAlgorithm->set_lambda_increment( mLambdaK );
        mNonlinearProblem->get_linearized_problem()->assemble_residual();
        mNonlinearProblem->get_linearized_problem()->assemble_jacobian();

        tMaxAssemblyTime = this->calculate_time_needed( tStartAssemblyTime );

        tHardBreak = false;
        //------------------------------------------------------------------------------
        /* Iteration Steps:
         * 1) determine partial derivatives
         * 2) run static condensation method
         * 3) update values and check convergence
         */

        /*
         * (1) Partials
         * --------------------------------------------
         */
        mDFArcDDeltaD->vec_put_scalar( 0 );    // reset the partials
        mDFArcDDeltaLambda = 0;

        for ( sint i = 0; i < tSize; i++ )
        {
            mDFArcDDeltaD->get_values_pointer()[ i ] += mDK->get_values_pointer()[ i ] * mJacVal0->get_values_pointer()[ i ];
        }
        moris::real tScaleVal = ( 1 - mB ) / ( mFArc * mArcDenom );
        mDFArcDDeltaD->scale_vector( tScaleVal );

        mDFArcDDeltaLambda = ( mB * mDeltaLambda ) / mFArc;

        //------------------------------------------------------------------------------
        /*
         * (2) Static Condensation Method
         * --------------------------------------------
         * 2.1) build denominator
         * 2.2) build numerator
         * 2.3) calculate delta_lambda
         * 2.4) use delta_lambda to calculate delta_d
         * --------------------------------------------
         */
        // reset vectors and scalars
        mDelLamNum->vec_put_scalar( 0 );
        mDelLambdaNum = 0;
        mDelLamDen->vec_put_scalar( 0 );
        mDelLambdaDen = 0;
        //-------------------------------------------------
        // solve linear system for denominator
        mGlobalRHS = mNonlinearProblem->get_linearized_problem()->get_solver_RHS();
        mGlobalRHS->vec_plus_vec( 1, *mFext, 0 );

        this->solve_linear_system( tIter, tHardBreak );

        mDelLamDen->vec_plus_vec( 1, *mNonlinearProblem->get_linearized_problem()->get_full_solver_LHS(), 0 );
        //-------------------------------------------------

        // solve linear system for numerator
        mNonlinearProblem->get_linearized_problem()->assemble_residual();
        this->solve_linear_system( tIter, tHardBreak );

        mDelLamNum->vec_plus_vec( 1, *mNonlinearProblem->get_linearized_problem()->get_full_solver_LHS(), 0 );

        //-------------------------------------------------
        for ( sint i = 0; i < tSize; i++ )
        {
            mDelLambdaNum += mDFArcDDeltaD->get_values_pointer()[ i ] * mDelLamNum->get_values_pointer()[ i ];
            mDelLambdaDen += mDFArcDDeltaD->get_values_pointer()[ i ] * mDelLamDen->get_values_pointer()[ i ];
        }

        mdeltaLambda = ( mDeltaA - mFArc - mDelLambdaNum ) / ( mDelLambdaDen + mDFArcDDeltaLambda );

        mdeltaD->vec_plus_vec( 1, *mFext, 0 );
        mdeltaD->scale_vector( mdeltaLambda );

        mGlobalRHS->vec_plus_vec( 1, *mdeltaD, 1 );
        this->solve_linear_system( tIter, tHardBreak );
        mdeltaD->vec_plus_vec( 1, *mNonlinearProblem->get_linearized_problem()->get_full_solver_LHS(), 0 );

        mNonlinearProblem->get_linearized_problem()->assemble_residual();
        //------------------------------------------------------------------------------
        /*
         * (3) Updates
         * --------------------------------------------
         * 3.1) update tD_k          = tD_k + tdelta_d
         * 3.2) update tLambda_k     = tLambda_k + tdelta_lambda
         * 3.3) update tDelta_d      = tDelta_d + tdelta_d
         * 3.4) update tDelta_lambda = tDelta_lambda + tdelta_lambda
         *
         * 3.5) update Fint --- N/A
         * 3.6) update f_arc
         * 3.7) update residual
         * 3.8) check convergence
         */
        //------------------------------------------------------------------------------
        mDK->vec_plus_vec( 1, *mdeltaD, 1 );

        mLambdaK = mLambdaK + mdeltaLambda;

        mDeltaD->vec_plus_vec( 1, *mdeltaD, 1 );

        mDeltaLambda = mDeltaLambda + mdeltaLambda;

        mArcNumer = 0.0;    // reset arc numerator
        for ( sint i = 0; i < tSize; i++ )
        {
            mArcNumer += mDeltaD->get_values_pointer()[ i ] * mJacVal0->get_values_pointer()[ i ] * mDeltaD->get_values_pointer()[ i ];
        }
        mFArc = std::sqrt( ( 1 - mB ) * ( mArcNumer / mArcDenom ) + mB * std::pow( mDeltaLambda, 2 ) );

        mMyTimeSolverAlgorithm->set_lambda_increment( mLambdaK );
        mNonlinearProblem->get_linearized_problem()->assemble_residual();
        //------------------------------------------------------------------------------

        Convergence tConvergence;

        tIsConverged = tConvergence.check_for_convergence( this, tIter, mMyNonLinSolverManager->get_ref_norm(), mMyNonLinSolverManager->get_residual_norm(), tMaxAssemblyTime, tMaxNewTime, tHardBreak );

        if ( tIsConverged )
        {
            if ( tHardBreak )
            {
                continue;
            }
            break;
        }

        // Solve linear system
        this->solve_linear_system( tIter, tHardBreak );

        // PreconTime
        // SolveTime
        ( mNonlinearProblem->get_full_vector() )->vec_plus_vec( -tRelaxation, *mNonlinearProblem->get_linearized_problem()->get_full_solver_LHS(), 1.0 );

        tMaxNewTime = this->calculate_time_needed( tArcLengthLoopStart );

        //            mNonlinearProblem->print_sol_vec( tIter );

        tIter++;
    }    // end iteration loop
    //------------------------------------------------------------------------------
    // store previously converged values
    if ( tTimeIter > 0 )
    {
        mDSolveNMinus2->vec_plus_vec( 1, *mDSolveNMinus1, 0 );
        mLambdaSolveNMinus2 = mLambdaSolveNMinus1;

        mDSolveNMinus1->vec_plus_vec( 1, *mDSolve, 0 );
        mLambdaSolveNMinus1 = mMyTimeSolverAlgorithm->get_new_lambda();
    }
    // update converged values
    mDSolve->vec_plus_vec( 1, *mDK, 0 );
    mMyTimeSolverAlgorithm->set_lambda_increment( mLambdaK );
    std::cout << "=================================================" << std::endl;
    std::cout << "lambda: " << mLambdaK << std::endl;
    std::cout << "=================================================" << std::endl;
    // clear iteration loop variables
    mDK->vec_put_scalar( 0 );
    mdeltaD->vec_put_scalar( 0 );
    mLambdaK     = 0.0;
    mDeltaLambda = 0.0;

    //------------------------------------------------------------------------------
    //------------------------------------------------------------------------------
    // save solutions into vector to print to screen
    // this is temporary and is only being used for debugging purposes
    mDis( tTimeIter, 0 ) = mDSolve->get_values_pointer()[ 1 ];
    mFor( tTimeIter, 0 ) = mFext->get_values_pointer()[ 1 ] * mMyTimeSolverAlgorithm->get_new_lambda();
    // tFor(tTimeIter,0) = tFext->get_values_pointer()[1]*tLambdaSolve;
    //------------------------------------------------------------------------------
    //------------------------------------------------------------------------------
    //  printing solutions
    //  this is temporary and is only being used for debugging purposes
    std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    std::cout << "displacement:   " << std::endl;
    std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    for ( uint i = 0; i < 50; i++ )
    {
        std::cout << mDis( i, 0 ) << std::endl;
    }
    std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    std::cout << "external force:  " << std::endl;
    std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    for ( uint i = 0; i < 50; i++ )
    {
        std::cout << mFor( i, 0 ) << std::endl;
    }
    //------------------------------------------------------------------------------
    //------------------------------------------------------------------------------

}    // end arc-length algorithm

//--------------------------------------------------------------------------------------------------------------------------
void Arc_Length_Solver::solve_linear_system( moris::sint &aIter,
        bool                                             &aHardBreak )
{
    // Solve linear system
    mLinSolverManager->solver_linear_system( mNonlinearProblem->get_linearized_problem(), aIter );
}
//--------------------------------------------------------------------------------------------------------------------------
void Arc_Length_Solver::get_full_solution( moris::Matrix< DDRMat > &LHSValues )
{
    mNonlinearProblem->get_full_vector()->extract_copy( LHSValues );
}

//--------------------------------------------------------------------------------------------------------------------------
void Arc_Length_Solver::get_solution( moris::Matrix< DDRMat > &LHSValues )
{
    mNonlinearProblem->get_full_vector()->extract_copy( LHSValues );
}
//--------------------------------------------------------------------------------------------------------------------------
void Arc_Length_Solver::extract_my_values( const moris::uint &aNumIndices,
        const moris::Matrix< DDSMat >                        &aGlobalBlockRows,
        const moris::uint                                    &aBlockRowOffsets,
        moris::Matrix< DDRMat >                              &LHSValues )
{
    mNonlinearProblem->get_full_vector()->extract_my_values( aNumIndices, aGlobalBlockRows, aBlockRowOffsets, LHSValues );
}
//--------------------------------------------------------------------------------------------------------------------------
void Arc_Length_Solver::set_arc_params( moris::real aBParam,
        moris::real                                 aDeltaAParam,
        moris::real                                 aForTol,
        moris::real                                 aResTol )
{
    mB      = aBParam;
    mDeltaA = aDeltaAParam;
    mForTol = aForTol;
    mResTol = aResTol;
}
//--------------------------------------------------------------------------------------------------------------------------
void Arc_Length_Solver::initialize_variables( Nonlinear_Problem *aNonlinearProblem )
{
    mNonlinearProblem = aNonlinearProblem;

    //------------------------------------------------------------------------------
    // temporary vectors for printing solutions (this is being used only for debugging/validating purposes)
    mDis.set_size( 50, 1, 0.0 );
    mFor.set_size( 50, 1, 0.0 );
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    //    Dist_Vector* tFext = mNonlinearProblem->get_f_ext();      //=============== the external force vector needs to be defined outside the algorithm and passed in
    //    moris::uint tNumDof = 1;
    //    moris::Matrix< DDSMat > tEleConn(1,1);
    //    tEleConn(0,0) = 0;
    //    moris::Matrix< DDRMat > tFextVal(1,1);
    //    tFextVal(0,0) = 8.0;
    //    tFext->sum_into_global_values( tNumDof,tEleConn,tFextVal );

    mFext                           = mNonlinearProblem->get_f_ext();    //=============== the external force vector needs to be defined outside the algorithm and passed in
    moris::uint             tNumDof = 2;
    moris::Matrix< DDSMat > tEleConn( 2, 1 );
    tEleConn( 0, 0 ) = 0;
    tEleConn( 1, 0 ) = 1;
    moris::Matrix< DDRMat > tFextVal( 2, 1, 0.0 );
    tFextVal( 0, 0 ) = 1.5;
    tFextVal( 1, 0 ) = 300.75;
    mFext->sum_into_global_values( tNumDof, tEleConn, tFextVal );
    //------------------------------------------------------------------------------

    mR0 = 1.0;

    mArcNumer = 0;
    mArcDenom = 0;
    mF_tilde  = 0;
    mFArc     = 0;

    mDeltaLambda = 0.0;
    mLambdaK     = 0.0;

    mLambdaSolveNMinus1 = 0;
    mLambdaSolveNMinus2 = 0;

    mDFArcDDeltaLambda = 0;

    mDelLambdaNum = 0;
    mDelLambdaDen = 0;
    mdeltaLambda  = 0;

    //------------------------------------------------------------------------------
    //--------------------get all vectors and matrices------------------------------
    //------------------------------------------------------------------------------
    mJac     = mNonlinearProblem->get_full_for_jacobian();    // full consistent tangent matrix (jacobian)
    mJacVal  = mNonlinearProblem->get_jacobian_diag();        // vector of diagonal values
    mD_tilde = mNonlinearProblem->get_d_tilde();              // displacement vector

    mJacVal0  = mNonlinearProblem->get_jacobian_diag_0();
    mD_tilde0 = mNonlinearProblem->get_d_tilde0();
    mDK       = mNonlinearProblem->get_d_k();

    mDSolve = mNonlinearProblem->get_d_solve();
    mDSolve->vec_put_scalar( 0 );
    mDSolveNMinus1 = mNonlinearProblem->get_d_solve_n_minus_1();
    mDSolveNMinus1->vec_put_scalar( 0 );
    mDSolveNMinus2 = mNonlinearProblem->get_d_solve_n_minus_2();
    mDSolveNMinus2->vec_put_scalar( 0 );

    mGlobalRHS = mNonlinearProblem->get_global_rhs();

    mDFArcDDeltaD = mNonlinearProblem->get_df_dDeltaD();
    mDFArcDDeltaD->vec_put_scalar( 0 );

    mDelLamNum = mNonlinearProblem->get_del_lam_num();
    mDelLamDen = mNonlinearProblem->get_del_lam_den();
    mDeltaD    = mNonlinearProblem->get_del_d_upper();
    mdeltaD    = mNonlinearProblem->get_del_d();
    //------------------------------------------------------------------------------
    bool tRebuildJacobian = true;
    bool tHardBreak       = false;
    sint tDummy           = 1;

    mNonlinearProblem->build_linearized_problem( tRebuildJacobian, tDummy );    // build the linearized problem
    this->solve_linear_system( tDummy, tHardBreak );                            // solve linearized problem

    mGlobalRHS = mNonlinearProblem->get_linearized_problem()->get_solver_RHS();    // set pointer to RHS
    mJac       = mNonlinearProblem->get_linearized_problem()->get_matrix();        // set pointer to jacobian matrix
}
//--------------------------------------------------------------------------------------------------------------------------
void Arc_Length_Solver::set_my_time_solver_algorithm( const std::shared_ptr< tsa::Time_Solver_Algorithm > &aMyTimeSolverAlgorithm )
{
    mMyTimeSolverAlgorithm = aMyTimeSolverAlgorithm;
}
