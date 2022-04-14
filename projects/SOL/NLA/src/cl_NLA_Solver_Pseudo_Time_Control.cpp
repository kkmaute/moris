/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_NLA_Solver_Pseudo_Time_Control.cpp
 *
 */
#include "cl_NLA_Solver_Pseudo_Time_Control.hpp"

#include "typedefs.hpp"

#include "cl_Communication_Tools.hpp"

#include "cl_DLA_Solver_Interface.hpp"
#include "cl_NLA_Nonlinear_Solver.hpp"

#include "cl_SOL_Dist_Vector.hpp"
#include "cl_SOL_Warehouse.hpp"
#include "cl_SOL_Matrix_Vector_Factory.hpp"

#include "cl_Logger.hpp"
#include "cl_Tracer.hpp"

#include "fn_norm.hpp"

namespace moris
{
    namespace NLA
    {
        //--------------------------------------------------------------------------------------------------------------------------

        Solver_Pseudo_Time_Control::Solver_Pseudo_Time_Control(
                ParameterList&    aParameterListNonlinearSolver,
                sol::Dist_Vector* aCurrentSolution,
                Nonlinear_Solver* aNonLinSolverManager )
        {
            // get relaxation strategy
            mTimeStepStrategy = static_cast< sol::SolverPseudoTimeControlType >(
                    aParameterListNonlinearSolver.get< uint >( "NLA_pseudo_time_control_strategy" ) );

            // skip setting remaining parameters if no pseudo time step control is used
            if ( mTimeStepStrategy == sol::SolverPseudoTimeControlType::None )
            {
                return;
            }

            // get required pseudo time step size needed for convergence
            mRequiredStepSize = aParameterListNonlinearSolver.get< real >( "NLA_pseudo_time_required_step_size" );

            // maximum time step size
            mMaxStepSize = aParameterListNonlinearSolver.get< real >( "NLA_pseudo_time_max_step_size" );

            // get maximum number of time steps
            mMaxNumTimeSteps = aParameterListNonlinearSolver.get< sint >( "NLA_pseudo_time_max_num_steps" );

            // get required relative residual drop for updating "previous" solution and time step
            mRelativeResidualDropThreshold = aParameterListNonlinearSolver.get< real >( "NLA_pseudo_time_relres" );

            // get time offsets for outputting pseudo time steps; if offset is zero no output is written
            mTimeOffSet = aParameterListNonlinearSolver.get< real >( "NLA_pseudo_time_offset" );

            // get time step size used once maximum number of step has been reached
            mSteadyStateStepSize = aParameterListNonlinearSolver.get< real >( "NLA_pseudo_time_steady_state_step_size" );

            // strategy depending parameters
            switch ( mTimeStepStrategy )
            {
                case sol::SolverPseudoTimeControlType::Polynomial:
                {
                    // get constant time step size
                    mConstantStepSize = aParameterListNonlinearSolver.get< real >( "NLA_pseudo_time_constant" );

                    // get pre-factor for time step index-based increase of time step
                    mTimeStepIndexFactor = aParameterListNonlinearSolver.get< real >( "NLA_pseudo_time_step_index_factor" );

                    // get exponent for time step index-based increase
                    mTimeStepExponent = aParameterListNonlinearSolver.get< real >( "NLA_pseudo_time_step_index_exponent" );

                    // set initial time step size
                    mInitialStepSize = mConstantStepSize;

                    break;
                }
                case sol::SolverPseudoTimeControlType::Exponential:
                {
                    // get constant time step size
                    mConstantStepSize = aParameterListNonlinearSolver.get< real >( "NLA_pseudo_time_constant" );

                    // get pre-factor for time step index-based increase of time step
                    mTimeStepIndexFactor = aParameterListNonlinearSolver.get< real >( "NLA_pseudo_time_step_index_factor" );

                    // get exponent for time step index-based increase
                    mTimeStepExponent = aParameterListNonlinearSolver.get< real >( "NLA_pseudo_time_step_index_exponent" );

                    // set initial time step size
                    mInitialStepSize = mConstantStepSize;

                    break;
                }
                case sol::SolverPseudoTimeControlType::InvResNorm:
                {
                    // get pre-factor for residual-based increase of time step
                    mResidualFactor = aParameterListNonlinearSolver.get< real >( "NLA_pseudo_time_residual_factor" );

                    // get exponent for residual-based increase
                    mResidualExponent = aParameterListNonlinearSolver.get< real >( "NLA_pseudo_time_residual_exponent" );

                    // set initial time step size
                    mInitialStepSize = 1.0;

                    break;
                }
                case sol::SolverPseudoTimeControlType::Hybrid:
                {
                    // get constant time step size
                    mConstantStepSize = aParameterListNonlinearSolver.get< real >( "NLA_pseudo_time_constant" );

                    // get pre-factor for time step index-based increase of time step
                    mTimeStepIndexFactor = aParameterListNonlinearSolver.get< real >( "NLA_pseudo_time_step_index_factor" );

                    // get exponent for time step index-based increase
                    mTimeStepExponent = aParameterListNonlinearSolver.get< real >( "NLA_pseudo_time_step_index_exponent" );

                    // get pre-factor for residual-based increase of time step
                    mResidualFactor = aParameterListNonlinearSolver.get< real >( "NLA_pseudo_time_residual_factor" );

                    // get exponent for residual-based increase
                    mResidualExponent = aParameterListNonlinearSolver.get< real >( "NLA_pseudo_time_residual_exponent" );

                    // set initial time step size
                    mInitialStepSize = mConstantStepSize;

                    break;
                }
                case sol::SolverPseudoTimeControlType::Comsol:
                {
                    // get constant time step size
                    mComsolParameter_1 = aParameterListNonlinearSolver.get< real >( "NLA_pseudo_time_comsol_1" );
                    mComsolParameter_2 = aParameterListNonlinearSolver.get< real >( "NLA_pseudo_time_comsol_2" );

                    // set initial time step size
                    mInitialStepSize = 1.0;

                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "Solver_Pseudo_Time_Control::Solver_Pseudo_Time_Control - Strategy not implemented yet." );
                }
            }

            // store solver interface
            mSolverInterface = aNonLinSolverManager->get_solver_interface();

            // create vector to store "previous solution"

            // get num RHS
            uint tNumRHMS = mSolverInterface->get_num_rhs();

            // get local map
            sol::Dist_Map* tLocalMap = aCurrentSolution->get_map();

            // create map object
            sol::Matrix_Vector_Factory tMatFactory( aNonLinSolverManager->get_solver_warehouse()->get_tpl_type() );

            // create vector
            mPreviousSolution = tMatFactory.create_vector( mSolverInterface, tLocalMap, tNumRHMS );

            // copy current solution onto "previous" solution
            mPreviousSolution->vec_plus_vec( 1.0, *( aCurrentSolution ), 0.0 );

            // save current solution vector set by time solver
            mOldPrevTimeStepVector = mSolverInterface->get_solution_vector_prev_time_step();

            // create distributed vectors to compute change in solution
            sol::Dist_Map* tFullMapCurrent  = tMatFactory.create_map( mSolverInterface->get_my_local_global_map() );
            sol::Dist_Map* tFullMapPrevious = tMatFactory.create_map( mSolverInterface->get_my_local_global_map() );

            mFullCurrentSolution  = tMatFactory.create_vector( mSolverInterface, tFullMapCurrent, tNumRHMS );
            mFullPreviousSolution = tMatFactory.create_vector( mSolverInterface, tFullMapPrevious, tNumRHMS );

            // set flag to compute static residual
            aNonLinSolverManager->set_compute_static_residual_flag( true );

            // save time frames set by time solver
            mTimeFrameCurrent  = mSolverInterface->get_time();
            mTimeFramePrevious = mSolverInterface->get_previous_time();

            // reset true previous solution
            mSolverInterface->set_solution_vector_prev_time_step( mPreviousSolution );

            // set pseudo time frames
            Matrix< DDRMat > tCurrent  = { { 1 }, { 2 } };
            Matrix< DDRMat > tPrevious = { { 0 }, { 1 } };

            mSolverInterface->set_time( tCurrent );
            mSolverInterface->set_previous_time( tPrevious );
        }

        //--------------------------------------------------------------------------------------------------------------------------

        Solver_Pseudo_Time_Control::~Solver_Pseudo_Time_Control()
        {
            // check if pseudo time stepping is used
            if ( mTimeStepStrategy == sol::SolverPseudoTimeControlType::None )
            {
                return;
            }

            // delete vector for "previous solution"
            delete mPreviousSolution;

            // delete vectors used to compute change in solution norm
            delete mFullCurrentSolution;
            delete mFullPreviousSolution;

            // reset true previous solution
            mSolverInterface->set_solution_vector_prev_time_step( mOldPrevTimeStepVector );

            // reset true time frames
            mSolverInterface->set_time( mTimeFrameCurrent );
            mSolverInterface->set_previous_time( mTimeFramePrevious );
        }

        //--------------------------------------------------------------------------------------------------------------------------

        bool
        Solver_Pseudo_Time_Control::get_initial_step_size( real& aInitialStepSize )
        {
            // set initial time step
            aInitialStepSize = mInitialStepSize;

            // check if pseudo time stepping is used - if not return true,i.e. converged
            if ( mTimeStepStrategy == sol::SolverPseudoTimeControlType::None )
            {
                return true;
            }

            // return not converged
            return false;
        }

        //--------------------------------------------------------------------------------------------------------------------------

        bool
        Solver_Pseudo_Time_Control::compute_time_step_size(
                Nonlinear_Solver* aNonLinSolverManager,
                sol::Dist_Vector* aCurrentSolution,
                real&             aTimeStep,
                real&             aTotalTime )
        {

            // skip setting remaining parameters if no pseudo time step control is used
            if ( mTimeStepStrategy == sol::SolverPseudoTimeControlType::None )
            {
                return true;
            }

            // initialize flag whether to perform update of "previous" solution
            bool tPerformUpdate = false;

            // get static residuals
            const real aRefNorm = aNonLinSolverManager->get_static_ref_norm();
            const real aResNorm = aNonLinSolverManager->get_static_residual_norm();

            // compute relative residual
            const real tRelResNorm = aResNorm / aRefNorm;

            MORIS_LOG_INFO( "Norm of static reference residual: %e", aRefNorm );
            MORIS_LOG_INFO( "Norm of static residual:           %e", aResNorm );
            MORIS_LOG_INFO( "Ratio of static residual norms:    %e", tRelResNorm );

            // compute pseudo time step and return convergence status
            switch ( mTimeStepStrategy )
            {
                // constant relaxation parameter
                case sol::SolverPseudoTimeControlType::Polynomial:
                {
                    // update increase time step only if requirement on relative residual is satisfied
                    if ( tRelResNorm < mRelativeResidualDropThreshold )
                    {
                        // set update flag to true
                        tPerformUpdate = true;

                        // compute new time step
                        aTimeStep = std::max( mConstantStepSize, mTimeStepIndexFactor * std::pow( mTimeStepCounter, mTimeStepExponent ) );

                        // increase time step counter
                        mTimeStepCounter++;
                    }
                    break;
                }
                case sol::SolverPseudoTimeControlType::Exponential:
                {
                    // update increase time step only if requirement on relative residual is satisfied
                    if ( tRelResNorm < mRelativeResidualDropThreshold )
                    {
                        // set update flag to true
                        tPerformUpdate = true;

                        // compute new time step
                        aTimeStep = mConstantStepSize * std::pow( mTimeStepIndexFactor, std::floor( mTimeStepCounter / mTimeStepExponent ) * mTimeStepExponent );

                        // increase time step counter
                        mTimeStepCounter++;
                    }
                    break;
                }
                case sol::SolverPseudoTimeControlType::InvResNorm:
                {
                    // compute new time step
                    aTimeStep = mResidualFactor / std::pow( tRelResNorm, mResidualExponent );

                    // update increase time step only if requirement on relative residual is satisfied
                    if ( tRelResNorm < mRelativeResidualDropThreshold )
                    {
                        // increase time step counter
                        mTimeStepCounter++;

                        // set update flag to true
                        tPerformUpdate = true;
                    }
                    break;
                }
                case sol::SolverPseudoTimeControlType::Hybrid:
                {
                    // compute new time step
                    real tIndexBased = std::max( mConstantStepSize, mTimeStepIndexFactor * std::pow( mTimeStepCounter, mTimeStepExponent ) );

                    real tResBased = mResidualFactor / std::pow( tRelResNorm, mResidualExponent );

                    aTimeStep = std::max( tIndexBased, tResBased );

                    // update increase time step only if requirement on relative residual is satisfied
                    if ( tRelResNorm < mRelativeResidualDropThreshold )
                    {
                        // increase time step counter
                        mTimeStepCounter++;

                        // set update flag to true
                        tPerformUpdate = true;
                    }
                    break;
                }
                case sol::SolverPseudoTimeControlType::SwitchedRelaxation:
                {
                    // update increase time step only if requirement on relative residual is satisfied
                    if ( tRelResNorm < mRelativeResidualDropThreshold )
                    {
                        // set update flag to true
                        tPerformUpdate = true;

                        if ( mTimeStepCounter < 2 )
                        {
                            aTimeStep = mConstantStepSize;
                        }
                        else
                        {
                            aTimeStep = mPrevStepSize * mPrevRelResNorm / tRelResNorm;
                        }

                        // save previous time step size
                        mPrevStepSize = aTimeStep;

                        // save previous relative residual
                        mPrevRelResNorm = tRelResNorm;

                        // increase time step counter
                        mTimeStepCounter++;
                    }
                    break;
                }
                case sol::SolverPseudoTimeControlType::Comsol:
                {
                    // update increase time step only if requirement on relative residual is satisfied
                    if ( tRelResNorm < mRelativeResidualDropThreshold )
                    {
                        // set update flag to true
                        tPerformUpdate = true;

                        // compute new time step
                        aTimeStep = std::pow( 1.3, std::min( (real)mTimeStepCounter, 9.0 ) );

                        if ( mTimeStepCounter > mComsolParameter_1 )
                        {
                            aTimeStep += 9.0 * std::pow( 1.3, std::min( mTimeStepCounter - mComsolParameter_1, 9.0 ) );
                        }

                        if ( mTimeStepCounter > mComsolParameter_2 )
                        {
                            aTimeStep += 90.0 * std::pow( 1.3, std::min( mTimeStepCounter - mComsolParameter_2, 9.0 ) );
                        }

                        // increase time step counter
                        mTimeStepCounter++;
                    }
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "Solver_Pseudo_Time_Control::eval - strategy not implemented.\n" );
                }
            }

            // set steady state step size
            if ( mTimeStepCounter > mMaxNumTimeSteps && mSteadyStateStepSize > 0.0 )
            {
                aTimeStep = mSteadyStateStepSize;
            }

            // enforce maximum time step size
            aTimeStep = std::min( aTimeStep, mMaxStepSize );

            // update total time
            mTotalTime += aTimeStep;
            aTotalTime = mTotalTime;

            // output pseudo time step
            if ( mTimeOffSet > 0.0 )
            {
                // increment pseudo time for output
                mOutputTime += mTimeOffSet;

                // write current solution to output 0
                mSolverInterface->initiate_output( 0, mOutputTime, false );
            }

            // initialize convergence flag
            bool tIsConverged = false;

            // update "previous" solution only if requirement on relative residual is satisfied
            if ( tPerformUpdate )
            {
                // log load factor
                MORIS_LOG_INFO( "Updated pseudo time step - updated previous time step in time step %d", mTimeStepCounter );

                // compute norms for previous and current solutions
                mFullPreviousSolution->import_local_to_global( *mPreviousSolution );
                real tPreviousNorm = mFullPreviousSolution->vec_norm2()( 0 );

                mFullCurrentSolution->import_local_to_global( *aCurrentSolution );
                real tCurrentNorm = mFullCurrentSolution->vec_norm2()( 0 );

                mFullPreviousSolution->vec_plus_vec( -1.0, *( mFullCurrentSolution ), 1.0 );
                real tDeltaNorm = mFullCurrentSolution->vec_norm2()( 0 );

                MORIS_LOG_INFO( "Solution norms: previous = %e  current %e  change %e", tPreviousNorm, tCurrentNorm, tDeltaNorm );

                // copy current solution onto "previous" solution
                mPreviousSolution->vec_plus_vec( 1.0, *( aCurrentSolution ), 0.0 );

                // check if time step size or number of time steps meets convergence criterion
                if ( aTimeStep > mRequiredStepSize || mTimeStepCounter > mMaxNumTimeSteps )
                {
                    tIsConverged = true;
                }
            }

            return tIsConverged;
        }
    }    // namespace NLA
}    // namespace moris
