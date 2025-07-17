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

#include "moris_typedefs.hpp"

#include "cl_Communication_Tools.hpp"

#include "cl_DLA_Solver_Interface.hpp"
#include "cl_NLA_Nonlinear_Solver.hpp"

#include "cl_SOL_Dist_Vector.hpp"
#include "cl_SOL_Warehouse.hpp"
#include "cl_SOL_Matrix_Vector_Factory.hpp"

#include "cl_Logger.hpp"
#include "cl_Tracer.hpp"

#include "fn_norm.hpp"

namespace moris::NLA
{
    //--------------------------------------------------------------------------------------------------------------------------

    Solver_Pseudo_Time_Control::Solver_Pseudo_Time_Control(
            Parameter_List&   aParameterListNonlinearSolver,
            sol::Dist_Vector* aCurrentSolution,
            Nonlinear_Solver* aNonLinSolverManager )
    {
        // get relaxation strategy
        mTimeStepStrategy = aParameterListNonlinearSolver.get< sol::SolverPseudoTimeControlType >( "NLA_pseudo_time_control_strategy" );

        // skip setting remaining parameters if no pseudo time step control is used
        if ( mTimeStepStrategy == sol::SolverPseudoTimeControlType::None )
        {
            return;
        }

        // maximum time step size
        mMaxStepSize = aParameterListNonlinearSolver.get< real >( "NLA_pseudo_time_max_step_size" );

        // get maximum number of time steps
        mMaxNumTimeSteps = aParameterListNonlinearSolver.get< sint >( "NLA_pseudo_time_max_num_steps" );

        // get required drop in relative static residual norm
        mRelResNormDrop = aParameterListNonlinearSolver.get< real >( "NLA_pseudo_time_rel_res_norm_drop" );

        // get required relative residual drop for updating "previous" solution and time step1
        mRelResUpdate = aParameterListNonlinearSolver.get< real >( "NLA_pseudo_time_rel_res_norm_update" );

        // get relative static residual norm for switching to steady state computation
        mSteadyStateRelRes = aParameterListNonlinearSolver.get< real >( "NLA_pseudo_time_steady_rel_res_norm" );

        // get time step size used once maximum number of step has been reached
        mSteadyStateStepSize = aParameterListNonlinearSolver.get< real >( "NLA_pseudo_time_steady_step_size" );

        // get time offsets for outputting pseudo time steps; if offset is zero no output is written
        mTimeOffSet = aParameterListNonlinearSolver.get< real >( "NLA_pseudo_time_offset" );

        // number of initial iterations without time step size control at constant time step
        mInitialIterations = aParameterListNonlinearSolver.get< sint >( "NLA_pseudo_time_initial_steps" );

        // iteration index at which reference norm is set
        mRefIterationID = aParameterListNonlinearSolver.get< sint >( "NLA_ref_iter" );

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
            case sol::SolverPseudoTimeControlType::SwitchedRelaxation:
            {
                // get initial time step size
                mConstantStepSize = aParameterListNonlinearSolver.get< real >( "NLA_pseudo_time_constant" );

                // set initial time step size
                mInitialStepSize = mConstantStepSize;

                break;
            }
            case sol::SolverPseudoTimeControlType::ResidualDifference:
            {
                // get initial time step size
                mConstantStepSize = aParameterListNonlinearSolver.get< real >( "NLA_pseudo_time_constant" );

                // get pre-factor
                mResidualFactor = aParameterListNonlinearSolver.get< real >( "NLA_pseudo_time_residual_factor" );

                // get scaling on previous residual
                mResidualExponent = aParameterListNonlinearSolver.get< real >( "NLA_pseudo_time_residual_exponent" );

                // set initial time step size
                mInitialStepSize = mConstantStepSize;

                break;
            }
            case sol::SolverPseudoTimeControlType::Expur:
            {
                // get initial time step size
                mConstantStepSize = aParameterListNonlinearSolver.get< real >( "NLA_pseudo_time_constant" );

                // get pre-factor
                mResidualFactor = aParameterListNonlinearSolver.get< real >( "NLA_pseudo_time_residual_factor" );

                // get scaling on previous residual
                mResidualExponent = aParameterListNonlinearSolver.get< real >( "NLA_pseudo_time_residual_exponent" );

                MORIS_ERROR( mResidualFactor >= 1.0,
                        "Solver_Pseudo_Time_Control::Solver_Pseudo_Time_Control - Expur strategy: NLA_pseudo_time_residual_factor needs to be larger 1.0" );

                MORIS_ERROR( mResidualExponent <= 1.0,
                        "Solver_Pseudo_Time_Control::Solver_Pseudo_Time_Control - Expur strategy: NLA_pseudo_time_residual_exponent needs to be smaller 1.0" );

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

        // overwrite initial time step size if initialization phase is used, i.e. mInitialIterations > 1
        if ( mInitialIterations > 1 )
        {
            mInitialStepSize = aParameterListNonlinearSolver.get< real >( "NLA_pseudo_time_initial" );
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
        sol::Dist_Map* tMapCurrent  = tMatFactory.create_map( mSolverInterface->get_my_local_global_map() );
        sol::Dist_Map* tMapPrevious = tMatFactory.create_map( mSolverInterface->get_my_local_global_map() );

        mFullCurrentSolution  = tMatFactory.create_vector( mSolverInterface, tMapCurrent, tNumRHMS, false, true );
        mFullPreviousSolution = tMatFactory.create_vector( mSolverInterface, tMapPrevious, tNumRHMS, false, true );

        // set flag to compute static residual
        aNonLinSolverManager->set_compute_static_residual_flag( true, true );

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
            real&             aTotalTime,
            real&             aRelResNorm,
            const uint        aIterationId )
    {
        // skip setting remaining parameters if no pseudo time step control is used
        if ( mTimeStepStrategy == sol::SolverPseudoTimeControlType::None )
        {
            return true;
        }

        // skip time control if in initial phase
        // note (a) the internal time step counter is not increased and total time remains at zero
        // note (b) since the new time step is computed for the next iteration the condition is
        //          aIterationId < mInitialIterations instead of aIterationId <= mInitialIterations
        if ( aIterationId < mInitialIterations )
        {
            aTimeStep   = mInitialStepSize;
            aTotalTime  = 0;
            aRelResNorm = 1.0;

            MORIS_LOG_INFO( "In initialization phase - using pseudo time step: %e", aTimeStep );

            // output pseudo time step
            if ( mTimeOffSet > 0.0 )
            {
                // increment pseudo time for output
                mOutputTime += mTimeOffSet;

                // write current solution to output 0
                mSolverInterface->initiate_output( 0, mOutputTime, false );
            }

            // copy current solution onto "previous" solution
            mPreviousSolution->vec_plus_vec( 1.0, *( aCurrentSolution ), 0.0 );

            // return that time continuation has not finished yet
            return false;
        }

        // initialize flag whether to perform update of "previous" solution
        bool tPerformUpdate = false;

        // get static residuals
        const real tRefNorm = aNonLinSolverManager->get_static_ref_norm();
        const real tResNorm = aNonLinSolverManager->get_static_residual_norm();

        // compute relative residual
        aRelResNorm = tResNorm / tRefNorm;

        MORIS_LOG_INFO( "Norm of static reference residual: %e", tRefNorm );
        MORIS_LOG_INFO( "Norm of static residual          : %e", tResNorm );
        MORIS_LOG_INFO( "Ratio of static residual norms   : %e", aRelResNorm );

        // get maximum of relative number of iterations
        const real tRelNumIter = aNonLinSolverManager->get_relative_number_iterations();

        MORIS_LOG_INFO( "Relative number of iterations    : %e", tRelNumIter );

        // compute pseudo time step and return convergence status
        switch ( mTimeStepStrategy )
        {
            // constant relaxation parameter
            case sol::SolverPseudoTimeControlType::Polynomial:
            {
                // update increase time step only if requirement on relative residual is satisfied
                if ( aRelResNorm < mRelResUpdate )
                {
                    // set update flag to true
                    tPerformUpdate = true;

                    // compute new time step
                    aTimeStep = std::max( mConstantStepSize, mTimeStepIndexFactor * std::pow( mTimeStepCounter, mTimeStepExponent ) );

                    // enforce maximum time step size
                    aTimeStep = std::min( aTimeStep, mMaxStepSize );

                    // increase time step counter
                    mTimeStepCounter++;
                }
                break;
            }
            case sol::SolverPseudoTimeControlType::Exponential:
            {
                // update increase time step only if requirement on relative residual is satisfied
                if ( aRelResNorm < mRelResUpdate )
                {
                    // set update flag to true
                    tPerformUpdate = true;

                    // compute new time step
                    aTimeStep = mConstantStepSize * std::pow( mTimeStepIndexFactor, std::floor( mTimeStepCounter / mTimeStepExponent ) * mTimeStepExponent );

                    // enforce maximum time step size
                    aTimeStep = std::min( aTimeStep, mMaxStepSize );

                    // increase time step counter
                    mTimeStepCounter++;
                }
                break;
            }
            case sol::SolverPseudoTimeControlType::InvResNorm:
            {
                // compute new time step
                aTimeStep = mResidualFactor / std::pow( aRelResNorm, mResidualExponent );

                // enforce maximum time step size
                aTimeStep = std::min( aTimeStep, mMaxStepSize );

                // update increase time step only if requirement on relative residual is satisfied
                if ( aRelResNorm < mRelResUpdate )
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

                real tResBased = mResidualFactor / std::pow( aRelResNorm, mResidualExponent );

                aTimeStep = std::max( tIndexBased, tResBased );

                // enforce maximum time step size
                aTimeStep = std::min( aTimeStep, mMaxStepSize );

                // update increase time step only if requirement on relative residual is satisfied
                if ( aRelResNorm < mRelResUpdate )
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
                if ( aRelResNorm < mRelResUpdate )
                {
                    // set update flag to true
                    tPerformUpdate = true;

                    // initialize previous step parameters
                    if ( mTimeStepCounter < 2 )
                    {
                        mPrevStepSize   = mConstantStepSize;
                        mPrevRelResNorm = 1.0;
                    }

                    // compute new time step
                    aTimeStep = mPrevStepSize * std::max( 1.0, mResidualExponent * mPrevRelResNorm / aRelResNorm );

                    // enforce maximum time step size
                    aTimeStep = std::min( aTimeStep, mMaxStepSize );

                    // save previous time step size
                    mPrevStepSize = aTimeStep;

                    // save previous relative residual
                    mPrevRelResNorm = aRelResNorm;

                    // increase time step counter
                    mTimeStepCounter++;
                }
                break;
            }
            case sol::SolverPseudoTimeControlType::ResidualDifference:
            {
                // update increase time step only if requirement on relative residual is satisfied
                if ( aRelResNorm < mRelResUpdate )
                {
                    // set update flag to true
                    tPerformUpdate = true;

                    // initialize previous step parameters
                    if ( mTimeStepCounter < 2 )
                    {
                        mPrevStepSize   = mConstantStepSize;
                        mPrevRelResNorm = 1.0;
                    }

                    // compute exponent
                    real tExponent = std::max( 0.0, ( mResidualExponent * mPrevRelResNorm - aRelResNorm ) / ( mPrevRelResNorm + MORIS_REAL_EPS ) );

                    // compute new time step
                    aTimeStep = mPrevStepSize * std::pow( mResidualFactor, tExponent );

                    // enforce maximum time step size
                    aTimeStep = std::min( aTimeStep, mMaxStepSize );

                    // save previous time step size
                    mPrevStepSize = aTimeStep;

                    // save previous relative residual
                    mPrevRelResNorm = aRelResNorm;

                    // increase time step counter
                    mTimeStepCounter++;
                }
                break;
            }
            case sol::SolverPseudoTimeControlType::Expur:
            {
                // update increase time step only if requirement on relative residual is satisfied
                if ( aRelResNorm < mRelResUpdate )
                {
                    // set update flag to true
                    tPerformUpdate = true;

                    // initialize previous step parameters
                    if ( mTimeStepCounter < 2 )
                    {
                        mPrevStepSize   = mConstantStepSize;
                        mPrevRelResNorm = 1.0;
                    }

                    // initialize new time step
                    aTimeStep = mPrevStepSize;

                    MORIS_LOG_INFO( "tRelNumIter %f", tRelNumIter );

                    // increase time step
                    if ( tRelNumIter < 0.7 )
                    {
                        aTimeStep = mPrevStepSize * mResidualFactor;
                    }

                    // reduce time step
                    if ( tRelNumIter > 0.8 )
                    {
                        aTimeStep = mPrevStepSize * mResidualExponent;
                    }

                    // perform update only if static residual has not increased by more than 25%
                    if ( aRelResNorm > 1.25 * mPrevRelResNorm && aIterationId > mRefIterationID )
                    {
                        tPerformUpdate = false;

                        MORIS_LOG_INFO( "Increase of static residual: %f (RelNumIter = %f)", aRelResNorm / mPrevRelResNorm, tRelNumIter );
                    }

                    // check that time step does not drop below initial one
                    if ( aTimeStep < mConstantStepSize )
                    {
                        // set time step size to initial one
                        aTimeStep = mConstantStepSize;

                        // force update
                        tPerformUpdate = true;
                    }

                    // enforce maximum time step size
                    aTimeStep = std::min( aTimeStep, mMaxStepSize );

                    MORIS_LOG_INFO( "mPrevStepSize %f  aTimeStep %f", mPrevStepSize, aTimeStep );

                    // save previous time step size
                    mPrevStepSize = aTimeStep;

                    // save previous relative residual only if update will be performed
                    if ( tPerformUpdate )
                    {
                        mPrevRelResNorm = aRelResNorm;
                    }

                    // increase time step counter
                    mTimeStepCounter++;
                }
                break;
            }
            case sol::SolverPseudoTimeControlType::Comsol:
            {
                // update increase time step only if requirement on relative residual is satisfied
                if ( aRelResNorm < mRelResUpdate )
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

                    // enforce maximum time step size
                    aTimeStep = std::min( aTimeStep, mMaxStepSize );

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
        if ( aRelResNorm < mSteadyStateRelRes || mSteadyStateMode )
        {
            // set time step to steady state time step
            aTimeStep = mSteadyStateStepSize;

            // activate steady state mode
            mSteadyStateMode = true;

            MORIS_LOG_INFO( "Switched to steady state mode - time step %e", aTimeStep );
        }

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
            real tDeltaNorm = mFullPreviousSolution->vec_norm2()( 0 );

            MORIS_LOG_INFO( "Solution norms: previous = %e  current %e  change %e", tPreviousNorm, tCurrentNorm, tDeltaNorm );

            // copy current solution onto "previous" solution
            mPreviousSolution->vec_plus_vec( 1.0, *( aCurrentSolution ), 0.0 );
        }
        else
        {
            // reset current solution to the one of the previous time step
            aCurrentSolution->vec_plus_vec( 1.0, *( mPreviousSolution ), 0.0 );
        }

        // check if time step size or number of time steps meets convergence criterion
        if ( aRelResNorm < mRelResNormDrop || mTimeStepCounter > mMaxNumTimeSteps )
        {
            MORIS_LOG_INFO( "Pseudo-transient continuation converged: %e < %e || %d > %d",
                    aRelResNorm,
                    mRelResNormDrop,
                    mTimeStepCounter,
                    mMaxNumTimeSteps );

            tIsConverged = true;
        }

        return tIsConverged;
    }
}    // namespace moris::NLA
