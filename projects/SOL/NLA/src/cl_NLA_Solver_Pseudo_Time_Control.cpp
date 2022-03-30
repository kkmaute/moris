/*
 * cl_NLA_Solver_Pseudo_Time_Control.cpp
 */

#include "cl_NLA_Solver_Pseudo_Time_Control.hpp"

#include "typedefs.hpp"

#include "cl_Communication_Tools.hpp"

#include "cl_DLA_Solver_Interface.hpp"
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
                ParameterList&      aParameterListNonlinearSolver,
                sol::Dist_Vector*   aCurrentSolution,
                Solver_Interface*   aSolverInterface,
                sol::SOL_Warehouse* aSolverWarehouse )
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
            mSolverInterface = aSolverInterface;

            // create vector to store "previous solution"

            // get num RHS
            uint tNumRHMS = mSolverInterface->get_num_rhs();

            // get full map
            sol::Dist_Map* tFullMap = aCurrentSolution->get_map();

            // create map object
            sol::Matrix_Vector_Factory tMatFactory( aSolverWarehouse->get_tpl_type() );

            // create vector
            mPreviousSolution = tMatFactory.create_vector( mSolverInterface, tFullMap, tNumRHMS );

            // copy current solution onto "previous" solution
            mPreviousSolution->vec_plus_vec( 1.0, *( aCurrentSolution ), 0.0 );

            // save current solution vector set by time solver
            mOldPrevTimeStepVector = mSolverInterface->get_solution_vector_prev_time_step();

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
                const real&       aRefNorm,
                const real&       aResNorm,
                sol::Dist_Vector* aCurrentSolution,
                real&             aTimeStep,
                real&             aTotalTime )
        {
            // initialize flag whether to perform update of "previous" solution
            bool tPerformUpdate = false;

            // compute relative residual
            real tRelResNorm = aResNorm / aRefNorm;

            // compute pseudo time step and return convergence status
            switch ( mTimeStepStrategy )
            {
                // no pseudo time stepping
                case sol::SolverPseudoTimeControlType::None:
                {
                    return true;
                }
                // constant relaxation parameter
                case sol::SolverPseudoTimeControlType::Exponential:
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
                // FIXME: should use vec_norm2 but this does not work in parallel
                Matrix< DDRMat > tVector;
                mPreviousSolution->extract_copy( tVector );
                real tPreviousNorm = sum_all( norm( tVector ) );

                aCurrentSolution->extract_copy( tVector );
                real tCurrentNorm = sum_all( norm( tVector ) );

                // compute change in full vector
                mPreviousSolution->vec_plus_vec( -1.0, *( aCurrentSolution ), 1.0 );
                mPreviousSolution->extract_copy( tVector );
                real tDeltaNorm = sum_all( norm( tVector ) );

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
