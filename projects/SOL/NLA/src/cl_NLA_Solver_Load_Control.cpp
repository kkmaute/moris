/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_NLA_Solver_Load_Control.cpp

 */
#include "cl_NLA_Solver_Load_Control.hpp"

#include "cl_SOL_Dist_Vector.hpp"

#include "cl_NLA_Nonlinear_Solver.hpp"

#include "moris_typedefs.hpp"

#include "cl_Logger.hpp"
#include "cl_Tracer.hpp"

namespace moris::NLA
{
    //--------------------------------------------------------------------------------------------------------------------------

    Solver_Load_Control::Solver_Load_Control( Parameter_List& aParameterListNonlinearSolver )
    {
        // get relaxation strategy
        mLoadControlStrategy = aParameterListNonlinearSolver.get< sol::SolverLoadControlType >( "NLA_load_control_strategy" );

        // get initial relaxation parameter
        mInitialLoadFactor = aParameterListNonlinearSolver.get< real >( "NLA_load_control_factor" );

        // get number of growth steps
        mNumLoadSteps = aParameterListNonlinearSolver.get< sint >( "NLA_load_control_steps" );

        // get required relative residual drop for triggering growth
        mRelativeResidualDropThreshold = aParameterListNonlinearSolver.get< real >( "NLA_load_control_relres" );

        // strategy depending parameters
        switch ( mLoadControlStrategy )
        {
            case sol::SolverLoadControlType::Exponential:
            {
                mExponent = aParameterListNonlinearSolver.get< real >( "NLA_load_control_exponent" );
                break;
            }
            case sol::SolverLoadControlType::UserDefined:
            {
                MORIS_ERROR( false, "Solver_Load_Control::Solver_Load_Control - User Defined Strategy not implemented yet." );
                break;
            }
            default:
            {
            }
        }

        // initialize load step counter
        mLoadStepCounter = 0;
    }

    //--------------------------------------------------------------------------------------------------------------------------

    real
    Solver_Load_Control::get_initial_load_factor()
    {
        return mInitialLoadFactor;
    }

    //--------------------------------------------------------------------------------------------------------------------------

    void
    Solver_Load_Control::eval(
            sint              aIter,
            Nonlinear_Solver* aNonLinSolverManager,
            real&             aLoadFactor )
    {
        // compute relaxation value depending on strategy
        switch ( mLoadControlStrategy )
        {
            // constant relaxation parameter
            case sol::SolverLoadControlType::Constant:
            {
                aLoadFactor = mInitialLoadFactor;

                break;
            }
            case sol::SolverLoadControlType::Linear:
            {
                if ( check_load_step_requirement( aNonLinSolverManager ) )
                {
                    mLoadStepCounter++;
                    aLoadFactor = mInitialLoadFactor + ( 1.0 - mInitialLoadFactor ) * std::min( 1.0, static_cast< real >( mLoadStepCounter ) / static_cast< real >( mNumLoadSteps ) );

                    // log load factor
                    MORIS_LOG_INFO( "Updated load factor (Linear): %7.5e in load step %d", aLoadFactor, mLoadStepCounter );
                }

                break;
            }
            case sol::SolverLoadControlType::Exponential:
            {
                if ( check_load_step_requirement( aNonLinSolverManager ) )
                {
                    mLoadStepCounter++;
                    aLoadFactor = mInitialLoadFactor
                                + ( 1.0 - mInitialLoadFactor ) * std::min( 1.0, std::pow( ( (real)mLoadStepCounter ) / mNumLoadSteps, mExponent ) );

                    // log load factor
                    MORIS_LOG_INFO( "Updated load factor (Exponential): %7.5e in load step %d", aLoadFactor, mLoadStepCounter );
                }

                break;
            }
            case sol::SolverLoadControlType::UserDefined:
            {
                MORIS_ERROR( false, "Solver_Load_Control::Solver_Load_Control - User Defined Strategy not implemented yet." );
                break;
            }
            default:
            {
                MORIS_ERROR( false, "Solver_Load_Control::eval - strategy not implemented.\n" );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------------------

    bool Solver_Load_Control::check_load_step_requirement( Nonlinear_Solver* aNonLinSolverManager )
    {
        real tRefNorm = 1.0;
        real tResNorm = 1.0;

        // use static residual if available
        if ( aNonLinSolverManager->get_use_static_residual_flag() )
        {
            tRefNorm = aNonLinSolverManager->get_static_ref_norm();
            tResNorm = aNonLinSolverManager->get_static_residual_norm();
        }
        else
        {
            tRefNorm = aNonLinSolverManager->get_ref_norm();
            tResNorm = aNonLinSolverManager->get_residual_norm();
        }

        // compute relative residual
        real tRelResNorm = tResNorm / tRefNorm;

        // update load factor if requirement on relative residual is satisfied
        return tRelResNorm < mRelativeResidualDropThreshold;
    }
}    // namespace moris::NLA
