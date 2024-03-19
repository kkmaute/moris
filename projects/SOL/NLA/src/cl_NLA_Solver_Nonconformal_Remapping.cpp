/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_NLA_Solver_Nonconformal_Remapping.cpp
 *
 */
#include "cl_NLA_Nonlinear_Solver.hpp"
#include "cl_NLA_Solver_Nonconformal_Remapping.hpp"

#include <math.h>
#include "cl_Param_List.hpp"
#include <cstdlib>

namespace moris::NLA
{
    Solver_Nonconformal_Remapping::Solver_Nonconformal_Remapping( ParameterList &aParameterListNonlinearSolver )
            : mStrategy( static_cast< sol::SolverRaytracingStrategy >( aParameterListNonlinearSolver.get< uint >( "NLA_remap_strategy" ) ) )
            , mLoadStepFrequency( aParameterListNonlinearSolver.get< int >( "NLA_remap_load_step_frequency" ) )
            , mIterationFrequency( aParameterListNonlinearSolver.get< int >( "NLA_remap_iteration_frequency" ) )
            , mResidualChangeTolerance( aParameterListNonlinearSolver.get< real >( "NLA_remap_residual_change_tolerance" ) )
            , mPreviousLoadFactor( 0.0 )
            , mLoadStepCounter( 0 )
    {
    }

    bool Solver_Nonconformal_Remapping::requires_remapping( uint aIter, Nonlinear_Solver *aNonLinSolverManager, real &aLoadFactor )
    {
        switch ( mStrategy )
        {
            case sol::SolverRaytracingStrategy::EveryNthLoadStep:
            {
                return check_every_nth_load_step( aLoadFactor );
            }
            case sol::SolverRaytracingStrategy::EveryNthIteration:
            {
                return check_every_nth_iteration( aIter, mIterationFrequency );
            }
            case sol::SolverRaytracingStrategy::ResidualChange:
            {
                return check_relative_residual_change( aNonLinSolverManager );
            }
            case sol::SolverRaytracingStrategy::EveryNthLoadStepOrNthIteration:
            {
                return check_every_nth_load_step( aLoadFactor ) || check_every_nth_iteration( aIter, mIterationFrequency );
            }
            case sol::SolverRaytracingStrategy::MixedNthLoadStepAndResidualChange:
            {
                // use relative residual drop if load factor is 1
                if ( std::abs( 1.0 - aLoadFactor ) < 1e-10 )
                {
                    return check_relative_residual_change( aNonLinSolverManager );
                }
                return check_every_nth_load_step( aLoadFactor );
            }
            case sol::SolverRaytracingStrategy::MixedNthLoadStepAndNthIteration:
            {
                // use iteration if load factor is 1
                if ( std::abs( 1.0 - aLoadFactor ) < 1e-10 )
                {
                    return check_every_nth_iteration( aIter, mIterationFrequency );
                }
                return check_every_nth_load_step( aLoadFactor );
            }
            default:
            {
                return false;
            }
        }
    }

    bool Solver_Nonconformal_Remapping::check_every_nth_load_step( real &aLoadFactor )
    {
        if ( std::abs( aLoadFactor - mPreviousLoadFactor ) > 1e-6 )
        {
            mLoadStepCounter++;
            mPreviousLoadFactor = aLoadFactor;
            return mLoadStepCounter % mLoadStepFrequency == 0;
        }
        return false;
    }

    bool Solver_Nonconformal_Remapping::check_every_nth_iteration( uint aIter, uint aFrequency )
    {
        return ( aIter % aFrequency ) == 0 || aIter == 1;    // always remap at first iteration or every nth iteration
    }

    bool Solver_Nonconformal_Remapping::check_relative_residual_change( Nonlinear_Solver *aNonLinSolverManager )
    {
        real tResidualNorm = NAN;
        if ( aNonLinSolverManager->get_compute_static_residual_flag() )
        {
            tResidualNorm = aNonLinSolverManager->get_static_residual_norm();
        }
        else
        {
            tResidualNorm = aNonLinSolverManager->get_residual_norm();
        }

        if ( mPreviousResidual < 0.0 )    // first iteration
        {
            mPreviousResidual = aNonLinSolverManager->get_residual_norm();
            return false;
        }

        real const tResidualChange = tResidualNorm - mPreviousResidual;    // compute residual change
        mPreviousResidual          = tResidualNorm;
        return tResidualChange > mResidualChangeTolerance;
    }
}    // namespace moris::NLA