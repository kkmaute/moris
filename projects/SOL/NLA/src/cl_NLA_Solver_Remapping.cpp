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
#include "cl_Param_List.hpp"

namespace moris::NLA
{
    Solver_Nonconformal_Remapping::Solver_Nonconformal_Remapping( ParameterList &aParameterListNonlinearSolver )
            : mRaytracingStrategy( static_cast< sol::SolverRaytracingStrategy >( aParameterListNonlinearSolver.get< uint >( "NLA_remap_strategy" ) ) )
            , mRaytracingFrequency( aParameterListNonlinearSolver.get< int >( "NLA_remap_frequency" ) )
            , mRaytracingRelResidualDrop( aParameterListNonlinearSolver.get< real >( "NLA_remap_relative_residual_drop" ) )
            , mPreviousLoadFactor( 0.0 )
            , mLoadStepCounter( 0 )
    {
    }

    bool Solver_Nonconformal_Remapping::requires_remapping( uint aIter, Nonlinear_Solver *aNonLinSolverManager, real &aLoadFactor )
    {
        switch ( mRaytracingStrategy )
        {
            case sol::SolverRaytracingStrategy::EveryNthLoadStep:
            {
                return check_every_nth_load_step( aLoadFactor );
            }
            case sol::SolverRaytracingStrategy::EveryNthIteration:
            {
                return check_every_nth_iteration( aIter );
            }
            case sol::SolverRaytracingStrategy::RelativeResidualDrop:
            {
                return check_relative_residual_drop( aNonLinSolverManager );
            }
            case sol::SolverRaytracingStrategy::MixedNthLoadStepAndResidualDrop:
            {
                // use relative residual drop if load factor is 1
                if ( std::abs( 1.0 - aLoadFactor ) < 1e-10 )
                {
                    return check_relative_residual_drop( aNonLinSolverManager );
                }
                return check_every_nth_load_step( aLoadFactor );
            }
            case sol::SolverRaytracingStrategy::MixedNthLoadStepAndNthIteration:
            {
                // use iteration if load factor is 1
                if ( std::abs( 1.0 - aLoadFactor ) < 1e-10 )
                {
                    return check_every_nth_iteration( aIter );
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
        if ( std::abs(aLoadFactor - mPreviousLoadFactor) > 1e-6 )
        {
            mLoadStepCounter++;
            mPreviousLoadFactor = aLoadFactor;
            return mLoadStepCounter % mRaytracingFrequency == 0;
        }
        return false;
    }

    bool Solver_Nonconformal_Remapping::check_every_nth_iteration( uint aIter ) const
    {
        return aIter % mRaytracingFrequency == 0;
    }

    bool Solver_Nonconformal_Remapping::check_relative_residual_drop( Nonlinear_Solver *aNonLinSolverManager ) const
    {
        real tRefNorm;
        real tResNorm;

        // use static residual if available
        if ( aNonLinSolverManager->get_compute_static_residual_flag() )
        {
            tRefNorm = aNonLinSolverManager->get_static_ref_norm();
            tResNorm = aNonLinSolverManager->get_static_residual_norm();
        }
        else
        {
            tRefNorm = aNonLinSolverManager->get_ref_norm();
            tResNorm = aNonLinSolverManager->get_residual_norm();
        }
        real const tRelResNorm = tResNorm / tRefNorm;    // compute relative residual
        return tRelResNorm < mRaytracingRelResidualDrop;
    }
}    // namespace moris::NLA