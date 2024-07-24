/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_NLA_Convergence.cpp
 *
 */
#include <ctime>

#include "cl_SOL_Enums.hpp"

#include "cl_SOL_Dist_Vector.hpp"
#include "cl_NLA_Convergence.hpp"

#include "cl_DLA_Linear_Problem.hpp"

#include "cl_Communication_Tools.hpp"

#include "cl_Logger.hpp"
#include "cl_Tracer.hpp"

namespace moris
{
    namespace NLA
    {
        //--------------------------------------------------------------------------------------------------------------------------
        bool
        Convergence::check_for_convergence(
                Nonlinear_Algorithm* tNonLinSolver,
                moris::sint          aIt,
                moris::sint          aMaxIter,
                bool&                aHartBreak )
        {
            // initialize return value
            bool tIsConverged = false;

            // compute norm of solution(s)
            Vector< moris::real > solNorm = tNonLinSolver->mNonlinearProblem->get_linearized_problem()->get_free_solver_LHS()->vec_norm2();

            // compute residual norm which is the norm of the RHS of the linearized problem
            real tResNorm = tNonLinSolver->mNonlinearProblem->get_linearized_problem()->get_solver_RHS()->vec_norm2()( 0 );

            // initialize reference residual norm
            real tRefNorm = 0.0;

            // set the residual norm to the reference norm for the first iteration
            if ( aIt <= mRefIterationID )
            {
                // set reference norm
                tRefNorm      = tResNorm;
                mPreviousNorm = tRefNorm;

                // store reference with solver manager of solver algorithm
                tNonLinSolver->mMyNonLinSolverManager->set_ref_norm( tRefNorm );

                // compute static residual norm
                if ( tNonLinSolver->mMyNonLinSolverManager->get_compute_static_residual_flag() )
                {
                    // compute and store static residual
                    real tStaticResNorm = tNonLinSolver->mNonlinearProblem->get_static_residual_norm();

                    // store reference with solver manager of solver algorithm
                    tNonLinSolver->mMyNonLinSolverManager->set_static_ref_norm( tStaticResNorm );

                    // log reference norm for static residual
                    MORIS_LOG_SPEC( "StaticReferenceNorm", tStaticResNorm );
                }

                // log reference norm
                MORIS_LOG_SPEC( "ReferenceNorm", tRefNorm );
            }
            else
            {
                // get reference norm
                tRefNorm = tNonLinSolver->mMyNonLinSolverManager->get_ref_norm();
            }

            // log residual and solution norms
            MORIS_LOG_SPEC( "ResidualNorm", tResNorm );
            MORIS_LOG_SPEC( "SolutionNorm", solNorm( 0 ) );

            // store current residual with solver manager of solver algorithm
            tNonLinSolver->mMyNonLinSolverManager->set_residual_norm( tResNorm );

            // log relative drop of residual
            if ( aIt > mRefIterationID )
            {
                MORIS_LOG_SPEC( "RelResidualDrop", tResNorm / tRefNorm );
                MORIS_LOG_SPEC( "RelResidualChange", ( tResNorm - mPreviousNorm ) / tRefNorm );
                mPreviousNorm = tResNorm;
            }

            // check that residual is valid
            MORIS_ERROR( !( std::isnan( tResNorm ) || std::isinf( tResNorm ) ),
                    "Convergence::check_for_convergence(): Residual contains NAN or INF, exiting!" );

            // check that residual does not exceed upper limit (should defined in MORIS global definitions)
            MORIS_ERROR( !( tResNorm > 1e20 ),
                    "Convergence::check_for_convergence(): Residual Norm has exceeded 1e20" );

            // Check for convergence
            if ( ( aIt > mRefIterationID ) && ( tResNorm < tRefNorm * tNonLinSolver->mParameterListNonlinearSolver.get< moris::real >( "NLA_rel_res_norm_drop" ) ) )
            {
                MORIS_LOG_INFO( "NlinResDrop < %6.1e", tNonLinSolver->mParameterListNonlinearSolver.get< moris::real >( "NLA_rel_res_norm_drop" ) );

                tIsConverged = true;
            }
            else if ( ( aIt > mRefIterationID ) && ( tResNorm < tNonLinSolver->mParameterListNonlinearSolver.get< moris::real >( "NLA_tot_res_norm" ) ) )
            {
                MORIS_LOG_INFO( "NlinResDrop < %6.1e", tNonLinSolver->mParameterListNonlinearSolver.get< moris::real >( "NLA_rel_res_norm_drop" ) );

                tIsConverged = true;
            }
            else if ( ( aIt > mRefIterationID ) && ( tResNorm > tRefNorm * tNonLinSolver->mParameterListNonlinearSolver.get< moris::real >( "NLA_max_rel_res_norm" ) ) )
            {
                MORIS_LOG_INFO( "MaxRelResNorm > %6.1e", tNonLinSolver->mParameterListNonlinearSolver.get< moris::real >( "NLA_max_rel_res_norm" ) );

                tIsConverged = false;

                if ( tNonLinSolver->mParameterListNonlinearSolver.get< bool >( "NLA_hard_break" ) )
                {
                    aHartBreak = true;
                }
            }

            // Store data at final iteration; skip if aMaxIter equals mRefIterationID
            if ( tIsConverged or aHartBreak or ( aIt == aMaxIter && aMaxIter > mRefIterationID ) )
            {
                // store relative number of iterations
                tNonLinSolver->mMyNonLinSolverManager->set_relative_number_iterations( ( (real)aIt ) / aMaxIter );

                // compute static residual norm
                if ( tNonLinSolver->mMyNonLinSolverManager->get_compute_static_residual_flag() )
                {
                    // compute and store static residual
                    real tStaticResNorm = tNonLinSolver->mNonlinearProblem->get_static_residual_norm();

                    // store reference with solver manager of solver algorithm
                    tNonLinSolver->mMyNonLinSolverManager->set_static_residual_norm( tStaticResNorm );

                    // log norm pf static residual
                    MORIS_LOG_SPEC( "StaticResidualNorm", tStaticResNorm );
                }
            }

            return tIsConverged;
        }

        //--------------------------------------------------------------------------------------------------------------------------

        bool
        Convergence::check_for_convergence(
                Nonlinear_Algorithm* tNonLinSolver,
                moris::sint          aIt,
                bool&                aHartBreak )
        {
            // initialize return value
            bool tIsConverged = false;

            // get current and reference residual norms
            real tResNorm = tNonLinSolver->mMyNonLinSolverManager->get_residual_norm();
            real tRefNorm = tNonLinSolver->mMyNonLinSolverManager->get_ref_norm();

            // log reference norm
            if ( aIt <= mRefIterationID )
            {
                MORIS_LOG_SPEC( "ReferenceNorm", tRefNorm );

                // log static reference residual norm
                if ( tNonLinSolver->mMyNonLinSolverManager->get_compute_static_residual_flag() )
                {
                    const real tStaticRefNorm = tNonLinSolver->mMyNonLinSolverManager->get_static_ref_norm();

                    MORIS_LOG_SPEC( "StaticReferenceNorm", tStaticRefNorm );
                }
            }

            // log residual and solution norms
            MORIS_LOG_SPEC( "ResidualNorm", tResNorm );

            // log relative drop of residual
            if ( aIt > mRefIterationID )
            {
                MORIS_LOG_SPEC( "RelResidualDrop", tResNorm / ( tRefNorm + MORIS_REAL_EPS ) );
            }

            // check that residual is valid
            MORIS_ERROR( !( std::isnan( tResNorm ) || std::isinf( tResNorm ) ),
                    "Convergence::check_for_convergence(): Residual contains NAN or INF, exiting!" );

            // check that residual does not exceed upper limit (should defined in MORIS global definitions)
            MORIS_ERROR( !( tResNorm > 1e20 ),
                    "Convergence::check_for_convergence(): Residual Norm has exceeded 1e20" );

            // log static residual norm and relative drop in static residual
            if ( tNonLinSolver->mMyNonLinSolverManager->get_compute_static_residual_flag() )
            {
                const real tStaticRefNorm = tNonLinSolver->mMyNonLinSolverManager->get_static_ref_norm();
                const real tStaticResNorm = tNonLinSolver->mMyNonLinSolverManager->get_static_residual_norm();

                MORIS_LOG_SPEC( "StaticResidualNorm", tStaticResNorm );

                if ( aIt > mRefIterationID )
                {
                    MORIS_LOG_SPEC( "StaticResDrop", tStaticResNorm / ( tStaticRefNorm + MORIS_REAL_EPS ) );
                }
            }

            // Check for convergence
            if ( ( aIt >= mRefIterationID ) && ( tResNorm < tRefNorm * tNonLinSolver->mParameterListNonlinearSolver.get< moris::real >( "NLA_rel_res_norm_drop" ) ) )
            {
                MORIS_LOG_INFO( "NlinResDrop < %6.1e", tNonLinSolver->mParameterListNonlinearSolver.get< moris::real >( "NLA_rel_res_norm_drop" ) );

                tIsConverged = true;
            }
            else if ( ( aIt >= mRefIterationID ) && ( tResNorm < tNonLinSolver->mParameterListNonlinearSolver.get< moris::real >( "NLA_tot_res_norm" ) ) )
            {
                MORIS_LOG_INFO( "NlinResNorm < %6.1e", tNonLinSolver->mParameterListNonlinearSolver.get< moris::real >( "NLA_tot_res_norm" ) );

                tIsConverged = true;
            }
            else if ( ( aIt >= mRefIterationID )
                      && ( tResNorm > tRefNorm * tNonLinSolver->mParameterListNonlinearSolver.get< moris::real >( "NLA_max_rel_res_norm" ) ) )
            {
                MORIS_LOG_INFO( "MaxRelResNorm > %6.1e", tNonLinSolver->mParameterListNonlinearSolver.get< moris::real >( "NLA_max_rel_res_norm" ) );

                tIsConverged = false;

                if ( tNonLinSolver->mParameterListNonlinearSolver.get< bool >( "NLA_hard_break" ) )
                {
                    aHartBreak = true;
                }
            }

            return tIsConverged;
        }
    }    // namespace NLA
}    // namespace moris
