/*
 * cl_NLA_Convergence.cpp
 *
 *  Created on: Sep 21, 2018
 *      Author: schmidt
 */
#include <ctime>

#include "cl_NLA_Convergence.hpp"

#include "cl_SOL_Enums.hpp"
#include "cl_SOL_Dist_Vector.hpp"

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
                moris::sint&         aIt,
                moris::real&         aRefNorm,
                moris::real&         aResNorm,
                const moris::real&   aAssemblyTime,
                const moris::real&   aSolvTime,
                bool&                aHartBreak )
        {
            // initialize return value
            bool tIsConverged = false;

            // compute norm of solution(s)
            Cell< moris::real > solNorm = tNonLinSolver->mNonlinearProblem->get_linearized_problem()->get_free_solver_LHS()->vec_norm2();

            // compute residual norm which is the norm of the RHS of the linearized problem
            aResNorm = tNonLinSolver->mNonlinearProblem->get_linearized_problem()->get_solver_RHS()->vec_norm2()( 0 );

            // set the residual norm to the reference norm for the first iteration
            if ( aIt <= 1 )
            {
                aRefNorm = aResNorm;

                MORIS_LOG_SPEC( "ReferenceNorm", aRefNorm );
            }

            // log residual and solution norms
            MORIS_LOG_SPEC( "ResidualNorm", aResNorm );
            MORIS_LOG_SPEC( "SolutionNorm", solNorm( 0 ) );

            // log relative drop of residual
            if ( aIt > 1 )
            {
                MORIS_LOG_SPEC( "NlinResDrop", aResNorm / aRefNorm );
            }

            // check that residual is valid
            MORIS_ERROR( !( std::isnan( aResNorm ) || std::isinf( aResNorm ) ),
                    "Convergence::check_for_convergence(): Residual contains NAN or INF, exiting!" );

            // check that residual does not exceed upper limit (should defined in MORIS global definitions)
            MORIS_ERROR( !( aResNorm > 1e20 ),
                    "Convergence::check_for_convergence(): Residual Norm has exceeded 1e20" );

            // Check for convergence
            if ( ( aIt > 1 ) && ( aResNorm < aRefNorm * tNonLinSolver->mParameterListNonlinearSolver.get< moris::real >( "NLA_rel_res_norm_drop" ) ) )
            {
                MORIS_LOG_INFO( "NlinResDrop < %6.1e", tNonLinSolver->mParameterListNonlinearSolver.get< moris::real >( "NLA_rel_res_norm_drop" ) );

                tIsConverged = true;
            }
            else if ( ( aIt > 1 ) && ( aResNorm < tNonLinSolver->mParameterListNonlinearSolver.get< moris::real >( "NLA_tot_res_norm" ) ) )
            {
                MORIS_LOG_INFO( "NlinResDrop < %6.1e", tNonLinSolver->mParameterListNonlinearSolver.get< moris::real >( "NLA_rel_res_norm_drop" ) );

                tIsConverged = true;
            }
            else if ( ( aIt > 1 ) && ( aResNorm > aRefNorm * tNonLinSolver->mParameterListNonlinearSolver.get< moris::real >( "NLA_max_rel_res_norm" ) ) )
            {
                MORIS_LOG_INFO( "MaxRelResNorm > %6.1e", tNonLinSolver->mParameterListNonlinearSolver.get< moris::real >( "NLA_max_rel_res_norm" ) );

                tIsConverged = false;

                if ( tNonLinSolver->mParameterListNonlinearSolver.get< bool >( "NLA_hard_break" ) )
                {
                    aIt = tNonLinSolver->mParameterListNonlinearSolver.get< moris::sint >( "NLA_max_iter" );

                    aHartBreak = true;
                }
            }

            return tIsConverged;
        }

        //--------------------------------------------------------------------------------------------------------------------------

        bool
        Convergence::check_for_convergence(
                Nonlinear_Algorithm* tNonLinSolver,
                moris::sint&         aIt,
                moris::real&         aRefNorm,
                moris::real&         aResNorm,
                const moris::real&   aSolvTime,
                bool&                aHartBreak )
        {
            // initialize return value
            bool tIsConverged = false;

            // log reference norm
            if ( aIt <= 1 )
            {
                MORIS_LOG_SPEC( "ReferenceNorm", aRefNorm );
            }

            // log residual and solution norms
            MORIS_LOG_SPEC( "ResidualNorm", aResNorm );

            // log relative drop of residual
            if ( aIt > 1 )
            {
                MORIS_LOG_SPEC( "NlinResDrop", aResNorm / aRefNorm );
            }

            // check that residual is valid
            MORIS_ERROR( !( std::isnan( aResNorm ) || std::isinf( aResNorm ) ),
                    "Convergence::check_for_convergence(): Residual contains NAN or INF, exiting!" );

            // check that residual does not exceed upper limit (should defined in MORIS global definitions)
            MORIS_ERROR( !( aResNorm > 1e20 ),
                    "Convergence::check_for_convergence(): Residual Norm has exceeded 1e20" );

            // Check for convergence
            if ( ( aIt >= 1 ) && ( aResNorm < aRefNorm * tNonLinSolver->mParameterListNonlinearSolver.get< moris::real >( "NLA_rel_res_norm_drop" ) ) )
            {
                MORIS_LOG_INFO( "NlinResDrop < %6.1e", tNonLinSolver->mParameterListNonlinearSolver.get< moris::real >( "NLA_rel_res_norm_drop" ) );

                tIsConverged = true;
            }
            else if ( ( aIt >= 1 ) && ( aResNorm < tNonLinSolver->mParameterListNonlinearSolver.get< moris::real >( "NLA_tot_res_norm" ) ) )
            {
                MORIS_LOG_INFO( "NlinResNorm < %6.1e", tNonLinSolver->mParameterListNonlinearSolver.get< moris::real >( "NLA_tot_res_norm" ) );

                tIsConverged = true;
            }
            else if ( ( aIt >= 1 )
                      && ( aResNorm > aRefNorm * tNonLinSolver->mParameterListNonlinearSolver.get< moris::real >( "NLA_max_rel_res_norm" ) ) )
            {
                MORIS_LOG_INFO( "MaxRelResNorm > %6.1e", tNonLinSolver->mParameterListNonlinearSolver.get< moris::real >( "NLA_max_rel_res_norm" ) );

                tIsConverged = false;

                if ( tNonLinSolver->mParameterListNonlinearSolver.get< bool >( "NLA_hard_break" ) )
                {
                    aIt        = tNonLinSolver->mParameterListNonlinearSolver.get< moris::sint >( "NLA_max_iter" );
                    aHartBreak = true;
                }
            }

            return tIsConverged;
        }
    }    // namespace NLA
}    // namespace moris
