/*
 * cl_NLA_Convergence.cpp
 *
 *  Created on: Sep 21, 2018
 *      Author: schmidt
 */
#include <ctime>

#include "cl_NLA_Convergence.hpp"

#include "cl_DLA_Enums.hpp"
#include "cl_Vector.hpp"

#include "cl_Communication_Tools.hpp"

#include "cl_Logger.hpp"

namespace moris
{
    namespace NLA
    {
//--------------------------------------------------------------------------------------------------------------------------
    bool Convergence::check_for_convergence(         Nonlinear_Solver * tNonLinSolver,
                                                     moris::sint & aIt,
                                                     moris::real & aRefNorm,
                                                     moris::real & aResNorm,
                                               const moris::real & aAssemblyTime,
                                               const moris::real & aSolvTime,
                                                     bool        & aHartBreak )
    {
        bool tIsConverged = false;

        moris::real solNorm = tNonLinSolver->mNonlinearProblem->get_full_vector()->vec_norm2();

        aResNorm = tNonLinSolver->mNonlinearProblem->get_linearized_problem()->get_solver_RHS()->vec_norm2();
        if ( aIt <= 1)
        {
            aRefNorm = aResNorm;
            MORIS_LOG( "--------------------------------------------------------------------------------\n");
            MORIS_LOG( " Newton ... refNorm for pseudo-time step is %+1.15e\n", aRefNorm );
            MORIS_LOG( "--------------------------------------------------------------------------------\n" );

            if ( par_rank() == 0 )
            {
                MORIS_LOG( "        NlinIt  |  NlinResNorm            |  NlinResDrop  |  SolVecNorm             ||  LinAsmTime  |  NewItrTime\n" );
                MORIS_LOG( "         %-5i  |  %-15.15e  |  %-11.5e  |  %-10.15e  |", 1, aResNorm, 0.0, solNorm );
                MORIS_LOG( "|  %9.4e  |  %9.4e \n", aAssemblyTime, 0.0 );
            }
        }

        MORIS_ERROR( !( std::isnan( aResNorm ) || std::isinf( aResNorm )), "Convergence::check_for_convergence(): Residual contains NAN or INF, exiting!");

        MORIS_ERROR( !( aResNorm > 1e20 ), "Convergence::check_for_convergence(): Residual Norm has exceeded 1e20");

        // Check for convergence
        if ( ( aIt > 1 ) && ( ( aResNorm/aRefNorm ) < tNonLinSolver->mParameterListNonlinearSolver.get< moris::real >( "NLA_tot_res_norm_drop" ) ) )
        {
            if ( par_rank() == 0 )
            {
                MORIS_LOG( "         %-5i  |  %-15.15e  |  %-11.5e  |  %-10.15e  |", aIt, aResNorm, ( aResNorm/aRefNorm ), solNorm );
                MORIS_LOG( "|  %-10.4e  |  (NlinResDrop < %6.1e)\n", aAssemblyTime, tNonLinSolver->mParameterListNonlinearSolver.get< moris::real >( "NLA_tot_res_norm_drop" ) );
            }

            tIsConverged = true;
        }
        else if ( ( aIt > 1 ) && ( aResNorm < tNonLinSolver->mParameterListNonlinearSolver.get< moris::real >( "NLA_tot_res_norm" ) ) )
        {
            if ( par_rank() == 0 )
            {
                MORIS_LOG( "         %-5i  |  %-15.15e  |  %-11.5e  |  %-10.15e  |", aIt, aResNorm, (aResNorm/aRefNorm), solNorm );
                MORIS_LOG( "|  %-10.4e  |  (NlinResNorm < %6.1e)\n", aAssemblyTime, tNonLinSolver->mParameterListNonlinearSolver.get< moris::real >( "NLA_tot_res_norm" ) );
            }

            tIsConverged = true;
        }
        else if ( ( aIt > 1 ) && ( ( aResNorm/aRefNorm ) > tNonLinSolver->mParameterListNonlinearSolver.get< moris::real >( "NLA_max_res_norm_drop" ) ) )
        {
            // case for residual drop getting too big, not converged, need to retry
            if ( par_rank() == 0 )
            {
                MORIS_LOG( "         %-5i  |  %-15.15e  |  %-11.5e  |  %-10.15e  |", aIt, aResNorm, ( aResNorm/aRefNorm ), solNorm );
                MORIS_LOG( "|  %-10.4e  |  (MaxResNormDrop > %6.1e)\n", aAssemblyTime, tNonLinSolver->mParameterListNonlinearSolver.get< moris::real >( "NLA_max_res_norm_drop" ) );
            }

            tIsConverged = false;

            if ( tNonLinSolver->mParameterListNonlinearSolver.get< bool >( "NLA_hard_break" ) )
            {
                aIt = tNonLinSolver->mParameterListNonlinearSolver.get< moris::sint >( "NLA_max_iter" );
                aHartBreak = true;
            }
        }
        else if( ( aIt > 1 ) )
        {
            if ( par_rank() == 0 )
            {
                MORIS_LOG( "         %-5i  |  %-15.15e  |  %-11.5e  |  %-10.15e  |", aIt, aResNorm, (aResNorm/aRefNorm), solNorm  );
                MORIS_LOG( "|  %9.4e  |  %9.4e \n", aAssemblyTime, aSolvTime );
            }
        }

        return tIsConverged;
    }
//--------------------------------------------------------------------------------------------------------------------------

    bool Convergence::check_for_convergence(         Nonlinear_Solver * tNonLinSolver,
                                                     moris::sint & aIt,
                                                     moris::real & aRefNorm,
                                                     moris::real & aResNorm,
                                               const moris::real & aSolvTime,
                                                     bool        & aHartBreak )
    {
        bool tIsConverged = false;

        if ( aIt <= 1)
        {
            aRefNorm = aResNorm;
            MORIS_LOG( "--------------------------------------------------------------------------------\n");
            MORIS_LOG( " NLBGS ... refNorm for pseudo-time step is %+1.15e\n", aRefNorm );
            MORIS_LOG( "--------------------------------------------------------------------------------\n" );

            if ( par_rank() == 0 )
            {
                MORIS_LOG( "        NlinIt  |  NlinResNorm            |  NlinResDrop  ||  NewItrTime\n" );
                MORIS_LOG( "         %-5i  |  %-15.15e  |  %-11.5e  |", 1, aResNorm, 0.0 );
                MORIS_LOG( "|  %9.4e \n",  0.0 );
            }
        }

        MORIS_ERROR( !( std::isnan( aResNorm ) || std::isinf( aResNorm )), "Convergence::check_for_convergence(): Residual contains NAN or INF, exiting!");

        MORIS_ERROR( !( aResNorm > 1e20 ), "Convergence::check_for_convergence(): Residual Norm has exceeded 1e20");

        // Check for convergence
        if ( ( aIt > 1 ) && ( ( aResNorm/aRefNorm ) < tNonLinSolver->mParameterListNonlinearSolver.get< moris::real >( "NLA_tot_res_norm_drop" ) ) )
        {
            if ( par_rank() == 0 )
            {
                MORIS_LOG( "         %-5i  |  %-15.15e  |  %-11.5e  |", aIt, aResNorm, ( aResNorm/aRefNorm ) );
                MORIS_LOG( "|  %-10.4e  |  (NlinResDrop < %6.1e)\n", aSolvTime, tNonLinSolver->mParameterListNonlinearSolver.get< moris::real >( "NLA_tot_res_norm_drop" ) );
            }

            tIsConverged = true;
        }
        else if ( ( aIt > 1 ) && ( aResNorm < tNonLinSolver->mParameterListNonlinearSolver.get< moris::real >( "NLA_tot_res_norm" ) ) )
        {
            if ( par_rank() == 0 )
            {
                MORIS_LOG( "         %-5i  |  %-15.15e  |  %-11.5e  |", aIt, aResNorm, (aResNorm/aRefNorm) );
                MORIS_LOG( "|  %-10.4e  |  (NlinResNorm < %6.1e)\n", aSolvTime, tNonLinSolver->mParameterListNonlinearSolver.get< moris::real >( "NLA_tot_res_norm" ) );
            }

            tIsConverged = true;
        }
        else if ( ( aIt > 1 ) && ( ( aResNorm/aRefNorm ) > tNonLinSolver->mParameterListNonlinearSolver.get< moris::real >( "NLA_max_res_norm_drop" ) ) )
        {
            // case for residual drop getting too big, not converged, need to retry
            if ( par_rank() == 0 )
            {
                MORIS_LOG( "         %-5i  |  %-15.15e  |  %-11.5e  |", aIt, aResNorm, ( aResNorm/aRefNorm ) );
                MORIS_LOG( "|  %-10.4e  |  (MaxResNormDrop > %6.1e)\n", aSolvTime, tNonLinSolver->mParameterListNonlinearSolver.get< moris::real >( "NLA_max_res_norm_drop" ) );
            }

            tIsConverged = false;

            if ( tNonLinSolver->mParameterListNonlinearSolver.get< bool >( "NLA_hard_break" ) )
            {
                aIt = tNonLinSolver->mParameterListNonlinearSolver.get< moris::sint >( "NLA_max_iter" );
                aHartBreak = true;
            }
        }
        else if( ( aIt > 1 ) )
        {
            if ( par_rank() == 0 )
            {
                MORIS_LOG( "         %-5i  |  %-15.15e  |  %-11.5e  |", aIt, aResNorm, (aResNorm/aRefNorm) );
                MORIS_LOG( "|  %9.4e  |  %9.4e \n", aSolvTime, aSolvTime );
            }
        }

        return tIsConverged;
    }

    }
}
