/*
 * cl_NLA_Convergence.cpp
 *
 *  Created on: Sep 21, 2018
 *      Author: schmidt
 */
#include <ctime>

#include "cl_NLA_Convergence.hpp"

#include "cl_DistLinAlg_Enums.hpp"
#include "cl_Vector.hpp"

#include "cl_Communication_Tools.hpp"

namespace moris
{
    namespace NLA
    {
//--------------------------------------------------------------------------------------------------------------------------
//    bool Newton_Solver::check_for_convergence(       moris::sint & aIt,
//                                                     moris::real & aRefNorm,
//                                               const moris::real & aAssemblyTime,
//                                                     bool        & aHartBreak )
//    {
//        bool tIsConverged = false;
//        moris::real resNorm = 0.0;
//        moris::real solNorm = mVectorFreeSol->vec_norm2();
//
//        resNorm = mLinearSolver->get_solver_RHS()->vec_norm2();
//
//        if ( aIt <= 1)
//        {
//            aRefNorm = resNorm;
//            std::printf( " ... refNorm for pseudo-time step is %+1.15e\n", aRefNorm );
//
//            if ( par_rank() == 0 )
//            {
//                fprintf( stdout, "        NlinIt  |  NlinResNorm            |  NlinResDrop  |  SolVecNorm             |" );
//
//                fprintf( stdout, "|  LinAsmTime  |  NewItrTime\n" );
//
//                // Print solution vector norm solNorm before the first Newton solve
//                fprintf( stdout, "         %-5i  |  %-15.15e  |  %-11.5e  |  %-10.15e  |", 1, resNorm, 0.0, solNorm );
//
//                fprintf( stdout, "|  %9.4e  |  %9.4e \n", aAssemblyTime, 0.0 );
//            }
//        }
//
//        MORIS_ERROR( !( std::isnan( resNorm ) || std::isinf( resNorm )), "Newton_Solver::check_for_convergence(): Residual contains NAN or INF, exiting!");
//
//        MORIS_ERROR( !( resNorm > 1e20 ), "Newton_Solver::check_for_convergence(): Residual Norm has exceeded 1e20");
//
//        // Check for convergence
//        if ( ( aIt > 1 ) && ( ( resNorm/aRefNorm ) < mParameterListNonlinearSolver.get< moris::real >( "NLA_tot_res_norm_drop" ) ) )
//        {
//            if ( par_rank() == 0 )
//            {
//                fprintf( stdout, "         %-5i  |  %-15.15e  |  %-11.5e  |  %-10.15e  |", aIt, resNorm, ( resNorm/aRefNorm ), solNorm );
//
//                fprintf( stdout, "|  %-10.4e  |  (NlinResDrop < %6.1e)\n", aAssemblyTime, mParameterListNonlinearSolver.get< moris::real >( "NLA_tot_res_norm_drop" ) );
//            }
//
//            tIsConverged = true;
//        }
//        else if ( ( aIt > 1 ) && ( resNorm < mParameterListNonlinearSolver.get< moris::real >( "NLA_tot_res_norm" ) ) )
//        {
//            if ( par_rank() == 0 )
//            {
//                fprintf( stdout, "         %-5i  |  %-15.15e  |  %-11.5e  |  %-10.15e  |", aIt, resNorm, (resNorm/aRefNorm), solNorm );
//
//                fprintf( stdout, "|  %-10.4e  |  (NlinResNorm < %6.1e)\n", aAssemblyTime, mParameterListNonlinearSolver.get< moris::real >( "NLA_tot_res_norm" ) );
//            }
//
//            tIsConverged = true;
//        }
//        else if ( ( aIt > 1 ) && ( ( resNorm/aRefNorm ) > mParameterListNonlinearSolver.get< moris::real >( "NLA_max_res_norm_drop" ) ) )
//        {
//            // case for residual drop getting too big, not converged, need to retry
//            if ( par_rank() == 0 )
//            {
//                fprintf( stdout, "         %-5i  |  %-15.15e  |  %-11.5e  |  %-10.15e  |", aIt, resNorm, ( resNorm/aRefNorm ), solNorm );
//
//                fprintf( stdout, "|  %-10.4e  |  (MaxResNormDrop > %6.1e)\n", aAssemblyTime, mParameterListNonlinearSolver.get< moris::real >( "NLA_max_res_norm_drop" ) );
//            }
//
//            tIsConverged = false;
//
//            if ( mParameterListNonlinearSolver.get< bool >( "NLA_hard_break" ) )
//            {
//                aIt = mParameterListNonlinearSolver.get< moris::sint >( "NLA_max_iter" );
//                aHartBreak = true;
//            }
//        }
//
//        return tIsConverged;
//    }
//--------------------------------------------------------------------------------------------------------------------------

    }
}
