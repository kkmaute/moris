/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_TSA_Staggered_Time_Solver.cpp
 *
 */

#include "cl_TSA_Staggered_Time_Solver.hpp"
#include "cl_TSA_Time_Solver.hpp"

#include "cl_SOL_Dist_Vector.hpp"
#include "cl_DLA_Solver_Interface.hpp"

// Logging package
#include "cl_Logger.hpp"
#include "cl_Tracer.hpp"

using namespace moris;
using namespace tsa;
//-------------------------------------------------------------------------------

void Staggered_Time_Solver::solve_staggered_time_system( Vector< sol::Dist_Vector * > & aFullVector )
{
    // trace this solve
    Tracer tTracer( "TimeSolverAlgorithm", "Staggered", "Solve" );

    this->finalize();

    moris::sint tMaxIts = 1;
    moris::uint tNumTimeSystems = 2;
    //bool tIsConverged            = false;

    // time loop
    for ( sint It = 1; It <= tMaxIts; ++It )
    {
        // log number of time steps
        MORIS_LOG_ITERATION();

        moris::real tMaxNewTime      = 0.0;

        //get_time_problem()
        clock_t tLoopStartTime = clock();

        // Loop over all time systems
        for (uint Ik = 0 ; Ik < tNumTimeSystems; Ik++)
        {
            mMyTimeSolver->get_sub_time_solver( Ik )->solve( aFullVector );

        } // end loop over all time sub-systems

        tMaxNewTime = this->calculate_time_needed( tLoopStartTime );

        std::cout<< "Time spend on time system: " << tMaxNewTime <<std::endl;

        mMyTimeSolver->prepare_sol_vec_for_next_time_step();

//                bool tHartBreak = false;
//
//                this->compute_norms( It );
//
//                Convergence tConvergence;
//
//                tIsConverged = tConvergence.check_for_convergence( this,
//                                                                   It,
//                                                                   mMyNonLinSolverManager->get_ref_norm(),
//                                                                   mMyNonLinSolverManager->get_residual_norm(),
//                                                                   tMaxNewTime,
//                                                                   tHartBreak);
//
//                if ( tIsConverged )
//                {
//                    if ( tHartBreak )
//                    {
//                        continue;
//                    }
//                    break;
//                }
    } // end time loop
}

//-------------------------------------------------------------------------------

void Staggered_Time_Solver::solve( Vector< sol::Dist_Vector * > & aFullVector )
 {
     this->solve_staggered_time_system( aFullVector );
 }

