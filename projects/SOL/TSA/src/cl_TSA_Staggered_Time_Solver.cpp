/*
 * cl_TSA_Staggered_Time_Solver.cpp
 *
 *  Created on: Feb 11, 2019
 *      Author: schmidt
 */

#include "cl_TSA_Staggered_Time_Solver.hpp"

#include "cl_TSA_Time_Solver.hpp"


using namespace moris;
using namespace tsa;
//-------------------------------------------------------------------------------

void Staggered_Time_Solver::solve_staggered_time_system()
{
    this->finalize();

    moris::sint tMaxIts = 1;
    moris::uint tNumTimeSystems = 2;
    //bool tIsConverged            = false;

    // time loop
    for ( sint It = 1; It <= tMaxIts; ++It )
    {
        moris::real tMaxNewTime      = 0.0;

        //get_time_problem()
        clock_t tLoopStartTime = clock();

        // Loop over all time systems
        for (uint Ik = 0 ; Ik < tNumTimeSystems; Ik++)
        {
            mMyTimeSolver->get_sub_time_solver( Ik )->solve( mFullVector );

             //mTimeSolverList( Ik )->solve( mFullVector );
        } // end loop over all time sub-systems

        tMaxNewTime = this->calculate_time_needed( tLoopStartTime );

        std::cout<< "Time spend on time system: " << tMaxNewTime <<std::endl;

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

void Staggered_Time_Solver::solve( Dist_Vector * aFullVector )
 {
     mFullVector = aFullVector;

     this->solve_staggered_time_system();
 }

 //-------------------------------------------------------------------------------

 void Staggered_Time_Solver::solve()
 {
     mIsMasterTimeSolver =  true;

     this->solve_staggered_time_system();
 }
