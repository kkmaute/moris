/*
 * cl_TSA_Monolithic_Time_Solver.cpp
 *
 *  Created on: Feb 11, 2019
 *      Author: schmidt
 */

#include "cl_TSA_Monolithic_Time_Solver.hpp"

using namespace moris;
using namespace tsa;
//-------------------------------------------------------------------------------

void Monolithic_Time_Solver::solve_monolytic_time_system()
{
    this->finalize();

    moris::real tTime_Scalar = 0;
    sint tTimeSteps = mParameterListTimeSolver.get< moris::sint >( "TSA_Num_Time_Steps" );
    moris::real tTimeFrame = mParameterListTimeSolver.get< moris::real >( "TSA_Time_Frame" );
    moris::real tTimeIncrements = tTimeFrame / tTimeSteps;


    for ( sint Ik = 0; Ik < tTimeSteps; Ik++ )
    {
        Matrix< DDRMat > tTime( 2, 1, tTime_Scalar );
        tTime_Scalar = tTime_Scalar + tTimeIncrements;
        tTime( 1, 0 ) = tTime_Scalar;
        
        mSolverInterface->set_time( tTime );

        mSolverInterface->set_solution_vector_prev_time_step( mPrevFullVector );

        mNonlinearSolver->solve( mFullVector );

        mPrevFullVector->vec_plus_vec( 1.0, *mFullVector, 0.0);

        mSolverInterface->perform_mapping();
    }
}

//-------------------------------------------------------------------------------

void Monolithic_Time_Solver::solve( Dist_Vector * aFullVector )
{
    mFullVector = aFullVector;

    this->solve_monolytic_time_system();
}

//-------------------------------------------------------------------------------

void Monolithic_Time_Solver::solve()
{
    mIsMasterTimeSolver =  true;

    this->solve_monolytic_time_system();
}
