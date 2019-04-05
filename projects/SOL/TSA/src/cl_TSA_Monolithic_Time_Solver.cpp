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

    uint tTimeSteps = 1000;
    moris::real tTime = 0;

    real tStepSize = 0.01;


    for ( uint Ik = 0; Ik < tTimeSteps; Ik++ )
    {
        tTime = tTime + tStepSize;
        mSolverInterface->set_time( tTime );

        mSolverInterface->set_solution_vector_prev_time_step( mPrevFullVector );

        mNonlinearSolver->solve( mFullVector );

        mPrevFullVector->vec_plus_vec( 1.0, *mFullVector, 0.0);

        this->perform_mapping();
    }
}

//-------------------------------------------------------------------------------

void Monolithic_Time_Solver::perform_mapping()
{
    Matrix< DDRMat > tMat;
    Matrix< DDSMat > tMatRows1 = mSolverInterface->get_time_level_Ids_minus();
    Matrix< DDSMat > tMatRows2 = mSolverInterface->get_time_level_Ids_plus();

    mPrevFullVector->extract_my_values( 1, tMatRows1, 0 , tMat );

    mPrevFullVector->sum_into_global_values( 1, tMatRows2, tMat );

    mPrevFullVector->vector_global_asembly();
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
