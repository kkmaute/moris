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

        bool tBreaker = false;
        Matrix< DDRMat > tTime( 2, 1, tTime_Scalar );

        if (mNonlinearSolver->get_nonlin_solver_type() == NLA::NonlinearSolverType::ARC_LENGTH_SOLVER )
        {
            mNonlinearSolver->set_time_solver_type( this );

            // update time increment (lambda) via the arc length function
            tTime( 1, 0 ) = mLambdaInc;

            if ( mLambdaInc >= 1)
            {
                mLambdaInc = 1.0;
                tBreaker = true;
            }

        }
        else
        {
            tTime_Scalar = tTime_Scalar + tTimeIncrements;
            tTime( 1, 0 ) = tTime_Scalar;

        }
        mSolverInterface->set_time( tTime );

        mSolverInterface->set_solution_vector_prev_time_step( mPrevFullVector );

        mNonlinearSolver->set_time_step_iter( Ik );

        mNonlinearSolver->solve( mFullVector );

        mPrevFullVector->vec_plus_vec( 1.0, *mFullVector, 0.0);

        if (tBreaker)
        {
            MORIS_ASSERT( false, "solve_monolytic_time_system(): lambda value greater than 1 detected...exiting monolithic time loop " );
            break;
        }

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

//-------------------------------------------------------------------------------

void Monolithic_Time_Solver::set_lambda_increment( moris::real aLambdaInc )
{
    mNonlinearSolver->get_my_nonlin_problem()->set_time_value( aLambdaInc );
    mLambdaInc = aLambdaInc;
}

moris::real Monolithic_Time_Solver::get_new_lambda()
{
    return mLambdaInc;
}
