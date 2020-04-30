/*
 * cl_TSA_Monolithic_Time_Solver.cpp
 *
 *  Created on: Feb 11, 2019
 *      Author: schmidt
 */

#include "cl_TSA_Monolithic_Time_Solver.hpp"
#include "cl_SOL_Dist_Vector.hpp"
#include "cl_DLA_Solver_Interface.hpp"
#include "cl_Matrix_Vector_Factory.hpp"
#include "cl_NLA_Nonlinear_Solver.hpp"
#include "cl_NLA_Nonlinear_Problem.hpp"

// for detailed logging
#include "cl_Logger.hpp"
#include "cl_Tracer.hpp"
#include "cl_Tracer_Enums.hpp"


using namespace moris;
using namespace tsa;
//-------------------------------------------------------------------------------

void Monolithic_Time_Solver::solve_monolytic_time_system()
{
    // trace this solve
    Tracer tTracer(EntityBase::TimeSolver, EntityType::Monolythic, EntityAction::Solve);

    this->finalize();

    moris::real tTime_Scalar = 0;
    sint tTimeSteps = mParameterListTimeSolver.get< moris::sint >( "TSA_Num_Time_Steps" );
    moris::real tTimeFrame = mParameterListTimeSolver.get< moris::real >( "TSA_Time_Frame" );
    moris::real tTimeIncrements = tTimeFrame / tTimeSteps;

    // init time for time slab
    Matrix< DDRMat > tTime( 2, 1, tTime_Scalar );

    for ( sint Ik = 0; Ik < tTimeSteps; Ik++ )
    {

        // log number of time steps
        MORIS_LOG_SPEC( OutputSpecifier::Iteration, (Ik+1) );
//        gLogger.log_specific(OutputSpecifier::Step, (Ik+1) );

        // set time for previous time slab
        mSolverInterface->set_previous_time( tTime );

        bool tBreaker = false;
        //Matrix< DDRMat > tTime( 2, 1, tTime_Scalar );

        if ( mNonlinearSolver->get_nonlin_solver_type() == NLA::NonlinearSolverType::ARC_LENGTH_SOLVER )
        {
            mNonlinearSolver->set_time_solver_type( this );

            // update time increment (lambda) via the arc length function
            tTime( 0, 0 ) = tTime( 1, 0 );
            tTime( 1, 0 ) = mLambdaInc;

            if ( mLambdaInc >= 1)
            {
                mLambdaInc = 1.0;
                tBreaker = true;
            }

        }
        else
        {
            tTime_Scalar += tTimeIncrements;
            tTime( 0, 0 ) = tTime( 1, 0 );
            tTime( 1, 0 ) = tTime_Scalar;
        }

        // set time for current time slab
        mSolverInterface->set_time( tTime );

        // set solution vector for previous time slab
        mSolverInterface->set_solution_vector_prev_time_step( mPrevFullVector );

        mNonlinearSolver->set_time_step_iter( Ik );

        mNonlinearSolver->solve( mFullVector );

        mPrevFullVector->vec_plus_vec( 1.0, *mFullVector, 0.0 );

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

//-------------------------------------------------------------------------------

moris::real Monolithic_Time_Solver::get_new_lambda()
{
    return mLambdaInc;
}
