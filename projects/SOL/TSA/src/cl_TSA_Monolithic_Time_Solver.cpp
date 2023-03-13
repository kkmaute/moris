/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_TSA_Monolithic_Time_Solver.cpp
 *
 */

#include "cl_TSA_Monolithic_Time_Solver.hpp"
#include "cl_TSA_Time_Solver.hpp"
#include "cl_SOL_Dist_Vector.hpp"
#include "cl_DLA_Solver_Interface.hpp"
#include "cl_SOL_Matrix_Vector_Factory.hpp"
#include "cl_NLA_Nonlinear_Solver.hpp"
#include "cl_NLA_Nonlinear_Problem.hpp"

// Logging package
#include "cl_Logger.hpp"
#include "cl_Tracer.hpp"

#include "HDF5_Tools.hpp"

using namespace moris;
using namespace tsa;
//-------------------------------------------------------------------------------

void
Monolithic_Time_Solver::solve_monolithic_time_system( moris::Cell< sol::Dist_Vector* >& aFullVector )
{
    // trace this solve
    Tracer tTracer( "Time Solver Algorithm", "Monolithic", "Solve" );

    this->finalize();

    moris::real tTime_Scalar    = 0.0;
    sint        tTimeSteps      = mParameterListTimeSolver.get< moris::sint >( "TSA_Num_Time_Steps" );
    moris::real tTimeFrame      = mParameterListTimeSolver.get< moris::real >( "TSA_Time_Frame" );
    moris::real tTimeIncrements = tTimeFrame / tTimeSteps;

    bool tMaxTimeIterationReached = false;

    // get list of time frames
    moris::Cell< Matrix< DDRMat > >& tTimeFrames = mMyTimeSolver->get_time_frames();

    Matrix< DDRMat > tTimeInitial( 2, 1, 0.0 );
    tTimeFrames.push_back( tTimeInitial );

    for ( sint Ik = 0; Ik < tTimeSteps; Ik++ )
    {
        // get solvec and prev solvec index
        uint tSolVecIndex     = Ik + 1;
        uint tPrevSolVecIndex = Ik;

        // initialize time for time slab
        Matrix< DDRMat > tTime( 2, 1, tTime_Scalar );
        tTime_Scalar += tTimeIncrements;
        tTime( 1, 0 ) = tTime_Scalar;

        tTimeFrames.push_back( tTime );

        // log number of time steps
        MORIS_LOG_ITERATION();

        mSolverInterface->set_solution_vector( aFullVector( tSolVecIndex ) );
        mSolverInterface->set_solution_vector_prev_time_step( aFullVector( tPrevSolVecIndex ) );

        // set time for previous time slab
        mSolverInterface->set_previous_time( tTimeFrames( tPrevSolVecIndex ) );

        // set time for current time slab
        mSolverInterface->set_time( tTimeFrames( tSolVecIndex ) );

        // log time slap
        MORIS_LOG_SPEC( "Forward Solve Time Slab", Ik + 1 );
        MORIS_LOG_SPEC( "Time Slap Start Time", tTimeFrames( tSolVecIndex )( 0, 0 ) );
        MORIS_LOG_SPEC( "Time Slap End Time", tTimeFrames( tSolVecIndex )( 1, 0 ) );

        mNonlinearSolver->set_time_step_iter( Ik );

        mNonlinearSolver->solve( aFullVector( tSolVecIndex ) );

        if ( Ik == tTimeSteps - 1 )
        {
            tMaxTimeIterationReached = true;
        }

        mSolverInterface->compute_IQI();

        // if separate output of solution to file is requested
        if ( mOutputSolVecFileName.size() > 0 )
        {
            // get iteration of time solver
            uint tTSAIter = gLogger.get_iteration( "TimeSolverAlgorithm", LOGGER_ARBITRARY_DESCRIPTOR, LOGGER_ARBITRARY_DESCRIPTOR );

            // construct string from output file name
            std::string tStrOutputFile = mOutputSolVecFileName + "." + ios::stringify( tTSAIter ) + ".hdf5";

            // log/print that the initial guess is read from file
            MORIS_LOG_INFO( "Saving solution vector to file: ", tStrOutputFile.c_str() );

            // FIXME: this option doesn't work in parallel, only for serial debugging purposes
            // MORIS_ERROR( par_size() == 1, "Monolithic_Time_Solver::solve_monolithic_time_system() - "
            //         "Writing solutions to hdf5 file only possible in serial." );

            // convert distributed vector to moris mat
            Matrix< DDRMat > tSolVec;
            aFullVector( tSolVecIndex )->extract_copy( tSolVec );

            // read HDF5 file to moris matrix
            hid_t  tFileID = create_hdf5_file( tStrOutputFile );
            herr_t tStatus = 0;
            save_matrix_to_hdf5_file( tFileID, "SolVec", tSolVec, tStatus );
        }

        // input second time slap value for output
        mMyTimeSolver->check_for_outputs( tTime( 1 ), tMaxTimeIterationReached );

        mMyTimeSolver->prepare_sol_vec_for_next_time_step();
    }
}

//-------------------------------------------------------------------------------

void
Monolithic_Time_Solver::solve_implicit_DqDs( moris::Cell< sol::Dist_Vector* >& aFullAdjointVector )
{
    // trace this solve
    Tracer tTracer( "TimeSolver", "Monolithic", "Solve" );

    sint tTimeSteps = mParameterListTimeSolver.get< moris::sint >( "TSA_Num_Time_Steps" );

    // initialize time for time slab
    moris::Cell< Matrix< DDRMat > >& tTimeFrames = mMyTimeSolver->get_time_frames();

    sint tStopTimeStepIndex = 0;    // Only consider last time step

    moris::Cell< sol::Dist_Vector* >& tSolVec = mMyTimeSolver->get_solution_vectors();

    // Loop over all time iterations backwards
    for ( sint Ik = tTimeSteps; Ik > tStopTimeStepIndex; --Ik )
    {
        // get solvec and prev solvec index
        uint tSolVecIndex     = Ik;
        uint tPrevSolVecIndex = Ik - 1;

        // log number of time steps
        MORIS_LOG_ITERATION();

        mSolverInterface->set_solution_vector( tSolVec( tSolVecIndex ) );
        mSolverInterface->set_solution_vector_prev_time_step( tSolVec( tPrevSolVecIndex ) );

        mSolverInterface->set_adjoint_solution_vector( aFullAdjointVector( 0 ) );
        mSolverInterface->set_previous_adjoint_solution_vector( aFullAdjointVector( 1 ) );

        // set time for current time slab ( since off-diagonal is computed on same time level for implicit DqDs)
        mSolverInterface->set_previous_time( tTimeFrames( tPrevSolVecIndex ) );
        mSolverInterface->set_time( tTimeFrames( tSolVecIndex ) );

        // log time slap
        MORIS_LOG_SPEC( "Adjoint Solve Time Slab", Ik );
        MORIS_LOG_SPEC( "Time Slab Start Time", tTimeFrames( tSolVecIndex )( 0, 0 ) );
        MORIS_LOG_SPEC( "Time Slab End Time", tTimeFrames( tSolVecIndex )( 1, 0 ) );

        mNonlinearSolverForAdjoint->set_time_step_iter( Ik );

        mNonlinearSolverForAdjoint->solve( aFullAdjointVector( 0 ) );

        moris::Cell< enum MSI::Dof_Type > tDofTypeUnion = mMyTimeSolver->get_dof_type_union();

        mSolverInterface->set_requested_dof_types( tDofTypeUnion );

        mSolverInterface->postmultiply_implicit_dQds();

        aFullAdjointVector( 1 )->vec_plus_vec( 1.0, *( aFullAdjointVector( 0 ) ), 0.0 );
    }
}

//-------------------------------------------------------------------------------

void
Monolithic_Time_Solver::solve( moris::Cell< sol::Dist_Vector* >& aFullVector )
{
    if ( mMyTimeSolver->get_is_forward_analysis() )
    {
        this->solve_monolithic_time_system( aFullVector );
    }
    else
    {
        this->solve_implicit_DqDs( aFullVector );
    }
}

//-------------------------------------------------------------------------------

void
Monolithic_Time_Solver::set_lambda_increment( moris::real aLambdaInc )
{
    mNonlinearSolver->get_my_nonlin_problem()->set_time_value( aLambdaInc );
    mLambdaInc = aLambdaInc;
}

//-------------------------------------------------------------------------------

moris::real
Monolithic_Time_Solver::get_new_lambda()
{
    return mLambdaInc;
}
