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

void Monolithic_Time_Solver::solve_monolithic_time_system( Vector< sol::Dist_Vector* >& aFullVector )
{
    // trace this solve
    Tracer tTracer( "Time Solver Algorithm", "Monolithic", "Solve" );

    this->finalize();

    moris::real tTime_Scalar = 0.0;

    bool tMaxTimeIterationReached = false;

    // get list of time frames
    Vector< Matrix< DDRMat > >& tTimeFrames = mMyTimeSolver->get_time_frames();

    Matrix< DDRMat > tTimeInitial( 2, 1, 0.0 );
    tTimeFrames.push_back( tTimeInitial );

    for ( uint iTimeStep = 0; iTimeStep < mTimeSteps; iTimeStep++ )
    {
        // get solvec and prev solvec index
        uint tSolVecIndex     = iTimeStep + 1;
        uint tPrevSolVecIndex = iTimeStep;

        // initialize time for time slab
        Matrix< DDRMat > tTime( 2, 1, tTime_Scalar );
        tTime_Scalar += mTimeIncrements;
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
        MORIS_LOG_SPEC( "Forward Solve Time Slab", iTimeStep + 1 );
        MORIS_LOG_SPEC( "Time Slab Start Time", tTimeFrames( tSolVecIndex )( 0, 0 ) );
        MORIS_LOG_SPEC( "Time Slab End Time", tTimeFrames( tSolVecIndex )( 1, 0 ) );

        mNonlinearSolver->set_time_step_iter( iTimeStep );

        mNonlinearSolver->solve( aFullVector( tSolVecIndex ) );

        if ( iTimeStep == mTimeSteps - 1 )
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
            MORIS_LOG_INFO( "Saving solution vector to file: %s", tStrOutputFile.c_str() );

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

        // print the solution vector
        Matrix< DDRMat > tSolVec;
        aFullVector( tSolVecIndex )->extract_copy( tSolVec );
        print( tSolVec, "Solution Vector" );

        // output state at end of time slab
        mMyTimeSolver->check_for_outputs( tTime( 1 ), tMaxTimeIterationReached, true );

        mMyTimeSolver->prepare_sol_vec_for_next_time_step();
    }
}

//-------------------------------------------------------------------------------

void Monolithic_Time_Solver::solve_implicit_DqDs( Vector< sol::Dist_Vector* >& aFullAdjointVector )
{
    // trace this solve
    Tracer tTracer( "TimeSolver", "Monolithic", "Solve" );

    // initialize time for time slab
    Vector< Matrix< DDRMat > >& tTimeFrames = mMyTimeSolver->get_time_frames();

    uint tStopTimeStepIndex = 0;    // Only consider last time step

    Vector< sol::Dist_Vector* >& tSolVec = mMyTimeSolver->get_solution_vectors();

    // Loop over all time iterations backwards
    for ( uint Ik = mTimeSteps; Ik > tStopTimeStepIndex; --Ik )
    {
        // get solvec and prev solvec index
        uint tSolVecIndex     = Ik;
        uint tPrevSolVecIndex = Ik - 1;

        // log number of time steps
        MORIS_LOG_ITERATION();

        mSolverInterface->set_solution_vector( tSolVec( tSolVecIndex ) );
        mSolverInterface->set_solution_vector_prev_time_step( tSolVec( tPrevSolVecIndex ) );

        std::cout << "need fix in Monolithic_Time_Solver::solve_implicit_DqDs \n";

        mSolverInterface->set_adjoint_solution_vector( aFullAdjointVector( 0 ) );
        mSolverInterface->set_previous_adjoint_solution_vector( aFullAdjointVector( 1 ) );

        // set time for current time slab ( since off-diagonal is computed on same time level for implicit DqDs)
        mSolverInterface->set_previous_time( tTimeFrames( tPrevSolVecIndex ) );
        mSolverInterface->set_time( tTimeFrames( tSolVecIndex ) );

        // log time slap
        MORIS_LOG_SPEC( "Adjoint Solve Time Slab", Ik );
        MORIS_LOG_SPEC( "Time Slab Start Time", tTimeFrames( tSolVecIndex )( 0, 0 ) );
        MORIS_LOG_SPEC( "Time Slab End Time", tTimeFrames( tSolVecIndex )( 1, 0 ) );

        mNonlinearSolverForSensitivityAnalysis->set_time_step_iter( Ik );

        mNonlinearSolverForSensitivityAnalysis->solve( aFullAdjointVector( 0 ) );

        Vector< enum MSI::Dof_Type > tDofTypeUnion = mMyTimeSolver->get_dof_type_union();

        mSolverInterface->set_requested_dof_types( tDofTypeUnion );

        mSolverInterface->postmultiply_implicit_dQds();

        aFullAdjointVector( 1 )->vec_plus_vec( 1.0, *( aFullAdjointVector( 0 ) ), 0.0 );

        // print the adjoint solution vector
        //        Matrix< DDRMat > tSolVec;
        //        aFullAdjointVector( 1 )->extract_copy( tSolVec );
        //        print( tSolVec, "aFullAdjointVector" );

        // output state at end of time slab
        mMyTimeSolver->check_for_outputs(
                tTimeFrames( tSolVecIndex )( 1 ),
                tPrevSolVecIndex == 0,
                false );
    }
}

//-------------------------------------------------------------------------------

void Monolithic_Time_Solver::solve( Vector< sol::Dist_Vector* >& aFullVector )
{
    // switch between forward and sensitivity analysis
    if ( mMyTimeSolver->is_forward_analysis() )
    {
        this->solve_monolithic_time_system( aFullVector );
    }
    else
    {
        this->solve_implicit_DqDs( aFullVector );
    }
}

//-------------------------------------------------------------------------------

void Monolithic_Time_Solver::set_lambda_increment( moris::real aLambdaInc )
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
