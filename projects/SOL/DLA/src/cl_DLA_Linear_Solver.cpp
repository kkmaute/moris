/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_DLA_Linear_Solver.cpp
 *
 */
#include "cl_DLA_Linear_Solver.hpp"

#include "cl_DLA_Solver_Factory.hpp"
#include "cl_DLA_Linear_Solver_Aztec.hpp"
#include "cl_DLA_Linear_Problem.hpp"
#include "cl_SOL_Enums.hpp"
#include "cl_SOL_Dist_Vector.hpp"

#include "cl_Communication_Tools.hpp"

#include "HDF5_Tools.hpp"

#include "cl_Logger.hpp"
#include "cl_Tracer.hpp"

#include "cl_Stopwatch.hpp"    //CHR/src

using namespace moris;
using namespace dla;

Linear_Solver::Linear_Solver()
{
    // create solver factory
    Solver_Factory tSolFactory;

    // create solver object
    std::shared_ptr< Linear_Solver_Algorithm > tLinSolver = tSolFactory.create_solver( sol::SolverType::AMESOS_IMPL );

    //    tLinSolver->set_param("AZ_diagnostics") = AZ_none;
    //    tLinSolver->set_param("AZ_output") = AZ_none;

    mLinearSolverList.clear();

    mLinearSolverList.resize( 1, nullptr );

    mLinearSolverList( 0 ) = tLinSolver;

    this->set_linear_solver_manager_parameters();
}

//--------------------------------------------------------------------------------------------------

Linear_Solver::Linear_Solver( const moris::ParameterList aParameterlist )
        : mParameterListLinearSolver( aParameterlist )
{
    // create solver factory
    Solver_Factory tSolFactory;

    // create solver object
    std::shared_ptr< Linear_Solver_Algorithm > tLinSolver = tSolFactory.create_solver( sol::SolverType::AMESOS_IMPL );

    //    tLinSolver->set_param("AZ_diagnostics") = AZ_none;
    //    tLinSolver->set_param("AZ_output") = AZ_none;

    mLinearSolverList.clear();

    mLinearSolverList.resize( 1, nullptr );

    mLinearSolverList( 0 ) = tLinSolver;

    this->set_linear_solver_manager_parameters();
}

// Linear_Solver::~Linear_Solver()
//{}

//--------------------------------------------------------------------------------------------------
void
Linear_Solver::set_linear_algorithm( std::shared_ptr< Linear_Solver_Algorithm > aLinSolverAlgorithm )
{
    if ( mCallCounter == 0 )
    {
        // removes all elements from the Cell and destroy them
        mLinearSolverList.clear();

        // Resize the Cell to size = 1
        mLinearSolverList.resize( 1 );

        // Set linear solver on first entry
        mLinearSolverList( 0 ) = aLinSolverAlgorithm;
    }
    else
    {
        // set nonlinear solver on next entry
        mLinearSolverList.push_back( aLinSolverAlgorithm );
    }

    mCallCounter = mCallCounter + 1;
}

//-------------------------------------------------------------------------------------------------------
void
Linear_Solver::set_linear_algorithm(
        const moris::uint                          aListEntry,
        std::shared_ptr< Linear_Solver_Algorithm > aLinSolverAlgorithm )
{
    // Check if list is smaller than given entry
    if ( mLinearSolverList.size() >= aListEntry )
    {
        // Resize to new entry value and set nullptr on new entries
        mLinearSolverList.resize( aListEntry + 1, nullptr );
    }
    // Set linear solver on entry
    mLinearSolverList( aListEntry ) = aLinSolverAlgorithm;
}

//-------------------------------------------------------------------------------------------------------
void
Linear_Solver::solver_linear_system(
        dla::Linear_Problem* aLinearProblem,
        const moris::sint    aIter )
{
    Tracer tTracer( "LinearSolver", LOGGER_NON_SPECIFIC_ENTITY_TYPE, "Solve" );

    moris::sint tErrorStatus        = 0;
    moris::sint tMaxNumLinRestarts  = mParameterListLinearSolver.get< moris::sint >( "DLA_max_lin_solver_restarts" );
    moris::sint tTryRestartOnFailIt = 1;

    // if printing of LHS requested through input file, initialize hdf5 files here
    // and save LHS before and after solve
    if ( !this->get_LHS_output_filename().empty() )
    {
        // get current solution vector
        Matrix< DDRMat > tLHS;
        aLinearProblem->get_free_solver_LHS()->extract_copy( tLHS );

        // generate .hdf5 file
        std::string tHdf5FilePath = this->get_LHS_output_filename() + ".hdf5";
        MORIS_LOG_INFO( "Save LHS to file: ", tHdf5FilePath.c_str() );
        hid_t  tFileID = create_hdf5_file( tHdf5FilePath );
        herr_t tStatus = 0;

        // write initial guess to hdf5
        save_matrix_to_hdf5_file( tFileID, "init_guess", tLHS, tStatus );

        // solve system
        tErrorStatus = mLinearSolverList( 0 )->solve_linear_system( aLinearProblem, aIter );

        // get update after solve
        aLinearProblem->get_free_solver_LHS()->extract_copy( tLHS );

        // solve update to hdf5
        save_matrix_to_hdf5_file( tFileID, "update", tLHS, tStatus );

        // close hdf5
        close_hdf5_file( tFileID );
    }
    // otherwise, don't save file and simply solve system
    else
    {
        // solve system
        tErrorStatus = mLinearSolverList( 0 )->solve_linear_system( aLinearProblem, aIter );
    }

    // Restart the linear solver using the current solution as an initial guess if the previous linear solve failed
    while ( tErrorStatus != 0 && tTryRestartOnFailIt <= tMaxNumLinRestarts && (moris::sint)mLinearSolverList.size() <= tMaxNumLinRestarts )
    {
        if ( par_rank() == 0 )
        {
            // Compute current solution vector norm
            Cell< moris::real > tSolVecNorm = aLinearProblem->get_free_solver_LHS()->vec_norm2();

            MORIS_LOG( " ... Previous linear solve failed. Trying restart %i of %i, using current solution with SolVecNorm = %5.15e as an initial guess. ",
                    tTryRestartOnFailIt,
                    tMaxNumLinRestarts,
                    tSolVecNorm( 0 ) );
        }

        // Re-solve scaled linear system with current solution as an initial guess
        tErrorStatus = mLinearSolverList( tTryRestartOnFailIt )->solve_linear_system( aLinearProblem, aIter );

        // Iterate TryRestartOnFailIt counter
        tTryRestartOnFailIt = tTryRestartOnFailIt + 1;
    }

    if ( ( tErrorStatus != 0 ) )
    {
        if ( par_rank() == 0 )
        {
            MORIS_LOG( " " );
            MORIS_LOG( "Linear Solver status absolute value = %i", tErrorStatus );
            MORIS_LOG( "Linear Solver did not exit with status 0!" );
        }
    }
}

//--------------------------------------------------------------------------------------------------------------------------
void
Linear_Solver::set_linear_solver_manager_parameters()
{
    // Maximal number of linear solver restarts on fail
    mParameterListLinearSolver.insert( "DLA_max_lin_solver_restarts", 0 );

    // Maximal number of linear solver restarts on fail
    mParameterListLinearSolver.insert( "DLA_hard_break", true );

    // Determines if lin solve should restart on fail
    mParameterListLinearSolver.insert( "DLA_rebuild_lin_solver_on_fail", false );
}
