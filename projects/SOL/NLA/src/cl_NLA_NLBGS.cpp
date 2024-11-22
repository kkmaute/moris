/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_NLA_NLBGS.cpp
 *
 */
#include <ctime>

#include "cl_NLA_NLBGS.hpp"

#include "cl_NLA_Convergence.hpp"
#include "cl_NLA_Nonlinear_Solver.hpp"
#include "cl_NLA_Solver_Load_Control.hpp"
#include "cl_NLA_Solver_Pseudo_Time_Control.hpp"

#include "cl_DLA_Linear_Solver_Algorithm.hpp"
#include "cl_DLA_Solver_Interface.hpp"

#include "cl_SOL_Enums.hpp"
#include "cl_SOL_Matrix_Vector_Factory.hpp"
#include "cl_SOL_Dist_Vector.hpp"
#include "cl_SOL_Warehouse.hpp"

#include "cl_Communication_Tools.hpp"

#include "cl_Logger.hpp"
#include "cl_Tracer.hpp"

using namespace moris;
using namespace NLA;
using namespace dla;

//--------------------------------------------------------------------------------------------------------------------------

NonLinBlockGaussSeidel::NonLinBlockGaussSeidel( const Parameter_List& aParameterList )
        : Nonlinear_Algorithm( aParameterList )
{
}

//--------------------------------------------------------------------------------------------------------------------------

NonLinBlockGaussSeidel::~NonLinBlockGaussSeidel()
{
}

//--------------------------------------------------------------------------------------------------------------------------

void NonLinBlockGaussSeidel::solver_nonlinear_system( Nonlinear_Problem* aNonlinearProblem )
{
    Tracer tTracer( "NonLinearAlgorithm", "NLBGS", "Solve" );

    sint tMaxIts           = mParameterListNonlinearSolver.get< sint >( "NLA_max_iter" );
    sint tRefIts           = mParameterListNonlinearSolver.get< sint >( "NLA_ref_iter" );
    uint tNonLinSysStartIt = 0;
    uint tNumNonLinSystems = mMyNonLinSolverManager->get_dof_type_list().size();

    // set solver load control strategy
    Solver_Load_Control tLoadControlStrategy( mParameterListNonlinearSolver );

    // set pseudo time control strategy
    Solver_Pseudo_Time_Control tPseudoTimeControl(
            mParameterListNonlinearSolver,
            aNonlinearProblem->get_full_vector(),
            mMyNonLinSolverManager );

    // initialize load control parameter
    real tLoadFactor = tLoadControlStrategy.get_initial_load_factor();

    // initialize pseudo time step
    real tPseudoTimeStep;
    real tPseudoTotalTime = 0.0;
    real tRelStaticRes    = 1.0;

    // initialize convergence monitoring
    bool tTimeStepIsConverged = tPseudoTimeControl.get_initial_step_size( tPseudoTimeStep );

    Convergence tConvergence( tRefIts );

    // NLBGS loop
    for ( sint It = 1; It <= tMaxIts; ++It )
    {
        // log iterations
        MORIS_LOG_ITERATION();

        // print and store pseudo time step, total time step, and relative static residual in logger
        MORIS_LOG_SPEC( "PseudoTimeStep", tPseudoTimeStep );
        MORIS_LOG_SPEC( "PseudoTotalTime", tPseudoTotalTime );
        MORIS_LOG_SPEC( "RelStaticResidual", tRelStaticRes );

        gLogger.set_action_data(
                "NonLinearAlgorithm",
                "NLBGS",
                "Solve",
                "PseudoTimeStep",
                tPseudoTimeStep );

        gLogger.set_action_data(
                "NonLinearAlgorithm",
                "NLBGS",
                "Solve",
                "PseudoTotalTime",
                tPseudoTotalTime );

        gLogger.set_action_data(
                "NonLinearAlgorithm",
                "NLBGS",
                "Solve",
                "RelativeStaticResidual",
                tRelStaticRes );

        // print and store load factor in logger
        MORIS_LOG_SPEC( "LoadFactor", tLoadFactor );

        gLogger.set_action_data(
                "NonLinearAlgorithm",
                "NLBGS",
                "Solve",
                "LoadFactor",
                tLoadFactor );

        // switch between forward and backward system
        if ( mMyNonLinSolverManager->get_solver_interface()->is_forward_analysis() )
        {
            // Pause if needed
            mForwardPauseFunction();
            
            // Loop over all non-linear systems
            for ( uint Ik = tNonLinSysStartIt; Ik < tNumNonLinSystems; Ik++ )
            {
                // Log/print which NL system is being solved
                MORIS_LOG_SPEC( "Forward Analysis Nonlinear System", Ik );

                // Set the nonlinear system index
                gLogger.set_iteration( "NonLinearSolver", LOGGER_ARBITRARY_DESCRIPTOR, LOGGER_ARBITRARY_DESCRIPTOR, Ik );

                // Get subsolver
                Nonlinear_Solver* tSubSolver = mMyNonLinSolverManager->get_sub_nonlinear_solver( Ik );

                MORIS_ERROR( tSubSolver != nullptr,
                        "NonLinBlockGaussSeidel::solver_nonlinear_system - forward analysis: sub-solver %d not defined.",
                        Ik );

                // Solver sub-system
                tSubSolver->solve( aNonlinearProblem->get_full_vector() );
            }
        }
        else
        {
            // Pause if needed
            mSensitivityPauseFunction();

            //            std::cout << "need fix in NonLinBlockGaussSeidel::solver_nonlinear_system\n";

            // Loop over all non-linear systems backwards
            for ( sint Ik = tNumNonLinSystems; Ik > (sint)tNonLinSysStartIt; Ik-- )
            {
                // Log/print which NL system is being solved
                MORIS_LOG_SPEC( "Sensitivity Analysis Nonlinear System", Ik );

                // Get subsolver
                Nonlinear_Solver* tSubSolver = mMyNonLinSolverManager->get_sub_nonlinear_solver( Ik - 1 );

                MORIS_ERROR( tSubSolver != nullptr,
                        "NonLinBlockGaussSeidel::solver_nonlinear_system - sensitivity analysis: sub-solver %d not defined.",
                        Ik - 1 );

                tSubSolver->solve( aNonlinearProblem->get_full_vector() );
            }
        }

        bool tHartBreak = false;

        // compute residual norms using residual from sub-solvers
        this->compute_norms( It, tRefIts );

        // check for convergence
        bool tIsConverged = tConvergence.check_for_convergence(
                this,
                It,
                tHartBreak );

        // exit if convergence criterion is met
        if ( tIsConverged and tLoadFactor >= 1.0 and tTimeStepIsConverged )
        {
            MORIS_LOG_INFO( "Number of Iterations (Convergence): %d", It );

            break;
        }

        // check if hard break is triggered
        if ( tHartBreak or It == tMaxIts )
        {
            MORIS_LOG_INFO( "Number of Iterations (Hard Stop): %d", It );

            break;
        }

        // compute new time step size and check for convergence of time stepping
        tTimeStepIsConverged = tPseudoTimeControl.compute_time_step_size(
                mMyNonLinSolverManager,
                aNonlinearProblem->get_full_vector(),
                tPseudoTimeStep,
                tPseudoTotalTime,
                tRelStaticRes,
                It );

        // Determine load factor
        tLoadControlStrategy.eval(
                It,
                mMyNonLinSolverManager,
                tLoadFactor );
    }    // end loop for NLBGS iterations
}

//--------------------------------------------------------------------------------------------------------------------------

void NonLinBlockGaussSeidel::solve_linear_system(
        sint& aIter,
        bool& aHardBreak )
{
    // Solve linear system
}

//--------------------------------------------------------------------------------------------------------------------------

void NonLinBlockGaussSeidel::compute_norms(
        const sint aIter,
        const sint aRefIts )
{
    uint tNumNonLinSystems = mMyNonLinSolverManager->get_dof_type_list().size();
    uint tNonLinSysStartIt = 0;

    // check whether static residual need to be evaluated
    bool tComputeStaticResidual = mMyNonLinSolverManager->get_compute_static_residual_flag();

    // compute reference residual norms
    if ( aIter <= aRefIts )
    {
        real tSqrtRefNorm       = 0.0;
        real tSqrtStaticRefNorm = 0.0;

        // Loop over all non-linear systems
        for ( uint Ik = tNonLinSysStartIt; Ik < tNumNonLinSystems; Ik++ )
        {
            real tSubSolverRefNorm = mMyNonLinSolverManager->get_sub_nonlinear_solver( Ik )->get_ref_norm();

            tSqrtRefNorm += std::pow( tSubSolverRefNorm, 2 );

            if ( tComputeStaticResidual )
            {
                real tSubSolverStaticRefNorm = mMyNonLinSolverManager->get_sub_nonlinear_solver( Ik )->get_static_ref_norm();

                tSqrtStaticRefNorm += std::pow( tSubSolverStaticRefNorm, 2 );
            }
        }

        // store reference residual norms
        mMyNonLinSolverManager->set_ref_norm( std::sqrt( tSqrtRefNorm ) );

        if ( tComputeStaticResidual )
        {
            mMyNonLinSolverManager->set_static_ref_norm( std::sqrt( tSqrtStaticRefNorm ) );
        }
    }

    // compute current residual norms
    real tSqrtResNorm       = 0.0;
    real tSqrtStaticResNorm = 0.0;

    // Loop over all non-linear systems
    for ( uint Ik = tNonLinSysStartIt; Ik < tNumNonLinSystems; Ik++ )
    {
        real tSubSolverResNorm = mMyNonLinSolverManager->get_sub_nonlinear_solver( Ik )->get_residual_norm();

        tSqrtResNorm += std::pow( tSubSolverResNorm, 2 );

        if ( tComputeStaticResidual )
        {
            // get static residual norm
            real tSubSolverStaticResNorm = mMyNonLinSolverManager->get_sub_nonlinear_solver( Ik )->get_static_residual_norm();

            // check that static residual is valid
            MORIS_ASSERT( tSubSolverStaticResNorm >= 0.0,
                    "NonLinBlockGaussSeidel::compute_norms - invalid static residual." );

            // square static residual
            tSqrtStaticResNorm += std::pow( tSubSolverStaticResNorm, 2 );
        }
    }

    // store current residual norms
    mMyNonLinSolverManager->set_residual_norm( std::sqrt( tSqrtResNorm ) );

    if ( tComputeStaticResidual )
    {
        mMyNonLinSolverManager->set_static_residual_norm( std::sqrt( tSqrtStaticResNorm ) );
    }

    // compute maximum relative number of iterations
    real tRelNumIterations = 0.0;

    // Loop over all non-linear systems
    for ( uint Ik = tNonLinSysStartIt; Ik < tNumNonLinSystems; Ik++ )
    {
        // get relative number of iterations from subsolver
        real tSubSolverRelNumIter = mMyNonLinSolverManager->get_sub_nonlinear_solver( Ik )->get_relative_number_iterations();

        // compute maximum
        tRelNumIterations = std::max( tRelNumIterations, tSubSolverRelNumIter );
    }

    // store maximum relative number of iterations
    mMyNonLinSolverManager->set_relative_number_iterations( tRelNumIterations );
}
