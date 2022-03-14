/*
 * cl_NLA_NLBGS.cpp
 *
 *  Created on: Jan 18, 2018
 *      Author: schmidt
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

NonLinBlockGaussSeidel::NonLinBlockGaussSeidel()
{
    // Set default parameters in parameter list for nonlinear solver
    this->set_nonlinear_solver_parameters();
}

//--------------------------------------------------------------------------------------------------------------------------

NonLinBlockGaussSeidel::NonLinBlockGaussSeidel( const ParameterList aParameterlist )
        : Nonlinear_Algorithm( aParameterlist )
{
}

//--------------------------------------------------------------------------------------------------------------------------

NonLinBlockGaussSeidel::~NonLinBlockGaussSeidel()
{
}

//--------------------------------------------------------------------------------------------------------------------------

void
NonLinBlockGaussSeidel::solver_nonlinear_system( Nonlinear_Problem* aNonlinearProblem )
{
    Tracer tTracer( "NonLinearAlgorithm", "NLBGS", "Solve" );

    moris::sint tMaxIts           = mParameterListNonlinearSolver.get< moris::sint >( "NLA_max_iter" );
    moris::uint tNonLinSysStartIt = 0;
    moris::uint tNumNonLinSystems = mMyNonLinSolverManager->get_dof_type_list().size();

    // set solver load control strategy
    Solver_Load_Control tLoadControlStrategy( mParameterListNonlinearSolver );

    // set pseudo time control strategy
    Solver_Pseudo_Time_Control tPseudoTimeControl(
            mParameterListNonlinearSolver,
            aNonlinearProblem->get_full_vector(),
            mMyNonLinSolverManager->get_solver_interface(),
            mMyNonLinSolverManager->get_solver_warehouse() );

    // initialize load control parameter
    real tLoadFactor = tLoadControlStrategy.get_initial_load_factor();

    // initialize pseudo time step
    real tPseudoTimeStep;
    bool tTimeStepIsConverged = tPseudoTimeControl.get_initial_step_size( tPseudoTimeStep );

    // NLBGS loop
    for ( sint It = 1; It <= tMaxIts; ++It )
    {
        // log iterations
        MORIS_LOG_ITERATION();

        moris::real tMaxNewTime = 0.0;

        // get_nonlinear_problem()
        clock_t tNewtonLoopStartTime = clock();

        // print and store pseudo time step in logger
        MORIS_LOG_SPEC( "PseudoTimeStep", tPseudoTimeStep );

        gLogger.set_action_data(
                "NonLinearAlgorithm",
                "NLBGS",
                "Solve",
                "PseudoTimeStep",
                tPseudoTimeStep );

        // print and store load factor in logger
        MORIS_LOG_SPEC( "LoadFactor", tLoadFactor );

        gLogger.set_action_data(
                "NonLinearAlgorithm",
                "NLBGS",
                "Solve",
                "LoadFactor",
                tLoadFactor );

        // switch between forward and backward system
        if ( mMyNonLinSolverManager->get_solver_interface()->get_is_forward_analysis() )
        {
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
            // Loop over all non-linear systems backwards
            for ( sint Ik = tNumNonLinSystems; Ik > (sint)tNonLinSysStartIt; Ik-- )
            {
                // Log/print which NL system is being solved
                MORIS_LOG_SPEC( "Sensitivity Analysis Nonlinear System", Ik );

                // Get subsolver
                Nonlinear_Solver* tSubSolver = mMyNonLinSolverManager->get_sub_nonlinear_solver( Ik - 1 );

                MORIS_ERROR( tSubSolver != nullptr,
                        "NonLinBlockGaussSeidel::solver_nonlinear_system - sensitivty analysis: sub-solver %d not defined.",
                        Ik - 1 );

                tSubSolver->solve( aNonlinearProblem->get_full_vector() );
            }
        }

        tMaxNewTime     = this->calculate_time_needed( tNewtonLoopStartTime );
        bool tHartBreak = false;

        // compute residual norms using residual from sub-solvers
        this->compute_norms( It );

        // check for convergence
        Convergence tConvergence;

        bool tIsConverged = tConvergence.check_for_convergence(
                this,
                It,
                mMyNonLinSolverManager->get_ref_norm(),
                mMyNonLinSolverManager->get_residual_norm(),
                tMaxNewTime,
                tHartBreak );

        // check if hard break is triggered
        if ( tHartBreak )
        {
            break;
        }

        // exit if convergence criterion is met
        if ( tIsConverged and tLoadFactor >= 1.0 and tTimeStepIsConverged )
        {
            break;
        }

        // compute new time step size and check for convergence of time stepping
        tTimeStepIsConverged = tPseudoTimeControl.compute_time_step_size(
                mMyNonLinSolverManager->get_ref_norm(),
                mMyNonLinSolverManager->get_residual_norm(),
                aNonlinearProblem->get_full_vector(),
                tPseudoTimeStep );

        // Determine load factor
        tLoadControlStrategy.eval(
                mMyNonLinSolverManager->get_ref_norm(),
                mMyNonLinSolverManager->get_residual_norm(),
                tLoadFactor );
    }    // end loop for NLBGS iterations
}

//--------------------------------------------------------------------------------------------------------------------------

void
NonLinBlockGaussSeidel::solve_linear_system(
        moris::sint& aIter,
        bool&        aHardBreak )
{
    // Solve linear system
}

//--------------------------------------------------------------------------------------------------------------------------

void
NonLinBlockGaussSeidel::compute_norms( const moris::sint aIter )
{
    moris::uint tNumNonLinSystems = mMyNonLinSolverManager->get_dof_type_list().size();
    moris::uint tNonLinSysStartIt = 0;

    if ( aIter == 1 )
    {
        moris::real tSqrtRefNorm = 0.0;

        // Loop over all non-linear systems
        for ( uint Ik = tNonLinSysStartIt; Ik < tNumNonLinSystems; Ik++ )
        {
            moris::real tSubSolverRefNorm = mMyNonLinSolverManager->get_sub_nonlinear_solver( Ik )->get_ref_norm();

            tSqrtRefNorm = tSqrtRefNorm + std::pow( tSubSolverRefNorm, 2 );
        }

        mMyNonLinSolverManager->get_ref_norm() = std::sqrt( tSqrtRefNorm );
    }

    moris::real tSqrtResNorm = 0.0;

    // Loop over all non-linear systems
    for ( uint Ik = tNonLinSysStartIt; Ik < tNumNonLinSystems; Ik++ )
    {
        moris::real tSubSolverResNorm = mMyNonLinSolverManager->get_sub_nonlinear_solver( Ik )->get_residual_norm();

        tSqrtResNorm = tSqrtResNorm + std::pow( tSubSolverResNorm, 2 );
    }

    mMyNonLinSolverManager->get_residual_norm() = std::sqrt( tSqrtResNorm );
}
