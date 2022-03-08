/*
 * cl_NLA_Newton_Solver.cpp
 *
 *  Created on: Sep 21, 2018
 *      Author: schmidt
 */
#include <ctime>

#include "cl_NLA_Newton_Solver.hpp"

#include "cl_NLA_Convergence.hpp"
#include "cl_NLA_Nonlinear_Solver.hpp"
#include "cl_NLA_Solver_Relaxation.hpp"
#include "cl_NLA_Solver_Load_Control.hpp"

#include "cl_SOL_Matrix_Vector_Factory.hpp"
#include "cl_SOL_Enums.hpp"
#include "cl_SOL_Dist_Vector.hpp"

#include "cl_DLA_Linear_Solver_Algorithm.hpp"
#include "cl_DLA_Solver_Interface.hpp"
#include "cl_DLA_Linear_Solver.hpp"
#include "cl_DLA_Linear_Problem.hpp"

#include "cl_Communication_Tools.hpp"

// Logging package
#include "cl_Logger.hpp"
#include "cl_Tracer.hpp"

using namespace moris;
using namespace NLA;
using namespace dla;

//--------------------------------------------------------------------------------------------------------------------------

Newton_Solver::Newton_Solver()
{
    // mLinSolverManager = new dla::Linear_Solver();

    // Set default parameters in parameter list for nonlinear solver
    this->set_nonlinear_solver_parameters();
}

//--------------------------------------------------------------------------------------------------------------------------

Newton_Solver::Newton_Solver( const ParameterList aParameterlist )
        : Nonlinear_Algorithm( aParameterlist )
{
}

//--------------------------------------------------------------------------------------------------------------------------

Newton_Solver::Newton_Solver( dla::Linear_Solver* aLinSolver )
        : Nonlinear_Algorithm()
{
    mLinSolverManager = aLinSolver;
}

//--------------------------------------------------------------------------------------------------------------------------

Newton_Solver::~Newton_Solver()
{
    if ( mLinSolverOwned )
    {
        delete mLinSolverManager;
    }
}

//--------------------------------------------------------------------------------------------------------------------------
void
Newton_Solver::solver_nonlinear_system( Nonlinear_Problem* aNonlinearProblem )
{
    Tracer tTracer( "NonLinearAlgorithm", "Newton", "Solve" );

    // set nonlinear system
    mNonlinearProblem = aNonlinearProblem;

    // get algorithmic parameters
    moris::sint tMaxIts                 = mParameterListNonlinearSolver.get< moris::sint >( "NLA_max_iter" );
    bool        tCombinedResJacAssembly = mParameterListNonlinearSolver.get< bool >( "NLA_combined_res_jac_assembly" );
    //        moris::sint tRebuildIterations = mParameterListNonlinearSolver.get< moris::sint >( "NLA_num_nonlin_rebuild_iterations" );

    // set relaxation strategy
    Solver_Relaxation tRelaxationStrategy( mParameterListNonlinearSolver );

    // set solver load control strategy
    Solver_Load_Control tLoadControlStrategy( mParameterListNonlinearSolver );

    // initialize flags
    bool tIsConverged     = false;
    bool tRebuildJacobian = true;

    moris::real tRelaxationParameter = 0.0;
    moris::real tMaxNewTime          = 0.0;
    moris::real tMaxAssemblyTime     = 0.0;

    // initialize load control parameter
    real tLoadFactor = tLoadControlStrategy.get_initial_load_factor();

    // Newton loop
    for ( moris::sint It = 1; It <= tMaxIts; ++It )
    {
        // log solver iteration
        MORIS_LOG_ITERATION();

        // start timing - purpose not clear
        clock_t tNewtonLoopStartTime = clock();
        clock_t tStartAssemblyTime   = clock();

        // assemble RHS and Jac
        if ( It > 1 )
        {
            tRebuildJacobian = mParameterListNonlinearSolver.get< bool >( "NLA_rebuild_jacobian" );
        }

        // For sensitivity analysis only: set current solution to LHS of linear system as residual is defined by A x - b
        if ( !mMyNonLinSolverManager->get_solver_interface()->get_is_forward_analysis() )
        {
            mNonlinearProblem->get_linearized_problem()->set_free_solver_LHS( mNonlinearProblem->get_full_vector() );
        }

        // print and store relative residual in logger
        MORIS_LOG_SPEC( "LoadFactor", tLoadFactor );

        gLogger.set_action_data(
                "NonLinearAlgorithm",
                "Newton",
                "Solve",
                "LoadFactor",
                tLoadFactor );

        // build residual and jacobian
        // restart and switch not clear
        if ( It == 1 && mParameterListNonlinearSolver.get< sint >( "NLA_restart" ) != 0 )
        {
            sint tRestart = mParameterListNonlinearSolver.get< sint >( "NLA_restart" );

            mNonlinearProblem->build_linearized_problem( tRebuildJacobian, It, tRestart );
        }
        else
        {
            mNonlinearProblem->build_linearized_problem( tRebuildJacobian, tCombinedResJacAssembly, It );
        }

        // timing not clear
        tMaxAssemblyTime = this->calculate_time_needed( tStartAssemblyTime );

        // check for convergence
        bool        tHardBreak = false;
        Convergence tConvergence;
        tIsConverged = tConvergence.check_for_convergence(
                this,
                It,
                mMyNonLinSolverManager->get_ref_norm(),
                mMyNonLinSolverManager->get_residual_norm(),
                tMaxAssemblyTime,
                tMaxNewTime,
                tHardBreak );

        // check if hard break is triggered
        if ( tHardBreak )
        {
            break;
        }

        // exit if convergence criterion is met
        if ( tIsConverged and tLoadFactor >= 1.0 )
        {
            break;
        }

        // Determine if new search direction needs to be computed and compute relaxation value
        bool tComputeSearchDirection = tRelaxationStrategy.eval(
                mMyNonLinSolverManager->get_ref_norm(),
                mMyNonLinSolverManager->get_residual_norm(),
                tRelaxationParameter );

        // Check that constant relaxation strategy with relaxation equals 1 is used in case of sensitivity analysis
        if ( !mMyNonLinSolverManager->get_solver_interface()->get_is_forward_analysis() )
        {
            MORIS_ERROR( tComputeSearchDirection && std::abs( tRelaxationParameter - 1.0 ) < 1e-12,
                    "Newton_Solver::solver_nonlinear_system - Incorrect relaxation strategy used for sensitivity analysis.\n" );
        }

        // Solve linear system to compute new search direction
        if ( tComputeSearchDirection )
        {
            // For forward analysis only: Use current solution as initial guess of linear solver
            if ( mMyNonLinSolverManager->get_solver_interface()->get_is_forward_analysis() )
            {
                mNonlinearProblem->get_linearized_problem()->set_free_solver_LHS( mNonlinearProblem->get_full_vector() );
            }

            // Solve linear system
            this->solve_linear_system( It, tHardBreak );

            // Determine load factor
            tLoadControlStrategy.eval(
                    mMyNonLinSolverManager->get_ref_norm(),
                    mMyNonLinSolverManager->get_residual_norm(),
                    tLoadFactor );
        }

        // Update solution
        ( mNonlinearProblem->get_full_vector() )->vec_plus_vec( -tRelaxationParameter, *mNonlinearProblem->get_linearized_problem()->get_full_solver_LHS(), 1.0 );

        tMaxNewTime = this->calculate_time_needed( tNewtonLoopStartTime );
    }
}

//--------------------------------------------------------------------------------------------------------------------------

void
Newton_Solver::solve_linear_system(
        moris::sint& aIter,
        bool&        aHardBreak )
{
    if ( !( mMyNonLinSolverManager->get_solver_interface()->get_is_forward_analysis() ) )
    {
        // Solve linear system
        mLinSolverManagerForAdjoint->solver_linear_system( mNonlinearProblem->get_linearized_problem(), aIter );
    }
    else
    {
        // Solve linear system
        mLinSolverManager->solver_linear_system( mNonlinearProblem->get_linearized_problem(), aIter );
    }
}

//--------------------------------------------------------------------------------------------------------------------------

void
Newton_Solver::get_full_solution( moris::Matrix< DDRMat >& LHSValues )
{
    mNonlinearProblem->get_full_vector()->extract_copy( LHSValues );
}

//--------------------------------------------------------------------------------------------------------------------------

void
Newton_Solver::get_solution( moris::Matrix< DDRMat >& LHSValues )
{
    mNonlinearProblem->get_full_vector()->extract_copy( LHSValues );
}

//--------------------------------------------------------------------------------------------------------------------------

void
Newton_Solver::extract_my_values(
        const moris::uint&                      aNumIndices,
        const moris::Matrix< DDSMat >&          aGlobalBlockRows,
        const moris::uint&                      aBlockRowOffsets,
        moris::Cell< moris::Matrix< DDRMat > >& LHSValues )
{
    mNonlinearProblem->get_full_vector()->extract_my_values( aNumIndices, aGlobalBlockRows, aBlockRowOffsets, LHSValues );
}

//--------------------------------------------------------------------------------------------------------------------------
