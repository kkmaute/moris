/*
 * cl_NLA_Nonlinear_Solver.cpp
 *
 *  Created on: Sep 21, 2018
 *      Author: schmidt
 */
#include "cl_NLA_Nonlinear_Solver.hpp"

#include <ctime>

#include "cl_DLA_Linear_Solver.hpp"
#include "cl_DLA_Solver_Interface.hpp"

//#include "cl_DLA_Linear_Solver_Manager.hpp"
#include "cl_Vector.hpp"

#include "cl_Communication_Tools.hpp"

using namespace moris;
using namespace NLA;
using namespace dla;

void Nonlinear_Solver::set_linear_solver_manager( dla::Linear_Solver_Manager * aLinSolverManager  )
{
    // Check if nullptr. If not delete liner solver manager
    if( mLinSolverManager != nullptr )
    {
        delete( mLinSolverManager );
    }

    // Set liner solver manager
    mLinSolverManager = aLinSolverManager;
}

//--------------------------------------------------------------------------------------------------------------------------
moris::real Nonlinear_Solver::calculate_time_needed( const clock_t aTime )
{
    moris::real tDeltaTime = (moris::real) ( clock() - aTime ) / CLOCKS_PER_SEC;

    moris::real tDeltaTimeMax   = tDeltaTime;

    MPI_Allreduce( &tDeltaTime, &tDeltaTimeMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );

    return tDeltaTimeMax;
}

//--------------------------------------------------------------------------------------------------------------------------
void Nonlinear_Solver::set_nonlinear_solver_manager( Nonlinear_Solver_Manager* aNonlinSolverManager )
{
    mMyNonLinSolverManager = aNonlinSolverManager;
}

//--------------------------------------------------------------------------------------------------------------------------
void Nonlinear_Solver::set_nonlinear_solver_parameters()
{
    // Allowable Newton solver iterations
    mParameterListNonlinearSolver.insert( "NLA_max_iter", 10 );

    // Allowable Newton solver iterations
    mParameterListNonlinearSolver.insert( "NLA_restart", 0 );

    // Allowable Newton irelative residual
    mParameterListNonlinearSolver.insert( "NLA_rel_residual" , 1e-08 );

    // Desired total residual norm drop
    mParameterListNonlinearSolver.insert( "NLA_tot_res_norm_drop" , 1e-08 );

    // Desired total residual norm
    mParameterListNonlinearSolver.insert( "NLA_tot_res_norm" , 1e-9 );

    // Maximal residual norm drop
    mParameterListNonlinearSolver.insert( "NLA_max_res_norm_drop" , 1e-6 );

    // Maximal number of linear solver restarts on fail
    mParameterListNonlinearSolver.insert( "NLA_max_lin_solver_restarts" , 0 );

    // Maximal number of linear solver restarts on fail
    mParameterListNonlinearSolver.insert( "NLA_relaxation_parameter" , 1.0 );

    // Maximal number of linear solver restarts on fail
    mParameterListNonlinearSolver.insert( "NLA_hard_break" , false );

    // Determines if lin solve should restart on fail
    mParameterListNonlinearSolver.insert( "NLA_rebuild_lin_solv_on_fail" , false );

    // Determines if lin solve should restart on fail
    mParameterListNonlinearSolver.insert( "NLA_rebuild_jacobian" , true );

    // Determines if newton should restart on fail
    mParameterListNonlinearSolver.insert( "NLA_rebuild_nonlin_solv_on_fail" , false );

    // Specifying the number of newton retries
    mParameterListNonlinearSolver.insert( "NLA_num_nonlin_rebuild_iterations" , 1 );

    // Determines relaxation multiplier
    mParameterListNonlinearSolver.insert( "NLA_relaxation_multiplier_on_fail" , 0.5 );

    // Determines newton maxits multiplier
    mParameterListNonlinearSolver.insert( "NLA_maxits_multiplier_on_fail" , 2 );
}


