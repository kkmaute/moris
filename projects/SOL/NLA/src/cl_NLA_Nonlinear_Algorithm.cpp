/*
 * cl_NLA_Nonlinear_Algorithm.cpp
 *
 *  Created on: Sep 21, 2018
 *      Author: schmidt
 */
#include <ctime>

#include "cl_NLA_Nonlinear_Algorithm.hpp"

#include "cl_DLA_Solver_Interface.hpp"

#include "cl_Vector.hpp"

#include "cl_Communication_Tools.hpp"

extern moris::Comm_Manager gMorisComm;

using namespace moris;
using namespace NLA;
using namespace dla;

void Nonlinear_Algorithm::set_linear_solver( dla::Linear_Solver * aLinSolver  )
{
    // Check if nullptr. If not delete liner solver manager
    if( mLinSolverManager != nullptr )
    {
        delete( mLinSolverManager );
    }

    // Set liner solver manager
    mLinSolverManager = aLinSolver;
}

//--------------------------------------------------------------------------------------------------------------------------
moris::real Nonlinear_Algorithm::calculate_time_needed( const clock_t aTime )
{
    moris::real tDeltaTime = (moris::real) ( clock() - aTime ) / CLOCKS_PER_SEC;

    moris::real tDeltaTimeMax   = tDeltaTime;

    max_all( tDeltaTime, tDeltaTimeMax );

    return tDeltaTimeMax;
}

//--------------------------------------------------------------------------------------------------------------------------
void Nonlinear_Algorithm::set_nonlinear_solver_manager( Nonlinear_Solver* aNonlinSolverManager )
{
    mMyNonLinSolverManager = aNonlinSolverManager;
}

//--------------------------------------------------------------------------------------------------------------------------
void Nonlinear_Algorithm::set_nonlinear_solver_parameters()
{
    // Allowable Newton solver iterations
    mParameterListNonlinearSolver.insert( "NLA_max_iter", 10 );

    // Allowable Newton solver iterations
    mParameterListNonlinearSolver.insert( "NLA_restart", 0 );

    // Allowable Newton irelative residual
//    mParameterListNonlinearSolver.insert( "NLA_rel_residual" , 1e-02 );
    mParameterListNonlinearSolver.insert( "NLA_rel_residual" , 1e-08 );

    // Desired total residual norm drop
//    mParameterListNonlinearSolver.insert( "NLA_tot_res_norm_drop" , 1e-02 );
    mParameterListNonlinearSolver.insert( "NLA_tot_res_norm_drop" , 1e-08 );

    // Desired total residual norm
//    mParameterListNonlinearSolver.insert( "NLA_tot_res_norm" , 1e-2 );
    mParameterListNonlinearSolver.insert( "NLA_tot_res_norm" , 1e-9 );

    // Maximal residual norm drop
//    mParameterListNonlinearSolver.insert( "NLA_max_res_norm_drop" , 1e-2 );
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


