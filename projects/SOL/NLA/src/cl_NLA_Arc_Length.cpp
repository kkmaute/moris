/*
 * cl_NLA_Arc_Length.cpp
 *
 *  Created on: Mar 26, 2019
 *      Author: sonne
 */
// moris includes
#include "cl_Communication_Tools.hpp"
//------------------------------------------------------------------------------
// nla includes
#include "cl_NLA_Arc_Length.hpp"
#include "cl_NLA_Convergence.hpp"
#include "cl_NLA_Nonlinear_Solver.hpp"
//------------------------------------------------------------------------------
// dla includes
#include "cl_Matrix_Vector_Factory.hpp"
#include "cl_DLA_Linear_Solver_Algorithm.hpp"
#include "cl_DLA_Solver_Interface.hpp"
#include "cl_DLA_Enums.hpp"
#include "cl_Sparse_Matrix.hpp"
#include "cl_Vector.hpp"
//------------------------------------------------------------------------------
// linalg includes
#include "cl_Matrix.hpp"
//------------------------------------------------------------------------------
// std includes
#include <ctime>
#include <cmath>
//------------------------------------------------------------------------------

using namespace moris;
using namespace NLA;
using namespace dla;

//--------------------------------------------------------------------------------------------------------------------------

Arc_Length_Solver::Arc_Length_Solver()
{
    mLinSolverManager = new dla::Linear_Solver();

    // Set default parameters in parameter list for nonlinear solver
    this->set_nonlinear_solver_parameters();
}

//--------------------------------------------------------------------------------------------------------------------------

Arc_Length_Solver::Arc_Length_Solver( dla::Linear_Solver * aLinSolver )
{
    mLinSolverManager = aLinSolver;

    // Set default parameters in parameter list for nonlinear solver
    this->set_nonlinear_solver_parameters();
}

//--------------------------------------------------------------------------------------------------------------------------

Arc_Length_Solver::~Arc_Length_Solver()
{
}

//--------------------------------------------------------------------------------------------------------------------------
// this is the work horse of the solution algorithm (where stuff actually gets done)
void Arc_Length_Solver::solver_nonlinear_system( Nonlinear_Problem * aNonlinearProblem )
{
    //------------------------------------------------------------------------------
    moris::real tB  = 0.50;
    moris::real tDa = 0.50;

    moris::real tF_tilde;
    moris::real tArcNumer,tArcDenom;

    Dist_Vector* tFext = mNonlinearProblem->get_f_ext();
    //------------------------------------------------------------------------------

    mNonlinearProblem = aNonlinearProblem;

    moris::sint tMaxIts  = mParameterListNonlinearSolver.get< moris::sint >( "NLA_max_iter" );
    //moris::real tRelRes = mParameterListNonlinearSolver.get< moris::real >( "NLA_rel_residual" );
    moris::real tRelaxation = mParameterListNonlinearSolver.get< moris::real >( "NLA_relaxation_parameter" );
    //moris::sint tRebuildIterations = mParameterListNonlinearSolver.get< moris::sint >( "NLA_num_nonlin_rebuild_iterations" );

    bool tIsConverged            = false;
    bool tRebuildJacobian        = true;
    moris::real tMaxNewTime      = 0.0;
    moris::real tMaxAssemblyTime = 0.0;
    //moris::real tErrorStatus     = 0;

    for ( sint timeStep = 1; timeStep < 10; timeStep++ )
    {   // temporary time step loop
        //------------------------------------------------------------------------------
        /*
         * 1) build K_tilde (Jacobian)
         * 2) use K_tilde to solve for d_tilde
         * 3) store K_tilde0 and d_tilde0
         */
        mNonlinearProblem->build_linearized_problem( tRebuildJacobian, timeStep ); // rebuild the linearized problem
        bool tHartBreak = false;
        this->solve_linear_system( timeStep, tHartBreak );                         // solve linearized problem
        //------------------------------------------------------------------------------

        mNonlinearProblem->get_linearized_problem()->assemble_residual_and_jacobian();  // assemble linear problem

        Sparse_Matrix* tK_tilde = mNonlinearProblem->get_full_for_consistent_tangent(); // full consistent tangent matrix
        tK_tilde->get_matrix();                                                         // get the Jacobian matrix
        Dist_Vector* tK_tildeVal = mNonlinearProblem->get_k_diag();                     // vector of diagonal values
        tK_tilde->get_diagonal( *tK_tildeVal );                                         // fill vector with diagonal values

        Dist_Vector* tD_tilde = mNonlinearProblem->get_d_tilde();                       // displacement vector
        tD_tilde = mNonlinearProblem->get_linearized_problem()->get_solver_RHS();
        tD_tilde->vec_plus_vec(1,*tFext,0);
        mNonlinearProblem->get_linearized_problem()->solve_linear_system();


        Dist_Vector* tK_tilde0Val = mNonlinearProblem->get_k0_diag();
        Dist_Vector* tD_tilde0 = mNonlinearProblem->get_d_tilde0();
        if ( timeStep==1 )
        {   // store K_tilde0 diagonal values and d_tilde0 vector
            tK_tilde0Val = tK_tildeVal;
            tD_tilde0    = tD_tilde;
        }

        uint k = 1;     // iteration counter
        // procedure 1
        //------------------------------------------------------------------------------
        sint tSize = tD_tilde->vec_local_length();
        if ( timeStep < 3 )
        {
            for ( sint i=0; i<tSize; i++ )
            {

            }
        }
        else
        {

        }
        //------------------------------------------------------------------------------
        // procedure 2

        //------------------------------------------------------------------------------

//------------------------------------------------------------------------------
        // Arc Length loop
//------------------------------------------------------------------------------
        // while ( (std::abs(tResidual/tResidualO) > tResTolerance) || ( (tF_arc-tDeltaA)/tDeltaA > tArcTol ) ) && (  iter < tMaxIts)
        for ( sint iter=1; iter<=tMaxIts; iter++ )
        {
            clock_t tArcLengthLoopStart = clock();
            clock_t tStartAssemblyTime  = clock();

            // assemble RHS and Jacobian
            if ( iter > 1 )
            {
                tRebuildJacobian = mParameterListNonlinearSolver.get< bool >( "NLA_rebuild_jacobian" );
            }
            if ( iter == 1 && mParameterListNonlinearSolver.get< sint >( "NLA_restart" ) != 0 )
            {
                sint tRestart = mParameterListNonlinearSolver.get< sint >( "NLA_restart" );
                mNonlinearProblem->build_linearized_problem( tRebuildJacobian, iter, tRestart );
            }
            else
            {
                mNonlinearProblem->build_linearized_problem( tRebuildJacobian, iter );
            }

            tMaxAssemblyTime = this->calculate_time_needed( tStartAssemblyTime );

            tHartBreak = false;
//------------------------------------------------------------------------------
            /*
             * 1) update Fint value --- N/A
             * 2) calculate residual
             * 3) determine consistent tangent (Jacobian) [K_tilde]
             * 4) calculate f_arc value
             * 5) calculate df_arc/dDelta_d
             * 6) calculate df_arc/dDelta_lambda
             */
            //------------------------------------------------------------------------------
            /*
             * --------------------------------------------
             * Static Condensation Method
             * 1) build numerator
             * 2) build denominator
             * 3) calculate delta_lambda
             * 4) use delta_lambda to calculate delta_d
             * --------------------------------------------
             */


            //------------------------------------------------------------------------------
            /*
             * Updates
             * 1) update tD_k          = tD_k + tdelta_d
             * 2) update tLambda_k     = tLambda_k + tdelta_lambda
             * 3) update tDelta_d      = tDelta_d + tdelta_d
             * 4) update tDelta_lambda = tDelta_lambda + tdelta_lambda
             *
             * 5) update Fint --- N/A
             * 6) update residual
             * 7) update f_arc
             */
            //------------------------------------------------------------------------------

            Convergence tConvergence;

            tIsConverged = tConvergence.check_for_convergence( this,
                    iter,
                    mMyNonLinSolverManager->get_ref_norm(),
                    mMyNonLinSolverManager->get_residual_norm(),
                    tMaxAssemblyTime,
                    tMaxNewTime,
                    tHartBreak ) ;

            if ( tIsConverged )
            {
                if ( tHartBreak )
                {
                    continue;
                }
                break;
            }

            // Solve linear system
            this->solve_linear_system( iter, tHartBreak );

            //PreconTime
            //SolveTime

            ( mNonlinearProblem->get_full_vector())->vec_plus_vec( -tRelaxation, *mNonlinearProblem->get_linearized_problem()
                                                   ->get_full_solver_LHS(), 1.0 );

            tMaxNewTime = this->calculate_time_needed( tArcLengthLoopStart );

        }//end iter
        //------------------------------------------------------------------------------
    }//end timeStep

}

//--------------------------------------------------------------------------------------------------------------------------
void Arc_Length_Solver::solve_linear_system( moris::sint & aIter,
                                             bool        & aHardBreak )
{
    // Solve linear system
    mLinSolverManager->solver_linear_system( mNonlinearProblem->get_linearized_problem(), aIter );
}

//--------------------------------------------------------------------------------------------------------------------------
void Arc_Length_Solver::get_full_solution( moris::Matrix< DDRMat > & LHSValues )
{
    mNonlinearProblem->get_full_vector()->extract_copy( LHSValues );
}

//--------------------------------------------------------------------------------------------------------------------------
void Arc_Length_Solver::get_solution( moris::Matrix< DDRMat > & LHSValues )
{
    mNonlinearProblem->get_full_vector()->extract_copy( LHSValues );
}
//--------------------------------------------------------------------------------------------------------------------------
void Arc_Length_Solver::extract_my_values( const moris::uint             & aNumIndices,
                                           const moris::Matrix< DDSMat > & aGlobalBlockRows,
                                           const moris::uint             & aBlockRowOffsets,
                                                 moris::Matrix< DDRMat > & LHSValues )
{
    mNonlinearProblem->get_full_vector()->extract_my_values( aNumIndices, aGlobalBlockRows, aBlockRowOffsets, LHSValues );
}
//--------------------------------------------------------------------------------------------------------------------------



