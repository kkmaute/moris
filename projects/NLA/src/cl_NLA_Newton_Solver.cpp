/*
 * cl_NLA_Newton_Solver.cpp
 *
 *  Created on: Sep 21, 2018
 *      Author: schmidt
 */
#include <ctime>

#include "cl_NLA_Newton_Solver.hpp"

#include "cl_NLA_Convergence.hpp"

#include "cl_Matrix_Vector_Factory.hpp"
#include "cl_DLA_Linear_Solver.hpp"
#include "cl_DLA_Solver_Interface.hpp"
#include "cl_DLA_Enums.hpp"
#include "cl_Vector.hpp"

#include "cl_Communication_Tools.hpp"

using namespace moris;
using namespace NLA;
using namespace dla;

    Newton_Solver::Newton_Solver()
    {
        // Set default parameters in parameter list for nonlinear solver
        this->set_nonlinear_solver_parameters();
    }

    Newton_Solver::Newton_Solver( Solver_Interface * aSolverInterface ) : Nonlinear_Solver( aSolverInterface )
    {
        // Set default parameters in parameter list for nonlinear solver
        this->set_nonlinear_solver_parameters();
    }

    Newton_Solver::~Newton_Solver()
    {
    }

//--------------------------------------------------------------------------------------------------------------------------
    void Newton_Solver::solver_nonlinear_system( Nonlinear_Problem * aNonlinearProblem )
    {
        mNonlinearProblem = aNonlinearProblem;

        moris::sint tMaxIts  = mParameterListNonlinearSolver.get< moris::sint >( "NLA_max_iter" );
        //moris::real tRelRes = mParameterListNonlinearSolver.get< moris::real >( "NLA_rel_residual" );
        moris::real tRelaxation = mParameterListNonlinearSolver.get< moris::real >( "NLA_relaxation_parameter" );
//        moris::sint tRebuildIterations = mParameterListNonlinearSolver.get< moris::sint >( "NLA_num_nonlin_rebuild_iterations" );

        bool tIsConverged            = false;
        moris::real refNorm          = 0.0;
        moris::real tMaxNewTime      = 0.0;
        moris::real tMaxAssemblyTime = 0.0;
        //moris::real tErrorStatus     = 0;

        // Newton retry loop
//        for ( moris::sint Ir = 0; Ir < tRebuildIterations; ++ Ir)
//        {
//            // Reset the old state
//            if ( Ir > 0 )
//            {
//                tRelaxation *= mParameterListNonlinearSolver.get< moris::real >( "NLA_relaxation_multiplier_on_fail" );
//                tMaxIts      = (moris::sint) ( tMaxIts * mParameterListNonlinearSolver.get< moris::sint >( "NLA_maxits_multiplier_on_fail" ));
//
//                mVectorFullSol->vec_plus_vec( 1.0, *mPrevVectorFreeSol, 0.0 );
//
//                if ( par_rank() == 0 )
//                {
//                    fprintf( stdout,"\n ... Previous nonlinear solve failed ... Retrying with: max nonlinear iterations = %i; relaxation = %f\n", tMaxIts, tRelaxation );
//                }
//            }

            // Newton loop
            for ( moris::sint It = 1; It <= tMaxIts; ++It )
            {
                  //get_nonlinear_problem()
                clock_t tNewtonLoopStartTime = clock();
                clock_t tStartAssemblyTime = clock();

                // assemble RHS and Jac
                mNonlinearProblem->build_linearized_problem();

                tMaxAssemblyTime = get_time_needed( tStartAssemblyTime );

                bool tHartBreak = false;
                Convergence tConvergence;
                tIsConverged = tConvergence.check_for_convergence( this, It, refNorm, tMaxAssemblyTime, tHartBreak );

                if ( tIsConverged )
                {
                    if ( tHartBreak )
                    {
                        continue;
                    }
                    break;
                }

                // Solve linear system
                this->solve_linear_system( It, tHartBreak );

//                if ( tHartBreak )
//                 {
//                     continue;
//                 }
//                 break;

                //PreconTime
                //SolveTime

                (mNonlinearProblem->get_full_vector())->vec_plus_vec( tRelaxation, *mNonlinearProblem->get_linearized_problem()->get_full_solver_LHS(), 1.0 );

                // Update the SolVecNorm
                // solNorm = mVectorFreeSol.Norm2();

                //mVectorFullSol->import_local_to_global( *mVectorFreeSol );

                tMaxNewTime = get_time_needed( tNewtonLoopStartTime );

                std::cout<<tMaxNewTime<<std::endl;
            }
//
//            // Check if Newton converged. Break retry loop if its converged
//            if ( tIsConverged )
//            {
//                break;
//            }
//        }
    }

//--------------------------------------------------------------------------------------------------------------------------
    void Newton_Solver::solver_nonlinear_system()
    {
        moris::sint tMaxIts  = mParameterListNonlinearSolver.get< moris::sint >( "NLA_max_iter" );
        //moris::real tRelRes = mParameterListNonlinearSolver.get< moris::real >( "NLA_rel_residual" );
        moris::real tRelaxation = mParameterListNonlinearSolver.get< moris::real >( "NLA_relaxation_parameter" );
//        moris::sint tRebuildIterations = mParameterListNonlinearSolver.get< moris::sint >( "NLA_num_nonlin_rebuild_iterations" );

        bool tIsConverged            = false;
        moris::real refNorm          = 0.0;
        moris::real tMaxNewTime      = 0.0;
        moris::real tMaxAssemblyTime = 0.0;
        //moris::real tErrorStatus     = 0;

        // Newton retry loop
//        for ( moris::sint Ir = 0; Ir < tRebuildIterations; ++ Ir)
//        {
//            // Reset the old state
//            if ( Ir > 0 )
//            {
//                tRelaxation *= mParameterListNonlinearSolver.get< moris::real >( "NLA_relaxation_multiplier_on_fail" );
//                tMaxIts      = (moris::sint) ( tMaxIts * mParameterListNonlinearSolver.get< moris::sint >( "NLA_maxits_multiplier_on_fail" ));
//
//                mVectorFullSol->vec_plus_vec( 1.0, *mPrevVectorFreeSol, 0.0 );
//
//                if ( par_rank() == 0 )
//                {
//                    fprintf( stdout,"\n ... Previous nonlinear solve failed ... Retrying with: max nonlinear iterations = %i; relaxation = %f\n", tMaxIts, tRelaxation );
//                }
//            }

            // Newton loop
            for ( moris::sint It = 1; It <= tMaxIts; ++It )
            {
                  //get_nonlinear_problem()
                clock_t tNewtonLoopStartTime = clock();
                clock_t tStartAssemblyTime = clock();

                // assemble RHS and Jac
                mNonlinearProblem->build_linearized_problem();

                tMaxAssemblyTime = get_time_needed( tStartAssemblyTime );

                bool tHartBreak = false;
                Convergence tConvergence;
                tIsConverged = tConvergence.check_for_convergence( this, It, refNorm, tMaxAssemblyTime, tHartBreak );

                if ( tIsConverged )
                {
                    if ( tHartBreak )
                    {
                        continue;
                    }
                    break;
                }

                // Solve linear system
                this->solve_linear_system( It, tHartBreak );

//                if ( tHartBreak )
//                 {
//                     continue;
//                 }
//                 break;

                //PreconTime
                //SolveTime

                (mNonlinearProblem->get_full_vector())->vec_plus_vec( tRelaxation, *mNonlinearProblem->get_linearized_problem()->get_full_solver_LHS(), 1.0 );

                // Update the SolVecNorm
                // solNorm = mVectorFreeSol.Norm2();

                //mVectorFullSol->import_local_to_global( *mVectorFreeSol );

                tMaxNewTime = get_time_needed( tNewtonLoopStartTime );

                std::cout<<tMaxNewTime<<std::endl;
            }
//
//            // Check if Newton converged. Break retry loop if its converged
//            if ( tIsConverged )
//            {
//                break;
//            }
//        }
    }

//--------------------------------------------------------------------------------------------------------------------------
    void Newton_Solver::solve_linear_system( moris::sint & aIter,
                                             bool        & aHardBreak )
    {
        moris::sint tErrorStatus = 0;
        moris::sint tMaxNumLinRestarts  = mParameterListNonlinearSolver.get< moris::sint >( "NLA_max_lin_solver_restarts" );
        moris::sint tTryRestartOnFailIt = 1;

        // Solve linear system
        ( mLinSolverManager.get_solver_list())( 0 )->set_linear_problem( mNonlinearProblem->get_linearized_problem() );

        tErrorStatus = ( mLinSolverManager.get_solver_list() )( 0 )->solve_linear_system();


        // Restart the linear solver using the current solution as an initial guess if the previous linear solve failed
        while ( tErrorStatus !=0 && tTryRestartOnFailIt <= tMaxNumLinRestarts && ( moris::sint )mLinSolverManager.get_solver_list().size() <= tMaxNumLinRestarts )
        {
            if ( par_rank() == 0 )
            {
                // Compute current solution vector norm
                moris::real tSolVecNorm = mNonlinearProblem->get_linearized_problem()->get_free_solver_LHS()->vec_norm2();

                fprintf( stdout, " ... Previous linear solve failed. Trying restart %i of %i, using current solution with SolVecNorm = %5.15e as an initial guess. \n",
                                       tTryRestartOnFailIt, tMaxNumLinRestarts, tSolVecNorm);
            }

            mLinSolverManager.get_solver_list()( tTryRestartOnFailIt )->set_linear_problem( mNonlinearProblem->get_linearized_problem() );

            // Re-solve scaled linear system with current solution as an initial guess
            tErrorStatus = ( mLinSolverManager.get_solver_list() )( tTryRestartOnFailIt )->solve_linear_system();

            // Iterate TryRestartOnFailIt counter
            tTryRestartOnFailIt = tTryRestartOnFailIt + 1;
        }

        if ( ( tErrorStatus != 0 && mParameterListNonlinearSolver.get< bool >( "NLA_hard_break" ) ) && !mParameterListNonlinearSolver.get< bool >( "NLA_rebuild_on_fail" ) )
        {
            if( par_rank() == 0 )
            {
                fprintf( stdout, "\n Linear Solver status absolute value = %i\n", tErrorStatus );
                MORIS_ERROR( false, "Linear Solver did not exit with status 0!\n" );
            }
        }
        else if( ( tErrorStatus != 0 ) )
        {
            if( par_rank() == 0)
            {
                fprintf( stdout, "\n Linear Solver status absolute value = %i\n", tErrorStatus );
                fprintf( stdout, "Linear Solver did not exit with status 0!\n" );
            }
        }

        if ( tErrorStatus != 0 && mParameterListNonlinearSolver.get< bool >( "NLA_hard_break" ) )
        {
            aIter = mParameterListNonlinearSolver.get< moris::sint >( "NLA_max_iter" );
            aHardBreak = true;
        }
    }

//--------------------------------------------------------------------------------------------------------------------------
    moris::real Newton_Solver::get_time_needed( const clock_t aTime )
    {
        moris::real tDeltaTime = (moris::real) ( clock() - aTime ) / CLOCKS_PER_SEC;

        moris::real tDeltaTimeMax   = tDeltaTime;

        MPI_Allreduce( &tDeltaTime, &tDeltaTimeMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );

        return tDeltaTimeMax;
    }

//--------------------------------------------------------------------------------------------------------------------------
    void Newton_Solver::get_full_solution( moris::Matrix< DDRMat > & LHSValues )
    {
        //mNonlinearProblem->get_full_vector()->extract_copy( LHSValues );
    }

    //--------------------------------------------------------------------------------------------------------------------------
    void Newton_Solver::get_solution( moris::Matrix< DDRMat > & LHSValues )
    {
        mNonlinearProblem->get_full_vector()->extract_copy( LHSValues );
    }
//--------------------------------------------------------------------------------------------------------------------------
    void Newton_Solver::extract_my_values( const moris::uint             & aNumIndices,
                                           const moris::Matrix< DDSMat > & aGlobalBlockRows,
                                           const moris::uint             & aBlockRowOffsets,
                                                 moris::Matrix< DDRMat > & LHSValues )
    {
        mNonlinearProblem->get_full_vector()->extract_my_values( aNumIndices, aGlobalBlockRows, aBlockRowOffsets, LHSValues );
    }
//--------------------------------------------------------------------------------------------------------------------------
    void Newton_Solver::set_nonlinear_solver_parameters()
    {
        // Allowable Newton solver iterations
        mParameterListNonlinearSolver.insert( "NLA_max_iter", 100 );

        // Allowable Newton irelative residual
        mParameterListNonlinearSolver.insert( "NLA_rel_residual" , 1e-08 );

        // Desired total residual norm drop
        mParameterListNonlinearSolver.insert( "NLA_tot_res_norm_drop" , 1e-08 );

        // Desired total residual norm
        mParameterListNonlinearSolver.insert( "NLA_tot_res_norm" , 1e-12 );

        // Maximal residual norm drop
        mParameterListNonlinearSolver.insert( "NLA_max_res_norm_drop" , 1e-6 );

        // Maximal number of linear solver restarts on fail
        mParameterListNonlinearSolver.insert( "NLA_max_lin_solver_restarts" , 0 );

        // Maximal number of linear solver restarts on fail
        mParameterListNonlinearSolver.insert( "NLA_relaxation_parameter" , 1.0 );

        // Maximal number of linear solver restarts on fail
        mParameterListNonlinearSolver.insert( "NLA_hard_break" , true );

        // Determines if lin solve should restart on fail
        mParameterListNonlinearSolver.insert( "NLA_rebuild_lin_solv_on_fail" , false );

        // Determines if newton should restart on fail
        mParameterListNonlinearSolver.insert( "NLA_rebuild_nonlin_solv_on_fail" , false );

        // Specifying the number of newton retries
        mParameterListNonlinearSolver.insert( "NLA_num_nonlin_rebuild_iterations" , 1 );

        // Determines relaxation multiplier
        mParameterListNonlinearSolver.insert( "NLA_relaxation_multiplier_on_fail" , 0.5 );

        // Determines newton maxits multiplier
        mParameterListNonlinearSolver.insert( "NLA_maxits_multiplier_on_fail" , 2 );

    }

