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
        bool tRebuildJacobian        = true;
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
                if ( It > 1 )
                {
                    tRebuildJacobian = mParameterListNonlinearSolver.get< bool >( "NLA_rebuild_jacobian" );
                }

                if ( It == 1 && mParameterListNonlinearSolver.get< sint >( "NLA_restart" ) != 0 )
                {
                    sint tRestart = mParameterListNonlinearSolver.get< sint >( "NLA_restart" );
                    mNonlinearProblem->build_linearized_problem( tRebuildJacobian, It, tRestart );
                }
                else
                {
                    mNonlinearProblem->build_linearized_problem( tRebuildJacobian, It );
                }

                tMaxAssemblyTime = get_time_needed( tStartAssemblyTime );

                bool tHartBreak = false;
                Convergence tConvergence;
                tIsConverged = tConvergence.check_for_convergence( this, It, refNorm, tMaxAssemblyTime, tMaxNewTime, tHartBreak );

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

                (mNonlinearProblem->get_full_vector())->vec_plus_vec( -tRelaxation, *mNonlinearProblem->get_linearized_problem()->get_full_solver_LHS(), 1.0 );

                // Update the SolVecNorm
                // solNorm = mVectorFreeSol.Norm2();

                //mVectorFullSol->import_local_to_global( *mVectorFreeSol );

                tMaxNewTime = get_time_needed( tNewtonLoopStartTime );

                //std::cout<<"Total iter time "<<tMaxNewTime<<std::endl;
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
//        moris::sint tMaxIts  = mParameterListNonlinearSolver.get< moris::sint >( "NLA_max_iter" );
//        //moris::real tRelRes = mParameterListNonlinearSolver.get< moris::real >( "NLA_rel_residual" );
//        moris::real tRelaxation = mParameterListNonlinearSolver.get< moris::real >( "NLA_relaxation_parameter" );
////        moris::sint tRebuildIterations = mParameterListNonlinearSolver.get< moris::sint >( "NLA_num_nonlin_rebuild_iterations" );
//
//        bool tIsConverged            = false;
//        bool tRebuildJacobian        = true;
//        moris::real refNorm          = 0.0;
//        moris::real tMaxNewTime      = 0.0;
//        moris::real tMaxAssemblyTime = 0.0;
//        //moris::real tErrorStatus     = 0;
//
//        // Newton retry loop
////        for ( moris::sint Ir = 0; Ir < tRebuildIterations; ++ Ir)
////        {
////            // Reset the old state
////            if ( Ir > 0 )
////            {
////                tRelaxation *= mParameterListNonlinearSolver.get< moris::real >( "NLA_relaxation_multiplier_on_fail" );
////                tMaxIts      = (moris::sint) ( tMaxIts * mParameterListNonlinearSolver.get< moris::sint >( "NLA_maxits_multiplier_on_fail" ));
////
////                mVectorFullSol->vec_plus_vec( 1.0, *mPrevVectorFreeSol, 0.0 );
////
////                if ( par_rank() == 0 )
////                {
////                    fprintf( stdout,"\n ... Previous nonlinear solve failed ... Retrying with: max nonlinear iterations = %i; relaxation = %f\n", tMaxIts, tRelaxation );
////                }
////            }
//
//            // Newton loop
//            for ( moris::sint It = 1; It <= tMaxIts; ++It )
//            {
//                std::cout<<"----1-----1-----1-----1-----1-----1----"<<std::endl;
//                  //get_nonlinear_problem()
//                clock_t tNewtonLoopStartTime = clock();
//                clock_t tStartAssemblyTime = clock();
//
//                // assemble RHS and Jac
//                if ( It > 1 )
//                {
//                    tRebuildJacobian = mParameterListNonlinearSolver.get< bool >( "NLA_rebuild_jacobian" );
//                }
//                mNonlinearProblem->build_linearized_problem( tRebuildJacobian );
//
//                tMaxAssemblyTime = get_time_needed( tStartAssemblyTime );
//
//                bool tHartBreak = false;
//                Convergence tConvergence;
//                tIsConverged = tConvergence.check_for_convergence( this, It, refNorm, tMaxAssemblyTime, tHartBreak );
//
//                if ( tIsConverged )
//                {
//                    if ( tHartBreak )
//                    {
//                        continue;
//                    }
//                    break;
//                }
//
//                // Solve linear system
//                this->solve_linear_system( It, tHartBreak );
//
////                if ( tHartBreak )
////                 {
////                     continue;
////                 }
////                 break;
//
//                //PreconTime
//                //SolveTime
//
//                (mNonlinearProblem->get_full_vector())->vec_plus_vec( tRelaxation, *mNonlinearProblem->get_linearized_problem()->get_full_solver_LHS(), 1.0 );
//
//                // Update the SolVecNorm
//                // solNorm = mVectorFreeSol.Norm2();
//
//                //mVectorFullSol->import_local_to_global( *mVectorFreeSol );
//
//                tMaxNewTime = get_time_needed( tNewtonLoopStartTime );
//
//                std::cout<<"Total iteration time "<<tMaxNewTime<<std::endl;
//            }
////
////            // Check if Newton converged. Break retry loop if its converged
////            if ( tIsConverged )
////            {
////                break;
////            }
////        }
    }

//--------------------------------------------------------------------------------------------------------------------------
    void Newton_Solver::solve_linear_system( moris::sint & aIter,
                                             bool        & aHardBreak )
    {
        // Solve linear system
        mLinSolverManager->solver_linear_system( mNonlinearProblem->get_linearized_problem(), aIter );
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
        mNonlinearProblem->get_full_vector()->extract_copy( LHSValues );
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

        // Allowable Newton solver iterations
        mParameterListNonlinearSolver.insert( "NLA_restart", 0 );

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

