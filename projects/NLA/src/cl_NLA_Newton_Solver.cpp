/*
 * cl_NLA_Newton_Solver.cpp
 *
 *  Created on: Sep 21, 2018
 *      Author: schmidt
 */
#include <ctime>

#include "cl_NLA_Newton_Solver.hpp"

#include "cl_Matrix_Vector_Factory.hpp"
#include "cl_Linear_Solver.hpp"
#include "cl_Solver_Input.hpp"
#include "cl_DistLinAlg_Enums.hpp"
#include "cl_Vector.hpp"

#include "cl_Communication_Tools.hpp"

namespace moris
{
    namespace NLA
    {
///    Newton_Solver::Newton_Solver( std::shared_ptr< Linear_Solver > aLinearSolver ) : mLinearSolver( aLinearSolver )
///    {
///        Solver_Input * tInput = mLinearSolver->get_solver_input();
///
///        Matrix_Vector_Factory    tMatFactory;
///
///        // create map object
///        mMap = tMatFactory.create_map( tInput->get_num_my_dofs(),
///                                       tInput->get_my_local_global_map(),
///                                       tInput->get_constr_dof(),
///                                       tInput->get_my_local_global_overlapping_map());
///
///        // Build RHS/LHS vector
///        mVectorFreeSol = tMatFactory.create_vector( tInput, mMap, VectorType::FREE );
///        mVectorFullSol = tMatFactory.create_vector( tInput, mMap, VectorType::FREE );
///
///        mPrevVectorFreeSol = tMatFactory.create_vector( tInput, mMap, VectorType::FREE );
///        mPrevVectorFullSol = tMatFactory.create_vector( tInput, mMap, VectorType::FREE );
///        //mVectorFullSol = tMatFactory.create_vector( tInput, mMap, VectorType::FULL_OVERLAPPING );
///
///        this->set_nonlinear_solver_parameters();
///    }

//--------------------------------------------------------------------------------------------------------------------------
    Newton_Solver::~Newton_Solver()
    {
        delete( mVectorFullSol );
        delete( mVectorFreeSol );
        delete( mPrevVectorFreeSol );
        delete( mPrevVectorFullSol );
        delete( mMap );
    }

//--------------------------------------------------------------------------------------------------------------------------
    void Newton_Solver::set_linear_solver( std::shared_ptr< Linear_Solver > aLinearSolver )
    {
        // Save pointer to linear solver as member variable
        mLinearSolver = aLinearSolver;

        // Get linear solver interface
        Solver_Input * tInput = mLinearSolver->get_solver_input();

        // Build Matrix vector factory
        Matrix_Vector_Factory    tMatFactory;

        // create map object
        mMap = tMatFactory.create_map( tInput->get_num_my_dofs(),
                                       tInput->get_my_local_global_map(),
                                       tInput->get_constr_dof(),
                                       tInput->get_my_local_global_overlapping_map());

        // Build free and full vector
        mVectorFreeSol = tMatFactory.create_vector( tInput, mMap, VectorType::FREE );
        mVectorFullSol = tMatFactory.create_vector( tInput, mMap, VectorType::FULL_OVERLAPPING );

        mPrevVectorFreeSol = tMatFactory.create_vector( tInput, mMap, VectorType::FREE );
        mPrevVectorFullSol = tMatFactory.create_vector( tInput, mMap, VectorType::FULL_OVERLAPPING );
        //mVectorFullSol = tMatFactory.create_vector( tInput, mMap, VectorType::FULL_OVERLAPPING );

        mVectorFullSol->vec_put_scalar( 0.0 );

        // Set default parameters in parameter list for nonlinear solver
        this->set_nonlinear_solver_parameters();
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

        moris::sint tRebuildIterations = mParameterListNonlinearSolver.get< moris::sint >( "NLA_num_nonlin_rebuild_iterations" );

        // Newton retry loop
        for ( moris::sint Ir = 0; Ir < tRebuildIterations; ++ Ir)
        {
            // Reset the old state
            if ( Ir > 0 )
            {
                tRelaxation *= mParameterListNonlinearSolver.get< moris::real >( "NLA_relaxation_multiplier_on_fail" );
                tMaxIts      = (moris::sint) ( tMaxIts * mParameterListNonlinearSolver.get< moris::sint >( "NLA_maxits_multiplier_on_fail" ));

                mVectorFullSol->vec_plus_vec( 1.0, *mPrevVectorFreeSol, 0.0 );

                if ( par_rank() == 0 )
                {
                    fprintf( stdout,"\n ... Previous nonlinear solve failed ... Retrying with: max nonlinear iterations = %i; relaxation = %f\n", tMaxIts, tRelaxation );
                }
            }

            // Newton loop
            for ( moris::sint It = 1; It <= tMaxIts; ++It )
            {
                mPrevVectorFreeSol->vec_plus_vec( 1.0, *mVectorFreeSol, 0.0 );
                mPrevVectorFullSol->import_local_to_global( *mPrevVectorFreeSol );

                clock_t tNewtonLoopStartTime = clock();
                clock_t tStartAssemblyTime = clock();

                // Set VectorFreeSol and LHS
                mVectorFreeSol->import_local_to_global( *mVectorFullSol );
                mLinearSolver->get_solver_LHS()->vec_plus_vec( 1.0, *mVectorFreeSol, 0.0 );

                // assemble RHS and Jac
                mLinearSolver->assemble_residual_and_jacobian( mVectorFullSol );


                tMaxAssemblyTime = get_time_needed( tStartAssemblyTime );

                bool tHartBreak = false;
                tIsConverged = this->check_for_convergence( It, refNorm, tMaxAssemblyTime, tHartBreak );

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

                mVectorFreeSol->vec_plus_vec( tRelaxation, *mLinearSolver->get_solver_LHS(), 1.0 );

                // Update the SolVecNorm
                // solNorm = mVectorFreeSol.Norm2();

                mVectorFullSol->import_local_to_global( *mVectorFreeSol );

                tMaxNewTime = get_time_needed( tNewtonLoopStartTime );

                std::cout<<tMaxNewTime<<std::endl;
            }

            // Check if Newton converged. Break retry loop if its converged
            if ( tIsConverged )
            {
                break;
            }
        }
    }

//--------------------------------------------------------------------------------------------------------------------------
    void Newton_Solver::solve_linear_system( moris::sint & aIter,
                                             bool        & aHardBreak )
    {
        moris::sint tErrorStatus = 0;
        moris::sint tMaxNumLinRestarts  = mParameterListNonlinearSolver.get< moris::sint >( "NLA_max_lin_solver_restarts" );
        moris::sint tTryRestartOnFailIt = 1;

        // Solve linear system
        tErrorStatus = mLinearSolver->solve_linear_system();

        // Restart the linear solver using the current solution as an initial guess if the previous linear solve failed
        while ( tErrorStatus !=0 && tTryRestartOnFailIt <= tMaxNumLinRestarts )
        {
            if ( par_rank() == 0 )
            {
                // Compute current solution vector norm
                moris::real tSolVecNorm = mLinearSolver->get_solver_LHS()->vec_norm2();

                fprintf( stdout, " ... Previous linear solve failed. Trying restart %i of %i, using current solution with SolVecNorm = %5.15e as an initial guess. \n",
                                       tTryRestartOnFailIt, tMaxNumLinRestarts, tSolVecNorm);
            }

            // Re-solve scaled linear system with current solution as an initial guess
            tErrorStatus = mLinearSolver->solve_linear_system();

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
    bool Newton_Solver::check_for_convergence(       moris::sint & aIt,
                                                     moris::real & aRefNorm,
                                               const moris::real & aAssemblyTime,
                                                     bool        & aHartBreak )
    {
        bool tIsConverged = false;
        moris::real resNorm = 0.0;
        moris::real solNorm = mVectorFreeSol->vec_norm2();

        resNorm = mLinearSolver->get_solver_RHS()->vec_norm2();

        if ( aIt <= 1)
        {
            aRefNorm = resNorm;
            std::printf( " ... refNorm for pseudo-time step is %+1.15e\n", aRefNorm );

            if ( par_rank() == 0 )
            {
                fprintf( stdout, "        NlinIt  |  NlinResNorm            |  NlinResDrop  |  SolVecNorm             |" );

                fprintf( stdout, "|  LinAsmTime  |  NewItrTime\n" );

                // Print solution vector norm solNorm before the first Newton solve
                fprintf( stdout, "         %-5i  |  %-15.15e  |  %-11.5e  |  %-10.15e  |", 1, resNorm, 0.0, solNorm );

                fprintf( stdout, "|  %9.4e  |  %9.4e \n", aAssemblyTime, 0.0 );
            }
        }

        MORIS_ERROR( !( std::isnan( resNorm ) || std::isinf( resNorm )), "Newton_Solver::check_for_convergence(): Residual contains NAN or INF, exiting!");

        MORIS_ERROR( !( resNorm > 1e20 ), "Newton_Solver::check_for_convergence(): Residual Norm has exceeded 1e20");

        // Check for convergence
        if ( ( aIt > 1 ) && ( ( resNorm/aRefNorm ) < mParameterListNonlinearSolver.get< moris::real >( "NLA_tot_res_norm_drop" ) ) )
        {
            if ( par_rank() == 0 )
            {
                fprintf( stdout, "         %-5i  |  %-15.15e  |  %-11.5e  |  %-10.15e  |", aIt, resNorm, ( resNorm/aRefNorm ), solNorm );

                fprintf( stdout, "|  %-10.4e  |  (NlinResDrop < %6.1e)\n", aAssemblyTime, mParameterListNonlinearSolver.get< moris::real >( "NLA_tot_res_norm_drop" ) );
            }

            tIsConverged = true;
        }
        else if ( ( aIt > 1 ) && ( resNorm < mParameterListNonlinearSolver.get< moris::real >( "NLA_tot_res_norm" ) ) )
        {
            if ( par_rank() == 0 )
            {
                fprintf( stdout, "         %-5i  |  %-15.15e  |  %-11.5e  |  %-10.15e  |", aIt, resNorm, (resNorm/aRefNorm), solNorm );

                fprintf( stdout, "|  %-10.4e  |  (NlinResNorm < %6.1e)\n", aAssemblyTime, mParameterListNonlinearSolver.get< moris::real >( "NLA_tot_res_norm" ) );
            }

            tIsConverged = true;
        }
        else if ( ( aIt > 1 ) && ( ( resNorm/aRefNorm ) > mParameterListNonlinearSolver.get< moris::real >( "NLA_max_res_norm_drop" ) ) )
        {
            // case for residual drop getting too big, not converged, need to retry
            if ( par_rank() == 0 )
            {
                fprintf( stdout, "         %-5i  |  %-15.15e  |  %-11.5e  |  %-10.15e  |", aIt, resNorm, ( resNorm/aRefNorm ), solNorm );

                fprintf( stdout, "|  %-10.4e  |  (MaxResNormDrop > %6.1e)\n", aAssemblyTime, mParameterListNonlinearSolver.get< moris::real >( "NLA_max_res_norm_drop" ) );
            }

            tIsConverged = false;

            if ( mParameterListNonlinearSolver.get< bool >( "NLA_hard_break" ) )
            {
                aIt = mParameterListNonlinearSolver.get< moris::sint >( "NLA_max_iter" );
                aHartBreak = true;
            }
        }

        return tIsConverged;
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
        mVectorFullSol->extract_copy( LHSValues );
    }

    //--------------------------------------------------------------------------------------------------------------------------
    void Newton_Solver::get_solution( moris::Matrix< DDRMat > & LHSValues )
    {
        mVectorFreeSol->extract_copy( LHSValues );
    }
//--------------------------------------------------------------------------------------------------------------------------
    void Newton_Solver::extract_my_values( const moris::uint             & aNumIndices,
                                           const moris::Matrix< DDSMat > & aGlobalBlockRows,
                                           const moris::uint             & aBlockRowOffsets,
                                                 moris::Matrix< DDRMat > & LHSValues )
    {
        mVectorFullSol->extract_my_values( aNumIndices, aGlobalBlockRows, aBlockRowOffsets, LHSValues );
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

    }
}
