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

#include "cl_Matrix_Vector_Factory.hpp"
#include "cl_DLA_Linear_Solver_Algorithm.hpp"
#include "cl_DLA_Solver_Interface.hpp"
#include "cl_DLA_Linear_Solver.hpp"
#include "cl_SOL_Enums.hpp"
#include "cl_SOL_Dist_Vector.hpp"
#include "cl_DLA_Linear_Problem.hpp"

#include "cl_Communication_Tools.hpp"

using namespace moris;
using namespace NLA;
using namespace dla;

//--------------------------------------------------------------------------------------------------------------------------

    Newton_Solver::Newton_Solver()
    {
        mLinSolverManager = new dla::Linear_Solver();

        // Set default parameters in parameter list for nonlinear solver
        this->set_nonlinear_solver_parameters();
    }

//--------------------------------------------------------------------------------------------------------------------------

    Newton_Solver::Newton_Solver( const ParameterList aParameterlist ) : Nonlinear_Algorithm( aParameterlist )
    {
        mLinSolverManager = new dla::Linear_Solver();
    }


//--------------------------------------------------------------------------------------------------------------------------

    Newton_Solver::Newton_Solver( dla::Linear_Solver * aLinSolver ) : Nonlinear_Algorithm()
    {
        mLinSolverManager = aLinSolver;

    }

//--------------------------------------------------------------------------------------------------------------------------

    Newton_Solver::~Newton_Solver()
    {
        if(mLinSolverOwned)
        {
            delete mLinSolverManager;
        }
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
        moris::real tMaxNewTime      = 0.0;
        moris::real tMaxAssemblyTime = 0.0;
        //moris::real tErrorStatus     = 0;



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

                tMaxAssemblyTime = this->calculate_time_needed( tStartAssemblyTime );

                bool tHartBreak = false;
                Convergence tConvergence;
                tIsConverged = tConvergence.check_for_convergence( this,
                                                                   It,
                                                                   mMyNonLinSolverManager->get_ref_norm(),
                                                                   mMyNonLinSolverManager->get_residual_norm(),
                                                                   tMaxAssemblyTime,
                                                                   tMaxNewTime,
                                                                   tHartBreak );

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

                ( mNonlinearProblem->get_full_vector() )
                                   ->vec_plus_vec( -tRelaxation, *mNonlinearProblem->get_linearized_problem()
                                                                                   ->get_full_solver_LHS(), 1.0 );

                // Update the SolVecNorm
                // solNorm = mVectorFreeSol.Norm2();

                //mFullVector->import_local_to_global( *mVectorFreeSol );

                tMaxNewTime = this->calculate_time_needed( tNewtonLoopStartTime );

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
    void Newton_Solver::solve_linear_system( moris::sint & aIter,
                                             bool        & aHardBreak )
    {
        // Solve linear system
        mLinSolverManager->solver_linear_system( mNonlinearProblem->get_linearized_problem(), aIter );
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
    void Newton_Solver::extract_my_values( const moris::uint                            & aNumIndices,
                                           const moris::Matrix< DDSMat >                & aGlobalBlockRows,
                                           const moris::uint                            & aBlockRowOffsets,
                                                 moris::Cell< moris::Matrix< DDRMat > > & LHSValues )
    {
        mNonlinearProblem->get_full_vector()->extract_my_values( aNumIndices, aGlobalBlockRows, aBlockRowOffsets, LHSValues );
    }
//--------------------------------------------------------------------------------------------------------------------------


