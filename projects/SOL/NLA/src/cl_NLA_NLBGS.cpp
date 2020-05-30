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
#include "cl_Matrix_Vector_Factory.hpp"
#include "cl_DLA_Linear_Solver_Algorithm.hpp"
#include "cl_DLA_Solver_Interface.hpp"
#include "cl_SOL_Enums.hpp"
#include "cl_SOL_Dist_Vector.hpp"

#include "cl_Communication_Tools.hpp"
#include "cl_SOL_Warehouse.hpp"

using namespace moris;
using namespace NLA;
using namespace dla;

    NonLinBlockGaussSeidel::NonLinBlockGaussSeidel()
    {
        // Set default parameters in parameter list for nonlinear solver
        this->set_nonlinear_solver_parameters();
    }

    NonLinBlockGaussSeidel::NonLinBlockGaussSeidel( const ParameterList aParameterlist ) : Nonlinear_Algorithm( aParameterlist )
    {

    }

//    NonLinBlockGaussSeidel::NonLinBlockGaussSeidel( Solver_Interface * aSolverInterface ) : Nonlinear_Solver( aSolverInterface )
//    {
//        // Set default parameters in parameter list for nonlinear solver
//        this->set_nonlinear_solver_parameters();
//    }

    NonLinBlockGaussSeidel::~NonLinBlockGaussSeidel()
    {
    }

//--------------------------------------------------------------------------------------------------------------------------
    void NonLinBlockGaussSeidel::solver_nonlinear_system( Nonlinear_Problem * aNonlinearProblem )
    {
        moris::sint tMaxIts = 10;
        moris::uint tNonLinSysStartIt = 0;
        moris::uint tNumNonLinSystems = mMyNonLinSolverManager->get_dof_type_list().size();
        bool tIsConverged             = false;

        // NLBGS loop
        for ( sint It = 1; It <= tMaxIts; ++It )
        {
            moris::real tMaxNewTime      = 0.0;

            //get_nonlinear_problem()
            clock_t tNewtonLoopStartTime = clock();

            // Loop over all non-linear systems
            for (uint Ik = tNonLinSysStartIt ; Ik < tNumNonLinSystems;Ik++)
            {
//                  moris::sint tNonlinSolverManagerIndex = mMyNonLinSolverManager->get_solver_warehouse()
//                                                                                ->get_nonlinear_solver_manager_index( mMyNonLinSolverManager->get_sonlinear_solver_manager_index(), Ik );
//
//                  Nonlinear_Solver * tMySubSolverManager = mMyNonLinSolverManager->get_solver_warehouse()
//                                                                                         ->get_nonliner_solver_manager_list()( tNonlinSolverManagerIndex );

                mMyNonLinSolverManager->get_sub_nonlinear_solver( Ik )->solve( aNonlinearProblem->get_full_vector() );
            } // end loop over all non-linear sub-systems

            tMaxNewTime = this->calculate_time_needed( tNewtonLoopStartTime );
            bool tHartBreak = false;

            this->compute_norms( It );

            Convergence tConvergence;

            tIsConverged = tConvergence.check_for_convergence( this,
                                                               It,
                                                               mMyNonLinSolverManager->get_ref_norm(),
                                                               mMyNonLinSolverManager->get_residual_norm(),
                                                               tMaxNewTime,
                                                               tHartBreak);

            if ( tIsConverged )
            {
                if ( tHartBreak )
                {
                    continue;
                }
                break;
            }
        } // end loop for NLBGS iterations
    }

//--------------------------------------------------------------------------------------------------------------------------
    void NonLinBlockGaussSeidel::solve_linear_system( moris::sint & aIter,
                                                      bool        & aHardBreak )
    {
        // Solve linear system
    }

    void NonLinBlockGaussSeidel::compute_norms( const moris::sint aIter )
    {
        moris::uint tNumNonLinSystems = mMyNonLinSolverManager->get_dof_type_list().size();
        moris::uint tNonLinSysStartIt = 0;

        if( aIter == 1 )
        {
            moris::real tSqrtRefNorm = 0.0;

            // Loop over all non-linear systems
            for (uint Ik = tNonLinSysStartIt ; Ik < tNumNonLinSystems; Ik++)
            {
//                moris::sint tNonlinSolverManagerIndex = mMyNonLinSolverManager->get_solver_warehouse()
//                                                                              ->get_nonlinear_solver_manager_index( mMyNonLinSolverManager->get_sonlinear_solver_manager_index(), Ik );
//
//                moris::real tSubSolverRefNorm = mMyNonLinSolverManager->get_solver_warehouse()
//                                                                      ->get_nonliner_solver_manager_list()( tNonlinSolverManagerIndex )
//                                                                      ->get_ref_norm();

                 moris::real tSubSolverRefNorm = mMyNonLinSolverManager->get_sub_nonlinear_solver( Ik )->get_ref_norm();

                tSqrtRefNorm = tSqrtRefNorm + std::pow( tSubSolverRefNorm, 2 );
            }

            mMyNonLinSolverManager->get_ref_norm() = std::sqrt( tSqrtRefNorm );
        }

        moris::real tSqrtResNorm = 0.0;

        // Loop over all non-linear systems
        for (uint Ik = tNonLinSysStartIt ; Ik < tNumNonLinSystems; Ik++)
        {
//            moris::sint tNonlinSolverManagerIndex = mMyNonLinSolverManager->get_solver_warehouse()
//                                                                          ->get_nonlinear_solver_manager_index( mMyNonLinSolverManager->get_sonlinear_solver_manager_index(), Ik );
//
//            moris::real tSubSolverFirstResNorm = mMyNonLinSolverManager->get_solver_warehouse()
//                                                                  ->get_nonliner_solver_manager_list()( tNonlinSolverManagerIndex )
//                                                                  ->get_ref_norm();

            moris::real tSubSolverResNorm = mMyNonLinSolverManager->get_sub_nonlinear_solver( Ik )->get_residual_norm();

            tSqrtResNorm = tSqrtResNorm + std::pow( tSubSolverResNorm, 2 );
        }

        mMyNonLinSolverManager->get_residual_norm() = std::sqrt( tSqrtResNorm );
    }


