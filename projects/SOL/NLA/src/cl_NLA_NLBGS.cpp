/*
 * cl_NLA_NLBGS.cpp
 *
 *  Created on: Jan 18, 2018
 *      Author: schmidt
 */
#include <ctime>

#include "cl_NLA_NLBGS.hpp"

#include "cl_NLA_Convergence.hpp"
#include "cl_NLA_Nonlinear_Solver_Manager.hpp"
#include "cl_Matrix_Vector_Factory.hpp"
#include "cl_DLA_Linear_Solver.hpp"
#include "cl_DLA_Solver_Interface.hpp"
#include "cl_DLA_Enums.hpp"
#include "cl_Vector.hpp"

#include "cl_Communication_Tools.hpp"
#include "cl_NLA_Nonlinear_Database.hpp"

using namespace moris;
using namespace NLA;
using namespace dla;

    NonLinBlockGaussSeidel::NonLinBlockGaussSeidel()
    {
        // Set default parameters in parameter list for nonlinear solver
        this->set_nonlinear_solver_parameters();
    }

    NonLinBlockGaussSeidel::NonLinBlockGaussSeidel( Solver_Interface * aSolverInterface ) : Nonlinear_Solver( aSolverInterface )
    {
        // Set default parameters in parameter list for nonlinear solver
        this->set_nonlinear_solver_parameters();
    }

    NonLinBlockGaussSeidel::~NonLinBlockGaussSeidel()
    {
    }

//--------------------------------------------------------------------------------------------------------------------------
    void NonLinBlockGaussSeidel::solver_nonlinear_system( Nonlinear_Problem * aNonlinearProblem )
    {
        moris::uint tMaxIts = 2;
        moris::uint tNonLinSysStartIt = 0;
        moris::uint tNumNonLinSystems = 2;

        // NLBGS loop
        for ( uint it = 1; it <= tMaxIts; ++it )
        {
            // Loop over all non-linear systems
            for (uint Ik = tNonLinSysStartIt ; Ik < tNumNonLinSystems;Ik++)
            {
                  moris::sint tNonlinSolverManagerIndex = mMyNonLinSolverManager->get_nonlinear_database()
                                                                                ->get_nonlinear_solver_manager_index( mMyNonLinSolverManager->get_sonlinear_solver_manager_index(), Ik );


                  Nonlinear_Solver_Manager * tMySubSolverManager = mMyNonLinSolverManager->get_nonlinear_database()->get_nonliner_solver_manager_list()(tNonlinSolverManagerIndex);

                  tMySubSolverManager->solve();

                  // moris::Cell< moris::Cell< enum MSI::Dof_Type > > tDofTypeList = mMyNonLinSolverManager->get_dof_type_list();

                  // FIXME call nonliner solver manager
//                // Solve the non-linear sub-problem, if an error occurred or the nonlinear problem did not converge, a NEGATIVE status is returned
//                NonLinStatus = mNonlinearSolverList[ik]->SolveNonlinearSystem( aMasterSolVec, aFullSolVec, aModel, ik );
//
//                // Update the total delta solution vector with the solution vector of the current nonlinear problem
//                mTimeSolver->GetMasterDeltaSolVec()->Update(1.0,*mNonlinearSolverList[ik]->GetLinearSolver()->GetSol(),0.0);
//
//                // Update the global = full delta solution vector by importing the local master delta solution vector into the global delta solution vector
//                mTimeSolver->GetDeltaSolVec()->Import(*mTimeSolver->GetMasterDeltaSolVec());
                //aNonlinearProblem->get_full_vector()
//
            } // end loop over all non-linear sub-systems
        } // end loop for NLBGS iterations
    }

//--------------------------------------------------------------------------------------------------------------------------
    void NonLinBlockGaussSeidel::solve_linear_system( moris::sint & aIter,
                                             bool        & aHardBreak )
    {
        // Solve linear system
    }

//--------------------------------------------------------------------------------------------------------------------------
    moris::real NonLinBlockGaussSeidel::get_time_needed( const clock_t aTime )
    {
        moris::real tDeltaTime = (moris::real) ( clock() - aTime ) / CLOCKS_PER_SEC;

        moris::real tDeltaTimeMax   = tDeltaTime;

        MPI_Allreduce( &tDeltaTime, &tDeltaTimeMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );

        return tDeltaTimeMax;
    }

//--------------------------------------------------------------------------------------------------------------------------
    void NonLinBlockGaussSeidel::set_nonlinear_solver_parameters()
    {

    }

