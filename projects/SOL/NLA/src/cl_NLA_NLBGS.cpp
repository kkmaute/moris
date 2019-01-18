/*
 * cl_NLA_NLBGS.cpp
 *
 *  Created on: Jan 18, 2018
 *      Author: schmidt
 */
#include <ctime>

#include "cl_NLA_NLBGS.hpp"

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
    void NonLinBlockGaussSeidel::solver_nonlinear_system( )
    {
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

                ( mNonlinearProblem->get_full_vector())->vec_plus_vec( -tRelaxation, *mNonlinearProblem->get_linearized_problem()->get_full_solver_LHS(), 1.0 );

                tMaxNewTime = get_time_needed( tNewtonLoopStartTime );
            }

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

