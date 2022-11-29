/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Linear_Solver_Test.cpp
 *
 */

#include "catch.hpp"
#include "fn_equal_to.hpp"    // ALG/src
#include "typedefs.hpp"       // COR/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "Epetra_FECrsMatrix.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_FECrsGraph.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_FEVector.h"
#include "Epetra_IntVector.h"

#include "cl_Communication_Manager.hpp"      // COM/src/
#include "cl_Communication_Tools.hpp"        // COM/src/
#include "cl_DLA_Linear_Solver_Aztec.hpp"    // DLA/src/

#include "cl_SOL_Matrix_Vector_Factory.hpp"    // DLA/src/
#include "cl_Solver_Interface_Proxy.hpp"       // DLA/src/
#include "cl_DLA_Solver_Factory.hpp"           // DLA/src/

#include "cl_DLA_Linear_System_Trilinos.hpp"    // DLA/src/

extern moris::Comm_Manager gMorisComm;
namespace moris
{
    namespace dla
    {
        TEST_CASE( "Linear Solver Trilinos", "[Linear Solver],[DistLinAlg],[Linear Solver test]" )
        {
            if ( par_size() == 4 )
            {
                /*!
                 * Create solver interface with Solver_Interface_Proxy
                 *
                 * \code{.cpp}
                 * Solver_Interface * tSolverInterface = new Solver_Interface_Proxy( );
                 * \endcode
                 */
                Solver_Interface* tSolverInterface = new Solver_Interface_Proxy();

                /*!
                 * Create solver factory
                 *
                 * \code{.cpp}
                 * Solver_Factory tSolFactory;
                 * \endcode
                 */
                Solver_Factory tSolFactory;

                // create solver object
                Linear_Problem* tLinProblem = tSolFactory.create_linear_system( tSolverInterface, sol::MapType::Epetra );

                tLinProblem->assemble_residual_and_jacobian();

                // call solve
                tLinProblem->solve_linear_system();

                // Set solution vector
                moris::Matrix< DDRMat > tSol;
                tLinProblem->get_solution( tSol );

                print( tSol, "tSol" );

                // Check if solution corresponds to given solution
                if ( par_rank() == 0 )
                {
                    CHECK( equal_to( tSol( 0, 0 ), -0.0138889, 1.0e+08 ) );
                    CHECK( equal_to( tSol( 5, 0 ), -0.00694444, 1.0e+08 ) );
                }
                if ( par_rank() == 3 )
                {
                    CHECK( equal_to( tSol( 3, 0 ), -0.0138889, 1.0e+08 ) );
                }

                // delete tEpetraComm;
                delete ( tSolverInterface );
                delete ( tLinProblem );
            }
        }

        TEST_CASE( "Linear Solver Aztec", "[Linear Solver Aztec],[Linear Solver],[DistLinAlg]" )
        {
            if ( par_size() == 4 )
            {
                /*!
                 * Create solver interface with Solver_Interface_Proxy
                 *
                 * \code{.cpp}
                 * Solver_Interface * tSolverInterface = new Solver_Interface_Proxy( );
                 * \endcode
                 */
                Solver_Interface* tSolverInterface = new Solver_Interface_Proxy();

                /*!
                 * Create solver factory
                 *
                 * \code{.cpp}
                 * Solver_Factory tSolFactory;
                 * \endcode
                 */
                Solver_Factory tSolFactory;

                /*!
                 * Create linear problem and linear solver
                 *
                 * \code{.cpp}
                 * Linear_Problem * tLinProblem = tSolFactory.create_linear_system( tSolverInterface, sol::MapType::Epetra );
                 * std::shared_ptr< Linear_Solver_Algorithm > tLinSolver = tSolFactory.create_solver( sol::SolverType::AZTEC_IMPL );
                 * \endcode
                 */
                Linear_Problem* tLinProblem = tSolFactory.create_linear_system( tSolverInterface, sol::MapType::Epetra );

                std::shared_ptr< Linear_Solver_Algorithm > tLinSolver = tSolFactory.create_solver( sol::SolverType::AZTEC_IMPL );

                /*!
                 * Assemble linear problem.
                 *
                 * \code{.cpp}
                 * tLinProblem->assemble_residual_and_jacobian();
                 * \endcode
                 */
                tLinProblem->assemble_residual_and_jacobian();

                /*!
                 * Set linear solver parameters.
                 *
                 * \code{.cpp}
                 * tLinSolver->set_param("AZ_precond") = AZ_dom_decomp;
                 * tLinSolver->set_param("AZ_max_iter") = 200;
                 * tLinSolver->set_param("AZ_diagnostics") = AZ_none;
                 * tLinSolver->set_param("AZ_output") = AZ_none;
                 * \endcode
                 */
                tLinSolver->set_param( "AZ_precond" )     = AZ_dom_decomp;
                tLinSolver->set_param( "AZ_max_iter" )    = 200;
                tLinSolver->set_param( "AZ_diagnostics" ) = AZ_none;
                tLinSolver->set_param( "AZ_output" )      = AZ_none;

                /*!
                 * Solver linear system
                 *
                 * \code{.cpp}
                 * tLinSolver->solve_linear_system();
                 * \endcode
                 */
                tLinSolver->solve_linear_system( tLinProblem );

                /*!
                 * extract solution
                 *
                 * \code{.cpp}
                 * moris::Matrix< DDRMat > tSol;
                 * tLinSystem->get_solution( tSol );
                 * \endcode
                 */
                moris::Matrix< DDRMat > tSol;
                tLinProblem->get_solution( tSol );

                // Check if solution corresponds to given solution
                if ( par_rank() == 0 )
                {
                    CHECK( equal_to( tSol( 0, 0 ), -0.0138889, 1.0e+08 ) );
                    CHECK( equal_to( tSol( 5, 0 ), -0.00694444, 1.0e+08 ) );
                }
                if ( par_rank() == 3 )
                {
                    CHECK( equal_to( tSol( 3, 0 ), -0.0138889, 1.0e+08 ) );
                }

                // delete tEpetraComm;
                delete ( tSolverInterface );
                delete ( tLinProblem );
            }
        }

        TEST_CASE( "Linear Solver Belos multiple RHS", "[Linear Solver multiple RHS],[Linear Solver],[DistLinAlg]" )
        {
            if ( par_size() == 1 )
            {
                Solver_Interface* tSolverInterface = new Solver_Interface_Proxy( 2 );

                Solver_Factory tSolFactory;

                Linear_Problem* tLinProblem = tSolFactory.create_linear_system( tSolverInterface, sol::MapType::Epetra );

                std::shared_ptr< Linear_Solver_Algorithm > tLinSolver = tSolFactory.create_solver( sol::SolverType::BELOS_IMPL );

                //        tLinProblem->assemble_residual_and_jacobian();
                tLinProblem->assemble_jacobian();
                tLinProblem->assemble_residual();

                tLinSolver->solve_linear_system( tLinProblem );

                // get solution vector (here: solution vector has only unconstrained dofs)
                moris::Matrix< DDRMat > tSol;
                tLinProblem->get_solution( tSol );

                // Check if solution corresponds to given solution
                CHECK( equal_to( tSol( 5, 0 ), -0.0138889, 1.0e+08 ) );
                CHECK( equal_to( tSol( 12, 0 ), -0.00694444, 1.0e+08 ) );

                CHECK( equal_to( tSol( 5, 1 ), -0.0138889, 1.0e+08 ) );
                CHECK( equal_to( tSol( 12, 1 ), -0.00694444, 1.0e+08 ) );

                // delete local variables
                delete ( tSolverInterface );
                delete ( tLinProblem );
            }
        }

#ifdef MORIS_HAVE_PETSC
        TEST_CASE( "Linear System PETSc single RHS", "[Linear Solver single RHS],[Linear Solver],[DistLinAlg]" )
        {
            if ( par_size() == 1 )
            {
                Solver_Interface* tSolverInterface = new Solver_Interface_Proxy( 1 );

                Solver_Factory tSolFactory;

                // create solver object
                Linear_Problem* tLinProblem = tSolFactory.create_linear_system( tSolverInterface, sol::MapType::Petsc, true );

                tLinProblem->assemble_residual_and_jacobian();

                // call solve
                tLinProblem->solve_linear_system();

                // get solution vector (here: solution vector includes the 3 constrained dofs)
                moris::Matrix< DDRMat > tSol;
                tLinProblem->get_solution( tSol );

                CHECK( equal_to( tSol( 5 + 3, 0 ), -0.0138889, 1.0e+08 ) );
                CHECK( equal_to( tSol( 12 + 3, 0 ), -0.00694444, 1.0e+08 ) );

                // delete local variables
                delete ( tSolverInterface );
                delete ( tLinProblem );
            }
        }

        TEST_CASE( "Linear Solver Petsc", "[Linear Solver Petsc],[Linear Solver],[DistLinAlg]" )
        {
            if ( par_size() == 4 )
            {
                Solver_Interface* tSolverInterface = new Solver_Interface_Proxy();

                Solver_Factory tSolFactory;

                Linear_Problem* tLinProblem = tSolFactory.create_linear_system( tSolverInterface, sol::MapType::Petsc, true );

                std::shared_ptr< Linear_Solver_Algorithm > tLinSolver = tSolFactory.create_solver( sol::SolverType::PETSC );

                tLinProblem->assemble_residual_and_jacobian();

                tLinSolver->set_param( "KSPType" ) = std::string( "fgmres" );
                tLinSolver->set_param( "PCType" )  = std::string( "none" );
                tLinSolver->set_param( "ILUFill" ) = 3;

                tLinSolver->solve_linear_system( tLinProblem );

                moris::Matrix< DDRMat > tSol;
                tLinProblem->get_solution( tSol );

                // Check if solution corresponds to given solution
                if ( par_rank() == 0 )
                {
                    print( tSol, "tSol" );
                    CHECK( equal_to( tSol( 2, 0 ), -0.0138889, 1.0e+08 ) );
                    CHECK( equal_to( tSol( 7, 0 ), -0.00694444, 1.0e+08 ) );
                }
                if ( par_rank() == 3 )
                {
                    print( tSol, "tSol" );
                    CHECK( equal_to( tSol( 3, 0 ), -0.0138889, 1.0e+08 ) );
                }

                delete ( tLinProblem );
                delete ( tSolverInterface );
            }
        }
#endif

    }    // namespace dla
}    // namespace moris
