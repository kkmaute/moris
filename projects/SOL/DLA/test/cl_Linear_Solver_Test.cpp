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
#include "fn_equal_to.hpp"       // ALG/src
#include "moris_typedefs.hpp"    // COR/src
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
#include "cl_DLA_Solver_Factory.hpp"           // DLA/src/
#include "cl_SOL_Warehouse.hpp"

#include "cl_DLA_Linear_System_Trilinos.hpp"    // DLA/src/

#ifdef MORIS_HAVE_PETSC
#include "cl_DLA_Preconditioner_PETSc.hpp"
#ifdef MORIS_HAVE_SLEPC
#include "cl_DLA_Eigen_Solver_SLEPc.hpp"
#endif
#endif

#define private public
#include "cl_Solver_Interface_Proxy.hpp"    // DLA/src/
#undef private

#include "fn_PRM_SOL_Parameters.hpp"
extern moris::Comm_Manager gMorisComm;
namespace moris::dla
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


            tLinProblem->get_free_solver_LHS()->print();
            tLinProblem->get_matrix()->save_matrix_to_matlab_file( "matrixnew" );
            tLinProblem->get_matrix()->print();
            tLinProblem->get_solver_RHS()->print();

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

            Parameter_List                             tLinearSolverParameterList = prm::create_linear_algorithm_parameter_list_aztec();
            std::shared_ptr< Linear_Solver_Algorithm > tLinSolver                 = tSolFactory.create_solver( tLinearSolverParameterList );

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
             * tLinearSolverParameterList.set( "AZ_precond", AZ_dom_decomp );
             * tLinearSolverParameterList.set( "AZ_max_iter", 200 );
             * tLinearSolverParameterList.set( "AZ_diagnostics", AZ_none );
             * tLinearSolverParameterList.set( "AZ_output", AZ_none );
             * \endcode
             */
            tLinearSolverParameterList.set( "AZ_precond", AZ_dom_decomp );
            tLinearSolverParameterList.set( "AZ_max_iter", 200 );
            tLinearSolverParameterList.set( "AZ_diagnostics", AZ_none );
            tLinearSolverParameterList.set( "AZ_output", AZ_none );

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

            Parameter_List                             tLinearSolverParameterList = prm::create_linear_algorithm_parameter_list_belos();
            std::shared_ptr< Linear_Solver_Algorithm > tLinSolver                 = tSolFactory.create_solver( tLinearSolverParameterList );

            Parameter_List tParamList;
            tParamList.insert( "ifpack_prec_type", std::string( "ILU" ) );
            tParamList.insert( "ml_prec_type", "" );
            tParamList.insert( "fact: level-of-fill", 1 );
            tParamList.insert( "fact: absolute threshold", 0.0 );
            tParamList.insert( "fact: relative threshold", 1.0 );
            tParamList.insert( "fact: relax value", 0.0 );

            tParamList.insert( "schwarz: combine mode", std::string( "add" ) );
            tParamList.insert( "schwarz: compute condest", true );
            tParamList.insert( "schwarz: filter singletons", false );
            tParamList.insert( "schwarz: reordering type", "rcm" );
            tParamList.insert( "overlap-level", 0 );
            tParamList.insert( "prec_reuse", false );


            // create preconditioner
            Preconditioner_Trilinos tPreconditioner( tParamList );
            tLinSolver->set_preconditioner( &tPreconditioner );

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
            PRINT( tSol );

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

            Parameter_List tLinearSolverParameterList = prm::create_linear_algorithm_parameter_list_petsc();
            tLinearSolverParameterList.set( "KSPType", std::string( "fgmres" ) );
            tLinearSolverParameterList.set( "PCType", std::string( "none" ) );
            tLinearSolverParameterList.set( "ILUFill", 3 );
            tLinearSolverParameterList.set( "ouput_eigenspectrum", (uint)1 );
            std::shared_ptr< Linear_Solver_Algorithm > tLinSolver = tSolFactory.create_solver( tLinearSolverParameterList );

            tLinProblem->assemble_residual_and_jacobian();

            // create preconditioner
            Parameter_List tParamList;
            tParamList.insert( "PCType", std::string( "none" ) );
            Preconditioner_PETSc tPreconditioner( tParamList );
            tLinSolver->set_preconditioner( &tPreconditioner );

            tLinSolver->solve_linear_system( tLinProblem );

            moris::Matrix< DDRMat > tSol;
            tLinProblem->get_solution( tSol );

            tLinProblem->get_free_solver_LHS()->print();
            tLinProblem->get_matrix()->save_matrix_to_matlab_file( "matrixnew" );
            tLinProblem->get_matrix()->print();
            tLinProblem->get_solver_RHS()->print();
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

    TEST_CASE( "Eigen Solver Block Davidson", "[Eigen Solver Block Davidson],[Eigen Solver], [EigSolve]" )
    {
        if ( par_size() == 1 )
        {
            std::string tProblem = "Block_Davidson_Eigen";
            /*!
             * Create solver interface with Solver_Interface_Proxy
             *
             * \code{.cpp}
             * Solver_Interface * tSolverInterface = new Solver_Interface_Proxy( );
             * \endcode
             */
            Solver_Interface_Proxy* tSolverInterface = new Solver_Interface_Proxy( tProblem );

            /*!
             * Create solver factory
             *
             * \code{.cpp}
             * Solver_Factory tSolFactory;
             * \endcode
             */
            Solver_Factory tSolFactory;

            sol::SOL_Warehouse tSolverWarehouse( tSolverInterface );

            // set RHS matrix type
            std::string tRHSMatType = std::string( "MassMat" );
            // ( &tSolverWarehouse )->set_RHS_mat_type( tRHSMatType );   // fix this this function is deleted

            const enum sol::MapType tMapType = sol::MapType::Epetra;

            // Build Matrix vector factory
            sol::Matrix_Vector_Factory tMatFactory( tMapType );

            // create map object FIXME ask linear problem for map
            sol::Dist_Map* tMap = tMatFactory.create_map( tSolverInterface->get_my_local_global_map() );

            // create map object FIXME ask linear problem for map
            sol::Dist_Map* tMapFull = tMatFactory.create_full_map(
                    tSolverInterface->get_my_local_global_map(),
                    tSolverInterface->get_my_local_global_overlapping_map() );

            // create solver object
            Linear_Problem* tEigProblem = tSolFactory.create_linear_system( tSolverInterface, &tSolverWarehouse, tMap, tMapFull, tMapType );

            // set eigen algorithm parameters
            Parameter_List tLinearSolverParameterList = prm::create_eigen_algorithm_parameter_list();
            tLinearSolverParameterList.set( "Eigen_Algorithm", std::string( "EIGALG_BLOCK_DAVIDSON" ) );
            tLinearSolverParameterList.set( "Which", std::string( "LM" ) );
            tLinearSolverParameterList.set( "Verbosity", false );
            tLinearSolverParameterList.set( "Block_Size", 1 );
            tLinearSolverParameterList.set( "Num_Blocks", 3 );
            tLinearSolverParameterList.set( "NumFreeDofs", 8 );
            tLinearSolverParameterList.set( "Num_Eig_Vals", 1 );
            tLinearSolverParameterList.set( "MaxSubSpaceDims", 6 );
            tLinearSolverParameterList.set( "MaxRestarts", 20 );
            tLinearSolverParameterList.set( "Initial_Guess", 0 );
            tLinearSolverParameterList.set( "Convergence_Tolerance", 1e-05 );
            tLinearSolverParameterList.set( "Relative_Convergence_Tolerance", true );
            tLinearSolverParameterList.set( "Update_Flag", false );    // false flag is set only for unit test. Default: True

            // create eigen solver
            std::shared_ptr< Linear_Solver_Algorithm > tEigSolver = tSolFactory.create_solver( tLinearSolverParameterList );

            Parameter_List tParamList;
            tParamList.insert( "ifpack_prec_type", std::string( "Amesos" ) );
            tParamList.insert( "ml_prec_type", "" );
            tParamList.insert( "amesos: solver type", std::string( "Amesos_Pardiso" ) );
            tParamList.insert( "overlap-level", 0 );
            tParamList.insert( "prec_reuse", false );

            tParamList.insert( "schwarz: combine mode", std::string( "add" ) );
            tParamList.insert( "schwarz: compute condest", true );
            tParamList.insert( "schwarz: filter singletons", false );
            tParamList.insert( "schwarz: reordering type", "rcm" );
            tParamList.insert( "Amesos_Direct_Solver_Type", "Amesos_Klu" );

            // create preconditioner
            Preconditioner_Trilinos tPreconditioner( tParamList );
            tEigSolver->set_preconditioner( &tPreconditioner );

            tEigProblem->set_rhs_matrix_type( tRHSMatType );

            // assemble jacobian
            tEigProblem->assemble_jacobian();

            moris::sint aIter = 0;

            // call solve
            tEigSolver->solve_linear_system( tEigProblem, aIter );

            // set eigenvalue solution
            std::vector< Anasazi::Value< double > > tSol;
            tSol        = dynamic_cast< Eigen_Solver* >( tEigSolver.get() )->get_eigen_values();
            double tVal = tSol[ 0 ].realpart;

            // Check if solution corresponds to given solution
            if ( par_rank() == 0 )
            {
                CHECK( equal_to( tVal, 0.5, 1.0e+08 ) );
            }

            // delete tEpetraComm;
            delete ( tSolverInterface );
            delete ( tEigProblem );
        }
    }

//----------------------------------------------------------------------------------------
   
}    // namespace moris::dla
