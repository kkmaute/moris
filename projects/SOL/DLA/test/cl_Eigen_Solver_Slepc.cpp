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
#include "cl_DLA_Eigen_Solver_SLEPc.hpp"
#endif

#define private public
#include "cl_Solver_Interface_Proxy.hpp"    // DLA/src/
#undef private

#include "fn_PRM_SOL_Parameters.hpp"
extern moris::Comm_Manager gMorisComm;
namespace moris::dla
{


    TEST_CASE( "Eigen Solver SLEPCv1", "[Eigen Solver SLEPC],[Eigen Solver], [EigSolveSLEPC]" )
    {
        if ( par_size() == 1 )
        {
            Solver_Interface_Proxy* tSolverInterface = new Solver_Interface_Proxy( 1 );
            tSolverInterface->mMyConstraintDofs.set_size( 0, 0 );

            Solver_Factory tSolFactory;

            Linear_Problem* tEigProblem = tSolFactory.create_linear_system( tSolverInterface, sol::MapType::Petsc, true );

            Parameter_List tLinearSolverParameterListLM = prm::create_slepc_algorithm_parameter_list();
            tLinearSolverParameterListLM.set( "Eigen_Algorithm", std::string( "power" ) );
            tLinearSolverParameterListLM.set( "Which", std::string( "LM" ) );
            tLinearSolverParameterListLM.set( "Num_Eig_Vals", 1 );    // 2; 1
            tLinearSolverParameterListLM.set( "Update_Flag", false );
            tLinearSolverParameterListLM.set( "Verbosity", true );


            Parameter_List tLinearSolverParameterListSM = prm::create_slepc_algorithm_parameter_list();
            tLinearSolverParameterListSM.set( "Eigen_Algorithm", std::string( "power" ) );
            tLinearSolverParameterListSM.set( "Which", std::string( "LM" ) );
            tLinearSolverParameterListSM.set( "Num_Eig_Vals", 1 );           // 2; 1
            tLinearSolverParameterListSM.set( "STType", "shift_invert" );    // 10 shift_invert
            tLinearSolverParameterListSM.set( "Update_Flag", false );
            tLinearSolverParameterListSM.set( "Verbosity", true );

            std::shared_ptr< Linear_Solver_Algorithm > tEigSolverLM = tSolFactory.create_solver( tLinearSolverParameterListLM );
            std::shared_ptr< Linear_Solver_Algorithm > tEigSolverSM = tSolFactory.create_solver( tLinearSolverParameterListSM );

            Parameter_List tParamListKSP = prm::create_linear_algorithm_parameter_list_petsc();
            tParamListKSP.set( "KSPType", std::string( "preonly" ) );

            Parameter_List tParamListPCNone = prm::create_preconditioner_parameter_list( moris::sol::PreconditionerType::PETSC );
            tParamListPCNone.set( "PCType", std::string( "none" ) );
            Parameter_List tParamListPCMumps = prm::create_preconditioner_parameter_list( moris::sol::PreconditionerType::PETSC );
            tParamListPCMumps.set( "PCType", std::string( "mumps" ) );

            tEigSolverLM->set_sublinear_solver_options( &tParamListKSP, &tParamListPCNone );
            tEigSolverSM->set_sublinear_solver_options( &tParamListKSP, &tParamListPCMumps );

            std::string tRHSMatType = std::string( "IdentityMat" );
            tEigProblem->set_rhs_matrix_type( tRHSMatType );

            // assemble jacobian
            tEigProblem->assemble_jacobian();

            moris::sint aIter = 0;

            // call solve
            tEigSolverLM->solve_linear_system( tEigProblem, aIter );
            tEigSolverSM->solve_linear_system( tEigProblem, aIter );

            tEigProblem->get_matrix()->save_matrix_to_matlab_file( "JacobianCheck.txt" );

            Vector<real> const & tLargeEigenValues = std::dynamic_pointer_cast<Eigen_Solver_SLEPc>(tEigSolverLM)->get_eigenvalues();
            Vector<real> const & tSmallEigenValues = std::dynamic_pointer_cast<Eigen_Solver_SLEPc>(tEigSolverSM)->get_eigenvalues();

            CHECK( equal_to( tLargeEigenValues(0), 60.864058, 1.0e+08 ) );
            CHECK( equal_to( tSmallEigenValues(0), 0.0, 1.0e+08 ) );

            delete tEigProblem;
            delete tSolverInterface;
            
        }
    }

    //----------------------------------------------------------------------------------------

    TEST_CASE( "Eigen Solver SLEPC2v1", "[Eigen Solver SLEPC],[Eigen Solver], [EigSolveSLEPC]" )
    {
        if ( par_size() == 1 )
        {
            Solver_Interface_Proxy* tSolverInterface = new Solver_Interface_Proxy( 1 );
            tSolverInterface->mMyConstraintDofs.set_size( 0, 0 );

            Solver_Factory tSolFactory;

            Linear_Problem* tEigProblem = tSolFactory.create_linear_system( tSolverInterface, sol::MapType::Petsc, true );

            Parameter_List tLinearSolverParameterListLM = prm::create_slepc_algorithm_parameter_list();
            tLinearSolverParameterListLM.set( "Eigen_Algorithm", std::string( "gd" ) );
            tLinearSolverParameterListLM.set( "Which", std::string( "LM" ) );
            tLinearSolverParameterListLM.set( "Num_Eig_Vals", 1 );      // 2; 1
            tLinearSolverParameterListLM.set( "STType", "precond" );    // 10 shift_invert
            tLinearSolverParameterListLM.set( "Update_Flag", false );


            Parameter_List tLinearSolverParameterListSM = prm::create_slepc_algorithm_parameter_list();
            tLinearSolverParameterListSM.set( "Eigen_Algorithm", std::string( "gd" ) );
            tLinearSolverParameterListSM.set( "Which", std::string( "SM" ) );
            tLinearSolverParameterListSM.set( "Num_Eig_Vals", 1 );      // 2; 1
            tLinearSolverParameterListSM.set( "STType", "precond" );    // 10 shift_invert
            tLinearSolverParameterListSM.set( "Update_Flag", false );

            std::shared_ptr< Linear_Solver_Algorithm > tEigSolverLM = tSolFactory.create_solver( tLinearSolverParameterListLM );
            std::shared_ptr< Linear_Solver_Algorithm > tEigSolverSM = tSolFactory.create_solver( tLinearSolverParameterListSM );

            Parameter_List tParamListKSP = prm::create_linear_algorithm_parameter_list_petsc();
            tParamListKSP.set( "KSPType", std::string( "preonly" ) );

            Parameter_List tParamListPCNone = prm::create_preconditioner_parameter_list( moris::sol::PreconditionerType::PETSC );
            tParamListPCNone.set( "PCType", std::string( "none" ) );
            Parameter_List tParamListPCMumps = prm::create_preconditioner_parameter_list( moris::sol::PreconditionerType::PETSC );
            tParamListPCMumps.set( "PCType", std::string( "superlu-dist" ) );

            tEigSolverLM->set_sublinear_solver_options( &tParamListKSP, &tParamListPCNone );
            tEigSolverSM->set_sublinear_solver_options( &tParamListKSP, &tParamListPCMumps );

            std::string tRHSMatType = std::string( "IdentityMat" );
            tEigProblem->set_rhs_matrix_type( tRHSMatType );

            // assemble jacobian
            tEigProblem->assemble_jacobian();

            moris::sint aIter = 0;

            // call solve
            tEigSolverLM->solve_linear_system( tEigProblem, aIter );
            //tEigSolverSM->solve_linear_system( tEigProblem, aIter );

            tEigProblem->get_matrix()->save_matrix_to_matlab_file( "JacobianCheck.txt" );

            Vector<real> const & tLargeEigenValues = std::dynamic_pointer_cast<Eigen_Solver_SLEPc>(tEigSolverLM)->get_eigenvalues();
            // Vector<real> const & tSmallEigenValues = std::dynamic_pointer_cast<Eigen_Solver_SLEPc>(tEigSolverSM)->get_eigenvalues();

            CHECK( equal_to( tLargeEigenValues(0), 60.864058, 1.0e+08 ) );
            // CHECK( equal_to( tSmallEigenValues(0), 0.0, 1.0e+08 ) );

            delete tEigProblem;
            delete tSolverInterface;
        }
    }


    TEST_CASE( "Eigen Solver SLEPC3v1", "[Eigen Solver SLEPC],[Eigen Solver], [EigSolveSLEPC]" )
    {
        if ( par_size() == 1 )
        {
            Solver_Interface_Proxy* tSolverInterface = new Solver_Interface_Proxy( 1 );
            tSolverInterface->mMyConstraintDofs.set_size( 0, 0 );

            Solver_Factory tSolFactory;

            Linear_Problem* tEigProblem = tSolFactory.create_linear_system( tSolverInterface, sol::MapType::Petsc, true );

            Parameter_List tLinearSolverParameterListLM = prm::create_slepc_algorithm_parameter_list();
            tLinearSolverParameterListLM.set( "Eigen_Algorithm", std::string( "krylovschur" ) );
            tLinearSolverParameterListLM.set( "Which", std::string( "LM" ) );
            tLinearSolverParameterListLM.set( "Num_Eig_Vals", 1 );    // 2; 1
            tLinearSolverParameterListLM.set( "STType", "shift" );    // 10 shift_invert
            tLinearSolverParameterListLM.set( "Update_Flag", false );


            Parameter_List tLinearSolverParameterListSM = prm::create_slepc_algorithm_parameter_list();
            tLinearSolverParameterListSM.set( "Eigen_Algorithm", std::string( "krylovschur" ) );
            tLinearSolverParameterListSM.set( "Which", std::string( "LM" ) );
            tLinearSolverParameterListSM.set( "Num_Eig_Vals", 1 );           // 2; 1
            tLinearSolverParameterListSM.set( "STType", "shift_invert" );    // 10 shift_invert
            tLinearSolverParameterListSM.set( "Update_Flag", false );

            std::shared_ptr< Linear_Solver_Algorithm > tEigSolverLM = tSolFactory.create_solver( tLinearSolverParameterListLM );
            std::shared_ptr< Linear_Solver_Algorithm > tEigSolverSM = tSolFactory.create_solver( tLinearSolverParameterListSM );

            Parameter_List tParamListKSP = prm::create_linear_algorithm_parameter_list_petsc();
            tParamListKSP.set( "KSPType", std::string( "gmres" ) );

            Parameter_List tParamListPCNone = prm::create_preconditioner_parameter_list( moris::sol::PreconditionerType::PETSC );
            tParamListPCNone.set( "PCType", std::string( "none" ) );
            Parameter_List tParamListPCMumps = prm::create_preconditioner_parameter_list( moris::sol::PreconditionerType::PETSC );
            tParamListPCMumps.set( "PCType", std::string( "none" ) );

            tEigSolverLM->set_sublinear_solver_options( &tParamListKSP, &tParamListPCNone );
            tEigSolverSM->set_sublinear_solver_options( &tParamListKSP, &tParamListPCMumps );

            std::string tRHSMatType = std::string( "IdentityMat" );
            tEigProblem->set_rhs_matrix_type( tRHSMatType );

            // assemble jacobian
            tEigProblem->assemble_jacobian();

            moris::sint aIter = 0;

            // call solve
            tEigSolverLM->solve_linear_system( tEigProblem, aIter );
            tEigSolverSM->solve_linear_system( tEigProblem, aIter );

            tEigProblem->get_matrix()->save_matrix_to_matlab_file( "JacobianCheck.txt" );

            Vector<real> const & tLargeEigenValues = std::dynamic_pointer_cast<Eigen_Solver_SLEPc>(tEigSolverLM)->get_eigenvalues();
            Vector<real> const & tSmallEigenValues = std::dynamic_pointer_cast<Eigen_Solver_SLEPc>(tEigSolverSM)->get_eigenvalues();


            CHECK( equal_to( tLargeEigenValues(0), 60.864058, 1.0e+08 ) );
            CHECK( equal_to( tSmallEigenValues(0), 0.0, 1.0e+08 ) );

            delete tEigProblem;
            delete tSolverInterface;
 
        }
    }
}    // namespace moris::dla