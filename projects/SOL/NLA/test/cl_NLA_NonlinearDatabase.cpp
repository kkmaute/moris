/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_NLA_NonlinearDatabase.cpp
 *
 */

#include "catch.hpp"
#include "fn_equal_to.hpp"
#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Communication_Tools.hpp"

#include "cl_MSI_Dof_Type_Enums.hpp"
#include "cl_SOL_Dist_Vector.hpp"

#define protected public
#define private   public
#include "cl_NLA_Nonlinear_Solver.hpp"
#include "cl_SOL_Warehouse.hpp"
#include "cl_NLA_Solver_Interface_Proxy2.hpp"
#undef protected
#undef private

#include "cl_SOL_Matrix_Vector_Factory.hpp"

namespace moris
{

namespace NLA
{
TEST_CASE("NonlinearDatabase3","[NLA],[NLA_Database3]")
{
    //FIXME Mathias commented this test out 05.30.2020 - I will fix it soon
//    if ( par_size() == 1 )
//    {
//        // Create and fill dof type lists
//        Vector< enum MSI::Dof_Type > tDofTypes1( 2 );
//        Vector< enum MSI::Dof_Type > tDofTypes2( 1 );
//
//        tDofTypes1( 0 ) = MSI::Dof_Type::UX;
//        tDofTypes1( 1 ) = MSI::Dof_Type::UY;
//        tDofTypes2( 0 ) = MSI::Dof_Type::TEMP;
//
//        // Create nonlinear solver manager
//        Nonlinear_Solver tNonlinearSolverManager1( NLA::NonlinearSolverType::NLBGS_SOLVER );
//        Nonlinear_Solver tNonlinearSolverManager2( NLA::NonlinearSolverType::NEWTON_SOLVER );
//        Nonlinear_Solver tNonlinearSolverManager3( NLA::NonlinearSolverType::NEWTON_SOLVER );
//
//        // Set dof type lists to nonlinear solver manager
//        tNonlinearSolverManager1.set_dof_type_list( tDofTypes2 );
//        tNonlinearSolverManager1.set_dof_type_list( tDofTypes1 );
//        tNonlinearSolverManager2.set_dof_type_list( tDofTypes1 );
//        tNonlinearSolverManager3.set_dof_type_list( tDofTypes2 );
//
//        tNonlinearSolverManager1.set_sub_nonlinear_solver( &tNonlinearSolverManager2 );
//        tNonlinearSolverManager1.set_sub_nonlinear_solver( &tNonlinearSolverManager3 );
//
//        // Create solver interface
//        Solver_Interface * tSolverInput = new NLA_Solver_Interface_Proxy_II();
//
//        // Build matrix vector factory
//        sol::Matrix_Vector_Factory    tMatFactory( sol::MapType::Epetra );
//
//        sol::Dist_Map*  tMap = tMatFactory.create_map( tSolverInput->get_my_local_global_overlapping_map() );
//
//        // Create Full Vector
//        sol::Dist_Vector * tFullVector = tMatFactory.create_vector( tSolverInput, tMap, 1 );
//
//        tSolverInput->set_solution_vector( tFullVector );
//
//        // Initilaze full vector with zeros
//        tFullVector->vec_put_scalar( 0.0 );
//
//        // Create solver database
//        sol::SOL_Warehouse tSolverWarehouse( tSolverInput );
//
//        tNonlinearSolverManager1.set_solver_warehouse( &tSolverWarehouse );
//        tNonlinearSolverManager2.set_solver_warehouse( &tSolverWarehouse );
//        tNonlinearSolverManager3.set_solver_warehouse( &tSolverWarehouse );
//
//        // Solve
//        tNonlinearSolverManager1.solve( tFullVector );
//
//        Matrix< DDRMat > tSol;
//        tNonlinearSolverManager1.get_full_solution( tSol );
//
//        print(tSol,"tSol");
//
//        CHECK( equal_to( tSol( 0, 0 ), 0.03510531645, 1.0e+08 ) );
//        CHECK( equal_to( tSol( 1, 0 ), 0.011710521925, 1.0e+08 ) );
//        CHECK( equal_to( tSol( 2, 0 ), 0.036574625191, 1.0e+08 ) );
//        CHECK( equal_to( tSol( 3, 0 ), 0.013057249537, 1.0e+08 ) );
//    }
}

    TEST_CASE("NonlinearDatabase2","[NLA],[NLA_Database2]")
    {
        if ( par_size() == 1 )
        {
        Vector< enum MSI::Dof_Type > tDofTypes1( 2 );
        Vector< enum MSI::Dof_Type > tDofTypes2( 1 );

        tDofTypes1( 0 ) = MSI::Dof_Type::UX;
        tDofTypes1( 1 ) = MSI::Dof_Type::UY;

        tDofTypes2( 0 ) = MSI::Dof_Type::TEMP;

        Nonlinear_Solver tNonlinearSolverManager;

        tNonlinearSolverManager.set_dof_type_list( tDofTypes1 );
        tNonlinearSolverManager.set_dof_type_list( tDofTypes2 );

        CHECK( equal_to( static_cast< int >( tNonlinearSolverManager.get_dof_type_list()( 0 )( 0 ) ), static_cast< int >( MSI::Dof_Type::UX ) ) );
        CHECK( equal_to( static_cast< int >( tNonlinearSolverManager.get_dof_type_list()( 0 )( 1 ) ), static_cast< int >( MSI::Dof_Type::UY ) ) );
        CHECK( equal_to( static_cast< int >( tNonlinearSolverManager.get_dof_type_list()( 1 )( 0 ) ), static_cast< int >( MSI::Dof_Type::TEMP ) ) );
        }
    }

//    TEST_CASE("NonlinearDatabase4","[NLA],[NLA_Database4]")
//    {
//        if ( par_size() == 1 )
//        {
//        Vector< enum MSI::Dof_Type > tDofTypes1( 1 );
//        Vector< enum MSI::Dof_Type > tDofTypes2( 3 );
//        Vector< enum MSI::Dof_Type > tDofTypes3( 2 );
//        Vector< enum MSI::Dof_Type > tDofTypes4( 1 );
//
//
//        tDofTypes1( 0 ) = MSI::Dof_Type::TEMP;
//
//        tDofTypes2( 0 ) = MSI::Dof_Type::UX;
//        tDofTypes2( 1 ) = MSI::Dof_Type::UY;
//        tDofTypes2( 2 ) = MSI::Dof_Type::UZ;
//
//        tDofTypes3( 0 ) = MSI::Dof_Type::UX;
//        tDofTypes3( 1 ) = MSI::Dof_Type::UY;
//
//        tDofTypes4( 0 ) = MSI::Dof_Type::UZ;
//
//        Nonlinear_Solver tNonlinearSolverManager1( NLA::NonlinearSolverType::NLBGS_SOLVER );
//        Nonlinear_Solver tNonlinearSolverManager2( NLA::NonlinearSolverType::NEWTON_SOLVER );
//        Nonlinear_Solver tNonlinearSolverManager3( NLA::NonlinearSolverType::NLBGS_SOLVER );
//        Nonlinear_Solver tNonlinearSolverManager4( NLA::NonlinearSolverType::NEWTON_SOLVER );
//        Nonlinear_Solver tNonlinearSolverManager5( NLA::NonlinearSolverType::NEWTON_SOLVER );
//
//        tNonlinearSolverManager1.set_dof_type_list( tDofTypes1 );
//        tNonlinearSolverManager1.set_dof_type_list( tDofTypes2 );
//
//        tNonlinearSolverManager2.set_dof_type_list( tDofTypes1 );
//
//        tNonlinearSolverManager3.set_dof_type_list( tDofTypes3 );
//        tNonlinearSolverManager3.set_dof_type_list( tDofTypes4 );
//
//        tNonlinearSolverManager4.set_dof_type_list( tDofTypes3 );
//
//        tNonlinearSolverManager5.set_dof_type_list( tDofTypes4 );
//
//        Solver_Interface * tSolverInput = new NLA_Solver_Interface_Proxy_II();
//
//        SOL_Warehouse tSolverWarehouse( tSolverInput );
//
//        tSolverWarehouse.set_nonlinear_solver_managers( & tNonlinearSolverManager1 );
//        tSolverWarehouse.set_nonlinear_solver_managers( & tNonlinearSolverManager2 );
//        tSolverWarehouse.set_nonlinear_solver_managers( & tNonlinearSolverManager3 );
//        tSolverWarehouse.set_nonlinear_solver_managers( & tNonlinearSolverManager4 );
//        tSolverWarehouse.set_nonlinear_solver_managers( & tNonlinearSolverManager5 );
//
//        tSolverWarehouse.create_solver_manager_dependencies();
//
//        CHECK( equal_to( tSolverWarehouse.get_nonlinear_solver_manager_index( 0, 0 ),  1 ) );
//        CHECK( equal_to( tSolverWarehouse.get_nonlinear_solver_manager_index( 0, 1 ),  2 ) );
//        CHECK( equal_to( tSolverWarehouse.get_nonlinear_solver_manager_index( 2, 0 ),  3 ) );
//        CHECK( equal_to( tSolverWarehouse.get_nonlinear_solver_manager_index( 2, 1 ),  4 ) );
//        }
//    }

}
}

