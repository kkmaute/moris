/*
 * cl_NLA_NonlinearDatabase.cpp
 *
 *  Created on: Jan 21, 2018
 *      Author: schmidt
 */
#ifdef MORIS_HAVE_PARALLEL
 #include <mpi.h>
#endif

#include "catch.hpp"
#include "fn_equal_to.hpp"
#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "fn_reshape.hpp"

#include "linalg_typedefs.hpp"
#include "cl_Communication_Tools.hpp"
#include "cl_Communication_Tools.hpp"

#include "cl_MSI_Dof_Type_Enums.hpp"

#define protected public
#define private   public
#include "cl_NLA_NonlinearDatabase.hpp"
#include "cl_NLA_Nonlinear_Solver_Manager.hpp"
#include "cl_NLA_Nonlinear_Manager.hpp"
#include "cl_NLA_Solver_Interface_Proxy2.hpp"
#undef protected
#undef private

namespace moris
{
namespace NLA
{
    TEST_CASE("NonlinearDatabase","[NLA],[NLA_Database]")
    {
        moris::Cell< enum MSI::Dof_Type > tDofTypes1( 2 );
        moris::Cell< enum MSI::Dof_Type > tDofTypes2( 1 );

        tDofTypes1( 0 ) = MSI::Dof_Type::UX;
        tDofTypes1( 1 ) = MSI::Dof_Type::UY;

        tDofTypes2( 0 ) = MSI::Dof_Type::TEMP;

        NonLinDatabase tNonlinearDatabase;

        tNonlinearDatabase.set_dof_type_list( tDofTypes1, 0 );
        tNonlinearDatabase.set_dof_type_list( tDofTypes2, 1 );

        CHECK( equal_to( tNonlinearDatabase.get_dof_type_list_level()( 0 ), 0 ) );
        CHECK( equal_to( tNonlinearDatabase.get_dof_type_list_level()( 1 ), 1 ) );

        CHECK( equal_to( static_cast< int >( tNonlinearDatabase.get_dof_type_list()( 0 )( 0 ) ), static_cast< int >( MSI::Dof_Type::UX ) ) );
        CHECK( equal_to( static_cast< int >( tNonlinearDatabase.get_dof_type_list()( 0 )( 1 ) ), static_cast< int >( MSI::Dof_Type::UY ) ) );
        CHECK( equal_to( static_cast< int >( tNonlinearDatabase.get_dof_type_list()( 1 )( 0 ) ), static_cast< int >( MSI::Dof_Type::TEMP ) ) );

        CHECK( equal_to( tNonlinearDatabase.get_dof_type_list().size() , 2 ) );
        CHECK( equal_to( tNonlinearDatabase.get_dof_type_list().capacity() , 2 ) );

        CHECK( equal_to( static_cast< int >( tNonlinearDatabase.get_dof_type_union()( 0 ) ), static_cast< int >( MSI::Dof_Type::UX ) ) );
        CHECK( equal_to( static_cast< int >( tNonlinearDatabase.get_dof_type_union()( 1 ) ), static_cast< int >( MSI::Dof_Type::UY ) ) );
        CHECK( equal_to( static_cast< int >( tNonlinearDatabase.get_dof_type_union()( 2 ) ), static_cast< int >( MSI::Dof_Type::TEMP ) ) );

    }

    TEST_CASE("NonlinearDatabase2","[NLA],[NLA_Database2]")
    {
        moris::Cell< enum MSI::Dof_Type > tDofTypes1( 2 );
        moris::Cell< enum MSI::Dof_Type > tDofTypes2( 1 );

        tDofTypes1( 0 ) = MSI::Dof_Type::UX;
        tDofTypes1( 1 ) = MSI::Dof_Type::UY;

        tDofTypes2( 0 ) = MSI::Dof_Type::TEMP;

        Nonlinear_Solver_Manager tNonlinearSolverManager;

        tNonlinearSolverManager.set_dof_type_list( tDofTypes1 );
        tNonlinearSolverManager.set_dof_type_list( tDofTypes2 );

        CHECK( equal_to( static_cast< int >( tNonlinearSolverManager.get_dof_type_list()( 0 )( 0 ) ), static_cast< int >( MSI::Dof_Type::UX ) ) );
        CHECK( equal_to( static_cast< int >( tNonlinearSolverManager.get_dof_type_list()( 0 )( 1 ) ), static_cast< int >( MSI::Dof_Type::UY ) ) );
        CHECK( equal_to( static_cast< int >( tNonlinearSolverManager.get_dof_type_list()( 1 )( 0 ) ), static_cast< int >( MSI::Dof_Type::TEMP ) ) );
    }

    TEST_CASE("NonlinearDatabase3","[NLA],[NLA_Database3]")
    {
        moris::Cell< enum MSI::Dof_Type > tDofTypes1( 2 );
        moris::Cell< enum MSI::Dof_Type > tDofTypes2( 1 );

        tDofTypes1( 0 ) = MSI::Dof_Type::UX;
        tDofTypes1( 1 ) = MSI::Dof_Type::UY;

        tDofTypes2( 0 ) = MSI::Dof_Type::TEMP;



        Nonlinear_Solver_Manager tNonlinearSolverManager1( NLA::NonlinearSolverType::NLBGS_SOLVER );
        Nonlinear_Solver_Manager tNonlinearSolverManager2( NLA::NonlinearSolverType::NEWTON_SOLVER );
        Nonlinear_Solver_Manager tNonlinearSolverManager3( NLA::NonlinearSolverType::NEWTON_SOLVER );
        tNonlinearSolverManager1.set_dof_type_list( tDofTypes2 );
        tNonlinearSolverManager1.set_dof_type_list( tDofTypes1 );
        tNonlinearSolverManager2.set_dof_type_list( tDofTypes1 );
        tNonlinearSolverManager3.set_dof_type_list( tDofTypes2 );

        CHECK( equal_to( static_cast< int >( tNonlinearSolverManager1.get_dof_type_list()( 1 )( 0 ) ), static_cast< int >( MSI::Dof_Type::UX ) ) );
        CHECK( equal_to( static_cast< int >( tNonlinearSolverManager1.get_dof_type_list()( 1 )( 1 ) ), static_cast< int >( MSI::Dof_Type::UY ) ) );
        CHECK( equal_to( static_cast< int >( tNonlinearSolverManager1.get_dof_type_list()( 0 )( 0 ) ), static_cast< int >( MSI::Dof_Type::TEMP ) ) );

        Solver_Interface * tSolverInput = new NLA_Solver_Interface_Proxy_II();

        Nonlinear_Manager tNonlinearDatabase( tSolverInput );

        tNonlinearDatabase.set_nonliner_solver_managers( & tNonlinearSolverManager1 );
        tNonlinearDatabase.set_nonliner_solver_managers( & tNonlinearSolverManager2 );
        tNonlinearDatabase.set_nonliner_solver_managers( & tNonlinearSolverManager3 );

        tNonlinearDatabase.solve();
    }
}
}

