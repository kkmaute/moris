/*
 * cl_TSA_Time_Solver_Test.cpp
 *
 *  Created on: Jan 21, 2018
 *      Author: schmidt
 */
#include "catch.hpp"
#include "fn_equal_to.hpp"
#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "fn_reshape.hpp"

#include "linalg_typedefs.hpp"
#include "cl_Communication_Tools.hpp"

#include "cl_MSI_Dof_Type_Enums.hpp"

#define protected public
#define private   public
#include "cl_TSA_Staggered_Time_Solver.hpp"
#include "cl_TSA_Monolithic_Time_Solver.hpp"
#include "cl_TSA_Solver_Interface_Proxy2.hpp"
#undef protected
#undef private

namespace moris
{

namespace tsa
{
    TEST_CASE("TimeSolverStaggered","[TSA],[TimeSolverStaggered]")
    {
        Time_Solver * tTimesolver1 = new Staggered_Time_Solver();
        Time_Solver * tTimesolver2 = new Monolithic_Time_Solver();
        Time_Solver * tTimesolver3 = new Monolithic_Time_Solver();

        // Create solver interface
        Solver_Interface * tSolverInput = new TSA_Solver_Interface_Proxy_II();

        // Create solver database
        NLA::Nonlinear_Database tNonlinearDatabase( tSolverInput );

        NLA::Nonlinear_Solver tNonlinearSolverManager1( NLA::NonlinearSolverType::NEWTON_SOLVER );
        NLA::Nonlinear_Solver tNonlinearSolverManager2( NLA::NonlinearSolverType::NEWTON_SOLVER );

        moris::Cell< enum MSI::Dof_Type > tDofTypes1( 1 );
        tDofTypes1( 0 ) = MSI::Dof_Type::TEMP;

        moris::Cell< enum MSI::Dof_Type > tDofTypes2( 1 );
        tDofTypes2( 0 ) = MSI::Dof_Type::UX;

        tNonlinearSolverManager1.set_dof_type_list( tDofTypes1 );
        tNonlinearSolverManager2.set_dof_type_list( tDofTypes2 );

        tTimesolver1->set_dof_type_list( tDofTypes1 );
        tTimesolver1->set_dof_type_list( tDofTypes2 );
        tTimesolver2->set_dof_type_list( tDofTypes1 );
        tTimesolver3->set_dof_type_list( tDofTypes2 );

        //------------------------------------------------------------------------------------------

        //tNonlinearSolverManager.set_dof_type_list( tDofTypes );

        tNonlinearDatabase.set_nonliner_solver_managers( & tNonlinearSolverManager1 );
        tNonlinearDatabase.set_nonliner_solver_managers( & tNonlinearSolverManager2 );

        tTimesolver2->set_nonlinear_solver( & tNonlinearSolverManager1 );
        tTimesolver3->set_nonlinear_solver( & tNonlinearSolverManager2 );

        tTimesolver1->set_time_solver( tTimesolver2 );
        tTimesolver1->set_time_solver( tTimesolver3 );

        tTimesolver1->set_database( &tNonlinearDatabase );
        tTimesolver2->set_database( &tNonlinearDatabase );
        tTimesolver3->set_database( &tNonlinearDatabase );

        tNonlinearDatabase.finalize();

        tTimesolver1 -> solve();

        Matrix< DDRMat > tSol;
        tTimesolver1 ->get_full_solution( tSol );

        print( tSol, "tSol");
        CHECK( equal_to( tSol( 0, 0 ), -8.869937049794211e-01, 1.0e+08 ) );
        CHECK( equal_to( tSol( 1, 0 ), -8.869937049794211e-01, 1.0e+08 ) );
    }
}
}

