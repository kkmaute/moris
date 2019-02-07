/*
 * cl_TSA_Time_Solver_Test.cpp
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
//#include "cl_TSA_Time_Solver.hpp"
#include "cl_TSA_Monolithic_Time_Solver.hpp"
#include "cl_TSA_Solver_Interface_Proxy.hpp"
#undef protected
#undef private

namespace moris
{

namespace tsa
{
    TEST_CASE("TimeSolverRest","[TSA],[TimeSolver]")
    {
        Time_Solver * tTimesolver = new Monolithic_Time_Solver();

        // Create solver interface
        Solver_Interface * tSolverInput = new TSA_Solver_Interface_Proxy();

        // Create solver database
        NLA::Nonlinear_Database tNonlinearDatabase( tSolverInput );

        NLA::Nonlinear_Solver tNonlinearSolverManager( NLA::NonlinearSolverType::NEWTON_SOLVER );

        moris::Cell< enum MSI::Dof_Type > tDofTypes( 1 );
        tDofTypes( 0 ) = MSI::Dof_Type::TEMP;
        tNonlinearSolverManager.set_dof_type_list( tDofTypes );

        tNonlinearDatabase.set_nonliner_solver_managers( & tNonlinearSolverManager );

        tTimesolver->set_database( &tNonlinearDatabase );

        tNonlinearDatabase.finalize();

        tTimesolver -> finalize();

        tTimesolver -> solve();

        Matrix< DDRMat > tSol;
        tTimesolver ->get_full_solution( tSol );

        CHECK( equal_to( tSol( 0, 0 ), -8.869937049794211e-01, 1.0e+08 ) );

    }
}
}

