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

#include "linalg_typedefs.hpp"
#include "cl_Communication_Tools.hpp"

#include "cl_MSI_Dof_Type_Enums.hpp"

#define protected public
#define private   public
#include "cl_TSA_Monolithic_Time_Solver.hpp"
#include "cl_TSA_Solver_Interface_Proxy.hpp"
#include "cl_TSA_Time_Solver.hpp"
#undef protected
#undef private

#include "cl_SOL_Warehouse.hpp"
#include "cl_NLA_Nonlinear_Solver.hpp"

namespace moris
{
    namespace tsa
    {
        TEST_CASE("TimeSolverRest","[TSA],[TimeSolver]")
            {
            if ( par_size() == 1 )
            {
                std::shared_ptr< Time_Solver_Algorithm > tTimesolverAlgorithm = std::make_shared< Monolithic_Time_Solver >();

                // Create solver interface
                Solver_Interface * tSolverInput = new TSA_Solver_Interface_Proxy();

                NLA::Nonlinear_Solver tNonlinearSolverManager( NLA::NonlinearSolverType::NEWTON_SOLVER );

                moris::Cell< enum MSI::Dof_Type > tDofTypes( 1 );
                tDofTypes( 0 ) = MSI::Dof_Type::TEMP;
                tNonlinearSolverManager.set_dof_type_list( tDofTypes );

                tTimesolverAlgorithm->set_nonlinear_solver( & tNonlinearSolverManager );

                tTimesolverAlgorithm->set_param("TSA_Num_Time_Steps")   = 1000;
                tTimesolverAlgorithm->set_param("TSA_Time_Frame")       = 10.0;

                Time_Solver tTimeSolver;

                tTimeSolver.set_time_solver_algorithm( tTimesolverAlgorithm );

                sol::SOL_Warehouse tSolverWarehouse( tSolverInput );

                tNonlinearSolverManager.set_solver_warehouse( &tSolverWarehouse );
                tTimeSolver.set_solver_warehouse( &tSolverWarehouse );

                tTimeSolver.set_dof_type_list( tDofTypes );

                tTimeSolver.solve();

                Matrix< DDRMat > tSol;
                tTimeSolver.get_full_solution( tSol );

                CHECK( equal_to( tSol( 0, 0 ), -8.869937049794211e-01, 1.0e+08 ) );
            }
            }
    }
}

