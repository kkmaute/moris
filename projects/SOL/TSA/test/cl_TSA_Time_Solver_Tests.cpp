/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_TSA_Time_Solver_Tests.cpp
 *
 */

#include "catch.hpp"
#include "fn_equal_to.hpp"
#include "moris_typedefs.hpp"
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
#include "cl_NLA_Nonlinear_Algorithm.hpp"
#include "cl_NLA_Nonlinear_Solver_Factory.hpp"

#include "cl_DLA_Solver_Factory.hpp"
#include "cl_DLA_Linear_Solver_Aztec.hpp"
#include "cl_DLA_Linear_Solver.hpp"

namespace moris::tsa
{
    TEST_CASE( "TimeSolverRest", "[TSA],[TimeSolver]" )
    {
        if ( par_size() == 1 )
        {
            Parameter_List tTimeSolverParameterList = prm::create_time_solver_algorithm_parameter_list();
            tTimeSolverParameterList.set( "TSA_Num_Time_Steps", 1000 );
            tTimeSolverParameterList.set( "TSA_Time_Frame", 10.0 );
            std::shared_ptr< Time_Solver_Algorithm > tTimesolverAlgorithm = std::make_shared< Monolithic_Time_Solver >(
                    tTimeSolverParameterList );

            // Create solver interface
            Solver_Interface *tSolverInput = new TSA_Solver_Interface_Proxy();

            dla::Solver_Factory tSolFactory;
            Parameter_List      tLinSolverParameterList = prm::create_linear_algorithm_parameter_list_aztec();
            tLinSolverParameterList.set( "AZ_diagnostics", AZ_none );
            tLinSolverParameterList.set( "AZ_output", AZ_none );
            tLinSolverParameterList.set( "AZ_solver", AZ_gmres );
            tLinSolverParameterList.set( "AZ_precond", AZ_dom_decomp );
            std::shared_ptr< dla::Linear_Solver_Algorithm > tLinSolverAlgorithm = tSolFactory.create_solver( tLinSolverParameterList );

            auto tLinSolManager = new dla::Linear_Solver();
            tLinSolManager->set_linear_algorithm( 0, tLinSolverAlgorithm );

            NLA::Nonlinear_Solver_Factory               tNonlinFactory;
            std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlLinSolverAlgorithm = tNonlinFactory.create_nonlinear_solver();
            tNonlLinSolverAlgorithm->set_linear_solver( tLinSolManager );

            NLA::Nonlinear_Solver tNonlinearSolverManager;
            tNonlinearSolverManager.set_nonlinear_algorithm( tNonlLinSolverAlgorithm, 0 );

            Vector< enum MSI::Dof_Type > tDofTypes( 1 );
            tDofTypes( 0 ) = MSI::Dof_Type::TEMP;
            tNonlinearSolverManager.set_dof_type_list( tDofTypes );

            tTimesolverAlgorithm->set_nonlinear_solver( &tNonlinearSolverManager );

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

            delete ( tLinSolManager );
            delete ( tSolverInput );
        }
    }
    }
