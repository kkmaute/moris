/*
 * cl_NLA_Newton_Solver.cpp
 *
 *  Created on: Sep 21, 2018
 *      Author: schmidt
 */
#ifdef MORIS_HAVE_PARALLEL
 #include <mpi.h>
#endif

#include "catch.hpp"
#include "fn_equal_to.hpp"
#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Communication_Tools.hpp"

#include "cl_Solver_Factory.hpp" // DLA/src/
#include "cl_Solver_Input_Test.cpp" // DLA/src
#include "cl_Vector.hpp" // DLA/src

#define protected public
#define private   public
#include "cl_NLA_Newton_Solver.hpp"
#undef protected
#undef private

#include "cl_NLA_Solver_Input_Test.cpp" // DLA/src

namespace moris
{
    namespace NLA
    {
    TEST_CASE("Newton Solver Test 1","[NLA],[NLA_Test1]")
    {
        Newton_Solver tNewton;

        // Build Input Class
        Solver_Input* tSolverInput = new NLA_Solver_Input_Test( &tNewton );

        // create solver factory
        Solver_Factory  tSolFactory;

        // create solver object
        std::shared_ptr< Linear_Solver > tLin = tSolFactory.create_solver( tSolverInput, SolverType::TRILINOSTEST );

        tNewton.set_linear_solver( tLin );

        tNewton.set_param("NLA_max_iter") = 10;
        tNewton.set_param("NLA_hard_break") = false;

        tSolverInput->set_test_problem();

        tNewton.solver_nonlinear_system();

        Matrix< DDSMat > tGlobalIndExtract( 2, 1, 0);
        tGlobalIndExtract( 1, 0 ) = 1;
        Matrix< DDRMat > tMyValues;

        tNewton.get_full_sol_vec()->extract_my_values( 2, tGlobalIndExtract, 0, tMyValues);

        CHECK( equal_to( tMyValues( 0, 0 ), 0.04011965, 1.0e+08 ) );
        CHECK( equal_to( tMyValues( 1, 0 ), 0.0154803, 1.0e+08 ) );
    }
    }
}


