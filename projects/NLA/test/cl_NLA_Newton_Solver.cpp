/*
 * cl_NLA_Newton_Solver.cpp
 *
 *  Created on: Sep 21, 2018
 *      Author: schmidt
 */
#ifdef MORIS_HAVE_PARALLEL
 #include "Epetra_MpiComm.h"
 #include <mpi.h>
#endif

#include "catch.hpp"
#include "fn_equal_to.hpp"
#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Communication_Tools.hpp"

#include "cl_Solver_Factory.hpp" // DLA/src/
#include "cl_Solver_Input_Test.hpp" // DLA/src

#define protected public
#define private   public
#include "cl_NLA_Newton_Solver.hpp"
#undef protected
#undef private

namespace moris
{
    namespace NLA
    {
    TEST_CASE("Newton Solver Test 1","[NLA],[NLA_Test1]")
    {
        // Build Input Class
        Solver_Input* tSolverInput = new Solver_Input_Test( );

        // create solver factory
        Solver_Factory  tSolFactory;

        // create solver object
        std::shared_ptr< Linear_Solver > tLin = tSolFactory.create_solver( tSolverInput, SolverType::TRILINOSTEST );

        Newton_Solver tNewton;
        tNewton.mA = 55;

        tNewton.devide( 5 );
        Matrix<DDRMat> tA(1,1,1.0);
        CHECK( equal_to( tA(0,0), 1.0 ) );
        CHECK( equal_to( tNewton.mA, 11 ) );
    }
    }
}


