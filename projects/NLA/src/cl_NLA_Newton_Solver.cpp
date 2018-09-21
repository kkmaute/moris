/*
 * cl_NLA_Newton_Solver.cpp
 *
 *  Created on: Sep 21, 2018
 *      Author: schmidt
 */
#include "cl_NLA_Newton_Solver.hpp"
#include "cl_Linear_Solver.hpp"
#include "cl_Solver_Input.hpp"
#include "cl_DistLinAlg_Enums.hpp"


namespace moris
{
    namespace NLA
    {
    Newton_Solver::Newton_Solver( std::shared_ptr< Linear_Solver > aLinearSolver ): mLinearSolver( aLinearSolver )
    {
        Solver_Input * tInput = mLinearSolver->get_solver_input();

        Matrix_Vector_Factory    tMatFactory;

        // create map object
        mMap = tMatFactory.create_map( tInput->get_num_my_dofs(),
                                       tInput->get_my_local_global_map(),
                                       tInput->get_constr_dof(),
                                       tInput->get_my_local_global_overlapping_map());

        // Build RHS/LHS vector
        mVectorMaster = tMatFactory.create_vector( tInput, mMap, VectorType::FREE );

    }

    }
}
