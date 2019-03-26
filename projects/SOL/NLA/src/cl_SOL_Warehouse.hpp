/*
 * cl_SOL_Warehouse.hpp
 *
 *  Created on: Okt 6, 2018
 *      Author: schmidt
 */
#ifndef MORIS_DISTLINALG_CL_SOL_WAREHOUSE_HPP_
#define MORIS_DISTLINALG_CL_SOL_WAREHOUSE_HPP_

// MORIS header files.
#include "typedefs.hpp" // CON/src
#include "cl_Cell.hpp"
#include <memory>
#include "cl_Param_List.hpp"
#include "cl_MSI_Dof_Type_Enums.hpp"

#include "cl_NLA_Nonlinear_Solver.hpp"

#include "cl_Map_Class.hpp"

namespace moris
{
namespace NLA
{
    class Nonlinear_Problem;
    class Nonlinear_Algorithm;
    class SOL_Warehouse
    {
    private:
        //! Pointer to the solver interface
        Solver_Interface * mSolverInterface;

    public:
        SOL_Warehouse( Solver_Interface * aSolverInterface ) : mSolverInterface( aSolverInterface ){};

        ~SOL_Warehouse(){};

        Solver_Interface * get_solver_interface(){ return mSolverInterface; };
    };
}
}
#endif /* MORIS_DISTLINALG_CL_SOL_WAREHOUSE_HPP_ */

