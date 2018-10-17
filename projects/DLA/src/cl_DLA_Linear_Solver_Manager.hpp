/*
 * cl_DLA_Linear_System_Manager.hpp
 *
 *  Created on: Okt 6, 2018
 *      Author: schmidt
 */
#ifndef MORIS_DISTLINALG_CL_DLA_Linear_Solver_MANAGER_HPP_
#define MORIS_DISTLINALG_CL_DLA_Linear_Solver_MANAGER_HPP_

// MORIS header files.
#ifdef MORIS_HAVE_PARALLEL
 #include <mpi.h>
#endif

#include "typedefs.hpp" // CON/src
#include "cl_Cell.hpp"
#include <memory>
//#include "cl_Param_List.hpp" // CON/src

namespace moris
{
namespace dla
{
    class Linear_Solver;
    class Linear_Solver_Manager
    {
    private:
        moris::Cell< std::shared_ptr< Linear_Solver > > mLinearSolverList;

        moris::uint mCallCounter = 0;

    protected:

    public:
        Linear_Solver_Manager();

        ~Linear_Solver_Manager();

        void set_linear_solver( std::shared_ptr< Linear_Solver > aLinSolver );

        void set_linear_solver( const moris::uint aListEntry,
                                      std::shared_ptr< Linear_Solver > aLinSolver );



        moris::Cell< std::shared_ptr< Linear_Solver > > get_solver_list(){ return mLinearSolverList; };

    };
}
}
#endif /* MORIS_DISTLINALG_CL_DLA_Linear_Solver_MANAGER_HPP_ */

