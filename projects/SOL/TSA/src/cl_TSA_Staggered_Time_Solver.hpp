/*
 * cl_TSA_Staggered_Time_Solver.hpp
 *
 *  Created on: Feb 02, 2019
 *      Author: schmidt
 */
#ifndef MORIS_DISTLINALG_CL_TSA_STAGGERED_TIME_SOLVER_HPP_
#define MORIS_DISTLINALG_CL_TSA_STAGGERED_TIME_SOLVER_HPP_

#include <iostream>

// MORIS header files.

namespace moris
{
class Map_Class;
class Dist_Vector;
class Solver_Interface;

namespace tsa
{
    class Staggered_Time_Solver : public Time_Solver
    {
    private:

    protected:
        //! Pointer to my nonlinear solver manager
        //Nonlinear_Solver * mMyNonLinSolverManager = nullptr;



    public:
        //-------------------------------------------------------------------------------

        Monolithic_Time_Solver();

        //-------------------------------------------------------------------------------

        ~Monolithic_Time_Solver();

        //-------------------------------------------------------------------------------

        void solve(){};

    };
}
}
#endif /* MORIS_DISTLINALG_CL_TSA_STAGGERED_TIME_SOLVER_HPP_ */
