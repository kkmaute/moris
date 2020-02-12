/*
 * cl_TSA_Staggered_Time_Solver.hpp
 *
 *  Created on: Feb 02, 2019
 *      Author: schmidt
 */
#ifndef MORIS_DISTLINALG_CL_TSA_STAGGERED_TIME_SOLVER_HPP_
#define MORIS_DISTLINALG_CL_TSA_STAGGERED_TIME_SOLVER_HPP_

#include "cl_TSA_Time_Solver_Algorithm.hpp"
#include "cl_SOL_Dist_Vector.hpp"
#include "cl_DLA_Solver_Interface.hpp"

namespace moris
{
class Dist_Vector;
class Solver_Interface;

namespace tsa
{
    class Staggered_Time_Solver : public Time_Solver_Algorithm
    {
    private:

        void solve_staggered_time_system();

    public:
        //-------------------------------------------------------------------------------

        Staggered_Time_Solver(){};

        //-------------------------------------------------------------------------------

        ~Staggered_Time_Solver(){};

        //-------------------------------------------------------------------------------

        void solve( Dist_Vector * aFullVector );

        //-------------------------------------------------------------------------------

        void solve();

        //-------------------------------------------------------------------------------
    };
}
}

#endif /* MORIS_DISTLINALG_CL_TSA_STAGGERED_TIME_SOLVER_HPP_ */
