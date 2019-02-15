/*
 * cl_TSA_Monolithic_Time_Solver.hpp
 *
 *  Created on: Feb 02, 2019
 *      Author: schmidt
 */
#ifndef MORIS_DISTLINALG_CL_TSA_MONOLITHIC_TIME_SOLVER_HPP_
#define MORIS_DISTLINALG_CL_TSA_MONOLITHIC_TIME_SOLVER_HPP_

#include "cl_TSA_Time_Solver_Algorithm.hpp"
#include "cl_Vector.hpp"

#include "cl_DLA_Solver_Interface.hpp"
#include "cl_Matrix_Vector_Factory.hpp"

namespace moris
{
class Dist_Vector;
class Solver_Interface;

namespace tsa
{
    class Monolithic_Time_Solver : public Time_Solver_Algorithm
    {
    private:

        void solve_monolytic_time_system();

        void perform_mapping();

    public:
        //-------------------------------------------------------------------------------

        Monolithic_Time_Solver( )
        {};

        //-------------------------------------------------------------------------------

        ~Monolithic_Time_Solver(){};

        //-------------------------------------------------------------------------------

        void solve( Dist_Vector * aFullVector );

        //-------------------------------------------------------------------------------

        void solve();

        //-------------------------------------------------------------------------------
    };
}
}
#endif /* MORIS_DISTLINALG_CL_TSA_MONOLITHIC_TIME_SOLVER_HPP_ */
