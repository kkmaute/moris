/*
 * cl_TSA_Time_Solver.hpp
 *
 *  Created on: Feb 02, 2019
 *      Author: schmidt
 */
#ifndef MORIS_DISTLINALG_CL_TSA_TIME_SOLVER_HPP_
#define MORIS_DISTLINALG_CL_TSA_TIME_SOLVER_HPP_

#include <iostream>

#include "cl_NLA_Nonlinear_Database.hpp"

// MORIS header files.

namespace moris
{
class Map_Class;
class Dist_Vector;
class Solver_Interface;

namespace tsa
{
    class Time_Solver
    {
    private:

    protected:
        //! Pointer to database
        NLA::Nonlinear_Database * mDatabase = nullptr;

        //! Pointer to solver interface
        Solver_Interface * mSolverInterface = nullptr;

    public:
        //-------------------------------------------------------------------------------

        Time_Solver();

        //-------------------------------------------------------------------------------

        ~Time_Solver(){};

        //-------------------------------------------------------------------------------

        virtual void solve(){};
        virtual void finalize(){};

        void set_database( NLA::Nonlinear_Database * aDatabase )
        {
            mDatabase = aDatabase;
        };

    };
}
}
#endif /* MORIS_DISTLINALG_CL_TSA_TIME_SOLVER_HPP_ */
