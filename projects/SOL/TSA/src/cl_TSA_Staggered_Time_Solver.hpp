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
        /**
         * @brief default constructor
         *
         */
        Staggered_Time_Solver(){};

        //-------------------------------------------------------------------------------
        /**
         * @brief Constructor using a given parameter list
         *
         * @param[in] rSolverDatabase Poiner to the solver database
         */
        Staggered_Time_Solver( const ParameterList aParameterlist ) : Time_Solver_Algorithm( aParameterlist )
        {};

        //-------------------------------------------------------------------------------

        ~Staggered_Time_Solver(){};

        //-------------------------------------------------------------------------------
        /**
         * @brief Solve call using a given soltion vector
         *
         * @param[in] aFullVector     Solution Vector
         */
        void solve( Dist_Vector * aFullVector );

        //-------------------------------------------------------------------------------
        /**
         * @brief Solve call
         *
         * @param[in] aFullVector
         */
        void solve();

        //-------------------------------------------------------------------------------
    };
}
}

#endif /* MORIS_DISTLINALG_CL_TSA_STAGGERED_TIME_SOLVER_HPP_ */
