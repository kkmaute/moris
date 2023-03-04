/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_TSA_Staggered_Time_Solver.hpp
 *
 */

#ifndef MORIS_DISTLINALG_CL_TSA_STAGGERED_TIME_SOLVER_HPP_
#define MORIS_DISTLINALG_CL_TSA_STAGGERED_TIME_SOLVER_HPP_

#include "cl_TSA_Time_Solver_Algorithm.hpp"

namespace moris
{
    namespace sol
    {
        class Dist_Vector;
    }
class Solver_Interface;

namespace tsa
{
    class Staggered_Time_Solver : public Time_Solver_Algorithm
    {
    private:

        void solve_staggered_time_system( moris::Cell< sol::Dist_Vector * > & aFullVector );

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
         * @param[in] rSolverDatabase Pointer to the solver database
         */
        Staggered_Time_Solver( const ParameterList aParameterlist ) : Time_Solver_Algorithm( aParameterlist )
        {};

        //-------------------------------------------------------------------------------

//        ~Staggered_Time_Solver(){};

        //-------------------------------------------------------------------------------
        /**
         * @brief Solve call using a given solution vector
         *
         * @param[in] aFullVector     Solution Vector
         */
        void solve( moris::Cell< sol::Dist_Vector * > & aFullVector );

        //-------------------------------------------------------------------------------
    };
}
}

#endif /* MORIS_DISTLINALG_CL_TSA_STAGGERED_TIME_SOLVER_HPP_ */

