/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_TSA_Monolithic_Time_Solver.hpp
 *
 */

#ifndef MORIS_DISTLINALG_CL_TSA_MONOLITHIC_TIME_SOLVER_HPP_
#define MORIS_DISTLINALG_CL_TSA_MONOLITHIC_TIME_SOLVER_HPP_

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
        class Monolithic_Time_Solver : public Time_Solver_Algorithm
        {
          private:
            void solve_monolithic_time_system( Vector< sol::Dist_Vector* >& aFullVector );

            void solve_implicit_DqDs( Vector< sol::Dist_Vector* >& aFullAdjointVector );

            moris::real mLambdaInc = 0;

          public:
            //-------------------------------------------------------------------------------
            /**
             * @brief default constructor
             *
             * @param[in] rSolverDatabase Pointer to the solver database
             */
            Monolithic_Time_Solver(){};

            //-------------------------------------------------------------------------------
            /**
             * @brief Constructor using a given parameter list
             *
             * @param[in] aParameterlist     User defined parameter list
             */
            Monolithic_Time_Solver( const Parameter_List aParameterlist )
                    : Time_Solver_Algorithm( aParameterlist ){};

            //-------------------------------------------------------------------------------

            //        ~Monolithic_Time_Solver(){};

            //-------------------------------------------------------------------------------
            /**
             * @brief Solve call using a given solution vector
             *
             * @param[in] aFullVector     Solution Vector
             */
            void solve( Vector< sol::Dist_Vector* >& aFullVector );

            //-------------------------------------------------------------------------------

            void set_lambda_increment( moris::real aLambdaInc );

            //-------------------------------------------------------------------------------

            moris::real get_new_lambda();
        };
    }    // namespace tsa
}    // namespace moris
#endif /* MORIS_DISTLINALG_CL_TSA_MONOLITHIC_TIME_SOLVER_HPP_ */
