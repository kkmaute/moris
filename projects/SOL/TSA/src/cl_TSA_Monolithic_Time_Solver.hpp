/*
 * cl_TSA_Monolithic_Time_Solver.hpp
 *
 *  Created on: Feb 02, 2019
 *      Author: schmidt
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

        void solve_monolytic_time_system();

        moris::real mLambdaInc = 0;

    public:
        //-------------------------------------------------------------------------------
        /**
         * @brief default constructor
         *
         * @param[in] rSolverDatabase Poiner to the solver database
         */
        Monolithic_Time_Solver( )
        {};

        //-------------------------------------------------------------------------------
        /**
         * @brief Constructor using a given parameter list
         *
         * @param[in] aParameterlist     User defined parameter list
         */
        Monolithic_Time_Solver( const ParameterList aParameterlist ) : Time_Solver_Algorithm( aParameterlist )
        {};

        //-------------------------------------------------------------------------------

//        ~Monolithic_Time_Solver(){};

        //-------------------------------------------------------------------------------
        /**
         * @brief Solve call using a given soltion vector
         *
         * @param[in] aFullVector     Solution Vector
         */
        void solve( sol::Dist_Vector * aFullVector );

        //-------------------------------------------------------------------------------
        /**
         * @brief Solve call
         *
         * @param[in] aFullVector
         */
        void solve();

        //-------------------------------------------------------------------------------

        void set_lambda_increment( moris::real aLambdaInc );

        //-------------------------------------------------------------------------------

        moris::real get_new_lambda();
    };
}
}
#endif /* MORIS_DISTLINALG_CL_TSA_MONOLITHIC_TIME_SOLVER_HPP_ */
