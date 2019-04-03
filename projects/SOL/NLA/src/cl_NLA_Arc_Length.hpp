/*
 * cl_NLA_Arc_Length.hpp
 *
 *  Created on: Mar 26, 2019
 *      Author: sonne
 */
#ifndef PROJECTS_SOL_NLA_SRC_CL_NLA_ARC_LENGTH_HPP_
#define PROJECTS_SOL_NLA_SRC_CL_NLA_ARC_LENGTH_HPP_

#include "typedefs.hpp"
#include "cl_NLA_Nonlinear_Algorithm.hpp"

namespace moris
{
class Dist_Vector;
namespace dla
{
    class Linear_Solver_Algorithm;
}
namespace NLA
{
    class Arc_Length_Solver : public Nonlinear_Algorithm
    {

    public:
        /**
         * @brief Constructor for Arc Length
         *
         */
        Arc_Length_Solver();

        Arc_Length_Solver( dla::Linear_Solver * aLinSolver );

        ~Arc_Length_Solver();

        /**
         * @brief Call to solve the nonlinear system
         *
         * @param[in] aNonlinearProblem Nonlinear problem
         */

        void solver_nonlinear_system( Nonlinear_Problem * aNonlinearProblem );

        void get_full_solution( moris::Matrix< DDRMat > & LHSValues );

        void get_solution( moris::Matrix< DDRMat > & LHSValues );

        void extract_my_values( const moris::uint             & aNumIndices,
                                const moris::Matrix< DDSMat > & aGlobalBlockRows,
                                const moris::uint             & aBlockRowOffsets,
                                      moris::Matrix< DDRMat > & LHSValues );

        /**
         * @brief Accessor to set a value in the parameter list of the Arc Length solver
         *
         * @param[in] aKey Key corresponding to the mapped value that
         *            needs to be accessed
         */
        boost::variant< bool, sint, real, const char* > &  set_param( char const* aKey )
        {
            return mParameterListNonlinearSolver( aKey );
        }

    private:
        /**
         * @brief Call for solve of linear system
         *
         * @param[in] aIter       Number of newton iterations
         * @param[in] aHardBreak  Flag for HartBreak
         */
        void solve_linear_system( moris::sint & aIter,
                                  bool        & aHardBreak);


    protected:

//        friend class Nonlinear_Problem;

    };
}
}



#endif /* PROJECTS_SOL_NLA_SRC_CL_NLA_ARC_LENGTH_HPP_ */
