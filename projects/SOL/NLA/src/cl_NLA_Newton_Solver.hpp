/*
 * cl_NLA_Newton_Solver.hpp
 *
 *  Created on: Sep 21, 2018
 *      Author: schmidt
 */
#ifndef SRC_FEM_CL_NEWTON_SOLVER_HPP_
#define SRC_FEM_CL_NEWTON_SOLVER_HPP_

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
    class Newton_Solver : public Nonlinear_Algorithm
    {
    private:
        /**
         * @brief Call for solve of linear system
         *
         * @param[in] aIter       Number of newton iterations
         * @param[in] aHardBreak  Flag for HartBreak
         */
        void solve_linear_system( moris::sint & aIter,
                                  bool        & aHardBreak);

    public:
        /**
         * @brief Constructor for Newton
         *
         */
        Newton_Solver();

        Newton_Solver( dla::Linear_Solver * aLinSolver );

        ~Newton_Solver();

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
         * @brief Accessor to set a value in the parameter list of the Newton solver
         *
         * @param[in] aKey Key corresponding to the mapped value that
         *            needs to be accessed
         */
        boost::variant< bool, sint, real, const char* > &  set_param( char const* aKey )
        {
            return mParameterListNonlinearSolver( aKey );
        }

    };
}
}

#endif /* SRC_FEM_CL_NEWTON_SOLVER_HPP_ */
