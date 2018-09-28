/*
 * cl_NLA_Newton_Solver.hpp
 *
 *  Created on: Sep 21, 2018
 *      Author: schmidt
 */
#ifndef SRC_FEM_CL_NEWTON_SOLVER_HPP_
#define SRC_FEM_CL_NEWTON_SOLVER_HPP_

#include "typedefs.hpp"
#include "cl_NLA_Nonlinear_Solver.hpp"

namespace moris
{
class Dist_Vector;
class Linear_Solver;
namespace NLA
{
    class Newton_Solver : public Nonlinear_Solver
    {
    private:
        void solve_linear_system( moris::sint & aIter,
                                  bool        & aHardBreak);

        bool check_for_convergence(       moris::sint & aIt,
                                          moris::real & aRefNorm,
                                    const moris::real & aAssemblyTime,
                                          bool        & aHartBreak);

        void set_nonlinear_solver_parameters();

        moris::real get_time_needed( const clock_t aTime );

    public:
        Newton_Solver()
        {};

//        Newton_Solver( std::shared_ptr< Linear_Solver > aLinearSolver );

        ~Newton_Solver();

        void set_linear_solver( std::shared_ptr< Linear_Solver > aLinearSolver );

        void solver_nonlinear_system();

        Dist_Vector * get_full_sol_vec()
        {
            return mVectorFullSol;
        };

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
