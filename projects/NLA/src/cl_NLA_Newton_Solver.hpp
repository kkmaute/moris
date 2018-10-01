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
        /**
         * @brief Call for solve of linear system
         *
         * @param[in] aIter       Number of newton iterations
         * @param[in] aHardBreak  Flag for HartBreak
         */
        void solve_linear_system( moris::sint & aIter,
                                  bool        & aHardBreak);

        /**
         * @brief Convergence check. Will be pushed to extra class
         *
         */
        bool check_for_convergence(       moris::sint & aIt,
                                          moris::real & aRefNorm,
                                    const moris::real & aAssemblyTime,
                                          bool        & aHartBreak);

        /**
         * @brief Set the parameters in the nonlinear solver parameter list
         *
         */
        void set_nonlinear_solver_parameters();

        /**
         * @brief Member function which keeps track of used time for a particular purpose.
         *
         */
        moris::real get_time_needed( const clock_t aTime );

    public:
        /**
         * @brief Constructor for Newton
         *
         */
        Newton_Solver()
        {};

        ~Newton_Solver();

        /**
         * @brief Set function to set linear solver. All distributed vectores used by the newton are created here.
         *
         * @param[in] aLinearSolver   Pointer to linear solver
         *
         */
        void set_linear_solver( std::shared_ptr< Linear_Solver > aLinearSolver );

        /**
         * @brief Call for solve of nonlinear system
         *
         */
        void solver_nonlinear_system();

        Dist_Vector * get_full_sol_vec()
        {
            return mVectorFullSol;
        };

        void get_full_solution( moris::Matrix< DDRMat > & LHSValues );

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
