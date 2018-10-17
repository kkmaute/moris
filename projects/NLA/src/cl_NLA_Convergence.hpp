/*
 * cl_NLA_Convergence.hpp
 *
 *  Created on: Sep 21, 2018
 *      Author: schmidt
 */
#ifndef SRC_FEM_CL_NLA_CONVERGENCE_HPP_
#define SRC_FEM_CL_NLA_CONVERGENCE_HPP_

#include "typedefs.hpp"
#include "cl_NLA_Nonlinear_Solver.hpp"

namespace moris
{
class Dist_Vector;
namespace NLA
{
    class Convergence
    {
    private:

    public:
        Convergence()
        {};

        ~Convergence(){};

        bool check_for_convergence(       Nonlinear_Solver * tNonLinSolver,
                                          moris::sint & aIt,
                                          moris::real & aRefNorm,
                                    const moris::real & aAssemblyTime,
                                          bool        & aHartBreak);
    };
}
}

#endif /* SRC_FEM_CL_NLA_CONVERGENCE_HPP_ */
