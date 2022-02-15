/*
 * cl_NLA_Convergence.hpp
 *
 *  Created on: Sep 21, 2018
 *      Author: schmidt
 */
#ifndef SRC_FEM_CL_NLA_CONVERGENCE_HPP_
#define SRC_FEM_CL_NLA_CONVERGENCE_HPP_

#include "typedefs.hpp"
#include "cl_NLA_Nonlinear_Algorithm.hpp"

namespace moris
{
    namespace sol
    {
        class Dist_Vector;
    }
    namespace NLA
    {
        class Convergence
        {
          private:

          public:
            Convergence(){};

            ~Convergence(){};

            bool check_for_convergence(
                    Nonlinear_Algorithm* tNonLinSolver,
                    moris::sint&         aIt,
                    moris::real&         aRefNorm,
                    moris::real&         aResNorm,
                    const moris::real&   aAssemblyTime,
                    const moris::real&   aSolvTime,
                    bool&                aHartBreak );

            bool check_for_convergence(
                    Nonlinear_Algorithm* tNonLinSolver,
                    moris::sint&         aIt,
                    moris::real&         aRefNorm,
                    moris::real&         aResNorm,
                    const moris::real&   aSolvTime,
                    bool&                aHartBreak );
        };
    }    // namespace NLA
}    // namespace moris

#endif /* SRC_FEM_CL_NLA_CONVERGENCE_HPP_ */
