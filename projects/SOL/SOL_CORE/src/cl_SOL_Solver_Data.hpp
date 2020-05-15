/*
 * cl_NLA_Convergence.hpp
 *
 *  Created on: Sep 21, 2018
 *      Author: schmidt
 */
#ifndef SRC_FEM_CL_SOL_SOLVER_DATA_HPP_
#define SRC_FEM_CL_SOL_SOLVER_DATA_HPP_

#include "typedefs.hpp"

namespace moris
{
    namespace sol
    {
        class Dist_Vector;

        class Solver_Data
        {
            private:
                bool        mIsConverged = false;

                moris::uint mNumIterationsToConvergence = 0;

                moris::real mResidualNorm = MORIS_REAL_MAX;


            public:
                Solver_Data(){};

                ~Solver_Data(){};


        };
    }
}

#endif /* SRC_FEM_CL_SOL_SOLVER_DATA_HPP_ */
