/*
 * cl_Solver_Factory.hpp
 *
 *  Created on: Mar 28, 2018
 *      Author: schmidt
 */
#include <memory>
//#include "cl_Linear_Solver_Trilinos.hpp"
//#include "cl_Linear_Solver_PETSc.hpp"
#include "cl_Linear_Solver.hpp"

#ifndef SRC_DISTLINALG_CL_SOLVER_FACTORY_HPP_
#define SRC_DISTLINALG_CL_SOLVER_FACTORY_HPP_

namespace moris
{
    class Solver_Input;

    class Solver_Factory
    {
    private:

    protected:

    public:
        Solver_Factory();

        ~Solver_Factory();

        std::shared_ptr< Linear_Solver > create_solver( Solver_Input*  aInput );
    };
}

#endif /* SRC_DISTLINALG_CL_SOLVER_FACTORY_HPP_ */
