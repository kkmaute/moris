/*
 * cl_DLA_Solver_Factory.hpp
 *
 *  Created on: Mar 28, 2018
 *      Author: schmidt
 */
#ifndef SRC_DISTLINALG_CL_DLA_SOLVER_FACTORY_HPP_
#define SRC_DISTLINALG_CL_DLA_SOLVER_FACTORY_HPP_

#include <memory>

#include "cl_DLA_Linear_Solver.hpp"

namespace moris
{
    class Solver_Interface;
    namespace dla
    {
    class Linear_Problem;

    class Solver_Factory
    {
    private:

    protected:

    public:
        Solver_Factory();

        ~Solver_Factory();

        std::shared_ptr< Linear_Solver > create_solver( const enum SolverType    aSolverType = SolverType::AZTEC_IMPL );

        std::shared_ptr< Linear_Problem > create_linear_system(      moris::Solver_Interface * aSolverInterface,
                                                               const enum MapType              aLinSysType = MapType::Epetra,
                                                               const bool                      aCreatedByNonLinSolver = false);
    };
    }
}

#endif /* SRC_DISTLINALG_CL_DLA_SOLVER_FACTORY_HPP_ */
