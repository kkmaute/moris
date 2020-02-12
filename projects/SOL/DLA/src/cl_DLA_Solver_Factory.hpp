/*
 * cl_DLA_Solver_Factory.hpp
 *
 *  Created on: Mar 28, 2018
 *      Author: schmidt
 */
#ifndef SRC_DISTLINALG_CL_DLA_SOLVER_FACTORY_HPP_
#define SRC_DISTLINALG_CL_DLA_SOLVER_FACTORY_HPP_

#include <memory>

#include "cl_DLA_Linear_Solver_Algorithm.hpp"

#include "cl_DLA_Linear_Solver_Aztec.hpp"
#include "cl_DLA_Linear_Solver_Amesos.hpp"
#include "cl_DLA_Linear_Solver_Amesos2.hpp"
#include "cl_DLA_Linear_Solver_PETSc.hpp"

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

        std::shared_ptr< Linear_Solver_Algorithm > create_solver( const enum SolverType    aSolverType = SolverType::AZTEC_IMPL );

        Linear_Problem * create_linear_system(       moris::Solver_Interface * aSolverInterface,
                                               const enum MapType              aLinSysType = MapType::Epetra,
                                               const bool                      aNotCreatedByNonLinSolver = false);

        Linear_Problem * create_linear_system(       moris::Solver_Interface * aSolverInterface,
        		Dist_Map               * aMap,
				Dist_Map               * aFullMap,
                                               const enum MapType              aLinSysType = MapType::Epetra,
                                               const bool                      aNotCreatedByNonLinSolver = false);
    };
    }
}

#endif /* SRC_DISTLINALG_CL_DLA_SOLVER_FACTORY_HPP_ */
