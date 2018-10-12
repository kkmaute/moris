/*
 * cl_Linear_Solver_Amesos.hpp
 *
 *  Created on: May 16, 2018
 *      Author: schmidt
 */
#ifndef SRC_DISTLINALG_CL_LINEAR_SOLVER_AMESOS_HPP_
#define SRC_DISTLINALG_CL_LINEAR_SOLVER_AMESOS_HPP_

// TPL header files
#include "Epetra_ConfigDefs.h"
#include "Amesos_ConfigDefs.h"

#include "cl_Linear_Solver_Trilinos.hpp"

#include "Amesos.h"
#include "Amesos_BaseSolver.h"

namespace moris
{
class Linear_Solver_Amesos : public Linear_Solver_Trilinos
{
private:

    Amesos_BaseSolver *mAmesosSolver;
    Amesos            mAmesosFactory;

    bool              mIsPastFirstSolve;

protected:
public:
    Linear_Solver_Amesos(Solver_Interface * aInput);

    ~Linear_Solver_Amesos();

    void set_solver_parameters();

    //int SetSystemMatrix ( bool aUseTranspose );

    moris::sint solve_linear_system();

    void set_solver_internal_parameters();
};
}

#endif /* SRC_DISTLINALG_CL_LINEAR_SOLVER_AMESOS_HPP_ */
