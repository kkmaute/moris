/*
 * cl_Linear_Solver_Aztec.hpp
 *
 *  Created on: May 14, 2018
 *      Author: schmidt
 */
#ifndef SRC_DISTLINALG_CL_LINEAR_SOLVER_AZTEC_HPP_
#define SRC_DISTLINALG_CL_LINEAR_SOLVER_AZTEC_HPP_

#include "cl_Linear_Solver_Trilinos.hpp"

namespace moris
{
class Linear_Solver_Aztec : public Linear_Solver_Trilinos
{
private:

    AztecOO                    mAztecSolver;

    Teuchos::ParameterList                mlParams;
    ML_Epetra::MultiLevelPreconditioner * mMlPrec;

protected:
public:
    Linear_Solver_Aztec( Solver_Input*   aInput );

    ~Linear_Solver_Aztec();

    void set_solver_parameters();

    moris::sint solve_linear_system();

    void set_solver_internal_parameters();
};
}

#endif /* SRC_DISTLINALG_CL_LINEAR_SOLVER_AZTEC_HPP_ */
