/*
 * cl_DLA_Linear_Solver_PETSc.hpp
 *
 *  Created on: Dez 11, 2018
 *      Author: schmidt
 */
#ifndef SRC_DISTLINALG_CL_LINEAR_SOLVER_PETSC_HPP_
#define SRC_DISTLINALG_CL_LINEAR_SOLVER_PETSC_HPP_

#include "core.hpp"
#include "cl_DLA_Linear_Solver.hpp"
#include "cl_VectorPETSc.hpp"
#include "cl_MatrixPETSc.hpp"

#include "cl_Matrix_Vector_Factory.hpp"
#include "cl_DLA_Solver_Interface.hpp"

#include "cl_DLA_Linear_Problem.hpp"

namespace moris
{
class Dist_Vector;
class Sparse_Matrix;
namespace dla
{
class Linear_Solver_PETSc : public moris::dla::Linear_Solver
{
    private:
        Linear_Problem * mLinearSystem;

        KSP mPetscKSPProblem;

        PC mpc;

    protected:

    public:
    Linear_Solver_PETSc();

    Linear_Solver_PETSc( moris::Solver_Interface * aInput );

    Linear_Solver_PETSc( Linear_Problem * aLinearSystem );

    ~Linear_Solver_PETSc();

    void set_linear_problem( Linear_Problem * aLinearSystem );

    void set_solver_parameters();

    void set_solver_internal_parameters();

    moris::sint solve_linear_system();

    moris::sint solve_linear_system(       Linear_Problem * aLinearSystem,
                                     const moris::sint      aIter );

    void build_multigrid_preconditioner( Linear_Problem * aLinearSystem );

//    void solve_eigenvalues(){};
//
//    void get_solution( moris::Matrix< DDRMat > & LHSValues );

};
}
}

#endif /* SRC_DISTLINALG_CL_LINEAR_SOLVER_PETSC_HPP_ */
