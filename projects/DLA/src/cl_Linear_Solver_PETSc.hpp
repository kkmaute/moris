/*
 * cl_Linear_Solver_PETSc.hpp
 *
 *  Created on: Mar 25, 2018
 *      Author: schmidt
 */

#ifndef SRC_DISTLINALG_CL_LINEAR_SOLVER_PETSC_HPP_
#define SRC_DISTLINALG_CL_LINEAR_SOLVER_PETSC_HPP_

#include "core.hpp"
#include "cl_Linear_Solver.hpp"
#include "cl_VectorPETSc.hpp"
#include "cl_MatrixPETSc.hpp"

#include "cl_Matrix_Vector_Factory.hpp"
#include "cl_Solver_Input.hpp"

#include "cl_Model_Solver_Interface_Solver.hpp"

namespace moris
{
class Dist_Vector;
class Sparse_Matrix;
class Linear_Solver_PETSc : public moris::Linear_Solver
{
private:
    KSP                  mksp;
    PC                   mpc;

protected:

    moris::Param_List< boost::variant< bool, moris::sint, moris::real, const char* > > mParameterList; // The Algorithm specific parameter list

public:

//    Linear_Solver_PETSc(Mat          aPETScMat,
//                        Vec          aPETScVector_x,
//                        Vec          aPETScVector_b);

    Linear_Solver_PETSc( moris::Solver_Input * aInput );

    ~Linear_Solver_PETSc();

//    void build_linear_system( Epetra_FECrsMatrix*       aEpetraMat,
//                              Epetra_FEVector*          aEpetraVector_x,
//                              Epetra_FEVector*          aEpetraVector_b ){};

    void assemble_residual_and_jacobian( Dist_Vector * aFullSolutionVector ){MORIS_ERROR( false, "not implemented in Petsc yet");};

    void assemble_residual_and_jacobian(){MORIS_ERROR( false, "not implemented in Petsc yet");};

    void build_linear_system();

    moris::sint solve_linear_system();

    void solve_eigenvalues(){};

    void get_solution( moris::Matrix< DDRMat > & LHSValues );

    //Vector_Epetra* GetVec()       { return mEpetraVector; };

    //void solve_linear_system();

    boost::variant< bool, moris::sint, moris::real, const char* > & set_param( char const* aKey )
    {
        return mParameterList(aKey);
    }

};
}


#endif /* SRC_DISTLINALG_CL_LINEAR_SOLVER_PETSC_HPP_ */
