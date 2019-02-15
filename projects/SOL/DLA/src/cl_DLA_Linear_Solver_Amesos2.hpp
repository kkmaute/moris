/*
 * cl_DLA_Linear_Solver_Amesos2.hpp
 *
 *  Created on: May 16, 2018
 *      Author: schmidt
 */
#ifndef SRC_DISTLINALG_CL_LINEAR_SOLVER_AMESOS2_HPP_
#define SRC_DISTLINALG_CL_LINEAR_SOLVER_AMESOS2_HPP_

#include "cl_DLA_Linear_Solver_Algorithm.hpp"

#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_CommandLineProcessor.hpp>

#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include "Amesos2.hpp"

namespace moris
{
namespace dla
{
//class Linear_Solver_Amesos2 : public Linear_Solver_Algorithm
//{
//private:
//    Teuchos::RCP< Amesos2::Solver< Epetra_CrsMatrix,Epetra_MultiVector > > mAmesos2Solver;
//
//    Epetra_LinearProblem      mEpetraProblem;
//
//    Teuchos::RCP<Epetra_CrsMatrix> mAMatrix;
//	Teuchos::RCP<Epetra_MultiVector> mRHS;
//    Teuchos::RCP<Epetra_MultiVector> mLHS;
//
//    bool              mIsPastFirstSolve;
//
//protected:
//public:
//    Linear_Solver_Amesos2(Solver_Interface * aInput);
//
//    ~Linear_Solver_Amesos2();
//
//    void set_solver_parameters();
//
//    //int SetSystemMatrix ( bool aUseTranspose );
//
//    moris::sint solve_linear_system();
//
//    void set_solver_internal_parameters();
//
//};
}
}



#endif /* SRC_DISTLINALG_CL_LINEAR_SOLVER_AMESOS2_HPP_ */
