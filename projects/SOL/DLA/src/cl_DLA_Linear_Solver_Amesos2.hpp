/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_DLA_Linear_Solver_Amesos2.hpp
 *
 */

#ifndef SRC_DISTLINALG_CL_LINEAR_SOLVER_AMESOS2_HPP_
#define SRC_DISTLINALG_CL_LINEAR_SOLVER_AMESOS2_HPP_

#include "cl_DLA_Linear_Solver_Algorithm.hpp"

//#include "Amesos2.hpp"
//#include "Amesos2_Version.hpp"

//#include <Teuchos_ScalarTraits.hpp>
//#include <Teuchos_RCP.hpp>
//#include <Teuchos_GlobalMPISession.hpp>
//#include <Teuchos_VerboseObject.hpp>
//#include <Teuchos_CommandLineProcessor.hpp>
//
////#include <Tpetra_DefaultPlatform.hpp>
//#include <Tpetra_Map.hpp>
//#include <Tpetra_MultiVector.hpp>
//#include <Tpetra_CrsMatrix.hpp>
//
//#include "Amesos2.hpp"

namespace moris
{
namespace dla
{
class Linear_Solver_Amesos2 : public Linear_Solver_Algorithm
{
private:
//    Teuchos::RCP< Amesos2::Solver< Epetra_CrsMatrix,Epetra_MultiVector > > mAmesos2Solver;

//    Epetra_LinearProblem      mEpetraProblem;

    bool              mIsPastFirstSolve;

    Linear_Problem   * mLinearSystem =  nullptr;

protected:
public:
    Linear_Solver_Amesos2(Solver_Interface * aInput);

    Linear_Solver_Amesos2();

    ~Linear_Solver_Amesos2();

    void set_solver_parameters();

    //int SetSystemMatrix ( bool aUseTranspose );

    moris::sint solve_linear_system();

    moris::sint solve_linear_system(      Linear_Problem * aLinearSystem,
                                     const moris::sint     aIter );

    void set_solver_internal_parameters();

};
}
}

#endif /* SRC_DISTLINALG_CL_LINEAR_SOLVER_AMESOS2_HPP_ */

