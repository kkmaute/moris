/*
 * LinearSolverTrilinos.hpp
 *
 *  Created on: Dec 6, 2017
 *      Author: schmidt
 */
#ifndef SRC_DISTLINALG_LINEARSOLVERTRILINOS_HPP_
#define SRC_DISTLINALG_LINEARSOLVERTRILINOS_HPP_

#include <map>
#include <vector>

// TPL header files
//#include "AnasaziBlockKrylovSchurSolMgr.hpp"
//#include "AnasaziBasicEigenproblem.hpp"
//#include "AnasaziEpetraAdapter.hpp"
//#include "AnasaziLOBPCGSolMgr.hpp"

#include "AztecOO.h"
//#include "AztecOO_Operator.h"

// ML
#include "ml_include.h"
#include "ml_epetra_utils.h"
#include "ml_epetra_preconditioner.h"

//#include "Ifpack.h"
#include "Teuchos_SerialDenseMatrix.hpp"

#include "Epetra_Operator.h"
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MultiVectorOut.h"
#include "Epetra_ConfigDefs.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_FEVector.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"

#include "cl_Vector_Epetra.hpp"
#include "cl_Sparse_Matrix_EpetraFECrs.hpp"
#include "cl_Matrix_Vector_Factory.hpp"
#include "cl_Linear_Solver.hpp"
#include "cl_Model_Solver_Interface_Solver.hpp"
#include "cl_Map_Class.hpp"

#include "cl_Param_List.hpp" // CON/src

namespace moris
{
class Linear_Solver_Trilinos : public Linear_Solver
{
private:

protected:
    Epetra_LinearProblem      mEpetraProblem;
    moris::Matrix< DDRMat >   mSolution;
    moris::Matrix< DDUMat >   mAdofIndMap;      //FIXME added to map MSI ind to HMR ind. will be replaced

    Param_List< boost::variant< bool, sint, real, const char* > > mParameterList; // The Algorithm specific parameter list

    //Epetra_CrsMatrix*   mMatFromMatrixMatket;
    //Epetra_MultiVector* mVecRHSFromMatrixMatket;
    //Epetra_MultiVector* mVecLHSFromMatrixMatket;

public:

//    Linear_Solver_Trilinos(Epetra_FECrsMatrix*       aEpetraMat,
//                           Epetra_FEVector*          aEpetraVector_x,
//                           Epetra_FEVector*          aEpetraVector_b);

    Linear_Solver_Trilinos( Solver_Input* aInput );

    Linear_Solver_Trilinos( const char* aString );

    ~Linear_Solver_Trilinos();

    void assemble_residual_and_jacobian();

    void build_linear_system();

    virtual moris::sint solve_linear_system();

    void solve_eigenvalues();

    void get_solution(Matrix< DDRMat > & LHSValues);

    void get_solution_full( Matrix< DDRMat > & LHSValues );

    void extract_my_values( const moris::uint      & aNumIndices,
                            const Matrix< DDSMat > & aGlobalBlockRows,
                            const moris::uint      & aBlockRowOffsets,
                                  Matrix< DDRMat > & LHSValues );

    void import( );

    /**
     * @brief Accessor for the parameter list of the LinearSolver
     */
//    template< typename Variant = boost::variant< bool, sint, real, const char* > >
//    Param_List< Variant > &
//    params()
//    {
//        return mParameterList;
//    }

    /**
     * @brief Accessor to set a value in the parameter list of the LinearSolver
     *
     * @param[in] aKey Key corresponding to the mapped value that
     *            needs to be accessed
     */
    boost::variant< bool, sint, real, const char* > &  set_param( char const* aKey)
    {
        return mParameterList(aKey);
    }
};
}
#endif /* SRC_DISTLINALG_LINEARSOLVERTRILINOS_HPP_ */
