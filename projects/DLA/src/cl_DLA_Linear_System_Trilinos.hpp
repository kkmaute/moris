/*
 * cl_DLA_Linear_System_Trilinos.hpp
 *
 *  Created on: Dec 6, 2017
 *      Author: schmidt
 */
#ifndef SRC_DISTLINALG_CL_DLA_LINEAR_SYSTEM_TRILISOS_HPP_
#define SRC_DISTLINALG_CL_DLA_LINEAR_SYSTEM_TRILISOS_HPP_

#include <map>
#include <vector>

#include "AztecOO.h"

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
#include "cl_DLA_Linear_Problem.hpp"
#include "cl_Map_Class.hpp"

#include "cl_Param_List.hpp" // CON/src

namespace moris
{
namespace dla
{
    class Linear_System_Trilinos : public Linear_Problem
    {
    private:

    protected:

    public:
        Linear_System_Trilinos( Solver_Interface * aInput );

        Linear_System_Trilinos( const char* aString );

        ~Linear_System_Trilinos();

        void assemble_residual_and_jacobian( Dist_Vector * aFullSolutionVector );
        void assemble_residual( Dist_Vector * aFullSolutionVector );
        void assemble_jacobian( Dist_Vector * aFullSolutionVector );

        void assemble_residual_and_jacobian( );

        void build_linear_system();

        virtual moris::sint solve_linear_system();

        void get_solution( Matrix< DDRMat > & LHSValues );
    };
}
}
#endif /* SRC_DISTLINALG_CL_DLA_LINEAR_SYSTEM_TRILISOS_HPP_ */
