/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_DLA_Linear_System_PETSc.hpp
 *
 */

#ifndef SRC_DISTLINALG_CL_DLA_LINEAR_SYSTEM_PETSC_HPP_
#define SRC_DISTLINALG_CL_DLA_LINEAR_SYSTEM_PETSC_HPP_

#include <map>
#include <vector>
#include "moris_typedefs.hpp"

#include "AztecOO.h"

// #include "Ifpack.h"
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
#include "cl_SOL_Matrix_Vector_Factory.hpp"
#include "cl_DLA_Linear_Problem.hpp"
#include "cl_SOL_Dist_Map.hpp"

#include "cl_Parameter_List.hpp"    // CON/src

namespace moris
{
    namespace dla
    {
        class Linear_System_PETSc : public Linear_Problem
        {
          private:
            // Flag for deconstructor. If PetscFinalize should be called in linear solver or in nonlinear
            bool mNotCreatedByNonLinearSolver = false;

          protected:

          public:
            //--------------------------------------------------------------------------

            /**
             * @brief Construct a new Linear_System_PETSc object
             * 
             */
            Linear_System_PETSc(){};

            //--------------------------------------------------------------------------
            Linear_System_PETSc(
                    Solver_Interface* aInput,
                    const bool        aNotCreatedByNonLinSolver = false );

            Linear_System_PETSc(
                    Solver_Interface*   aInput,
                    sol::SOL_Warehouse* aSolverWarehouse,
                    sol::Dist_Map*      aFreeMap,
                    sol::Dist_Map*      aFullMap,
                    const bool          aNotCreatedByNonLinSolver = false );

            Linear_System_PETSc( const char* aString );

            ~Linear_System_PETSc();

            moris::sint solve_linear_system();

            void get_solution( Matrix< DDRMat >& LHSValues );
        };
    }    // namespace dla
}    // namespace moris
#endif /* SRC_DISTLINALG_CL_DLA_LINEAR_SYSTEM_PETSC_HPP_ */
