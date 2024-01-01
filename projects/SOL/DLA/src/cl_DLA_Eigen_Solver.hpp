/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
The Eigen Solver class.

*/
/* ------------------------------------------------------------------------- */

#ifndef SOLVERS_EIGEN_SOLVER_H_
#define SOLVERS_EIGEN_SOLVER_H_

// C system files

// C++ system files
#include <cstddef>
#include <assert.h>
#include "moris_typedefs.hpp"

// Project header files
#include "fn_PRM_SOL_Parameters.hpp"
#include "cl_DLA_Linear_Solver_Algorithm_Trilinos.hpp"
#include "cl_DLA_Preconditioner_Trilinos.hpp"
#include "cl_SOL_Amesos_GenOp.hpp"
#include "cl_Sparse_Matrix_EpetraFECrs.hpp"
#include "cl_Vector_Epetra.hpp"

// TPL header files

// Project header files
#include "AnasaziConfigDefs.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziBlockDavidsonSolMgr.hpp"
#include "AnasaziBasicOutputManager.hpp"
#include "AnasaziEpetraAdapter.hpp"

#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MultiVectorOut.h"

#include "Epetra_FEVector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_MultiVector.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_Vector.h"

#include "Amesos.h"
#include "Ifpack.h"
#include "Ifpack_Preconditioner.h"
#include "Epetra_InvOperator.h"

// ML
#include "ml_include.h"
#include "ml_epetra_utils.h"
#include "ml_epetra_preconditioner.h"

// Class forward declarations
class Vector_Epetra;
class SparseMatrix;
class Epetra_Vector;

/* ------------------------------------------------------------------------- */
namespace moris
{
    namespace dla
    {
        class Linear_Problem;
        class Eigen_Solver : public Linear_Solver_Algorithm_Trilinos
        {
          private:
            Preconditioner_Trilinos* mLeftPreconditioner = nullptr;

          protected:
            typedef Epetra_MultiVector                                    MV;
            typedef Epetra_Operator                                       OP;
            typedef Anasazi::MultiVecTraits< double, Epetra_MultiVector > MVT;
            typedef Anasazi::OperatorTraits< double, MV, OP >             OPT;

            Epetra_FECrsMatrix* mMat;
            Epetra_FECrsMatrix* mMassMat;

            sol::Dist_Matrix* mNewMat;

            real mSolTime;

            Anasazi::BasicOutputManager< double > mPrinter;
            Teuchos::RCP< Epetra_MultiVector >    mIvec;

            Teuchos::RCP< Epetra_Operator > mSPmat;        // K matrix sparse
            Teuchos::RCP< Epetra_Operator > mSPmassmat;    // M matrix sparse

            sol::Dist_Map* mMap;

            Vector_Epetra* mFreeSolVec = nullptr;    // Solution vector containing a certain eigenvector (ony free DoFs)

            std::vector< Anasazi::Value< double > > mevals;

            Teuchos::RCP< Anasazi::BasicEigenproblem< double, MV, OP > > mMyEigProblem;

            Anasazi::Eigensolution< double, MV > mSol;

            int mNumReturnedEigVals;

            std::vector< double > mNormEvecsMassMatEvecs;    // Vector if phi^T * M * phi norms. One for each eigenpair

            sol::EigSolMethod mEigSolMethod;

            Epetra_LinearProblem mEpetraProblem;

            bool mComputePrecondionerOpertor = false;

          public:
            //-----------------------------------------------------------------------
            /**
             * @brief Constructor. Creates a default eigen solver.
             */
            Eigen_Solver();

            Eigen_Solver( const ParameterList* aParameterList );

            //-----------------------------------------------------------------------
            /**
             * @brief Destructor.
             */

            // virtual ~Eigen_Solver(){};
            ~Eigen_Solver();

            //-----------------------------------------------------------------------
            /**
             * @brief extracts mass and stiffness matrix from linear system
             */

            void build_linearized_system( Linear_Problem* aLinearSystem );

            //-----------------------------------------------------------------------
            /**
             * @brief solve eigen system through eigen algorithm
             */

            moris::sint
            solve_linear_system()
            {
                return 0;
            };

            //-----------------------------------------------------------------------

            /**
             * @brief Solve linear system
             *
             * @param[in] aLinearSystem Pointer to linear system.
             * @param[in] aLinearSystem Iteration number.
             */
            moris::sint solve_linear_system( Linear_Problem* aLinearSystem,
                    const moris::sint                        aIter );

            //-----------------------------------------------------------------------
            /**
             * @brief sets parameters for eigen solver
             */

            void set_eigen_solver_manager_parameters();

            //-----------------------------------------------------------------------

            /**
             * @brief Eigenvalue solver using Anasazi's implementation of the Block Davidson method
             * @param[in] aLinearSystem Pointer to linear system.
             */
            int solve_block_davidson_system( Linear_Problem* aLinearSystem );

            //-----------------------------------------------------------------------

            /**
             * @brief Eigenvalue solver using Anasazi's implementation of the Generalized Davidson method an Ifpack ILUT factorization of K used as a preconditioner
             * @param[in] aLinearSystem Pointer to linear system.
             */
            int solve_generalized_davidson_system( Linear_Problem* aLinearSystem );

            //-----------------------------------------------------------------------

            /**
             * @brief Eigenvalue solver using Anasazi's implementation of the Block Krylov Schur method with a GMRES and a preconditioner (Aztec)
             * @param[in] aLinearSystem Pointer to linear system.
             */
            int solve_block_krylov_schur_system( Linear_Problem* aLinearSystem );

            //-----------------------------------------------------------------------

            /**
             * @brief Eigenvalue solver using Anasazi's implementation of the Block Krylov Schur Amesos method
             * @param[in] aLinearSystem Pointer to linear system.
             */
            int solve_block_krylov_schur_amesos_system( Linear_Problem* aLinearSystem );

            //-----------------------------------------------------------------------
            /**
             * @brief get eigenvector solution and save as hdf5 file
             * @param[in] aEigValIndex  Index of eigenvalue
             * @param[in] aSolVec       Current eigenvector
             * @param[in] aLeaderSolVec Leader solution vector
             * @param[in] aEigValReal   Real eigenvalue
             * @param[in] aEigValImag   Imaginary eigenvalue
             */

            int get_solution(
                    uint           aEigValIndex,
                    Vector_Epetra* aSolVec,
                    Vector_Epetra* aLeaderSolVec,
                    real&          aEigValReal,
                    real&          aEigValImag );

            //-----------------------------------------------------------------------
            /**
             * @brief returns eigenvector as distributed vector
             */

            sol::Dist_Vector*
            get_eigen_vector()
            {
                return mFreeSolVec;
            }

            //-------------------------------------------------------------------------------------------------
            /**
             * @brief return eigenvalues
             */

            std::vector< Anasazi::Value< double > >
            get_eigen_values()
            {
                return mevals;
            }

            //-------------------------------------------------------------------------------------------------

            /**
             * set the left hand side preocndioer, is useful for eigenanlaysis of the P^{-1}A system instead of the A system
             *
             */

            virtual void
            set_left_hand_side_preconditioner( Preconditioner_Trilinos* aPreconditioner )
            {
                mLeftPreconditioner = aPreconditioner;
            }

            //-------------------------------------------------------------------------------------------------

            /**
             * @brief Compute the preconditioned operator
             * 
             */

            void 
            compute_preconditioned_operator();
        };

    }    // namespace dla

}    // namespace moris
#endif
