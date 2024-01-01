/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MatrixPETSc.hpp
 *
 */

#ifndef SRC_DISTLINALG_CL_MATRIXPETSC_HPP_
#define SRC_DISTLINALG_CL_MATRIXPETSC_HPP_

// C++ system files
#include <cstdio>
#include <iostream>

// MORIS header files.
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

// Project header files
#include "cl_Map_PETSc.hpp"
#include "cl_Vector_PETSc.hpp"
#include "cl_SOL_Dist_Matrix.hpp"

// TPL header files
#include <petsc.h>
#include <petscsys.h>

namespace moris
{
    class Matrix_PETSc : public sol::Dist_Matrix
    {
      private:
        moris::Matrix< DDUMat > mDirichletBCVec;

        void dirichlet_BC_vector(
                moris::Matrix< DDUMat >&       aDirichletBCVec,
                const moris::Matrix< DDUMat >& aMyConstraintDofs );

      protected:

      public:
        /** Default constructor */
        Matrix_PETSc(
                moris::Solver_Interface* aInput,
                sol::Dist_Map*           aMap );

        Matrix_PETSc( const moris::uint aRows,
                const moris::uint       aCols );

        /** Destructor */
        virtual ~Matrix_PETSc();

        void fill_matrix(
                const moris::uint&             aNumMyDofs,
                const moris::Matrix< DDRMat >& aA_val,
                const moris::Matrix< DDSMat >& aEleDofConnectivity );

        /**
         * Inserts values into the matrix at locations corresponding to the given row and column IDs.
         *
         * @param aRowIDs Row IDs
         * @param aColumnIDs Column IDs
         * @param aMatrixValues Values to be inserted
         */
        void insert_values(
                const Matrix< DDSMat >& aRowIDs,
                const Matrix< DDSMat >& aColumnIDs,
                const Matrix< DDRMat >& aMatrixValues );

        /**
         * Sums values into the matrix at locations corresponding to the given row and column IDs.
         *
         * @param aRowIDs Row IDs
         * @param aColumnIDs Column IDs
         * @param aMatrixValues Values to be summed into the existing matrix values
         */
        void sum_into_values(
                const Matrix< DDSMat >& aRowIDs,
                const Matrix< DDSMat >& aColumnIDs,
                const Matrix< DDRMat >& aMatrixValues );

        void get_matrix_values(
                const moris::Matrix< DDSMat >& aRequestedIds,
                moris::Matrix< DDRMat >&       aValues );

        void matrix_global_assembly();

        void build_graph(
                const moris::uint&             aNumMyDof,
                const moris::Matrix< DDSMat >& aElementTopology );

        void get_diagonal( moris::sol::Dist_Vector& aDiagVec ) const {};

        // FIXME mat_put_scalar only implemented for zeros with petsc. has to be changed
        void
        mat_put_scalar( const moris::real& aValue )
        {
            MORIS_ERROR( std::abs( aValue ) < MORIS_REAL_EPS,
                    "mat_put_scalar only implemented for zeros with petsc. " );

            MatZeroEntries( mPETScMat );
        }

        void
        sparse_mat_left_scale( const moris::sol::Dist_Vector& aScaleVector )
        {
            MORIS_ERROR( false, "not yet implemented for petsc" );
        };

        void
        sparse_mat_right_scale( const moris::sol::Dist_Vector& aScaleVector )
        {
            MORIS_ERROR( false, "not yet implemented for petsc" );
        };

        void
        replace_diagonal_values( const moris::sol::Dist_Vector& aDiagVec )
        {
            MORIS_ERROR( false, "not yet implemented for petsc" );
        };

        virtual void
        mat_vec_product(
                const moris::sol::Dist_Vector& aInputVec,
                moris::sol::Dist_Vector&       aResult,
                const bool                     aUseTranspose )
        {
            MORIS_ERROR( false, "not yet implemented for petsc" );
        };

        void print() const;

        void save_matrix_to_matlab_file( const char* aFilename );

        void save_matrix_to_matrix_market_file( const char* aFilename ){};

        void save_matrix_map_to_matrix_market_file( const char* aFilename ){};

        // ----------------------------------------------------------------------------

        /**
         * @brief build the graph of matrix, allocate enough space in CSR format
         *
         * @param aNumMyDof
         * @param aElementTopology
         */

        virtual void build_graph(
                Vector< moris_id >& aNonZeroDiagonal,
                Vector< moris_id >& aNonZeroOffDiagonal ) override;

        // void BuildSparseGraph(int numElements = 5);

        // Mat get_petsc_matrix()       { return mPETScMat; }
    };

}    // namespace moris

#endif /* SRC_DISTLINALG_CL_MATRIXPETSC_HPP_ */
