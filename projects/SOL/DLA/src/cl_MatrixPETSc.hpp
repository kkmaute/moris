/*
 * MatrixPETSc.hpp
 *
 *  Created on: Dec 5, 2018
 *      Author: schmidt
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
#include "cl_VectorPETSc.hpp"
#include "cl_Sparse_Matrix.hpp"

// TPL header files
#include <petsc.h>
#include <petscsys.h>

namespace moris
{
class Matrix_PETSc : public Sparse_Matrix
{
private:
    moris::Matrix< DDUMat >   mDirichletBCVec;

    void dirichlet_BC_vector(       moris::Matrix< DDUMat > & aDirichletBCVec,
                              const moris::Matrix< DDUMat > & aMyConstraintDofs);

protected:

public:
    /** Default contructor */
    Matrix_PETSc(       moris::Solver_Interface * aInput,
                  const moris::Map_Class        * aMap );

    Matrix_PETSc( const moris::uint aRows,
                  const moris::uint aCols );

    /** Destructor */
    ~Matrix_PETSc();

    void fill_matrix( const moris::uint             & aNumMyDofs,
                      const moris::Matrix< DDRMat > & aA_val,
                      const moris::Matrix< DDSMat > & aEleDofConectivity );

    void fill_matrix_row( const moris::Matrix< DDRMat > & aA_val,
                          const moris::Matrix< DDSMat > & aRow,
                          const moris::Matrix< DDSMat > & aCols );

    void matrix_global_assembly();

    void build_graph( const moris::uint             & aNumMyDof,
                      const moris::Matrix< DDSMat > & aElementTopology );

    void get_diagonal( moris::Dist_Vector & aDiagVec ) const{};

    //FIXME mat_put_scalar only implemented for zeros with petsc. has to be changed
    void mat_put_scalar( const moris::real & aValue )
    {
        MatZeroEntries( mPETScMat );
//      MORIS_ERROR(false, "mat_put_scalar only implemented for zeros with petsc. has to be changed.");
    }

    void sparse_mat_left_scale( const moris::Dist_Vector & aScaleVector ){};

    void sparse_mat_right_scale( const moris::Dist_Vector & aScaleVector ){};

    void replace_diagonal_values( const moris::Dist_Vector & aDiagVec ){};

    void print() const;

    void save_matrix_to_matlab_file( const char* aFilename ){};

    void save_matrix_to_matrix_market_file( const char* aFilename ){};

    void save_matrix_map_to_matrix_market_file( const char* aFilename ){};

    //void BuildSparseGraph(int numElements = 5);

    //Mat get_petsc_matrix()       { return mPETScMat; }
};

}

#endif /* SRC_DISTLINALG_CL_MATRIXPETSC_HPP_ */
