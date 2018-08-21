/*
 * cl_Sparse_Matrix.hpp
 *
 *  Created on: Mar 27, 2018
 *      Author: schmidt
 */
#ifndef SRC_DISTLINALG_CL_SPARSE_MATRIX_HPP_
#define SRC_DISTLINALG_CL_SPARSE_MATRIX_HPP_

// MORIS header files.
#ifdef MORIS_HAVE_PARALLEL
 #include <mpi.h>
#endif

#include "cl_Vector_Epetra.hpp"
#include "cl_Solver_Input.hpp"
#include "linalg.hpp"

// TPL header files
#include "Epetra_FECrsMatrix.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_FECrsGraph.h"

#include "EpetraExt_OperatorOut.h"
#include "EpetraExt_CrsMatrixIn.h"
#include "Epetra_InvOperator.h"
#include "EpetraExt_BlockMapIn.h"
#include "EpetraExt_BlockMapOut.h"
#include "EpetraExt_RowMatrixOut.h"

#include <petsc.h>
#include <petscis.h>
#include <petscao.h>
#include <petscsys.h>

#include "cl_DistLinAlg_Enums.hpp"

class Sparse_Matrix
{
private:
protected:
          Epetra_FECrsMatrix   * mEpetraMat;
    const moris::Map_Class     * mMap;
          Mat                  mPETScMat;

public:
    Sparse_Matrix( const moris::Map_Class  * aMap ) : mEpetraMat( NULL ),
                                                      mMap( aMap ),
                                                      mPETScMat( NULL )
    {
    };

    virtual ~Sparse_Matrix(){};

    virtual void fill_matrix(const moris::uint               & anumDofs,
                             const moris::Mat< moris::real > & aA_val,
                             const moris::Mat< int >         & aEleDofConectivity) = 0;

    virtual void matrix_global_asembly() = 0;

    virtual void dirichlet_BC_vector(      moris::Mat< moris::uint > & aDirichletBCVec,
                                     const moris::Mat< uint >        & aMyConstraintDofs) = 0;

    virtual void build_graph(const moris::uint       & anumDofs,
                             const moris::Mat< int > & aEleDofConectivity) = 0;

    virtual void get_diagonal( moris::Dist_Vector & aDiagVec ) const = 0;

    virtual void sparse_mat_left_scale( const moris::Dist_Vector & aScaleVector ) = 0;

    virtual void sparse_mat_right_scale( const moris::Dist_Vector & aScaleVector ) = 0;

    virtual void replace_diagonal_values( const moris::Dist_Vector & aDiagVec ) = 0;

    virtual void print_matrix_to_screen() const = 0;

    virtual void save_matrix_to_matlab_file( char* aFilename ) = 0;

    virtual void save_matrix_to_matrix_market_file( const char* aFilename ) = 0;

    virtual void save_matrix_map_to_matrix_market_file( const char* aFilename ) = 0;

    //---------------------------------------------------------------------------------

    Epetra_FECrsMatrix* get_matrix()       { return mEpetraMat; }

    Mat get_petsc_matrix()       { return mPETScMat; }
};

#endif /* SRC_DISTLINALG_CL_SPARSE_MATRIX_HPP_ */
