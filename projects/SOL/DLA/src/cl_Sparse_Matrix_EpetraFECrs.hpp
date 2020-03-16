/*
 * SparseMatrixEpetra.hpp
 *
 *  Created on: Dec 6, 2017
 *      Author: schmidt
 */
#ifndef SRC_DISTLINALG_SPARSEMATRIXEPETRAFECRS_HPP_
#define SRC_DISTLINALG_SPARSEMATRIXEPETRAFECRS_HPP_

// MORIS header files.
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_Map_Epetra.hpp"
#include "cl_SOL_Dist_Matrix.hpp"
#include "cl_Vector_Epetra.hpp"

// C system files
#include <cstdio>


namespace moris
{
// Project header files
class Sparse_Matrix_EpetraFECrs : public Dist_Matrix
{
private:
    moris::Matrix< DDUMat > mDirichletBCVec;

    void dirichlet_BC_vector(       moris::Matrix< DDUMat > & aDirichletBCVec,
                              const moris::Matrix< DDUMat > & aMyConstraintDofs );

protected:

public:
    Sparse_Matrix_EpetraFECrs(       moris::Solver_Interface * aInput,
                               const moris::Dist_Map        * aMap );

    Sparse_Matrix_EpetraFECrs( const moris::uint aRows,
                               const moris::uint aCols )
    { MORIS_ERROR( false, "Sparse_Matrix_EpetraFECrs::Sparse_Matrix_EpetraFECrs: not set yet with epetra"); };

    /** Destructor */
    ~Sparse_Matrix_EpetraFECrs();

    void fill_matrix( const moris::uint             & aNumMyDofs,
                      const moris::Matrix< DDRMat > & aA_val,
                      const moris::Matrix< DDSMat > & aEleDofConectivity );

	void fill_matrix_row( const moris::Matrix< DDRMat > & aA_val,
                          const moris::Matrix< DDSMat > & aRow,
                          const moris::Matrix< DDSMat > & aCols )
	{ MORIS_ERROR( false, "Sparse_Matrix_EpetraFECrs::fill_matrix_row: not set yet with epetra"); };

    void matrix_global_assembly();

    void build_graph( const moris::uint             & aNumMyDof,
                      const moris::Matrix< DDSMat > & aElementTopology );

    void get_diagonal( moris::Dist_Vector & aDiagVec ) const;

    void mat_put_scalar( const moris::real & aValue );

    void sparse_mat_left_scale( const moris::Dist_Vector & aScaleVector );

    void sparse_mat_right_scale( const moris::Dist_Vector & aScaleVector );

    void replace_diagonal_values( const moris::Dist_Vector & aDiagVec );

    void  print() const;

    void save_matrix_to_matlab_file( const char* aFilename );

    void save_matrix_to_matrix_market_file( const char* aFilename );

    void save_matrix_map_to_matrix_market_file( const char* aFilename );

};
}

#endif /* SRC_DISTLINALG_SPARSEMATRIXEPETRAFECRS_HPP_ */
