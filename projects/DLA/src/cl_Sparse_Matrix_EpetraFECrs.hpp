/*
 * SparseMatrixEpetra.hpp
 *
 *  Created on: Dec 6, 2017
 *      Author: schmidt
 */
#ifndef SRC_DISTLINALG_SPARSEMATRIXEPETRAFECRS_HPP_
#define SRC_DISTLINALG_SPARSEMATRIXEPETRAFECRS_HPP_

// MORIS header files.
#ifdef MORIS_HAVE_PARALLEL
 #include <mpi.h>
#endif

// C system files
#include <cstdio>
#include <iostream>

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_Map_Epetra.hpp"
#include "cl_Sparse_Matrix.hpp"
#include "cl_Vector_Epetra.hpp"
namespace moris
{

// Project header files
class Sparse_Matrix_EpetraFECrs : public Sparse_Matrix
{
private:
    //Epetra_FECrsMatrix *     mEpetraMat;
    //Epetra_FECrsGraph*       mEpetraGraph;
    //Map_Epetra *             mEpetraMap;
    //const Map_Class * mMap;
    //const Map_Epetra * mMap;

    moris::Matrix< DDUMat > DirichletBCVec;


    void dirichlet_BC_vector(       moris::Matrix< DDUMat > & aDirichletBCVec,
                              const moris::Matrix< DDUMat > & aMyConstraintDofs );

protected:

public:
    Sparse_Matrix_EpetraFECrs(       moris::Solver_Input * aInput,
                               const moris::Map_Class    * aMap );
    /** Destructor */
    ~Sparse_Matrix_EpetraFECrs();

    void fill_matrix( const moris::uint               & aNumMyDofs,
                      const moris::Matrix< DDRMat > & aA_val,
                      const moris::Matrix< DDSMat >         & aEleDofConectivity );

    void matrix_global_asembly();

    void build_graph( const moris::uint       & aNumMyDof,
                      const moris::Matrix< DDSMat > & aElementTopology );

    void get_diagonal( moris::Dist_Vector & aDiagVec ) const;

    void sparse_mat_left_scale( const moris::Dist_Vector & aScaleVector );

    void sparse_mat_right_scale( const moris::Dist_Vector & aScaleVector );

    void replace_diagonal_values( const moris::Dist_Vector & aDiagVec );

    void  print_matrix_to_screen() const;

    void save_matrix_to_matlab_file( char* aFilename );

    void save_matrix_to_matrix_market_file( const char* aFilename );

    void save_matrix_map_to_matrix_market_file( const char* aFilename );

    //const MapClass* GetMap() const { return ( Map_Class * ) mEpetraMap; }
    //MapClass* GetMap()       { return ( Map_Class * ) mEpetraMap; }
    //
    //const MapEpetra* GetFreeMap() const { return mEpetraMap; }
    //MapEpetra* get_epetra_free_map()       { return mEpetraMap; }
};
}

#endif /* SRC_DISTLINALG_SPARSEMATRIXEPETRAFECRS_HPP_ */
