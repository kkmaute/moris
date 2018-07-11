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

#include "linalg.hpp"

#include "cl_Map_Epetra.hpp"
#include "cl_Sparse_Matrix.hpp"
#include "cl_Vector_Epetra.hpp"

// Project header files
class Sparse_Matrix_EpetraFECrs : public Sparse_Matrix
{
private:
    //Epetra_FECrsMatrix *     mEpetraMat;
    //Epetra_FECrsGraph*       mEpetraGraph;
    //Map_Epetra *             mEpetraMap;
    //const Map_Class * mMap;
    //const Map_Epetra * mMap;

    moris::Mat< moris::uint > DirichletBCVec;


    void dirichlet_BC_vector(       moris::Mat< moris::uint > & aDirichletBCVec,
                              const moris::Mat< uint >        & aMyConstraintDofs );

protected:

public:
    Sparse_Matrix_EpetraFECrs(       moris::Solver_Input * aInput,
                               const moris::Map_Class    * aMap );
    /** Destructor */
    ~Sparse_Matrix_EpetraFECrs();

    void fill_matrix( const moris::uint               & aNumMyDofs,
                      const moris::Mat< moris::real > & aA_val,
                      const moris::Mat< int >         & aEleDofConectivity );

    void matrix_global_asembly();

    void build_graph( const moris::uint       & aNumMyDof,
                      const moris::Mat< int > & aElementTopology );

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

#endif /* SRC_DISTLINALG_SPARSEMATRIXEPETRAFECRS_HPP_ */
