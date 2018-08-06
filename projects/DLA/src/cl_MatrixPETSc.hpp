/*
 * MatrixPETSc.hpp
 *
 *  Created on: Jan 10, 2018
 *      Author: schmidt
 */
#ifndef SRC_DISTLINALG_CL_MATRIXPETSC_HPP_
#define SRC_DISTLINALG_CL_MATRIXPETSC_HPP_

// C system files

// C++ system files
#include <cstdio>
#include <iostream>

// MORIS header files.
#ifdef MORIS_HAVE_PARALLEL
 #include <mpi.h>
#endif

#include "linalg.hpp"



// Project header files
#include "cl_Map_PETSc.hpp"
#include "cl_VectorPETSc.hpp"
#include "cl_Sparse_Matrix.hpp"

/*
// TPL header files
//#include <petscksp.h>
#include <petsc.h>
//#include <petscis.h>
#include <petscsys.h>

class Matrix_PETSc : public Sparse_Matrix
{
private:
    //const moris::Map_Class*       mMap;
    //Mat                     mPETScMat;

    moris::Mat< moris::uint >   DirichletBCVec;

    void dirichlet_BC_vector(       moris::Mat< moris::uint > & aDirichletBCVec,
                              const moris::Mat< uint >        & aMyConstraintDofs);

protected:

public:

    Matrix_PETSc(       moris::Solver_Input * aInput,
                  const moris::Map_Class    * aMap );


    ~Matrix_PETSc();

    void fill_matrix( const moris::uint               & aNumMyDofs,
                      const moris::Mat< moris::real > & aA_val,
                      const moris::Mat< int >         & aEleDofConectivity );

    void matrix_global_asembly();

    void build_graph( const moris::uint       & aNumMyDof,
                      const moris::Mat< int > & aElementTopology ){};

    void get_diagonal( moris::Dist_Vector & aDiagVec ) const{};

    void sparse_mat_left_scale( const moris::Dist_Vector & aScaleVector ){};

    void sparse_mat_right_scale( const moris::Dist_Vector & aScaleVector ){};

    void replace_diagonal_values( const moris::Dist_Vector & aDiagVec ){};

    void  print_matrix_to_screen() const{};

    void save_matrix_to_matlab_file( char* aFilename ){};

    void save_matrix_to_matrix_market_file( const char* aFilename ){};

    void save_matrix_map_to_matrix_market_file( const char* aFilename ){};

    //void BuildSparseGraph(int numElements = 5);

    //Mat get_petsc_matrix()       { return mPETScMat; }
};*/

#endif /* SRC_DISTLINALG_CL_MATRIXPETSC_HPP_ */
