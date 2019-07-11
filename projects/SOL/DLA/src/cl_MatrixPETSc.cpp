/*
 * MatrixPETSc.cpp
 *
 *  Created on: Dec 5, 2018
 *      Author: schmidt
 */
#include "cl_MatrixPETSc.hpp"

#include <cstddef>
#include <cassert>
#include <algorithm>
#include <iostream>

#include "petscmat.h"

#include "fn_trans.hpp"

extern moris::Comm_Manager gMorisComm;

// TPL header files
using namespace moris;

Matrix_PETSc::Matrix_PETSc(       moris::Solver_Interface * aInput,
                            const moris::Map_Class        * aMap ) : Sparse_Matrix( aMap )
{
    //moris::uint             aNumMyDofs          = aInput->get_num_my_dofs();
    moris::uint aNumMyDofs = aInput->get_my_local_global_map().n_rows();
    //moris::Matrix< DDSMat > aMyLocaltoGlobalMap = aInput->get_my_local_global_map();
    moris::Matrix< DDUMat > aMyConstraintDofs   = aInput->get_constr_dof();

    // Fixme Implement nonzero algorithm
    PetscInt    tNonzeros =16;
    moris::uint    tNumMyDofs = aNumMyDofs;
    moris::uint tNumGlobalDofs=  aNumMyDofs;

    // sum up all distributed dofs
    sum_all( tNumMyDofs, tNumGlobalDofs );

    //FIXME insert boolian array for BC-- insert NumGlobalElements-- size
    mDirichletBCVec.set_size( tNumGlobalDofs, 1, 0 );

    // build BC vector
    this->dirichlet_BC_vector( mDirichletBCVec, aMyConstraintDofs );

    // Create and set Matrix
    MatCreate( PETSC_COMM_WORLD, &mPETScMat );

    MatSetSizes( mPETScMat, tNumMyDofs, tNumMyDofs, PETSC_DETERMINE, PETSC_DETERMINE );
    MatSetFromOptions( mPETScMat );
    MatMPIAIJSetPreallocation( mPETScMat, tNonzeros, NULL, tNonzeros, NULL );

    //FIXME extra matrix for serial (performance)
    //MatSeqAIJSetPreallocation(mPETScMat, tNonzeros, NULL);

    MatSetUp( mPETScMat );
}

Matrix_PETSc::Matrix_PETSc( const moris::uint aRows,
                            const moris::uint aCols )
{
    // Fixme Implement nonzero algorithm
    PetscInt    tNonzeros =16;

    // Create and set Matrix
    MatCreate( PETSC_COMM_WORLD, &mPETScMat );



    MatSetSizes( mPETScMat, aCols, aRows, PETSC_DETERMINE, PETSC_DETERMINE );
    MatSetFromOptions( mPETScMat );
    MatMPIAIJSetPreallocation( mPETScMat, tNonzeros, NULL, tNonzeros, NULL );

//    MatSetOption( mPETScMat, MAT_COLUMN_ORIENTED, PETSC_TRUE );

    //FIXME extra matrix for serial (performance)
    //MatSeqAIJSetPreallocation(mPETScMat, tNonzeros, NULL);

    MatSetUp( mPETScMat );
}

// ----------------------------------------------------------------------------
Matrix_PETSc::~Matrix_PETSc()
{
    MatDestroy( &mPETScMat );
}

void Matrix_PETSc::build_graph( const moris::uint             & aNumMyDof,
                                const moris::Matrix< DDSMat > & aElementTopology )
{
    moris::Matrix< DDSMat >tTempElemDofs( aNumMyDof, 1 );
    tTempElemDofs = aElementTopology;

    // Build Zero matrix and matrix for element free dof id
    moris::Matrix< DDRMat > tZeros (aNumMyDof*aNumMyDof, 1, 0.0);

    //loop over elemental dofs
    for ( moris::uint Ij=0; Ij< aNumMyDof; Ij++ )
    {
        //set constrDof to neg value
        if ( aElementTopology(Ij,0) < 0)
        {
            tTempElemDofs( Ij, 0) = -1;
        }
        else if ( aElementTopology(Ij,0) > (sint)(mDirichletBCVec.length()-1) )
        {
            tTempElemDofs( Ij, 0) = -1;
        }
        //set constrDof to neg value
        if ( mDirichletBCVec( aElementTopology(Ij,0),   0) == 1 )
        {
            tTempElemDofs( Ij, 0 ) = -1;
        }
     }



    // Applying Petsc map AO
    AOApplicationToPetsc( mMap->get_petsc_map(), aNumMyDof, tTempElemDofs.data() );

    MatSetValues( mPETScMat, aNumMyDof, tTempElemDofs.data(), aNumMyDof, tTempElemDofs.data(), tZeros.data(), ADD_VALUES );
    //MatSetValuesBlocked();                                                  //important+
}

void Matrix_PETSc::fill_matrix( const moris::uint             & aNumMyDof,
                                const moris::Matrix< DDRMat > & aA_val,
                                const moris::Matrix< DDSMat > & aEleDofConectivity)
{
    moris::Matrix< DDSMat >tTempElemDofs( aNumMyDof, 1 );
    moris::Matrix< DDRMat >tTempVal( aNumMyDof, aNumMyDof );
    tTempElemDofs = aEleDofConectivity;

    //loop over elemental dofs
    for ( moris::uint Ij=0; Ij< aNumMyDof; Ij++ )
    {
        //set constrDof to neg value
        if ( mDirichletBCVec( aEleDofConectivity(Ij,0),   0) == 1 )
        {
            tTempElemDofs( Ij, 0 ) = -1;
        }
     }

    tTempVal=trans(aA_val);

    // Applying Petsc map AO
    AOApplicationToPetsc( mMap->get_petsc_map(), aNumMyDof, tTempElemDofs.data() );

    MatSetValues( mPETScMat, aNumMyDof, tTempElemDofs.data(), aNumMyDof, tTempElemDofs.data(), tTempVal.data(), ADD_VALUES );
    //MatSetValuesBlocked();                                                  //important+
}

void Matrix_PETSc::fill_matrix_row( const moris::Matrix< DDRMat > & aA_val,
                                    const moris::Matrix< DDSMat > & aRow,
                                    const moris::Matrix< DDSMat > & aCols )
{
    MatSetValues( mPETScMat, 1, aRow.data(), aCols.length(), aCols.data(), aA_val.data(), INSERT_VALUES );
}

void Matrix_PETSc::matrix_global_assembly()
{
    MatAssemblyBegin( mPETScMat, MAT_FINAL_ASSEMBLY );
    MatAssemblyEnd( mPETScMat, MAT_FINAL_ASSEMBLY );

    //MatView(mPETScMat, PETSC_VIEWER_STDOUT_(PETSC_COMM_WORLD) );
}

void Matrix_PETSc::dirichlet_BC_vector(       moris::Matrix< DDUMat > & aDirichletBCVec,
                                        const moris::Matrix< DDUMat > & aMyConstraintDofs )
{
    //build vector with constraint values. unconstraint=0 constraint =1. change this to true/false
    for ( moris::uint Ik=0; Ik< aMyConstraintDofs.n_rows(); Ik++ )
    {
        aDirichletBCVec( aMyConstraintDofs( Ik, 0 )     ,0 )  = 1;
    }
}

void Matrix_PETSc::print() const
{
    MatView(mPETScMat, PETSC_VIEWER_STDOUT_(PETSC_COMM_WORLD) );
}


