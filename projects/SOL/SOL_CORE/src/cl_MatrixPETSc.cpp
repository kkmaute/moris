/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MatrixPETSc.cpp
 *
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

// ----------------------------------------------------------------------------

Matrix_PETSc::Matrix_PETSc(
        moris::Solver_Interface* aInput,
        sol::Dist_Map*           aMap )
        : sol::Dist_Matrix( aMap )
{
    moris::uint aNumMyDofs = aInput->get_my_local_global_map().numel();

    moris::Matrix< DDUMat > aMyConstraintDofs = aInput->get_constrained_Ids();

    // Fixme Implement nonzero algorithm
    PetscInt    tNonzeros  = 16;
    moris::uint tNumMyDofs = aNumMyDofs;

    // sum up all distributed dofs
    moris::uint tNumGlobalDofs = sum_all( tNumMyDofs );

    // FIXME insert boolean array for BC-- insert NumGlobalElements-- size
    mDirichletBCVec.set_size( tNumGlobalDofs, 1, 0 );

    // build BC vector
    this->dirichlet_BC_vector( mDirichletBCVec, aMyConstraintDofs );

    // Create and set Matrix
    MatCreate( PETSC_COMM_WORLD, &mPETScMat );

    MatSetSizes( mPETScMat, tNumMyDofs, tNumMyDofs, PETSC_DETERMINE, PETSC_DETERMINE );
    MatSetFromOptions( mPETScMat );
    MatMPIAIJSetPreallocation( mPETScMat, tNonzeros, NULL, tNonzeros, NULL );

    // FIXME extra matrix for serial (performance)
    // MatSeqAIJSetPreallocation(mPETScMat, tNonzeros, NULL);

    MatSetUp( mPETScMat );

    // allow for column based inputs
    MatSetOption( mPETScMat, MAT_ROW_ORIENTED, PETSC_FALSE );
}

// ----------------------------------------------------------------------------

Matrix_PETSc::Matrix_PETSc(
        const moris::uint aRows,
        const moris::uint aCols )
{
    // FIXME Implement nonzero algorithm
    PetscInt tNonzeros = 16;

    // Create and set Matrix
    MatCreate( PETSC_COMM_WORLD, &mPETScMat );

    mDirichletBCVec.set_size( aRows, 1, 0 );

    MatSetSizes( mPETScMat, aCols, aRows, PETSC_DETERMINE, PETSC_DETERMINE );
    MatSetFromOptions( mPETScMat );
    MatMPIAIJSetPreallocation( mPETScMat, tNonzeros, NULL, tNonzeros, NULL );

    //    MatSetOption( mPETScMat, MAT_COLUMN_ORIENTED, PETSC_TRUE );

    // FIXME extra matrix for serial (performance)
    // MatSeqAIJSetPreallocation(mPETScMat, tNonzeros, NULL);

    MatSetUp( mPETScMat );
    // allow for column based inputs
    MatSetOption( mPETScMat, MAT_ROW_ORIENTED, PETSC_FALSE );
}

// ----------------------------------------------------------------------------

Matrix_PETSc::~Matrix_PETSc()
{
    MatDestroy( &mPETScMat );
}

// ----------------------------------------------------------------------------

void
Matrix_PETSc::build_graph(
        const moris::uint&             aNumMyDof,
        const moris::Matrix< DDSMat >& aElementTopology )
{
    moris::Matrix< DDSMat > tTempElemDofs = aElementTopology;

    // Build Zero matrix and matrix for element free dof id
    moris::Matrix< DDRMat > tZeros( aNumMyDof * aNumMyDof, 1, 0.0 );

    // loop over elemental dofs
    for ( moris::uint Ij = 0; Ij < aNumMyDof; Ij++ )
    {
        // set constrDof to neg value
        if ( aElementTopology( Ij ) < 0 )
        {
            tTempElemDofs( Ij, 0 ) = -1;
        }
        else if ( aElementTopology( Ij ) > (sint)( mDirichletBCVec.length() - 1 ) )
        {
            tTempElemDofs( Ij, 0 ) = -1;
        }
        // set constrDof to neg value
        if ( mDirichletBCVec( aElementTopology( Ij ), 0 ) == 1 )
        {
            tTempElemDofs( Ij, 0 ) = -1;
        }
    }

    // Applying Petsc map AO
    AOApplicationToPetsc( mMap->get_petsc_map(), aNumMyDof, tTempElemDofs.data() );

    MatSetValues( mPETScMat, aNumMyDof, tTempElemDofs.data(), aNumMyDof, tTempElemDofs.data(), tZeros.data(), ADD_VALUES );
    //    MatSetValues( mPETScMat, aNumMyDof, aElementTopology.data(), aNumMyDof, aElementTopology.data(), tZeros.data(), ADD_VALUES );
    // MatSetValuesBlocked();                                                  //important+
}

// ----------------------------------------------------------------------------

void
Matrix_PETSc::fill_matrix(
        const moris::uint&             aNumMyDof,
        const moris::Matrix< DDRMat >& aA_val,
        const moris::Matrix< DDSMat >& aEleDofConectivity )
{
    moris::Matrix< DDSMat > tTempElemDofs( aNumMyDof, 1 );
    tTempElemDofs = aEleDofConectivity;

    // loop over elemental dofs
    for ( moris::uint Ij = 0; Ij < aNumMyDof; Ij++ )
    {
        // set constrDof to neg value
        if ( mDirichletBCVec( aEleDofConectivity( Ij ), 0 ) == 1 )
        {
            tTempElemDofs( Ij, 0 ) = -1;
        }
    }

    // Applying Petsc map AO
    AOApplicationToPetsc( mMap->get_petsc_map(), aNumMyDof, tTempElemDofs.data() );

    MatSetValues( mPETScMat, aNumMyDof, tTempElemDofs.data(), aNumMyDof, tTempElemDofs.data(), aA_val.data(), ADD_VALUES );
    // MatSetValuesBlocked();                                                  //important+
}

// ----------------------------------------------------------------------------

void
Matrix_PETSc::insert_values(
        const Matrix< DDSMat >& aRowIDs,
        const Matrix< DDSMat >& aColumnIDs,
        const Matrix< DDRMat >& aMatrixValues )
{
    MatSetValues( mPETScMat,
            aRowIDs.numel(),
            aRowIDs.data(),
            aColumnIDs.numel(),
            aColumnIDs.data(),
            aMatrixValues.data(),
            INSERT_VALUES );
}

// ----------------------------------------------------------------------------

void
Matrix_PETSc::sum_into_values(
        const Matrix< DDSMat >& aRowIDs,
        const Matrix< DDSMat >& aColumnIDs,
        const Matrix< DDRMat >& aMatrixValues )
{
    MatSetValues( mPETScMat,
            aRowIDs.numel(),
            aRowIDs.data(),
            aColumnIDs.numel(),
            aColumnIDs.data(),
            aMatrixValues.data(),
            ADD_VALUES );
}

// ----------------------------------------------------------------------------

void
Matrix_PETSc::get_matrix_values(
        const moris::Matrix< DDSMat >& aRequestedIds,
        moris::Matrix< DDRMat >&       aValues )
{
    // get values in row based format. There is no other way
    MatGetValues( mPETScMat,
            aRequestedIds.numel(),
            aRequestedIds.data(),
            aRequestedIds.numel(),
            aRequestedIds.data(),
            aValues.data() );

    // moris is column based.
    aValues = trans( aValues );
}

// ----------------------------------------------------------------------------

void
Matrix_PETSc::matrix_global_assembly()
{
    MatAssemblyBegin( mPETScMat, MAT_FINAL_ASSEMBLY );
    MatAssemblyEnd( mPETScMat, MAT_FINAL_ASSEMBLY );

    // MatView(mPETScMat, PETSC_VIEWER_STDOUT_(PETSC_COMM_WORLD) );
}

// ----------------------------------------------------------------------------

void
Matrix_PETSc::dirichlet_BC_vector(
        moris::Matrix< DDUMat >&       aDirichletBCVec,
        const moris::Matrix< DDUMat >& aMyConstraintDofs )
{
    // build vector with constraint values. unconstraint=0 constraint =1. change this to true/false
    for ( moris::uint Ik = 0; Ik < aMyConstraintDofs.n_rows(); Ik++ )
    {
        aDirichletBCVec( aMyConstraintDofs( Ik, 0 ), 0 ) = 1;
    }
}

// ----------------------------------------------------------------------------

void
Matrix_PETSc::print() const
{
    MatView( mPETScMat, PETSC_VIEWER_STDOUT_( PETSC_COMM_WORLD ) );
}

// ----------------------------------------------------------------------------

void
Matrix_PETSc::save_matrix_to_matlab_file( const char* aFilename )
{
    PetscViewer tViewer;

    PetscViewerCreate( PETSC_COMM_WORLD, &tViewer );
    PetscViewerSetType( tViewer, PETSCVIEWERASCII );
    PetscViewerPushFormat( tViewer, PETSC_VIEWER_ASCII_MATLAB );
    PetscViewerFileSetName( tViewer, aFilename );

    MatView( mPETScMat, tViewer );

    PetscViewerDestroy( &tViewer );
}
