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
    // get number of owned dofs
    moris::uint tMyNumOwnedDofs = aInput->get_my_local_global_map().numel();

    // get constrained dofs
    moris::Matrix< DDUMat > tMyConstraintDofs = aInput->get_constrained_Ids();

    // get total number of dofs
    moris::uint tNumGlobalDofs = sum_all( tMyNumOwnedDofs );

    // FIXME insert boolean array for BC or - better - create map
    mDirichletBCVec.set_size( tNumGlobalDofs, 1, 0 );

    // build BC vector
    this->dirichlet_BC_vector( mDirichletBCVec, tMyConstraintDofs );

    // Create matrix
    MatCreate( PETSC_COMM_WORLD, &mPETScMat );

    // Set size of matrix
    MatSetSizes( mPETScMat, tMyNumOwnedDofs, tMyNumOwnedDofs, PETSC_DETERMINE, PETSC_DETERMINE );

    // Set options
    MatSetFromOptions( mPETScMat );

       // Finalize setup of matrix
    MatSetUp( mPETScMat );
}

// ----------------------------------------------------------------------------

Matrix_PETSc::Matrix_PETSc(
        const moris::uint aRows,
        const moris::uint aCols )
{
    // build BC vector (here: no BCs)
    mDirichletBCVec.set_size( aRows, 1, 0 );

    // Create matrix
    MatCreate( PETSC_COMM_WORLD, &mPETScMat );

    // Set size of matrix
    MatSetSizes( mPETScMat, aCols, aRows, PETSC_DETERMINE, PETSC_DETERMINE );

    // Set options
    MatSetFromOptions( mPETScMat );

    // allow for column based inputs
    MatSetOption( mPETScMat, MAT_ROW_ORIENTED, PETSC_FALSE );

    // Fixme Implement sparsity algorithm
    PetscInt tNonZeros = 16;

    // Define sparsity structure
    MatMPIAIJSetPreallocation( mPETScMat, tNonZeros, NULL, tNonZeros, NULL );

    // Finalize setup of matrix
    MatSetUp( mPETScMat );
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

    // add values into matrix
    MatSetValues( mPETScMat, aNumMyDof, tTempElemDofs.data(), aNumMyDof, tTempElemDofs.data(), tZeros.data(), ADD_VALUES );
}

// ----------------------------------------------------------------------------

void
Matrix_PETSc::fill_matrix(
        const moris::uint&             aNumMyDof,
        const moris::Matrix< DDRMat >& aA_val,
        const moris::Matrix< DDSMat >& aEleDofConnectivity )
{
    // check for consistent sizes of vectors of IDs and values
    MORIS_ASSERT( aEleDofConnectivity.numel() == aNumMyDof,
            "Matrix_PETSc::fill_matrix - inconsistent sizes of ID and value vectors" );

    // create copy of vector with moris IDs; will be overwritten in AOApplicationToPetsc
    Matrix< DDSMat > tTempElemDofs = aEleDofConnectivity;

    // loop over elemental dofs
    for ( moris::uint Ij = 0; Ij < aNumMyDof; Ij++ )
    {
        // set constrDof to neg value
        if ( mDirichletBCVec( aEleDofConnectivity( Ij ), 0 ) == 1 )
        {
            tTempElemDofs( Ij, 0 ) = -1;
        }
    }

    // map moris IDs into petsc IDs
    AOApplicationToPetsc(
            mMap->get_petsc_map(),
            aNumMyDof,
            tTempElemDofs.data() );

    // add values into matrix
    MatSetValues(
            mPETScMat,
            aNumMyDof,
            tTempElemDofs.data(),
            aNumMyDof,
            tTempElemDofs.data(),
            aA_val.data(),
            ADD_VALUES );
}

// ----------------------------------------------------------------------------

void
Matrix_PETSc::insert_values(
        const Matrix< DDSMat >& aRowIDs,
        const Matrix< DDSMat >& aColumnIDs,
        const Matrix< DDRMat >& aMatrixValues )
{
    // create copies of moris IDs
    Matrix< DDSMat > tTempRowIDs    = aRowIDs;
    Matrix< DDSMat > tTempColumnIDs = aColumnIDs;

    // map moris IDs into petsc IDs if map exists
    if ( mMap )
    {
        AOApplicationToPetsc(
                mMap->get_petsc_map(),
                tTempRowIDs.numel(),
                tTempRowIDs.data() );

        AOApplicationToPetsc(
                mMap->get_petsc_map(),
                tTempColumnIDs.numel(),
                tTempColumnIDs.data() );
    }

    // insert values into matrix
    MatSetValues( mPETScMat,
            tTempRowIDs.numel(),
            tTempRowIDs.data(),
            tTempColumnIDs.numel(),
            tTempColumnIDs.data(),
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
    // create copies of moris IDs
    Matrix< DDSMat > tTempRowIDs    = aRowIDs;
    Matrix< DDSMat > tTempColumnIDs = aColumnIDs;

    // map moris IDs into petsc IDs if map exists
    if ( mMap )
    {
        AOApplicationToPetsc(
                mMap->get_petsc_map(),
                tTempRowIDs.numel(),
                tTempRowIDs.data() );

        AOApplicationToPetsc(
                mMap->get_petsc_map(),
                tTempColumnIDs.numel(),
                tTempColumnIDs.data() );
    }

    // insert values into matrix
    MatSetValues( mPETScMat,
            tTempRowIDs.numel(),
            tTempRowIDs.data(),
            tTempColumnIDs.numel(),
            tTempColumnIDs.data(),
            aMatrixValues.data(),
            ADD_VALUES );
}

// ----------------------------------------------------------------------------

void
Matrix_PETSc::get_matrix_values(
        const moris::Matrix< DDSMat >& aRequestedIds,
        moris::Matrix< DDRMat >&       aValues )
{
    // create copies of requested moris IDs
    Matrix< DDSMat > tTempIDs = aRequestedIds;

    // map moris IDs into petsc IDs if map exists
    if ( mMap )
    {
        AOApplicationToPetsc(
                mMap->get_petsc_map(),
                tTempIDs.numel(),
                tTempIDs.data() );
    }

    // get values in row based format; column based format not available
    MatGetValues( mPETScMat,
            tTempIDs.numel(),
            tTempIDs.data(),
            tTempIDs.numel(),
            tTempIDs.data(),
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
}

// ----------------------------------------------------------------------------

void
Matrix_PETSc::dirichlet_BC_vector(
        moris::Matrix< DDUMat >&       aDirichletBCVec,
        const moris::Matrix< DDUMat >& aMyConstraintDofs )
{
    // build vector with constrained dofs: unconstrained=0; constrained =1
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

// ----------------------------------------------------------------------------

void
Matrix_PETSc::build_graph(
        Vector< moris_id >& aNonZeroDiagonal,
        Vector< moris_id >& aNonZeroOffDiagonal )
{
    // Define sparsity structure
    MatMPIAIJSetPreallocation( mPETScMat, 0, aNonZeroDiagonal.memptr(), 0, aNonZeroOffDiagonal.memptr() );

    // ask the matrix to keep the sparsity pattern
    MatSetOption( mPETScMat, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE );

    // Allow colum wise element insertion
    MatSetOption( mPETScMat, MAT_ROW_ORIENTED, PETSC_FALSE );

    // Finalize setup of matrix
    MatSetUp( mPETScMat );

    // note: if there is a problem with this routine turn this option off at the cost of performance
    // MatSetOption(mPETScMat,  MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
}
