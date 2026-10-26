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
#include "cl_Vector_PETSc_Multi.hpp"

#include <cstddef>
#include <cassert>
#include <algorithm>
#include <iostream>

#include "petscmat.h"
#include "cl_SOL_Dist_Map_Custom.hpp"

#include "fn_trans.hpp"

extern moris::Comm_Manager gMorisComm;

// TPL header files
using namespace moris;

// ----------------------------------------------------------------------------

Matrix_PETSc::Matrix_PETSc(
        moris::Solver_Interface* aInput,
        sol::Dist_Map*           aMap,
        const std::string&       aMatrixType )
        : sol::Dist_Matrix( aMap )
{
    // cast the map to the custom map type
    Dist_Map_Custom* tCustomMap = dynamic_cast< Dist_Map_Custom* >( aMap );

    const Matrix< DDSMat >& tOwnedIds = tCustomMap->get_moris_ids_owned();

    // get number of owned dofs
    moris::uint tMyNumOwnedDofs = tOwnedIds.numel();

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

    // Preserve PETSc options overrides while keeping AIJ as the default.
    PetscCallAbort( PETSC_COMM_WORLD, MatSetType( mPETScMat, aMatrixType.c_str() ) );
    PetscCallAbort( PETSC_COMM_WORLD, MatSetFromOptions( mPETScMat ) );

    PetscBool tIsMATIS = PETSC_FALSE;
    PetscCallAbort( PETSC_COMM_WORLD, PetscObjectTypeCompare( (PetscObject)mPETScMat, MATIS, &tIsMATIS ) );

    if ( tIsMATIS )
    {
        // Interface constraints must be known on every contributing subdomain.
        MPI_Allreduce( MPI_IN_PLACE, mDirichletBCVec.data(), tNumGlobalDofs, get_comm_datatype( uint( 0 ) ), MPI_MAX, PETSC_COMM_WORLD );
        mHasConstrainedDofs = std::any_of( mDirichletBCVec.begin(), mDirichletBCVec.end(), []( uint aFlag ) { return aFlag != 0; } );

        // MATIS stores element contributions in the owned-and-shared ordering.
        // Its map translates these subdomain indices to PETSc global IDs.
        Matrix< DDSMat > tGlobalIds = tCustomMap->get_moris_ids_owned_and_shared();
        tCustomMap->map_from_moris_ids_to_petsc_ids( tGlobalIds );
        PetscCallAbort( PETSC_COMM_WORLD, ISLocalToGlobalMappingCreate( PETSC_COMM_WORLD, 1, tGlobalIds.numel(), tGlobalIds.data(), PETSC_COPY_VALUES, &mSubdomainMap ) );
        PetscCallAbort( PETSC_COMM_WORLD, MatSetLocalToGlobalMapping( mPETScMat, mSubdomainMap, mSubdomainMap ) );

        // MORIS element matrices use column-major storage.
        PetscCallAbort( PETSC_COMM_WORLD, MatSetOption( mPETScMat, MAT_ROW_ORIENTED, PETSC_FALSE ) );

        // The element-graph path grows the local graph dynamically. The sparsity
        // path below replaces this with preallocation based on global row counts.
        Mat tLocal;
        PetscCallAbort( PETSC_COMM_WORLD, MatISGetLocalMat( mPETScMat, &tLocal ) );
        PetscCallAbort( PETSC_COMM_WORLD, MatSetUp( tLocal ) );
        PetscCallAbort( PETSC_COMM_WORLD, MatSetOption( tLocal, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE ) );
    }

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
    MatMPIAIJSetPreallocation( mPETScMat, tNonZeros, nullptr, tNonZeros, nullptr );

    // Finalize setup of matrix
    MatSetUp( mPETScMat );
}

// ----------------------------------------------------------------------------

Matrix_PETSc::~Matrix_PETSc()
{
    MatDestroy( &mPETScMat );
    ISLocalToGlobalMappingDestroy( &mSubdomainMap );
}

// ----------------------------------------------------------------------------

void Matrix_PETSc::build_graph(
        const moris::uint&             aNumMyDof,
        const moris::Matrix< DDSMat >& aElementTopology )
{
    // MORIS_ASSERT( false, "Matrix_PETSc::build_graph - not implemented yet" );
    moris::Matrix< DDSMat > tTempElemDofs = aElementTopology;

    // Build Zero matrix and matrix for element free dof id
    moris::Matrix< DDRMat > tZeros( aNumMyDof * aNumMyDof, 1, 0.0 );

    Dist_Map_Custom* tCustomMap = dynamic_cast< Dist_Map_Custom* >( mMap );
    tCustomMap->map_from_moris_ids_to_petsc_ids( tTempElemDofs );

    // Mask constrained or invalid PETSc IDs before constructing the graph.
    for ( auto& tId : tTempElemDofs )
    {
        if ( tId < 0 || tId >= (sint)mDirichletBCVec.numel() || mDirichletBCVec( tId ) )
        {
            tId = -1;
        }
    }

    if ( mSubdomainMap )
    {
        // MATIS assembles into the overlapping subdomain matrix.
        this->map_to_subdomain( tTempElemDofs );
        PetscCallAbort( PETSC_COMM_WORLD, MatSetValuesLocal( mPETScMat, aNumMyDof, tTempElemDofs.data(), aNumMyDof, tTempElemDofs.data(), tZeros.data(), ADD_VALUES ) );
    }
    else
    {
        MatSetValues( mPETScMat, aNumMyDof, tTempElemDofs.data(), aNumMyDof, tTempElemDofs.data(), tZeros.data(), ADD_VALUES );
    }
}

// ----------------------------------------------------------------------------

void Matrix_PETSc::fill_matrix(
        const moris::uint&             aNumMyDof,
        const moris::Matrix< DDRMat >& aA_val,
        const moris::Matrix< DDSMat >& aEleDofConnectivity )
{

    // check for consistent sizes of vectors of IDs and values
    MORIS_ASSERT( aEleDofConnectivity.numel() == aNumMyDof,
            "Matrix_PETSc::fill_matrix - inconsistent sizes of ID and value vectors" );

    // create temporary matrix for element dofs
    Matrix< DDSMat > tTempElemDofs = aEleDofConnectivity;

    Dist_Map_Custom* tCustomMap = dynamic_cast< Dist_Map_Custom* >( mMap );
    tCustomMap->map_from_moris_ids_to_petsc_ids( tTempElemDofs );

    // loop over elemental dofs
    // FIXME: This is inficent, the correct way to do it would be
    // to set the dirichelet values after the assembely is done, it would decrase FLOP
    for ( moris::uint Ij = 0; Ij < aNumMyDof; Ij++ )
    {
        // set constrDof to neg value
        if ( tTempElemDofs( Ij ) < 0 || mDirichletBCVec( tTempElemDofs( Ij ), 0 ) == 1 )
        {
            tTempElemDofs( Ij, 0 ) = -1;
        }
    }

    if ( mSubdomainMap )
    {
        // MATIS assembles into the overlapping subdomain matrix.
        this->map_to_subdomain( tTempElemDofs );
        PetscCallAbort( PETSC_COMM_WORLD, MatSetValuesLocal( mPETScMat, aNumMyDof, tTempElemDofs.data(), aNumMyDof, tTempElemDofs.data(), aA_val.data(), ADD_VALUES ) );
        return;
    }

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

void Matrix_PETSc::insert_values(
        const Matrix< DDSMat >& aRowIDs,
        const Matrix< DDSMat >& aColumnIDs,
        const Matrix< DDRMat >& aMatrixValues )
{
    //    MORIS_ASSERT(false, "Matrix_PETSc::insert_values - not implemented yet");
    // create copies of moris IDs
    Matrix< DDSMat > tTempRowIDs    = aRowIDs;
    Matrix< DDSMat > tTempColumnIDs = aColumnIDs;

    // map moris IDs into petsc IDs if map exists
    if ( mMap )
    {

        Dist_Map_Custom* tCustomMap = dynamic_cast< Dist_Map_Custom* >( mMap );
        tCustomMap->map_from_moris_ids_to_petsc_ids( tTempRowIDs );

        tCustomMap->map_from_moris_ids_to_petsc_ids( tTempColumnIDs );
    }

    if ( mSubdomainMap )
    {
        // Convert both index sets from PETSc global to subdomain numbering.
        this->map_to_subdomain( tTempRowIDs );
        this->map_to_subdomain( tTempColumnIDs );
        PetscCallAbort( PETSC_COMM_WORLD, MatSetValuesLocal( mPETScMat, tTempRowIDs.numel(), tTempRowIDs.data(), tTempColumnIDs.numel(), tTempColumnIDs.data(), aMatrixValues.data(), INSERT_VALUES ) );
        return;
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

void Matrix_PETSc::sum_into_values(
        const Matrix< DDSMat >& aRowIDs,
        const Matrix< DDSMat >& aColumnIDs,
        const Matrix< DDRMat >& aMatrixValues )
{
    MORIS_ERROR( false, "Matrix_PETSc::sum_into_values - not implemented yet" );
    // create copies of moris IDs
    //     Matrix< DDSMat > tTempRowIDs    = aRowIDs;
    //     Matrix< DDSMat > tTempColumnIDs = aColumnIDs;

    //     // map moris IDs into petsc IDs if map exists
    //     if ( mMap )
    //     {
    //         AOApplicationToPetsc(
    //                 mMap->get_petsc_map(),
    //                 tTempRowIDs.numel(),
    //                 tTempRowIDs.data() );

    //         AOApplicationToPetsc(
    //                 mMap->get_petsc_map(),
    //                 tTempColumnIDs.numel(),
    //                 tTempColumnIDs.data() );
    //     }

    //     // insert values into matrix
    //     MatSetValues( mPETScMat,
    //             tTempRowIDs.numel(),
    //             tTempRowIDs.data(),
    //             tTempColumnIDs.numel(),
    //             tTempColumnIDs.data(),
    //             aMatrixValues.data(),
    //             ADD_VALUES );
}

// ----------------------------------------------------------------------------

void Matrix_PETSc::get_matrix_values(
        const moris::Matrix< DDSMat >& aRequestedIds,
        moris::Matrix< DDRMat >&       aValues )
{
    // MORIS_ERROR(false, "Matrix_PETSc::get_matrix_values - not implemented yet");
    //     // create copies of requested moris IDs
    Matrix< DDSMat > tTempIDs = aRequestedIds;

    //     // map moris IDs into petsc IDs if map exists
    if ( mMap )
    {
        Dist_Map_Custom* tCustomMap = dynamic_cast< Dist_Map_Custom* >( mMap );
        tCustomMap->map_from_moris_ids_to_petsc_ids( tTempIDs );
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

void Matrix_PETSc::matrix_global_assembly()
{
    MatAssemblyBegin( mPETScMat, MAT_FINAL_ASSEMBLY );
    MatAssemblyEnd( mPETScMat, MAT_FINAL_ASSEMBLY );

    if ( mSubdomainMap && mHasConstrainedDofs )
    {
        // PETSc's row/column elimination requires an explicit diagonal even
        // where MORIS omitted every constrained element entry. Insert structural
        // zeros after assembly so INSERT_VALUES callers do not mix insert modes.
        auto tLocalConstraints = this->get_constrained_subdomain_dofs();
        if ( !tLocalConstraints.empty() )
        {
            Mat tLocal;
            PetscCallAbort( PETSC_COMM_WORLD, MatISGetLocalMat( mPETScMat, &tLocal ) );
            PetscCallAbort( PETSC_COMM_WORLD, MatSetOption( tLocal, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE ) );

            for ( PetscInt tIndex : tLocalConstraints )
            {
                PetscScalar tZero = 0.0;
                PetscCallAbort( PETSC_COMM_WORLD, MatSetValuesLocal( mPETScMat, 1, &tIndex, 1, &tIndex, &tZero, ADD_VALUES ) );
            }
        }

        // Finish diagonal insertion collectively, including ranks without constraints.
        PetscCallAbort( PETSC_COMM_WORLD, MatAssemblyBegin( mPETScMat, MAT_FINAL_ASSEMBLY ) );
        PetscCallAbort( PETSC_COMM_WORLD, MatAssemblyEnd( mPETScMat, MAT_FINAL_ASSEMBLY ) );

        // Supply each retained constrained global row from its owning rank.
        PetscInt tStart, tEnd;
        PetscCallAbort( PETSC_COMM_WORLD, MatGetOwnershipRange( mPETScMat, &tStart, &tEnd ) );

        std::vector< PetscInt > tConstrained;

        for ( PetscInt i = tStart; i < tEnd; ++i )
        {
            if ( mDirichletBCVec( i ) )
            {
                tConstrained.push_back( i );
            }
        }

        // MORIS skips constrained rows/columns during element assembly. Supply
        // unit rows for these retained global DOFs rather than singular zero rows.
        PetscCallAbort( PETSC_COMM_WORLD, MatZeroRowsColumns( mPETScMat, tConstrained.size(), tConstrained.data(), 1.0, nullptr, nullptr ) );
    }
}

// ----------------------------------------------------------------------------

void Matrix_PETSc::dirichlet_BC_vector(
        moris::Matrix< DDUMat >&       aDirichletBCVec,
        const moris::Matrix< DDUMat >& aMyConstraintDofs )
{

    // create temporary matrix for element dofs
    Matrix< DDUMat > tTempIDs = aMyConstraintDofs;

    Dist_Map_Custom* tCustomMap = dynamic_cast< Dist_Map_Custom* >( mMap );
    tCustomMap->map_from_moris_ids_to_petsc_ids( tTempIDs );

    // build vector with constrained dofs: unconstrained=0; constrained =1
    for ( moris::uint Ik = 0; Ik < aMyConstraintDofs.n_rows(); Ik++ )
    {
        aDirichletBCVec( tTempIDs( Ik, 0 ), 0 ) = 1;
    }
}

// ----------------------------------------------------------------------------

void Matrix_PETSc::print() const
{
    MatView( mPETScMat, PETSC_VIEWER_STDOUT_( PETSC_COMM_WORLD ) );
}

// ----------------------------------------------------------------------------

void Matrix_PETSc::save_matrix_to_matlab_file( const char* aFilename )
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

void Matrix_PETSc::build_graph(
        Vector< moris_id >& aNonZeroDiagonal,
        Vector< moris_id >& aNonZeroOffDiagonal )
{
    // Define sparsity structure
    if ( mSubdomainMap )
    {
        // MATIS scatters owned-row counts to the overlapping local matrices.
        PetscCallAbort( PETSC_COMM_WORLD, MatISSetPreallocation( mPETScMat, 0, aNonZeroDiagonal.memptr(), 0, aNonZeroOffDiagonal.memptr() ) );
    }
    else
    {
        MatMPIAIJSetPreallocation( mPETScMat, 0, aNonZeroDiagonal.memptr(), 0, aNonZeroOffDiagonal.memptr() );
    }

    // ask the matrix to keep the sparsity pattern
    MatSetOption( mPETScMat, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE );

    // Allow colum wise element insertion
    MatSetOption( mPETScMat, MAT_ROW_ORIENTED, PETSC_FALSE );

    // Finalize setup of matrix
    MatSetUp( mPETScMat );

    // note: if there is a problem with this routine turn this option off at the cost of performance
    // MatSetOption(mPETScMat,  MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
}

// ----------------------------------------------------------------------------

void Matrix_PETSc::mat_vec_product(
        const moris::sol::Dist_Vector& aInputVec,
        moris::sol::Dist_Vector&       aResult,
        const bool                     aUseTranspose )
{
    // check if the bool value if false
    MORIS_ASSERT( aUseTranspose == false, "Sparse_Matrix_EpetraFECrs::mat_vec_product: aUseTranspose must be false" );

    // case the input and output vector into the petsc multivector(Mat) object
    //  cast source vector to MultiVector_PETSc
    const MultiVector_PETSc& tPetscSourceVec      = dynamic_cast< const MultiVector_PETSc& >( aInputVec );
    Mat                      tPetscSourceMultiVec = tPetscSourceVec.get_petsc_vector();

    MultiVector_PETSc& tPetscResultVec      = dynamic_cast< MultiVector_PETSc& >( aResult );
    Mat                tPetscResultMultiVec = tPetscResultVec.get_petsc_vector();

    if ( mSubdomainMap )
    {
        // MATIS implements MatMult; dense matrix products are not available in
        // all supported PETSc versions. Apply each RHS column independently.
        PetscInt tNumColumns;
        PetscCallAbort( PETSC_COMM_WORLD, MatGetSize( tPetscSourceMultiVec, nullptr, &tNumColumns ) );

        for ( PetscInt i = 0; i < tNumColumns; ++i )
        {
            Vec tSource, tResult;
            PetscCallAbort( PETSC_COMM_WORLD, MatDenseGetColumnVecRead( tPetscSourceMultiVec, i, &tSource ) );
            PetscCallAbort( PETSC_COMM_WORLD, MatDenseGetColumnVecWrite( tPetscResultMultiVec, i, &tResult ) );

            PetscCallAbort( PETSC_COMM_WORLD, MatMult( mPETScMat, tSource, tResult ) );

            PetscCallAbort( PETSC_COMM_WORLD, MatDenseRestoreColumnVecRead( tPetscSourceMultiVec, i, &tSource ) );
            PetscCallAbort( PETSC_COMM_WORLD, MatDenseRestoreColumnVecWrite( tPetscResultMultiVec, i, &tResult ) );
        }
        return;
    }

    // perform the matrix vector multiplication
    MatMatMult( mPETScMat, tPetscSourceMultiVec, MAT_REUSE_MATRIX, PETSC_DEFAULT, &tPetscResultMultiVec );
}

// ----------------------------------------------------------------------------

void Matrix_PETSc::map_to_subdomain( Matrix< DDSMat >& aIds ) const
{
    // Keep the input shape and negative sentinels used to skip constrained entries.
    Matrix< DDSMat > tLocalIds( aIds.n_rows(), aIds.n_cols() );
    PetscCallAbort( PETSC_COMM_WORLD, ISGlobalToLocalMappingApply( mSubdomainMap, IS_GTOLM_MASK, aIds.numel(), aIds.data(), nullptr, tLocalIds.data() ) );

    // A nonnegative global ID must belong to this rank's overlap map.
    for ( uint i = 0; i < aIds.numel(); ++i )
    {
        MORIS_ERROR( aIds( i ) < 0 || tLocalIds( i ) >= 0,
                "MATIS assembly references a DOF absent from the owned-and-shared map." );
    }

    aIds = tLocalIds;
}

// ----------------------------------------------------------------------------

std::vector< PetscInt > Matrix_PETSc::get_constrained_subdomain_dofs() const
{
    std::vector< PetscInt > tDofs;

    if ( !mSubdomainMap )
    {
        return tDofs;
    }

    // Boundary flags use PETSc global IDs; BDDC needs subdomain indices.
    PetscInt        tSize;
    const PetscInt* tIds;
    PetscCallAbort( PETSC_COMM_WORLD, ISLocalToGlobalMappingGetSize( mSubdomainMap, &tSize ) );
    PetscCallAbort( PETSC_COMM_WORLD, ISLocalToGlobalMappingGetIndices( mSubdomainMap, &tIds ) );

    for ( PetscInt i = 0; i < tSize; ++i )
    {
        if ( mDirichletBCVec( tIds[ i ] ) )
        {
            tDofs.push_back( i );
        }
    }

    PetscCallAbort( PETSC_COMM_WORLD, ISLocalToGlobalMappingRestoreIndices( mSubdomainMap, &tIds ) );

    return tDofs;
}
