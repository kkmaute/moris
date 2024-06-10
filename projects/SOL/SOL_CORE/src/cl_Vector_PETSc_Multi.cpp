/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MultiVector_PETSc.cpp
 *
 */

//VecSetOption( tPetscSingleVec, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE );
//if a vector is extracted ignorenegidx = tPetscSingleVec->stash.ignorenegidx;
// #include <petsc/private/vecimpl.h>

#include "cl_Vector_PETSc_Multi.hpp"
#include <petscviewerhdf5.h>

extern moris::Comm_Manager gMorisComm;

using namespace moris;

// ----------------------------------------------------------------------------

MultiVector_PETSc::MultiVector_PETSc(
        Solver_Interface* aInput,
        sol::Dist_Map*    aMap,
        const sint        aNumVectors,
        bool              aManageMap )
        : sol::Dist_Vector( aManageMap ), mNumVectors( aNumVectors )
{
    // store map as PETSc map
    mMap = reinterpret_cast< Map_PETSc* >( aMap );

    // build either vector of only owned or vector of owned and shared DOFs, i.e., full vector
    if ( mMap->is_full_map() )
    {
        // get number of owned and shared DOFs on this processor
        uint tMyNumOwnedAndSharedDofs = aInput->get_my_local_global_overlapping_map().n_rows();

        MatCreate( PETSC_COMM_SELF, &mPetscVector );
        MatSetSizes( mPetscVector, tMyNumOwnedAndSharedDofs, aNumVectors, PETSC_DECIDE, PETSC_DECIDE );
        MatSetType(mPetscVector, MATDENSE);
        // Allow colum wise element insertion
        MatSetOption( mPetscVector, MAT_ROW_ORIENTED, PETSC_FALSE );
        MatSetUp( mPetscVector);
    }
    else
    {
        // get number of owned Dofs on this processor
        uint tMyNumOwnedDofs = aInput->get_my_local_global_map().n_rows();

        // get owned Dof Ids on this processor
        Matrix< DDSMat > tMyDofIds = aInput->get_my_local_global_map();

        // get constrained Dof Ids on this processor
        Matrix< DDUMat > tMyConstrainedDofIds = aInput->get_constrained_Ids();

        // sum up all owned Dofs
        uint tNumGlobalDofs = sum_all( tMyNumOwnedDofs );

        // FIXME insert boolean array for BC or - better - map
        mDirichletBCVec.set_size( tNumGlobalDofs, 1, 0 );

        // build BC vector
        this->dirichlet_BC_vector( mDirichletBCVec, tMyConstrainedDofIds );


        MatCreate( PETSC_COMM_WORLD, &mPetscVector );
        MatSetSizes( mPetscVector, tMyNumOwnedDofs, aNumVectors, tNumGlobalDofs, aNumVectors );
        MatSetType(mPetscVector, MATDENSE);
         // Allow colum wise element insertion
        MatSetOption( mPetscVector, MAT_ROW_ORIENTED, PETSC_FALSE );
        MatSetUp( mPetscVector);
    }
}

//-----------------------------------------------------------------------------

MultiVector_PETSc::~MultiVector_PETSc()
{
    MatDestroy( &mPetscVector );

    if ( mManageMap )
    {
        delete mMap;
    }
}

//-----------------------------------------------------------------------------

real&
MultiVector_PETSc::operator()( sint aGlobalId, uint aVectorIndex )
{
    MORIS_ERROR( false, "operator() not implemented for PETSc vector." );
    Matrix< DDRMat > tLHSValues;
    extract_copy( tLHSValues );
    return tLHSValues( 0 );
}

//-----------------------------------------------------------------------------

void
MultiVector_PETSc::sum_into_global_values(
        const Matrix< DDSMat >& aGlobalIds,
        const Matrix< DDRMat >& aValues,
        const uint&             aVectorIndex )
{
    // get number of dofs to be summed into
    uint tNumMyDofs = aGlobalIds.numel();

    // check for consistent sizes of vectors of IDs and values
    MORIS_ASSERT( aValues.numel() == tNumMyDofs,
            "MultiVector_PETSc::sum_into_global_values - inconsistent sizes of ID and value vectors" );

    // create copy of vector with moris IDs; will be overwritten in AOApplicationToPetsc
    Matrix< DDSMat > tTempElemDofs = aGlobalIds;

    // loop over elemental dofs
    for ( uint Ij = 0; Ij < tNumMyDofs; Ij++ )
    {
        // set constrDof to neg value
        if ( mDirichletBCVec( aGlobalIds( Ij, 0 ), 0 ) == 1 )
        {
            tTempElemDofs( Ij, 0 ) = -1;
        }
    }

    // map moris IDs into petsc IDs
    AOApplicationToPetsc( mMap->get_petsc_map(), tNumMyDofs, tTempElemDofs.data() );

    // create vector of column indices
    std::vector<sint> tColumnIndices(tNumMyDofs,aVectorIndex);
    
    // add values into matrix, negative values are ignored
    MatSetValues( mPetscVector, tNumMyDofs, tTempElemDofs.data(), 1 , tColumnIndices.data(), aValues.data(), ADD_VALUES );

}

//-----------------------------------------------------------------------------

void
MultiVector_PETSc::sum_into_global_values(
        const Vector< sint >&   aGlobalIds,
        const Matrix< DDRMat >& aValues,
        const uint&             aVectorIndex )
{
    MORIS_ERROR( false, "sum_into_global_values()v2, not implemented for petsc" );
}

//-----------------------------------------------------------------------------


void
MultiVector_PETSc::replace_global_values(
        const moris::Matrix< DDSMat >& aGlobalIds,
        const moris::Matrix< DDRMat >& aValues,
        const uint&                    aVectorIndex )
{
    // get number of dofs to be summed into
    uint tNumMyDofs = aGlobalIds.numel();

    // check for consistent sizes of vectors of IDs and values
    MORIS_ASSERT( aValues.numel() == tNumMyDofs,
            "MultiVector_PETSc::sum_into_global_values - inconsistent sizes of ID and value vectors" );

    // create copy of vector with moris IDs; will be overwritten in AOApplicationToPetsc
    Matrix< DDSMat > tTempElemDofs = aGlobalIds;

    // loop over elemental dofs
    for ( uint Ij = 0; Ij < tNumMyDofs; Ij++ )
    {
        // set constrDof to neg value
        if ( mDirichletBCVec( aGlobalIds( Ij, 0 ), 0 ) == 1 )
        {
            tTempElemDofs( Ij, 0 ) = -1;
        }
    }

    // map moris IDs into petsc IDs
    AOApplicationToPetsc( mMap->get_petsc_map(), tNumMyDofs, tTempElemDofs.data() );

    // create vector of column indices
    std::vector<sint> tColumnIndices(tNumMyDofs,aVectorIndex);
    
    // add values into matrix, negative values are ignored
    MatSetValues( mPetscVector, tNumMyDofs, tTempElemDofs.data(), 1 , tColumnIndices.data(), aValues.data(), INSERT_VALUES );
}

//-----------------------------------------------------------------------------

void MultiVector_PETSc::replace_global_values(
        const Vector< sint >& aGlobalIds,
        const Vector< real >& aValues )
{
    MORIS_ERROR( false, "replace_global_values()v2, not implemented for petsc" );
}

//-----------------------------------------------------------------------------

void
MultiVector_PETSc::dirichlet_BC_vector(
        Matrix< DDUMat >&       aDirichletBCVec,
        const Matrix< DDUMat >& aMyConstraintDofs )
{
    // build vector with constraint values: unconstrained =0 constrained =1
    for ( uint Ik = 0; Ik < aMyConstraintDofs.n_rows(); Ik++ )
    {
        aDirichletBCVec( aMyConstraintDofs( Ik, 0 ), 0 ) = 1;
    }
}

//-----------------------------------------------------------------------------

void
MultiVector_PETSc::vector_global_assembly()
{
    MatAssemblyBegin( mPetscVector, MAT_FINAL_ASSEMBLY );
    MatAssemblyEnd( mPetscVector, MAT_FINAL_ASSEMBLY );
}

//-----------------------------------------------------------------------------

void
MultiVector_PETSc::vec_plus_vec(
        const real&       aScaleA,
        sol::Dist_Vector& aVecA,
        const real&       aScaleThis )
{
    // set scaling value of given vector
    PetscScalar tValueA = aScaleA;

    // set scaling value for member vector
    PetscScalar tValueThis = aScaleThis;

    // if( tValueThis not_eq 1.0 ) 
    MatScale( mPetscVector, tValueThis );

    // FIXME: check if the nonzero parttern can be improved
    MatAXPY( mPetscVector, tValueA, dynamic_cast< MultiVector_PETSc& >( aVecA ).get_petsc_vector(), UNKNOWN_NONZERO_PATTERN );
}

//-----------------------------------------------------------------------------

void
MultiVector_PETSc::scale_vector(
        const real& aValue,
        const uint& aVecIndex )
{
    PetscScalar tValue = aValue;

    MatScale( mPetscVector, tValue );

}

//-----------------------------------------------------------------------------

void
MultiVector_PETSc::vec_put_scalar( const real& aValue )
{
    PetscScalar tValue = aValue;

    MatZeroEntries( mPetscVector );

    if( tValue != 0.0 )
    {
        MatShift( mPetscVector, tValue );
    }
}

//-----------------------------------------------------------------------------

void
MultiVector_PETSc::random()
{
        // Assume that rctx is your PetscRandom context
    PetscRandom rctx;

    // Create the random number context
    PetscRandomCreate(PETSC_COMM_WORLD, &rctx);

    // Set the type of the random number generator
    PetscRandomSetType(rctx, PETSCRAND);

    // Fill the matrix with random numbers
    MatSetRandom(mPetscVector, rctx);

    // Destroy the random number context
    PetscRandomDestroy(&rctx);
}

//-----------------------------------------------------------------------------

sint
MultiVector_PETSc::vec_local_length() const
{
    sint tVecLocSize;
    MatGetLocalSize( mPetscVector, &tVecLocSize, NULL );
    return tVecLocSize;
}

//-----------------------------------------------------------------------------

sint
MultiVector_PETSc::vec_global_length() const
{
    sint tVecSize;
    MatGetSize( mPetscVector, &tVecSize, NULL );
    return tVecSize;
}

//-----------------------------------------------------------------------------

Vector< real >
MultiVector_PETSc::vec_norm2()
{
    Vector< real > tVecNorm( mNumVectors, 0.0 );
    MatGetColumnNorms( mPetscVector, NORM_2, tVecNorm.memptr() );
    return tVecNorm;
}

//-----------------------------------------------------------------------------

void
MultiVector_PETSc::check_vector()
{
    MORIS_ASSERT( false, "epetra vector should not have any input on the petsc vector" );
}

//-----------------------------------------------------------------------------

void
MultiVector_PETSc::extract_copy( Matrix< DDRMat >& LHSValues )
{
    moris::sint tVectorLength = this->vec_local_length();
    LHSValues.set_size( tVectorLength, mNumVectors );

    // loop over the number of multi-vectors and copy each individual vector
    for(uint iVecIndex = 0 ; iVecIndex < mNumVectors ; iVecIndex++ )
    {
        // get a pointer to 
        PetscScalar *Aa;
        MatDenseGetColumn( mPetscVector,iVecIndex,  &Aa );

        // copy the data from Aa to LHSValues
        std::copy( Aa, Aa + tVectorLength, LHSValues.colptr( iVecIndex ) );

        // restore the column
        MatDenseRestoreColumn( mPetscVector, &Aa );
    }

}

//-----------------------------------------------------------------------------

void MultiVector_PETSc::extract_copy( Vector< real >& aVector )
{
    MORIS_ERROR( false, "extract_copy(), not implemented for petsc" );
}

//-----------------------------------------------------------------------------

void
MultiVector_PETSc::import_local_to_global( sol::Dist_Vector& aSourceVec )
{
    // cast source vector to MultiVector_PETSc
    MultiVector_PETSc& tPetscSourceVec = dynamic_cast< MultiVector_PETSc& >( aSourceVec );

    // get raw multi vector of source vector
    Mat tSourceVec = tPetscSourceVec.get_petsc_vector();

    // get petsc map of source vector
    Map_PETSc* tSourceMap = tPetscSourceVec.get_petsc_map();

    // get list of source petsc ids of from map
    IS tPetscSourceIds = tSourceMap->get_petsc_ids();

    // get list of target petsc ids of from map
    IS tPetscTargetIds = mMap->get_petsc_ids();

    // if source vector is full vector and local vector is vector of only owned dofs
    if ( tSourceMap->is_full_map() )
    {
        // check that target vector is owned vector
        MORIS_ERROR( !mMap->is_full_map(),
                "MultiVector_PETSc::import_local_to_global - target and source vectors are both full vectors: case not implemented." );

        // get number of source dofs
        PetscInt tNumSourceIds;
        ISGetLocalSize( tPetscSourceIds, &tNumSourceIds );

        // get list of source petsc IDs
        const PetscInt* tSourceIdList;
        ISGetIndices( tPetscSourceIds, &tSourceIdList );

        // loop over every multi-vector and set values in target vector
        for ( uint iVecIndex = 0; iVecIndex < mNumVectors; iVecIndex++ )
        {
            Vec tSourceVecSingle;
            MatDenseGetColumnVec( tSourceVec, iVecIndex, &tSourceVecSingle );

            Vec mPetscVectorSingle;
            MatDenseGetColumnVec( mPetscVector, iVecIndex, &mPetscVectorSingle );

            // get array from source petsc vector
            PetscScalar* tSourceValues;
            VecGetArray( tSourceVecSingle, &tSourceValues );

            // set values in target vector
            // VecSetValues( mPetscVectorSingle, tNumSourceIds, tSourceIdList, tSourceValues, INSERT_VALUES );
            // create vector of column indices
            std::vector<sint> tColumnIndices(tNumSourceIds,iVecIndex);
            
            // add values into matrix, negative values are ignored
            MatSetValues( mPetscVector, tNumSourceIds, tSourceIdList, 1 , tColumnIndices.data(), tSourceValues, INSERT_VALUES );

            // free memory allocated in petsc calls
            VecRestoreArray( tSourceVecSingle, &tSourceValues );
            MatDenseRestoreColumnVec( tSourceVec, iVecIndex, &tSourceVecSingle );
            MatDenseRestoreColumnVec( mPetscVector, iVecIndex, &mPetscVectorSingle );

            // Flush the assembly
            MatAssemblyBegin( mPetscVector, MAT_FLUSH_ASSEMBLY );
            MatAssemblyEnd( mPetscVector, MAT_FLUSH_ASSEMBLY );
        }

        ISRestoreIndices( tPetscSourceIds, &tSourceIdList );
    }
    else
    {
        // loop over the multivector indices and import the local vector to the global vector
        for ( uint iVecIndex = 0; iVecIndex < mNumVectors; iVecIndex++ )
        {
            Vec tSourceVecSingle;
            MatDenseGetColumnVec( tSourceVec, iVecIndex, &tSourceVecSingle );

            Vec mPetscVectorSingle;
            MatDenseGetColumnVec( mPetscVector, iVecIndex, &mPetscVectorSingle );

            VecScatter tVecScatter;
            VecScatterCreate( tSourceVecSingle, tPetscTargetIds, mPetscVectorSingle, NULL, &tVecScatter );

            // perform scattering
            VecScatterBegin( tVecScatter, tSourceVecSingle, mPetscVectorSingle, INSERT_VALUES, SCATTER_FORWARD );
            VecScatterEnd( tVecScatter, tSourceVecSingle, mPetscVectorSingle, INSERT_VALUES, SCATTER_FORWARD );

            // destroy the scatter vector object
            VecScatterDestroy( &tVecScatter );

            // free the petsc memoery objects
            MatDenseRestoreColumnVec( tSourceVec, iVecIndex, &tSourceVecSingle );
            MatDenseRestoreColumnVec( mPetscVector, iVecIndex, &mPetscVectorSingle );
        }
    }

    this->vector_global_assembly();
}

//-----------------------------------------------------------------------------
// this is only used if the local vector is full vector and source vector is vector of only owned dofs

void
MultiVector_PETSc::import_local_to_global( Vec aSourceVec, uint aVecIndex )
{
    Vec mPetscVectorSingle;
    MatDenseGetColumnVec( mPetscVector, aVecIndex, &mPetscVectorSingle );

    // get list of target petsc ids of from map
    IS tPetscTargetIds = mMap->get_petsc_ids();

    // create scatter object
    VecScatter tVecScatter;
    VecScatterCreate( aSourceVec, tPetscTargetIds, mPetscVectorSingle, NULL, &tVecScatter );

    // perform scattering
    VecScatterBegin( tVecScatter, aSourceVec, mPetscVectorSingle, INSERT_VALUES, SCATTER_FORWARD );
    VecScatterEnd( tVecScatter, aSourceVec, mPetscVectorSingle, INSERT_VALUES, SCATTER_FORWARD );

    // destroy the scatter vector object
    VecScatterDestroy( &tVecScatter );

    // free the petsc memoery objects
    MatDenseRestoreColumnVec( mPetscVector, aVecIndex, &mPetscVectorSingle );
}

//-----------------------------------------------------------------------------

void
MultiVector_PETSc::extract_my_values(
        const uint&               aNumIndices,
        const Matrix< DDSMat >&   aGlobalBlockRows,
        const uint&               aBlockRowOffsets,
        Vector< Matrix< DDRMat > >& ExtractedValues )
{

    ExtractedValues.resize( mNumVectors );

    for ( moris::uint Ik = 0; Ik < mNumVectors; ++Ik )
    {
        ExtractedValues( Ik ).set_size( aNumIndices, 1 );
    }

    // moris::sint tVecLength = this->vec_local_length();

    // get map from moris id to indices in vector
    Matrix< DDSMat > tIndices = mMap->map_from_moris_ids_to_indices( aGlobalBlockRows );

    for ( moris::uint Ik = 0; Ik < mNumVectors; ++Ik )
    {
        Vec mPetscVectorSingle;
        MatDenseGetColumnVec( mPetscVector, Ik, &mPetscVectorSingle );

        // get values from petsc vector
        VecGetValues( mPetscVectorSingle, tIndices.numel(), tIndices.data(), ExtractedValues( Ik ).data() );

        // free the petsc memoery objects
        MatDenseRestoreColumnVec( mPetscVector, Ik, &mPetscVectorSingle );
    }

    // check that vector index  is zero
    MORIS_ASSERT( aBlockRowOffsets == 0,
            "MultiVector_PETSc::extract_my_values - petsc not implemented yet for aBlockRowOffsets neq 0" );
}

//-----------------------------------------------------------------------------

void
MultiVector_PETSc::print() const
{
    MatView( mPetscVector, PETSC_VIEWER_STDOUT_WORLD );
}

//-----------------------------------------------------------------------------

void
MultiVector_PETSc::save_vector_to_HDF5( const char* aFilename )
{
    PetscViewer tViewer;

    PetscViewerHDF5Open( PETSC_COMM_WORLD, aFilename, FILE_MODE_WRITE, &tViewer );

    PetscObjectSetName( (PetscObject)mPetscVector, "Res_Vec" );

    MatView( mPetscVector, tViewer );

    PetscViewerDestroy( &tViewer );
}

//-----------------------------------------------------------------------------

void
MultiVector_PETSc::read_vector_from_HDF5(
        const char* aFilename,
        std::string aGroupName,
        sint        aVectorindex )
{
    PetscViewer tViewer;

    PetscViewerHDF5Open( PETSC_COMM_WORLD, aFilename, FILE_MODE_READ, &tViewer );

    PetscObjectSetName( (PetscObject)mPetscVector, "Res_Vec" );

    MatLoad( mPetscVector, tViewer );

    PetscViewerDestroy( &tViewer );
}

//-----------------------------------------------------------------------------

real*
MultiVector_PETSc::get_values_pointer()
{
    MORIS_ERROR( false, "get_values_pointer() not implemented yet for a PETSc distributed vector." );
    return nullptr;
}

// ----------------------------------------------------------------------------

void
MultiVector_PETSc::save_vector_to_matlab_file( const char* aFilename )
{
    PetscViewer tViewer;

    PetscViewerCreate( PETSC_COMM_WORLD, &tViewer );
    PetscViewerSetType( tViewer, PETSCVIEWERASCII );
    PetscViewerPushFormat( tViewer, PETSC_VIEWER_ASCII_MATLAB );
    PetscViewerFileSetName( tViewer, aFilename );

    MatView( mPetscVector, tViewer );

    PetscViewerDestroy( &tViewer );
}
