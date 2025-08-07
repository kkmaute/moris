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
    mMap = reinterpret_cast< Dist_Map_Custom* >( aMap );

    // build either vector of only owned or vector of owned and shared DOFs, i.e., full vector
    if ( mMap->is_full_map() )
    {
        // get number of owned and shared DOFs on this processor
        uint tMyNumOwnedAndSharedDofs = mMap->get_moris_ids_owned_and_shared().n_rows();

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
        uint tMyNumOwnedDofs = mMap->get_moris_ids_owned().n_rows();

        // get owned Dof Ids on this processor
        Matrix< DDSMat > tMyDofIds = mMap->get_moris_ids_owned();

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

    Dist_Map_Custom* tCustomMap = mMap;
    tCustomMap->map_from_moris_ids_to_petsc_ids( tTempElemDofs);

    // loop over elemental dofs
    for ( uint Ij = 0; Ij < tNumMyDofs; Ij++ )
    {
        
        // set constrDof to neg value
        if ( mDirichletBCVec( tTempElemDofs( Ij, 0 ), 0 ) == 1 )
        {
            tTempElemDofs( Ij, 0 ) = -1;
        }
    }

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

    Dist_Map_Custom* tCustomMap = mMap;
    tCustomMap->map_from_moris_ids_to_petsc_ids( tTempElemDofs);

    // loop over elemental dofs
    for ( uint Ij = 0; Ij < tNumMyDofs; Ij++ )
    {
        // set constrDof to neg value
        if ( mDirichletBCVec( aGlobalIds( Ij, 0 ), 0 ) == 1 )
        {
            tTempElemDofs( Ij, 0 ) = -1;
        }
    }

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
    MatGetLocalSize( mPetscVector, &tVecLocSize, nullptr );
    return tVecLocSize;
}

//-----------------------------------------------------------------------------

sint
MultiVector_PETSc::vec_global_length() const
{
    sint tVecSize;
    MatGetSize( mPetscVector, &tVecSize, nullptr );
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
        // MORIS_ASSERT( false, "import_local_to_global() not implemented for petsc" );
    // cast source vector to MultiVector_PETSc
    MultiVector_PETSc& tPetscSourceVec = dynamic_cast< MultiVector_PETSc& >( aSourceVec );

    // get raw multi vector of source vector
    Mat tSourceVec = tPetscSourceVec.get_petsc_vector();

    // // get petsc map of source vector
    Dist_Map_Custom* tSourceMap = dynamic_cast<Dist_Map_Custom*>(tPetscSourceVec.get_map());

    // the case where one imnports a full vector(sequntail) to a paralle vector
    if ( tSourceMap->is_full_map() and !mMap->is_full_map() ) 
    {
        // create the IS objects to do the scatter
        IS tFrom,tTo; 

        // the data is transfred from the parallel owned and shared 
        Matrix< DDSMat > tOwnedMorisIds = mMap->get_moris_ids_owned();
        Matrix< DDSMat > tSequantalVectorIndex = tSourceMap->map_from_moris_ids_to_indices( tOwnedMorisIds );
        ISCreateGeneral( PETSC_COMM_WORLD, tSequantalVectorIndex.numel(), tSequantalVectorIndex.data(), PETSC_USE_POINTER, &tFrom );

        //
        Matrix< DDSMat > tOwnedMorisIdsDest    = mMap->get_moris_ids_owned();
        mMap->map_from_moris_ids_to_petsc_ids( tOwnedMorisIdsDest );
        ISCreateGeneral( PETSC_COMM_WORLD, tOwnedMorisIdsDest.numel(), tOwnedMorisIdsDest.data(), PETSC_USE_POINTER, &tTo );

        for ( uint iVecIndex = 0; iVecIndex < mNumVectors; iVecIndex++ )
        {
            Vec tSourceVecSingle;
            MatDenseGetColumnVec( tSourceVec, iVecIndex, &tSourceVecSingle );

            Vec mPetscVectorSingle;
            MatDenseGetColumnVec( mPetscVector, iVecIndex, &mPetscVectorSingle );

            VecScatter tVecScatter;
            VecScatterCreate( tSourceVecSingle, tFrom, mPetscVectorSingle, tTo, &tVecScatter );

            // perform scattering
            VecScatterBegin( tVecScatter, tSourceVecSingle, mPetscVectorSingle, INSERT_VALUES, SCATTER_FORWARD );
            VecScatterEnd( tVecScatter, tSourceVecSingle, mPetscVectorSingle, INSERT_VALUES, SCATTER_FORWARD );

            // destroy the scatter vector object
            VecScatterDestroy( &tVecScatter );

            // free the petsc memoery objects
            MatDenseRestoreColumnVec( tSourceVec, iVecIndex, &tSourceVecSingle );
            MatDenseRestoreColumnVec( mPetscVector, iVecIndex, &mPetscVectorSingle );
        }

        ISDestroy( &tFrom );
        ISDestroy( &tTo );

    }

    // the case where one imnports a parallel vector(source) to a full vector(sequntail)
    if ( mMap->is_full_map() and !tSourceMap->is_full_map() )
    {
        // create the IS objects to do the scatter
        IS tFrom,tTo; 

        // the data is transfred from the parallel owned and shared 
        Matrix< DDSMat > tOwnedMorisIds = tSourceMap->get_moris_ids_owned_and_shared();
        tSourceMap->map_from_moris_ids_to_petsc_ids( tOwnedMorisIds );
        ISCreateGeneral( PETSC_COMM_WORLD, tOwnedMorisIds.numel(), tOwnedMorisIds.data(), PETSC_USE_POINTER, &tFrom );

        //
        Matrix< DDSMat > tOwnedMorisIdsDest    = tSourceMap->get_moris_ids_owned_and_shared();  // changed from mMap
        Matrix< DDSMat > tSequantalVectorIndex = mMap->map_from_moris_ids_to_indices( tOwnedMorisIdsDest );
        ISCreateGeneral( PETSC_COMM_WORLD, tSequantalVectorIndex.numel(), tSequantalVectorIndex.data(), PETSC_USE_POINTER, &tTo );

        for ( uint iVecIndex = 0; iVecIndex < mNumVectors; iVecIndex++ )
        {
            Vec tSourceVecSingle;
            MatDenseGetColumnVec( tSourceVec, iVecIndex, &tSourceVecSingle );

            Vec mPetscVectorSingle;
            MatDenseGetColumnVec( mPetscVector, iVecIndex, &mPetscVectorSingle );

            VecScatter tVecScatter;
            VecScatterCreate( tSourceVecSingle, tFrom, mPetscVectorSingle, tTo, &tVecScatter );

            // perform scattering
            VecScatterBegin( tVecScatter, tSourceVecSingle, mPetscVectorSingle, INSERT_VALUES, SCATTER_FORWARD );
            VecScatterEnd( tVecScatter, tSourceVecSingle, mPetscVectorSingle, INSERT_VALUES, SCATTER_FORWARD );

            // destroy the scatter vector object
            VecScatterDestroy( &tVecScatter );

            // free the petsc memoery objects
            MatDenseRestoreColumnVec( tSourceVec, iVecIndex, &tSourceVecSingle );
            MatDenseRestoreColumnVec( mPetscVector, iVecIndex, &mPetscVectorSingle );
        }

        ISDestroy( &tFrom );
        ISDestroy( &tTo );
    }

}

//-----------------------------------------------------------------------------
// this is only used if the local vector is full vector and source vector is vector of only owned dofs

void MultiVector_PETSc::import_local_to_global( Vec aSourceVec, uint aVecIndex, Dist_Map_Custom* tSourceMap )
{

    // create the IS objects to do the scatter
    IS tFrom, tTo;

    // the data is transfred from the parallel owned and shared
    Matrix< DDSMat > tOwnedMorisIds = tSourceMap->get_moris_ids_owned_and_shared();
    tSourceMap->map_from_moris_ids_to_petsc_ids( tOwnedMorisIds );
    ISCreateGeneral( PETSC_COMM_WORLD, tOwnedMorisIds.numel(), tOwnedMorisIds.data(), PETSC_USE_POINTER, &tFrom );

    //
    Matrix< DDSMat > tOwnedMorisIdsDest    = tSourceMap->get_moris_ids_owned_and_shared();    // changed from mMap
    Matrix< DDSMat > tSequantalVectorIndex = mMap->map_from_moris_ids_to_indices( tOwnedMorisIdsDest );
    ISCreateGeneral( PETSC_COMM_WORLD, tSequantalVectorIndex.numel(), tSequantalVectorIndex.data(), PETSC_USE_POINTER, &tTo );

    Vec mPetscVectorSingle;
    MatDenseGetColumnVec( mPetscVector, aVecIndex, &mPetscVectorSingle );

    VecScatter tVecScatter;
    VecScatterCreate( aSourceVec, tFrom, mPetscVectorSingle, tTo, &tVecScatter );

    // perform scattering
    VecScatterBegin( tVecScatter, aSourceVec, mPetscVectorSingle, INSERT_VALUES, SCATTER_FORWARD );
    VecScatterEnd( tVecScatter, aSourceVec, mPetscVectorSingle, INSERT_VALUES, SCATTER_FORWARD );

    // destroy the scatter vector object
    VecScatterDestroy( &tVecScatter );

    // free the petsc memoery objects
    MatDenseRestoreColumnVec( mPetscVector, aVecIndex, &mPetscVectorSingle );

    ISDestroy( &tFrom );
    ISDestroy( &tTo );
}

//-----------------------------------------------------------------------------

void
MultiVector_PETSc::extract_my_values(
        const uint&               aNumIndices,
        const Matrix< DDSMat >&   aGlobalBlockRows,
        const uint&               aBlockRowOffsets,
        Vector< Matrix< DDRMat > >& ExtractedValues )
{
    // MORIS_ASSERT( false, "extract_my_values() not implemented for petsc" );
    ExtractedValues.resize( mNumVectors );

    for ( moris::uint Ik = 0; Ik < mNumVectors; ++Ik )
    {
        ExtractedValues( Ik ).set_size( aNumIndices, 1 );
    }

    // moris::sint tVecLength = this->vec_local_length();

    // get map from moris id to indices in vector
    MORIS_ASSERT( mMap->is_full_map(), "MultiVector_PETSc::extract_my_values - full map required" );

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
