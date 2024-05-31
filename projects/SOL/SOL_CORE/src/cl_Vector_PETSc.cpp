/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Vector_PETSc.cpp
 *
 */

#include "cl_Vector_PETSc.hpp"

#include <petscviewerhdf5.h>

extern moris::Comm_Manager gMorisComm;

using namespace moris;

// ----------------------------------------------------------------------------

Vector_PETSc::Vector_PETSc(
        Solver_Interface* aInput,
        sol::Dist_Map*    aMap,
        const sint        aNumVectors,
        bool              aManageMap )
        : sol::Dist_Vector( aManageMap )
{
    // store map as PETSc map
    mMap = reinterpret_cast< Map_PETSc* >( aMap );

    // build either vector of only owned or vector of owned and shared DOFs, i.e., full vector
    if ( mMap->is_full_map() )
    {
        // get number of owned and shared DOFs on this processor
        uint tMyNumOwnedAndSharedDofs = aInput->get_my_local_global_overlapping_map().n_rows();

        // build sequential vector
        VecCreateSeq( PETSC_COMM_SELF, tMyNumOwnedAndSharedDofs, &mPetscVector );
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

        // build distributed vector
        VecCreateMPI( PETSC_COMM_WORLD, tMyNumOwnedDofs, tNumGlobalDofs, &mPetscVector );

        // set option that negative IDs are ignored
        VecSetOption( mPetscVector, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE );
        VecSetFromOptions( mPetscVector );

        // finalize setup of vector
        VecSetUp( mPetscVector );
    }
}

//-----------------------------------------------------------------------------

Vector_PETSc::~Vector_PETSc()
{
    VecDestroy( &mPetscVector );

    if ( mManageMap )
    {
        delete mMap;
    }
}

//-----------------------------------------------------------------------------

real&
Vector_PETSc::operator()( sint aGlobalId, uint aVectorIndex )
{
    MORIS_ERROR( false, "operator() not implemented for PETSc vector." );
    Matrix< DDRMat > tLHSValues;
    extract_copy( tLHSValues );
    return tLHSValues( 0 );
}

//-----------------------------------------------------------------------------

void
Vector_PETSc::sum_into_global_values(
        const Matrix< DDSMat >& aGlobalIds,
        const Matrix< DDRMat >& aValues,
        const uint&             aVectorIndex )
{
    // Petsc implementation is for single column vector as of now
    MORIS_ASSERT( aVectorIndex == 0,
            "Vector_PETSc::sum_into_global_values - Petsc implementation is for single column vector only." );

    // get number of dofs to be summed into
    uint tNumMyDofs = aGlobalIds.numel();

    // check for consistent sizes of vectors of IDs and values
    MORIS_ASSERT( aValues.numel() == tNumMyDofs,
            "Vector_PETSc::sum_into_global_values - inconsistent sizes of ID and value vectors" );

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

    // add values into vector
    VecSetValues( mPetscVector, tNumMyDofs, tTempElemDofs.data(), aValues.data(), ADD_VALUES );
}

//-----------------------------------------------------------------------------

void
Vector_PETSc::sum_into_global_values(
        const Vector< sint >&   aGlobalIds,
        const Matrix< DDRMat >& aValues,
        const uint&             aVectorIndex )
{
    // Petsc implementation is for single column vector as of now
    MORIS_ASSERT( aVectorIndex == 0,
            "Vector_PETSc::sum_into_global_values - Petsc implementation is for single column vector only." );

    // get number of dofs to be summed into
    uint tNumMyDofs = aGlobalIds.size();

    // check for consistent sizes of vectors of IDs and values
    MORIS_ASSERT( aValues.numel() == tNumMyDofs,
            "Vector_PETSc::sum_into_global_values - inconsistent sizes of ID and value vectors" );

    // create copy of vector with moris IDs; will be overwritten in AOApplicationToPetsc
    Vector< sint > tTempElemDofs = aGlobalIds;

    // loop over elemental dofs
    for ( uint Ij = 0; Ij < tNumMyDofs; Ij++ )
    {
        // set constrDof to neg value
        if ( mDirichletBCVec( aGlobalIds( Ij ), 0 ) == 1 )
        {
            tTempElemDofs( Ij ) = -1;
        }
    }

    // map moris IDs into petsc IDs
    AOApplicationToPetsc( mMap->get_petsc_map(), tNumMyDofs, tTempElemDofs.memptr() );

    // add values into vector
    VecSetValues( mPetscVector, tNumMyDofs, tTempElemDofs.memptr(), aValues.data(), ADD_VALUES );
}

//-----------------------------------------------------------------------------


void
Vector_PETSc::replace_global_values(
        const moris::Matrix< DDSMat >& aGlobalIds,
        const moris::Matrix< DDRMat >& aValues,
        const uint&                    aVectorIndex )
{
    // check that Id and value vectors have same length
    MORIS_ASSERT( aGlobalIds.numel() == aValues.numel(),
            "Vector_PETSc::replace_global_values - inconsistent number of IDs and values" );

    // check that vector index  is zero
    MORIS_ASSERT( aVectorIndex == 0,
            "Vector_PETSc::replace_global_values - petsc not implemented for multi-vectors yet" );

    // get map from moris id to indices in vector
    Matrix< DDSMat > tIndices = mMap->map_from_moris_ids_to_indices( aGlobalIds );

    // set values in petsc vector
    VecSetValuesLocal( mPetscVector, tIndices.numel(), tIndices.data(), aValues.data(), INSERT_VALUES );
}

//-----------------------------------------------------------------------------

void Vector_PETSc::replace_global_values(
        const Vector< sint >& aGlobalIds,
        const Vector< real >& aValues )
{
    // check that Id and value vectors have same length
    MORIS_ASSERT( aGlobalIds.size() == aValues.size(),
            "Vector_PETSc::replace_global_values - inconsistent number of IDs and values" );

    // get map from moris id to indices in vector
    Vector< sint > tIndices = mMap->map_from_moris_ids_to_indices( aGlobalIds );

    // set values in petsc vector
    VecSetValuesLocal( mPetscVector, tIndices.size(), tIndices.memptr(), aValues.memptr(), INSERT_VALUES );
}

//-----------------------------------------------------------------------------

void
Vector_PETSc::dirichlet_BC_vector(
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
Vector_PETSc::vector_global_assembly()
{
    VecAssemblyBegin( mPetscVector );
    VecAssemblyEnd( mPetscVector );
}

//-----------------------------------------------------------------------------

void
Vector_PETSc::vec_plus_vec(
        const real&       aScaleA,
        sol::Dist_Vector& aVecA,
        const real&       aScaleThis )
{
    // set scaling value of given vector
    PetscScalar tValueA = aScaleA;

    // set scaling value for member vector
    PetscScalar tValueThis = aScaleThis;

    // perform operation
    VecAXPBY( mPetscVector, tValueA, tValueThis, dynamic_cast< Vector_PETSc& >( aVecA ).get_petsc_vector() );
}

//-----------------------------------------------------------------------------

void
Vector_PETSc::scale_vector(
        const real& aValue,
        const uint& aVecIndex )
{
    PetscScalar tValue = aValue;

    VecScale( mPetscVector, tValue );
}

//-----------------------------------------------------------------------------

void
Vector_PETSc::vec_put_scalar( const real& aValue )
{
    PetscScalar tValue = aValue;

    VecSet( mPetscVector, tValue );
}

//-----------------------------------------------------------------------------

void
Vector_PETSc::random()
{
    MORIS_ERROR( false, "random(), not implemented for petsc" );
}

//-----------------------------------------------------------------------------

sint
Vector_PETSc::vec_local_length() const
{
    sint tVecLocSize = 0;
    VecGetLocalSize( mPetscVector, &tVecLocSize );
    return tVecLocSize;
}

//-----------------------------------------------------------------------------

sint
Vector_PETSc::vec_global_length() const
{
    sint tVecSize = 0;
    VecGetSize( mPetscVector, &tVecSize );
    return tVecSize;
}

//-----------------------------------------------------------------------------

Vector< real >
Vector_PETSc::vec_norm2()
{
    Vector< real > tVecNorm( mNumVectors, 0.0 );

    VecNorm( mPetscVector, NORM_2, tVecNorm.data().data() );
    return tVecNorm;
}

//-----------------------------------------------------------------------------

void
Vector_PETSc::check_vector()
{
    MORIS_ASSERT( false, "epetra vector should not have any input on the petsc vector" );
}

//-----------------------------------------------------------------------------

void
Vector_PETSc::extract_copy( Matrix< DDRMat >& LHSValues )
{
    // VecGetArray (tSolution, &  LHSValues.data());

    sint tVecLocSize;
    VecGetLocalSize( mPetscVector, &tVecLocSize );

    // FIXME replace with VecGetArray()
    Matrix< DDSMat > tVal( tVecLocSize, 1, 0 );
    LHSValues.set_size( tVecLocSize, 1 );

    // Get list containing the number of owned adofs of each processor
    Matrix< DDUMat > tNumOwnedList;
    comm_gather_and_broadcast( tVecLocSize, tNumOwnedList );

    Matrix< DDUMat > tOwnedOffsetList( tNumOwnedList.length(), 1, 0 );

    // Loop over all entries to create the offsets. Starting with 1
    for ( uint Ij = 1; Ij < tOwnedOffsetList.length(); Ij++ )
    {
        // Add the number of owned adofs of the previous processor to the offset of the previous processor
        tOwnedOffsetList( Ij, 0 ) = tOwnedOffsetList( Ij - 1, 0 ) + tNumOwnedList( Ij - 1, 0 );
    }

    for ( sint Ik = 0; Ik < tVecLocSize; Ik++ )
    {
        tVal( Ik, 0 ) = tOwnedOffsetList( par_rank(), 0 ) + Ik;
    }

    VecGetValues( mPetscVector, tVecLocSize, tVal.data(), LHSValues.data() );
}

//-----------------------------------------------------------------------------

void Vector_PETSc::extract_copy( Vector< real >& aVector )
{
    sint tVecLocSize;
    VecGetLocalSize( mPetscVector, &tVecLocSize );

    // FIXME replace with VecGetArray()
    Matrix< DDSMat > tVal( tVecLocSize, 1, 0 );
    aVector.resize( tVecLocSize );

    // Get list containing the number of owned adofs of each processor
    Matrix< DDUMat > tNumOwnedList;
    comm_gather_and_broadcast( tVecLocSize, tNumOwnedList );

    Matrix< DDUMat > tOwnedOffsetList( tNumOwnedList.length(), 1, 0 );

    // Loop over all entries to create the offsets. Starting with 1
    for ( uint Ij = 1; Ij < tOwnedOffsetList.length(); Ij++ )
    {
        // Add the number of owned adofs of the previous processor to the offset of the previous processor
        tOwnedOffsetList( Ij, 0 ) = tOwnedOffsetList( Ij - 1, 0 ) + tNumOwnedList( Ij - 1, 0 );
    }

    for ( sint Ik = 0; Ik < tVecLocSize; Ik++ )
    {
        tVal( Ik, 0 ) = tOwnedOffsetList( par_rank(), 0 ) + Ik;
    }

    VecGetValues( mPetscVector, tVecLocSize, tVal.data(), aVector.memptr() );
}

//-----------------------------------------------------------------------------

void
Vector_PETSc::import_local_to_global( sol::Dist_Vector& aSourceVec )
{
    // cast source vector to vector_petsc
    Vector_PETSc& tPetscSourceVec = dynamic_cast< Vector_PETSc& >( aSourceVec );

    // get raw petsc vector of source vector
    Vec tSourceVec = tPetscSourceVec.get_petsc_vector();

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
                "Vector_PETSc::import_local_to_global - target and source vectors are both full vectors: case not implemented." );

        // get number of source dofs
        PetscInt tNumSourceIds;
        ISGetLocalSize( tPetscSourceIds, &tNumSourceIds );

        // get list of source petsc IDs
        const PetscInt* tSourceIdList;
        ISGetIndices( tPetscSourceIds, &tSourceIdList );

        // get array from source petsc vector
        PetscScalar* tSourceValues;
        VecGetArray( tSourceVec, &tSourceValues );

        // set values in target vector
        VecSetValues( mPetscVector, tNumSourceIds, tSourceIdList, tSourceValues, INSERT_VALUES );

        // free memory allocated in petsc calls
        ISRestoreIndices( tPetscSourceIds, &tSourceIdList );
        VecRestoreArray( tSourceVec, &tSourceValues );

        // assemble target vector
        VecAssemblyBegin( mPetscVector );
        VecAssemblyEnd( mPetscVector );
    }
    else
    {
        // create scatter object
        VecScatter tVecScatter;
        VecScatterCreate( tSourceVec, tPetscTargetIds, mPetscVector, NULL, &tVecScatter );

        // perform scattering
        VecScatterBegin( tVecScatter, tSourceVec, mPetscVector, INSERT_VALUES, SCATTER_FORWARD );
        VecScatterEnd( tVecScatter, tSourceVec, mPetscVector, INSERT_VALUES, SCATTER_FORWARD );

        // destroy the scatter vector object
        VecScatterDestroy(&tVecScatter);
    }
}

//-----------------------------------------------------------------------------
// this is only used if the local vector is full vector and source vector is vector of only owned dofs

void
Vector_PETSc::import_local_to_global( Vec aSourceVec )
{
    // get list of target petsc ids of from map
    IS tPetscTargetIds = mMap->get_petsc_ids();

    // create scatter object
    VecScatter tVecScatter;
    VecScatterCreate( aSourceVec, tPetscTargetIds, mPetscVector, NULL, &tVecScatter );

    // perform scattering
    VecScatterBegin( tVecScatter, aSourceVec, mPetscVector, INSERT_VALUES, SCATTER_FORWARD );
    VecScatterEnd( tVecScatter, aSourceVec, mPetscVector, INSERT_VALUES, SCATTER_FORWARD );

    // destroy the scatter vector object
    VecScatterDestroy( &tVecScatter );
}

//-----------------------------------------------------------------------------

void
Vector_PETSc::extract_my_values(
        const uint&               aNumIndices,
        const Matrix< DDSMat >&   aGlobalBlockRows,
        const uint&               aBlockRowOffsets,
        Vector< Matrix< DDRMat > >& ExtractedValues )
{
    // check that vector index  is zero
    MORIS_ASSERT( aBlockRowOffsets == 0,
            "Vector_PETSc::extract_my_values - petsc not implemented yet for aBlockRowOffsets neq 0" );

    // check that aNumIndices equals size of aGlobalBlockRows
    MORIS_ASSERT( aNumIndices == aGlobalBlockRows.numel(),
            "Vector_PETSc::extract_my_values - number of indices does not match size of aGlobalBlockRows" );

    // get map from moris id to indices in vector
    Matrix< DDSMat > tIndices = mMap->map_from_moris_ids_to_indices( aGlobalBlockRows );

    // check that aNumIndices equals size of tIndices
    MORIS_ASSERT( aNumIndices == tIndices.numel(),
            "Vector_PETSc::extract_my_values - number of indices does not match size of tIndices" );


    // allocate memory for extracted values
    ExtractedValues.resize( 1 );
    ExtractedValues( 0 ).set_size( aNumIndices, 1 );

    // get values from petsc vector
    VecGetValues( mPetscVector, tIndices.numel(), tIndices.data(), ExtractedValues( 0 ).data() );
}

//-----------------------------------------------------------------------------

void
Vector_PETSc::print() const
{
    VecView( mPetscVector, PETSC_VIEWER_STDOUT_WORLD );
}

//-----------------------------------------------------------------------------

void
Vector_PETSc::save_vector_to_HDF5( const char* aFilename )
{
    PetscViewer tViewer;

    PetscViewerHDF5Open( PETSC_COMM_WORLD, aFilename, FILE_MODE_WRITE, &tViewer );

    PetscObjectSetName( (PetscObject)mPetscVector, "Res_Vec" );

    VecView( mPetscVector, tViewer );

    PetscViewerDestroy( &tViewer );
}

//-----------------------------------------------------------------------------

void
Vector_PETSc::read_vector_from_HDF5(
        const char* aFilename,
        std::string aGroupName,
        sint        aVectorindex )
{
    PetscViewer tViewer;

    PetscViewerHDF5Open( PETSC_COMM_WORLD, aFilename, FILE_MODE_READ, &tViewer );

    PetscObjectSetName( (PetscObject)mPetscVector, "Res_Vec" );

    VecLoad( mPetscVector, tViewer );

    PetscViewerDestroy( &tViewer );
}

//-----------------------------------------------------------------------------

real*
Vector_PETSc::get_values_pointer()
{
    MORIS_ERROR( false, "get_values_pointer() not implemented yet for a PETSc distributed vector." );
    return nullptr;
}

// ----------------------------------------------------------------------------

void
Vector_PETSc::save_vector_to_matlab_file( const char* aFilename )
{
    PetscViewer tViewer;

    PetscViewerCreate( PETSC_COMM_WORLD, &tViewer );
    PetscViewerSetType( tViewer, PETSCVIEWERASCII );
    PetscViewerPushFormat( tViewer, PETSC_VIEWER_ASCII_MATLAB );
    PetscViewerFileSetName( tViewer, aFilename );

    VecView( mPetscVector, tViewer );

    PetscViewerDestroy( &tViewer );
}
