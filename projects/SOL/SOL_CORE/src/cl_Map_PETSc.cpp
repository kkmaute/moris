/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Map_PETSc.cpp
 *
 */

#include "cl_Map_PETSc.hpp"

extern moris::Comm_Manager gMorisComm;
using namespace moris;

Map_PETSc::Map_PETSc(
        const Matrix< DDSMat >& aMyGlobalOwnedIds,
        const Matrix< DDUMat >& aMyConstraintDofs )
        : moris::sol::Dist_Map()
{
    // clear member variables
    AODestroy( &mPETScMap );
    ISLocalToGlobalMappingDestroy( &mMorisIDtoIndexMap );
    ISDestroy( &mPETScIDs );

    // Get number of owned DOFs
    moris::uint tNumMyOwnedDofs = aMyGlobalOwnedIds.n_rows();

    // Check that IDs are numbered consecutively from 0 to max(tNumMyDofs)
    MORIS_ASSERT( (sint)( sum_all( tNumMyOwnedDofs ) - 1 ) == (sint)max_all( aMyGlobalOwnedIds.max() ),
            "Map_PETSc::Map_PETSc - IDs are not numbered consecutively starting from 0" );

    // Build PETSc AO map between moris IDs and petsc IDs based on owned free dofs
    AOCreateBasic(
            PETSC_COMM_WORLD,
            tNumMyOwnedDofs,
            aMyGlobalOwnedIds.data(),
            PETSC_NULL,
            &mPETScMap );

    // build map between moris IDs and indices of owned (free and constraint) vector
    ISLocalToGlobalMappingCreate(
            PETSC_COMM_WORLD,
            1,
            tNumMyOwnedDofs,
            aMyGlobalOwnedIds.data(),
            PETSC_COPY_VALUES,
            &mMorisIDtoIndexMap );

    ISLocalToGlobalMappingSetFromOptions( mMorisIDtoIndexMap );

    // build index list of petsc IDs of owned vector
    Matrix< DDSMat > tPetscIDs = aMyGlobalOwnedIds;

    // get petsc IDs for owned free moris IDs
    AOApplicationToPetsc(
            mPETScMap,
            tNumMyOwnedDofs,
            tPetscIDs.data() );

    // store petsc IDs in index set
    ISCreateGeneral(
            PETSC_COMM_SELF,
            tNumMyOwnedDofs,
            tPetscIDs.data(),
            PETSC_COPY_VALUES,
            &mPETScIDs );
}

// ----------------------------------------------------------------------------

Map_PETSc::Map_PETSc( const Matrix< DDSMat >& aMyGlobalOwnedIds )
        : moris::sol::Dist_Map()
{
    // clear member variables
    AODestroy( &mPETScMap );
    ISLocalToGlobalMappingDestroy( &mMorisIDtoIndexMap );
    ISDestroy( &mPETScIDs );

    // Get number of owned DOFs
    moris::uint tNumMyOwnedDofs = aMyGlobalOwnedIds.n_rows();

    // Check that IDs are numbered consecutively from 0 to max(tNumMyDofs)
    MORIS_ASSERT( (sint)( sum_all( tNumMyOwnedDofs ) - 1 ) == (sint)max_all( aMyGlobalOwnedIds.max() ),
            "Map_PETSc::Map_PETSc - IDs are not numbered consecutively starting from 0" );

    // Build PETSc AO map between moris IDs and petsc IDs based on owned dofs
    AOCreateBasic(
            PETSC_COMM_WORLD,
            tNumMyOwnedDofs,
            aMyGlobalOwnedIds.data(),
            PETSC_NULL,
            &mPETScMap );

    // build map between moris IDs and indices of owned vector
    ISLocalToGlobalMappingCreate(
            PETSC_COMM_WORLD,
            1,
            tNumMyOwnedDofs,
            aMyGlobalOwnedIds.data(),
            PETSC_COPY_VALUES,
            &mMorisIDtoIndexMap );

    ISLocalToGlobalMappingSetFromOptions( mMorisIDtoIndexMap );

    // build index list of petsc IDs of owned vector
    Matrix< DDSMat > tPetscIDs = aMyGlobalOwnedIds;

    // get petsc IDs for owned moris IDs
    AOApplicationToPetsc(
            mPETScMap,
            tNumMyOwnedDofs,
            tPetscIDs.data() );

    // store petsc IDs in index set
    ISCreateGeneral(
            PETSC_COMM_SELF,
            tNumMyOwnedDofs,
            tPetscIDs.data(),
            PETSC_COPY_VALUES,
            &mPETScIDs );
}

// ----------------------------------------------------------------------------

Map_PETSc::Map_PETSc(
        const Matrix< DDSMat >& aMyGlobalOwnedIds,
        const Matrix< DDSMat >& aMyGlobalOwnedAndSharedIds,
        bool                    aIsFullMap )
{
    // clear member variables
    AODestroy( &mPETScMap );
    ISLocalToGlobalMappingDestroy( &mMorisIDtoIndexMap );
    ISDestroy( &mPETScIDs );

    // set full map flag
    MORIS_ASSERT( aIsFullMap, "Map_PETSc::Map_PETSc - full map can only be constructed with aIsFullMap = true" );
    mIsFullMap = aIsFullMap;

    // Get number of owned DOFs
    moris::uint tNumMyOwnedDofs = aMyGlobalOwnedIds.n_rows();

    // Check that IDs are numbered consecutively from 0 to max(tNumMyDofs)
    MORIS_ASSERT( (sint)( sum_all( tNumMyOwnedDofs ) - 1 ) == (sint)max_all( aMyGlobalOwnedIds.max() ),
            "Map_PETSc::Map_PETSc - IDs are not numbered consecutively starting from 0" );

    // Build PETSc AO map between moris IDs and petsc IDs based on owned dofs
    AOCreateBasic(
            PETSC_COMM_WORLD,
            tNumMyOwnedDofs,
            aMyGlobalOwnedIds.data(),
            PETSC_NULL,
            &mPETScMap );

    // Get number of owned and shared DOFs
    moris::uint tNumMyOwnedAndSharedDofs = aMyGlobalOwnedAndSharedIds.n_rows();

    // build map between moris IDs and indices of full vector
    ISLocalToGlobalMappingCreate(
            PETSC_COMM_WORLD,
            1,
            tNumMyOwnedAndSharedDofs,
            aMyGlobalOwnedAndSharedIds.data(),
            PETSC_COPY_VALUES,
            &mMorisIDtoIndexMap );

    ISLocalToGlobalMappingSetFromOptions( mMorisIDtoIndexMap );

    // build index list of petsc IDs of full vector
    Matrix< DDSMat > tPetscIDs = aMyGlobalOwnedAndSharedIds;

    // get petsc IDs for owned and shared moris IDs
    AOApplicationToPetsc(
            mPETScMap,
            tNumMyOwnedAndSharedDofs,
            tPetscIDs.data() );

    // store petsc IDs in index set
    ISCreateGeneral(
            PETSC_COMM_SELF,
            tNumMyOwnedAndSharedDofs,
            tPetscIDs.data(),
            PETSC_COPY_VALUES,
            &mPETScIDs );
}

// ----------------------------------------------------------------------------
Map_PETSc::~Map_PETSc()
{
    // clear member variables
    AODestroy( &mPETScMap );
    ISLocalToGlobalMappingDestroy( &mMorisIDtoIndexMap );
    ISDestroy( &mPETScIDs );
}

// ----------------------------------------------------------------------------

void
Map_PETSc::translator(
        const moris::uint&      aMaxDofsId,
        const moris::uint&      aNumMyDofs,
        const Matrix< DDSMat >& aMyLocaltoGlobalMap,
        Matrix< DDSMat >&       aMyFreeDofs,
        const Matrix< DDUMat >& aMyConstraintDofs )
{
    // Set size of vector local constraint dofs
    aMyFreeDofs.set_size( aNumMyDofs, 1 );

    // Set Bitset with global dof size
    moris::BoostBitset tBitset( aMaxDofsId + 1 );

    // Set bitset entry to true if dof is constrained
    for ( moris::uint Ik = 0; Ik < aMyConstraintDofs.n_rows(); Ik++ )
    {
        tBitset.set( aMyConstraintDofs( Ik, 0 ) );
    }

    // if bitset entry = false, than add my global dofs to aMyFreeDofs vector
    moris::uint tCount = 0;
    for ( moris::uint k = 0; k < aNumMyDofs; ++k )
    {
        if ( !tBitset.test( aMyLocaltoGlobalMap( k, 0 ) ) )
        {
            aMyFreeDofs( tCount++, 0 ) = aMyLocaltoGlobalMap( k, 0 );
        }
    }

    aMyFreeDofs.resize( tCount, 1 );
}

// ----------------------------------------------------------------------------

Matrix< DDSMat >
Map_PETSc::map_from_moris_ids_to_indices( const Matrix< DDSMat >& aGlobalIds )
{
    // get number of IDs
    uint tNumIds = aGlobalIds.numel();

    // allocate vector for indices
    Matrix< DDSMat > tIndices( tNumIds, 1 );

    // determine indices for given IDs, if give ID does not exists a "-1" is inserted
    ISGlobalToLocalMappingApply( mMorisIDtoIndexMap, IS_GTOLM_MASK, tNumIds, aGlobalIds.data(), NULL, tIndices.data() );

    return tIndices;
}
