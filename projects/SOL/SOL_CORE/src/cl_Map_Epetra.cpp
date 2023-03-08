/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Map_Epetra.cpp
 *
 */

#include "cl_Map_Epetra.hpp"
#include "cl_Communication_Tools.hpp"    // COM/src
#include "cl_DLA_Solver_Interface.hpp"

#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_FEVector.h"
#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_Import.h"

extern moris::Comm_Manager gMorisComm;

using namespace moris;

// ----------------------------------------------------------------------------------------------------------------------

Map_Epetra::Map_Epetra(
        const Matrix< DDSMat >& aMyGlobalIds,
        const Matrix< DDUMat >& aMyConstraintDofs )
        : sol::Dist_Map()
{
    // Minimum index value used for arrays that use this map. Typically 0 for C/C++ and 1 for Fortran.
    moris::uint tIndexBase = 0;

    // Get necessary inputs for Epetra Maps
    moris::uint tNumMyDofs = aMyGlobalIds.numel();
    moris::uint tMaxDofId  = aMyGlobalIds.max();

    // vector constraint dofs
    Matrix< DDSMat > tMyGlobalConstraintDofs;

    this->translator( tMaxDofId, tNumMyDofs, aMyGlobalIds, tMyGlobalConstraintDofs, aMyConstraintDofs );

    // build maps
    mEpetraMap = new Epetra_Map(
            -1,
            tMyGlobalConstraintDofs.n_rows(),
            tMyGlobalConstraintDofs.data(),
            tIndexBase,
            *mEpetraComm.get_epetra_comm() );

    this->build_point_map();
}

// ----------------------------------------------------------------------------------------------------------------------

Map_Epetra::Map_Epetra( const Matrix< DDSMat >& aMyGlobalIds )
        : sol::Dist_Map()
{
    // Minimum index value used for arrays that use this map. Typically 0 for C/C++ and 1 for Fortran.
    moris::uint tIndexBase = 0;

    // build maps
    mEpetraMap = new Epetra_Map(
            -1,
            aMyGlobalIds.n_rows(),
            aMyGlobalIds.data(),
            tIndexBase,
            *mEpetraComm.get_epetra_comm() );

    this->build_point_map();
}

// ----------------------------------------------------------------------------------------------------------------------

Map_Epetra::~Map_Epetra()
{
    delete mEpetraMap;
    delete mEpetraPointMap;
    delete mFullOverlappingMap;

    delete mFullToFreePoint;
}

// ----------------------------------------------------------------------------------------------------------------------

void
Map_Epetra::build_point_map()
{
    // build point map
    sint tNumMyPoints = mEpetraMap->NumMyPoints();

    moris::uint tGlobID = get_processor_offset( tNumMyPoints );

    //    int* myDofIds = (int*) alloca (sizeof(int)*numMyPoints);
    moris::Matrix< IdMat > tMyDofIds( tNumMyPoints, 1 );

    for ( sint Ik = 0; Ik < tNumMyPoints; Ik++ )
    {
        tMyDofIds( Ik ) = tGlobID++;
    }

    mEpetraPointMap = new Epetra_Map(
            -1,
            tNumMyPoints,
            tMyDofIds.data(),
            0,
            *mEpetraComm.get_epetra_comm() );
}

// ----------------------------------------------------------------------------------------------------------------------

void
Map_Epetra::translator(
        const moris::uint&      aMaxDofsId,
        const moris::uint&      aNumMyDofs,
        const Matrix< DDSMat > & aMyLocalToGlobalMap,
        Matrix< DDSMat >&       aMyGlobalConstraintDofs,
        const Matrix< DDUMat >& aMyConstraintDofs )
{
    // Set size of vector local constraint dofs
    aMyGlobalConstraintDofs.set_size( aNumMyDofs, 1 );

    // Set Bitset with global dof size
    moris::BoostBitset tBitset( aMaxDofsId + 1 );

    // Set bitset entry to true if dof is constrained
    for ( moris::uint Ik = 0; Ik < aMyConstraintDofs.numel(); Ik++ )
    {
        tBitset.set( aMyConstraintDofs( Ik ) );
    }

    // if bitset entry = false, than add my global topology to aMyGlobalConstraintDofs
    moris::uint tCount = 0;
    for ( moris::uint k = 0; k < aNumMyDofs; ++k )
    {
        if ( !tBitset.test( aMyLocalToGlobalMap( k ) ) )
        {
            aMyGlobalConstraintDofs( tCount++ ) = aMyLocalToGlobalMap( k );
        }
    }

    aMyGlobalConstraintDofs.resize( tCount, 1 );
}

// ----------------------------------------------------------------------------------------------------------------------

void
Map_Epetra::build_dof_translator(
        const Matrix< IdMat >& aFullMap,
        const bool             aFlag )
{
    mFullOverlappingMap = new Epetra_Map(
            -1,
            aFullMap.n_rows(),
            aFullMap.data(),
            0,
            *mEpetraComm.get_epetra_comm() );

    mFullToFreePoint = new Epetra_MultiVector( *mFullOverlappingMap, 1 );

    // Initialize every point as constrained
    mFullToFreePoint->PutScalar( -1 );

    //    Epetra_BlockMap* masterMap     = aModel->GetDofHandler()->GetMasterDofMap()->GetEpetraMap();
    Epetra_MultiVector* tTempVec = new Epetra_MultiVector( *mEpetraMap, 1 );

    // Initialize every point as constrained
    tTempVec->PutScalar( -1 );

    sint* tPointToElementList = mEpetraMap->PointToElementList();

    // Get the number of free DoFs
    sint tNumFreeDoFs = mEpetraMap->NumMyPoints();

    // Make sure the number of free DoFs is consistent
    MORIS_ERROR( tNumFreeDoFs == mEpetraPointMap->NumMyPoints(), "Map_Epetra::build_dof_translator(), map size must be the same." );

    // Loop over all free DoFs
    for ( sint Ik = 0; Ik < tNumFreeDoFs; Ik++ )
    {
        sint GlobalFreeDofID = mEpetraMap->GID( tPointToElementList[ Ik ] );

        tTempVec->ReplaceGlobalValue(
                GlobalFreeDofID,
                0,
                (moris::real)mEpetraPointMap->GID( Ik ) );
    }

    Epetra_Import* tImporter = new Epetra_Import( *mFullOverlappingMap, *mEpetraMap );

    // Update mFullToFreePoint by importing the local master masterTemp into the global mFullToFreePoint
    sint tStatus = mFullToFreePoint->Import( *tTempVec, *tImporter, Insert );

    // std::cout<<*mFullToFreePoint<<std::endl;

    if ( tStatus != 0 )
    {
        MORIS_ERROR( false, "Status return error!\n" );
    }

    delete ( tImporter );
    delete ( tTempVec );
}

// ----------------------------------------------------------------------------------------------------------------------

void
Map_Epetra::translate_ids_to_free_point_ids(
        const moris::Matrix< IdMat >& aIdsIn,
        moris::Matrix< IdMat >&       aIdsOut,
        const bool&                   aIsBuildGraph )
{
    uint tNumIds = aIdsIn.numel();

    aIdsOut.set_size( tNumIds, 1, MORIS_ID_MAX );

    MORIS_ASSERT( mFullOverlappingMap != nullptr,
            "Map_Epetra::translate_ids_to_free_point_ids(), mFullOverlappingMap not set.\n");

    // Loop over all DoFs of the current element
    for ( uint Ik = 0; Ik < tNumIds; Ik++ )
    {
        // Get local index
        moris_index tLocalIndex = mFullOverlappingMap->LID( aIdsIn( Ik ) );

        MORIS_ASSERT( mFullToFreePoint->MyLength() > tLocalIndex,
                "Map_Epetra::translate_ids_to_free_point_ids(), tLocalIndex out of bounds.\n" );

        // FIXME temporary workaround for tLocalIndex < 0
        moris::real tIdOut;

        if ( tLocalIndex >= 0 )
        {
            tIdOut = mFullToFreePoint->Values()[ tLocalIndex ];
        }
        else
        {
            tIdOut = -1;
        }

        // FIXME: comment needed
        if ( !aIsBuildGraph and tIdOut == -1 )
        {
            aIdsOut( Ik ) = MORIS_ID_MAX;
        }
        else
        {
            aIdsOut( Ik ) = tIdOut;
        }
    }
}

// ----------------------------------------------------------------------------------------------------------------------

moris::sint
Map_Epetra::return_local_ind_of_global_Id( moris::uint aGlobalId ) const
{
    MORIS_ERROR( mEpetraMap != NULL, "Map_Epetra::return_local_ind_of_global_Id(), Map does not exist" );

    return mEpetraMap->LID( (int)aGlobalId );
}

// ----------------------------------------------------------------------------------------------------------------------

void
Map_Epetra::print()
{
    std::cout << *mEpetraMap << std::endl;
}
