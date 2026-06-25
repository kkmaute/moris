/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SOL_Dist_Map_Custom.cpp
 *
 */

#include "cl_SOL_Dist_Map_Custom.hpp"

namespace moris
{
    Dist_Map_Custom::Dist_Map_Custom( const Matrix< DDSMat >& aMyGlobalOwnedIds, const Matrix< DDSMat >& aMyGlobalOwnedAndSharedIds )
    {
        // Create a map from the application ids to the petsc ids
        // to this we transfer everything to all processors ( this is similiar to AO object implemenation in petsc)

        // STEP 1: Get the size on each processor
        Vector< sint > tRecvCounts( par_size() );
        sint           tSendCount = aMyGlobalOwnedIds.numel();
        MPI_Allgather( &tSendCount, 1, MPI_INT, tRecvCounts.memptr(), 1, MPI_INT, MPI_COMM_WORLD );

        // STEP 2: Create the displacement vector indicating offeset on each processor
        Vector< sint > tDispls( par_size(), 0 );
        for ( moris_id iProc = 1; iProc < par_size(); iProc++ )
        {
            tDispls( iProc ) = tDispls( iProc - 1 ) + tRecvCounts( iProc - 1 );
        }

        // STEP 3: Prpeare the receive buffer
        int            tTotalCount = tDispls( par_size() - 1 ) + tRecvCounts( par_size() - 1 );
        Vector< sint > tRecvBuf( tTotalCount, 1 );

        // STEP 4: Allgather the data
        MPI_Allgatherv( aMyGlobalOwnedIds.data(), tSendCount, MPI_INT, tRecvBuf.memptr(), tRecvCounts.memptr(), tDispls.memptr(), MPI_INT, MPI_COMM_WORLD );

        // now put the data in a map
        for ( sint i = 0; i < tTotalCount; i++ )
        {
            mApplicationToPetsc[ tRecvBuf( i ) ] = i;
        }

        mMorisIDsOwned          = aMyGlobalOwnedIds;
        mMorisIDsOwnedAndShared = aMyGlobalOwnedAndSharedIds;
    }

    //---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    //---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Dist_Map_Custom::Dist_Map_Custom( const Matrix< DDSMat >& aMyGlobalOwnedIds, const Matrix< DDSMat >& aMyGlobalOwnedAndSharedIds, bool aIsFullMap )
    {
        MORIS_ASSERT( aIsFullMap, "Map_PETSc::Map_PETSc - full map can only be constructed with aIsFullMap = true" );
        mIsFullMap = aIsFullMap;

        // Create a map from the application ids to the petsc ids
        // to this we transfer everything to all processors ( this is similiar to AO object implemenation in petsc)

        // STEP 1: Get the size on each processor
        Vector< sint > tRecvCounts( par_size() );
        sint           tSendCount = aMyGlobalOwnedIds.numel();
        MPI_Allgather( &tSendCount, 1, MPI_INT, tRecvCounts.memptr(), 1, MPI_INT, MPI_COMM_WORLD );

        // STEP 2: Create the displacement vector indicating offeset on each processor
        Vector< sint > tDispls( par_size() , 0);
        for ( moris_id iProc = 1; iProc < par_size(); iProc++ )
        {
            tDispls( iProc ) = tDispls( iProc - 1 ) + tRecvCounts( iProc - 1 );
        }

        // STEP 3: Prpeare the receive buffer
        int            tTotalCount = tDispls( par_size() - 1 ) + tRecvCounts( par_size() - 1 );
        Vector< sint > tRecvBuf( tTotalCount, 1 );

        // STEP 4: Allgather the data
        MPI_Allgatherv( aMyGlobalOwnedIds.data(), tSendCount, MPI_INT, tRecvBuf.memptr(), tRecvCounts.memptr(), tDispls.memptr(), MPI_INT, MPI_COMM_WORLD );

        // Create a map from the application ids to the petsc ids
        // to this we transfer everything to all processors ( this is similiar to AO object implemenation in petsc)
        // now put the data in a map
        for ( sint i = 0; i < tTotalCount; i++ )
        {
            mApplicationToPetsc[ tRecvBuf( i ) ] = i;
        }

        mMorisIDsOwned          = aMyGlobalOwnedIds;
        mMorisIDsOwnedAndShared = aMyGlobalOwnedAndSharedIds;

        // create a local to global mapping petsc object
        ISLocalToGlobalMappingCreate(
                PETSC_COMM_WORLD,
                1,
                aMyGlobalOwnedAndSharedIds.n_rows(),
                aMyGlobalOwnedAndSharedIds.data(),
                PETSC_COPY_VALUES,
                &mMorisIDtoIndexMap );
    }

    //---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Dist_Map_Custom::Dist_Map_Custom( const Matrix< DDSMat >& aMyGlobalOwnedIds, const Matrix< DDUMat >& aConstrainedDOFs )
    {
        // Create a map from the application ids to the petsc ids
        // to this we transfer everything to all processors ( this is similiar to AO object implemenation in petsc)

        // STEP 1: Get the size on each processor
        Vector< sint > tRecvCounts( par_size() );
        sint           tSendCount = aMyGlobalOwnedIds.numel();
        MPI_Allgather( &tSendCount, 1, MPI_INT, tRecvCounts.memptr(), 1, MPI_INT, MPI_COMM_WORLD );

        // STEP 2: Create the displacement vector indicating offeset on each processor
        Vector< sint > tDispls( par_size(),  0);
        for ( moris_id iProc = 1; iProc < par_size(); iProc++ )
        {
            tDispls( iProc ) = tDispls( iProc - 1 ) + tRecvCounts( iProc - 1 );
        }

        // STEP 3: Prpeare the receive buffer
        int            tTotalCount = tDispls( par_size() - 1 ) + tRecvCounts( par_size() - 1 );
        Vector< sint > tRecvBuf( tTotalCount, 1 );

        // STEP 4: Allgather the data
        MPI_Allgatherv( aMyGlobalOwnedIds.data(), tSendCount, MPI_INT, tRecvBuf.memptr(), tRecvCounts.memptr(), tDispls.memptr(), MPI_INT, MPI_COMM_WORLD );

        // now put the data in a map
        for ( sint i = 0; i < tTotalCount; i++ )
        {
            mApplicationToPetsc[ tRecvBuf( i ) ] = i;
        }

        mMorisIDsOwned = aMyGlobalOwnedIds;
        // mMorisIDsOwnedAndShared = aConstrainedDOFs;
    }

    //---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void
    Dist_Map_Custom::map_from_moris_ids_to_petsc_ids( Matrix< DDSMat >& aGlobalIds )
    {
        // Assuming Matrix<DDSMat> can be iterated and modified directly
        for ( auto& id : aGlobalIds )
        {
            if ( id != -1 )
            {
                id = mApplicationToPetsc[ id ];
            }
        }
    }

    //---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        void
    Dist_Map_Custom::map_from_moris_ids_to_petsc_ids( Matrix< DDUMat >& aGlobalIds )
    {
        // Assuming Matrix<DDSMat> can be iterated and modified directly
        for ( auto& id : aGlobalIds )
        {
            id = (uint)mApplicationToPetsc[ id ];
        }
    }

    Matrix< DDSMat >
    Dist_Map_Custom::map_from_moris_ids_to_indices( const Matrix< DDSMat >& aGlobalIds )
    {
        // get number of IDs
        uint tNumIds = aGlobalIds.numel();

        // allocate vector for indices
        Matrix< DDSMat > tIndices( tNumIds, 1 );

        // determine indices for given IDs, if give ID does not exists a "-1" is inserted
        ISGlobalToLocalMappingApply( mMorisIDtoIndexMap, IS_GTOLM_MASK, tNumIds, aGlobalIds.data(), nullptr, tIndices.data() );

        return tIndices;
    }

    //---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Matrix< DDSMat >
    Dist_Map_Custom::map_from_moris_ids_to_indices( const Vector< sint >& aGlobalIds )
    {
        // get number of IDs
        uint tNumIds = aGlobalIds.size();

        // allocate vector for indices
        Matrix< DDSMat > tIndices( tNumIds, 1 );

        // determine indices for given IDs, if give ID does not exists a "-1" is inserted
        ISGlobalToLocalMappingApply( mMorisIDtoIndexMap, IS_GTOLM_MASK, tNumIds, aGlobalIds.memptr(), nullptr, tIndices.data() );

        return tIndices;
    }

}    // namespace moris
