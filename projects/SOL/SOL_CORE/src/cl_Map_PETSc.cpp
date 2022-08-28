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
        const Matrix< DDSMat > & aMyGlobalIds,
        const Matrix< DDUMat > & aMyConstraintDofs)
: moris::sol::Dist_Map()
{
    AODestroy( &mPETScMap );
    //size_t rank = par_rank();
    //MPI_Comm_split(gMorisComm.get_global_comm(), rank%1, 0, &PETSC_COMM_WORLD);

    // Get necessary inputs for maps
    moris::uint tNumMyDofs =  aMyGlobalIds.numel();
    moris::uint tMaxDofId  =  aMyGlobalIds.max();

    // vector constraint dofs
    moris::Matrix< DDSMat > tMyGlobalConstraintDofs;

    this->translator( tMaxDofId, tNumMyDofs, aMyGlobalIds, tMyGlobalConstraintDofs, aMyConstraintDofs );

    // Build PETSc AO map
    AOCreateBasic( PETSC_COMM_WORLD, tNumMyDofs, aMyGlobalIds.data(), PETSC_NULL, &mPETScMap );              //PETSC_NULL for natural ordeing

    //AOView(mPETScMap,PETSC_VIEWER_STDOUT_WORLD);
}

// ----------------------------------------------------------------------------

Map_PETSc::Map_PETSc( const Matrix< DDSMat > & aMyGlobalIds )
: moris::sol::Dist_Map()
{
    AODestroy( &mPETScMap );
    //size_t rank = par_rank();
    //MPI_Comm_split(gMorisComm.get_global_comm(), rank%1, 0, &PETSC_COMM_WORLD);

    // Get necessary inputs for Epetra Maps
    moris::uint tNumMyDofs =  aMyGlobalIds.n_rows();

    // Build PETSc AO map
    AOCreateBasic( PETSC_COMM_WORLD, tNumMyDofs, aMyGlobalIds.data(), PETSC_NULL, &mPETScMap );              //PETSC_NULL for natural ordeing

    //AOView(mPETScMap,PETSC_VIEWER_STDOUT_WORLD);
}

// ----------------------------------------------------------------------------
Map_PETSc::~Map_PETSc()
{
    AODestroy( &mPETScMap );
}

// ----------------------------------------------------------------------------

void Map_PETSc::translator(
        const moris::uint       & aMaxDofsId,
        const moris::uint       & aNumMyDofs,
        const Matrix< DDSMat >  & aMyLocaltoGlobalMap,
        Matrix< DDSMat >        & aMyGlobalConstraintDofs,
        const Matrix< DDUMat >  & aMyConstraintDofs )
{
    // Set size of vector local constraint dofs
    aMyGlobalConstraintDofs.set_size( aNumMyDofs, 1 );

    // Set Bitset with global dof size
    moris::BoostBitset tBitset( aMaxDofsId + 1 );

    // Set bitset entry to true if dof is constrained
    for ( moris::uint Ik=0; Ik< aMyConstraintDofs.n_rows(); Ik++ )
    {
        tBitset.set( aMyConstraintDofs( Ik, 0 ) );
    }

    // if bitset entry = false, than add my global topology to aMyGlobalConstraintDofs
    moris::uint tCount = 0;
    for ( moris::uint k = 0; k < aNumMyDofs; ++k )
    {
        if ( !tBitset.test( aMyLocaltoGlobalMap( k, 0 ) ) )
        {
            aMyGlobalConstraintDofs( tCount++, 0 ) = aMyLocaltoGlobalMap( k, 0 );
        }
    }

    aMyGlobalConstraintDofs.resize( tCount, 1 );
}

