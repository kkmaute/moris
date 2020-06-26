/*
 * cl_Map_Epetra.cpp
 *
 *  Created on: Apr 10, 2018
 *      Author: schmidt
 */
#include "cl_Map_Epetra.hpp"
#include "cl_Communication_Tools.hpp" // COM/src

extern moris::Comm_Manager gMorisComm;

using namespace moris;

// ----------------------------------------------------------------------------------------------------------------------

Map_Epetra::Map_Epetra(
        const Matrix< DDSMat > & aMyGlobalIds,
        const Matrix< DDUMat > & aMyConstraintDofs )
: sol::Dist_Map()
{
    delete( mEpetraMap );

    // Minimum index value used for arrays that use this map. Typically 0 for C/C++ and 1 for Fortran.
    moris::uint tIndexBase = 0;

    // Get necessary inputs for Epetra Maps
    moris::uint tNumMyDofs =  aMyGlobalIds.numel();
    moris::uint tMaxDofId  =  aMyGlobalIds.max();

    // vector constraint dofs
    Matrix< DDSMat > tMyGlobalConstraintDofs;

    this->translator( tMaxDofId, tNumMyDofs,  aMyGlobalIds, tMyGlobalConstraintDofs, aMyConstraintDofs );

    // build maps
    mEpetraMap = new Epetra_Map( -1, tMyGlobalConstraintDofs.n_rows(), tMyGlobalConstraintDofs.data() , tIndexBase, *mEpetraComm.get_epetra_comm() );
}


// ----------------------------------------------------------------------------------------------------------------------

Map_Epetra::Map_Epetra( const Matrix< DDSMat > & aMyGlobalIds )
: sol::Dist_Map()
{
    delete( mEpetraMap );

    // Minimum index value used for arrays that use this map. Typically 0 for C/C++ and 1 for Fortran.
    moris::uint tIndexBase = 0;

    // build maps
    mEpetraMap = new Epetra_Map( -1, aMyGlobalIds.n_rows(), aMyGlobalIds.data() , tIndexBase, *mEpetraComm.get_epetra_comm() );
}

// ----------------------------------------------------------------------------------------------------------------------

Map_Epetra::~Map_Epetra()
{
    delete( mEpetraMap );
}

// ----------------------------------------------------------------------------------------------------------------------
void Map_Epetra::translator(
        const moris::uint      & aMaxDofsId,
        const moris::uint      & aNumMyDofs,
        const Matrix< DDSMat > & aMyLocaltoGlobalMap,
        Matrix< DDSMat >       & aMyGlobalConstraintDofs,
        const Matrix< DDUMat > & aMyConstraintDofs )
{
    // Set size of vector local constraint dofs
    aMyGlobalConstraintDofs.set_size( aNumMyDofs, 1 );

    // Set Bitset with global dof size
    moris::BoostBitset tBitset( aMaxDofsId + 1 );

    // Set bitset entry to true if dof is constrained
    for ( moris::uint Ik=0; Ik< aMyConstraintDofs.numel(); Ik++ )
    {
        tBitset.set( aMyConstraintDofs( Ik ) );
    }

    // if bitset entry = false, than add my global topology to aMyGlobalConstraintDofs
    moris::uint tCount = 0;
    for ( moris::uint k = 0; k < aNumMyDofs; ++k )
    {
        if ( !tBitset.test( aMyLocaltoGlobalMap( k ) ) )
        {
            aMyGlobalConstraintDofs( tCount++ ) = aMyLocaltoGlobalMap( k );
        }
    }

    aMyGlobalConstraintDofs.resize( tCount, 1 );
}
// ----------------------------------------------------------------------------------------------------------------------
moris::sint Map_Epetra::return_local_ind_of_global_Id( moris::uint aGlobalId ) const
{
    MORIS_ERROR( mEpetraMap != NULL, "Map_Epetra::return_local_ind_of_global_Id(), Map does not exist");

    return mEpetraMap->LID( ( int ) aGlobalId );
}
