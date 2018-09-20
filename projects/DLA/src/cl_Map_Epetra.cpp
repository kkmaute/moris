
#include "cl_Map_Epetra.hpp"

using namespace moris;

Map_Epetra::Map_Epetra( const moris::uint      & aNumMyDofs,
                        const Matrix< DDSMat > & aMyLocaltoGlobalMap,
                        const Matrix< DDUMat > & aMyConstraintDofs,
                        const Matrix< DDSMat > & aOverlappingLocaltoGlobalMap ) :  Map_Class()
{
    delete( mFreeEpetraMap );
    delete( mFullEpetraMap );
    delete( mFullOverlappingEpetraMap );

    // Minimum index value used for arrays that use this map. Typically 0 for C/C++ and 1 for Fortran.
    moris::uint tIndexBase = 0;

    // Get necessary inputs for Epetra Maps
    moris::uint tNumMyDofs        =  aNumMyDofs;
    moris::uint tNumGlobalDofs    =  aNumMyDofs;

    // sum up all distributed dofs
#ifdef MORIS_HAVE_PARALLEL
        MPI_Allreduce(&tNumMyDofs,&tNumGlobalDofs,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#endif

    // vector constraint dofs
    Matrix< DDSMat > tMyGlobalConstraintDofs;

    this->translator(aNumMyDofs, tNumGlobalDofs,  aMyLocaltoGlobalMap, tMyGlobalConstraintDofs, aMyConstraintDofs);

    // build maps
    mFreeEpetraMap = new Epetra_Map( -1, tMyGlobalConstraintDofs.n_rows(), tMyGlobalConstraintDofs.data() , tIndexBase, *mEpetraComm.get_epetra_comm() );

    mFullEpetraMap = new Epetra_Map( -1, aMyLocaltoGlobalMap.n_rows(), aMyLocaltoGlobalMap.data() , tIndexBase, *mEpetraComm.get_epetra_comm() );

    mFullOverlappingEpetraMap = new Epetra_Map( -1, aOverlappingLocaltoGlobalMap.n_rows(), aOverlappingLocaltoGlobalMap.data() , tIndexBase, *mEpetraComm.get_epetra_comm() );
}

// ----------------------------------------------------------------------------------------------------------------------

Map_Epetra::~Map_Epetra()
{
}

// ----------------------------------------------------------------------------------------------------------------------
void Map_Epetra::translator( const moris::uint      & aNumMyDofs,
                             const moris::uint      & aNumGlobalDofs,
                             const Matrix< DDSMat > & aMyLocaltoGlobalMap,
                                   Matrix< DDSMat > & aMyGlobalConstraintDofs,
                             const Matrix< DDUMat > & aMyConstraintDofs )
{
    // Set size of vector local constraint dofs
    aMyGlobalConstraintDofs.set_size( aNumMyDofs - aMyConstraintDofs.n_rows(), 1 );

    // Set Bitset with global dof size
    moris::BoostBitset tBitset( aNumGlobalDofs );

    // Set bitset entry to true if dof is constrained
    for ( moris::uint Ik=0; Ik< aMyConstraintDofs.n_rows(); Ik++ )
       {
           tBitset.set( aMyConstraintDofs( Ik, 0 ) );
       }

    // if bitset entry = false, than add my global topology to aMyGlobalConstraintDofs
    moris::uint tCount = 0;
    for ( moris::uint k=0; k<aNumMyDofs; ++k )
    {
        if ( !tBitset.test( aMyLocaltoGlobalMap( k, 0 ) ) )
        {
            aMyGlobalConstraintDofs( tCount++, 0 ) = aMyLocaltoGlobalMap( k, 0 );
        }
    }
}
// ----------------------------------------------------------------------------------------------------------------------
const moris::sint Map_Epetra::return_local_ind_of_global_Id( moris::uint aGlobalId ) const
{
    // FIXME only work for the full maps right now
     if( mFullOverlappingEpetraMap != NULL )
     {
         return mFullOverlappingEpetraMap->LID( ( int ) aGlobalId );
     }
//    else if( mFreeEpetraMap != NULL )
//    {
//        return mFreeEpetraMap->LID( ( int ) aGlobalId );
//    }
    MORIS_ERROR( false, "Map_Epetra: no free nor full map exists");

    return -1;
}



