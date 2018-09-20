/*
 * MapPETSc.cpp
 *
 *  Created on: Jan 10, 2018
 *      Author: schmidt
 */
#include "cl_Map_PETSc.hpp"

extern moris::Comm_Manager gMorisComm;
using namespace moris;

Map_PETSc::Map_PETSc(const moris::uint        & aNumMyDofs,
                     const Matrix< DDSMat >    & aMyLocaltoGlobalMap,
                     const Matrix< DDUMat > & aMyConstraintDofs) : moris::Map_Class()
{
    AODestroy( &mPETScMap );
    //size_t rank = par_rank();
    //MPI_Comm_split(gMorisComm.get_global_comm(), rank%1, 0, &PETSC_COMM_WORLD);

    // Get necessary inputs for Epetra Maps
    moris::uint tNumMyDofs        =  aNumMyDofs;
    moris::uint tNumGlobalDofs    =  aNumMyDofs;

    // sum up all distributed dofs
    Sum_All_Local_Int( tNumMyDofs, tNumGlobalDofs );

    // vector constraint dofs
    moris::Matrix< DDSMat > tMyGlobalConstraintDofs;

    this->translator( aNumMyDofs, tNumGlobalDofs,  aMyLocaltoGlobalMap, tMyGlobalConstraintDofs, aMyConstraintDofs );

    // Build PETSc AO map
    AOCreateBasic( PETSC_COMM_WORLD, tNumMyDofs, aMyLocaltoGlobalMap.data(), PETSC_NULL, &mPETScMap );              //PETSC_NULL for natural ordeing

    //AOView(mPETScMap,PETSC_VIEWER_STDOUT_WORLD);
}

// ----------------------------------------------------------------------------
Map_PETSc::~Map_PETSc()
{
}

// ----------------------------------------------------------------------------
void Map_PETSc::translator(const moris::uint        & aNumMyDofs,
                           const moris::uint        & aNumGlobalDofs,
                           const Matrix< DDSMat >  & aMyLocaltoGlobalMap,
                                 Matrix< DDSMat >  & aMyGlobalConstraintDofs,
                           const Matrix< DDUMat > & aMyConstraintDofs)
{
    // Set size of vector local constraint dofs
    aMyGlobalConstraintDofs.set_size( aNumMyDofs - aMyConstraintDofs.n_rows(), 1 );

    // Set Bitset with global dof size
    moris::BoostBitset tBitset( aNumGlobalDofs );

    // Set bitset entry to true if dof is constrained
    for ( moris::uint Ik=0; Ik< aMyConstraintDofs.n_rows(); Ik++ )
       {
           tBitset.set( aMyConstraintDofs( Ik,0) );
       }

    // if bitset entry = false, than add my global topology to aMyGlobalConstraintDofs
    moris::uint tCount = 0;
    for ( moris::uint k=0; k<aNumMyDofs; ++k )
    {
        if (!tBitset.test( aMyLocaltoGlobalMap(k,0) ))
        {
            aMyGlobalConstraintDofs( tCount++, 0 ) = aMyLocaltoGlobalMap( k, 0 );
        }
    }
}

