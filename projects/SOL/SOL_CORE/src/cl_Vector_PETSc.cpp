/*
 * cl_VectorPETSc.cpp
 *
 *  Created on: Dec 5, 2018
 *      Author: schmidt
 */
#include "cl_Vector_PETSc.hpp"

#include <petscviewerhdf5.h>

extern moris::Comm_Manager gMorisComm;

using namespace moris;

// ----------------------------------------------------------------------------

Vector_PETSc::Vector_PETSc(
        moris::Solver_Interface * aInput,
        sol::Dist_Map*            aMap,
        const sint                aNumVectors,
        bool                      aManageMap)
: moris::sol::Dist_Vector( aMap, aManageMap )
{
    //PetscScalar    tZero = 0;
    //moris::uint             aNumMyDofs          = aInput->get_num_my_dofs();
    moris::uint aNumMyDofs                      = aInput->get_my_local_global_map().n_rows();
    moris::Matrix< DDSMat > aMyLocaltoGlobalMap = aInput->get_my_local_global_map();
    moris::Matrix< DDUMat > aMyConstraintDofs   = aInput->get_constrained_Ids();

    // Get PETSc communicator
    //    PetscMPIInt                rank;
    //    MPI_Comm_rank(mComm->GetPETScComm(), &rank);
    //    MPI_Comm_split(mComm->GetPETScComm(), rank%1, 0, &PETSC_COMM_WORLD);

    // sum up all distributed dofs
    moris::uint tNumGlobalDofs = sum_all( aNumMyDofs );

    //FIXME insert boolean array for BC-- insert NumGlobalElements-- size
    mDirichletBCVec.set_size( tNumGlobalDofs, 1, 0 );

    // build BC vector
    this->dirichlet_BC_vector( mDirichletBCVec, aMyConstraintDofs );

    // Set up RHS b
    VecCreateMPI( PETSC_COMM_WORLD, aNumMyDofs, PETSC_DETERMINE, &mPetscVector );
    VecSetFromOptions( mPetscVector );
    VecSetOption(mPetscVector, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);
    VecSetUp( mPetscVector );
}

//-----------------------------------------------------------------------------

Vector_PETSc::~Vector_PETSc()
{
    VecDestroy( &mPetscVector );
}

//-----------------------------------------------------------------------------

real& Vector_PETSc::operator()( sint aGlobalId, uint aVectorIndex )
{
    MORIS_ERROR(false, "operator() not implemented for PETSc vector.");
    Matrix<DDRMat> tLHSValues;
    extract_copy( tLHSValues );
    return tLHSValues(0);
}

//-----------------------------------------------------------------------------

void Vector_PETSc::sum_into_global_values(
        const moris::Matrix< DDSMat > & aGlobalIds,
        const moris::Matrix< DDRMat > & aValues,
        const uint                    & aVectorIndex )
{
    uint tNumMyDofs = aGlobalIds.numel();
    moris::Matrix< DDSMat >tTempElemDofs ( tNumMyDofs, 1 );
    tTempElemDofs = aGlobalIds;

    //loop over elemental dofs
    for ( moris::uint Ij=0; Ij< tNumMyDofs; Ij++ )
    {
        //set constrDof to neg value
        if (mDirichletBCVec( aGlobalIds( Ij, 0 ), 0 ) == 1 )
        {
            tTempElemDofs( Ij, 0) = -1;
        }
    }

    // Apply PETSc map AO
    AOApplicationToPetsc( mMap->get_petsc_map(), tNumMyDofs, tTempElemDofs.data() );

    // Insert values into vector
    VecSetValues( mPetscVector, tNumMyDofs, tTempElemDofs.data(), aValues.data(), ADD_VALUES );
}

//-----------------------------------------------------------------------------

void Vector_PETSc::dirichlet_BC_vector(
        moris::Matrix< DDUMat >       & aDirichletBCVec,
        const moris::Matrix< DDUMat > & aMyConstraintDofs )
{
    //build vector with constraint values. unconstrained =0 constrained =1. change this to true/false
    for ( moris::uint Ik=0; Ik< aMyConstraintDofs.n_rows(); Ik++ )
    {
        aDirichletBCVec( aMyConstraintDofs( Ik, 0 ), 0 )  = 1;
    }
}

//-----------------------------------------------------------------------------

void Vector_PETSc::vector_global_assembly()
{
    VecAssemblyBegin( mPetscVector );
    VecAssemblyEnd  ( mPetscVector );

    //VecView(mPetscVector, PETSC_VIEWER_STDOUT_(PETSC_COMM_WORLD));
}

//-----------------------------------------------------------------------------

void Vector_PETSc::vec_plus_vec(
        const moris::real & aScaleA,
        sol::Dist_Vector  & aVecA,
        const moris::real & aScaleThis )
{
    PetscScalar tValueA = aScaleA;

    PetscScalar tValueThis = aScaleThis;

    VecAXPBY( mPetscVector, tValueA, tValueThis, dynamic_cast<Vector_PETSc&>(aVecA).get_petsc_vector() );
}

//-----------------------------------------------------------------------------

void Vector_PETSc::scale_vector(
        const moris::real & aValue,
        const moris::uint & aVecIndex )
{
    PetscScalar tValue = aValue;

    VecScale( mPetscVector, tValue );
}

//-----------------------------------------------------------------------------

void Vector_PETSc::vec_put_scalar( const moris::real & aValue )
{
    PetscScalar tValue = aValue;

    VecSet( mPetscVector, tValue);
}

//-----------------------------------------------------------------------------

void Vector_PETSc::random()
{
    MORIS_ERROR(false, "random(), not implemented for petsc");
}

//-----------------------------------------------------------------------------

moris::sint Vector_PETSc::vec_local_length() const
{
    moris::sint tVecLocSize = 0;
    VecGetLocalSize( mPetscVector, &tVecLocSize);
    return tVecLocSize;
}

//-----------------------------------------------------------------------------

moris::sint Vector_PETSc::vec_global_length() const
{
    moris::sint tVecSize = 0;
    VecGetSize( mPetscVector, &tVecSize );
    return  tVecSize;
}

//-----------------------------------------------------------------------------

Cell< moris::real > Vector_PETSc::vec_norm2()
{
    Cell< moris::real > tVecNorm( mNumVectors, 0.0);

    VecNorm( mPetscVector, NORM_2, tVecNorm.data().data() );
    return tVecNorm;
}

//-----------------------------------------------------------------------------

void Vector_PETSc::check_vector( )
{
    MORIS_ASSERT( false, "epetra vector should not have any input on the petsc vector" );
}

//-----------------------------------------------------------------------------

void Vector_PETSc::extract_copy( moris::Matrix< DDRMat > & LHSValues )
{
    //VecGetArray (tSolution, &  LHSValues.data());

    moris::sint tVecLocSize;
    VecGetLocalSize( mPetscVector, &tVecLocSize );

    // FIXME replace with VecGetArray()
    moris::Matrix< DDSMat > tVal ( tVecLocSize, 1 , 0 );
    LHSValues.set_size( tVecLocSize, 1 );

    // Get list containing the number of owned adofs of each processor
    Matrix< DDUMat > tNumOwnedList;
    comm_gather_and_broadcast( tVecLocSize, tNumOwnedList );

    Matrix< DDUMat > tOwnedOffsetList( tNumOwnedList.length(), 1, 0 );

    // Loop over all entries to create the offsets. Starting with 1
    for ( moris::uint Ij = 1; Ij < tOwnedOffsetList.length(); Ij++ )
    {
        // Add the number of owned adofs of the previous processor to the offset of the previous processor
        tOwnedOffsetList( Ij, 0 ) = tOwnedOffsetList( Ij-1, 0 ) + tNumOwnedList( Ij-1, 0 );
    }

    for ( moris::sint Ik=0; Ik< tVecLocSize; Ik++ )
    {
        tVal( Ik, 0 ) = tOwnedOffsetList( par_rank(), 0)+Ik;
    }

    VecGetValues( mPetscVector, tVecLocSize, tVal.data(), LHSValues.data() );
}

//-----------------------------------------------------------------------------

void Vector_PETSc::import_local_to_global( sol::Dist_Vector & aSourceVec )
{
    // FIXME change this to scatter thus that it works better in parallel
    PetscScalar tValueA = 1;
    PetscScalar tValueThis = 0;
    VecAXPBY( mPetscVector, tValueA, tValueThis, dynamic_cast<Vector_PETSc&>(aSourceVec).get_petsc_vector() );
}

//-----------------------------------------------------------------------------

void Vector_PETSc::extract_my_values(
        const moris::uint                       & aNumIndices,
        const moris::Matrix< DDSMat >           & aGlobalBlockRows,
        const moris::uint                       & aBlockRowOffsets,
        moris::Cell< moris::Matrix< DDRMat > >  & ExtractedValues )
{
    ExtractedValues.resize( 1 );
    ExtractedValues( 0 ).set_size( aNumIndices, 1 );

    VecGetValues( mPetscVector, aNumIndices, aGlobalBlockRows.data(), ExtractedValues( 0 ).data() );
}

//-----------------------------------------------------------------------------

void Vector_PETSc::print() const
{
    VecView( mPetscVector, PETSC_VIEWER_STDOUT_WORLD );
}

//-----------------------------------------------------------------------------

void Vector_PETSc::save_vector_to_HDF5( const char* aFilename )
{
    PetscViewer tViewer;

    PetscViewerHDF5Open( PETSC_COMM_WORLD, aFilename, FILE_MODE_WRITE, &tViewer );

    PetscObjectSetName( (PetscObject) mPetscVector, "Res_Vec");

    VecView( mPetscVector, tViewer );

    PetscViewerDestroy( &tViewer );
}

//-----------------------------------------------------------------------------

void Vector_PETSc::read_vector_from_HDF5( const char* aFilename )
{
    PetscViewer tViewer;

    PetscViewerHDF5Open( PETSC_COMM_WORLD, aFilename, FILE_MODE_READ, &tViewer);

    PetscObjectSetName( (PetscObject) mPetscVector, "Res_Vec");

    VecLoad( mPetscVector, tViewer );

    PetscViewerDestroy( &tViewer );
}

//-----------------------------------------------------------------------------

moris::real* Vector_PETSc::get_values_pointer()
{
    MORIS_ERROR(false, "get_values_pointer() not implemented yet for a PETSc distributed vector.");
    return nullptr;
}

// ----------------------------------------------------------------------------

void Vector_PETSc::save_vector_to_matlab_file( const char* aFilename )
{
    PetscViewer tViewer;

    PetscViewerCreate( PETSC_COMM_WORLD, &tViewer );
    PetscViewerSetType( tViewer, PETSCVIEWERASCII );
    PetscViewerPushFormat( tViewer, PETSC_VIEWER_ASCII_MATLAB );
    PetscViewerFileSetName( tViewer, aFilename );

    VecView( mPetscVector, tViewer );

    PetscViewerDestroy( &tViewer );
}
