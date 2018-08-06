/*
 * cl_VectorPETSc.cpp
 *
 *  Created on: Mar 25, 2018
 *      Author: schmidt
 */
#include "cl_VectorPETSc.hpp"

Vector_PETSc::Vector_PETSc(       moris::Solver_Input * aInput,
                            const moris::Map_Class    * aMap,
                            const enum moris::VectorType       aVectorType ) : moris::Dist_Vector( aMap )
{
    //PetscScalar    tZero = 0;
    moris::uint               aNumMyDofs          = aInput->get_num_my_dofs();
    moris::Mat< int >         aMyLocaltoGlobalMap = aInput->get_my_local_global_map();
    moris::Mat< moris::uint > aMyConstraintDofs   = aInput->get_constr_dof();
    // Get PETSc communicator
//    PetscMPIInt                rank;
//    MPI_Comm_rank(mComm->GetPETScComm(), &rank);
//    MPI_Comm_split(mComm->GetPETScComm(), rank%1, 0, &PETSC_COMM_WORLD);

    moris::uint tNumGlobalDofs=   aNumMyDofs;

    // sum up all distributed dofs
#ifdef MORIS_HAVE_PARALLEL
        MPI_Allreduce(&aNumMyDofs,&tNumGlobalDofs,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#endif

    //FIXME insert boolian array for BC-- insert NumGlobalElements-- size
    DirichletBCVec.set_size( tNumGlobalDofs, 1, 0 );

    // build BC vector
    this->dirichlet_BC_vector( DirichletBCVec, aMyConstraintDofs );

    // Set up RHS b
    VecCreateMPI( PETSC_COMM_WORLD, aNumMyDofs, PETSC_DETERMINE, &mPetscVector );
    VecSetFromOptions( mPetscVector );
    VecSetUp( mPetscVector );
}

Vector_PETSc::~Vector_PETSc()
{
    VecDestroy( &mPetscVector );
}

void Vector_PETSc::sum_into_global_values(const moris::uint               & aNumMyDof,
                                          const moris::Mat< int >         & aEleDofConectivity,
                                          const moris::Mat< moris::real > & aRHSVal)
{
    moris::Mat< int >tTempElemDofs       ( aNumMyDof,    1 );
    tTempElemDofs = aEleDofConectivity;

    //loop over elemental dofs
        for ( moris::uint Ij=0; Ij< aNumMyDof; Ij++ )
        {
            //set constrDof to neg value
            if (DirichletBCVec( aEleDofConectivity( Ij, 0 ), 0 ) == 1 )
             {
                 tTempElemDofs( Ij, 0) = -1;
             }
        }

    // Apply PETSc map AO
    AOApplicationToPetsc( mMap->get_petsc_map(), aNumMyDof, mem_pointer( tTempElemDofs ) );

    // Insert values into vector
    VecSetValues( mPetscVector, aNumMyDof, mem_pointer( tTempElemDofs ), mem_pointer( aRHSVal ), ADD_VALUES );
}

void Vector_PETSc::dirichlet_BC_vector(       moris::Mat< moris::uint > & aDirichletBCVec,
                                        const moris::Mat< uint >        & aMyConstraintDofs )
{
//build vector with constraint values. unconstraint=0 constraint =1. change this to true/false
for ( moris::uint Ik=0; Ik< aMyConstraintDofs.n_rows(); Ik++ )
    {
    aDirichletBCVec( aMyConstraintDofs( Ik, 0 ), 0 )  = 1;
    }
}

void Vector_PETSc::vector_global_asembly()
{
    VecAssemblyBegin( mPetscVector );
    VecAssemblyEnd  ( mPetscVector );
}

void Vector_PETSc::vec_plus_vec( const moris::real & aScaleA,
                                       Dist_Vector & aVecA,
                                 const moris::real & aScaleThis )
{
    PetscScalar tValueA = aScaleA;
    PetscScalar tValueThis = aScaleThis;
    VecAXPBY( mPetscVector, tValueA, tValueThis, aVecA.get_petsc_vector() );
}

void Vector_PETSc::scale_vector( const moris::real & aValue,
                                 const moris::uint & aVecIndex )
{
    PetscScalar tValue = aValue;
    VecScale( mPetscVector, tValue );
}

void Vector_PETSc::vec_put_scalar( const moris::real & aValue )
{
    PetscScalar tValue = aValue;
    VecSet( mPetscVector, tValue);
}

moris::sint Vector_PETSc::vec_local_length() const
{
    moris::sint * tVecLocSize = 0;
    VecGetLocalSize( mPetscVector, tVecLocSize);
    return *tVecLocSize;
}

moris::sint Vector_PETSc::vec_global_length() const
{
    moris::sint * tVecSize = 0;
    VecGetSize( mPetscVector, tVecSize );
    return  *tVecSize;
}

moris::real Vector_PETSc::vec_norm2()
{
    moris::real * tVecNorm = 0 ;
    VecNorm( mPetscVector, NORM_2, tVecNorm );
    return *tVecNorm;
}

void Vector_PETSc::check_vector( )
{
    if ( mEpetraVector != NULL )
    {
        MORIS_ASSERT( false, "epetra vector should not have any input on the petsc vector" );
    }
}

