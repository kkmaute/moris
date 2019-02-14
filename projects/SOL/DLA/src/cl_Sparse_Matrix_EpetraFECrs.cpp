/*
 * cl_Sparse_Matrix_EpetraFECrs.cpp
 *
 *  Created on: Jun 28, 2018
 *      Author: schmidt
 */
#include "cl_Sparse_Matrix_EpetraFECrs.hpp"

extern moris::Comm_Manager gMorisComm;

using namespace moris;

Sparse_Matrix_EpetraFECrs::Sparse_Matrix_EpetraFECrs(       Solver_Interface * aInput,
                                                      const Map_Class        * aMap ) : Sparse_Matrix( aMap )
{
    // Fixme implement get function for nonzero rows
    //BSpline_Mesh_Base::get_number_of_basis_connected_to_basis( const moris_index aIndex )
    moris::uint nonzerosRow =2;

    moris::uint aNumMyDofs = aInput->get_my_local_global_map().n_rows();

    moris::uint tNumMyDofs     = aNumMyDofs;
    moris::uint tNumGlobalDofs = aNumMyDofs;

    // sum up all distributed dofs
    sum_all( tNumMyDofs, tNumGlobalDofs );

    //FIXME insert boolian array for BC-- insert NumGlobalElements-- size
    mDirichletBCVec.set_size  ( aInput->get_max_num_global_dofs(), 1, 0 );

    // build BC vector
    this->dirichlet_BC_vector( mDirichletBCVec, aInput->get_constr_dof() );

    // create matrix class
    mEpetraMat = new Epetra_FECrsMatrix( Copy, *aMap->mFreeEpetraMap, nonzerosRow );

    //mEpetraMap = ( Map_Epetra * ) aMap;
}

// ----------------------------------------------------------------------------------------------------------------------
/** Destructor */
Sparse_Matrix_EpetraFECrs::~Sparse_Matrix_EpetraFECrs()
{
    delete( mEpetraMat);
    //delete( mEpetraMap);
    //delete( mEpetraGraph);
}

// ----------------------------------------------------------------------------------------------------------------------
void Sparse_Matrix_EpetraFECrs::dirichlet_BC_vector(       moris::Matrix< DDUMat > & aDirichletBCVec,
                                                     const moris::Matrix< DDUMat > & aMyConstraintDofs)
 {
     //build vector with constraint values. unconstraint=0 constraint =1. change this to true/false
     for (moris::uint Ik=0; Ik< aMyConstraintDofs.n_rows(); Ik++)
     {
         aDirichletBCVec( aMyConstraintDofs( Ik,0)     ,0 )  = 1;
     }
 }

// ----------------------------------------------------------------------------------------------------------------------

void Sparse_Matrix_EpetraFECrs::fill_matrix( const moris::uint             & aNumMyDofs,
                                             const moris::Matrix< DDRMat > & aA_val,
                                             const moris::Matrix< DDSMat > & aEleDofConectivity )
{
    // insert values to matrix
    //mEpetraMat->SumIntoGlobalValues(aNumMyDofs, mem_pointer( aEleDofConectivity ), mem_pointer( aA_val ), Epetra_FECrsMatrix::ROW_MAJOR);
    mEpetraMat->SumIntoGlobalValues(aNumMyDofs, aEleDofConectivity.data(), aA_val.data(), Epetra_FECrsMatrix::COLUMN_MAJOR);
}

// ----------------------------------------------------------------------------------------------------------------------
void Sparse_Matrix_EpetraFECrs::matrix_global_asembly()
{
    // Assemble matrix
    mEpetraMat->GlobalAssemble();

    //std::cout<<*mEpetraMat<<std::endl;
}

// ----------------------------------------------------------------------------------------------------------------------
void Sparse_Matrix_EpetraFECrs::build_graph( const moris::uint             & aNumMyDof,
                                             const moris::Matrix< DDSMat > & aElementTopology )
{
   // Build temporary matrix FIXME
   moris::Matrix< DDSMat >TempElemDofs (aNumMyDof, 1);
   TempElemDofs = aElementTopology;

   // Build Zero matrix and matrix for element free dof id
   moris::Matrix< DDRMat > tZeros (aNumMyDof*aNumMyDof, 1, 0.0);
   moris::Matrix< DDSMat > tFreeDofIds (aNumMyDof, 1, -1.0);

   //loop over elemental dofs
   for (moris::uint Ij=0; Ij< aNumMyDof; Ij++)
   {
        //set constrDof to neg value
       if ( aElementTopology(Ij,0) < 0)
       {
           TempElemDofs( Ij, 0) = -1;
       }
       else if ( aElementTopology(Ij,0) > (sint)(mDirichletBCVec.length()-1) )
       {
           TempElemDofs( Ij, 0) = -1;
       }
       else if ( mDirichletBCVec( aElementTopology(Ij,0), 0) == 1)          //FIXME
       {
           TempElemDofs( Ij, 0) = -1;
       }
   }

   // Set counter of number free dofs to 0
   moris::uint tNumFreeDofs = 0;

   for(moris::uint Ik=0 ; Ik< aNumMyDof ; Ik++)
   {
       //if (!GMultigrid==true)
       //{
           if ( TempElemDofs(Ik,0) < 0 ) continue;                   //elemDofs
       //}
       tFreeDofIds(tNumFreeDofs,0) = TempElemDofs(Ik,0);
       tNumFreeDofs++;
   }
   // Fill matrix with zeros to initialize
   mEpetraMat->InsertGlobalValues(tNumFreeDofs, tFreeDofIds.data(), tZeros.data(), Epetra_FECrsMatrix::COLUMN_MAJOR);
}

// ----------------------------------------------------------------------------------------------------------------------
void Sparse_Matrix_EpetraFECrs::mat_put_scalar( const moris::real & aValue )
{
    mEpetraMat->PutScalar( aValue );
}

// ----------------------------------------------------------------------------------------------------------------------

void Sparse_Matrix_EpetraFECrs::get_diagonal( Dist_Vector & aDiagVec ) const
{
    // check if matrix is filled
    if ( mEpetraMat->Filled() == false )
    {
        MORIS_ASSERT( false, "Matrix not filled, cannot extract diagonal \n" );
    }

    // extract diagonal values into vector
    int error = mEpetraMat->ExtractDiagonalCopy( ( Epetra_Vector & ) *(aDiagVec.get_vector()) );

    if ( error != 0 )
    {
        MORIS_ASSERT( false, "SparseMatrixEpetraFECrs::Get_diagonal - extracting diagonal copy failed \n" );
    }
}

// ----------------------------------------------------------------------------------------------------------------------

void Sparse_Matrix_EpetraFECrs::sparse_mat_left_scale( const Dist_Vector & aScaleVector )
{
    // scale matrix with vector from the left
    int error = mEpetraMat->LeftScale( ( Epetra_Vector & ) *aScaleVector.get_vector() );

    if ( error != 0 )
    {
        MORIS_ASSERT( false, "SparseMatrixEpetraFECrs::sparse_mat_left_scale - Scaling matrix from left failed \n" );
    }
}

// ----------------------------------------------------------------------------------------------------------------------

void Sparse_Matrix_EpetraFECrs::sparse_mat_right_scale( const Dist_Vector & aScaleVector )
{
    // scale matrix with vector from the right
    int error = mEpetraMat->RightScale( ( Epetra_Vector & ) *aScaleVector.get_vector() );

    if ( error != 0 )
    {
        MORIS_ASSERT( false, "SparseMatrixEpetraFECrs::sparse_mat_right_scale - Scaling matrix from right failed  \n" );
    }
}

void Sparse_Matrix_EpetraFECrs::replace_diagonal_values( const Dist_Vector & aDiagVec )
{
    // check if matrix is filled
    if ( mEpetraMat->Filled() == false )
    {
        MORIS_ASSERT( false, "Matrix not filled, cannot replace diagonal values \n" );
    }

    // replace diagonal matrix values with vector values
    int error = mEpetraMat->ReplaceDiagonalValues( ( const Epetra_Vector & ) *(aDiagVec.get_vector()) );

    if ( error != 0 )
    {
        MORIS_ASSERT( false, "SparseMatrixEpetraFECrs::replace_diagonal_values - replacing diagonal was failed!!! \n" );

    }
}

// ----------------------------------------------------------------------------------------------------------------------

void  Sparse_Matrix_EpetraFECrs::print_matrix_to_screen() const
{
    std::cout << *mEpetraMat <<std::endl;
}

// ----------------------------------------------------------------------------------------------------------------------

void Sparse_Matrix_EpetraFECrs::save_matrix_to_matlab_file( char* aFilename )
{
    EpetraExt::RowMatrixToMatlabFile( aFilename, *mEpetraMat );
}

// ----------------------------------------------------------------------------------------------------------------------

void Sparse_Matrix_EpetraFECrs::save_matrix_to_matrix_market_file( const char* aFilename )
{
    EpetraExt::RowMatrixToMatrixMarketFile( aFilename, *mEpetraMat );
}

// ----------------------------------------------------------------------------------------------------------------------

void Sparse_Matrix_EpetraFECrs::save_matrix_map_to_matrix_market_file( const char* aFilename )
{
    EpetraExt::BlockMapToMatrixMarketFile( aFilename, *mMap->mFreeEpetraMap );
}

