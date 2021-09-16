/*
 * cl_Sparse_Matrix_EpetraFECrs.cpp
 *
 *  Created on: Jun 28, 2018
 *      Author: schmidt
 */
#include "cl_Sparse_Matrix_EpetraFECrs.hpp"

extern moris::Comm_Manager gMorisComm;

using namespace moris;

// ----------------------------------------------------------------------------

Sparse_Matrix_EpetraFECrs::Sparse_Matrix_EpetraFECrs(
        Solver_Interface* aInput,
        sol::Dist_Map*    aMap,
        bool aPointMap,
        bool aBuildGraph)
: sol::Dist_Matrix( aMap ),
  mMatBuildWithPointMap( aPointMap ),
  mBuildGraph( aBuildGraph )
{
    // Fixme implement get function for nonzero rows
    //BSpline_Mesh_Base::get_number_of_basis_connected_to_basis( const moris_index aIndex )
    moris::uint nonzerosRow =5;

    //FIXME insert boolean array for BC-- insert NumGlobalElements-- size
    mDirichletBCVec.set_size  ( aInput->get_max_num_global_dofs(), 1, 0 );

    // build BC vector
    this->dirichlet_BC_vector( mDirichletBCVec, aInput->get_constrained_Ids() );


    if( mBuildGraph )
    {
        // create matrix class
        if( aPointMap )
        {
            mEpetraGraph = new Epetra_FECrsGraph( Copy, *aMap->get_epetra_point_map(), nonzerosRow );
        }
        else
        {
            mEpetraGraph = new Epetra_FECrsGraph( Copy, *aMap->get_epetra_map(), nonzerosRow );
        }
    }
    else
    {
        // create matrix class
        if( aPointMap )
        {
            mEpetraMat = new Epetra_FECrsMatrix( Copy, *aMap->get_epetra_point_map(), nonzerosRow );
        }
        else
        {
            mEpetraMat = new Epetra_FECrsMatrix( Copy, *aMap->get_epetra_map(), nonzerosRow );
        }
    }
}

Sparse_Matrix_EpetraFECrs::Sparse_Matrix_EpetraFECrs(
        const sol::Dist_Map*  aRowMap,
        const sol::Dist_Map*  aColMap  )
{
    // Fixme implement get function for nonzero rows
    //BSpline_Mesh_Base::get_number_of_basis_connected_to_basis( const moris_index aIndex )
    moris::uint nonzerosRow =5;

    // create matrix class
    mEpetraMat = new Epetra_FECrsMatrix( Copy, *aRowMap->get_epetra_map(), *aColMap->get_epetra_map(), nonzerosRow );
}

// ----------------------------------------------------------------------------------------------------------------------
/** Destructor */
Sparse_Matrix_EpetraFECrs::~Sparse_Matrix_EpetraFECrs()
{
    delete( mEpetraMat);
    //delete( mEpetraMap);
    if( mBuildGraph && mEpetraGraph!= nullptr)
    {
        delete( mEpetraGraph);
    }
}

// ----------------------------------------------------------------------------------------------------------------------
void Sparse_Matrix_EpetraFECrs::dirichlet_BC_vector(
        moris::Matrix< DDUMat >       & aDirichletBCVec,
        const moris::Matrix< DDUMat > & aMyConstraintDofs)
{
    //build vector with constraint values. unconstraint=0 constraint =1. change this to true/false
    for (moris::uint Ik=0; Ik< aMyConstraintDofs.n_rows(); Ik++)
    {
        aDirichletBCVec( aMyConstraintDofs( Ik ) )  = 1;
    }
}

// ----------------------------------------------------------------------------------------------------------------------

void Sparse_Matrix_EpetraFECrs::fill_matrix(
        const moris::uint             & aNumMyDofs,
        const moris::Matrix< DDRMat > & aA_val,
        const moris::Matrix< DDSMat > & aEleDofConectivity )
{
    if( mMatBuildWithPointMap )
    {
        Matrix< IdMat > tPointFreeIds;
        mMap->translate_ids_to_free_point_ids( aEleDofConectivity, tPointFreeIds, false );

        // insert values to matrix
        //mEpetraMat->SumIntoGlobalValues(aNumMyDofs, mem_pointer( aEleDofConectivity ), mem_pointer( aA_val ), Epetra_FECrsMatrix::ROW_MAJOR);
        mEpetraMat->SumIntoGlobalValues(
                aNumMyDofs,
                tPointFreeIds.data(),
                aA_val.data(),
                Epetra_FECrsMatrix::COLUMN_MAJOR);
    }
    else
    {
        // insert values to matrix
        //mEpetraMat->SumIntoGlobalValues(aNumMyDofs, mem_pointer( aEleDofConectivity ), mem_pointer( aA_val ), Epetra_FECrsMatrix::ROW_MAJOR);
        mEpetraMat->SumIntoGlobalValues(
                aNumMyDofs,
                aEleDofConectivity.data(),
                aA_val.data(),
                Epetra_FECrsMatrix::COLUMN_MAJOR);
    }
}

// ----------------------------------------------------------------------------------------------------------------------

void Sparse_Matrix_EpetraFECrs::insert_values(
        const Matrix<DDSMat>& aRowIDs,
        const Matrix<DDSMat>& aColumnIDs,
        const Matrix<DDRMat>& aMatrixValues)
{
    // insert values to matrix
    mEpetraMat->InsertGlobalValues(
            aRowIDs.numel(),
            aRowIDs.data(),
            aColumnIDs.numel(),
            aColumnIDs.data(),
            aMatrixValues.data(),
            Epetra_FECrsMatrix::COLUMN_MAJOR);
}

// ----------------------------------------------------------------------------------------------------------------------

void Sparse_Matrix_EpetraFECrs::sum_into_values(
        const Matrix<DDSMat>& aRowIDs,
        const Matrix<DDSMat>& aColumnIDs,
        const Matrix<DDRMat>& aMatrixValues)
{
    // insert values to matrix
    mEpetraMat->SumIntoGlobalValues(
            aRowIDs.numel(),
            aRowIDs.data(),
            aColumnIDs.numel(),
            aColumnIDs.data(),
            aMatrixValues.data(),
            Epetra_FECrsMatrix::COLUMN_MAJOR);
}

// ----------------------------------------------------------------------------------------------------------------------

void Sparse_Matrix_EpetraFECrs::matrix_global_assembly()
{
    // Assemble matrix
    mEpetraMat->GlobalAssemble();
}

// ----------------------------------------------------------------------------------------------------------------------

void Sparse_Matrix_EpetraFECrs::initial_matrix_global_assembly()
{
    if( mBuildGraph )
    {
        mEpetraGraph->GlobalAssemble();
        mEpetraMat = new Epetra_FECrsMatrix( Copy, *mEpetraGraph );
        delete( mEpetraGraph);
        mEpetraGraph = nullptr;
    }
    else
    {
        // Assemble matrix
        mEpetraMat->GlobalAssemble();
    }
}

// ----------------------------------------------------------------------------------------------------------------------

void Sparse_Matrix_EpetraFECrs::build_graph(
        const moris::uint             & aNumMyDof,
        const moris::Matrix< DDSMat > & aElementTopology )
{
    // Build temporary matrix FIXME
    //moris::Matrix< DDSMat >TempElemDofs (aNumMyDof, 1);
    //TempElemDofs = aElementTopology;

    //moris::Matrix< DDSMat > tFreeDofIds (aNumMyDof, 1, -1.0);

    //loop over elemental dofs
    //    for (moris::uint Ij=0; Ij< aNumMyDof; Ij++)
    //    {
    //        if ( aElementTopology( Ij ) < 0)
    //        {
    //            TempElemDofs( Ij ) = -1;
    //        }
    //        else if ( aElementTopology( Ij ) > (sint)(mDirichletBCVec.numel()-1) )
    //        {
    //            TempElemDofs( Ij ) = -1;
    //        }
    //        else if ( mDirichletBCVec( aElementTopology( Ij ) ) == 1)          //FIXME
    //        {
    //            TempElemDofs( Ij ) = -1;
    //        }
    //    }

    if( mMatBuildWithPointMap )
    {
        moris::Matrix< DDSMat > tFreeDofIds (aNumMyDof, 1, -1.0);

        Matrix< IdMat > tPointFreeIds;
        mMap->translate_ids_to_free_point_ids( aElementTopology, tPointFreeIds );

        // Set counter of number free dofs to 0
        moris::uint tNumFreeDofs = 0;

        for(moris::uint Ik=0 ; Ik< aNumMyDof ; Ik++)
        {
            if ( tPointFreeIds( Ik ) < 0 ) continue;                   //elemDofs

            tFreeDofIds( tNumFreeDofs ) = tPointFreeIds( Ik );
            tNumFreeDofs++;
        }

        if( mBuildGraph )
        {
            mEpetraGraph->InsertGlobalIndices( tNumFreeDofs, tFreeDofIds.data(), tNumFreeDofs, tFreeDofIds.data() );
        }
        else
        {
            // Build Zero matrix and matrix for element free dof id
            moris::Matrix< DDRMat > tZeros (aNumMyDof*aNumMyDof, 1, 0.0);

            // Fill matrix with zeros to initialize
            mEpetraMat->InsertGlobalValues(tNumFreeDofs, tFreeDofIds.data(), tZeros.data(), Epetra_FECrsMatrix::COLUMN_MAJOR);
        }
    }
    else
    {
        moris::uint tNumFreeDofs = aElementTopology.numel();

        if( mBuildGraph )
        {
            mEpetraGraph->InsertGlobalIndices( tNumFreeDofs, aElementTopology.data(), tNumFreeDofs, aElementTopology.data() );
        }
        else
        {
            // Build Zero matrix and matrix for element free dof id
            moris::Matrix< DDRMat > tZeros (aNumMyDof*aNumMyDof, 1, 0.0);

            // Fill matrix with zeros to initialize
            mEpetraMat->InsertGlobalValues(tNumFreeDofs, aElementTopology.data(), tZeros.data(), Epetra_FECrsMatrix::COLUMN_MAJOR);
        }
    }
}

// ----------------------------------------------------------------------------------------------------------------------

void Sparse_Matrix_EpetraFECrs::mat_put_scalar( const moris::real & aValue )
{
    mEpetraMat->PutScalar( aValue );
}

// ----------------------------------------------------------------------------------------------------------------------

void Sparse_Matrix_EpetraFECrs::get_diagonal( sol::Dist_Vector & aDiagVec ) const
{
    // check if matrix is filled
    if ( mEpetraMat->Filled() == false )
    {
        MORIS_ASSERT( false, "Matrix not filled, cannot extract diagonal \n" );
    }

    // extract diagonal values into vector
    int error = mEpetraMat->ExtractDiagonalCopy( *static_cast<Epetra_Vector*>(dynamic_cast<Vector_Epetra&>(aDiagVec).get_epetra_vector()) );

    if ( error != 0 )
    {
        MORIS_ASSERT( false, "SparseMatrixEpetraFECrs::Get_diagonal - extracting diagonal copy failed \n" );
    }
}

// ----------------------------------------------------------------------------------------------------------------------

void Sparse_Matrix_EpetraFECrs::sparse_mat_left_scale( const sol::Dist_Vector & aScaleVector )
{
    // scale matrix with vector from the left
    int error = mEpetraMat->LeftScale( *static_cast<Epetra_Vector*>( dynamic_cast<const Vector_Epetra&>(aScaleVector).get_epetra_vector() ) );

    if ( error != 0 )
    {
        MORIS_ASSERT( false, "SparseMatrixEpetraFECrs::sparse_mat_left_scale - Scaling matrix from left failed \n" );
    }
}

// ----------------------------------------------------------------------------------------------------------------------

void Sparse_Matrix_EpetraFECrs::sparse_mat_right_scale( const sol::Dist_Vector & aScaleVector )
{
    // scale matrix with vector from the right
    int error = mEpetraMat->RightScale( *static_cast<Epetra_Vector*>( dynamic_cast<const Vector_Epetra&>(aScaleVector).get_epetra_vector() ) );

    if ( error != 0 )
    {
        MORIS_ASSERT( false, "SparseMatrixEpetraFECrs::sparse_mat_right_scale - Scaling matrix from right failed  \n" );
    }
}

// ----------------------------------------------------------------------------

void Sparse_Matrix_EpetraFECrs::replace_diagonal_values( const sol::Dist_Vector & aDiagVec )
{
    // check if matrix is filled
    if ( mEpetraMat->Filled() == false )
    {
        MORIS_ASSERT( false, "Matrix not filled, cannot replace diagonal values \n" );
    }

    // replace diagonal matrix values with vector values
    int error = mEpetraMat->ReplaceDiagonalValues( *static_cast<Epetra_Vector*>( dynamic_cast<const Vector_Epetra&>(aDiagVec).get_epetra_vector() ) );

    if ( error != 0 )
    {
        MORIS_ASSERT( false, "SparseMatrixEpetraFECrs::replace_diagonal_values - replacing diagonal was failed!!! \n" );
    }
}

// ----------------------------------------------------------------------------------------------------------------------

void Sparse_Matrix_EpetraFECrs::mat_vec_product(
        const moris::sol::Dist_Vector & aInputVec,
        moris::sol::Dist_Vector       & aResult,
        const bool                      aUseTranspose )
{
    mEpetraMat->Multiply(
            aUseTranspose,
            *dynamic_cast<const Vector_Epetra&>(aInputVec).get_epetra_vector(),
            *dynamic_cast<Vector_Epetra&>(aResult).get_epetra_vector() );
}

// ----------------------------------------------------------------------------------------------------------------------

void  Sparse_Matrix_EpetraFECrs::print() const
{
    std::cout << *mEpetraMat <<std::endl;
}

// ----------------------------------------------------------------------------------------------------------------------

void Sparse_Matrix_EpetraFECrs::save_matrix_to_matlab_file( const char* aFilename )
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
    EpetraExt::BlockMapToMatrixMarketFile( aFilename, *mMap->get_epetra_map() );
}

