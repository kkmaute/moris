/*
 * cl_Vector_Epetra.cpp
 *
 *  Created on: Jun 28, 2018
 *      Author: schmidt
 */
#include "cl_Vector_Epetra.hpp"

using namespace moris;

Vector_Epetra::Vector_Epetra( const Map_Class       * aMapClass,
                              const enum VectorType   aVectorType) : Dist_Vector( aMapClass )
{
    // Build Epetra Vector
    if ( aVectorType == VectorType::FREE )
    {
        mEpetraVector = new Epetra_FEVector( *aMapClass->get_epetra_free_map(), true );

        // Get pointer to epetra free map
        mEpetraMap = aMapClass->get_epetra_free_map();
    }
    else if ( aVectorType == VectorType::FULL )
    {
        mEpetraVector = new Epetra_FEVector( *aMapClass->get_epetra_full_map(), true );

        // Get pointer to epetra free map
        mEpetraMap = aMapClass->get_epetra_full_map();
    }
    else if ( aVectorType == VectorType::FULL_OVERLAPPING )
    {
        mEpetraVector = new Epetra_FEVector( *aMapClass->get_epetra_full_overlapping_map(), true );

        // Get pointer to epetra free map
        mEpetraMap = aMapClass->get_epetra_full_overlapping_map();
    }
    else
    {
        MORIS_ERROR( false, "Dist_Vector type not implemented. Use VectorType::FREE or VectorType::FULL" );
    }

//    // Get pointer to epetra free map
//    mEpetraMap = aMapClass->get_epetra_free_map();

    // Get pointer to MultiVector values
    mValuesPtr = mEpetraVector->Values();
}

//----------------------------------------------------------------------------------------------
Vector_Epetra::~Vector_Epetra()
{
    delete( mEpetraVector );
}

//----------------------------------------------------------------------------------------------
void Vector_Epetra::replace_global_values()
{
    MORIS_ASSERT( false, "replace_global_values not implemented yet" );
}

//----------------------------------------------------------------------------------------------
void Vector_Epetra::sum_into_global_values( const moris::uint               & aNumMyDofs,
                                            const moris::Mat< int >         & aElementTopology,
                                            const moris::Mat< moris::real > & aRHSVal )
{
    // sum a nuber (aNumMyDofs)  of values (mem_pointer( aRHSVal )) into given positions (mem_pointer( aElementTopology )) of the vector
    mEpetraVector->SumIntoGlobalValues( aNumMyDofs, mem_pointer( aElementTopology ), mem_pointer( aRHSVal ) );
}

//----------------------------------------------------------------------------------------------
void Vector_Epetra::vector_global_asembly()
{
    // Gather any overlapping/shared data into the non-overlapping partitioning defined by the Map.
    mEpetraVector->GlobalAssemble();
}

//----------------------------------------------------------------------------------------------
void Vector_Epetra::vec_plus_vec( const moris::real & aScaleA,
                                        Dist_Vector & aVecA,
                                  const moris::real & aScaleThis )
{
    // check if both vectors are build with the same map
    const Epetra_Map* tMap = aVecA.get_vector_map();

    if ( mEpetraMap->PointSameAs( *tMap ) )
    {
        //currently guessing Epetra update is smart enough to switch to replace if aScaleThis is 0.0
        mEpetraVector->Update( aScaleA, *aVecA.get_vector(), aScaleThis );
        return;
    }
    else
    {
        MORIS_ASSERT( false, "Update option not implemented for vectors with different maps, yet" );
    }
}

//----------------------------------------------------------------------------------------------
void Vector_Epetra::scale_vector( const moris::real & aValue,
                                  const moris::uint & aVecIndex )
{
    // check if index of vector is 0. might not be zero for a multivector
    if ( aVecIndex==0 )
    {
        // scale this vector with the aValue
        mEpetraVector->Scale( aValue );
        return;
    }

    // if multivector and aVecIndex != 0 than get the pointer to all these vectors.
    double **pointers = mEpetraVector->Pointers();

    // get local length of these vectors
    moris::uint length = this->vec_local_length();
    for ( moris::uint i=0; i < length; ++i )
    {
        // scale all values of vector number aVecIndex
        pointers[ aVecIndex ][i] *= aValue;
    }
}

//----------------------------------------------------------------------------------------------
// Import the local vector into the global vector object mEpetraVec
void Vector_Epetra::import_local_to_global( const Dist_Vector & aSourceVec )
{
    // check if both vectores have the same map
    const Epetra_Map* tMap = aSourceVec.get_vector_map();

    //std::cout<<*tMap<<std::endl;
    //std::cout<<*mEpetraMap<<std::endl;
    if ( mEpetraMap->PointSameAs( *tMap ) )
    {
        MORIS_ERROR( false, "Both vectors have the same map. Use vec_plus_vec() instead" );
    }

    // Build importer oject
    if ( !mImporter )
    {
        mImporter = new Epetra_Import( *mEpetraMap, *aSourceVec.get_vector_map() );
    }

    int status = mEpetraVector->Import( *aSourceVec.get_vector(), *mImporter, Insert );
    if ( status!=0 )
    {
        MORIS_ERROR( false, "failed to import local to global vector" );
    }
}

//----------------------------------------------------------------------------------------------
void Vector_Epetra::vec_put_scalar( const moris::real & aValue )
{
    // set all entries of this vector to aValue
    mEpetraVector->PutScalar( aValue );
}

//----------------------------------------------------------------------------------------------
moris::sint Vector_Epetra::vec_local_length() const
{
    // get local length
    return ( moris::sint ) mEpetraVector->MyLength();
}

//----------------------------------------------------------------------------------------------
moris::sint Vector_Epetra::vec_global_length() const
{
    // get global lengt of this vector
    return ( moris::sint ) mEpetraVector->GlobalLength();
}

//----------------------------------------------------------------------------------------------
moris::real Vector_Epetra::vec_norm2()
{
    moris::real tNorm = 0.0;

    // get the norm2 of this vector
    mEpetraVector->Norm2( &tNorm );

    return tNorm;
}

//----------------------------------------------------------------------------------------------
void Vector_Epetra::extract_copy( moris::Mat< moris::real > & LHSValues )
{
    //std::cout<<*mEpetraVector<<std::endl;

    LHSValues.set_size( this->vec_local_length(), 1 );

    // needed as offset parameter for Epetra. =0
    sint tMyLDA = 0;

    // Get solution and output it in moris::Mat LHSValues
    mEpetraVector->ExtractCopy( mem_pointer( LHSValues ), tMyLDA );

}

//----------------------------------------------------------------------------------------------
void Vector_Epetra::extract_my_values( const moris::uint               & aNumIndices,
                                       const moris::Mat< moris::sint > & aGlobalRows,
                                       const moris::uint               & aRowOffsets,
                                             moris::Mat< moris::real > & LHSValues )
{
    LHSValues.set_size( aNumIndices, 1 );

    for ( moris::uint Ii = 0; Ii < aNumIndices; ++Ii )
    {
        const int tLocIndex =  mMap->return_local_ind_of_global_Id( aGlobalRows( Ii, 0) );

        MORIS_ASSERT( !( tLocIndex < 0 ), "Vector_Epetra::extract_my_values: local index < 0. this is not allowed");

//            if (!offsets)
//            {
//                LHSValues[i] = mValuesPtr[locIndex];
//                continue;
//            }

        MORIS_ASSERT( !( aRowOffsets < 0 ), "Vector_Epetra::extract_my_values: offset < 0. this is not allowed");

        LHSValues( Ii, 0) = mValuesPtr[ tLocIndex + aRowOffsets ];
    }
}

//----------------------------------------------------------------------------------------------
void Vector_Epetra::save_vector_to_matrix_market_file( const char* aFilename )
{
    EpetraExt::MultiVectorToMatrixMarketFile( aFilename, *mEpetraVector );
}

//----------------------------------------------------------------------------------------------
void Vector_Epetra::check_vector( )
{
    if ( mPetscVector != NULL )
    {
        MORIS_ASSERT( false, "epetra vector should not have any input on the petsc vector" );
    }
}
