/*
 * cl_Vector_Epetra.cpp
 *
 *  Created on: Jun 28, 2018
 *      Author: schmidt
 */
#include "cl_Vector_Epetra.hpp"
#include <string>

using namespace moris;

//----------------------------------------------------------------------------------------------

Vector_Epetra::Vector_Epetra(       sol::Dist_Map   * aMapClass,
                              const sint              aNumVectors ) : sol::Dist_Vector( aMapClass )
{
    mNumVectors = aNumVectors;

    // Build Epetra Vector
    mEpetraVector = new Epetra_FEVector( *aMapClass->get_epetra_map(), aNumVectors );

    // Get pointer to MultiVector values
    mValuesPtr = mEpetraVector->Values();
}
//----------------------------------------------------------------------------------------------

Vector_Epetra::~Vector_Epetra()
{
    delete( mEpetraVector );
}

//----------------------------------------------------------------------------------------------

void Vector_Epetra::replace_global_values( const moris::Matrix< DDSMat > & aGlobalIds,
                                           const moris::Matrix< DDRMat > & aValues,
                                           const uint                    & aVectorIndex )
{
    reinterpret_cast< Epetra_FEVector* >( mEpetraVector )->ReplaceGlobalValues( aGlobalIds.numel(),
                                                                                aGlobalIds.data(),
                                                                                aValues.data(),
                                                                                aVectorIndex );
}

//----------------------------------------------------------------------------------------------

void Vector_Epetra::sum_into_global_values( const moris::Matrix< DDSMat > & aGlobalIds,
                                            const moris::Matrix< DDRMat > & aValues,
                                            const uint                    & aVectorIndex )
{
    // sum a number (aNumMyDofs) of values (mem_pointer( aRHSVal )) into given positions (mem_pointer( aElementTopology )) of the vector
    reinterpret_cast< Epetra_FEVector* >( mEpetraVector )->SumIntoGlobalValues( aGlobalIds.numel(),
                                                                                aGlobalIds.data(),
                                                                                aValues.data(),
                                                                                aVectorIndex );
}

//----------------------------------------------------------------------------------------------

void Vector_Epetra::vector_global_asembly()
{
    // Gather any overlapping/shared data into the non-overlapping partitioning defined by the Map.
    reinterpret_cast< Epetra_FEVector* >( mEpetraVector )->GlobalAssemble();
}

//----------------------------------------------------------------------------------------------

void Vector_Epetra::vec_plus_vec( const moris::real      & aScaleA,
                                        sol::Dist_Vector & aVecA,
                                  const moris::real      & aScaleThis )
{
    // check if both vectors are build with the same map
    const Epetra_Map* tMap = aVecA.get_map()->get_epetra_map();

    if ( mMap->get_epetra_map()->PointSameAs( *tMap ) )
    {
        //currently guessing Epetra update is smart enough to switch to replace if aScaleThis is 0.0
        mEpetraVector->Update( aScaleA, *aVecA.get_epetra_vector(), aScaleThis );
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
void Vector_Epetra::import_local_to_global( sol::Dist_Vector & aSourceVec )
{
    // check if both vectores have the same map
    const Epetra_Map* tMap = aSourceVec.get_map()->get_epetra_map();

    if ( mMap->get_epetra_map()->PointSameAs( *tMap ) )
    {
        mEpetraVector->Update( 1.0, *aSourceVec.get_epetra_vector(), 0.0 );
        //MORIS_ERROR( false, "Both vectors have the same map. Use vec_plus_vec() instead" );
    }
    else
    {
        // Build importer oject
        if ( !mImporter )
        {
            mImporter = new Epetra_Import( *mMap->get_epetra_map(), *aSourceVec.get_map()->get_epetra_map() );
        }

        int status = mEpetraVector->Import( *aSourceVec.get_epetra_vector(), *mImporter, Insert );

        if ( status!=0 )
        {
            MORIS_ERROR( false, "failed to import local to global vector" );
        }
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

Cell< moris::real > Vector_Epetra::vec_norm2()
{
    Cell< moris::real > tNorm( mNumVectors, 0.0);

    // get the norm2 of this vector
    mEpetraVector->Norm2( tNorm.data().data() );

    return tNorm;
}

//----------------------------------------------------------------------------------------------

void Vector_Epetra::extract_copy( moris::Matrix< DDRMat > & LHSValues )
{
    //std::cout<<*mEpetraVector<<std::endl;

    moris::sint tVectorLenght = this->vec_local_length();

    LHSValues.set_size( tVectorLenght, mNumVectors );

    // needed as offset parameter for Epetra. =tVectorLenght
    sint tMyLDA = tVectorLenght;

    // Get solution and output it in moris::Mat LHSValues
    mEpetraVector->ExtractCopy( LHSValues.data(), tMyLDA );

}

//----------------------------------------------------------------------------------------------

void Vector_Epetra::extract_my_values( const moris::uint                            & aNumIndices,
                                       const moris::Matrix< DDSMat >                & aGlobalRows,
                                       const moris::uint                            & aRowOffsets,
                                             moris::Cell< moris::Matrix< DDRMat > > & ExtractedValues )
{
    ExtractedValues.resize( mNumVectors );

    for ( moris::sint Ik = 0; Ik < mNumVectors; ++Ik )
    {
        ExtractedValues( Ik ).set_size( aNumIndices, 1 );
    }

    moris::sint tVecLenght = this->vec_local_length();

    for ( moris::sint Ik = 0; Ik < mNumVectors; ++Ik )
    {
        moris::sint tOffset = tVecLenght * Ik;

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

            ExtractedValues( Ik )( Ii ) = mValuesPtr[ tLocIndex + tOffset ];
        }
    }
}

//----------------------------------------------------------------------------------------------

void Vector_Epetra::print() const
{
    std::cout << *mEpetraVector <<std::endl;
}

//----------------------------------------------------------------------------------------------

void Vector_Epetra::save_vector_to_matrix_market_file( const char* aFilename )
{
    EpetraExt::MultiVectorToMatrixMarketFile( aFilename, *mEpetraVector );
}

//----------------------------------------------------------------------------------------------

void Vector_Epetra::save_vector_to_HDF5( const char* aFilename )
{
    EpetraExt::HDF5 HDF5( mMap->get_epetra_map()->Comm() );
    HDF5.Create( aFilename );

    HDF5.Write("map-" + std::to_string( mMap->get_epetra_map()->Comm().NumProc() ), *mMap->get_epetra_map());
    HDF5.Write( "LHS", *mEpetraVector );

    HDF5.Close( );
}

//-----------------------------------------------------------------------------

void Vector_Epetra::read_vector_from_HDF5( const char* aFilename )
{
    Communicator_Epetra   mEpetraComm;

    EpetraExt::HDF5 HDF5( *mEpetraComm.get_epetra_comm() );
    HDF5.Open( aFilename );

    Epetra_Map * NewMap;
    HDF5.Read("map-" + std::to_string( (mEpetraComm.get_epetra_comm())->NumProc() ), NewMap );

    Epetra_MultiVector * NewVector = NULL;

    HDF5.Read( "LHS", *NewMap ,NewVector );
    HDF5.Close( );

    mEpetraVector = NewVector;
    //FIXME
//    mMap->get_epetra_map() = NewMap;
    //mMap->get_epetra_full_overlapping_map() = NewMap;
}

//----------------------------------------------------------------------------------------------

void Vector_Epetra::check_vector( )
{

}
