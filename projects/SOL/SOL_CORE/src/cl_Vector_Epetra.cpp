/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Vector_Epetra.cpp
 *
 */

#include "cl_Vector_Epetra.hpp"
#include <string>
#include "fn_print.hpp"

using namespace moris;

//----------------------------------------------------------------------------------------------

Vector_Epetra::Vector_Epetra(
        sol::Dist_Map* aMapClass,
        const sint     aNumVectors,
        bool           aPointMap,
        bool           aManageMap )
        : sol::Dist_Vector( aManageMap )
        , mVecBuildWithPointMap( aPointMap )
{
    // store map as Epetra map
    mMap = reinterpret_cast< Map_Epetra* >( aMapClass );

    // store number of columns for multi-column vectors
    mNumVectors = aNumVectors;

    // Build Epetra Vector
    if ( aPointMap )
    {
        mEpetraVector = new Epetra_FEVector( *aMapClass->get_epetra_point_map(), aNumVectors );
    }
    else
    {
        mEpetraVector = new Epetra_FEVector( *aMapClass->get_epetra_map(), aNumVectors );
    }

    // Get pointer to MultiVector values
    mValuesPtr = mEpetraVector->Values();
}
//----------------------------------------------------------------------------------------------

Vector_Epetra::~Vector_Epetra()
{
    delete mEpetraVector;

    delete mImporter;
    mImporter = nullptr;

    if ( mManageMap )
    {
        delete mMap;
    }
}

//-----------------------------------------------------------------------------

real&
Vector_Epetra::operator()( sint aGlobalId, uint aVectorIndex )
{
    // Get offset for this vector
    uint tOffset = this->vec_local_length() * aVectorIndex;

    // Get local index
    uint tLocIndex = mMap->return_local_ind_of_global_Id( aGlobalId );

    // Return value
    return mValuesPtr[ tLocIndex + tOffset ];
}

//----------------------------------------------------------------------------------------------

void
Vector_Epetra::replace_global_values(
        const moris::Matrix< DDSMat >& aGlobalIds,
        const moris::Matrix< DDRMat >& aValues,
        const uint&                    aVectorIndex )
{
    // check for empty vector
    if ( aGlobalIds.numel() == 0 )
    {
        return;
    }

    // check for valid IDs
    MORIS_ASSERT( aGlobalIds.min() >= 0 && aGlobalIds.max() < MORIS_SINT_MAX,
            "Vector_Epetra::sum_into_global_values - invalid ID provided" );

    // call native epetra function
    int error = reinterpret_cast< Epetra_FEVector* >( mEpetraVector )->    //
                ReplaceGlobalValues(                                       //
                        aGlobalIds.numel(),
                        aGlobalIds.data(),
                        aValues.data(),
                        aVectorIndex );

    MORIS_ERROR( error == 0,
            "Vector_Epetra::replace_global_values - failed" );
}

//----------------------------------------------------------------------------------------------

void
Vector_Epetra::sum_into_global_values(
        const moris::Matrix< DDSMat >& aGlobalIds,
        const moris::Matrix< DDRMat >& aValues,
        const uint&                    aVectorIndex )
{
    // check for empty vector
    if ( aGlobalIds.numel() == 0 )
    {
        return;
    }

    // check for valid IDs
    MORIS_ASSERT( aGlobalIds.min() >= 0 && aGlobalIds.max() < MORIS_SINT_MAX,
            "Vector_Epetra::sum_into_global_values - invalid ID provided" );

    // call native epetra function
    if ( mVecBuildWithPointMap )
    {
        Matrix< IdMat > tPointFreeIds;
        mMap->translate_ids_to_free_point_ids( aGlobalIds, tPointFreeIds, false );

        // sum a number (aNumMyDofs) of values (mem_pointer( aRHSVal )) into given positions (mem_pointer( aElementTopology )) of the vector
        int error = reinterpret_cast< Epetra_FEVector* >( mEpetraVector )->    //
                    SumIntoGlobalValues(                                       //
                            tPointFreeIds.numel(),
                            tPointFreeIds.data(),
                            aValues.data(),
                            aVectorIndex );

        MORIS_ERROR( error == 0,
                "Vector_Epetra::sum_into_global_values - failed" );
    }
    else
    {
        // sum a number (aNumMyDofs) of values (mem_pointer( aRHSVal )) into given positions (mem_pointer( aElementTopology )) of the vector
        int error = reinterpret_cast< Epetra_FEVector* >( mEpetraVector )->    //
                    SumIntoGlobalValues(                                       //
                            aGlobalIds.numel(),
                            aGlobalIds.data(),
                            aValues.data(),
                            aVectorIndex );

        MORIS_ERROR( error == 0,
                "Vector_Epetra::sum_into_global_values - failed" );
    }
}

//----------------------------------------------------------------------------------------------

void
Vector_Epetra::vector_global_assembly()
{
    // Gather any overlapping/shared data into the non-overlapping partitioning defined by the Map.
    int error = reinterpret_cast< Epetra_FEVector* >( mEpetraVector )->GlobalAssemble();

    MORIS_ERROR( error == 0,
            "Vector_Epetra::vector_global_assembly - failed" );
}

//----------------------------------------------------------------------------------------------

void
Vector_Epetra::vec_plus_vec(
        const moris::real& aScaleA,
        sol::Dist_Vector&  aVecA,
        const moris::real& aScaleThis )
{
    // check if both vectors are build with the same map
    const Epetra_BlockMap* tMap = aVecA.get_map()->get_epetra_map();

    if ( mMap->get_epetra_map()->PointSameAs( *tMap ) )
    {
        // currently guessing Epetra update is smart enough to switch to replace if aScaleThis is 0.0
        mEpetraVector->Update( aScaleA, *dynamic_cast< Vector_Epetra& >( aVecA ).get_epetra_vector(), aScaleThis );
        return;
    }
    else
    {
        if ( !mVecBuildWithPointMap )
        {
            mEpetraVector->Scale( aScaleThis );

            sint tNumElements = mEpetraVector->MyLength();

            Matrix< IdMat > tIdMat( tNumElements, 1, MORIS_ID_MAX );

            for ( sint Ik = 0; Ik < tNumElements; Ik++ )
            {
                tIdMat( Ik ) = mMap->get_epetra_map()->GID( Ik );
            }

            // FIXME adjust for multivector
            Matrix< DDRMat > tValues;
            aVecA.extract_copy( tValues );

            this->sum_into_global_values(
                    tIdMat,
                    tValues,
                    0 );
        }
        else
        {
            mEpetraVector->Scale( aScaleThis );

            sint tNumElements = mEpetraVector->MyLength();

            Matrix< IdMat > tIdMat( tNumElements, 1, MORIS_ID_MAX );

            for ( sint Ik = 0; Ik < tNumElements; Ik++ )
            {
                tIdMat( Ik ) = mMap->get_epetra_map()->GID( Ik );
            }

            Matrix< IdMat > tPointIds;
            mMap->translate_ids_to_free_point_ids( tIdMat, tPointIds, false );

            // FIXME adjust for multivector
            Matrix< DDRMat > tValues;
            aVecA.extract_copy( tValues );

            this->sum_into_global_values(
                    tPointIds,
                    tValues,
                    0 );
        }
    }
}

//----------------------------------------------------------------------------------------------

void
Vector_Epetra::scale_vector(
        const moris::real& aValue,
        const moris::uint& aVecIndex )
{
    // check if index of vector is 0. might not be zero for a multivector
    if ( aVecIndex == 0 )
    {
        // scale this vector with the aValue
        mEpetraVector->Scale( aValue );
        return;
    }

    // if multivector and aVecIndex != 0 than get the pointer to all these vectors.
    double** pointers = mEpetraVector->Pointers();

    // get local length of these vectors
    moris::uint length = this->vec_local_length();
    for ( moris::uint i = 0; i < length; ++i )
    {
        // scale all values of vector number aVecIndex
        pointers[ aVecIndex ][ i ] *= aValue;
    }
}

//----------------------------------------------------------------------------------------------

// Import the local vector into the global vector object mEpetraVec
void
Vector_Epetra::import_local_to_global( sol::Dist_Vector& aSourceVec )
{
    // check if both vectors have the same map
    const Epetra_BlockMap* tMap = aSourceVec.get_map()->get_epetra_map();

    if ( mMap->get_epetra_map()->PointSameAs( *tMap ) )
    {
        mEpetraVector->Update( 1.0, *dynamic_cast< Vector_Epetra& >( aSourceVec ).get_epetra_vector(), 0.0 );
        // MORIS_ERROR( false, "Both vectors have the same map. Use vec_plus_vec() instead" );
    }
    else
    {
        // Build importer object
        if ( !mImporter )
        {
            mImporter = new Epetra_Import(
                    *mMap->get_epetra_map(),
                    *aSourceVec.get_map()->get_epetra_map() );
        }

        int status = mEpetraVector->Import(
                *dynamic_cast< Vector_Epetra& >( aSourceVec ).get_epetra_vector(),
                *mImporter,
                Insert );

        if ( status != 0 )
        {
            MORIS_ERROR( false, "failed to import local to global vector" );
        }
    }
}

//----------------------------------------------------------------------------------------------

void
Vector_Epetra::vec_put_scalar( const moris::real& aValue )
{
    // set all entries of this vector to aValue
    mEpetraVector->PutScalar( aValue );
}

//----------------------------------------------------------------------------------------------

void
Vector_Epetra::random()
{
    // set all entries to random values
    mEpetraVector->Random();
}

//----------------------------------------------------------------------------------------------

moris::sint
Vector_Epetra::vec_local_length() const
{
    // get local length
    return (moris::sint)mEpetraVector->MyLength();
}

//----------------------------------------------------------------------------------------------

moris::sint
Vector_Epetra::vec_global_length() const
{
    // get global length of this vector
    return (moris::sint)mEpetraVector->GlobalLength();
}

//----------------------------------------------------------------------------------------------

Cell< moris::real >
Vector_Epetra::vec_norm2()
{
    Cell< moris::real > tNorm( mNumVectors, 0.0 );

    // get the norm2 of this vector
    mEpetraVector->Norm2( tNorm.data().data() );

    return tNorm;
}

//----------------------------------------------------------------------------------------------

void
Vector_Epetra::extract_copy( moris::Matrix< DDRMat >& LHSValues )
{
    // std::cout<<*mEpetraVector<<std::endl;

    moris::sint tVectorLength = this->vec_local_length();

    LHSValues.set_size( tVectorLength, mNumVectors );

    // needed as offset parameter for Epetra. =tVectorLength
    sint tMyLDA = tVectorLength;

    // Get solution and output it in moris::Mat LHSValues
    mEpetraVector->ExtractCopy( LHSValues.data(), tMyLDA );
}

//----------------------------------------------------------------------------------------------

void
Vector_Epetra::extract_my_values(
        const moris::uint&                      aNumIndices,
        const moris::Matrix< DDSMat >&          aGlobalRows,
        const moris::uint&                      aRowOffsets,
        moris::Cell< moris::Matrix< DDRMat > >& ExtractedValues )
{
    ExtractedValues.resize( mNumVectors );

    for ( moris::sint Ik = 0; Ik < mNumVectors; ++Ik )
    {
        ExtractedValues( Ik ).set_size( aNumIndices, 1 );
    }

    moris::sint tVecLength = this->vec_local_length();

    for ( moris::sint Ik = 0; Ik < mNumVectors; ++Ik )
    {
        moris::sint tOffset = tVecLength * Ik;

        for ( moris::uint Ii = 0; Ii < aNumIndices; ++Ii )
        {
            const int tLocIndex = mMap->return_local_ind_of_global_Id( aGlobalRows( Ii ) );

            MORIS_ASSERT( !( tLocIndex < 0 ), "Vector_Epetra::extract_my_values: local index < 0. this is not allowed" );

            ExtractedValues( Ik )( Ii ) = mValuesPtr[ tLocIndex + tOffset ];
        }
    }
}

//----------------------------------------------------------------------------------------------

void
Vector_Epetra::print() const
{
    std::cout << *mEpetraVector << std::endl;
}

//----------------------------------------------------------------------------------------------

void
Vector_Epetra::save_vector_to_matrix_market_file( const char* aFilename )
{
    EpetraExt::MultiVectorToMatrixMarketFile( aFilename, *mEpetraVector );
}

//----------------------------------------------------------------------------------------------

void
Vector_Epetra::save_vector_to_matlab_file( const char* aFilename )
{
    EpetraExt::MultiVectorToMatlabFile( aFilename, *mEpetraVector );
}

//----------------------------------------------------------------------------------------------

void
Vector_Epetra::save_vector_to_HDF5( const char* aFilename )
{
    EpetraExt::HDF5 HDF5( mMap->get_epetra_map()->Comm() );
    HDF5.Create( aFilename );

    HDF5.Write( "map-" + std::to_string( mMap->get_epetra_map()->Comm().NumProc() ), *mMap->get_epetra_map() );
    HDF5.Write( "LHS", *mEpetraVector );

    HDF5.Close();
}

//-----------------------------------------------------------------------------

void
Vector_Epetra::read_vector_from_HDF5(
        const char* aFilename,
        std::string aGroupName,
        sint        aVectorindex )
{
    // Note: the following function uses the existing map of the vector;
    //       thus the map of the vector read from hdf5 file needs to be consistent
    //       wit the map of the current vector

    // delete old Epetra Vector
    if ( mEpetraVector )
    {
        delete mEpetraVector;
    }

    // open hdf5 file
    Communicator_Epetra mEpetraComm;

    EpetraExt::HDF5 HDF5( *mEpetraComm.get_epetra_comm() );
    HDF5.Open( aFilename );

    // read new map
    Epetra_Map* tNewMap;

    HDF5.Read( "map-" + std::to_string( ( mEpetraComm.get_epetra_comm() )->NumProc() ), tNewMap );

    // read new multi-vector
    Epetra_MultiVector* tNewVector = NULL;

    HDF5.Read( aGroupName, *tNewMap, tNewVector );

    // close hdf5 file
    HDF5.Close();

    // extract single vector from multi-vector
    if ( aVectorindex > 0 )
    {
        // get number of read vectors
        MORIS_ERROR( aVectorindex < tNewVector->NumVectors(),
                "Vector_Epetra::read_vector_from_HDF5 - requested vector does not exist" );

        // get length of vector
        moris::sint tVecLength = tNewVector->MyLength();

        // create new vector
        mEpetraVector = new Epetra_FEVector( *tNewMap, 1 );

        // start address in multi-vector
        real* tStart = tNewVector->Values() + tVecLength * aVectorindex;
        real* tEnd   = tStart + tVecLength - 1;

        // copy vector of requested index onto member vector
        std::copy( tStart, tEnd, mEpetraVector->Values() );
    }
    else
    {
        // set vector to vector read from file
        mEpetraVector = tNewVector;
    }

    // set member variables
    mNumVectors = mEpetraVector->NumVectors();
    mValuesPtr  = mEpetraVector->Values();

    // delete map read from hdf5 file
    delete tNewMap;
}

//----------------------------------------------------------------------------------------------

void
Vector_Epetra::check_vector()
{
}
