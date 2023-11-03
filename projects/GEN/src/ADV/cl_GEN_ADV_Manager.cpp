/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_ADV_Manager.cpp
 *
 */

#include "cl_GEN_ADV_Manager.hpp"
#include "cl_SOL_Matrix_Vector_Factory.hpp"
#include "cl_SOL_Dist_Map.hpp"

namespace moris::ge
{
    //--------------------------------------------------------------------------------------------------------------

    ADV_Manager::ADV_Manager(
            Matrix< DDRMat >&       aADVs,
            const Matrix< DDUMat >& aVariableIndices,
            const Matrix< DDUMat >& aADVIndices,
            const Matrix< DDRMat >& aConstants )
            : mDeterminingADVIds( aVariableIndices.length() + aConstants.length(), 1, gNoID )
            , mHasADVs( aADVIndices.length() )
    {
        // Check that the number of field variables indices equals the number of ADV indices
        MORIS_ERROR( aVariableIndices.length() == aADVIndices.length(),
                "Number of field variables indices must equal the number of ADV indices in a GEN ADV_Manager." );

        // Set ADV dependencies
        for ( uint tADVFillIndex = 0; tADVFillIndex < aVariableIndices.length(); tADVFillIndex++ )
        {
            mDeterminingADVIds( aVariableIndices( tADVFillIndex ) ) = aADVIndices( tADVFillIndex );
        }

        // Fill with pointers to ADVs
        this->create_advs( aADVs, aConstants );
    }

    //--------------------------------------------------------------------------------------------------------------

    ADV_Manager::ADV_Manager(
            const Matrix< DDRMat >& aConstants )
            : mDeterminingADVIds( aConstants.length(), 1, gNoID )
            , mHasADVs( false )
    {
        // Set ADVs
        Matrix< DDRMat > tDummyADVs( 0, 0 );
        this->create_advs( tDummyADVs, aConstants );
    }

    //--------------------------------------------------------------------------------------------------------------

    ADV_Manager::ADV_Manager( const Matrix< DDSMat >& aSharedADVIds )
            : mDeterminingADVIds( aSharedADVIds )
            , mHasADVs( true )
    {
        // Create shared distributed vector
        sol::Matrix_Vector_Factory tDistributedFactory;
        sol::Dist_Map* tSharedADVMap = tDistributedFactory.create_map( aSharedADVIds );
        mSharedADVs = tDistributedFactory.create_vector( tSharedADVMap, 1, false, true );

        // Set variables from ADVs
        uint tNumSharedADVs = aSharedADVIds.length();
        mADVs.reserve( tNumSharedADVs );
        for ( uint iVariableIndex = 0; iVariableIndex < tNumSharedADVs; iVariableIndex++ )
        {
            mADVs.emplace_back( mSharedADVs, aSharedADVIds( iVariableIndex ) );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    ADV_Manager::ADV_Manager(
            const ADV_Manager& aCopyADVManager,
            const Cell< uint >& aReplaceVariables,
            const Cell< real >& aNewConstants )
            : mADVs( aCopyADVManager.mADVs )
            , mDeterminingADVIds( aCopyADVManager.mDeterminingADVIds )
            , mHasADVs( aCopyADVManager.mHasADVs )
            , mSharedADVs( aCopyADVManager.mSharedADVs )
    {
        // Ensure the number of replacement variables equals the number of new constants
        MORIS_ERROR( aReplaceVariables.size() == aNewConstants.size(),
                "ADV copy constructor must be given same amount of variable indices to replace and new constants." );

        // Replace constant variables
        for ( uint iReplacementIndex = 0; iReplacementIndex < aReplaceVariables.size(); iReplacementIndex++ )
        {
            mADVs( aReplaceVariables( iReplacementIndex ) ).replace_constant( aNewConstants( iReplacementIndex ) );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    ADV_Manager::~ADV_Manager()
    {
        delete mSharedADVs;
    }

    //--------------------------------------------------------------------------------------------------------------

    template< typename Vector_Type >
    void
    ADV_Manager::set_advs( Vector_Type& aADVs )
    {
        for ( uint iVariableIndex = 0; iVariableIndex < mDeterminingADVIds.length(); iVariableIndex++ )
        {
            if ( mDeterminingADVIds( iVariableIndex ) >= 0 )
            {
                mADVs( iVariableIndex ) = ADV( aADVs, mDeterminingADVIds( iVariableIndex ) );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    real ADV_Manager::get_variable( uint aVariableIndex )
    {
        return mADVs( aVariableIndex ).get_value();
    }

    //--------------------------------------------------------------------------------------------------------------

    void ADV_Manager::import_advs( sol::Dist_Vector* aOwnedADVs )
    {
        if ( mSharedADVs )
        {
            mSharedADVs->import_local_to_global( *aOwnedADVs );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDSMat > ADV_Manager::get_determining_adv_ids()
    {
        return mDeterminingADVIds;
    }

    //--------------------------------------------------------------------------------------------------------------

    bool ADV_Manager::has_advs()
    {
        return mHasADVs;
    }

    //--------------------------------------------------------------------------------------------------------------

    void ADV_Manager::create_advs(
            Matrix< DDRMat >& aADVs,
            const Matrix< DDRMat >& aConstants )
    {
        // Reserve ADVs
        mADVs.reserve( mDeterminingADVIds.length() );

        // Create ADVs
        uint tConstantIndex = 0;
        for ( uint iVariableIndex = 0; iVariableIndex < mDeterminingADVIds.length(); iVariableIndex++ )
        {
            if ( mDeterminingADVIds( iVariableIndex ) >= 0 )
            {
                mADVs.emplace_back( aADVs, mDeterminingADVIds( iVariableIndex ) );
            }
            else
            {
                mADVs.emplace_back( aConstants( tConstantIndex++ ) );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------
    // Explicit template instantiation
    //--------------------------------------------------------------------------------------------------------------

    template void ADV_Manager::set_advs( Matrix< DDRMat >& aADVs );
    template void ADV_Manager::set_advs( sol::Dist_Vector*& aADVs );

    //--------------------------------------------------------------------------------------------------------------
}
