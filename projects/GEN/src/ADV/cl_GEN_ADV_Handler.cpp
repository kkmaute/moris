/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_ADV_Handler.cpp
 *
 */

#include "cl_GEN_ADV_Handler.hpp"
#include "cl_SOL_Matrix_Vector_Factory.hpp"
#include "cl_SOL_Dist_Map.hpp"

namespace moris::gen
{
    //--------------------------------------------------------------------------------------------------------------

    ADV_Handler::ADV_Handler(
            Vector< real >&       aADVs,
            const Vector< uint >& aVariableIndices,
            const Vector< uint >& aADVIndices,
            const Vector< real >& aConstants )
            : mDeterminingADVIds( aVariableIndices.size() + aConstants.size(), gNoID )
            , mHasADVs( aADVIndices.size() )
    {
        // Check that the number of field variables indices equals the number of ADV indices
        MORIS_ERROR( aVariableIndices.size() == aADVIndices.size(),
                "Number of field variables indices must equal the number of ADV indices in a GEN ADV_Manager." );

        // Set ADV dependencies
        for ( uint tADVFillIndex = 0; tADVFillIndex < aVariableIndices.size(); tADVFillIndex++ )
        {
            mDeterminingADVIds( aVariableIndices( tADVFillIndex ) ) = aADVIndices( tADVFillIndex );
        }

        // Fill with pointers to ADVs
        this->create_advs( aADVs, aConstants );
    }

    //--------------------------------------------------------------------------------------------------------------

    ADV_Handler::ADV_Handler( const Vector< ADV >& aADVs )
            : mADVs( aADVs )
            , mDeterminingADVIds( aADVs.size(), gNoID )
            , mHasADVs( false )
    {
        // Set ADV dependencies
        for ( uint iLocalADVIndex = 0; iLocalADVIndex < aADVs.size(); iLocalADVIndex++ )
        {
            mDeterminingADVIds( iLocalADVIndex ) = aADVs( iLocalADVIndex ).get_id();
            if ( mDeterminingADVIds( iLocalADVIndex ) != gNoID )
            {
                mHasADVs = true;
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    ADV_Handler::ADV_Handler( const Vector< sint >& aSharedADVIds )
            : mDeterminingADVIds( aSharedADVIds )
            , mHasADVs( true )
    {
        // Create shared distributed vector
        sol::Matrix_Vector_Factory tDistributedFactory;
        sol::Dist_Map*             tSharedADVMap = tDistributedFactory.create_map( aSharedADVIds );
        mSharedADVs                              = tDistributedFactory.create_vector( tSharedADVMap, 1, false, true );

        // Set variables from ADVs
        uint tNumSharedADVs = aSharedADVIds.size();
        mADVs.reserve( tNumSharedADVs );
        for ( uint iVariableIndex = 0; iVariableIndex < tNumSharedADVs; iVariableIndex++ )
        {
            mADVs.emplace_back( mSharedADVs, aSharedADVIds( iVariableIndex ) );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    ADV_Handler::ADV_Handler(
            const ADV_Handler&    aCopyADVManager,
            const Vector< uint >& aReplaceVariables,
            const Vector< real >& aNewConstants )
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

    ADV_Handler::~ADV_Handler()
    {
        delete mSharedADVs;
    }

    //--------------------------------------------------------------------------------------------------------------

    template< typename Vector_Type >
    void
    ADV_Handler::set_advs( Vector_Type& aADVs )
    {
        for ( uint iVariableIndex = 0; iVariableIndex < mDeterminingADVIds.size(); iVariableIndex++ )
        {
            if ( mDeterminingADVIds( iVariableIndex ) >= 0 )
            {
                mADVs( iVariableIndex ) = ADV( aADVs, mDeterminingADVIds( iVariableIndex ) );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    real ADV_Handler::get_variable( uint aVariableIndex )
    {
        return mADVs( aVariableIndex ).get_value();
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< real > ADV_Handler::get_values() const
    {
        Vector< real > tValues( mADVs.size() );
        for ( uint iADV = 0; iADV < mADVs.size(); iADV++ )
        {
            tValues( iADV ) = mADVs( iADV ).get_value();
        }

        return tValues;
    }

    //--------------------------------------------------------------------------------------------------------------

    void ADV_Handler::import_advs( sol::Dist_Vector* aOwnedADVs )
    {
        if ( mSharedADVs )
        {
            mSharedADVs->import_local_to_global( *aOwnedADVs );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< sint > ADV_Handler::get_determining_adv_ids() const
    {
        return mDeterminingADVIds;
    }

    //--------------------------------------------------------------------------------------------------------------

    bool ADV_Handler::has_advs()
    {
        return mHasADVs;
    }

    //--------------------------------------------------------------------------------------------------------------

    void ADV_Handler::create_advs(
            Vector< real >&       aADVs,
            const Vector< real >& aConstants )
    {
        // Reserve ADVs
        mADVs.reserve( mDeterminingADVIds.size() );

        // Create ADVs
        uint tConstantIndex = 0;
        for ( uint iVariableIndex = 0; iVariableIndex < mDeterminingADVIds.size(); iVariableIndex++ )
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

    template void ADV_Handler::set_advs( Vector< real >& aADVs );
    template void ADV_Handler::set_advs( sol::Dist_Vector*& aADVs );

    //--------------------------------------------------------------------------------------------------------------
}    // namespace moris::gen
