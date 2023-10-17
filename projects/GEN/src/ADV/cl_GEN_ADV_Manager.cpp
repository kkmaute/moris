/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Field.cpp
 *
 */

#include "cl_GEN_ADV_Manager.hpp"
#include "cl_SOL_Matrix_Vector_Factory.hpp"
#include "cl_SOL_Dist_Map.hpp"

namespace moris::ge
{
    //--------------------------------------------------------------------------------------------------------------
    // Get addresses inside of different vector types with overloaded function
    //--------------------------------------------------------------------------------------------------------------

    real*
    get_address( Matrix< DDRMat >& aVector, uint aIndex )
    {
        return &aVector( aIndex );
    }

    real*
    get_address( sol::Dist_Vector* aVector, uint aIndex )
    {
        return &( *aVector )( aIndex );
    }

    //--------------------------------------------------------------------------------------------------------------
    // Definitions
    //--------------------------------------------------------------------------------------------------------------

    ADV_Manager::ADV_Manager(
            Matrix< DDRMat >&       aADVs,
            const Matrix< DDUMat >& aVariableIndices,
            const Matrix< DDUMat >& aADVIndices,
            const Matrix< DDRMat >& aConstants )
            : mVariables( aVariableIndices.length() + aConstants.length() )
            , mSensitivities( 1, aVariableIndices.length() + aConstants.length() )
            , mConstants( aConstants )
            , mDeterminingADVIds( aVariableIndices.length() + aConstants.length(), 1, -1 )
            , mHasADVs( aADVIndices.length() )
    {
        // Check and assign ADV dependencies
        this->assign_adv_dependencies( aVariableIndices, aADVIndices );

        // Fill with pointers to ADVs
        this->set_advs( aADVs );

        // Fill constant parameters
        this->fill_constant_parameters();
    }

    //--------------------------------------------------------------------------------------------------------------

    ADV_Manager::ADV_Manager(
            const Matrix< DDRMat >& aConstants )
            : mVariables( aConstants.length() )
            , mSensitivities( 1, aConstants.length() )
            , mConstants( aConstants )
            , mDeterminingADVIds( aConstants.length(), 1, -1 )
            , mHasADVs( false )
    {
        // Fill constant parameters
        this->fill_constant_parameters();
    }

    //--------------------------------------------------------------------------------------------------------------

    ADV_Manager::ADV_Manager(
            const Matrix< DDUMat >& aFieldVariableIndices,
            const Matrix< DDSMat >& aSharedADVIds )
            : mVariables( aSharedADVIds.length() )
            , mSensitivities( 1, aSharedADVIds.length() )
            , mDeterminingADVIds( aSharedADVIds )
            , mHasADVs( true )
    {
        // Check that the field variable indices match the shared ADV Ids
        MORIS_ERROR( aFieldVariableIndices.length() == aSharedADVIds.length(),
                "Number of field variable indices must equal the number of ADV IDs in a GEN Field." );

        // Create shared distributed vector
        sol::Matrix_Vector_Factory tDistributedFactory;
        sol::Dist_Map* tSharedADVMap = tDistributedFactory.create_map( aSharedADVIds );
        mSharedADVs = tDistributedFactory.create_vector( tSharedADVMap, 1, false, true );

        // Set variables from ADVs
        uint tNumSharedADVs = aSharedADVIds.length();
        for ( uint tVariable = 0; tVariable < tNumSharedADVs; tVariable++ )
        {
            mVariables( aFieldVariableIndices( tVariable ) ) = &( *mSharedADVs )( aSharedADVIds( tVariable ) );
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
        for ( uint tVariableIndex = 0; tVariableIndex < mDeterminingADVIds.length(); tVariableIndex++ )
        {
            if ( mDeterminingADVIds( tVariableIndex ) >= 0 )
            {
                mVariables( tVariableIndex ) = get_address( aADVs, mDeterminingADVIds( tVariableIndex ) );
            }
        }
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

    Matrix< DDSMat > ADV_Manager::get_determining_adv_ids(
            uint                    aNodeIndex,
            const Matrix< DDRMat >& aCoordinates )
    {
        return mDeterminingADVIds;
    }

    //--------------------------------------------------------------------------------------------------------------

    bool ADV_Manager::has_advs()
    {
        return mHasADVs;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    ADV_Manager::assign_adv_dependencies(
            const Matrix< DDUMat >& aVariableIndices,
            const Matrix< DDUMat >& aADVIndices )
    {
        // Check that the number of field variables indices equals the number of ADV indices
        MORIS_ERROR( aVariableIndices.length() == aADVIndices.length(),
                "Number of field variables indices must equal the number of ADV indices in a GEN ADV entity." );

        // Set ADV dependencies
        for ( uint tADVFillIndex = 0; tADVFillIndex < aVariableIndices.length(); tADVFillIndex++ )
        {
            mDeterminingADVIds( aVariableIndices( tADVFillIndex ) ) = aADVIndices( tADVFillIndex );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    ADV_Manager::fill_constant_parameters()
    {
        uint tParameterIndex = 0;
        for ( uint tVariableIndex = 0; tVariableIndex < mVariables.size(); tVariableIndex++ )
        {
            if ( mVariables( tVariableIndex ) == nullptr )
            {
                mVariables( tVariableIndex ) = &( mConstants( tParameterIndex++ ) );
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
