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
            : mConstants( aConstants )
            , mDeterminingADVIds( aVariableIndices.length() + aConstants.length(), 1, -1 )
            , mHasADVs( aADVIndices.length() )
    {
        // Check and assign ADV dependencies
        this->assign_adv_dependencies( aVariableIndices, aADVIndices );

        // Reserve ADVs
        mADVs.reserve( mDeterminingADVIds.length() );

        // Fill with pointers to ADVs
        this->set_advs( aADVs );
    }

    //--------------------------------------------------------------------------------------------------------------

    ADV_Manager::ADV_Manager(
            const Matrix< DDRMat >& aConstants )
            : mConstants( aConstants )
            , mDeterminingADVIds( aConstants.length(), 1, -1 )
            , mHasADVs( false )
    {
        // Reserve ADVs
        mADVs.reserve( aConstants.length() );

        // Set ADVs
        Matrix< DDRMat > tDummyADVs( 0, 0 );
        this->set_advs( tDummyADVs );
    }

    //--------------------------------------------------------------------------------------------------------------

    ADV_Manager::ADV_Manager(
            const Matrix< DDUMat >& aFieldVariableIndices,
            const Matrix< DDSMat >& aSharedADVIds )
            : mDeterminingADVIds( aSharedADVIds )
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
        mADVs.reserve( tNumSharedADVs );
        for ( uint tVariable = 0; tVariable < tNumSharedADVs; tVariable++ )
        {
            mADVs.emplace_back( mSharedADVs, aSharedADVIds( tVariable ) );
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
        uint tConstantIndex = 0;
        mADVs.clear();
        for ( uint tVariableIndex = 0; tVariableIndex < mDeterminingADVIds.length(); tVariableIndex++ )
        {
            if ( mDeterminingADVIds( tVariableIndex ) >= 0 )
            {
                mADVs.emplace_back( aADVs, mDeterminingADVIds( tVariableIndex ) );
            }
            else
            {
                mADVs.emplace_back( mConstants( tConstantIndex++ ) );
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
    // Explicit template instantiation
    //--------------------------------------------------------------------------------------------------------------

    template void ADV_Manager::set_advs( Matrix< DDRMat >& aADVs );
    template void ADV_Manager::set_advs( sol::Dist_Vector*& aADVs );

    //--------------------------------------------------------------------------------------------------------------
}
