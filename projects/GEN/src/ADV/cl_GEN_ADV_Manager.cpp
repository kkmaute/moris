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

namespace moris::gen
{

    //--------------------------------------------------------------------------------------------------------------

    void ADV_Manager::register_parameter_ids( const Vector< char >& aParameterIDs )
    {
        // Check for finalization
        MORIS_ASSERT( not mParameterIDsFinalized, "The ADV manager cannot add additional parameter IDs after they have already been finalized." );

        // Loop over new parameter IDs
        for ( auto iParameterID : aParameterIDs )
        {
            // Add ID if it doesn't exist already
            auto tFindID = std::find( mParameterIDs.begin(), mParameterIDs.end(), iParameterID );
            if ( tFindID == mParameterIDs.end() )
            {
                mParameterIDs.push_back( iParameterID );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void ADV_Manager::finalize_parameter_ids()
    {
        // Check for finalization
        MORIS_ASSERT( not mParameterIDsFinalized, "The ADV manager can only finalize parameter IDs once." );

        // Reserve space
        mADVs.reserve( mParameterIDs.size() );
        mLowerBounds.reserve( mParameterIDs.size() );
        mUpperBounds.reserve( mParameterIDs.size() );

        // Set flag
        mParameterIDsFinalized = true;
    }

    //--------------------------------------------------------------------------------------------------------------

    ADV ADV_Manager::create_adv( const Design_Variable& aDesignVariable )
    {
        // Check for finalization
        MORIS_ASSERT( mParameterIDsFinalized, "ADVs can only be created if the parameter IDs are finalized." );

        if ( aDesignVariable.is_constant() )
        {
            // Constant ADV
            return ADV( aDesignVariable.get_value() );
        }
        else
        {
            // Find ID in the list
            auto tFindID = std::find( mParameterIDs.begin(), mParameterIDs.end(), aDesignVariable.get_id() );
            uint tFindIndex = tFindID - mParameterIDs.begin();

            // Determine what to do
            if ( tFindIndex == mADVs.size() )
            {
                // New ADV
                mADVs.push_back( aDesignVariable.get_value() );
                mLowerBounds.push_back( aDesignVariable.get_lower_bound() );
                mUpperBounds.push_back( aDesignVariable.get_upper_bound() );
                return ADV( mADVs, mADVs.size() - 1 );
            }
            else
            {
                // ADV that has the same value as a previously added ADV
                return ADV( mADVs, tFindIndex );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------
}
