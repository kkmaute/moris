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

    ADV ADV_Manager::create_adv( const Design_Variable& aDesignVariable )
    {
        if ( aDesignVariable.is_constant() )
        {
            return ADV( aDesignVariable.get_value() );
        }
        else
        {
            mADVs.push_back( aDesignVariable.get_value() );
            mLowerBounds.push_back( aDesignVariable.get_lower_bound() );
            mUpperBounds.push_back( aDesignVariable.get_upper_bound() );
            return ADV( mADVs, mADVs.size() - 1 );
        }
    }

    //--------------------------------------------------------------------------------------------------------------
}
