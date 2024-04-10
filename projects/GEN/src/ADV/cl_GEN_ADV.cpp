/*
* Copyright (c) 2022 University of Colorado
* Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
*
*------------------------------------------------------------------------------------
*
* cl_GEN_ADV.cpp
*
*/

#include "cl_GEN_ADV.hpp"
#include "cl_Matrix.hpp"
#include "cl_SOL_Dist_Vector.hpp"

namespace moris::gen
{
    //--------------------------------------------------------------------------------------------------------------

    ADV::ADV( Vector< real >& aADVs, sint aID )
            : mValue( &aADVs( aID ) )
            , mID( aID )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    ADV::ADV( sol::Dist_Vector* aADVs, sint aID )
            : mValue( &( *aADVs )( aID ) )
            , mID( aID )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    ADV::ADV( real aValue )
            : mValue( new real( aValue ) )
            , mID( gNoID )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    ADV::ADV( const ADV& aCopyADV )
            : mValue( aCopyADV.mID == gNoID ? new real( *aCopyADV.mValue ) : aCopyADV.mValue )
            , mID( aCopyADV.mID )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    ADV::~ADV()
    {
        if ( mID == gNoID )
        {
            delete mValue;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    real ADV::get_value()
    {
        return *mValue;
    }

    //--------------------------------------------------------------------------------------------------------------

    void ADV::replace_constant( real aNewValue )
    {
        if ( mID == gNoID )
        {
            *mValue = aNewValue;
        }
    }

    //--------------------------------------------------------------------------------------------------------------
}
