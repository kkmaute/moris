/*
* Copyright (c) 2022 University of Colorado
* Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
*
*------------------------------------------------------------------------------------
*
* cl_Design_Variable.cpp
*
*/

#include "cl_Design_Variable.hpp"

#include <utility>

namespace moris
{
    //--------------------------------------------------------------------------------------------------------------

    Design_Variable::Design_Variable( real aConstantValue )
            : mID( mCount++ )
            , mValue( aConstantValue )
            , mLowerBound( aConstantValue )
            , mUpperBound( aConstantValue )
            , mIsConstant( true )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    Design_Variable::Design_Variable(
            real        aLowerBound,
            real        aInitialValue,
            real        aUpperBound,
            std::string aName )
            : mID( mCount++ )
            , mValue ( aInitialValue )
            , mLowerBound( aLowerBound )
            , mUpperBound( aUpperBound )
            , mIsConstant( false )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    bool Design_Variable::is_constant() const
    {
        return mIsConstant;
    }

    //--------------------------------------------------------------------------------------------------------------

    char Design_Variable::get_id() const
    {
        return mID;
    }

    //--------------------------------------------------------------------------------------------------------------

    real Design_Variable::get_value() const
    {
        return mValue;
    }

    //--------------------------------------------------------------------------------------------------------------

    real Design_Variable::get_lower_bound() const
    {
        return mLowerBound;
    }

    //--------------------------------------------------------------------------------------------------------------

    real Design_Variable::get_upper_bound() const
    {
        return mUpperBound;
    }
    
    //--------------------------------------------------------------------------------------------------------------

    bool operator <( const Design_Variable& aLeft, const Design_Variable& aRight )
    {
        return aLeft.get_value() < aRight.get_value()
           and aLeft.get_lower_bound() < aRight.get_lower_bound()
           and aLeft.get_upper_bound() < aRight.get_upper_bound();
    }

    //--------------------------------------------------------------------------------------------------------------

    bool operator >( const Design_Variable& aLeft, const Design_Variable& aRight )
    {
        return aLeft.get_value() > aRight.get_value()
           and aLeft.get_lower_bound() > aRight.get_lower_bound()
           and aLeft.get_upper_bound() > aRight.get_upper_bound();
    }

    //--------------------------------------------------------------------------------------------------------------
}
