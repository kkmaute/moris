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

namespace moris
{
    //--------------------------------------------------------------------------------------------------------------

    Design_Variable::Design_Variable( real aConstantValue )
            : mValue( aConstantValue )
            , mLowerBound( aConstantValue )
            , mUpperBound( aConstantValue )
            , mIsConstant( true )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    Design_Variable::Design_Variable( real aLowerBound, real aInitialValue, real aUpperBound )
            : mValue ( aInitialValue )
            , mLowerBound( aLowerBound )
            , mUpperBound( aUpperBound )
            , mIsConstant( false )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    bool Design_Variable::is_constant()
    {
        return mIsConstant;
    }

    //--------------------------------------------------------------------------------------------------------------

    real Design_Variable::get_value()
    {
        return mValue;
    }

    //--------------------------------------------------------------------------------------------------------------

    real Design_Variable::get_lower_bound()
    {
        return mLowerBound;
    }

    //--------------------------------------------------------------------------------------------------------------

    real Design_Variable::get_upper_bound()
    {
        return mUpperBound;
    }
    
    //--------------------------------------------------------------------------------------------------------------
}
