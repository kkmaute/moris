/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MSI_QI.cpp
 *
 */

#include "cl_MSI_QI.hpp"


namespace moris::MSI
{
    //------------------------------------------------------------------------------

    QI::QI(
            Module_Type             aModule,
            std::function< real() > aValueFunction )
            : mModule( aModule )
            , mValueFunction( std::move( aValueFunction ) )
    {
    }

    //------------------------------------------------------------------------------

    QI::QI(
            Module_Type aModule,
            const real  aValue )
            : mModule( aModule )
            , mValue( std::move( aValue ) )
            , mIsEvaluated( true )
    {
    }

    //------------------------------------------------------------------------------

    real QI::val() const
    {
        if ( not mIsEvaluated )
        {
            MORIS_ERROR( mValueFunction != nullptr, "QI::val() - No function to compute QI value was provided." );
            mValue       = mValueFunction();
            mIsEvaluated = true;
        }
        return mValue;
    }

    //------------------------------------------------------------------------------

    void QI::set_val( real aValue )
    {
        mValue       = aValue;
        mIsEvaluated = true;
    }

    //------------------------------------------------------------------------------

    void QI::reset()
    {
        mIsEvaluated = false;
    }
}    // namespace moris::MSI