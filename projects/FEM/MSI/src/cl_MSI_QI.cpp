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
            Module_Type       aModule,
            const real        aValue,
            sol::Dist_Vector* adQI )
            : mModule( aModule )
            , mValue( std::move( aValue ) )
            , mIsEvaluated( true )
    {
        if ( adQI != nullptr )
        {
            mdQI = std::move( adQI );
        }
    }

    //------------------------------------------------------------------------------

    QI::QI(
            Module_Type                          aModule,
            std::function< real() >              aValueFunction,
            std::function< sol::Dist_Vector*() > adQIFunction )
            : mModule( aModule )
            , mValueFunction( std::move( aValueFunction ) )
            , mdQIFunction( std::move( adQIFunction ) )
    {
    }

    //------------------------------------------------------------------------------

    QI::QI(
            Module_Type                          aModule,
            const real                           aValue,
            std::function< sol::Dist_Vector*() > adQIFunction )
            : mModule( aModule )
            , mValue( std::move( aValue ) )
            , mIsEvaluated( true )
            , mdQIFunction( std::move( adQIFunction ) )
    {
    }

    //------------------------------------------------------------------------------

    QI::QI(
            Module_Type             aModule,
            std::function< real() > aValueFunction,
            sol::Dist_Vector*       adQI )
            : mModule( aModule )
            , mValueFunction( std::move( aValueFunction ) )
    {
        if ( adQI != nullptr )
        {
            mdQI = std::move( adQI );
        }
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

    bool QI::has_sensitivities() const
    {
        return ( mdQI.has_value() && mdQI.value() != nullptr ) || mdQIFunction != nullptr;
    }

    //------------------------------------------------------------------------------

    sol::Dist_Vector* QI::sensitivity() const
    {
        if ( not mdQI.has_value() )
        {
            MORIS_ERROR( mdQIFunction != nullptr, "QI::sensitivity() - No function to compute QI sensitivity was provided." );
            mdQI = mdQIFunction();
        }
        return mdQI.value();
    }

    //------------------------------------------------------------------------------

    void QI::set_sensitivity( sol::Dist_Vector* adQI )
    {
        mdQI = adQI;
    }

    //------------------------------------------------------------------------------

    void QI::set_sensitivity( uint aIndex, real aValue )
    {
        MORIS_ERROR( mdQI.has_value(), "QI::set_sensitivity - Sensitivity vector has not been set." );
        ( *mdQI.value() )( aIndex ) = aValue;
    }

    //------------------------------------------------------------------------------

    void QI::reset()
    {
        mIsEvaluated = false;
        mdQI.reset();
    }
}    // namespace moris::MSI