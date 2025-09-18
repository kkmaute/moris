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
    QI::QI(
            uint              aIndex,
            Module_Type       aModule,
            const real        aValue,
            sol::Dist_Vector* adQI )
            : mIndex( aIndex )
            , mModule( aModule )
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
            uint                                 aIndex,
            Module_Type                          aModule,
            std::function< real() >              aValueFunction,
            std::function< sol::Dist_Vector*() > adQIFunction )
            : mIndex( aIndex )
            , mModule( aModule )
            , mValueFunction( std::move( aValueFunction ) )
            , mdQIFunction( std::move( adQIFunction ) )
    {
    }

    //------------------------------------------------------------------------------

    QI::QI(
            uint                                 aIndex,
            Module_Type                          aModule,
            const real                           aValue,
            std::function< sol::Dist_Vector*() > adQIFunction )
            : mIndex( aIndex )
            , mModule( aModule )
            , mValue( std::move( aValue ) )
            , mIsEvaluated( true )
            , mdQIFunction( std::move( adQIFunction ) )
    {
    }

    //------------------------------------------------------------------------------

    QI::QI(
            uint                    aIndex,
            Module_Type             aModule,
            std::function< real() > aValueFunction,
            sol::Dist_Vector*       adQI )
            : mIndex( aIndex )
            , mModule( aModule )
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

    const Matrix< DDRMat > QI::dADV() const
    {
        Matrix< DDRMat > tdQI;
        if ( not mdQI.has_value() )
        {
            if ( mdQIFunction == nullptr )
            {
                return tdQI;    // return empty matrix
            }
            else
            {
                mdQI = mdQIFunction();
            }
        }

        mdQI.value()->extract_copy( tdQI );
        return tdQI;
    }

    //------------------------------------------------------------------------------

    void QI::set_dADV( sol::Dist_Vector* adQI )
    {
        mdQI = adQI;
    }

    //------------------------------------------------------------------------------

    void QI::reset()
    {
        mIsEvaluated = false;
        mdQI.reset();
    }
}    // namespace moris::MSI