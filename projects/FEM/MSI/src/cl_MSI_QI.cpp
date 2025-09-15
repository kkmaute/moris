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
            Module_Type             aModule,
            const real&             aValue,
            const Matrix< DDRMat >& aDQIdADV )
            : mModule( aModule )
            , mValue( std::move( aValue ) )
            , mIsEvaluated( true )
    {
        if ( aDQIdADV.numel() != 0 )
        {
            mdQIdADV = std::move( aDQIdADV );
        }
    }

    //------------------------------------------------------------------------------

    QI::QI(
            Module_Type                          aModule,
            std::function< real() >              aValueFunction,
            std::function< Matrix< DDRMat >&() > aDQIdADVFunction )
            : mModule( aModule )
            , mValueFunction( std::move( aValueFunction ) )
            , mDQIdADVFunction( std::move( aDQIdADVFunction ) )
    {
    }

    //------------------------------------------------------------------------------

    QI::QI(
            Module_Type                          aModule,
            const real&                          aValue,
            std::function< Matrix< DDRMat >&() > aDQIdADVFunction )
            : mModule( aModule )
            , mValue( std::move( aValue ) )
            , mIsEvaluated( true )
            , mDQIdADVFunction( std::move( aDQIdADVFunction ) )
    {
    }

    //------------------------------------------------------------------------------

    QI::QI(
            Module_Type             aModule,
            std::function< real() > aValueFunction,
            const Matrix< DDRMat >& aDQIdADV )
            : mModule( aModule )
            , mValueFunction( std::move( aValueFunction ) )
    {
        if ( aDQIdADV.numel() != 0 )
        {
            mdQIdADV = std::move( aDQIdADV );
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

    const Matrix< DDRMat >& QI::dADV() const
    {
        if ( not mdQIdADV.has_value() )
        {
            if ( mDQIdADVFunction == nullptr )
            {
                mdQIdADV = Matrix< DDRMat >( 0, 0 );    // brendan maybe set size
            }
            else
            {
                mdQIdADV = mDQIdADVFunction();
            }
        }
        return mdQIdADV.value();
    }

    //------------------------------------------------------------------------------

    void QI::set_dADV( const Matrix< DDRMat >& aDQIdADV )
    {
        mdQIdADV = aDQIdADV;
    }

    //------------------------------------------------------------------------------

    void QI::reset()
    {
        mIsEvaluated = false;
        mdQIdADV.reset();
    }
}    // namespace moris::MSI