/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MSI_Design_Variable_Interface.cpp
 *
 */

#include "cl_MSI_Design_Variable_Interface.hpp"
#include "cl_MSI_Equation_Model.hpp"
// SOL/src
#include "cl_SOL_Dist_Vector.hpp"
// LINALG/src
#include "fn_trans.hpp"

namespace moris::MSI
{

    //------------------------------------------------------------------------------

    void
    Design_Variable_Interface::set_requested_QIs(
            const Vector< std::string >& aRequestedIQIs )
    {
        mRequestedQIs = std::move( aRequestedIQIs );
    }

    //------------------------------------------------------------------------------

    const Vector< std::string >&
    Design_Variable_Interface::get_requested_QI_names() const
    {
        return mRequestedQIs;
    }

    //------------------------------------------------------------------------------

    const Vector< real >
    Design_Variable_Interface::get_requested_QI_values() const
    {
        Vector< real > tQIValues( mRequestedQIs.size(), 0.0 );

        for ( uint iQI = 0; iQI < mRequestedQIs.size(); iQI++ )
        {
            tQIValues( iQI ) = this->get_QI( mRequestedQIs( iQI ) );
        }

        return tQIValues;
    }

    //------------------------------------------------------------------------------

    const Vector< Matrix< DDRMat > >
    Design_Variable_Interface::get_requested_QI_values_mat() const
    {
        Vector< Matrix< DDRMat > > tQIValues( mRequestedQIs.size(), Matrix< DDRMat >( 1, 1, 0.0 ) );

        for ( uint iQI = 0; iQI < mRequestedQIs.size(); iQI++ )
        {
            tQIValues( iQI )( 0, 0 ) = this->get_QI( mRequestedQIs( iQI ) );
        }

        return tQIValues;
    }

    //------------------------------------------------------------------------------

    void
    Design_Variable_Interface::register_QI( QI aQI )
    {
        // check if the name already exists
        MORIS_ERROR( mQIs.find( aQI.name() ) == mQIs.end(),
                "Design_Variable_Interface::register_QI - Quantity of interest with name %s already exists. Please use a unique name for each quantity of interest.",
                aQI.name().c_str() );

        // add to map
        mQIs.emplace( aQI.name(), std::move( aQI ) );
    }

    //------------------------------------------------------------------------------

    real
    Design_Variable_Interface::get_QI( const std::string& aQIName ) const
    {
        return mQIs.at( aQIName ).val();
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    Design_Variable_Interface::get_dQIdADV( const std::string& aQIName ) const
    {
        MORIS_ASSERT( mQIs.find( aQIName ) != mQIs.end(),
                "Design_Variable_Interface::get_dQIdADV - Quantity of interest with name %s not found.",
                aQIName.c_str() );

        return mQIs.at( aQIName ).dADV();
    }

    //------------------------------------------------------------------------------

    void
    Design_Variable_Interface::update_QI( const std::string& aQIName, real aValue )
    {
        MORIS_ASSERT( mQIs.find( aQIName ) != mQIs.end(),
                "Design_Variable_Interface::update_QI - Quantity of interest with name %s not found.",
                aQIName.c_str() );

        mQIs.at( aQIName ).set_val( aValue );
    }

    //------------------------------------------------------------------------------

    void
    Design_Variable_Interface::update_QI( const std::string& aQIName, real aValue, const Matrix< DDRMat >& adQIdADV )
    {
        MORIS_ASSERT( mQIs.find( aQIName ) != mQIs.end(),
                "Design_Variable_Interface::update_QI - Quantity of interest with name %s not found.",
                aQIName.c_str() );

        mQIs.at( aQIName ).set_val( aValue );
        mQIs.at( aQIName ).set_dADV( adQIdADV );
    }

    //------------------------------------------------------------------------------

    void
    Design_Variable_Interface::update_QI( const std::string& aQIName, const Matrix< DDRMat >& adQIdADV )
    {
        MORIS_ASSERT( mQIs.find( aQIName ) != mQIs.end(),
                "Design_Variable_Interface::update_QI - Quantity of interest with name %s not found.",
                aQIName.c_str() );

        mQIs.at( aQIName ).set_dADV( adQIdADV );
    }

    //------------------------------------------------------------------------------

    sol::Dist_Vector*
    Design_Variable_Interface::get_dQIdp()
    {
        if ( !mdQIdpImported )
        {
            MORIS_ASSERT( mModel != nullptr,
                    "Design_Variable_Interface::get_dQIdp - mModel has not been set." );

            return mModel->get_dQIdp();
        }
        else
        {
            return mdQIdp;
        }
    }

    //------------------------------------------------------------------------------

    void
    Design_Variable_Interface::set_dQIdp_dist_vect( sol::Dist_Vector* adQIdp )
    {
        mdQIdpImported = true;
        mdQIdp         = adQIdp;
    }

}    // namespace moris::MSI
