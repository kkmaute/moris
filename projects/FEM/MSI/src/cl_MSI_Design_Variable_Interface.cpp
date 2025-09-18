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
#include "cl_SOL_Matrix_Vector_Factory.hpp"
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
    Design_Variable_Interface::get_all_QI_names() const
    {
        return mRequestedQIs;
    }

    //------------------------------------------------------------------------------

    const Vector< std::string >
    Design_Variable_Interface::get_QI_names( Module_Type aModule ) const
    {

        // Dry run to size output vector
        uint tNumQIsInModule = 0;
        for ( const auto& tQI : mQIs )
        {
            if ( tQI.second.mModule == aModule )
            {
                tNumQIsInModule++;
            }
        }

        // Fill the vector with the QI values
        Vector< std::string > tRequestedQINames( tNumQIsInModule );
        uint                  tIndex = 0;
        for ( const auto& tQI : mQIs )
        {
            if ( tQI.second.mModule == aModule )
            {
                tRequestedQINames( tIndex++ ) = tQI.first;
            }
        }

        return tRequestedQINames;
    }

    //------------------------------------------------------------------------------

    const Vector< real >
    Design_Variable_Interface::get_all_QI_values() const
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
    Design_Variable_Interface::get_all_QI_values_mat() const
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
    Design_Variable_Interface::register_QI( const std::string& aName, QI aQI )
    {
        // check if the name already exists
        MORIS_ERROR( mQIs.find( aName ) == mQIs.end(),
                "Design_Variable_Interface::register_QI - Quantity of interest with name %s already exists. Please use a unique name for each quantity of interest.",
                aName.c_str() );

        // add to map
        mQIs.emplace( aName, std::move( aQI ) );
    }

    //------------------------------------------------------------------------------

    real
    Design_Variable_Interface::get_QI( const std::string& aQIName ) const
    {
        return mQIs.at( aQIName ).val();
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat >
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
    Design_Variable_Interface::update_QI( const std::string& aQIName, real aValue, sol::Dist_Vector* adQI )
    {
        MORIS_ASSERT( mQIs.find( aQIName ) != mQIs.end(),
                "Design_Variable_Interface::update_QI - Quantity of interest with name %s not found.",
                aQIName.c_str() );

        mQIs.at( aQIName ).set_val( aValue );
        mQIs.at( aQIName ).set_dADV( adQI );
    }

    //------------------------------------------------------------------------------

    void
    Design_Variable_Interface::update_QI( const std::string& aQIName, sol::Dist_Vector* adQI )
    {
        MORIS_ASSERT( mQIs.find( aQIName ) != mQIs.end(),
                "Design_Variable_Interface::update_QI - Quantity of interest with name %s not found.",
                aQIName.c_str() );

        mQIs.at( aQIName ).set_dADV( adQI );
    }

    //------------------------------------------------------------------------------

    void
    Design_Variable_Interface::update_QIs( Module_Type aModule, sol::Dist_Vector* adQIdp )
    {
        // BRENDAN FIXME

        // Get the QI names for this module
        Vector< std::string > tQINames = this->get_QI_names( aModule );

        MORIS_ASSERT( tQINames.size() == (uint)adQIdp->get_num_vectors(),
                "Design_Variable_Interface::update_QIs - Number of QI names (%zu) for Module %s does not match number of vectors in adQIdp (%d).",
                tQINames.size(),
                convert_parameter_list_enum_to_string( aModule ).c_str(),
                adQIdp->get_num_vectors() );

        // Create factory for resulting distributed vector
        sol::Matrix_Vector_Factory tDistributedFactory;

        // Loop through the QIs and update them
        for ( uint iQI = 0; iQI < tQINames.size(); iQI++ )
        {
            // Make new dist vector that points to the correct vector in adQIdp
            sol::Dist_Vector* tdQIdp = tDistributedFactory.create_vector( adQIdp->get_map(), 1, false, true );    // brendan check flags

            // Copy the correct vector into the new dist vector
            tdQIdp->vec_plus_vec( 1.0, *adQIdp, 0.0 ); // wrong vector

            // Update the QI
            this->update_QI( tQINames( iQI ), tdQIdp );
        }
    }

    //------------------------------------------------------------------------------

    void
    Design_Variable_Interface::update_QIs( Vector< std::string >& aQINames, sol::Dist_Vector* adQIdp )
    {
        MORIS_ASSERT( aQINames.size() == (uint)adQIdp->get_num_vectors(),
                "Design_Variable_Interface::update_QIs - Number of QI names (%zu) does not match number of vectors in adQIdp (%d).",
                aQINames.size(),
                adQIdp->get_num_vectors() );
    }

    //------------------------------------------------------------------------------

    sol::Dist_Vector*
    Design_Variable_Interface::get_dQIdp()
    {
        return mdQIdp;
    }

    //------------------------------------------------------------------------------

    // void
    // Design_Variable_Interface::set_dQIdp_dist_vect( sol::Dist_Vector* adQIdp )
    // {
    //     mdQIdpImported = true;
    //     mdQIdp         = adQIdp;
    // }

}    // namespace moris::MSI
