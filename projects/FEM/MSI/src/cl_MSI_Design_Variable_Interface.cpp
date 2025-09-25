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

    template<>
    Vector< std::string >
    Design_Variable_Interface::get_requested_QIs< std::string >( Module_Type aModule ) const
    {
        // Check if we want all the requested QIs or from a specific module
        if ( aModule == Module_Type::END_ENUM )
        {
            return mRequestedQIs;
        }

        // Get the indices of all the requested QIs
        Vector< std::string > tRequestedQINames;
        tRequestedQINames.reserve( mRequestedQIs.size() );

        // Go through the list of requested QIs, add the index only if the module matches the requested module
        for ( uint iQI = 0; iQI < mRequestedQIs.size(); iQI++ )
        {
            // Get the name of the requested QI and check it exists
            const std::string& tQIName = mRequestedQIs( iQI );
            MORIS_ASSERT( mQINameToIndexMap.key_exists( tQIName ),
                    "Design_Variable_Interface::get_requested_QI_indices - Quantity of interest with name %s not found.",
                    tQIName.c_str() );

            // Get the index of the requested QI
            uint tQIIndex = mQINameToIndexMap.at( tQIName );

            // Add the index to the list if the module matches or we didn't get a module specified
            if ( ( mQIs( tQIIndex ).mModule == aModule ) )
            {
                tRequestedQINames.push_back( tQIName );
            }
        }

        return tRequestedQINames;
    }

    //------------------------------------------------------------------------------

    template<>
    Vector< uint >
    Design_Variable_Interface::get_requested_QIs< uint >( Module_Type aModule ) const
    {
        // Get the indices of all the requested QIs
        Vector< uint > tRequestedQIIndices;
        tRequestedQIIndices.reserve( mRequestedQIs.size() );

        // Go through the list of requested QIs, add the index only if the module matches the requested module
        for ( uint iQI = 0; iQI < mRequestedQIs.size(); iQI++ )
        {
            // Get the name of the requested QI and check it exists
            const std::string& tQIName = mRequestedQIs( iQI );
            MORIS_ASSERT( mQINameToIndexMap.key_exists( tQIName ),
                    "Design_Variable_Interface::get_requested_QI_indices - Quantity of interest with name %s not found.",
                    tQIName.c_str() );

            // Get the index of the requested QI
            uint tQIIndex = mQINameToIndexMap.at( tQIName );

            // Add the index to the list if the module matches
            if ( ( mQIs( tQIIndex ).mModule == aModule ) or aModule == Module_Type::END_ENUM )
            {
                tRequestedQIIndices.push_back( tQIIndex );
            }
        }

        return tRequestedQIIndices;
    }

    //------------------------------------------------------------------------------

    template<>
    const Vector< real >
    Design_Variable_Interface::get_all_QI_values< real >() const
    {
        Vector< real > tQIValues( mRequestedQIs.size(), 0.0 );

        for ( uint iQI = 0; iQI < mRequestedQIs.size(); iQI++ )
        {
            tQIValues( iQI ) = this->get_QI( mRequestedQIs( iQI ) );
        }

        return tQIValues;
    }

    //------------------------------------------------------------------------------

    template<>
    const Vector< Matrix< DDRMat > >
    Design_Variable_Interface::get_all_QI_values< Matrix< DDRMat > >() const
    {
        Vector< Matrix< DDRMat > > tQIValues( mRequestedQIs.size(), Matrix< DDRMat >( 1, 1, 0.0 ) );

        for ( uint iQI = 0; iQI < mRequestedQIs.size(); iQI++ )
        {
            tQIValues( iQI )( 0, 0 ) = this->get_QI( mRequestedQIs( iQI ) );
        }

        return tQIValues;
    }

    //------------------------------------------------------------------------------

    real Design_Variable_Interface::get_QI( const std::string& aQIName ) const
    {
        MORIS_ASSERT( mQINameToIndexMap.key_exists( aQIName ),
                "Design_Variable_Interface::get_QI - Quantity of interest with name %s not found.",
                aQIName.c_str() );

        return mQIs( mQINameToIndexMap.at( aQIName ) ).val();
    }

    //------------------------------------------------------------------------------

    real Design_Variable_Interface::get_QI( uint aIndex ) const
    {
        MORIS_ASSERT( aIndex < mQIs.size(),
                "Design_Variable_Interface::get_QI - Quantity of interest index %d out of bounds.",
                aIndex );

        return mQIs( aIndex ).val();
    }

    //------------------------------------------------------------------------------

    // const Matrix< DDRMat >
    // Design_Variable_Interface::get_dQIdADV( const std::string& aQIName ) const
    // {
    //     MORIS_ASSERT( mQINameToIndexMap.key_exists( aQIName ),
    //             "Design_Variable_Interface::get_dQIdADV - Quantity of interest with name %s not found.",
    //             aQIName.c_str() );

    //     return mQIs( mQINameToIndexMap[ aQIName ] ).sensitivity();
    // }

    //------------------------------------------------------------------------------

    void Design_Variable_Interface::update_QI( const std::string& aQIName, real aValue )
    {
        MORIS_ASSERT( mQINameToIndexMap.key_exists( aQIName ),
                "Design_Variable_Interface::update_QI - Quantity of interest with name %s not found.",
                aQIName.c_str() );

        mQIs( mQINameToIndexMap[ aQIName ] ).set_val( aValue );
    }

    //------------------------------------------------------------------------------

    void Design_Variable_Interface::update_QI( const std::string& aQIName, real aValue, sol::Dist_Vector* adQI )
    {
        MORIS_ASSERT( mQINameToIndexMap.key_exists( aQIName ),
                "Design_Variable_Interface::update_QI - Quantity of interest with name %s not found.",
                aQIName.c_str() );

        mQIs( mQINameToIndexMap[ aQIName ] ).set_val( aValue );
        mQIs( mQINameToIndexMap[ aQIName ] ).set_sensitivity( adQI );
    }

    //------------------------------------------------------------------------------

    void Design_Variable_Interface::update_QI( const std::string& aQIName, sol::Dist_Vector* adQI )
    {
        MORIS_ASSERT( mQINameToIndexMap.key_exists( aQIName ),
                "Design_Variable_Interface::update_QI - Quantity of interest with name %s not found.",
                aQIName.c_str() );

        mQIs( mQINameToIndexMap[ aQIName ] ).set_sensitivity( adQI );
    }

    //------------------------------------------------------------------------------

    void Design_Variable_Interface::update_QI_sensitivity( Module_Type aModule, sol::Dist_Vector* adQIdp )
    {
        switch ( aModule )
        {
            case Module_Type::FEM:
            {
                mdIQIdPDV = adQIdp;
                break;
            }
            case Module_Type::XTK:
            {
                mdXQIdPDV = adQIdp;
                break;
            }
            case Module_Type::GEN:
            {
                mdGQIdADV = adQIdp;
                break;
            }
            default:
            {
                MORIS_ASSERT( false,
                        "Design_Variable_Interface::update_QI_sensitivity - No design variable sensitivities for module type %d.",
                        static_cast< int >( aModule ) );
                break;
            }
        }
        // MORIS_ASSERT( mModuleToQIIndicesMap.key_exists( aModule ),
        //         "Design_Variable_Interface::update_QI_sensitivity - No QIs found for module type %d.",
        //         static_cast< int >( aModule ) );

        // // Get the indices of the QIs associated with the module
        // Vector< uint > tQIIndices = mModuleToQIIndicesMap[ aModule ];

        // MORIS_ASSERT( adQIdp->get_num_vectors() == (sint)tQIIndices.size(),
        //         "Design_Variable_Interface::update_QI_sensitivity - Number of QIs from module type %d (%lu) does not match the number of vectors in adQIdp (%d).",
        //         static_cast< int >( aModule ),
        //         tQIIndices.size(),
        //         adQIdp->get_num_vectors() );

        // // Loop over the QIs in adQIdp
        // for ( uint iQI = 0; iQI < tQIIndices.size(); iQI++ )
        // {
        //     uint tQIIndex = tQIIndices( iQI );

        //     // Loop over all entries in the Dist_Vector and set the sensitivity in the corresponding QI
        //     int tGlobalLength = adQIdp->vec_global_length();

        //     for ( int iDV = 0; iDV < tGlobalLength; iDV++ )
        //     {
        //         // Set the value in the sensitivity matrix
        //         mQIs( tQIIndex ).set_sensitivity( iDV, ( *adQIdp )( iDV, iQI ) );
        //     }
        // }
    }

    //------------------------------------------------------------------------------

    sol::Dist_Vector*
    Design_Variable_Interface::get_dQIdp( Module_Type aModule )
    {
        switch ( aModule )
        {
            case Module_Type::FEM:
            {
                return mdIQIdPDV;
            }
            case Module_Type::XTK:
            {
                return mdXQIdPDV;
            }
            case Module_Type::GEN:
            {
                return mdGQIdADV;
            }
            default:
            {
                MORIS_ASSERT( false,
                        "Design_Variable_Interface::get_dQIdp - No design variable sensitivities for module type %d.",
                        static_cast< int >( aModule ) );
                return nullptr;
            }
        }
    }
    //------------------------------------------------------------------------------

    // void
    // Design_Variable_Interface::set_dQIdp_dist_vect( sol::Dist_Vector* adQIdp )
    // {
    //     mdQIdpImported = true;
    //     mdQIdp         = adQIdp;
    // }

}    // namespace moris::MSI
