/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Parameter_List.cpp
 *
 */

#include "cl_Parameter_List.hpp"

namespace moris
{

    //--------------------------------------------------------------------------------------------------------------

    void
    Parameter_List::erase( const std::string& aKey )
    {
        mParamMap.erase( aKey );
    }

    //--------------------------------------------------------------------------------------------------------------

    void Parameter_List::copy_parameters( const Parameter_List& aParameterList )
    {
        for ( const auto& iCopyIterator : aParameterList )
        {
            // Get iterator from this map
            auto tFoundIterator = mParamMap.find( iCopyIterator.first );

            // Determine if parameter needs to be insert or set
            if ( tFoundIterator == mParamMap.end() )
            {
                // Insert parameter
                mParamMap.insert( { iCopyIterator.first, iCopyIterator.second } );
            }
            else
            {
                // Set parameter
                tFoundIterator->second = iCopyIterator.second;
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    [[nodiscard]] bool Parameter_List::exists( const std::string& aKey ) const
    {
        return mParamMap.find( aKey ) not_eq mParamMap.end();
    }

    //--------------------------------------------------------------------------------------------------------------

    moris::sint
    Parameter_List::which( const std::string& aKey )
    {
        auto tIterator = mParamMap.find( aKey );

        // throw error
        MORIS_ERROR( tIterator != mParamMap.end(),
                "The requested parameter %s does not exist.\n",
                aKey.c_str() );

        return tIterator->second->which();
    }

    //--------------------------------------------------------------------------------------------------------------

    auto
    Parameter_List::begin() const -> decltype( mParamMap.begin() )
    {
        return mParamMap.begin();
    }

    //--------------------------------------------------------------------------------------------------------------

    auto
    Parameter_List::end() const -> decltype( mParamMap.end() )
    {
        return mParamMap.end();
    }

    //--------------------------------------------------------------------------------------------------------------

    bool
    Parameter_List::is_empty()
    {
        return mParamMap.empty();
    }

    //--------------------------------------------------------------------------------------------------------------

    size_t
    Parameter_List::size()
    {
        return mParamMap.size();
    }

    //--------------------------------------------------------------------------------------------------------------

}
