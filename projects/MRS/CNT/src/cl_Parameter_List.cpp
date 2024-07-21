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

#include <utility>

namespace moris
{
    //--------------------------------------------------------------------------------------------------------------

    void Parameter_List::insert(
            const std::string&           aName,
            const std::string&           aDefaultValue,
            const Vector< std::string >& aValidSelections )
    {
        // Check for leading and trailing whitespaces in key
        std::string tKeyWithoutSpaces = aName;
        split_trim_string( tKeyWithoutSpaces, "" );
        MORIS_ERROR( aName == tKeyWithoutSpaces,
                "Param_List::insert - key contains whitespaces" );

        // Insert new value
        Parameter tParameter( aDefaultValue, aValidSelections );
        mParamMap.insert( { aName, tParameter } );
    }

    //--------------------------------------------------------------------------------------------------------------

    void Parameter_List::insert_enum(
            const std::string&           aName,
            const Vector< std::string >& aEnumStrings )
    {
        // Check for leading and trailing whitespaces in key
        std::string tKeyWithoutSpaces = aName;
        split_trim_string( tKeyWithoutSpaces, "" );
        MORIS_ERROR( aName == tKeyWithoutSpaces,
                "Param_List::insert - key contains whitespaces" );

        // Insert new value
        Parameter tParameter( aEnumStrings );
        mParamMap.insert( { aName, tParameter } );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Parameter_List::erase( const std::string& aName )
    {
        mParamMap.erase( aName );
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

    [[nodiscard]] bool Parameter_List::exists( const std::string& aName ) const
    {
        return mParamMap.find( aName ) not_eq mParamMap.end();
    }

    //--------------------------------------------------------------------------------------------------------------

    uint Parameter_List::index( const std::string& aName )
    {
        auto tIterator = mParamMap.find( aName );

        // throw error
        MORIS_ERROR( tIterator != mParamMap.end(),
                "The requested parameter %s does not exist.\n",
                aName.c_str() );

        return tIterator->second.index();
    }
    
    //--------------------------------------------------------------------------------------------------------------
    
    const Variant& Parameter_List::get( const std::string& aName ) const
    {
        auto tIterator = mParamMap.find( aName );
        
        // throw error
        MORIS_ERROR( tIterator != mParamMap.end(),
                "The requested parameter %s does not exist.\n",
                aName.c_str() );

        return tIterator->second.get_value();
    }

    //--------------------------------------------------------------------------------------------------------------

    auto
    Parameter_List::begin() -> decltype( mParamMap.begin() )
    {
        return mParamMap.begin();
    }

    //--------------------------------------------------------------------------------------------------------------

    auto
    Parameter_List::end() -> decltype( mParamMap.end() )
    {
        return mParamMap.end();
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
    Parameter_List::size() const
    {
        return mParamMap.size();
    }

    //--------------------------------------------------------------------------------------------------------------

}
