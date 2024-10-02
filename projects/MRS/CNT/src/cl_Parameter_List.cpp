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

    Parameter_Iterator::Parameter_Iterator(
            const std::map< std::string, Parameter >& aParameterMap,
            const Vector< std::string >&              aOrderedKeys,
            luint                                     aKeyIndex )
            : mParameterMap( aParameterMap )
            , mOrderedKeys( aOrderedKeys )
            , mKeyIndex( aKeyIndex )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    const Parameter_Iterator& Parameter_Iterator::operator*() const
    {
        return *this;
    };

    //--------------------------------------------------------------------------------------------------------------

    Parameter_Iterator& Parameter_Iterator::operator++()
    {
        // Increment key index
        mKeyIndex++;

        // Return this
        return *this;
    }

    //--------------------------------------------------------------------------------------------------------------

    bool Parameter_Iterator::operator!=( const Parameter_Iterator& aComparisonIterator ) const
    {
        return mKeyIndex != aComparisonIterator.mKeyIndex;
    }

    //--------------------------------------------------------------------------------------------------------------

    const std::string& Parameter_Iterator::get_name() const
    {
        return mOrderedKeys( mKeyIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    const Parameter& Parameter_Iterator::get_parameter() const
    {
        return mParameterMap.find( mOrderedKeys( mKeyIndex ) )->second;
    }

    //--------------------------------------------------------------------------------------------------------------

    void Parameter_List::insert(
            const std::string&             aName,
            const std::string&             aDefaultValue,
            const std::set< std::string >& aValidSelections )
    {
        // Check for leading and trailing whitespaces in key
        std::string tKeyWithoutSpaces = aName;
        split_trim_string( tKeyWithoutSpaces, "" );
        MORIS_ERROR( aName == tKeyWithoutSpaces,
                "Param_List::insert - key contains whitespaces" );

        // Insert new value
        Parameter tParameter( aDefaultValue, aValidSelections );
        mParameterMap.insert( { aName, tParameter } );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Parameter_List::erase( const std::string& aName )
    {
        mParameterMap.erase( aName );
    }

    //--------------------------------------------------------------------------------------------------------------

    void Parameter_List::copy_parameters( const Parameter_List& aParameterList )
    {
        for ( const Parameter_Iterator& iCopyIterator : aParameterList )
        {
            // Get iterator from this map
            auto tFoundIterator = mParameterMap.find( iCopyIterator.get_name() );

            // Determine if parameter needs to be insert or set
            if ( tFoundIterator == mParameterMap.end() )
            {
                // Insert parameter
                mParameterMap.insert( { iCopyIterator.get_name(), iCopyIterator.get_parameter() } );
            }
            else
            {
                // Set parameter
                tFoundIterator->second = iCopyIterator.get_parameter();
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    [[nodiscard]] bool Parameter_List::exists( const std::string& aName ) const
    {
        return mParameterMap.find( aName ) not_eq mParameterMap.end();
    }

    //--------------------------------------------------------------------------------------------------------------

    uint Parameter_List::index( const std::string& aName )
    {
        auto tIterator = mParameterMap.find( aName );

        // throw error
        MORIS_ERROR( tIterator != mParameterMap.end(),
                "The requested parameter %s does not exist.\n",
                aName.c_str() );

        return tIterator->second.index();
    }
    
    //--------------------------------------------------------------------------------------------------------------
    
    const Variant& Parameter_List::get_variant( const std::string& aName ) const
    {
        auto tIterator = mParameterMap.find( aName );
        
        // throw error
        MORIS_ERROR( tIterator != mParameterMap.end(),
                "The requested parameter %s does not exist.\n",
                aName.c_str() );

        return tIterator->second.get_value();
    }

    //--------------------------------------------------------------------------------------------------------------

    Parameter_Iterator Parameter_List::begin() const
    {
        return { mParameterMap, mOrderedKeys, 0 };
    }

    //--------------------------------------------------------------------------------------------------------------

    Parameter_Iterator Parameter_List::end() const
    {
        return { mParameterMap, mOrderedKeys, mOrderedKeys.size() };
    }

    //--------------------------------------------------------------------------------------------------------------

    bool
    Parameter_List::is_empty()
    {
        return mParameterMap.empty();
    }

    //--------------------------------------------------------------------------------------------------------------

    size_t
    Parameter_List::size() const
    {
        return mParameterMap.size();
    }

    //--------------------------------------------------------------------------------------------------------------

    std::string Parameter_List::register_key( const std::string& aName )
    {
        // Trim leading and trailing whitespaces from name to form key
        std::string tKey = aName;
        trim_string( tKey );

        // Add key to ordered list
        mOrderedKeys.push_back( tKey );

        // Return key
        return tKey;
    }

    //--------------------------------------------------------------------------------------------------------------

}
