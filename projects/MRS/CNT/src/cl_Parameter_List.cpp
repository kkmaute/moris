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

    template< typename Map_Type >
    Parameter_Iterator< Map_Type >::Parameter_Iterator(
            Map_Type aParameterMap,
            const Vector< std::string >&              aOrderedKeys,
            luint                                     aKeyIndex )
            : mParameterMap( aParameterMap )
            , mOrderedKeys( aOrderedKeys )
            , mKeyIndex( aKeyIndex )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    template< typename Map_Type >
    const Parameter_Iterator< Map_Type >& Parameter_Iterator< Map_Type >::operator*() const
    {
        return *this;
    };

    //--------------------------------------------------------------------------------------------------------------

    template< typename Map_Type >
    Parameter_Iterator< Map_Type >& Parameter_Iterator< Map_Type >::operator++()
    {
        // Increment key index
        mKeyIndex++;

        // Return this
        return *this;
    }

    //--------------------------------------------------------------------------------------------------------------

    template< typename Map_Type >
    bool Parameter_Iterator< Map_Type >::operator!=( const Parameter_Iterator< Map_Type >& aComparisonIterator ) const
    {
        return mKeyIndex != aComparisonIterator.mKeyIndex;
    }

    //--------------------------------------------------------------------------------------------------------------

    template< typename Map_Type >
    const std::string& Parameter_Iterator< Map_Type >::get_name() const
    {
        return mOrderedKeys( mKeyIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    template< typename Map_Type >
    const Parameter& Parameter_Iterator< Map_Type >::get_parameter() const
    {
        return mParameterMap.find( mOrderedKeys( mKeyIndex ) )->second;
    }

    //--------------------------------------------------------------------------------------------------------------

    Parameter_List::Parameter_List( std::string aName )
            : mName( std::move( aName ) )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    void Parameter_List::set_name( std::string aName )
    {
        mName = std::move( aName );
    }

    //--------------------------------------------------------------------------------------------------------------

    const std::string& Parameter_List::get_name()
    {
        return mName;
    }

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
        mParameterMap.insert( { aName, tParameter } );
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
        mParameterMap.insert( { aName, tParameter } );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Parameter_List::erase( const std::string& aName )
    {
        mParameterMap.erase( aName );
        mOrderedKeys.remove( aName );
    }

    //--------------------------------------------------------------------------------------------------------------

    void Parameter_List::copy_parameters( const Parameter_List& aParameterList )
    {
        for ( Parameter_List::const_iterator iCopyIterator : aParameterList )
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

    Parameter_List::iterator Parameter_List::begin()
    {
        return { mParameterMap, mOrderedKeys, 0 };
    }

    //--------------------------------------------------------------------------------------------------------------

    Parameter_List::iterator Parameter_List::end()
    {
        return { mParameterMap, mOrderedKeys, mOrderedKeys.size() };
    }

    //--------------------------------------------------------------------------------------------------------------

    Parameter_List::const_iterator Parameter_List::begin() const
    {
        return { mParameterMap, mOrderedKeys, 0 };
    }

    //--------------------------------------------------------------------------------------------------------------

    Parameter_List::const_iterator Parameter_List::end() const
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
