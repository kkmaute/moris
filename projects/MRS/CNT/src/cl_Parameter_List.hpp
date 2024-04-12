/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Parameter_List.hpp
 *
 */

#pragma once

// c++ header files.
#include <cstring>
#include <map>
#include <string>
#include <sstream>

// MORIS library header files.
#include "moris_typedefs.hpp"
#include "assert.hpp"
#include "core.hpp"
#include "ios.hpp"
#include "cl_Parameter.hpp"

namespace moris::containers
{
    struct strcmp
    {
        bool
        operator()( const std::string& a, const std::string& b ) const
        {
            return std::strcmp( a.c_str(), b.c_str() ) < 0;
        }
    };
}

namespace moris
{
    /**
     * The parameter list class.
     *
     * A parameter list can be created as follows:
     * @include CON/src/cl_param_list/cl_param_list_insert.inc
     */
    class Parameter_List
    {

      private:
        std::map< std::string, Parameter, moris::containers::strcmp > mParamMap;

      public:

        /**
         * Constructor
         */
        Parameter_List() = default;

        /**
         * Destructor
         */
        ~Parameter_List() = default;

        /**
         * Adds a new parameter to this parameter list. This parameter's type will be validated when being set.
         *
         * @param[in] aKey Key corresponding to the mapped value that
         *            needs to be accessed; should not contain trailing or leading whitespaces
         * @param[in] aValue Value corresponding to aKey; if value is a string,
         *            pair of strings, or character array leading and trailing whitespaces
         *            are trimmed off
         */
        template< typename T >
        void insert( const std::string& aKey, T aValue )
        {
            // Check for leading and trailing whitespaces in key
            std::string tKeyWithoutSpaces = aKey;
            split_trim_string( tKeyWithoutSpaces, "" );
            MORIS_ERROR( aKey == tKeyWithoutSpaces,
                    "Param_List::insert - key contains whitespaces" );

            // Insert new value
            Parameter tParameter( aValue );
            mParamMap.insert( { aKey, tParameter } );
        }

        /**
         * Adds a new parameter to this parameter list. This parameter must be set with a valid range.
         *
         * @param[in] aKey Key corresponding to the mapped value that
         *            needs to be accessed; should not contain trailing or leading whitespaces
         * @param[in] aValue Value corresponding to aKey; if value is a string,
         *            pair of strings, or character array leading and trailing whitespaces
         *            are trimmed off
         * @param aMinimumValue Maximum permitted parameter value
         * @param aMaximumValue Minimum permitted parameter value
         */
        template< typename T >
        void insert( const std::string& aKey, T aValue, T aMinimumValue, T aMaximumValue )
        {
            // Delegate to private implementation overload, depending on if T is an enum
            this->insert( aKey, aValue, aMinimumValue, aMaximumValue, std::is_enum< T >() );
        }

        /**
         * Removes the specified key and its associated value from the map
         *
         * @param aKey the key to be erased
         */
        void
        erase( const std::string& aKey );

        /**
         * Sets an element to a value if it exists, otherwise an error is thrown
         *        whitespaces in key and value string will be removed
         *
         * @param[in] aKey Key corresponding to the mapped value that
         *            needs to be accessed; spaces will be removed
         * @param[in] aValue Value corresponding to aKey; if values is a string, pair of
         *            strings or character array, whitespaces will be removed
         * @param aLockValue If the set value is to be locked, and unable to be set again.
         */
        template< typename T >
        void set( const std::string& aKey, T aValue, bool aLockValue = true )
        {
            // Delegate to private implementation overload, depending on if T is an enum
            this->set( aKey, aValue, aLockValue, std::is_enum< T >() );
        }

        /**
         * Sets an element to a moris vector using a std::initializer_list
         *
         * @param[in] aKey Key corresponding to the mapped value that
         *            needs to be accessed; spaces will be removed
         * @param[in] aValue Value corresponding to aKey; if values is a string, pair of
         *            strings or character array, whitespaces will be removed
         * @param aLockValue If the set value is to be locked, and unable to be set again.
         */
        template< typename T >
        void set( const std::string& aKey, std::initializer_list< T > aValue, bool aLockValue = true )
        {
            // Delegate to private implementation overload, T is guaranteed not an enum
            this->set( aKey, Vector< T >( aValue ), aLockValue, std::false_type() );
        }

        /**
         * Copies the parameters over from a different parameter list.
         * If the parameter exists already, the parameter is set to the copied value.
         * If the parameter does not exist, the parameter is inserted into this parameter list.
         *
         * @param aParameterList Parameter list to copy values over from.
         */
        void copy_parameters( const Parameter_List& aParameterList );

        /**
         * Checks whether or not the given key exists. Useful for when a parameter may or may not have been inserted.
         *
         * @param aKey Key in the map
         * @return whether or not this key is in the map with a corresponding value
         */
        [[nodiscard]] bool exists( const std::string& aKey ) const;

        /**
         * @brief Determine the underlying type of the entry.
         *
         * @param[in] aKey Key corresponding to the mapped value that
         *            needs to be accessed
         *
         * @return The index of the type as ordered in the variant entry
         *         of mParamMap
         */
        uint index( const std::string& aKey );

        /**
         * Gets the parameter value corresponding to a given key.
         *
         * @param[in] aKey Key corresponding to the mapped value that
         *            needs to be accessed
         *
         * @return The value corresponding to aKey.
         */
        template< typename T >
        const T& get( const std::string& aKey ) const
        {
            // Delegate to private implementation overload, depending on if T is an enum
            return this->get< T >( aKey, std::is_enum< T >() );
        }

        /**
         * Gets a cell from a paramter that is stored as a std::string
         *
         * @tparam T Cell data type
         * @param aKey Key corresponding to mapped string
         * @return Cell of values
         */
        template< typename T >
        Vector< T > get_cell( const std::string& aKey ) const
        {
            return string_to_cell< T >( this->get< std::string >( aKey ) );
        }

        /**
         * @brief Get the beginning iterator of the underlying map.
         *
         * @return An iterator pointing to the beginning of mParamMap.
         */
        auto
        begin() const -> decltype( mParamMap.begin() );

        /**
         * @brief Get the end iterator of the underlying map.
         *
         * @return An iterator pointing to the end of mParamMap.
         */
        auto
        end() const -> decltype( mParamMap.end() );

        /**
         * @brief check if parameter list is empty
         *
         * @return true if parameter list doesn't have any entries
         */
        bool is_empty();

        /**
         * @brief get the size of the parameter list
         *
         * @return size_t number of entries in the parameter list
         */
        size_t size();

      private:

        /**
         * Adds a new parameter to this parameter list (private implementation).
         *
         * @param[in] aKey Key corresponding to the mapped value that
         *            needs to be accessed; should not contain trailing or leading whitespaces
         * @param[in] aValue Value corresponding to aKey; if value is a string,
         *            pair of strings, or character array leading and trailing whitespaces
         *            are trimmed off
         */
        template< typename T >
        void insert( const std::string& aKey, T aValue, T aMinimumValue, T aMaximumValue, std::false_type )
        {
            // Check for leading and trailing whitespaces in key
            std::string tKeyWithoutSpaces = aKey;
            split_trim_string( tKeyWithoutSpaces, "" );
            MORIS_ERROR( aKey == tKeyWithoutSpaces,
                    "Param_List::insert - key contains whitespaces" );

            // Insert new value
            Parameter tParameter( aValue, aMinimumValue, aMaximumValue );
            mParamMap.insert( { aKey, tParameter } );
        }

        /**
         * Insert function overload, for static casting enums to a uint.
         */
        template< typename T >
        void insert( const std::string& aKey, T aValue, T aMinimumValue, T aMaximumValue, std::true_type )
        {
            this->insert(
                    aKey,
                    static_cast< uint >( aValue ),
                    static_cast< uint >( aMinimumValue ),
                    static_cast< uint >( aMaximumValue ),
                    std::false_type() );
        }

        /**
         * Sets an element to a value if it exists (private implementation)
         *
         * @param[in] aKey Key corresponding to the mapped value that
         *            needs to be accessed; spaces will be removed
         * @param[in] aValue Value corresponding to aKey; if values is a string, pair of
         *            strings or character array, whitespaces will be removed
         * @param aLockValue If the set value is to be locked, and unable to be set again.
         */
        template< typename T >
        void set( const std::string& aKey, T aValue, bool aLockValue, std::false_type )
        {
            // create copy of key string such that tIterator can be manipulated
            std::string tKey = aKey;

            // remove spurious whitespaces from key string
            split_trim_string( tKey, "" );

            // find key in map
            auto tIterator = mParamMap.find( tKey );

            // if key does not exist in map
            MORIS_ERROR( tIterator != mParamMap.end(),
                    "The requested parameter %s can not be set because it does not exist.\n",
                    tKey.c_str() );

            tIterator->second.set_value( aKey, aValue, aLockValue );
        }

        /**
         * Set function overload, for static casting enums to a uint.
         */
        template< typename T >
        void set( const std::string& aKey, T aValue, bool aLockValue, std::true_type )
        {
            this->set( aKey, static_cast< uint >( aValue ), aLockValue, std::false_type() );
        }

        /**
         * Gets the parameter value corresponding to a given key (private implementation)
         *
         * @param[in] aKey Key corresponding to the mapped value that
         *            needs to be accessed
         *
         * @return The value corresponding to aKey.
         */
        template< typename T >
        const T& get( const std::string& aKey, std::false_type ) const
        {
            auto tIterator = mParamMap.find( aKey );

            // check if key exists
            MORIS_ERROR( tIterator != mParamMap.end(),
                    "The requested parameter %s does not exist.\n",
                    aKey.c_str() );

            return tIterator->second.get_value< T >();
        }

        /**
         * Get function overload, for static casting uints to a requested enum.
         */
        template< typename T >
        const T& get( const std::string& aKey, std::true_type ) const
        {
            return ( const T& ) this->get< uint >( aKey, std::false_type() );
        }
    };

    //------------------------------------------------------------------------------

}
