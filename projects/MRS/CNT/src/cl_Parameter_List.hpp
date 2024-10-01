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
#include <utility>

// MORIS library header files.
#include "moris_typedefs.hpp"
#include "assert.hpp"
#include "core.hpp"
#include "ios.hpp"
#include "cl_Parameter.hpp"

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
        std::map< std::string, Parameter > mParamMap;

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
         * Adds a new parameter to this parameter list. This parameter may be validated against other parameters.
         * If no parameter list type is given, then only the parameter type will be validated upon being set.
         *
         * @tparam T Parameter type
         * @param aName Parameter name
         * @param aDefaultValue Default parameter value
         * @param aExternalValidationType Type of external validation to perform
         * @param aExternalParameterName Name of external parameter to validate with
         * @param aExternalParameterListType External parameter list type to validate with
         * @param aExternalParameterListIndex Index of given parameter list type to search
         */
        template< typename T >
        void insert(
                const std::string&  aName,
                const T&            aDefaultValue,
                Validation_Type     aExternalValidationType     = Validation_Type::NONE,
                std::string         aExternalParameterName      = "",
                Parameter_List_Type aExternalParameterListType  = Parameter_List_Type::END_ENUM,
                uint                aExternalParameterListIndex = 0 )
        {
            // Register new
            std::string tKey = this->register_key( aName );

            // Insert new value
            Parameter tParameter( std::move( aDefaultValue ), aExternalValidationType, std::move( aExternalParameterName ), aExternalParameterListType, aExternalParameterListIndex );
            mParamMap.insert( { tKey, tParameter } );
        }

        /**
         * Adds a new string parameter to this parameter list. The parameter must be set from a selection of values.
         *
         * @param aName Parameter name
         * @param aDefaultValue Default parameter value
         * @param aValidSelections Valid values the parameter can be set to
         */
        void insert(
                const std::string&             aName,
                const std::string&             aDefaultValue,
                const std::set< std::string >& aValidSelections );

        /**
         * Adds a new parameter to this parameter list. This parameter must be set with a valid range.
         *
         * @tparam T Parameter type
         * @param aName Parameter name
         * @param aDefaultValue Default parameter value
         * @param aMinimumValue Maximum permitted parameter value
         * @param aMaximumValue Minimum permitted parameter value
         */
        template< typename T >
        void insert(
                const std::string& aName,
                T                  aDefaultValue,
                T                  aMinimumValue,
                T                  aMaximumValue )
        {
            // Delegate to private implementation overload, depending on if T is an enum
            this->convert_and_insert( aName, aDefaultValue, aMinimumValue, aMaximumValue, std::is_enum< T >() );
        }

        /**
         * Removes the specified key and its associated value from the map
         *
         * @param aName the key to be erased
         */
        void
        erase( const std::string& aName );

        /**
         * Sets an element to a value if it exists, otherwise an error is thrown
         *        whitespaces in key and value string will be removed
         *
         * @param aName Parameter name
         * @param aValue Parameter value
         * @param aLockValue If the set value is to be locked, and unable to be set again.
         */
        template< typename T >
        void set(
                const std::string& aName,
                const T&           aValue,
                bool               aLockValue = true )
        {
            // Delegate to private implementation overload, depending on if T is an enum
            this->convert_and_set( aName, std::move( aValue ), aLockValue, std::is_enum< T >() );
        }

        /**
         * Sets an element to a moris vector using a parameter pack
         *
         * @param aName Parameter name
         * @param aFirstValue First value to put into the vector
         * @param aSecondValue Second value to put into the vector
         * @param aMoreValues Parameter pack of more values
         */
        template< typename T, typename... Arg_Types >
        void set(
                const std::string& aName,
                T                  aFirstValue,
                T                  aSecondValue,
                Arg_Types...       aMoreValues )
        {
            // Delegate to private implementation overload, with lock on
            this->convert_and_set( aName, Vector< T >( { aFirstValue, aSecondValue, aMoreValues... } ), true, std::false_type() );
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
         * @param aName Key in the map
         * @return whether or not this key is in the map with a corresponding value
         */
        [[nodiscard]] bool exists( const std::string& aName ) const;

        /**
         * @brief Determine the underlying type of the entry.
         *
         * @param[in] aName Key corresponding to the mapped value that
         *            needs to be accessed
         *
         * @return The index of the type as ordered in the variant entry
         *         of mParamMap
         */
        uint index( const std::string& aName );

        /**
         * Gets the parameter value corresponding to a given key.
         *
         * @param[in] aName Name of the parameter to get
         * @return The value corresponding to aName
         */
        template< typename T >
        const T& get( const std::string& aName ) const
        {
            // Delegate to private implementation overload, depending on if T is an enum
            return this->get_and_convert< T >( aName, std::is_enum< T >() );
        }

        /**
         * Gets the parameter value corresponding to a given name, as a variant.
         *
         * @param aName Name of the parameter to get
         * @return The variant corresponding to aName
         */
        const Variant& get_variant( const std::string& aName ) const;

        /**
         * Gets a cell from a paramter that is stored as a std::string
         *
         * @tparam T Cell data type
         * @param aName Key corresponding to mapped string
         * @return Cell of values
         */
        template< typename T >
        Vector< T > get_cell( const std::string& aName ) const
        {
            return string_to_vector< T >( this->get< std::string >( aName ) );
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
        size_t size() const;

      private:

        /**
         * Registers a new key to this parameter list's vector of names, based on the parameter name.
         *
         * @param aName Parameter name
         * @return Valid map key
         */
        std::string register_key( const std::string& aName )
        {
            // Trim leading and trailing whitespaces from name to form key
            std::string tKey = aName;
            trim_string( tKey );

            // Add key to list TODO

            // Return key
            return tKey;
        }

        /**
         * Adds a new parameter to this parameter list (private implementation).
         *
         * @param[in] aName Key corresponding to the mapped value that
         *            needs to be accessed; should not contain trailing or leading whitespaces
         * @param[in] aValue Value corresponding to aName; if value is a string,
         *            pair of strings, or character array leading and trailing whitespaces
         *            are trimmed off
         */
        template< typename T >
        void convert_and_insert( const std::string& aName, T aValue, T aMinimumValue, T aMaximumValue, std::false_type )
        {
            // Register new key
            std::string tKey = this->register_key( aName );

            // Insert new value
            Parameter tParameter( aValue, aMinimumValue, aMaximumValue );
            mParamMap.insert( { tKey, tParameter } );
        }

        /**
         * Insert function overload, for static casting enums to a uint.
         */
        template< typename T >
        void convert_and_insert( const std::string& aName, T aValue, T aMinimumValue, T aMaximumValue, std::true_type )
        {
            this->convert_and_insert(
                    aName,
                    static_cast< uint >( aValue ),
                    static_cast< uint >( aMinimumValue ),
                    static_cast< uint >( aMaximumValue ),
                    std::false_type() );
        }

        /**
         * Sets an element to a value if it exists (private implementation)
         *
         * @param aName Parameter name
         * @param aValue Parameter value
         * @param aLockValue If the set value is to be locked, and unable to be set again.
         */
        template< typename T >
        void convert_and_set(
                const std::string& aName,
                const T&           aValue,
                bool               aLockValue,
                std::false_type )
        {
            // Remove leading and trailing whitespaces in key
            std::string tKey = aName;
            trim_string( tKey );

            // find key in map
            auto tIterator = mParamMap.find( tKey );

            // if key does not exist in map
            MORIS_ERROR( tIterator != mParamMap.end(),
                    "The requested parameter %s can not be set because it does not exist.\n",
                    tKey.c_str() );

            tIterator->second.set_value( tKey, std::move( aValue ), aLockValue );
        }

        /**
         * Set function overload, for static casting enums to a uint.
         */
        template< typename T >
        void convert_and_set(
                const std::string& aName,
                T                  aValue,
                bool               aLockValue,
                std::true_type )
        {
            this->convert_and_set( aName, static_cast< uint >( aValue ), aLockValue, std::false_type() );
        }

        /**
         * Gets the parameter value corresponding to a given key (private implementation)
         *
         * @param[in] aName Key corresponding to the mapped value that
         *            needs to be accessed
         *
         * @return The value corresponding to aName.
         */
        template< typename T >
        const T& get_and_convert( const std::string& aName, std::false_type ) const
        {
            auto tIterator = mParamMap.find( aName );

            // check if key exists
            MORIS_ERROR( tIterator != mParamMap.end(),
                    "The requested parameter %s does not exist.\n",
                    aName.c_str() );

            return tIterator->second.get_value< T >();
        }

        /**
         * Get function overload, for static casting uints to a requested enum.
         */
        template< typename T >
        const T& get_and_convert( const std::string& aName, std::true_type ) const
        {
            return ( const T& ) this->get_and_convert< uint >( aName, std::false_type() );
        }
    };

    //------------------------------------------------------------------------------

}
