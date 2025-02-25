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
    template< typename Map_Type >
    class Parameter_Iterator
    {
      private:
        Map_Type                     mParameterMap;
        const Vector< std::string >& mOrderedKeys;
        luint                        mKeyIndex;

      public:
        /**
         * Parameter iterator constructor, gives all the information needed to increment this iterator.
         *
         * @param aParameterMap Map where the std::string to Parameter relationship is found
         * @param aOrderedKeys Vector of keys in the order to be iterated
         * @param aKeyIndex Key index to start iterating from
         */
        Parameter_Iterator(
                Map_Type                     aParameterMap,
                const Vector< std::string >& aOrderedKeys,
                luint                        aKeyIndex );

        /**
         * Indirection operator, used for range-based for loops.
         *
         * @return Parameter iterator by const reference
         */
        const Parameter_Iterator& operator*() const;

        /**
         * Increment operator
         *
         * @return Parameter iterator, with 1 added to the underlying key index.
         */
        Parameter_Iterator& operator++();

        /**
         * Comparison operator for checking if iterator is at the end of the parameter list container.
         *
         * @return If iterators are equal
         */
        bool operator!=( const Parameter_Iterator& ) const;

        /**
         * Gets the name of the parameter at the current key index.
         *
         * @return Parameter name
         */
        const std::string& get_name() const;

        /**
         * Gets the parameter associated with the key at the current index.
         *
         * @return Parameter currently mapped to by the key
         */
        template< bool not_const = not std::is_const< typename std::remove_reference< Map_Type >::type >::value >
        typename std::enable_if< not_const, Parameter& >::type
        get_parameter()
        {
            return mParameterMap.find( mOrderedKeys( mKeyIndex ) )->second;
        }

        /**
         * Gets a const reference to the parameter associated with the key at the current index.
         *
         * @return Parameter currently mapped to by the key
         */
        const Parameter& get_parameter() const;
    };

    // Explicit instantiation
    template class Parameter_Iterator< std::map< std::string, Parameter >& >;
    template class Parameter_Iterator< const std::map< std::string, Parameter >& >;

    class Parameter_List
    {
      private:
        std::string mName;
        std::map< std::string, Parameter > mParameterMap;
        Vector< std::string > mOrderedKeys;

        // Scoped iterator types
        typedef Parameter_Iterator< std::map< std::string, Parameter >& > iterator;
        typedef Parameter_Iterator< const std::map< std::string, Parameter >& > const_iterator;

      public:

        /**
         * Parameter list constructor
         *
         * @param aName Name of this parameter list
         */
        explicit Parameter_List( std::string aName );

        /**
         * Sets the name of this parameter list.
         *
         * @param aName Name of this collection of parameters
         */
        void set_name( std::string aName );

        /**
         * Gets the name of this parameter list.
         *
         * @return Name of this collection of parameters
         */
        const std::string& get_name();

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
                T                   aDefaultValue,
                Entry_Type          aExternalValidationType = Entry_Type::FREE,
                std::string         aExternalParameterName = "",
                Module_Type        aExternalParameterListType = Module_Type::END_ENUM,
                uint                aExternalParameterListIndex = 0 )
        {
            // Register new
            std::string tKey = this->register_key( aName );

            // Insert new value
            Parameter tParameter( std::move( aDefaultValue ), aExternalValidationType, std::move( aExternalParameterName ), aExternalParameterListType, aExternalParameterListIndex );
            mParameterMap.insert( { tKey, tParameter } );
        }

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
         * Adds a new string parameter to this parameter list. The parameter must be set from a selection of values.
         *
         * @param aName Parameter name
         * @param aDefaultValue Default parameter value
         * @param aValidSelections Valid values the parameter can be set to
         */
        void insert(
                const std::string&           aName,
                const std::string&           aDefaultValue,
                const Vector< std::string >& aValidSelections );

        /**
         * Adds a new enum parameter to this parameter list. The parameter can be set
         * as an enum directly or with one of the names provided.
         *
         * @param aName Parameter name
         * @param aEnumStrings Valid enum strings (must not be empty)
         */
        void insert_enum(
                const std::string&           aName,
                const Vector< std::string >& aEnumStrings );

        /**
         * Removes the specified key and its associated value from the map
         *
         * @param aName the key to be erased
         */
        void
        erase( const std::string& aName );

        /**
         * Sets an element to a value if it exists, otherwise an error is thrown
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
        Vector< T > get_vector( const std::string& aName ) const
        {
            return string_to_vector< T >( this->get< std::string >( aName ) );
        }

        /**
         * Gets a custom begin iterator for going through the parameter map in the order parameters were inserted.
         *
         * @return Beginning of the parameter map, by insertion
         */
        [[nodiscard]] iterator begin();

        /**
         * Gets a custom end iterator for going through the parameter map in the order parameters were inserted.
         *
         * @return End of the parameter map, by insertion
         */
        [[nodiscard]] iterator end();

        /**
         * Gets a custom begin iterator for going through the parameter map in the order parameters were inserted.
         *
         * @return Beginning of the parameter map, by insertion
         */
        [[nodiscard]] const_iterator begin() const;

        /**
         * Gets a custom end iterator for going through the parameter map in the order parameters were inserted.
         *
         * @return End of the parameter map, by insertion
         */
        [[nodiscard]] const_iterator end() const;

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
        std::string register_key( const std::string& aName );

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
            mParameterMap.insert( { tKey, tParameter } );
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
            auto tIterator = mParameterMap.find( tKey );

            // if key does not exist in map
            MORIS_ERROR( tIterator != mParameterMap.end(),
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
            auto tIterator = mParameterMap.find( aName );

            // check if key exists
            MORIS_ERROR( tIterator != mParameterMap.end(),
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
