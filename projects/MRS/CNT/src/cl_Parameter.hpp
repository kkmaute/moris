/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Parameter.hpp
 *
 */

#pragma once

#include <utility>

#include "moris_typedefs.hpp"
#include "cl_Validator.hpp"
#include "cl_Library_Enums.hpp"
#include "fn_Parsing_Tools.hpp"

namespace moris
{
    enum class Entry_Type
    {
        FREE,
        SELECTION,
        FIXED_SIZE_VECTOR
    };

    /**
     * For validation after all parameters have been set.
     */
    struct External_Validator
    {
        std::string         mParameterName;
        Parameter_List_Type mParameterListType  = Parameter_List_Type::END_ENUM;
        uint                mParameterListIndex = 0;
    };

    class Parameter
    {
      private:
        Variant mValue;
        Entry_Type mEntryType = Entry_Type::FREE;
        uint mNumberOfEntries = 1;
        Validator* mValidator;
        External_Validator mExternalValidator;

      public:
        /**
         * Constructor for a parameter with a type validator.
         *
         * @tparam T Input parameter type
         * @param aParameterValue Default value
         * @param aExternalValidationType Type of external validation to perform
         * @param aExternalParameterName Name of external parameter to validate with
         * @param aExternalParameterListType External parameter list type to validate with
         * @param aExternalParameterListIndex Index of given parameter list type to search
         */
        template< typename T >
        Parameter(
                T                   aParameterValue,
                Entry_Type          aExternalValidationType,
                std::string         aExternalParameterListName,
                Parameter_List_Type aExternalParameterListType,
                uint                aExternalParameterListIndex )
        {
            // Set default value without validation
            mValue = make_variant( aParameterValue );

            // Create type validator
            mValidator = new Type_Validator< T >();

            // Set type
            mEntryType = aExternalValidationType;

            // Set external validator
            mExternalValidator.mParameterName = std::move( aExternalParameterListName );
            mExternalValidator.mParameterListType = aExternalParameterListType;
            mExternalValidator.mParameterListIndex = aExternalParameterListIndex;
        }

        /**
         * Constructor for a vector parameter with a type validator.
         *
         * @tparam T Input parameter type
         * @param aParameterValue Default value
         * @param aExternalValidationType Type of external validation to perform
         * @param aExternalParameterName Name of external parameter to validate with
         * @param aExternalParameterListType External parameter list type to validate with
         * @param aExternalParameterListIndex Index of given parameter list type to search
         */
        template< typename T >
        Parameter(
                Vector< T >         aParameterValue,
                Entry_Type          aExternalValidationType,
                std::string         aExternalParameterListName,
                Parameter_List_Type aExternalParameterListType,
                uint                aExternalParameterListIndex )
                : mNumberOfEntries( 0 )
        {
            // Set default value without validation
            mValue = make_variant( aParameterValue );

            // Create type validator
            mValidator = new Vector_Validator< T >();

            // Set entry type
            mEntryType = aExternalValidationType;

            // Set external validator
            mExternalValidator.mParameterName = std::move( aExternalParameterListName );
            mExternalValidator.mParameterListType = aExternalParameterListType;
            mExternalValidator.mParameterListIndex = aExternalParameterListIndex;
        }
        
        /**
         * Constructor for a parameter with a range validator.
         *
         * @tparam T Input parameter type
         * @param aParameterValue Default value
         * @param aMinimumValue Maximum permitted parameter value
         * @param aMaximumValue Minimum permitted parameter value
         */
        template< typename T >
        Parameter( T aParameterValue, T aMinimumValue, T aMaximumValue )
        {
            // Set default value without validation
            mValue = make_variant( aParameterValue );
            
            // Create range validator
            mValidator = new Range_Validator( aMinimumValue, aMaximumValue );
        }

        /**
         * Constructor for a parameter with a selection validator.
         *
         * @tparam T Input parameter type
         * @param aParameterValue Default value
         * @param aValidSelections Set of valid values
         */
        template< typename T >
        Parameter( T aParameterValue, const std::set< T >& aValidSelections )
                : mEntryType( Entry_Type::SELECTION )
        {
            // Set default value without validation
            mValue = make_variant( aParameterValue );

            // Create selection validator
            mValidator = new Selection_Validator( aValidSelections );
        }

        /**
         * Custom copy constructor, used to ensure each parameter has a validator.
         *
         * @param aParameter Parameter to copy
         */
        Parameter( const Parameter& aParameter );

        /**
         * Parameter destructor, deletes the validator.
         */
        ~Parameter();

        /**
         * Sets the value of this parameter
         *
         * @tparam T Input parameter type
         * @param aParameterName Parameter name, for error reporting
         * @param aParameterValue Input value
         * @param aLockValue If this value should be locked after setting
         */
        template< typename T >
        void set_value(
                const std::string& aParameterName,
                T                  aParameterValue,
                bool               aLockValue = true )
        {
            // Make sure parameter is not locked
            MORIS_ERROR( mValidator,
                    "Parameter %s has already been set and locked, it cannot be set again.",
                    aParameterName.c_str() );

            // Make value into a variant
            Variant tParameterVariant = make_variant( aParameterValue );

            // Validate the variant
            MORIS_ERROR( mValidator->make_valid_parameter( tParameterVariant ),
                    "Parameter %s was set incorrectly as %s. Valid values are: %s.",
                    aParameterName.c_str(),
                    convert_variant_to_string( tParameterVariant ).c_str(),
                    mValidator->get_valid_values().c_str() );

            // Set the value
            mValue = tParameterVariant;

            // Lock the parameter by deleting the validator
            if ( aLockValue )
            {
                delete mValidator;
                mValidator = nullptr;
            }
        }

        /**
         * Gets the value of this parameter.
         *
         * @tparam T Requested type
         * @return Stored value
         */
        template< typename T >
        const T& get_value() const
        {
            return std::get< T >( mValue );
        }

        /**
         * Gets the value of this parameter, as a variant.
         *
         * @return Stored variant
         */
        const Variant& get_value() const;

        /**
         * Gets this parameter value as a string
         *
         * @return Parameter string
         */
        [[nodiscard]] std::string get_string() const;

        /**
         * Gets the underlying type index of this parameter variant.
         *
         * @return Variant index
         */
        [[nodiscard]] uint index() const;

        /**
         * Gets the type of entry for this parameter.
         *
         * @return Parameter type
         */
        Entry_Type get_entry_type() const;

        /**
         * Gets the number of entries this parameter requires.
         *
         * @return Number of entries, or 0 if resizable
         */
        uint get_number_of_entries() const;

        /**
         * Gets the external validator from this parameter, for validation after all parameter have been set.
         *
         * @return External validator
         */
        const External_Validator& get_external_validator() const;

        /**
         * Equality operator for parameters, comparing if their stored variants are the same.
         *
         * @param aOther Other parameter argument
         * @return If parameters are equal
         */
        bool operator==( const Parameter& aOther );
    };

    //--------------------------------------------------------------------------------------------------------------

    // Declare template specializations of the Parameter constructor
    template<> Parameter::Parameter( const char*, Entry_Type, std::string, Parameter_List_Type, uint );
}
