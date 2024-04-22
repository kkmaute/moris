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

#include "moris_typedefs.hpp"
#include "cl_Validator.hpp"
#include "fn_Parsing_Tools.hpp"

namespace moris
{
    class Parameter
    {
      private:
        Variant mValue;
        Validator* mValidator;

      public:
        /**
         * Constructor for a parameter with a type validator.
         *
         * @tparam T Input parameter type
         * @param aParameterValue Default value
         */
        template< typename T >
        explicit Parameter( T aParameterValue )
        {
            // Set default value without validation
            mValue = make_variant( aParameterValue );

            // Create type validator
            mValidator = new Type_Validator( mValue.index() );
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
            mValidator = new Range_Validator( mValue.index(), aMinimumValue, aMaximumValue );
        }

        /**
         * Constructor for a parameter with a selection validator.
         *
         * @tparam T Input parameter type
         * @param aParameterValue Default value
         * @param aValidValues Set of valid values
         */
        template< typename T >
        Parameter( T aParameterValue, const std::set< T >& aValidValues )
        {
            // Set default value without validation
            mValue = make_variant( aParameterValue );

            // Create selection validator
            mValidator = new Selection_Validator( mValue.index(), aValidValues );
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
            Variant tParameterVariant = this->make_valid_variant( aParameterValue );

            // Validate the variant
            MORIS_ERROR( mValidator->parameter_is_valid( tParameterVariant ),
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
         * Makes the given parameter into a variant to be checked
         *
         * @tparam T Input parameter type
         * @param aParameterValue Parameter value
         * @return Parameter variant
         */
        template< typename T >
        Variant make_valid_variant( T aParameterValue )
        {
            return make_variant( aParameterValue );
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
    };

    //--------------------------------------------------------------------------------------------------------------

    // Declare template specializations
    template<> Variant Parameter::make_valid_variant( uint );
    template<> Variant Parameter::make_valid_variant( real );
}
