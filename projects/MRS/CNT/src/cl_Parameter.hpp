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

#include <boost/variant.hpp>
#include "moris_typedefs.hpp"
#include "cl_Validator.hpp"
#include "fn_Parsing_Tools.hpp"

namespace moris
{
    /**
     * Gets a string of the value stored inside of a variant.
     *
     * @param aVariant Input variant
     * @return Value as a string
     */
    std::string convert_variant_to_string( Variant aVariant );

    class Parameter
    {
      private:
        Variant mValue;
        Validator* mValidator;

      public:
        /**
         * Constructor for general parameter type
         *
         * @tparam T Input parameter type
         * @param aParameterValue Default value
         */
        template< typename T >
        explicit Parameter( T aParameterValue )
        {
            // Set default value without validation
            mValue = this->make_variant( aParameterValue );

            // Create validator with default
            mValidator = new Type_Validator( mValue.index() );
        }
        
        /**
         * Constructor for general parameter type
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
            mValue = this->make_variant( aParameterValue );
            
            // Create validator with default
            mValidator = new Range_Validator( mValue.index(), aMinimumValue, aMaximumValue );
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
            Variant tParameterVariant = this->make_variant( aParameterValue );

            // Validate the variant
            MORIS_ERROR( mValidator->parameter_is_valid( tParameterVariant ),
                    "Parameter %s was set incorrectly. Valid values are: %s.",
                    aParameterName.c_str(),
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
         * Takes the input value and converts it to a variant type.
         *
         * @tparam T Input parameter type
         * @param aParameterValue Input parameter
         */
        template< typename T >
        Variant make_variant( T aParameterValue )
        {
            return aParameterValue;
        }

        /**
         * Gets this parameter value as a string
         *
         * @return Parameter string
         */
        [[nodiscard]] std::string get_string() const
        {
            return convert_variant_to_string( mValue );
        }

        /**
         * Gets the underlying type index of this parameter variant.
         *
         * @return Variant index
         */
        [[nodiscard]] uint index() const
        {
            return mValue.index();
        }
    };

    //--------------------------------------------------------------------------------------------------------------

    // Declare template specializations making variants
    template<> Variant Parameter::make_variant( std::string aParameter );
    template<> Variant Parameter::make_variant( std::pair< std::string, std::string > aParameterValue );
    template<> Variant Parameter::make_variant( const char* aParameterValue );

    //--------------------------------------------------------------------------------------------------------------
}
