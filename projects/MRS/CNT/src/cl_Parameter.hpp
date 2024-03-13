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
         * Constructor for general type
         *
         * @tparam T Input parameter type
         * @param aParameterName Parameter name, for error reporting
         * @param aParameterValue Default value
         */
        template< typename T >
        explicit Parameter( T aParameterValue )
        {
            // Set default value without validation
            mValue = this->make_variant( aParameterValue );

            // Create validator with default
            mValidator = new Validator( mValue );
        }

        /**
         * Parameter destructor, deletes the validator.
         */
        ~Parameter()
        {
            delete mValidator;
        }

        /**
         * Sets the value of this parameter
         *
         * @tparam T Input parameter type
         * @param aParameterName Parameter name, for error reporting
         * @param aParameterValue Input value
         */
        template< typename T >
        void set_value( const std::string& aParameterName, T aParameterValue )
        {
            // Make value into a variant
            Variant tParameterVariant = this->make_variant( aParameterValue );

            // Validate the variant
            mValidator->check_parameter_type( aParameterName, tParameterVariant );
            MORIS_ERROR( mValidator->is_parameter_valid( tParameterVariant ),
                    "Parameter %s was set with an invalid value. Valid values are: %s.",
                    aParameterName.c_str(),
                    mValidator->get_valid_values().c_str() );

            // Set the value
            mValue = tParameterVariant;
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
            return boost::get< T >( mValue );
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
        [[nodiscard]] sint which() const
        {
            return mValue.which();
        }
    };

    //--------------------------------------------------------------------------------------------------------------

    // Declare template specializations making variants
    template<> Variant Parameter::make_variant( std::string aParameter );
    template<> Variant Parameter::make_variant( std::pair< std::string, std::string > aParameterValue );
    template<> Variant Parameter::make_variant( const char* aParameterValue );

    //--------------------------------------------------------------------------------------------------------------
}
