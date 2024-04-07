/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Validator.hpp
 *
 */

#pragma once

#include "cl_Variant.hpp"

#define VALIDATOR_OVERRIDES                                                      \
    bool        parameter_is_valid( const Variant& aParameterVariant ) override; \
    std::string get_valid_values() override;                                     \
    Validator*  copy() override;

namespace moris
{
    class Validator
    {
      protected:
        uint mTypeIndex;

      public:

        /**
         * Validator constructor
         *
         * @param aTypeIndex Variant type index
         */
        explicit Validator( uint aTypeIndex );

        /**
         * Validator destructor
         */
        virtual ~Validator() = default;

        /**
         * Checks a parameter to make sure it is valid.
         *
         * @param aParameterVariant Input parameter
         * @return If the parameter passes this validator's checks
         */
        virtual bool parameter_is_valid( const Variant& aParameterVariant ) = 0;

        /**
         * Gets the valid values this validator is evaluating against.
         *
         * @return Valid values, represented with a string
         */
        virtual std::string get_valid_values() = 0;

        /**
         * Copies this validator and returns a new pointer.
         *
         * @return Validator pointer
         */
        virtual Validator* copy() = 0;

      protected:
        /**
         * Helper function that checks if a given variant has the same type index as the default parameter.
         *
         * @param aParameterVariant Input parameter variant
         * @return If the input parameter variant has the correct type
         */
        bool same_type_index( const Variant& aParameterVariant );
    };

    class Type_Validator : public Validator
    {
      public:
        /**
         * A type validator checks inputs to make sure the type is correct, but otherwise every value is valid.
         *
         * @param aTypeIndex Variant type index
         */
        explicit Type_Validator( uint aTypeIndex );

        VALIDATOR_OVERRIDES
    };

    template< typename T >
    class Range_Validator : public Validator
    {
      private:
        T mMinimumValue;
        T mMaximumValue;

      public:
        /**
         * A range validator checks inputs to ensure the value lies within a specific range.
         *
         * @param aTypeIndex Variant type index
         * @param aMinimumValue Maximum permitted parameter value
         * @param aMaximumValue Minimum permitted parameter value
         */
        Range_Validator( uint aTypeIndex, T aMinimumValue, T aMaximumValue )
                : Validator( aTypeIndex )
                , mMinimumValue( aMinimumValue )
                , mMaximumValue( aMaximumValue )
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        bool parameter_is_valid( const Variant& aParameterVariant ) override
        {
            return this->same_type_index( aParameterVariant )
               and std::get< T >( aParameterVariant ) >= mMinimumValue
               and std::get< T >( aParameterVariant ) <= mMaximumValue;
        }

        //--------------------------------------------------------------------------------------------------------------

        std::string get_valid_values() override
        {
            return get_variant_name( mTypeIndex ) + ", [" + std::to_string( mMinimumValue ) + ", " + std::to_string( mMaximumValue ) + "]";
        }

        //--------------------------------------------------------------------------------------------------------------

        Validator* copy() override
        {
            return new Range_Validator< T >( mTypeIndex, mMinimumValue, mMaximumValue );
        }

        //--------------------------------------------------------------------------------------------------------------
    };
}
