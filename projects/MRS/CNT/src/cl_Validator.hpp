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
#include <set>
#include <utility>

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

    //--------------------------------------------------------------------------------------------------------------

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

    //--------------------------------------------------------------------------------------------------------------

    /**
     * A range validator checks inputs to ensure the value lies within a specific range.
     *
     * Note: For some reason, the constructor must be defined outside of the class declaration.
     * Currently, if this is changed, it will cause the parameter to have a segmentation fault
     * in its destructor when LTO compiler optimizations are turned on.
     *
     * @tparam T Range validator type
     */
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
        Range_Validator( uint aTypeIndex, T aMinimumValue, T aMaximumValue );

        VALIDATOR_OVERRIDES
    };

    template< typename T >
    Range_Validator< T >::Range_Validator( uint aTypeIndex, T aMinimumValue, T aMaximumValue )
                : Validator( aTypeIndex )
                , mMinimumValue( aMinimumValue )
                , mMaximumValue( aMaximumValue )
    {
    }

    // Explicitly instantiate only the range validator types that make sense
    template class Range_Validator< uint >;
    template class Range_Validator< sint >;
    template class Range_Validator< real >;

    //--------------------------------------------------------------------------------------------------------------

    /**
     * A range validator checks inputs to ensure the value lies within a specific range.
     *
     * Note: For some reason, the constructor must be defined outside of the class declaration.
     * Currently, if this is changed, it will cause the parameter to have a segmentation fault
     * in its destructor when LTO compiler optimizations are turned on.
     *
     * @tparam T Range validator type
     */
    template< typename T >
    class Selection_Validator : public Validator
    {
      private:
        std::set< T > mValidValues;

      public:
        /**
         * A selection validator checks inputs to ensure the value is one of the given selections.
         *
         * @param aTypeIndex Variant type index
         * @param aValidValues Valid values for this parameter
         */
        Selection_Validator(
                uint                 aTypeIndex,
                const std::set< T >& aValidValues );

        VALIDATOR_OVERRIDES
    };

    template< typename T >
    Selection_Validator< T >::Selection_Validator(
            uint                 aTypeIndex,
            const std::set< T >& aValidValues )
            : Validator( aTypeIndex )
            , mValidValues( aValidValues )
    {
    }

    // Explicitly instantiate only the selection validator types that make sense
    template class Selection_Validator< std::string >;

    //--------------------------------------------------------------------------------------------------------------
}
