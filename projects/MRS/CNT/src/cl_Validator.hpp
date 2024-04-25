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

#define VALIDATOR_OVERRIDES                                         \
    bool        make_valid_parameter( Variant& aVariant ) override; \
    std::string get_valid_values() override;                        \
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
        Validator() = default;

        /**
         * Validator destructor
         */
        virtual ~Validator() = default;

        /**
         * Makes a parameter into a valid 
         *
         * @param aVariant Input variant parameter
         * @return If the validator could successfully make a valid parameter
         */
        virtual bool make_valid_parameter( Variant& aVariant ) = 0;

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
    };

    /**
     * A type validator checks to make sure that the input type is correct, but otherwise every value is valid.
     *
     * @tparam T Type to check against
     */
    template< typename T >
    class Type_Validator : public Validator
    {
      public:
        VALIDATOR_OVERRIDES
    };

    // Explicitly instantiate only the type validator types that make sense
    template class Type_Validator< bool >;
    template class Type_Validator< uint >;
    template class Type_Validator< sint >;
    template class Type_Validator< real >;
    template class Type_Validator< std::string >;
    template class Type_Validator< std::pair< std::string, std::string > >;

    /**
     * A vector validator checks to make sure that the input type is either the correct vector type,
     * or that there is a single value that would fit into this vector type
     *
     * @tparam T Vector type
     */
    template< typename T >
    class Vector_Validator : public Validator
    {
      public:
        VALIDATOR_OVERRIDES
    };

    // Type validators with vectors as their template arguments are defined as vector validators
    template class Vector_Validator< uint >;
    template class Vector_Validator< real >;
    template class Vector_Validator< std::string >;
    template<> class Type_Validator< Vector< uint > > : public Vector_Validator< uint >{};
    template<> class Type_Validator< Vector< real > > : public Vector_Validator< real >{};
    template<> class Type_Validator< Vector< std::string > > : public Vector_Validator< std::string >{};

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
        Range_Validator( T aMinimumValue, T aMaximumValue );

        VALIDATOR_OVERRIDES
    };

    template< typename T >
    Range_Validator< T >::Range_Validator( T aMinimumValue, T aMaximumValue )
            : mMinimumValue( aMinimumValue )
            , mMaximumValue( aMaximumValue )
    {
    }

    // Explicitly instantiate only the range validator types that make sense
    template class Range_Validator< uint >;
    template class Range_Validator< sint >;
    template class Range_Validator< real >;

    //--------------------------------------------------------------------------------------------------------------

    /**
     * A selection validator checks inputs to ensure the value lies within a specific range.
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
        Selection_Validator( const std::set< T >& aValidValues );

        VALIDATOR_OVERRIDES
    };

    template< typename T >
    Selection_Validator< T >::Selection_Validator(
            const std::set< T >& aValidValues )
            : mValidValues( aValidValues )
    {
    }

    // Explicitly instantiate only the selection validator types that make sense
    template class Selection_Validator< std::string >;

    //--------------------------------------------------------------------------------------------------------------

    /**
     * A design variable validator has special rules. It allows a parameter to be set as a design variable,
     * a single real value, or a vector of 3 real values (lower bound, initial value, and upper bound).
     */
    class Design_Variable_Validator : public Validator
    {
      public:
        VALIDATOR_OVERRIDES
    };

    //--------------------------------------------------------------------------------------------------------------

}
