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

#define VALIDATOR_OVERRIDES                                         \
    bool        make_valid_parameter( Variant& aVariant ) override; \
    std::string get_validation_message() override;                  \
    Validator*  copy() override;

namespace moris
{
    class Validator
    {
      public:

        /**
         * Validator constructor
         */
        Validator() = default;

        /**
         * Validator destructor
         */
        virtual ~Validator() = default;
        
        /**
         * Gets the valid selections available based on this validator.
         * 
         * @return Vector of valid selections
         */
        virtual const Vector< std::string >& get_selection_names();

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
        virtual std::string get_validation_message() = 0;

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

    // Specialize for uints to be able to also accept sints
    template<> bool Type_Validator< uint >::make_valid_parameter( Variant& aVariant );

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

    // Specialize for uint vectors to be able to also accept sints and sint vectors
    template<> bool Vector_Validator< uint >::make_valid_parameter( Variant& aVariant );

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
         * Range validator constructor, with given minimum and maximum values.
         *
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
     * @tparam T Selection validator type
     */
    template< typename T >
    class Selection_Validator : public Validator
    {
      private:
        Vector< T > mValidSelections;

      public:
        /**
         * Selection validator constructor, with set of valid selections.
         *
         * @param aValidSelections Valid values for this parameter
         */
        Selection_Validator( const Vector< T >& aValidSelections );
        
        /**
         * Gets the names of valid selections stored in the selection validator
         * 
         * @return Vector of valid selections
         */
        const Vector< std::string >& get_selection_names() override;

        VALIDATOR_OVERRIDES
    };

    template< typename T >
    Selection_Validator< T >::Selection_Validator(
            const Vector< T >& aValidSelections )
            : mValidSelections( aValidSelections )
    {
    }

    // Explicitly instantiate only the selection validator types that make sense
    template class Selection_Validator< std::string >;

    //--------------------------------------------------------------------------------------------------------------

    /**
     * An enum validator makes sure that the parameter variant stores a uint, but the set value can be
     * either a uint (via the enum itself) or a string generated by an ENUM_MACRO
     */
    class Enum_Validator : public Validator
    {
      private:
        Vector< std::string > mEnumStrings;

      public:
        /**
         * Enum validator constructor
         */
        Enum_Validator();

        /**
         * Gets the names of valid selections stored in the selection validator
         *
         * @return Vector of valid selections
         */
        const Vector< std::string >& get_selection_names() override;

        VALIDATOR_OVERRIDES
    };

    //--------------------------------------------------------------------------------------------------------------

    /**
     * Design variable validators have special rules. They allow a parameter to be set as a design variable,
     * a single real value, or a vector of 3 real values (lower bound, initial value, and upper bound).
     * This can be either a type validator or a range validator (where the range applies to all 3 values).
     */
    template<> bool Type_Validator< Design_Variable >::make_valid_parameter( Variant& aVariant );
    template<> bool Range_Validator< Design_Variable >::make_valid_parameter( Variant& aVariant );
    template<> std::string Type_Validator< Design_Variable >::get_validation_message();
    template<> std::string Range_Validator< Design_Variable >::get_validation_message();
    template class Type_Validator< Design_Variable >;
    template class Range_Validator< Design_Variable >;

    //--------------------------------------------------------------------------------------------------------------

}
