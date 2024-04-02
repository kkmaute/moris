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

#include <variant>
#include "moris_typedefs.hpp"
#include "fn_Parsing_Tools.hpp"

namespace moris
{
    // Variant typedef
    typedef std::variant< bool, uint, sint, real, std::string, std::pair< std::string, std::string > > Variant;

    class Validator
    {
      private:
        uint mTypeIndex;

      public:

        /**
         * Validator constructor
         *
         * @param aDefaultParameter Default variant parameter, for grabbing its type_index
         */
        explicit Validator( const Variant& aDefaultParameter );

        /**
         * Validator destructor
         */
        virtual ~Validator() = default;

        /**
         * Checks if a given variant is valid, by making sure its type matches what is expected.
         *
         * @param aParameterName Parameter name (for error output)
         * @param aParameterVariant Input parameter variant
         * @return If the input parameter variant has the correct type
         */
        void check_parameter_type( const std::string& aParameterName, const Variant& aParameterVariant );

        /**
         * Checks a parameter to make sure it is valid. For a default validator, this will always return true.
         *
         * @param aParameterVariant Input parameter
         * @return If the parameter passes this validator's checks
         */
        virtual bool is_parameter_valid( const Variant& aParameterVariant );

        /**
         * Gets the valid values this validator is evaluating against.
         *
         * @return Valid values, represented with a string
         */
        virtual std::string get_valid_values();
    };
}
