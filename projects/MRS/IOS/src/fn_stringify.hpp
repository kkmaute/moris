/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_stringify.hpp
 *
 */

#ifndef PROJECTS_MRS_IOS_SRC_FN_STRINGIFY_HPP_
#define PROJECTS_MRS_IOS_SRC_FN_STRINGIFY_HPP_

#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include <limits>

#include "Log_Constants.hpp"
#include "moris_typedefs.hpp"

namespace moris::ios
{
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------

    // converts all output values to string formatted according to type
    template< typename T >
    inline std::string
    stringify( const T &aValue )
    {
        std::ostringstream out;
        out << aValue;
        return out.str();
    }

    // ----------------------------------------------------------------------------

    // bool specialization
    template<>
    inline std::string
    stringify< bool >( const bool &aValue )
    {
        std::ostringstream out;
        out << std::boolalpha << aValue;
        return out.str();
    }

    // ----------------------------------------------------------------------------

    // float specializations
    template<>
    inline std::string
    stringify< double >( const double &aValue )
    {
        std::ostringstream out;
        out << std::setprecision( LOGGER_FLOAT_PRECISION ) << std::scientific << aValue;
        return out.str();
    }

    // ----------------------------------------------------------------------------

    template<>
    inline std::string
    stringify< float >( const float &aValue )
    {
        std::ostringstream out;
        out << std::setprecision( LOGGER_FLOAT_PRECISION ) << std::scientific << aValue;
        return out.str();
    }

    // ----------------------------------------------------------------------------

    template<>
    inline std::string
    stringify< long double >( const long double &aValue )
    {
        std::ostringstream out;
        out << std::setprecision( LOGGER_FLOAT_PRECISION ) << std::scientific << aValue;
        return out.str();
    }

    // ----------------------------------------------------------------------------

    template< typename T >
    inline std::string
    stringify_cell( const std::vector< T > &aCellOfValues )
    {
        // initialize string stream
        std::ostringstream out;
        out << "{";

        // go through elements of the cell
        for ( const T &iValue : aCellOfValues )
        {
            out << " " << iValue << ";";
        }

        // convert to string and close cell notation
        std::string tOutput = out.str();
        tOutput.pop_back();
        tOutput += " }";

        // return
        return tOutput;
    }

    // ----------------------------------------------------------------------------

}    // namespace moris::ios
// end namespace moris

#endif /* PROJECTS_MRS_IOS_SRC_FN_STRINGIFY_HPP_ */
