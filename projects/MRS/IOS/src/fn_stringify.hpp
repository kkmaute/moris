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
#include "typedefs.hpp"

namespace moris
{
    namespace ios
    {
        // converts all output values to string formated according to type
        template<typename T>
        inline std::string stringify(T aValue)
        {
            std::ostringstream out;
            out << aValue;
            return out.str();
        }

        template<>
        inline std::string stringify<bool>(bool aValue)
        {
            std::ostringstream out;
            out << std::boolalpha << aValue;
            return out.str();
        }

        template<>
        inline std::string stringify<double>(double aValue)
        {
            std::ostringstream out;
            out << std::setprecision(LOGGER_FLOAT_PRECISION) << std::scientific << aValue;
            return out.str();
        }

        template<>
        inline std::string stringify<long double>(long double aValue)
        {
            std::ostringstream out;
            out << std::setprecision(LOGGER_FLOAT_PRECISION) << std::scientific << aValue;
            return out.str();
        }

        template<>
        inline std::string stringify<float>(float aValue)
        {
            std::ostringstream out;
            out << std::setprecision(LOGGER_FLOAT_PRECISION) << std::scientific << aValue;
            return out.str();
        }

    } // end namespace ios
} // end namespace moris

#endif /* PROJECTS_MRS_IOS_SRC_FN_STRINGIFY_HPP_ */

