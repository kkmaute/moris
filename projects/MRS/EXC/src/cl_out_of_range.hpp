/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_out_of_range.hpp
 *
 */

#ifndef MORIS_EXCEPTIONS_CL_OUT_OF_RANGE_HPP_
#define MORIS_EXCEPTIONS_CL_OUT_OF_RANGE_HPP_

// C++ header files.
#include <stdexcept>
#include <string>

// moris header files.
#include "core.hpp"

#include <iostream>

namespace moris::exceptions
{
    /**
     * @brief Out-of-range exception.
     *
     * @include "snippets/exceptions/cl_out_of_range.inc" // snippets EXC
     */
    class out_of_range: public std::out_of_range
    {
    public:

        /**
         * moris::exceptions::out_of_range constructor.
         *
         * @note[chvillanuevap@gmail.com]
         * Does std::string throw an exception on construction?
         * If so, using a std::string is dangerous here.
         * Don't use any construct that might throw an exception in an exception.
         * @see [(not) using std::string in exceptions]
         * (http://stackoverflow.com/questions/15831029/not-using-stdstring-in-exceptions)
         */
        out_of_range(
                moris::size_t const & index)
            : std::out_of_range("")
            , index_m(index)
            , message_m("Index " + std::to_string(index_m) + " out of bounds.")
        {
        }

        /**
         * moris::exceptions::out_of_range constructor.
         */
        ~out_of_range() override = default;

        /**
         * Overrides std::exception::what.
         */
        const char *
        what() const noexcept override
        {
            return message_m.c_str();
        }

    private:

        /**
         * Out-of-range index.
         */
        moris::size_t index_m;

        /**
         * Error message.
         */
        std::string message_m;
    };

}    // namespace moris::exceptions

#endif /* MORIS_EXCEPTIONS_CL_OUT_OF_RANGE_HPP_ */

