/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_invalid_size.hpp
 *
 */

#ifndef MORIS_EXCEPTIONS_CL_INVALID_SIZE_HPP_
#define MORIS_EXCEPTIONS_CL_INVALID_SIZE_HPP_

// C++ header files.
#include <stdexcept>
#include <string>

// moris header files.
#include "core.hpp"

namespace moris
{
namespace exceptions
{
    /**
     * @brief Invalid size exception.
     *
     * @include "snippets/exceptions/cl_invalid_size.inc" // snippets EXC
     */
    class invalid_size: public std::logic_error
    {
    public:

        /**
         * moris::exceptions::invalid_size constructor.
         */
        invalid_size(
                moris::size_t const & index)
            : std::logic_error("")
            , index_m(index)
            , message_m("Invalid size " + std::to_string(index_m) + ".")
        {
        }

        /**
         * moris::exceptions::invalid_size constructor.
         */
        ~invalid_size() = default;

        /**
         * Overrides std::exception::what.
         */
        const char *
        what() const noexcept
        {
            return message_m.c_str();
        }

    private:

        /**
         * Invalid size.
         */
        moris::size_t index_m;

        /**
         * Error message.
         */
        std::string message_m;
    };

}    // namespace exceptions
}    // namespace moris

#endif /* MORIS_EXCEPTIONS_CL_INVALID_SIZE_HPP_ */

