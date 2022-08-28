/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_tic_toc_error.hpp
 *
 */

#ifndef MORIS_EXCEPTIONS_CL_TIC_TOC_ERROR_HPP_
#define MORIS_EXCEPTIONS_CL_TIC_TOC_ERROR_HPP_

// C++ header files.
#include <stdexcept>

namespace moris
{
namespace exceptions
{
    /**
     * @brief Tic toc exception.
     *
     * @include "snippets/exceptions/cl_tic_toc_error.inc" // snippets EXC
     */
    class tic_toc_error: public std::exception
    {
    public:

        /**
         * moris::exceptions::tic_toc_error constructor.
         */
        tic_toc_error()
            : std::exception()
        {
        }

        /**
         * moris::exceptions::tic_toc_error constructor.
         */
        ~tic_toc_error() = default;

        /**
         * Overrides std::exception::what.
         */
        const char *
        what() const noexcept
        {
            return
                    "A mismatch between the moris::tic and moris::toc occurred.\n"
                    "Review your code and assert that for every tic, there is a toc.\n"
                    "In the case that a tic is inside a for loop,\n"
                    "make sure that the matching toc is also inside.\n"
                    "Remember that tic and toc are scoped function calls.\n"
                    "Do not use moris::Tic() directly, rather use moris::tic().";
        }
    };

}    // namespace exceptions
}    // namespace moris

#endif    /* MORIS_EXCEPTIONS_CL_TIC_TOC_ERROR_HPP_ */

