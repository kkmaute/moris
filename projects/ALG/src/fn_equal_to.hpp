/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_equal_to.hpp
 *
 */

#ifndef MORIS_ALGORITHMS_FN_EQUAL_TO_HPP_
#define MORIS_ALGORITHMS_FN_EQUAL_TO_HPP_

// C++ header files.
#include <algorithm>

// MORIS header files.
#include "typedefs.hpp" // COR/src

namespace moris
{
/**
 * Floating point equality comparison. For further discussion, see this StackOverflow
 * <a href="http://www.stackoverflow.com/questions/4010240/comparing-doubles">
 *    Question</a>.
 *
 * @note Modeled after std::equal_to.
 *
 * @param[in] a            Floating point number.
 * @param[in] b            Floating point number.
 * @param[in] error_factor Number of machine epsilon differences
 *                         between a and b.
 *
 * @return Boolean true if equal.
 */
template< class T >
bool
equal_to(
        T           const & a,
        T           const & b,
        moris::real const & error_factor = 1.0e+06 );

template< class T, class A >
bool
equal_to(
        T           const & a,
        A           const & b,
        moris::real const & error_factor = 1.0e+06 );

    /**
     * Floating point equality comparison, specialized for complex numbers.
     *
     * @note Modeled after std::equal_to.
     *
     * @param[in] a            Floating point number.
     * @param[in] b            Floating point number.
     * @param[in] error_factor Number of machine epsilon differences
     *                         between a and b.
     *
     * @return Boolean true if equal.
     */
    template<>
    bool
    equal_to<moris::cplx>(
            moris::cplx const & a,
            moris::cplx const & b,
            moris::real const & error_factor );
}

// Template implementation file.
#include "fn_equal_to.tpp"

#endif /* MORIS_ALGORITHMS_FN_EQUAL_TO_HPP_ */

