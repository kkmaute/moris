/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_num_digits.hpp
 *
 */

#ifndef MORIS_ALGORITHMS_FN_NUM_DIGITS_HPP_
#define MORIS_ALGORITHMS_FN_NUM_DIGITS_HPP_

namespace moris
{
    /**
     * @brief Determine number of digits in an integer.
     *
     * @see [Efficient way to determine number of digits in an integer]
     * (http://stackoverflow.com/questions/1489830/efficient-way-to-determine-number-of-digits-in-an-integer)
     *
     * @param[in] number Number parameter.
     * @return Number of digits.
     *
     * @include "snippets/algorithms/fn_num_digits.inc" // snippets ALG
     */
    template<typename T>
    moris::size_t
    num_digits(
            T number)
    {
        moris::size_t digits = 0;

        if (number < 0)
            digits = 1;

        while (number)
        {
            number /= 10;
            digits++;
        }

        return digits;
    }

}    // namespace moris

#endif    /* MORIS_ALGORITHMS_FN_NUM_DIGITS_HPP_ */

