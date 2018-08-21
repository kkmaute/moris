/*
 * fn_num_digits.hpp
 *
 *  Created on: Jun 19, 2017
 *      Author: ktdoble
 */

#ifndef SRC_ALGORITHMS_FN_NUM_DIGITS_HPP_
#define SRC_ALGORITHMS_FN_NUM_DIGITS_HPP_

namespace xtk
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
     * @include "snippets/algorithms/fn_num_digits.inc"
     */
    template<typename T>
    size_t
    num_digits(
            T number)
    {
        size_t digits = 0;

        if (number < 0)
            digits = 1;

        while (number)
        {
            number /= 10;
            digits++;
        }

        return digits;
    }

}    // namespace xtk
#endif /* SRC_ALGORITHMS_FN_NUM_DIGITS_HPP_ */
