/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * fn_XTK_punch_card_sum.hpp
 *
 */
#ifndef SRC_fn_XTK_punch_card_sum
#define SRC_fn_XTK_punch_card_sum

#include "containers.hpp"

using namespace moris;

namespace moris::xtk
{

    //-------------------------------------------------------------------------------------

    /**
     * @brief counts the number of positive entries in a punch card
     *
     * @param aPunchCard cell of booleans
     * @return uint number of true entries in the cell
     */
    inline uint
    punch_card_sum( const Vector< bool >& aPunchCard )
    {
        uint tSum = 0;

        for ( uint i = 0; i < aPunchCard.size(); i++ )
        {
            tSum += aPunchCard( i );
        }

        return tSum;
    }

    //-------------------------------------------------------------------------------------

    /**
     * @brief generalized sum over all entries of a Vector
     *
     * @tparam T type of the Cell entries
     * @param aCell Vector to sum up
     * @return T sum over all entries
     */
    template< class T >
    inline T
    sum_over_cell( Vector< T > const & aCell )
    {
        // initialize sum with 0
        T tSum = (T)0;

        // sum up entries of cell
        for ( uint i = 0; i < aCell.size(); i++ )
        {
            tSum += aCell( i );
        }

        // return sum
        return tSum;
    }

    //-------------------------------------------------------------------------------------

}    // namespace moris::xtk

#endif /* fn_XTK_punch_card_sum.hpp */
