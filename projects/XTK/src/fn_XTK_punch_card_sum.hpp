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

#include "cl_Cell.hpp"

using namespace moris;

namespace xtk
{

    //-------------------------------------------------------------------------------------

    /**
     * @brief counts the number of positive entries in a punch card
     * 
     * @param aPunchCard cell of booleans
     * @return uint number of true entries in the cell
     */
    uint
    punch_card_sum( const Cell< bool >& aPunchCard )
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
     * @brief generalized sum over all entries of a moris::Cell
     * 
     * @tparam T type of the Cell entries
     * @param aCell moris::Cell to sum up
     * @return T sum over all entries
     */
    template< class T >
    T
    sum_over_cell( Cell< T > const & aCell )
    {
        // initialize sum with 0
        T tSum = (T)0;

        // sum up entries of cell
        for( uint i = 0; i < aCell.size(); i++ )
        {
            tSum += aCell( i );
        }

        // return sum
        return tSum;
    }

    //-------------------------------------------------------------------------------------
    
} // namespace xtk

#endif /* fn_XTK_punch_card_sum.hpp */