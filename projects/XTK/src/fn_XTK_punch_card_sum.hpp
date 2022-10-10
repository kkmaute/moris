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
    
} // namespace xtk

#endif /* fn_XTK_punch_card_sum.hpp */