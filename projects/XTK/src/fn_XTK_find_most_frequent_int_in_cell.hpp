/*
 * Copyright (c) 2022 University of Colorado 
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details. 
 * 
 * ------------------------------------------------------------------------------------ 
 * 
 * fn_XTK_find_most_frequent_int_in_cell.hpp  
 * 
 */
#ifndef SRC_fn_XTK_find_most_frequent_int_in_cell
#define SRC_fn_XTK_find_most_frequent_int_in_cell

#include "cl_Cell.hpp"

using namespace moris;

namespace xtk
{
    /**
     * @brief finds the most occuring index in a cell of indices
     * 
     * @param aCellOfIndices list of indices  
     * @param aCount returns the number this index occurs
     * @return moris_index returns the index that occurs the most frequent in the list
     */
    moris_index find_most_frequent_index_in_cell( 
            moris::Cell< moris_index >& aCellOfIndices,
            uint&                       aCount )
    {
        uint tNumElemsInArray = aCellOfIndices.size();

        // list all array elements in hash.
        std::unordered_map< moris_index, moris_index > tHash;
        for( uint i = 0; i < tNumElemsInArray; i++ )
        {
            tHash[ aCellOfIndices( i ) ]++;
        }

        // find the array element with the max frequency
        moris_index tMaxCount = 0;
        moris_index tResult = -1;
        for( auto i : tHash ) 
        {
            if ( tMaxCount < i.second ) 
            {
                tResult = i.first;
                tMaxCount = i.second;
            }
        }

        // return the array element that occurs the most often
        aCount = (moris_index) tMaxCount;
        return tResult;        
    }
}

#endif /* fn_XTK_find_most_frequent_int_in_cell.hpp */