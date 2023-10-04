/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_XTK_convert_cell_to_map.hpp
 *
 */

#ifndef SRC_fn_XTK_convert_cell_to_map
#define SRC_fn_XTK_convert_cell_to_map

#include "containers.hpp"

using namespace moris;
namespace xtk
{
    //-------------------------------------------------------------------------------------

    /**
     * @brief convert a cell of typename T to a moris::Mini_Map
     * where the keys are the contents in the cell, and they map to their index in the cell
     */
    template< class T >
    inline void
    convert_cell_to_map(
            moris::Cell< T > const &    aCellToConvert,
            Mini_Map< T, moris_index >& aIndexMapToFill )
    {
        for ( uint i = 0; i < aCellToConvert.size(); i++ )
        {
            aIndexMapToFill[ aCellToConvert( i ) ] = i;
        }
    }

    //-------------------------------------------------------------------------------------

//     /**
//      * @brief convert a cell of typename T to an std::unordered_map
//      * where the keys are the contents in the cell, and they map to their index in the cell
//      */
//     template< class T >
//     inline void
//     convert_cell_to_map(
//             moris::Cell< T > const &              aCellToConvert,
//             std::unordered_map< T, moris_index >& aIndexMapToFill )
//     {
//         for ( uint i = 0; i < aCellToConvert.size(); i++ )
//         {
//             aIndexMapToFill[ aCellToConvert( i ) ] = i;
//         }
//     }

    //-------------------------------------------------------------------------------------

}    // namespace xtk

#endif /* fn_XTK_convert_cell_to_map.hpp */
