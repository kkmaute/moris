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

#include "cl_Cell.hpp"

using namespace moris;
namespace xtk
{
    //-------------------------------------------------------------------------------------

    // define an Index-map
    typedef std::unordered_map< moris::moris_index, moris::moris_index > IndexMap;
    typedef std::unordered_map< moris::moris_id, moris::moris_id >       IdMap;
    typedef std::unordered_map< moris::moris_id, moris::moris_index >    IdToIndexMap;

    //-------------------------------------------------------------------------------------

    template< class T >
    inline void
    convert_cell_to_map(
            moris::Cell< T > const &              aCellToConvert,
            std::unordered_map< T, moris_index >& aIndexMapToFill )
    {
        for ( uint i = 0; i < aCellToConvert.size(); i++ )
        {
            aIndexMapToFill[ aCellToConvert( i ) ] = i;
        }
    }

    //-------------------------------------------------------------------------------------

}    // namespace xtk

#endif /* fn_XTK_convert_cell_to_map.hpp */
