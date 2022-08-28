/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_XTK_convert_index_cell_to_index_map.hpp
 *
 */

#ifndef SRC_fn_XTK_convert_index_cell_to_index_map
#define SRC_fn_XTK_convert_index_cell_to_index_map

#include "cl_Cell.hpp"

using namespace moris;
namespace xtk
{
    // define an Index-map
    typedef std::unordered_map< moris::moris_index, moris::moris_index > IndexMap;

    //-------------------------------------------------------------------------------------

    void
    convert_index_cell_to_index_map(
        moris::Cell< moris_index > const& aCellToConvert,
        IndexMap&                         aIndexMapToFill )
    {
        for( uint i = 0; i < aCellToConvert.size(); i++ )
        {
            aIndexMapToFill[ aCellToConvert( i ) ] = i;
        }
    }

    //-------------------------------------------------------------------------------------
}

#endif /* fn_XTK_convert_index_cell_to_index_map.hpp */
