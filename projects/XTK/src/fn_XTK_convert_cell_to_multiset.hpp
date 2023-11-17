/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * fn_XTK_convert_cell_to_multiset.hpp
 *
 */
#ifndef SRC_fn_XTK_convert_cell_to_multiset
#define SRC_fn_XTK_convert_cell_to_multiset

#include <set>
#include <algorithm>
#include "containers.hpp"

using namespace moris;

namespace xtk
{
    //-------------------------------------------------------------------------------------

    /**
     * @brief converts a cell of <T> into an ordered multiset of <T>
     *
     * @param aIndexCell moris::Cell<T>
     * @param aMultiSet std::multiset<T> to fill
     */
    template< class T >
    inline void
    convert_index_cell_to_index_multiset(
            Cell< T > const &   aIndexCell,
            std::multiset< T >& aMultiSet )
    {
        for ( uint i = 0; i < aIndexCell.size(); i++ )
        {
            aMultiSet.insert( aIndexCell( i ) );
        }
    }

    //-------------------------------------------------------------------------------------

}    // namespace xtk

#endif /* fn_XTK_convert_cell_to_multiset.hpp */
