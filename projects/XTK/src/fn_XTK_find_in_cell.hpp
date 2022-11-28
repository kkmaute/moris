/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * fn_XTK_find_in_cell.hpp
 *
 */
#ifndef SRC_fn_XTK_find_in_cell
#define SRC_fn_XTK_find_in_cell

#include "cl_Cell.hpp"

using namespace moris;
namespace xtk
{
    //-------------------------------------------------------------------------------------

    /**
     * @brief finds the first instance of a value in a moris cell
     *
     * @tparam T value type
     * @param aList cell to find the value in
     * @param aKeyToFind value to find
     * @return moris_index index at which the value is first found
     */
    template< class T >
    inline moris_index
    find_in_cell(
            Cell< T > const & aList,
            T                 aKeyToFind )
    {
        // get the size of the cell
        moris_index tSize = (moris_index)aList.size();

        // go over elements and find
        for ( moris_index i = 0; i < tSize; i++ )
        {
            // if found, stop search and return index
            if ( aList( i ) == aKeyToFind )
            {
                return i;
            }
        }

        // if not found, return -1
        return -1;
    }

    //-------------------------------------------------------------------------------------

    /**
     * @brief finds the first instance of a value in a moris cell that has not been found before
     *
     * @tparam T value type
     * @param aList cell to find the value in
     * @param aKeyToFind value to find
     * @param aPunchCard boolean array that marks objects that have already been found with true
     * @return moris_index index at which the value is first found and has not been found before
     */
    template< class T >
    inline moris_index
    find_in_cell(
            Cell< T > const & aList,
            T                 aKeyToFind,
            Cell< bool >&     aPunchCard )
    {
        // get the size of the cell
        moris_index tSize = (moris_index)aList.size();

        // check that the Punchcard makes sense
        MORIS_ASSERT( aPunchCard.size() == (uint)tSize,
                "xtk::find_in_cell() - The cell and punch card provided have different sizes" );

        // go over elements and find
        for ( moris_index i = 0; i < tSize; i++ )
        {
            // check if this is the value looked for
            bool tValueFound = ( aList( i ) == aKeyToFind );

            // check if found and not been found before
            if ( tValueFound && ( !aPunchCard( i ) ) )
            {
                // update the punch card
                aPunchCard( i ) = true;

                // return the current index
                return i;
            }
        }

        // if not found, return -1
        return -1;
    }

    //-------------------------------------------------------------------------------------

}    // namespace xtk

#endif /* fn_XTK_find_in_cell.hpp */
