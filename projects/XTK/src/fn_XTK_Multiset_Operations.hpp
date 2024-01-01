/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * fn_XTK_Multiset_Operations.hpp
 *
 */
#ifndef SRC_fn_XTK_Multiset_Operations
#define SRC_fn_XTK_Multiset_Operations

#include <set>
#include <algorithm>
#include "cl_Vector.hpp"
#include "fn_XTK_convert_cell_to_multiset.hpp"

using namespace moris;

namespace xtk
{
    //-------------------------------------------------------------------------------------

    // TODO: these operations can still be substantially optimized in how information is passed, sorted, and stored, but this is likely not a problem

    //-------------------------------------------------------------------------------------

    inline void
    multiset_union(
            const Vector< moris_index >& aMultiSet1,
            const Vector< moris_index >& aMultiSet2,
            Vector< moris_index >&       aMultiSetUnion )
    {
        // convert the first multiset into std::multisets
        std::multiset< moris_index > tFirstMultiSet;
        xtk::convert_index_cell_to_index_multiset( aMultiSet1, tFirstMultiSet );

        // convert the second multiset into std::multisets
        std::multiset< moris_index > tSecondMultiSet;
        xtk::convert_index_cell_to_index_multiset( aMultiSet2, tSecondMultiSet );

        // resize
        aMultiSetUnion.resize( aMultiSet1.size() + aMultiSet2.size() );

        // get access to the data of the underlying vector
        std::vector< moris_index >&          tVector = aMultiSetUnion.data();
        std::vector< moris_index >::iterator tIter;

        // perform union operation
        tIter = std::set_union( tFirstMultiSet.begin(), tFirstMultiSet.end(), tSecondMultiSet.begin(), tSecondMultiSet.end(), tVector.begin() );

        // get the used length
        uint tNumElems = (uint)( tIter - tVector.begin() );

        // resize out unused space
        aMultiSetUnion.resize( tNumElems );
    }

    //-------------------------------------------------------------------------------------

    void
    multiset_difference(
            const Vector< moris_index >& aMultiSet,
            const Vector< moris_index >& aMultiSetToSubtract,
            Vector< moris_index >&       aMultiSetDifference )
    {
        // convert the first multiset into std::multisets
        std::multiset< moris_index > tMultiSet;
        xtk::convert_index_cell_to_index_multiset( aMultiSet, tMultiSet );

        // convert the second multiset into std::multisets
        std::multiset< moris_index > tMultiSetToSubtract;
        xtk::convert_index_cell_to_index_multiset( aMultiSetToSubtract, tMultiSetToSubtract );

        // resize
        aMultiSetDifference.resize( aMultiSet.size() );

        // get access to the data of the underlying vector
        std::vector< moris_index >&          tVector = aMultiSetDifference.data();
        std::vector< moris_index >::iterator tIter;

        // perform difference operation on the two sets
        tIter = std::set_difference( tMultiSet.begin(), tMultiSet.end(), tMultiSetToSubtract.begin(), tMultiSetToSubtract.end(), tVector.begin() );

        // get the used length
        uint tNumElems = (uint)( tIter - tVector.begin() );

        // resize out unused space
        aMultiSetDifference.resize( tNumElems );
    }

    //-------------------------------------------------------------------------------------

}    // namespace xtk

#endif /* fn_XTK_Multiset_Operations.hpp */
