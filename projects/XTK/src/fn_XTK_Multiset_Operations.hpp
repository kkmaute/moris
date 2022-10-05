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
#include "cl_Cell.hpp"

using namespace moris;

namespace xtk
{
    //-------------------------------------------------------------------------------------

    // TODO: these operations can still be substantially optimized in how information is passed, sorted, and stored, but this is likely not a problem

    //-------------------------------------------------------------------------------------

    void
    multiset_union( 
            const Cell< moris_index >& aMultiSet1, 
            const Cell< moris_index >& aMultiSet2, 
                  Cell< moris_index >& aMultiSetUnion )
    {        
        // convert the first multiset into std::multisets
        std::multiset< moris_index > tFirstMultiSet;
        uint tFirstMultisetSize = aMultiSet1.size();
        for( uint i = 0; i < tFirstMultisetSize; i++ )
        {
            tFirstMultiSet.insert( aMultiSet1( i ) );
        }

        // convert the second multiset into std::multisets
        std::multiset< moris_index > tSecondMultiSet;
        uint tSecondMultiSetSize = aMultiSet2.size();
        for( uint i = 0; i < tSecondMultiSetSize; i++ )
        {
            tSecondMultiSet.insert( aMultiSet2( i ) );
        }

        // resize
        aMultiSetUnion.resize( tFirstMultisetSize + tSecondMultiSetSize );

        // get access to the data of the underlying vector
        std::vector< moris_index >& tVector = aMultiSetUnion.data();
        std::vector< moris_index >::iterator tIter;

        // perform union operation
        tIter = std::set_union( tFirstMultiSet.begin(), tFirstMultiSet.end(), tSecondMultiSet.begin(), tSecondMultiSet.end(), tVector.begin() );

        // get the used length
        uint tNumElems = ( uint )( tIter - tVector.begin() );

        // resize out unused space
        aMultiSetUnion.resize( tNumElems );
    }

    //-------------------------------------------------------------------------------------

    void
    multiset_difference( 
            const Cell< moris_index >& aMultiSet, 
            const Cell< moris_index >& aMultiSetToSubtract, 
                  Cell< moris_index >& aMultiSetDifference )
    {        
        // convert the first multiset into std::multisets
        std::multiset< moris_index > tMultiSet;
        uint tMultiSetSize = aMultiSet.size();
        for( uint i = 0; i < tMultiSetSize; i++ )
        {
            tMultiSet.insert( aMultiSet( i ) );
        }

        // convert the second multiset into std::multisets
        std::multiset< moris_index > tMultiSetToSubtract;
        uint tMultiSetToSubtractSize = aMultiSetToSubtract.size();
        for( uint i = 0; i < tMultiSetToSubtractSize; i++ )
        {
            tMultiSetToSubtract.insert( aMultiSetToSubtract( i ) );
        }

        // resize
        aMultiSetDifference.resize( tMultiSetSize );

        // get access to the data of the underlying vector
        std::vector< moris_index >& tVector = aMultiSetDifference.data();
        std::vector< moris_index >::iterator tIter;

        // perform union operation
        tIter = std::set_difference( tMultiSet.begin(), tMultiSet.end(), tMultiSetToSubtract.begin(), tMultiSetToSubtract.end(), tVector.begin() );

        // get the used length
        uint tNumElems = ( uint )( tIter - tVector.begin() );

        // resize out unused space
        aMultiSetDifference.resize( tNumElems );
    }

    //-------------------------------------------------------------------------------------

} // namespace xtk

#endif /* fn_XTK_Multiset_Operations.hpp */