/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_mesh_flood_fill.hpp
 *
 */

#pragma once

// XTKL: Linear Algebra Includes
#include "cl_Matrix.hpp"

// Unordered Map Include
#include <unordered_map>

namespace moris::mtk
{
    /* full_flood_fill runs a floodfill algorithm traversing the element to element connectivity returning subphases
     * values at each index on a subdomain indicated by aActiveElements
     *
     * @param[in]  aElementToElement   - Element to Element connectivity size < numElements > x < maxNeighbors >
     * @param[in]  aElementPhaseIndex  - Element Phase Index (row vector)  size < numElements >
     * @param[in]  aActiveElements     - Short list of active elements to limit elements needed to iterate over
     * @param[in]  aElementsToInclude  - Marks the elements to consider == 1 and the ones to ignore == 0 (could be expensive to create inside if large mesh)
     * @param[in]  aIncludeAllElements - Says to include all elements
     * @param[out] Element Subphase Indices
     *
     */
    inline Matrix< IndexMat >
    flood_fill(
            Matrix< IndexMat > const &    aElementToElement,
            Vector< moris_index > const & aElementPhaseIndex,
            Vector< moris_index > const & aActiveElements,
            Vector< moris_index > const & aElementsToInclude,
            moris::moris_index            aDummyValue,
            moris::moris_index&           aMaxValueAssigned,
            bool                          aIncludeAllElements = false )
    {
        // Active phase index
        moris::size_t tPhaseIndex = 0;

        // Number of elements in the flood fill
        moris::size_t tNumElements = aActiveElements.size();

        // Number of Elements with Set Phases (This allows for early termination of code if every element has been set)
        moris::size_t tNumPhasesSet = 0;

        // Maximum number of neighbors per element
        moris::size_t tMaxNumNeighbors = aElementToElement.n_cols();

        // Current Element Index
        moris::moris_index tElementIndex = 0;

        // Neighbor Element Phase
        moris::size_t tNeighborPhase = 0;

        // Neighbor Element local Index in the list aElementsToInclude
        moris::moris_index tNeighborIndex = 0;

        // Current Subphase
        moris::size_t tCurrentSubphase = 0;

        // Track which elements have their phase set
        Matrix< moris::DDBMat > tPhaseSet( 1, tNumElements, 0 );

        // Initialize element sub-phases
        Matrix< IndexMat > tElementSubphase( tNumElements, 1, aDummyValue );

        // Initialize Active Front
        moris::size_t      tActiveFrontCount   = 0;
        moris::size_t      tActiveFrontElement = 0;
        Matrix< IndexMat > tActiveFront( 1, 10 * tNumElements + 1, 0 );

        // Map between the active element indexes provided and their corresponding iE (Only needed if all elements are not included)
        // key   - Element Index
        // value - flood fill local index
        std::unordered_map< moris::moris_index, moris::moris_index > tElementToLocalIndex;
        if ( !aIncludeAllElements )
        {
            for ( moris::size_t iE = 0; iE < tNumElements; iE++ )
            {
                tElementToLocalIndex[ aActiveElements( iE ) ] = iE;
            }
        }

        // Loop over all elements
        for ( moris::size_t iE = 0; iE < tNumElements; iE++ )
        {
            tElementIndex = aActiveElements( iE );

            // If this element phase has not been set
            if ( !tPhaseSet( 0, iE ) )
            {
                // Phase Index of the element
                tPhaseIndex = aElementPhaseIndex( aActiveElements( iE ) );

                // Set the elements subphase value
                tElementSubphase( iE, 0 ) = tCurrentSubphase;

                // Mark this element as set
                tPhaseSet( 0, iE ) = 1;

                // Update active front
                for ( moris::size_t iN = 0; iN < tMaxNumNeighbors; iN++ )
                {
                    // Move on if we see a dummy value
                    if ( aElementToElement( tElementIndex, iN ) == aDummyValue )
                    {
                        break;
                    }

                    tNeighborPhase = aElementPhaseIndex( aElementToElement( tElementIndex, iN ) );

                    if ( !aIncludeAllElements )
                    {
                        tNeighborIndex = tElementToLocalIndex[ aElementToElement( tElementIndex, iN ) ];
                    }

                    // Otherwise the neighbor element index is located easily without a map
                    else
                    {
                        tNeighborIndex = aElementToElement( tElementIndex, iN );
                    }

                    // If this is a neighbor element to include in the subdomain, has not already
                    // been set and its phase matches the current element's phase then
                    // add it to the active front and increment the count
                    if ( aElementsToInclude( aElementToElement( tElementIndex, iN ) ) == 1 &&    //
                            tPhaseSet( 0, tNeighborIndex ) != 1 &&                               //
                            tNeighborPhase == tPhaseIndex )
                    {
                        tActiveFront( 0, tActiveFrontCount ) = aElementToElement( tElementIndex, iN );
                        tActiveFrontCount++;
                    }

                }    // end for: each neighbor of the current element

                // Iterate through active front until there are no more elements in the active front
                // We start at the end of the front and work backwards
                while ( tActiveFrontCount != 0 )
                {
                    // Current Element Index in the Active Front
                    tActiveFrontElement = tActiveFront( 0, tActiveFrontCount - 1 );
                    // Get Neighbor index from map if we're not considering the full domain
                    if ( !aIncludeAllElements )
                    {
                        tNeighborIndex = tElementToLocalIndex[ tActiveFront( 0, tActiveFrontCount - 1 ) ];
                    }

                    // Otherwise the neighbor element index is located easily without a map
                    else
                    {
                        tNeighborIndex = tActiveFrontElement;
                    }

                    // Get the neighbors phase
                    tNeighborPhase = aElementPhaseIndex( tActiveFrontElement );

                    // If the neighbor phase matches our phase, then we add it's neighbor to the active front
                    // Unless it has already been set
                    if ( tNeighborPhase == tPhaseIndex && tPhaseSet( 0, tNeighborIndex ) != 1 )
                    {
                        // Set the neighbor elements subphase value
                        tElementSubphase( tNeighborIndex, 0 ) = tCurrentSubphase;

                        // Mark element as set
                        tPhaseSet( 0, tNeighborIndex ) = 1;

                        // Increase the number of phases set
                        tNumPhasesSet++;

                        // Add the elements other neighbors to the active front
                        bool tReplaced = false;
                        for ( moris::size_t i = 0; i < tMaxNumNeighbors; i++ )
                        {
                            tElementIndex = aElementToElement( tActiveFrontElement, i );

                            // Get Neighbor index from map if we're not considering the full domain
                            if ( !aIncludeAllElements )
                            {
                                tNeighborIndex = tElementToLocalIndex[ tElementIndex ];
                            }

                            // Otherwise the neighbor element index is located easily without a map
                            else
                            {
                                tNeighborIndex = tElementIndex;
                            }

                            // If the elements neighbor is a dummy value
                            if ( tElementIndex == aDummyValue )
                            {
                                if ( !tReplaced )
                                {
                                    tActiveFrontCount--;
                                }
                                break;
                            }

                            // If this element is active and its phase hasn't been set
                            if ( aElementsToInclude( tElementIndex ) == 1 &&    //
                                    tPhaseSet( 0, tNeighborIndex ) != 1 )
                            {
                                // and the previous element hasn't been replaced, then replace it
                                // don't add to active front count
                                if ( !tReplaced )
                                {
                                    tActiveFront( 0, tActiveFrontCount - 1 ) = tElementIndex;

                                    tReplaced = true;
                                }

                                // Else add to end of active front and add to the count
                                else
                                {

                                    MORIS_ASSERT( tActiveFrontCount < tActiveFront.numel(),
                                            " Active front in flood fill not big enough" );

                                    tActiveFront( 0, tActiveFrontCount ) = tElementIndex;
                                    tActiveFrontCount++;
                                }
                            }

                            // If we reached the end and haven't replaced the element from before, then remove it.
                            if ( i == tMaxNumNeighbors - 1 && !tReplaced )
                            {
                                tActiveFront( 0, tActiveFrontCount - 1 ) = aDummyValue;
                                tActiveFrontCount--;
                            }
                        }
                    }

                    // Else if: the phase doesn't match we remove that element from the active front
                    else
                    {
                        tActiveFront( 0, tActiveFrontCount - 1 ) = aDummyValue;
                        tActiveFrontCount--;
                    }

                }    // end while: active front is not empty

                tCurrentSubphase++;

            }    // end if: enrichment level has not been determined yet
        }    // end for: all elements in graph

        aMaxValueAssigned = tCurrentSubphase - 1;

        return tElementSubphase;

    }    // end function: flood_fill()
}    // namespace moris::mtk