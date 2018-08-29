/*
 * fn_flood_fill.hpp
 *
 *  Created on: Feb 9, 2018
 *      Author: ktdoble
 */

#ifndef XTK_SRC_XTK_FN_MESH_FLOOD_FILL_HPP_
#define XTK_SRC_XTK_FN_MESH_FLOOD_FILL_HPP_

// XTKL: Linear Algebra Includes
#include "linalg/cl_XTK_Matrix_Base_Utilities.hpp"
#include "linalg/cl_XTK_Matrix_Base.hpp"


// Unordered Map Include
#include <unordered_map>


#include "cl_XTK_Child_Mesh.hpp"


namespace xtk
{


/* full_flood_fill runs a floodfill algorithm traversing the element to element connectivity returning subphases
 * values at each index on a subdomain indicated by aActiveElements
 *
 * @param[in]  aElementToElement   - Element to Element connectivity
 * @param[in]  aElementPhaseIndex  - Element Phase Index (row vector)
 * @param[in]  aActiveElements     - Short list of active elements to limit elements needed to iterate over
 * @param[in]  aElementsToInclude  - Marks the elements to consider == 1 and the ones to ignore == 0 (could be expensive to create inside if large mesh)
 * @param[in]  aNumPhases          - Number of Phases possible
 * @param[in]  aIncludeAllElements - Says to include all elements
 * @param[out] Element Subphase Indices
 *
 */
template<typename Integer, typename Integer_Matrix>
moris::Matrix<Integer, Integer_Matrix>
flood_fill( moris::Matrix<Integer, Integer_Matrix> const & aElementToElement,
            moris::Matrix<Integer, Integer_Matrix> const & aElementPhaseIndex,
            moris::Matrix<Integer, Integer_Matrix> const & aActiveElements,
            moris::Matrix<Integer, Integer_Matrix> const & aElementsToInclude,
            Integer                              aNumPhases,
            Integer                              aDummyValue,
            bool aIncludeAllElements = false)
            {

    // Active phase index
    Integer tPhaseIndex  = 0;

    // Number of elements in the flood fill
    Integer tNumElements = aActiveElements.n_cols();

    // Number of Elements with Set Phases (This allows for early termination of code if every element has been set)
    Integer tNumPhasesSet = 0;

    // Maximum number of neighbors per element
    Integer tMaxNumNeighbors = aElementToElement.n_cols();

    // Current Element Index
    Integer tElementIndex = 0;

    // Neighbor Element Phase
    Integer tNeighborPhase = 0;

    // Neighbor Element local Index in the list aElementsToInclude
    Integer tNeighborIndex = 0;

    // Current Subphase
    Integer tCurrentSubphase = 0;

    // Track which elements have their phase set
    moris::Matrix<Integer, Integer_Matrix> tPhaseSet(1,tNumElements,0);

    // Initialize element sub-phases
    moris::Matrix<Integer, Integer_Matrix> tElementSubphase(1,tNumElements,aDummyValue);

    // Initialize Active Front
    Integer tActiveFrontCount = 0;
    Integer tActiveFrontElement = 0;
    moris::Matrix<Integer, Integer_Matrix> tActiveFront(1,tNumElements,0);

    // Map between the active element indexes provided and their corresponding iE (Only needed if all elements are not included)
    // key   - Element Index
    // value - flood fill local index
    std::unordered_map<Integer,Integer> tElementToLocalIndex;
    if(!aIncludeAllElements)
    {
        for(Integer iE = 0; iE<tNumElements; iE++)
        {
            tElementToLocalIndex[aActiveElements(0,iE)] = iE;
        }
    }

    // Loop over all elements
    for(Integer iE = 0; iE<tNumElements; iE ++)
    {
        tElementIndex = aActiveElements(0,iE);

        // If this element phase has not been set
        if(!tPhaseSet(0,iE))
        {
            // Phase Index of the element
            tPhaseIndex = aElementPhaseIndex(0,aActiveElements(0,iE));


            // Set the elements subphase value
            tElementSubphase(0,iE) = tCurrentSubphase;

            // Mark this element as set
            tPhaseSet(0,iE) = 1;

            // Update active front
            for(Integer iN = 0; iN<tMaxNumNeighbors; iN++)
            {
                // Move on if we see a dummy value
                if(aElementToElement(tElementIndex,iN) ==  aDummyValue)
                {
                    break;
                }

                tNeighborPhase = aElementPhaseIndex(0,aElementToElement(tElementIndex,iN));



                if(!aIncludeAllElements)
                {
                    tNeighborIndex = tElementToLocalIndex[aElementToElement(tElementIndex,iN)];
                }

                // Otherwise the neighbor element index is located easily without a map
                else
                {
                    tNeighborIndex = aElementToElement(tElementIndex,iN);
                }


                // If this is an neighbor element to include in the subdomain, has not already
                // been set and its phase matches the current elements phase then
                // add it to the active front and increment the count
                if(aElementsToInclude(0,aElementToElement(tElementIndex,iN)) == 1 &&
                   tPhaseSet(0,tNeighborIndex)!=1 &&
                   tNeighborPhase == tPhaseIndex)
                {
                    tActiveFront(0,tActiveFrontCount) = aElementToElement(tElementIndex,iN);
                    tActiveFrontCount++;
                }
            }

            // Iterate through active front until there are no more elements in the active front
            // We start at the end of the front and work backwards
            while(tActiveFrontCount!=0)
            {

                // Current Element Index in the Active Front
                tActiveFrontElement = tActiveFront(0,tActiveFrontCount-1);
                // Get Neighbor index from map if we're not considering the full domain
                if(!aIncludeAllElements)
                {
                    tNeighborIndex = tElementToLocalIndex[tActiveFront(0,tActiveFrontCount-1)];
                }

                // Otherwise the neighbor element index is located easily without a map
                else
                {
                    tNeighborIndex = tActiveFrontElement;
                }

                // Get the neighbors phase
                tNeighborPhase = aElementPhaseIndex(0,tActiveFrontElement);

                // If the neighbor phase matches our phase, then we add it's neighbor to the active front
                // Unless it has already been set
                if(tNeighborPhase == tPhaseIndex &&
                  tPhaseSet(0,tNeighborIndex)!=1)
                {
                    // Set the neighbor elements subphase value
                    tElementSubphase(0,tNeighborIndex) = tCurrentSubphase;

                    // Mark element as set
                    tPhaseSet(0,tNeighborIndex) = 1;

                    // Increase the number of phases set
                    tNumPhasesSet++;


                    // Add the elements other neighbors to the active front
                    bool tReplaced = false;
                    for(Integer i = 0; i<tMaxNumNeighbors; i++)
                    {

                        tElementIndex = aElementToElement(tActiveFrontElement,i);

                        // Get Neighbor index from map if we're not considering the full domain
                        if(!aIncludeAllElements)
                        {
                            tNeighborIndex = tElementToLocalIndex[tElementIndex];
                        }

                        // Otherwise the neighbor element index is located easily without a map
                        else
                        {
                            tNeighborIndex = tElementIndex;
                        }


                        // If the elements neighbor is a dummy value
                        if(tElementIndex == aDummyValue)
                        {
                            if(!tReplaced)
                            {
                                tActiveFrontCount--;
                            }
                            break;
                        }

                        // If this element is active and its phase hasn't been set
                        if(   aElementsToInclude(0,tElementIndex) == 1 &&
                            tPhaseSet(0,tNeighborIndex) != 1 )
                        {
                            // and the previous element hasn't been replaced, then replace it
                            // don't add to active front count
                            if(!tReplaced)
                            {
                                tActiveFront(0,tActiveFrontCount-1) = tElementIndex;
                                tReplaced = true;
                            }

                            // Else add to end of active front and add to the count
                            else
                            {
                                tActiveFront(0,tActiveFrontCount) = tElementIndex;
                                tActiveFrontCount++;
                            }
                        }

                        // If we reached the end and haven't replaced the element from before, then remove it.
                        if( i == tMaxNumNeighbors-1 && !tReplaced)
                        {
                            tActiveFront(0,tActiveFrontCount-1) = aDummyValue;
                            tActiveFrontCount--;
                        }
                    }

                }

                // Else if the phase doesn't match we remove that element from the active front
                else
                {
                    tActiveFront(0,tActiveFrontCount-1) = aDummyValue;
                    tActiveFrontCount--;
                }
            }

            tCurrentSubphase++;
        }
    }

    return tElementSubphase;
            }
// Full Mesh Version

}


#endif /* XTK_SRC_XTK_FN_MESH_FLOOD_FILL_HPP_ */
