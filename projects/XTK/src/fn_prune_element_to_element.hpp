/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_prune_element_to_element.hpp
 *
 */

#ifndef SRC_XTK_FN_PRUNE_ELEMENT_TO_ELEMENT_HPP_
#define SRC_XTK_FN_PRUNE_ELEMENT_TO_ELEMENT_HPP_

// Linear Algebra Includes
// XTKL: Linear Algebra Includes
#include "cl_XTK_Matrix_Base_Utilities.hpp"
#include "cl_Matrix.hpp"

#include "cl_Cell.hpp"

// Unordered Map Include
#include <unordered_map>

namespace xtk
{

/*
 * Starting with a full element to element graph, the graph is pruned removing elements that are not desired in the pruned graph.
 * Also prunes the shared face lists
 *
 * @param[in] aElementToElement - Full element to element graph
 * @param[in] aElementsInPrunedGraph -  Elements to include in the pruned graph
 * @param[in] aSharedFaces - Shared face list with element neighbor
 * @param[in] aDummyValue - Dummy value needed since not all elements have the same number of neighbors
 * @param[out] Cell(0) Pruned element to element graph only including the elements found in the aElementsInPrunedGraph vector
 *             Cell(1) Pruned shared faces
 *
 */
inline
Cell<moris::Matrix< moris::IndexMat >>
prune_element_to_element(moris::Matrix< moris::IndexMat > const & aElementToElement,
                         moris::Matrix< moris::IndexMat > const & aElementsInPrunedGraph,
                         moris::Matrix< moris::IndexMat > const & aSharedFaces,
                         moris::moris_index                       aDummyValue = MORIS_INDEX_MAX)
{

    // Number of elements in to included in pruned graph
    size_t tNumIncludedElems = aElementsInPrunedGraph.numel();

    MORIS_ASSERT(aElementToElement.n_rows()== tNumIncludedElems,"Included elements and number of element neighbor relationships in element to element graph do not match");

    // Intialize pruned results where cell 0 is for element to element and cell 1 is the corresponding shared face
    moris::Matrix< moris::IndexMat > tPrunedElements(tNumIncludedElems, aElementToElement.n_cols(),aDummyValue);
    moris::Matrix< moris::IndexMat > tPrunedSharedFaces(tNumIncludedElems, aElementToElement.n_cols(),aDummyValue);

    // Generate Map
    std::unordered_map<moris::moris_index,moris::moris_index> tElementMap(tNumIncludedElems);
    for(size_t iE = 0; iE<tNumIncludedElems; iE++)
    {
        tElementMap[aElementsInPrunedGraph(iE)] = iE;
    }

    for(size_t iE = 0; iE< tNumIncludedElems; iE++)
    {
        size_t tCount = 0;
        for(size_t iNE = 0; iNE < aElementToElement.n_cols(); iNE++)
        {
            // If the element is in the map
            if(tElementMap.find(aElementToElement(iE,iNE)) != tElementMap.end())
            {
                (tPrunedElements)(iE,tCount) = aElementToElement(iE,iNE);
                (tPrunedSharedFaces)(iE,tCount) = aSharedFaces(iE,iNE);
                tCount++;
            }
        }
    }
    return {tPrunedElements,tPrunedSharedFaces};
}

}

#endif /* SRC_XTK_FN_PRUNE_ELEMENT_TO_ELEMENT_HPP_ */

