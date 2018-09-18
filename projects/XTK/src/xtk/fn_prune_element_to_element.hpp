/*
 * fn_prune_element_to_element.hpp
 *
 *  Created on: May 16, 2018
 *      Author: doble
 */

#ifndef SRC_XTK_FN_PRUNE_ELEMENT_TO_ELEMENT_HPP_
#define SRC_XTK_FN_PRUNE_ELEMENT_TO_ELEMENT_HPP_

// Linear Algebra Includes
// XTKL: Linear Algebra Includes
#include "linalg/cl_XTK_Matrix_Base_Utilities.hpp"
#include "linalg/cl_XTK_Matrix.hpp"

#include "containers/cl_XTK_Cell.hpp"

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
template<typename Integer, typename Integer_Matrix>
Cell<moris::Matrix< Integer_Matrix >>
prune_element_to_element(moris::Matrix< Integer_Matrix > const & aElementToElement,
                         moris::Matrix< Integer_Matrix > const & aElementsInPrunedGraph,
                         moris::Matrix< Integer_Matrix > const & aSharedFaces,
                         Integer aDummyValue = std::numeric_limits<Integer>::max())
{

    // Number of elements in to included in pruned graph
    Integer tNumIncludedElems = aElementsInPrunedGraph.n_cols();

    XTK_ASSERT(aElementToElement.n_rows()== tNumIncludedElems,"Included elements and number of element neighbor relationships in element to element graph do not match");

    // Intialize pruned results where cell 0 is for element to element and cell 1 is the corresponding shared face
    moris::Matrix< Integer_Matrix > tPrunedElements(tNumIncludedElems, aElementToElement.n_cols(),aDummyValue);
    moris::Matrix< Integer_Matrix > tPrunedSharedFaces(tNumIncludedElems, aElementToElement.n_cols(),aDummyValue);

    // Generate Map
    std::unordered_map<Integer,Integer> tElementMap(tNumIncludedElems);
    for(Integer iE = 0; iE<tNumIncludedElems; iE++)
    {
        tElementMap[aElementsInPrunedGraph(0,iE)] = iE;
    }

    for(Integer iE = 0; iE< tNumIncludedElems; iE++)
    {
        Integer tCount = 0;
        for(Integer iNE = 0; iNE < aElementToElement.n_cols(); iNE++)
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
