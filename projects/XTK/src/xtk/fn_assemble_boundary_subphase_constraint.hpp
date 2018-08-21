/*
 * fn_assemble_boundary_subphase_constraint.hpp
 *
 *  Created on: May 16, 2018
 *      Author: doble
 */

#ifndef SRC_XTK_FN_ASSEMBLE_BOUNDARY_SUBPHASE_CONSTRAINT_HPP_
#define SRC_XTK_FN_ASSEMBLE_BOUNDARY_SUBPHASE_CONSTRAINT_HPP_

// XTKL: Linalg Includes
#include "linalg/cl_XTK_Matrix_Base.hpp"

namespace xtk
{

/*
 * Given a element to element connectivity and element phase values, constructs the constraint matrix A
 * as seen in the system Ax = 0. x is a vector of elemental enrichment levels. For a given i in the system
 * x_i = x_j where j!= i. i.e
 *
 * A = [ 1  0 -1  0
 *       0  1  0 -1
 *      -1  0  1  0
 *       0 -1  0 -1];
 *
 *
 *  @param[in]  aElementToElement - Element to Element Connectivity Table (n_els x n_neighbors)
 *  @param[in]  aElementPhase     - Element Phase Value ( 1 x n_els )
 *  @param[in]  aDummyVal         - A value indicating an entry is not to be used (since not all elements have the same number of neighbors)
 *  @param[out] A - Constraint Matrix
 */
template<typename Real, typename Integer, typename Real_Matrix, typename Integer_Matrix>
static void
assemble_boundary_subphase_constraint(Matrix_Base<Integer, Integer_Matrix> const & aElementToElement,
                                      Matrix_Base<Integer, Integer_Matrix> const & aElementPhase,
                                      Integer const & aDummyVal,
                                      Matrix_Base<Real,Real_Matrix> & aA)
{

    // Number of elements
    Integer tNumElements  = aElementToElement.get_num_rows();
    Integer tNumNeighbors = aElementToElement.get_num_columns();

    // Initialize my phase and neighbor element phase
    Integer tMyPhase = aDummyVal;
    Integer tNeighborPhase = aDummyVal;

    // Initialize A such that there is the maximum number of constraint rows
    aA.resize(tNumElements*tNumNeighbors,tNumElements);
    aA.fill(0);


    // Count Number of rows used
    Integer tCount = 0;


    // Loop over all elements and their neighbors, check phase value and see if they are the same
    for( Integer i = 0;  i<tNumElements; i++)
    {

        // Get my phase value
        tMyPhase = aElementPhase(0,i);

        for( Integer j = 0; j<tNumNeighbors; j++)
        {
            // Move onto next element if this isn't a valid neighbor
            if(aElementToElement(i,j) == aDummyVal)
            {
                break;
            }

            // Get my neighbors phase
            tNeighborPhase = aElementPhase(0, aElementToElement(i,j) );

            // If my phase matches my neighbor's phase then constrain them to be equal
            if(tMyPhase == tNeighborPhase)
            {
                aA(tCount,i) = 1.0;
                aA(tCount,aElementToElement(i,j)) = -1.0;
                tCount++;
            }

        }
    }

    aA.resize(tCount,tNumElements);

}
}
#endif /* SRC_XTK_FN_ASSEMBLE_BOUNDARY_SUBPHASE_CONSTRAINT_HPP_ */
