/*
 * fn_bounding_box.hpp
 *
 *  Created on: Oct 22, 2019
 *      Author: sonne
 */

#ifndef PROJECTS_GEN_SRC_FN_BOUNDING_BOX_HPP_
#define PROJECTS_GEN_SRC_FN_BOUNDING_BOX_HPP_

#include "cl_Matrix.hpp"
#include "op_equal_equal.hpp"

namespace moris{
namespace ge{

inline
bool check_if_in_bounding_box( moris::Matrix<DDRMat> aEvalPoint,
                               moris::Matrix<DDRMat> aTestPointA,
                               moris::Matrix<DDRMat> aTestPointB,
                               moris::real aCriticalLength = 10.0 )
{
    MORIS_ASSERT(aEvalPoint.n_rows()<2,  "check_if_in_bounding_box() - evaluation point should be a row vector");
    MORIS_ASSERT(aTestPointA.n_rows()<2, "check_if_in_bounding_box() - test point should be a row vector");
    MORIS_ASSERT(aTestPointB.n_rows()<2, "check_if_in_bounding_box() - test point should be a row vector");
    bool tOutsideBox = false; // initially assume in box

    // check for 2D or 3D
    moris::uint tDim = aEvalPoint.n_cols();

    if ( tDim == 3 )
    {
        for ( uint k=0; k<3; ++k )
        {
            real tUpper = aEvalPoint(k) + aCriticalLength;
            if ( aTestPointA(k) > tUpper )
            {
                tOutsideBox = true;
                break;
            }
            if ( aTestPointB(k) > tUpper )
            {
                tOutsideBox = true;
                break;
            }

            real tLower = aEvalPoint(k) - aCriticalLength;
            if ( aTestPointA(k) < tLower )
            {
                tOutsideBox = true;
                break;
            }
            if ( aTestPointB(k) < tLower )
            {
                tOutsideBox = true;
                break;
            }
        }
    }

    if ( tDim == 2 )
    {
        for ( uint k=0; k<2; ++k )
        {
            real tUpper = aEvalPoint(k) + aCriticalLength;
            if ( aTestPointA(k) > tUpper )
            {
                tOutsideBox = true;
                break;
            }
            if ( aTestPointB(k) > tUpper )
            {
                tOutsideBox = true;
                break;
            }

            real tLower = aEvalPoint(k) - aCriticalLength;
            if ( aTestPointA(k) < tLower )
            {
                tOutsideBox = true;
                break;
            }
            if ( aTestPointB(k) < tLower )
            {
                tOutsideBox = true;
                break;
            }
        }
    }

    return tOutsideBox;
}


} // end ge namespace
} // end moris namespace


#endif /* PROJECTS_GEN_SRC_FN_BOUNDING_BOX_HPP_ */
