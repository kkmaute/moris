/*
 * cl_Test_Functions.hpp
 *
 *  Created on: Sep 1, 2017
 *      Author: doble
 */

#ifndef TEST_INCLUDE_CL_TEST_ALGORITHMS_HPP_
#define TEST_INCLUDE_CL_TEST_ALGORITHMS_HPP_

#include "catch.hpp"

#include "cl_Mat.hpp" // LNA/src


/*
 * For functions which are convenient to use for the test
 */
namespace moris
{
    template< typename T>
    bool equal_to(Mat<T> const & aMat1,
                  Mat<T> const & aMat2)
    {

        if( aMat1.n_cols() != aMat2.n_cols())
        {
            MORIS_LOG_ERROR<<"Unequal number of rows columns";
            return false;
        }

        if( aMat1.n_rows() != aMat2.n_rows())
        {
            MORIS_LOG_ERROR<<"Unequal  number of rows";
            return false;
        }

        uint tNumRows = aMat1.n_rows();
        uint tNumCols = aMat1.n_cols();

        for( uint iCol = 0; iCol <tNumCols; iCol++)
        {
            for (uint iRow = 0; iRow <tNumRows; iRow++)
            {
                if(aMat1(iRow,iCol) == Approx( aMat2(iRow,iCol)))
                {

                }
                else
                {
                    return false;
                }
            }
        }

        return true;
    }
}


#endif /* TEST_INCLUDE_CL_TEST_ALGORITHMS_HPP_ */
