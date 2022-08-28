/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_all_true.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_FN_ALL_TRUE_HPP_
#define PROJECTS_LINALG_SRC_FN_ALL_TRUE_HPP_

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

namespace moris
{

    /*
     * Checks whether all elements in a matrix evaluate to true.
     * For example if the result of an tMat1 == tMat2 holds
     * for all elements, this will return true.
     */
    inline
    bool
    all_true(Matrix< DDBMat > const & aMatrix)
    {
        bool tAllTrue = false;
        uint tNumRows = aMatrix.n_rows();
        uint tNumCols = aMatrix.n_cols();

        for(uint i = 0; i<tNumRows; i++)
        {
            for(uint j = 0; j<tNumCols; j++)
            {
                if(!aMatrix(i,j))
                {
                    return false;
                }
            }
        }

        tAllTrue = true;
        return tAllTrue;
    }

}

#endif /* PROJECTS_LINALG_SRC_FN_ALL_TRUE_HPP_ */

