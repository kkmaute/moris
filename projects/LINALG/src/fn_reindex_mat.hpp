/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_reindex_mat.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_FN_REINDEX_MAT_HPP_
#define PROJECTS_LINALG_SRC_FN_REINDEX_MAT_HPP_

#include "cl_Matrix.hpp"

/*
 *
 * I.e
 * aMatToReindex = [2,3,4,5;
 *                  1,2,3,4;
 *                  4,2,1,1];
 * aMatRow = 0;
 * aIndexMap = [0,1,2,3;
 *              3,2,1,0];
 *
 * would result in
 *
 * [2,3,4,5;
 *  5,4,3,2];
 *
 */

template<typename Matrix_Type, typename Integer_Matrix>
moris::Matrix< Matrix_Type >
reindex_mat(moris::Matrix<Integer_Matrix> const & aIndexMap,
            size_t                        const & aMatRow,
            moris::Matrix< Matrix_Type >  const & aMatToReindex)
{
    size_t tNumRow = aIndexMap.n_rows();
    size_t tNumCol = aIndexMap.n_cols();
    moris::Matrix< Matrix_Type > tReindexedMat(tNumRow,tNumCol);

    for( size_t i = 0; i<tNumRow; i++)
    {
        for(size_t j =0; j<tNumCol; j++)
        {
            tReindexedMat(i,j) = aMatToReindex(aMatRow,aIndexMap(i,j));
        }
    }

    return tReindexedMat;
}

#endif /* PROJECTS_LINALG_SRC_FN_REINDEX_MAT_HPP_ */

