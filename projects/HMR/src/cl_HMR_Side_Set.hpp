/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Side_Set.hpp
 *
 */

#ifndef PROJECTS_HMR_SRC_DOC_CL_HMR_SIDE_SET_HPP_
#define PROJECTS_HMR_SRC_DOC_CL_HMR_SIDE_SET_HPP_

#include <string>

#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_MTK_Side_Sets_Info.hpp"

namespace moris::hmr
{
    struct Side_Set
    {
        Matrix<IdMat>        mElemIdsAndSideOrds;
        Matrix<IndexMat>     mElemIndices;
        mtk::MtkSideSetInfo  mInfo;
//-------------------------------------------------------------------------------
    };

}

#endif /* PROJECTS_HMR_SRC_DOC_CL_HMR_SIDE_SET_HPP_ */

