/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Side_Set.hpp
 *
 */

#pragma once

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
