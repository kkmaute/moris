/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_MTK_QuadraturePointMapper.cpp
 *
 */
#include "cl_MTK_QuadraturePointMapper.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_Vector.hpp"
#include "cl_MTK_Side_Set.hpp"
#include "fn_assert.hpp"

namespace moris::mtk
{
    QuadraturePointMapper::QuadraturePointMapper(
            Integration_Mesh                                      *aIGMesh,
            Vector< Side_Set const * >                            &aSideSets,
            Vector< std::pair< moris_index, moris_index > > const &aCandidatePairs )
            : mIGMesh( aIGMesh )
            , mSideSets( aSideSets )
            , mCandidatePairs( aCandidatePairs )
    {
        MORIS_ASSERT( aSideSets.size() > 0, "The ContactMapper needs at least one SideSet!" );
    }
}    // namespace moris::mtk
