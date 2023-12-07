//
// Created by frank on 11/27/23.
//

#include "cl_MTK_ContactMapper.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_Cell.hpp"
#include "cl_MTK_Side_Set.hpp"
#include "fn_assert.hpp"

namespace moris::mtk
{
    ContactMapper::ContactMapper(
            Integration_Mesh                                           *aIGMesh,
            moris::Cell< Side_Set * >                                  &aSideSets,
            moris::Cell< std::pair< moris_index, moris_index > > const &aCandidatePairs )
            : mIntegrationMeshes( aIGMesh )
            , mSideSets( aSideSets )
            , mCandidatePairs( aCandidatePairs )
    {
        MORIS_ASSERT( aSideSets.size() > 0, "The ContactMapper needs at least one SideSet!" );
    }

    std::ostream &operator<<( std::ostream &aOs, ContactMapper::MappingResult const &aPoint )
    {
        //        aOs << "Point: (";
        //        uint tNumElements = aPoint.parametric_coordinate.numel();
        //        for ( uint tIndex = 0; tIndex < tNumElements; ++tIndex )
        //        {
        //            aOs << aPoint.parametric_coordinate( tIndex, 0 );
        //            if ( tIndex < tNumElements - 1 ) aOs << ", ";
        //        }
        //        aOs << ")"
        //            << ", SideSet-Index: " << aPoint.sideset_index
        //            << ", Cluster-Index: " << aPoint.cluster_index
        //            << ", Cell-Index: " << aPoint.cell_index
        //            << ", Side-Ordinal: " << aPoint.side_ordinal;
        return aOs;
    }
}    // namespace moris::mtk