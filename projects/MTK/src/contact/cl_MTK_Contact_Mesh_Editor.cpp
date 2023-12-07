/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_MTK_Contact_Mesh_Editor.cpp
 *
 */

#include "cl_MTK_Contact_Mesh_Editor.hpp"

namespace moris::mtk
{
    void Contact_Mesh_Editor::update_nonconformal_side_sets()
    {
        Matrix< DDRMat > tQPoints;
        mIntegrationRule.get_points( tQPoints );
        tQPoints = tQPoints.get_row( 0 ).eval();

        Matrix< DDRMat > tQWeights;
        mIntegrationRule.get_weights( tQWeights );

//        moris_index      tSourceCell        = 30;
//        moris_index      tSourceCluster     = 7;
//        moris_index      tSourceSideOrdinal = 0;
//        Matrix< DDRMat > tParametricCoord   = { { 0 } };
//
//        ContactMapper::Point tSource{ tParametricCoord, 0, tSourceCluster, tSourceCell, tSourceSideOrdinal };


        mPointMapper.map( tQPoints );


//        auto tTarget = mPointMapper.query(tSource);

//        std::cout << tTarget << std::endl;

    }

    void Contact_Mesh_Editor::update_displacements( Matrix< DDRMat >& aDisplacements )
    {
    }
}    // namespace moris::mtk
