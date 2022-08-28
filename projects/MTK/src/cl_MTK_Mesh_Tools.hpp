/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Mesh_Tools.hpp
 *
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_MESH_TOOLS_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_MESH_TOOLS_HPP_

#include "cl_Mesh_Enums.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Mesh_Core.hpp"
#include "cl_Matrix.hpp"

namespace moris
{
namespace mtk
{

    /*
     * Takes a matrix of entity indices and entity rank, and returns
     * a matrix of entity glb ids.
     */
    inline
    Matrix< IdMat >
    convert_entity_indices_to_ids(
                    Matrix< IndexMat > aEntityIndices,
                    EntityRank         aEntityRank,
                    Mesh*              aMesh)
    {
        Matrix< IdMat > tEntityIds ( aEntityIndices.n_rows(), aEntityIndices.n_cols());

        for( uint  i = 0; i < aEntityIndices.n_rows(); i++)
        {
            for( uint  j = 0; j < aEntityIndices.n_cols(); j++)
            {
                tEntityIds(i,j) = aMesh->get_glb_entity_id_from_entity_loc_index(aEntityIndices(i,j),aEntityRank);
            }
        }

        return tEntityIds;
    }

}
}

#endif /* PROJECTS_MTK_SRC_CL_MTK_MESH_TOOLS_HPP_ */

