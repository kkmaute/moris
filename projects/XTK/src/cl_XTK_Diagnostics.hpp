/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Diagnostics.hpp
 *
 */

#ifndef SRC_cl_XTK_Diagnostics
#define SRC_cl_XTK_Diagnostics

#include "cl_Matrix.hpp"
#include "cl_Vector.hpp"
#include "moris_typedefs.hpp"
#include <unordered_map>
namespace moris::xtk
{
    class Cut_Integration_Mesh;
    class Model;

    /**
     * @brief Collection of diagnostic free functions for XTK
     */

    /**
     * @brief Checks if the interpolation to a coordinate is as expected
     *
     * @param aCutMesh
     * @return true - all coordinates of integration vertices evaluate to their expected value
     * @return false
     */
    bool
    interpolated_coordinate_check( Cut_Integration_Mesh* aCutMesh );

    bool
    check_vertices(
            xtk::Model*                                                          aModel,
            Vector< moris::moris_index > const &                                 aGoldNumVerts,
            std::unordered_map< moris::moris_index, moris::moris_index > const & aGoldVertexMap,
            Matrix< DDRMat > const &                                             aGoldVertexCoords,
            moris::real                                                          aTolerance = 1e-12 );

    bool
    check_cells(
            xtk::Model*                                                          aModel,
            Vector< moris::moris_index > const &                                 aGoldNumCells,
            std::unordered_map< moris::moris_index, moris::moris_index > const & aGoldCellMap,
            Matrix< IndexMat > const &                                           aGoldCellConn );
}    // namespace moris::xtk

#endif /* cl_XTK_Diagnostics.hpp */
