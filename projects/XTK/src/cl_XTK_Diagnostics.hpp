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
namespace xtk
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

/**
 * @brief Checks that interface vertices are sufficiently close to the interface. Note, due to geometric approximation error this function may not work for functions very different from SDF
 *
 * @param aModel - XTK model
 * @param aIsocontourThreshold - Isocontour threshold
 * @param aIsocontourTolerance - Isocontour tolerance
 * @return true - all interface vertices are on the interface
 * @return true - all interface vertices are NOT on the interface
 */
bool
verify_interface_vertices(
    xtk::Model*                           aModel,
    moris::Matrix< moris::DDRMat > const& aIsocontourThreshold,
    moris::Matrix< moris::DDRMat > const& aIsocontourTolerance );

bool
check_vertices(
    xtk::Model*                                                         aModel,
    moris::Vector< moris::moris_index > const&                            aGoldNumVerts,
    std::unordered_map< moris::moris_index, moris::moris_index > const& aGoldVertexMap,
    moris::Matrix< moris::DDRMat > const&                               aGoldVertexCoords,
    moris::real                                                         aTolerance = 1e-12 );

bool
check_cells(
    xtk::Model*                                                         aModel,
    moris::Vector< moris::moris_index > const&                            aGoldNumCells,
    std::unordered_map< moris::moris_index, moris::moris_index > const& aGoldCellMap,
    moris::Matrix< moris::IndexMat > const&                             aGoldCellConn );
}// namespace xtk

#endif /* cl_XTK_Diagnostics.hpp */
