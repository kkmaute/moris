/*
 * fn_compute_interface_surface_area.hpp
 *
 *  Created on: Mar 21, 2019
 *      Author: doble
 */

#ifndef PROJECTS_XTK_SRC_XTK_FN_COMPUTE_INTERFACE_SURFACE_AREA_HPP_
#define PROJECTS_XTK_SRC_XTK_FN_COMPUTE_INTERFACE_SURFACE_AREA_HPP_


#include "cl_Matrix.hpp"
#include "cl_XTK_Background_Mesh.hpp"
#include "cl_XTK_Cut_Mesh.hpp"
#include "cl_XTK_Model.hpp"
#include "fn_tri_area.hpp"
#include "cl_MTK_Tetra4_Connectivity.hpp"

namespace xtk
{
inline
moris::real
compute_interface_surface_area(moris::Matrix<moris::DDRMat> const & aNodeCoordinates,
                               Model const &                        aXTKModel)
{
    moris::Matrix<moris::IdMat> tInterfaceElemIndandSideOrd = aXTKModel.get_cut_mesh().pack_interface_sides_loc_inds();

    Background_Mesh const & tBackgroundMesh = aXTKModel.get_background_mesh();

    moris::real tSurfArea = 0.0;

    moris::Matrix<moris::IndexMat> tNodeMapOnFace = moris::Tetra4_Connectivity::get_node_to_face_map();

    // iterate through tInterface sides
    for(moris::uint i = 0; i<tInterfaceElemIndandSideOrd.n_rows(); i++)
    {
        moris_index tElemInd = tInterfaceElemIndandSideOrd(i,0);
        moris_index tSideOrd = tInterfaceElemIndandSideOrd(i,1);

        if(tBackgroundMesh.get_element_phase_index(tElemInd) == 0)
        {
            // cell reference
            moris::mtk::Cell const & tCell = tBackgroundMesh.get_mtk_cell(tElemInd);

            // vertex indices
            moris::Matrix<IndexMat> tVertexInds = tCell.get_vertex_inds();

            moris::Matrix<moris::DDRMat> tTriCoords(3,3);
            moris::moris_index tIndex0 = tVertexInds( tNodeMapOnFace(tSideOrd,0) );
            moris::moris_index tIndex1 = tVertexInds( tNodeMapOnFace(tSideOrd,1) );
            moris::moris_index tIndex2 = tVertexInds( tNodeMapOnFace(tSideOrd,2) );

            tTriCoords.get_row(0) = aNodeCoordinates.get_row(tIndex0);
            tTriCoords.get_row(1) = aNodeCoordinates.get_row(tIndex1);
            tTriCoords.get_row(2) = aNodeCoordinates.get_row(tIndex2);

            tSurfArea = tSurfArea + area_tri(tTriCoords);
        }
    }

    return tSurfArea;
}
}

#endif /* PROJECTS_XTK_SRC_XTK_FN_COMPUTE_INTERFACE_SURFACE_AREA_HPP_ */
