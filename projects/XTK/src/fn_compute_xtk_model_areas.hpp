/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_compute_xtk_model_areas.hpp
 *
 */

#ifndef PROJECTS_XTK_SRC_XTK_FN_COMPUTE_XTK_MODEL_VOLUMES_HPP_
#define PROJECTS_XTK_SRC_XTK_FN_COMPUTE_XTK_MODEL_VOLUMES_HPP_

#include "cl_Matrix.hpp"
#include "cl_XTK_Background_Mesh.hpp"
#include "cl_XTK_Cut_Mesh.hpp"
#include "cl_XTK_Model.hpp"
#include "fn_tet_volume.hpp"
#include "fn_hex_8_volume.hpp"
namespace xtk
{

    /*
     * Given a global coordinate list (so we can just access rather than copy every time we encounter a node)
     * Compute the volume of all background element which are not intersected by the geometry
     *
     * @param[in] aCoords - Node coordinates ordered by index
     * @param[in] aModel  - An XTK Model
     */
    inline moris::real
    compute_non_intersected_parent_element_volume_by_phase(
            moris::moris_index                     aPhaseIndex,
            moris::Matrix< moris::DDRMat > const & aNodeCoordinates,
            Model const &                          aXTKModel )
    {
        // Proc Rank
        moris_index tParRank = par_rank();

        // Get a reference to the XTK Mesh from the Model
        Background_Mesh const & tXTKBMesh = aXTKModel.get_background_mesh();

        // Get the underlying background mesh data
        moris::mtk::Mesh const & tBMMeshData = tXTKBMesh.get_mesh_data();

        // Get the non intersected background elements
        moris::Matrix< moris::IndexMat > tUnintersectedElements = tXTKBMesh.get_all_non_intersected_elements_loc_inds();

        // Determine parent element topology (note: this assumes a uniform background mesh)
        enum CellTopology tParentTopo = tXTKBMesh.get_parent_cell_topology();

        moris::real tVolume = 0;
        for ( size_t i = 0; i < tUnintersectedElements.numel(); i++ )
        {
            // Get the nodes connected to this element
            if ( tXTKBMesh.get_element_phase_index( tUnintersectedElements( i ) ) == aPhaseIndex &&    //
                    tBMMeshData.get_entity_owner( tUnintersectedElements( i ), moris::EntityRank::ELEMENT ) == tParRank )
            {

                moris::Matrix< moris::IndexMat > tElementToNode =
                        tBMMeshData.get_entity_connected_to_entity_loc_inds( i, moris::EntityRank::ELEMENT, moris::EntityRank::NODE );

                if ( tParentTopo == CellTopology::HEX8 )
                {
                    tVolume += compute_hex_8_volume( aNodeCoordinates, tElementToNode );
                }

                else if ( tParentTopo == CellTopology::TET4 )
                {
                    tVolume += vol_tetrahedron( aNodeCoordinates, tElementToNode );
                }

                else
                {
                    MORIS_ERROR( 0, "ELEMENT TYPE VOLUME CALCULATION NOT IMPLEMENTED" );
                }
            }
        }

        return tVolume;
    }
}    // namespace xtk

#endif /* PROJECTS_XTK_SRC_XTK_FN_COMPUTE_XTK_MODEL_VOLUMES_HPP_ */
