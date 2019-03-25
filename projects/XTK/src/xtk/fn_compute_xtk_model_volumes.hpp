/*
 * fn_compute_xtk_mesh_volume.hpp
 *
 *  Created on: Oct 5, 2018
 *      Author: doble
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
moris::real
compute_non_intersected_parent_element_volume_by_phase(moris::moris_index                   aPhaseIndex,
                                                       moris::Matrix<moris::DDRMat> const & aNodeCoordinates,
                                                       Model const &                        aXTKModel)
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
    enum CellTopology tParentTopo = tXTKBMesh.get_XTK_mesh_element_topology();

    moris::real tVolume = 0;
    for(size_t i = 0; i < tUnintersectedElements.numel(); i++)
    {
        // Get the nodes connected to this element
        if(tXTKBMesh.get_element_phase_index(tUnintersectedElements(i)) == aPhaseIndex &&
                tBMMeshData.get_entity_owner(tUnintersectedElements(i),moris::EntityRank::ELEMENT) == tParRank)
        {

            moris::Matrix< moris::IndexMat > tElementToNode
            = tBMMeshData.get_entity_connected_to_entity_loc_inds(i,moris::EntityRank::ELEMENT,moris::EntityRank::NODE);
            if(tParentTopo == CellTopology::HEX8)
            {
                tVolume += compute_hex_8_volume(aNodeCoordinates, tElementToNode);
            }

            else if (tParentTopo == CellTopology::TET4)
            {
                tVolume+=vol_tetrahedron(aNodeCoordinates, tElementToNode);
            }

            else
            {
                MORIS_ERROR(0,"ELEMENT TYPE VOLUME CALCULATION NOT IMPLEMENTED");
            }
        }
    }


    return tVolume;
}

/*
 * Given a global coordinate list (so we can just access rather than copy every time we encounter a node)
 * Compute the volume of all background element which are not intersected by the geometry
 *
 * @param[in] aCoords - Node coordinates ordere by index
 * @param[in] aModel  - An XTK Model
 */
moris::real
compute_child_element_volume_by_phase(moris::moris_index                                     aPhaseIndex,
                                      moris::Matrix<moris::DDRMat> const &                     aNodeCoordinates,
                                      Model const & aXTKModel)
{
    // Proc Rank
    moris_index tParRank = par_rank();

    // Get a reference to the XTK Mesh from the Model
    Background_Mesh const & tXTKBMesh = aXTKModel.get_background_mesh();

    // Get the underlying background mesh data
    moris::mtk::Mesh const & tBMMeshData = tXTKBMesh.get_mesh_data();

    // Get a reference to the XTK Mesh from the Model
    Cut_Mesh const & tCutMesh = aXTKModel.get_cut_mesh();

    moris::real tVolume = 0;
    for(size_t i = 0; i < tCutMesh.get_num_child_meshes(); i++)
    {
        // Get reference to Child Mesh
        Child_Mesh const & tChildMesh = tCutMesh.get_child_mesh(i);

        // Get reference to nodes connected to elements
        moris::Matrix<moris::IndexMat> const & tElementToNode = tChildMesh.get_element_to_node();

        moris::Matrix< moris::IndexMat > const & tElementPhase = tChildMesh.get_element_phase_indices();

        moris::size_t tNumElems = tChildMesh.get_num_entities(EntityRank::ELEMENT);

        if(tBMMeshData.get_entity_owner(tChildMesh.get_parent_element_index(),moris::EntityRank::ELEMENT) == tParRank)
        {
        for(size_t j = 0; j <tNumElems; j++)
        {

            if(tElementPhase(j) == aPhaseIndex )
            {
                moris::Matrix< moris::IndexMat > tElementToNodeCM = tElementToNode.get_row(j);

                tVolume += vol_tetrahedron(aNodeCoordinates, tElementToNodeCM);

            }
        }
        }
    }

    return tVolume;
}


}




#endif /* PROJECTS_XTK_SRC_XTK_FN_COMPUTE_XTK_MODEL_VOLUMES_HPP_ */
