/*
 * fn_XTK_IO_Utilities.cpp
 *
 *  Created on: Jun 20, 2019
 *      Author: doble
 */



#include "fn_XTK_IO_Utilities.hpp"
#include "cl_XTK_Cut_Mesh.hpp"
#include "cl_XTK_Background_Mesh.hpp"
#include "cl_MTK_Mesh_Core.hpp"

#include <iomanip>      // std::setw

namespace xtk
{
void
print_cut_mesh(Model const & aModel)
{
    uint tWidth = 7;

    Cut_Mesh         const & tCutMesh        = aModel.get_cut_mesh();
    Background_Mesh  const & tBackgroundMesh = aModel.get_background_mesh();
    moris::mtk::Mesh const & tMeshData       = tBackgroundMesh.get_mesh_data();

    std::cout<<"\nCut Mesh Details:"<<std::endl;
    std::cout<<"Num Background Cells With Children Cells: "<<std::right<<std::setw(tWidth)<<tCutMesh.get_num_child_meshes()<<std::endl;
    std::cout<<"Background Cells with Children Cells:     ";
    for(moris::uint  i = 0; i <tCutMesh.get_num_child_meshes(); i++)
    {
        Child_Mesh const & tCM = tCutMesh.get_child_mesh((moris_index)i);
        moris_index tParentCellIndex = tCM.get_parent_element_index();
        moris_id tParentCellId = tMeshData.get_glb_entity_id_from_entity_loc_index(tParentCellIndex,EntityRank::ELEMENT);

        std::cout<<std::right<<std::setw(tWidth)<<tParentCellId;
    }
    std::cout<<"\n";

    for(moris::uint  i = 0; i <tCutMesh.get_num_child_meshes(); i++)
    {
        Child_Mesh const & tCM = tCutMesh.get_child_mesh((moris_index)i);
        moris_index tParentCellIndex = tCM.get_parent_element_index();
        moris_id tParentCellId = tMeshData.get_glb_entity_id_from_entity_loc_index(tParentCellIndex,EntityRank::ELEMENT);
        Matrix<IdMat> const & tNodeInds     = tCM.get_node_indices();
        Matrix<IdMat> const & tNodeIds      = tCM.get_node_ids();
        Matrix<IdMat> const & tCellIds      = tCM.get_element_ids();
        Matrix<DDRMat> const & tParamCoords = tCM.get_parametric_coordinates();
        std::cout<<"\nChild Mesh in Background Cell Id: "<<std::right<<std::setw(tWidth)<<tParentCellId;
        std::cout<<"\n     Num Vertices: "<<std::right<<std::setw(tWidth)<<tNodeInds.numel();
        std::cout<<"\n     Num Edges:    "<<std::right<<std::setw(tWidth)<<tCM.get_num_entities(EntityRank::EDGE);
        std::cout<<"\n     Num Edges:    "<<std::right<<std::setw(tWidth)<<tCM.get_num_entities(EntityRank::FACE);
        std::cout<<"\n     Num Cells:    "<<std::right<<std::setw(tWidth)<<tCM.get_num_entities(EntityRank::ELEMENT);

        std::cout<<"\n    Vertex Info:"<< std::endl;
        std::cout<<"\n     Vertex Parametric Coordinates:"<< std::endl;
        for(moris::uint k = 0; k < tNodeIds.numel(); k++)
        {
            std::cout<<"    "<<std::right<<std::setw(tWidth-4)<<" Vert Id: "<<std::right<<std::setw(tWidth)<<tNodeIds(k)<<" | Vert Ind:"<<std::right<<std::setw(tWidth)<<tNodeInds(k)<<" | CM Vert Ind:"<<std::right<<std::setw(tWidth)<<k<<std::endl;
        }

        std::cout<<"\n     Vertex Parametric Coordinates:"<< std::endl;
        for(moris::uint k = 0; k < tNodeIds.numel(); k++)
        {
            std::cout<<"    "<<std::right<<std::setw(tWidth-4)<<" Vert Id: "<<std::right<<std::setw(tWidth)<<tNodeIds(k)<<" | ";
            for(moris::uint j = 0; j <tParamCoords.n_cols(); j++)
            {
                std::cout<<std::scientific<<std::right<<std::setw(14)<<tParamCoords(k,j)<< "   ";
            }
            std::cout<<std::endl;
        }

        std::cout<<"\n     Vertex Physical Coordinates:"<< std::endl;
        Matrix<DDRMat> tVertCoord = tBackgroundMesh.get_selected_node_coordinates_loc_inds(tNodeInds);
        for(moris::uint k = 0; k < tNodeIds.numel(); k++)
        {
            std::cout<<"    "<<std::right<<std::setw(tWidth-4)<<" Vert Id: "<<std::right<<std::setw(tWidth)<<tNodeIds(k)<<" | ";

            for(moris::uint j = 0; j <tVertCoord.n_cols(); j++)
            {
                std::cout<<std::scientific<<std::right<<std::setw(14)<<tVertCoord(k,j)<< "   ";
            }
            std::cout<<std::endl;
        }

        std::cout<<"\n     Cell Info:"<< std::endl;
        for(moris::uint i = 0; i < tCellIds.numel(); i++)
        {
            std::cout<<"     Cell Id: "<<std::right<<std::setw(5)<<tCellIds(i)<<" | ";

            Matrix<IdMat> tCellToNode = tCM.get_element_to_node_glob_ids((moris_index)i);
            for(moris::uint j = 0; j < tCellToNode.numel(); j++)
            {
                std::cout<<std::right<<std::setw(10)<<(tCellToNode)(j)<< "   ";
            }

            std::cout<<std::endl;
        }

    }



}
}
