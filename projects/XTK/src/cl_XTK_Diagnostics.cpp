#include "cl_XTK_Diagnostics.hpp"

#include "cl_XTK_Cut_Integration_Mesh.hpp"
#include "cl_XTK_Model.hpp"
#include "cl_GEN_Geometry_Engine.hpp"
#include "cl_MTK_Vertex.hpp"

using namespace moris;
namespace xtk
{
    bool
    interpolated_coordinate_check(Cut_Integration_Mesh* aCutMesh)
    {
        moris::Cell<std::shared_ptr<Matrix<DDRMat>>>* tCoords = aCutMesh->get_all_vertex_coordinates_loc_inds();

        // Allocate a basis function weight matrix
        moris::Matrix< moris::DDRMat > tBasisWeights(1,8);

        // tolerance for difference between coordinates
        real tTol = 1e-12;

        // iterate through 
        moris::uint tNumIgCellGroups = aCutMesh->get_num_ig_cell_groups();

        // iterate over ig cell groups
        for( moris::uint iIgCellGroup = 0; iIgCellGroup < tNumIgCellGroups; iIgCellGroup++ )
        {
            // get the cell group
            std::shared_ptr<IG_Cell_Group> tIgCellGroup = aCutMesh->get_ig_cell_group((moris_index)iIgCellGroup);

            // get the vertex group
            std::shared_ptr<IG_Vertex_Group> tIgVertexGroup = aCutMesh->get_vertex_group((moris_index)iIgCellGroup);

            // the base cell 
            moris::mtk::Cell* tBaseCell = aCutMesh->get_ig_cell_group_parent_cell((moris_index)iIgCellGroup);

            moris::Matrix<moris::DDRMat> tBaseCellCoords = tBaseCell->get_vertex_coords();

            moris::mtk::Cell_Info const * tCellInfo = tBaseCell->get_cell_info();

            // iterate through the IgCellGroup
            for ( moris::uint iCell = 0; iCell < tIgCellGroup->mIgCellGroup.size(); iCell++)
            {
                moris::mtk::Cell* tIgCell = tIgCellGroup->mIgCellGroup(iCell);

                moris::Cell<moris::mtk::Vertex *> tVertices =  tIgCell->get_vertex_pointers();

                // iterate through vertices attached to the cell
                for( moris::uint iV = 0; iV < tVertices.size(); iV++ ) 
                {
                    // local coordinate of this vertex wrt the current cell group
                    std::shared_ptr<Matrix<DDRMat>> tVertexLocalCoords = tIgVertexGroup->get_vertex_local_coords(tVertices(iV)->get_index());

                    // evaluate the basis function
                    tCellInfo->eval_N(*tVertexLocalCoords,tBasisWeights);

                    // Evaluate the nodes global coordinate from the basis weights
                    moris::Matrix<moris::DDRMat> tInterpNodeCoord = tBasisWeights*tBaseCellCoords;

                    // node index
                    moris_index tNodeIndex = tVertices(iV)->get_index();

                    if(norm(tInterpNodeCoord - *((*tCoords)(tNodeIndex))) > tTol)
                    {
                        return false;
                    }
                }
            }
        }
        return true;
    }

    bool
    verify_interface_vertices(
        xtk::Model* aModel, 
        moris::Matrix<moris::DDRMat> const & aIsocontourThreshold,
        moris::Matrix<moris::DDRMat> const & aIsocontourTolerance)
    {
        // get the cut integration mesh
        Cut_Integration_Mesh* tCutMesh  = aModel->get_cut_integration_mesh();

        // get the geometry engine
        moris::ge::Geometry_Engine* tGeomEngine = aModel->get_geom_engine();

        // number of vertices
        moris::uint tNumVertices = tCutMesh->get_num_entities(EntityRank::NODE,0);

        // iterate through vertices in the integration mesh
        for(moris::uint iV = 0; iV < tNumVertices; iV++)
        {
            // iterate through geometries
            for(moris::uint iGeom = 0; iGeom < tGeomEngine->get_num_geometries(); iGeom++)
            {
                moris::mtk::Vertex* tVertex = &tCutMesh->get_mtk_vertex((moris_index)iV);
                if(tGeomEngine->is_interface_vertex(tVertex->get_index(), iGeom))
                {
                    moris::real tGeomVal = tGeomEngine->get_field_value(iGeom,tVertex->get_index(),tVertex->get_coords());
                    if( std::abs(tGeomVal - aIsocontourThreshold(iGeom)) > aIsocontourTolerance(iGeom) )
                    {
                        return false;
                    }

                }

            }
        }
        return true;
        
    }
}