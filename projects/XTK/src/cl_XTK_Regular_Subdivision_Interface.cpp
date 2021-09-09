#include "cl_XTK_Decomposition_Algorithm.hpp"
#include "cl_XTK_Regular_Subdivision_Interface.hpp"
#include "cl_XTK_Integration_Mesh_Generator.hpp"
#include "cl_XTK_Cut_Integration_Mesh.hpp"
#include "cl_MTK_Cell_Info_Factory.hpp"
#include "cl_MTK_Cell_Info.hpp"

namespace xtk
{
    bool 
    Regular_Subdivision_Interface::has_geometric_independent_vertices() const
    {
        return true;
    }

    void 
    Regular_Subdivision_Interface::perform_impl_vertex_requests(
        Integration_Mesh_Generation_Data*  aMeshGenerationData,
        Decomposition_Data*                 aDecompositionData,
        Cut_Integration_Mesh*               aCutIntegrationMesh,
        moris::mtk::Mesh*                   aBackgroundMesh,
        Integration_Mesh_Generator*         aMeshGenerator)
    {
        aDecompositionData->mHasSecondaryIdentifier = false;

        // allocate data in decomposition data 
        aDecompositionData->tCMNewNodeLoc        = Cell<Cell<moris_index>>(aMeshGenerationData->mIntersectedBackgroundCellIndexToChildMeshIndex.size(),this->get_num_new_nodes());
        aDecompositionData->tCMNewNodeParamCoord = Cell<Cell<Matrix<DDRMat>>>(aMeshGenerationData->mIntersectedBackgroundCellIndexToChildMeshIndex.size(),
                                                                        Cell<Matrix<DDRMat>>(this->get_num_new_nodes(),
                                                                        Matrix<DDRMat>(1,this->get_parametric_dimension())));
        
        // allocate a data struct to pass into functions - stores the matrix but doesnt require the method to store data
        Regular_Subdivision_Interface_Data tRegSubInterfaceData;

        // get the new nodes on each entity rank
        tRegSubInterfaceData.mNewNodesOnEdges = this->get_new_node_on_parent_edge();
        tRegSubInterfaceData.mNewNodesOnFaces = this->get_new_node_on_parent_face();
        tRegSubInterfaceData.mNewNodesOnCells = this->get_new_node_in_parent_cell();

        // get the new node edge,face,cell ordinal
        tRegSubInterfaceData.mNewNodesOnEdgesOrd = this->get_new_node_on_parent_edge_edge_ordinal();
        tRegSubInterfaceData.mNewNodesOnFacesOrd = this->get_new_node_on_parent_face_face_ordinal();
        tRegSubInterfaceData.mNewNodesOnCellsOrd = this->get_new_node_on_parent_cell_cell_ordinal();

        tRegSubInterfaceData.mNewNodeXi = this->get_new_vertex_parametric_coordinates_wrt_parent();

        // iterate through all intersected background cells and make vertex requests
        for (auto& it: aMeshGenerationData->mIntersectedBackgroundCellIndexToChildMeshIndex) 
        {
            std::shared_ptr<Child_Mesh_Experimental> tChildMesh = aCutIntegrationMesh->get_child_mesh(it.second);  

            // make new vertex requests
            this->make_new_vertex_requests(tChildMesh.get(), this, &tRegSubInterfaceData, aBackgroundMesh, aDecompositionData);

        }

        // handle the requests
        aDecompositionData->tSecondaryIdentifiers = Cell<moris_index>(aDecompositionData->tNewNodeParentIndex.size(), MORIS_INDEX_MAX);

    }

    void
    Regular_Subdivision_Interface::perform_impl_generate_mesh(
                               Integration_Mesh_Generation_Data* aMeshGenerationData,
                               Decomposition_Data*               aDecompositionData,
                               Cut_Integration_Mesh*             aCutIntegrationMesh,
                               moris::mtk::Mesh*                 aBackgroundMesh,
                               Integration_Mesh_Generator*       aMeshGenerator)
    {
        // how many cells am I constructing?
        mNumNewCells = aMeshGenerationData->tNumChildMeshes * this->get_num_ig_cells();

        // new cell info
        moris::mtk::Cell_Info_Factory tFactory;
        std::shared_ptr<moris::mtk::Cell_Info> tIgCellInfo = tFactory.create_cell_info_sp(this->get_ig_cell_topology());

        // this is going to be the cell infor pointer for all new cells
        mNewCellCellInfo = moris::Cell<std::shared_ptr<moris::mtk::Cell_Info>>(mNumNewCells,tIgCellInfo);

        // number of vertices per cell
        moris::uint tVerticesPerCell = tIgCellInfo->get_num_verts();

        // allocate data in the new ig cell data
        mNewCellToVertexConnectivity = moris::Cell<moris::Cell<moris::moris_index>>(mNumNewCells,tVerticesPerCell);
        mNewCellChildMeshIndex       = moris::Cell<moris::moris_index>(mNumNewCells);

        // for this method we are not going to replace any cells, this is because the cell we would replace corresponds to the background mesh cell
        mNewCellCellIndexToReplace   = moris::Cell<moris::moris_index>(mNumNewCells,MORIS_INDEX_MAX);

        // get the cell to vertex template
        moris::Cell<moris::Cell<moris::moris_index>> tIgCellToVertexTemplate = this->get_ig_cell_to_vertex_connectivity();

        // populate new cell data
        moris::moris_index tCurrentCellIndex = 0;
        for(moris::moris_index iCM = 0; iCM < aMeshGenerationData->tNumChildMeshes; iCM++)
        {
            std::shared_ptr<Child_Mesh_Experimental> tChildMesh = aCutIntegrationMesh->get_child_mesh(iCM);

            for(moris::moris_index iNewCell = 0; iNewCell < this->get_num_ig_cells(); iNewCell++)
            {
                mNewCellChildMeshIndex(tCurrentCellIndex) = iCM;

                for(moris::uint iV = 0; iV < tVerticesPerCell; iV++)
                {
                    moris_index tNewVertexCMOrdinal = tIgCellToVertexTemplate(iNewCell)(iV);
                    MORIS_ERROR(tNewVertexCMOrdinal < (moris::moris_index)tChildMesh->mIgVerts->size(),"Template ordinal out of bounds" );
                    mNewCellToVertexConnectivity(tCurrentCellIndex)(iV) = tChildMesh->mIgVerts->get_vertex(tNewVertexCMOrdinal)->get_index();
                }

                tCurrentCellIndex++;
            }
        }

    }                               

    enum Decomposition_Algorithm_Type
    Regular_Subdivision_Interface::get_algorithm_type() const 
    {
        return Decomposition_Algorithm_Type::REGULAR_TEMPLATE_NONCONFORMING;
    }

    void
    Regular_Subdivision_Interface::make_new_vertex_requests(
                             Child_Mesh_Experimental* aChildMesh,
                             Regular_Subdivision_Interface* aRegularSubdivisionInterface,
                             Regular_Subdivision_Interface_Data* aRegularSubdivisionInterfaceData,
                             moris::mtk::Mesh*   aBackgroundMesh,
                             Decomposition_Data* aDecompositionData)
    {
        moris::mtk::Cell* tParentCell = aChildMesh->get_parent_cell();

        moris::mtk::Cell_Info  const * tParentCellInfo = tParentCell->get_cell_info();

        // iterate through nodes on edges
        // moris::Matrix<moris::IndexMat> tEdgeIndices(0,0);
        // for(moris::uint iEdge =0 ;iEdge < aRegularSubdivisionInterfaceData->mNewNodesOnEdges.numel(); iEdge++)
        // {
        //     if(iEdge == 0)
        //     {
        //         moris::Matrix<moris::IndexMat> tEdgeIndices = tBackgroundMeshData.get_entity_connected_to_entity_loc_inds(
        //                 aChildMesh->get_parent_element_index(),
        //                 moris::EntityRank::ELEMENT,
        //                 moris::EntityRank::EDGE);
        //     }

        //     // moris_index tEdgeOwner = tMeshData.get_entity_owner(tRegSubInterfaceData->mNewNodesOnEdgesOrd(iEdge), EntityRank::EDGE);


        //     // evaluate the shape functions at this point relative to the background cell
        //     tParentCellInfo->eval_N(aRegularSubdivisionInterfaceData->mNewNodeXi(aRegularSubdivisionInterfaceData->mNewNodesOnEdgesOrd(iEdge)),aRegularSubdivisionInterfaceData->mNXi);

        //     moris::print(aRegularSubdivisionInterfaceData->mNXi,"tRegSubInterfaceData->mNXi");

        //     // Matrix<DDRMat> tNewNodeCoordinates = tParentCellInfo.get_
        // }

        // iterate through nodes on faces and make requests
        moris::Matrix<moris::IndexMat> tFaceIndices(0,0);
        for(moris::uint iFace = 0 ;iFace < aRegularSubdivisionInterfaceData->mNewNodesOnFaces.numel(); iFace++)
        {
            if(iFace == 0)
            {
                tFaceIndices = aBackgroundMesh->get_entity_connected_to_entity_loc_inds(
                        aChildMesh->get_parent_cell()->get_index(),
                        moris::EntityRank::ELEMENT,
                        moris::EntityRank::FACE);
            }

            moris_index tNewNodeFaceOrdinal = aRegularSubdivisionInterfaceData->mNewNodesOnFacesOrd(iFace);
            moris_index tRequestLoc = MORIS_INDEX_MAX;
            bool tRequestExists = aDecompositionData->request_exists(tFaceIndices(tNewNodeFaceOrdinal),EntityRank::FACE,tRequestLoc);

            moris_index tNewNodeTemplateOrd = aRegularSubdivisionInterfaceData->mNewNodesOnFaces(iFace);
            if(!tRequestExists)
            {
                moris_index tOwner = aBackgroundMesh->get_entity_owner(tFaceIndices(tNewNodeFaceOrdinal), EntityRank::FACE);

                // evaluate the shape functions at this point relative to the background cell
                tParentCellInfo->eval_N(aRegularSubdivisionInterfaceData->mNewNodeXi(tNewNodeTemplateOrd),aRegularSubdivisionInterfaceData->mNXi);

                std::shared_ptr<Matrix<DDRMat>> tNewNodeXi = std::make_shared<Matrix<DDRMat>>(aRegularSubdivisionInterfaceData->mNewNodeXi(tNewNodeTemplateOrd));

                Matrix<DDRMat> tNewNodeCoordinates = aRegularSubdivisionInterfaceData->mNXi * tParentCell->get_vertex_coords();

                moris_index tNewNodeIndexInSubdivision = aDecompositionData->register_new_request(tFaceIndices(tNewNodeFaceOrdinal),tOwner,EntityRank::FACE,tNewNodeCoordinates,tParentCell,tNewNodeXi);

                aDecompositionData->tCMNewNodeParamCoord(aChildMesh->get_child_mesh_index())(tNewNodeTemplateOrd) = aRegularSubdivisionInterfaceData->mNewNodeXi(tNewNodeTemplateOrd);
                aDecompositionData->tCMNewNodeLoc(aChildMesh->get_child_mesh_index())(tNewNodeTemplateOrd)        = tNewNodeIndexInSubdivision;
            }

            else
            {
                aDecompositionData->tCMNewNodeLoc(aChildMesh->get_child_mesh_index())(tNewNodeTemplateOrd) = tRequestLoc;
                aDecompositionData->tCMNewNodeParamCoord(aChildMesh->get_child_mesh_index())(tNewNodeTemplateOrd) = aRegularSubdivisionInterfaceData->mNewNodeXi(tNewNodeTemplateOrd);
            }
        }

        moris::Matrix<moris::IndexMat> tElementIndices(0,0);
        for(moris::uint iCell = 0 ; iCell < aRegularSubdivisionInterfaceData->mNewNodesOnCells.numel(); iCell++)
        {
            moris_index tRequestLoc = MORIS_INDEX_MAX;
            bool tRequestExists = aDecompositionData->request_exists(tParentCell->get_index(),EntityRank::ELEMENT,tRequestLoc);

            moris_index tNewNodeTemplateOrd = aRegularSubdivisionInterfaceData->mNewNodesOnCells(iCell);
            if(!tRequestExists)
            {
                moris_index tOwner = aBackgroundMesh->get_entity_owner(tParentCell->get_index(), EntityRank::ELEMENT);

                // evaluate the shape functions at this point relative to the background cell
                tParentCellInfo->eval_N(aRegularSubdivisionInterfaceData->mNewNodeXi(tNewNodeTemplateOrd),aRegularSubdivisionInterfaceData->mNXi);

                std::shared_ptr<Matrix<DDRMat>> tNewNodeXi = std::make_shared<Matrix<DDRMat>>(aRegularSubdivisionInterfaceData->mNewNodeXi(tNewNodeTemplateOrd));

                Matrix<DDRMat> tNewNodeCoordinates = aRegularSubdivisionInterfaceData->mNXi * tParentCell->get_vertex_coords();

                moris_index tNewNodeIndexInSubdivision = aDecompositionData->register_new_request(tParentCell->get_index(),tOwner,EntityRank::ELEMENT,tNewNodeCoordinates,tParentCell,tNewNodeXi);

                aDecompositionData->tCMNewNodeLoc(aChildMesh->get_child_mesh_index())(tNewNodeTemplateOrd)        = tNewNodeIndexInSubdivision;
                aDecompositionData->tCMNewNodeParamCoord(aChildMesh->get_child_mesh_index())(tNewNodeTemplateOrd) = aRegularSubdivisionInterfaceData->mNewNodeXi(tNewNodeTemplateOrd);
            }

            else
            {
                aDecompositionData->tCMNewNodeLoc(aChildMesh->get_child_mesh_index())(tNewNodeTemplateOrd) = tRequestLoc;
                aDecompositionData->tCMNewNodeParamCoord(aChildMesh->get_child_mesh_index())(tNewNodeTemplateOrd) = aRegularSubdivisionInterfaceData->mNewNodeXi(tNewNodeTemplateOrd);
            }
        }
                
    }
}
