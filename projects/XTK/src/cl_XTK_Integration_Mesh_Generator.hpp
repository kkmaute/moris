#ifndef MORIS_CL_XTK_Integration_Mesh_Generator_HPP_
#define MORIS_CL_XTK_Integration_Mesh_Generator_HPP_


#include "cl_XTK_Model.hpp"
#include "cl_GEN_Geometry_Engine.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_GEN_Geometric_Query_Interface.hpp"
#include "cl_Mesh_Enums.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Cell_Info.hpp"
#include "cl_XTK_Child_Mesh.hpp"
#include "cl_MTK_Mesh_Core.hpp"
#include "cl_XTK_Background_Mesh.hpp"
#include "cl_XTK_Decomposition_Data.hpp"
#include "cl_XTK_Cut_Integration_Mesh.hpp"
#include "cl_MTK_Cell_Info_Factory.hpp"
#include "cl_XTK_Cell_No_CM.hpp"


namespace xtk
{


struct Integration_Mesh_Generation_Data
{
    // number of child meshes
    moris_index tNumChildMeshes = 0;
    // outer cell geometry index (if the geometry is inactive the cell is empty)
    // inner cell active child mesh index by cut mesh index
    moris::Cell<moris::Cell<moris_index>> mIntersectedBackgroundCellIndex;

    // this tells me which geometry indices are associated with each intersected background cell. (transpose of mIntersectedBackgroundCellIndex)
    moris::Cell<moris::Cell<moris_index>> mBackgroundCellGeometryIndices;

    // this maps from the background cell index to the child mesh index
    std::unordered_map<moris_index,moris_index> mIntersectedBackgroundCellIndexToChildMeshIndex;


    // continuous matrix of coordinates to use (ensuring no copies here as this is quite large)
    // here to keep in scope until end of integration mesh generation
    Matrix<DDRMat> tVertexCoordinates;
};




class Geometric_Query_XTK : public moris::ge::Geometric_Query_Interface
{
    private:
    // tell the geometry engine what you are interested
    enum moris::ge::Query_Type mQueryType;

    // associated with a given child mesh
    enum EntityRank  mQueryEntityRank;
    Matrix<IndexMat> mQueryEntityToVertices;

    // child mesh pointer
    Child_Mesh* mChildMesh;

    // parent cell is the query object when checking for intersection
    moris::mtk::Cell* mQueryParentCell;

    // parent info 
    moris::mtk::Cell* mParentCell;

    // row-based and indexed based on the indexing in mQueryEntityToVertices.
    moris::Matrix<moris::DDRMat>* mQueryEntityIndexedCoordinates;    

    // geometric index
    moris_index mGeometricIndex;

    public:
    Geometric_Query_XTK():mChildMesh(nullptr),mQueryParentCell(nullptr),mParentCell(nullptr),mGeometricIndex(MORIS_INDEX_MAX)
    {

    }
    ~Geometric_Query_XTK(){}

    void
    set_query_type(enum moris::ge::Query_Type aQueryType)
    {
        mQueryType = aQueryType;
    }

    void
    set_parent_cell(moris::mtk::Cell* aParentCell)
    {
        mParentCell = aParentCell;
    }

    void
    set_query_parent_cell(moris::mtk::Cell* aQueryParentCell)
    {
        mQueryParentCell = aQueryParentCell;
        mQueryEntityToVertices = aQueryParentCell->get_vertex_inds();
        mChildMesh = nullptr;
    }

    void
    set_query_child_mesh( 
        Child_Mesh* aQueryChildMesh, 
        moris_index aEdgeIndex)
    {
        mQueryParentCell       = nullptr;
        mChildMesh             = aQueryChildMesh;
        mQueryEntityToVertices = mChildMesh->get_edge_to_node().get_row(aEdgeIndex);
    }

    void
    set_coordinates_matrix(Matrix<DDRMat> * aCoordinates)
    {
        mQueryEntityIndexedCoordinates = aCoordinates;
    }

    void
    set_geometric_index(moris_index aGeometricIndex)
    {
        mGeometricIndex = aGeometricIndex;
    }
    
    enum moris::ge::Query_Type get_query_type() const
    { 
        return mQueryType;
    }

    moris_index 
    get_geomtric_index() const
    {
        return mGeometricIndex;
    }

    enum EntityRank 
    get_query_entity_rank() const
    {
        return mQueryEntityRank;
    }

    Matrix<IndexMat> const & 
    get_query_entity_to_vertex_connectivity() const
    {
        return mQueryEntityToVertices;
    }

    Matrix<DDRMat> const * get_query_indexed_coordinates() const
    {
        return mQueryEntityIndexedCoordinates;
    }

    Matrix<DDRMat> 
    get_vertex_local_coord_wrt_parent_entity( moris_index aVertexIndex ) const
    {
        MORIS_ERROR(mChildMesh != nullptr ,"get_vertex_local_coord_wrt_parent_entity Only supported on child meshes queries.");
        return mChildMesh->get_parametric_coordinates(aVertexIndex);
    }

    enum EntityRank 
    get_query_parent_entity_rank() const 
    {
        return EntityRank::ELEMENT;
    }

    Matrix<IndexMat>
    get_query_parent_entity_connectivity() const 
    {
        return mParentCell->get_vertex_inds();
    }   

    moris_index 
    get_query_parent_entity_id() const 
    {
        return mParentCell->get_id();
    }

    
    moris_index 
    max_query_entity_intersection() const 
    {
        return 1;
    }
};


class Decomposition_Algorithm;
class Integration_Mesh_Generator
{
    private:
        Model*                         mXTKModel;
        moris::ge::Geometry_Engine*    mGeometryEngine;
        moris::Matrix<moris::IndexMat> mActiveGeometries;
        bool                           mAllActiveGeometries;
        enum Subdivision_Method        mRegularSubdivision; // Required to have a inheritance structure
        enum Subdivision_Method        mConformalSubdivision; // Recursive no inheritance

        

    public:
        Integration_Mesh_Generator( xtk::Model*                    aXTKModelPtr,
                                    Cell<enum Subdivision_Method>  aMethods,
                                    moris::Matrix<moris::IndexMat> aActiveGeometries);
        ~Integration_Mesh_Generator();

        bool
        perform();

        bool
        determine_intersected_background_cells( Integration_Mesh_Generation_Data & aMeshGenerationData, 
                                                Cut_Integration_Mesh*              aCutIntegrationMesh,
                                                moris::mtk::Mesh*                  aBackgroundMesh )
        {
            // set the global coordinate list to use throughout the mesh generation
            aMeshGenerationData.tVertexCoordinates = aCutIntegrationMesh->get_all_vertex_coordinates_loc_inds();
            
            uint tNumGeometries = mAllActiveGeometries? mXTKModel->mGeometryEngine->get_num_geometries() : mActiveGeometries.numel();

            uint tNumCells = aBackgroundMesh->get_num_elems();

            aMeshGenerationData.mIntersectedBackgroundCellIndex.resize(tNumGeometries);
            aMeshGenerationData.mIntersectedBackgroundCellIndex.reserve(tNumGeometries*tNumCells);

            aMeshGenerationData.mBackgroundCellGeometryIndices.resize(tNumCells);
            aMeshGenerationData.mBackgroundCellGeometryIndices.reserve(tNumGeometries*tNumCells);

            // Initialize geometric query
            Geometric_Query_XTK tGeometricQuery;

            // say I am just interested in a yes or no answer
            tGeometricQuery.set_query_type(moris::ge::Query_Type::INTERSECTION_NO_LOCATION);

            // large coord matrix that I want to keep in scope for a long time avoid copying coordinate all the time.
            tGeometricQuery.set_coordinates_matrix(&aMeshGenerationData.tVertexCoordinates);

            for(moris::uint iCell = 0; iCell < tNumCells; iCell++)
            {
                // setup geometric query with this current cell information
                tGeometricQuery.set_parent_cell(&aBackgroundMesh->get_mtk_cell((moris_index) iCell));
                tGeometricQuery.set_query_parent_cell(&aBackgroundMesh->get_mtk_cell((moris_index) iCell));

                for(moris::size_t iGeom = 0; iGeom<tNumGeometries; iGeom++)
                {
                    // current index for this geometry
                    moris_index tGeometryIndex = mAllActiveGeometries? iGeom : mActiveGeometries(iGeom);

                    // tell the query which geometric index we are working on
                    tGeometricQuery.set_geometric_index(tGeometryIndex);

                    if(mXTKModel->get_geom_engine()->geometric_query(&tGeometricQuery))
                    {
                        // add background cell to the list for iGEOM
                        aMeshGenerationData.mIntersectedBackgroundCellIndex(iGeom).push_back(iCell);

                        // add the geometry to the background cell list
                        aMeshGenerationData.mBackgroundCellGeometryIndices(iCell).push_back(iGeom);

                        // if this one is intersected for the first time give it a child mesh index
                        if(aMeshGenerationData.mIntersectedBackgroundCellIndexToChildMeshIndex.find(iCell) == aMeshGenerationData.mIntersectedBackgroundCellIndexToChildMeshIndex.end())
                        {
                            aMeshGenerationData.mIntersectedBackgroundCellIndexToChildMeshIndex[iCell] = aMeshGenerationData.tNumChildMeshes;
                            aMeshGenerationData.tNumChildMeshes++;
                        }
                    }
                }
            }

            // remove the excess space
            aMeshGenerationData.mIntersectedBackgroundCellIndex.shrink_to_fit();
            aMeshGenerationData.mBackgroundCellGeometryIndices.shrink_to_fit();

            return true;
        }

        bool
        allocate_child_meshes(Integration_Mesh_Generation_Data & aMeshGenerationData, 
                              Cut_Integration_Mesh* aCutIntegrationMesh,
                              moris::mtk::Mesh*     aBackgroundMesh  )
        {
            
            aCutIntegrationMesh->mChildMeshes.resize(aMeshGenerationData.tNumChildMeshes);
            aCutIntegrationMesh->mIntegrationCellGroups.resize(aMeshGenerationData.tNumChildMeshes,nullptr);
            aCutIntegrationMesh->mIntegrationVertexGroups.resize(aMeshGenerationData.tNumChildMeshes,nullptr);
            aCutIntegrationMesh->mIntegrationCellGroupsParentCell.resize(aMeshGenerationData.tNumChildMeshes,nullptr);

            // create the child meshes
            for (auto& it: aMeshGenerationData.mIntersectedBackgroundCellIndexToChildMeshIndex) 
            {
                moris_index tCMIndex = it.second;
                moris::mtk::Cell* tParentCell = &aBackgroundMesh->get_mtk_cell(it.first);
                aCutIntegrationMesh->mChildMeshes(tCMIndex)                     = std::make_shared<Child_Mesh_Experimental>();
                aCutIntegrationMesh->mIntegrationCellGroups(tCMIndex)           = std::make_shared<IG_Cell_Group>(0);
                aCutIntegrationMesh->mChildMeshes(tCMIndex)->mIgCells           = aCutIntegrationMesh->mIntegrationCellGroups(tCMIndex);
                aCutIntegrationMesh->mIntegrationCellGroupsParentCell(tCMIndex) = tParentCell;
                aCutIntegrationMesh->mChildMeshes(tCMIndex)->mParentCell        = tParentCell;
                aCutIntegrationMesh->mChildMeshes(tCMIndex)->mChildMeshIndex    = tCMIndex;


                moris_index tNumGeometricVertices = 8;
                moris::Cell<moris::mtk::Vertex*> tParentCellVerts = tParentCell->get_vertex_pointers();
                aCutIntegrationMesh->mIntegrationVertexGroups(tCMIndex) = std::make_shared<IG_Vertex_Group>(tNumGeometricVertices);
                aCutIntegrationMesh->mChildMeshes(tCMIndex)->mIgVerts = aCutIntegrationMesh->mIntegrationVertexGroups(tCMIndex);
                //FIXME: GET GEOMETRIC VERTICES FROM MTK CELL HARDCODED TO HEXFAMILY
                for(moris::moris_index i = 0 ; i < tNumGeometricVertices; i++)
                {
                    aCutIntegrationMesh->mIntegrationVertexGroups(tCMIndex)->mIgVertexGroup(i) = tParentCellVerts(i);
                }
                
            }

            return true;
        }

        bool
        determine_intersected_ig_edges( Integration_Mesh_Generation_Data & aMeshGenerationData,
                                        Cut_Integration_Mesh*              aCutIntegrationMesh,
                                        moris::mtk::Mesh*                  aBackgroundMesh )
        {
            // Generate edge to node connectivity for the entire cut integration mesh
            this->setup_edge_based_connectivity_of_ig_mesh(aMeshGenerationData, aCutIntegrationMesh);

        }

        void
        setup_child_mesh_edge_based_connectivity_of_ig_mesh(
            Integration_Mesh_Generation_Data & aMeshGenerationData,
            Cut_Integration_Mesh*              aCutIntegrationMesh)
        {
            // iterate through child meshes and compute their edge based connectivity
            
            
        }

        void
        commit_new_ig_vertices_to_cut_mesh(
            Integration_Mesh_Generation_Data*  aMeshGenerationData,
            Decomposition_Data *               aDecompositionData,
            Cut_Integration_Mesh *             aCutIntegrationMesh,
            moris::mtk::Mesh *                 aBackgroundMesh)
        {
            // current index
            moris_index tControlledVertexIndex = aCutIntegrationMesh->mControlledIgVerts.size();

            // allocate new vertices
            moris::uint tNumNewIgVertices = aDecompositionData->tNewNodeIndex.size();
            aCutIntegrationMesh->mControlledIgVerts.resize(aCutIntegrationMesh->mControlledIgVerts.size() + tNumNewIgVertices);
            aCutIntegrationMesh->mIntegrationVertices.resize(aCutIntegrationMesh->mIntegrationVertices.size() + tNumNewIgVertices);
            
            // iterate and create new vertices
            for(moris::uint iV = 0; iV < aDecompositionData->tNewNodeId.size(); iV++ )
            {
                // create a controlled vertex (meaning I need to manage memory of it)
                aCutIntegrationMesh->mControlledIgVerts(tControlledVertexIndex) = std::make_shared<moris::mtk::Vertex_XTK>( 
                    aDecompositionData->tNewNodeId(iV),
                    aDecompositionData->tNewNodeIndex(iV),
                    aDecompositionData->tNewNodeCoordinate(iV));
                    
                aCutIntegrationMesh->mIntegrationVertices(aDecompositionData->tNewNodeIndex(iV)) = aCutIntegrationMesh->mControlledIgVerts(tControlledVertexIndex).get();
                tControlledVertexIndex++;
            }
            

            // iterate through child meshes and commit the vertices to their respective vertex groups
            for(moris::uint iCM = 0; iCM < (moris::uint) aMeshGenerationData->tNumChildMeshes; iCM++)
            {
                MORIS_ERROR(aDecompositionData->tCMNewNodeLoc.size() == (moris::uint) aCutIntegrationMesh->mChildMeshes.size(),"Mismatch in child mesh sizes. All child meshes need to be present in the decomposition data");

                // add the vertices to child mesh groups
                moris_index tStartIndex = (moris_index) aCutIntegrationMesh->mIntegrationVertexGroups(iCM)->mIgVertexGroup.size();
                moris_index tNumNewVertices = (moris_index) aDecompositionData->tCMNewNodeLoc(iCM).size();

                // resize the vertices in the group
                aCutIntegrationMesh->mIntegrationVertexGroups(iCM)->mIgVertexGroup.resize(tNumNewVertices + tStartIndex,nullptr);

                for(moris::moris_index iCMVerts = 0; iCMVerts < tNumNewVertices; iCMVerts++)
                {
                    moris_index tNewNodeLocInDecomp = aDecompositionData->tCMNewNodeLoc(iCM)(iCMVerts);
                    moris_index tNewNodeIndex = aDecompositionData->tNewNodeIndex(tNewNodeLocInDecomp);

                    aCutIntegrationMesh->mIntegrationVertexGroups(iCM)->mIgVertexGroup(iCMVerts + tStartIndex) = aCutIntegrationMesh->mIntegrationVertices(tNewNodeIndex);
                }
            }
        }

        void
        commit_new_ig_cells_to_cut_mesh(
            Integration_Mesh_Generation_Data*  aMeshGenerationData,
            Decomposition_Data *               aDecompositionData,
            Cut_Integration_Mesh *             aCutIntegrationMesh,
            moris::mtk::Mesh *                 aBackgroundMesh,
            Decomposition_Algorithm*           aDecompositionAlgorithm);


    private:
    void
    setup_subdivision_methods(Cell<enum Subdivision_Method> aMethods)
    {
        MORIS_ERROR(aMethods.size() == 2,  "Invalid methods provided to integration mesh generation, to be extended to support more cases");

        if(aMethods.size() == 1)
        {
            if(aMethods(0) == Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4 
            || aMethods(0) == Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8)
            {
                mRegularSubdivision = aMethods(0);
                mConformalSubdivision = Subdivision_Method::NO_METHOD;
            }
            else
            {
                mRegularSubdivision = Subdivision_Method::NO_METHOD;
                mConformalSubdivision = aMethods(0);
            }
        }

        else if(aMethods.size() == 2)
        {
            MORIS_ERROR(aMethods(0) == Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4 || aMethods(0) == Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8,"Regular subdivision method must come first");
            MORIS_ERROR(aMethods(1) == Subdivision_Method::C_TRI3 || aMethods(1) == Subdivision_Method::C_HIERARCHY_TET4,"Recursive conformal method must come second");
            mRegularSubdivision   = aMethods(0);
            mConformalSubdivision = aMethods(1);
        }
    }
    // Parallel assignment of new integration vertex ids
    public:
    void
    assign_node_requests_identifiers( 
            Decomposition_Data &  aDecompData,
            Cut_Integration_Mesh* aCutIntegrationMesh,
            moris::mtk::Mesh*     aBackgroundMesh,
            moris::moris_index    aMPITag)
    {
        moris_index tNodeIndex = aCutIntegrationMesh->get_first_available_index(EntityRank::NODE);

        for(moris::uint i = 0; i < aDecompData.tNewNodeIndex.size(); i++)
        {
            // set the new node index
            aDecompData.tNewNodeIndex(i) = tNodeIndex;
            tNodeIndex++;
        }

        barrier();
        // asserts
        MORIS_ERROR(aDecompData.tNewNodeId.size() == aDecompData.tNewNodeIndex.size(),
                "Dimension mismatch in assign_node_requests_identifiers");
        MORIS_ERROR(aDecompData.tNewNodeId.size() == aDecompData.tNewNodeParentRank.size(),
                "Dimension mismatch in assign_node_requests_identifiers");
        MORIS_ERROR(aDecompData.tNewNodeId.size() == aDecompData.tNewNodeParentIndex.size(),
                "Dimension mismatch in assign_node_requests_identifiers");

        // owned requests and shared requests sorted by owning proc
        Cell<uint> tOwnedRequest;
        Cell<Cell<uint>> tNotOwnedRequests;
        Cell<uint> tProcRanks;
        std::unordered_map<moris_id,moris_id> tProcRankToDataIndex;
        this->sort_new_node_requests_by_owned_and_not_owned(
                aDecompData,
                aBackgroundMesh,
                tOwnedRequest,
                tNotOwnedRequests,
                tProcRanks,
                tProcRankToDataIndex);

        // allocate ids for nodes I own
        moris::moris_id tNodeId  = aCutIntegrationMesh->allocate_entity_ids(aDecompData.tNewNodeId.size(), EntityRank::NODE);

        // Assign owned request identifiers
        this->assign_owned_request_id(aDecompData, tOwnedRequest, tNodeId);

        // prepare node information request data
        Cell<Matrix<IndexMat>> tOutwardRequests;
        this->setup_outward_requests(aDecompData, aBackgroundMesh, tNotOwnedRequests, tProcRanks, tProcRankToDataIndex, tOutwardRequests);

        // send requests to owning processor
        mXTKModel->send_outward_requests(aMPITag,tProcRanks,tOutwardRequests);

        // hold on to make sure everyone has sent all their information
        barrier();

        // receive the requests
        Cell<Matrix<IndexMat>> tReceivedRequests;
        Cell<uint> tProcsReceivedFrom;
        mXTKModel->inward_receive_requests(aMPITag, 3, tReceivedRequests, tProcsReceivedFrom);

        // Prepare request answers
        Cell<Matrix<IndexMat>> tRequestAnwers;
        this->prepare_request_answers(aDecompData,aBackgroundMesh,tReceivedRequests,tRequestAnwers);

        // send the answers back
        mXTKModel->return_request_answers(aMPITag+1, tRequestAnwers, tProcsReceivedFrom);

        barrier();

        // receive the answers
        Cell<Matrix<IndexMat>> tReceivedRequestsAnswers;
        mXTKModel->inward_receive_request_answers(aMPITag+1,1,tProcRanks,tReceivedRequestsAnswers);

        // handle received information
        this->handle_received_request_answers(aDecompData, aBackgroundMesh, tOutwardRequests,tReceivedRequestsAnswers,tNodeId);

        MORIS_ERROR(mXTKModel->verify_successful_node_assignment(aDecompData),
                "Unsuccesssful node assignment detected.");

        barrier();
    }

   void
   sort_new_node_requests_by_owned_and_not_owned(
            Decomposition_Data                    & tDecompData,
            moris::mtk::Mesh*                       aBackgroundMesh,
            Cell<uint>                            & aOwnedRequests,
            Cell<Cell<uint>>                      & aNotOwnedRequests,
            Cell<uint>                            & aProcRanks,
            std::unordered_map<moris_id,moris_id> & aProcRankToIndexInData)
    {
        // access the communication
        Matrix<IdMat> tCommTable = aBackgroundMesh->get_communication_table();

        // number of new nodes
        moris::uint tNumNewNodes = tDecompData.tNewNodeParentIndex.size();

        // Par rank
        moris::moris_index tParRank = par_rank();

        // resize proc ranks and setup map to comm table
        aProcRanks.resize(tCommTable.numel());
        for(moris::uint i = 0; i <tCommTable.numel(); i++)
        {
            aProcRankToIndexInData[tCommTable(i)] = i;
            aProcRanks(i) = (tCommTable(i));
            aNotOwnedRequests.push_back(Cell<uint>(0));
        }

        // iterate through each node request and figure out the owner
        for(moris::uint i = 0; i <tNumNewNodes; i++)
        {
            // Parent Rank
            enum EntityRank    tParentRank  = tDecompData.tNewNodeParentRank(i);
            moris::moris_index tParentIndex = tDecompData.tNewNodeParentIndex(i);

            // get the owner processor
            moris::moris_index tOwnerProc = aBackgroundMesh->get_entity_owner(tParentIndex,tParentRank);

            // If i own the request keep track of the index
            if(tOwnerProc == tParRank)
            {
                aOwnedRequests.push_back(i);
            }
            else
            {
                moris_index tIndex = aProcRankToIndexInData[tOwnerProc];

                aNotOwnedRequests(tIndex).push_back(i);
            }
        }
    }

    void
    assign_owned_request_id(
            Decomposition_Data & aDecompData,
            Cell<uint> const &   aOwnedRequest,
            moris::moris_id &    aNodeId)
    {
        for(moris::uint i = 0; i < aOwnedRequest.size(); i++)
        {
            moris_index tRequestIndex = aOwnedRequest(i);

            // set the new node id
            aDecompData.tNewNodeId(tRequestIndex) = aNodeId;
            aNodeId++;

            // increment number of new nodes with set ids (for assertion purposes)
            aDecompData.mNumNewNodesWithIds++;
        }
    }

    void
    setup_outward_requests(
            Decomposition_Data              const & aDecompData,
            moris::mtk::Mesh*                       aBackgroundMesh,
            Cell<Cell<uint>>                const & aNotOwnedRequests,
            Cell<uint>                      const & aProcRanks,
            std::unordered_map<moris_id,moris_id> & aProcRankToIndexInData,
            Cell<Matrix<IndexMat>>                & aOutwardRequests)
    {
        // size data
        aOutwardRequests.resize(aProcRanks.size());

        // iterate through the processors we need information from and package the matrix
        for(moris::uint i = 0; i < aProcRanks.size(); i++)
        {
            uint tProcRank = aProcRanks(i);

            MORIS_ASSERT(aProcRankToIndexInData.find(tProcRank) != aProcRankToIndexInData.end(),"Proc rank not in map");
            uint tIndexInData = aProcRankToIndexInData[tProcRank];

            uint tNumRequests = aNotOwnedRequests(tIndexInData).size();

            // size the sending matrix
            // column - request
            //   r0 - parent entity id
            //   r1 - parent entity rank
            //   r2 - Secondary id
            if(tNumRequests > 0)
            {
                aOutwardRequests(i) = moris::Matrix<IndexMat>(3,tNumRequests);
            }

            else
            {
                aOutwardRequests(i) = moris::Matrix<IndexMat>(3,1,MORIS_INDEX_MAX);
            }

            // populate matrix to send;
            for(moris::uint j = 0; j < tNumRequests; j++)
            {
                moris_index     tRequestIndex = aNotOwnedRequests(tIndexInData)(j);
                moris_index     tParentIndex  = aDecompData.tNewNodeParentIndex(tRequestIndex);
                moris_index     tSecondaryId  = aDecompData.tSecondaryIdentifiers(tRequestIndex);
                enum EntityRank tParentRank   = aDecompData.tNewNodeParentRank(tRequestIndex);

                aOutwardRequests(i)(0,j) = aBackgroundMesh->get_glb_entity_id_from_entity_loc_index(tParentIndex,tParentRank);
                aOutwardRequests(i)(1,j) = (moris_index)tParentRank;
                aOutwardRequests(i)(2,j) = tSecondaryId;
            }
        }
    }

    void
    prepare_request_answers(
            Decomposition_Data           & aDecompData,
            moris::mtk::Mesh*              aBackgroundMesh,
            Cell<Matrix<IndexMat>> const & aReceiveData,
            Cell<Matrix<IndexMat>>       & aRequestAnswers)
    {
        // allocate answer size
        aRequestAnswers.resize(aReceiveData.size());

        // iterate through received data
        for(moris::uint i = 0; i < aReceiveData.size(); i++)
        {
            uint tNumReceivedReqs = aReceiveData(i).n_cols();

            aRequestAnswers(i).resize(1,tNumReceivedReqs);

            aRequestAnswers(i)(0) = MORIS_INDEX_MAX;

            // avoid the dummy message
            if(aReceiveData(i)(0,0) != MORIS_INDEX_MAX)
            {
                // iterate through received requests
                for(moris::uint j = 0; j < tNumReceivedReqs; j++)
                {
                    moris_id        tParentId      = aReceiveData(i)(0,j);
                    enum EntityRank tParentRank    = (enum EntityRank) aReceiveData(i)(1,j);
                    moris_id        tSecondaryId   = aReceiveData(i)(2,j);
                    moris_index     tParentInd     = aBackgroundMesh->get_loc_entity_ind_from_entity_glb_id(tParentId,tParentRank);
                    bool            tRequestExists = false;
                    moris_index     tRequestIndex  = MORIS_INDEX_MAX;

                    if(aDecompData.mHasSecondaryIdentifier)
                    {
                        tRequestExists = aDecompData.request_exists(
                                tParentInd,
                                tSecondaryId,
                                (EntityRank)tParentRank,
                                tRequestIndex);
                    }
                    else
                    {
                        tRequestExists = aDecompData.request_exists(
                                tParentInd,
                                (EntityRank)tParentRank,
                                tRequestIndex);
                    }

                    if(tRequestExists)
                    {
                        moris_id tNodeId =aDecompData.tNewNodeId(tRequestIndex);

                        aRequestAnswers(i)(j) = tNodeId;

                        if(tNodeId == MORIS_ID_MAX)
                        {
                            std::cout<<"tParentId = "<<tParentId<<" | Rank "<<(uint)tParentRank<<std::endl;
                            //                    MORIS_ERROR(0,"Max node");
                        }
                    }
                    else
                    {
                        aRequestAnswers(i)(j) = MORIS_ID_MAX;
                    }
                }
            }
        }
    }
    void
    handle_received_request_answers(
            Decomposition_Data           & aDecompData,
            moris::mtk::Mesh*              aBackgroundMesh,
            Cell<Matrix<IndexMat>> const & aRequests,
            Cell<Matrix<IndexMat>> const & aRequestAnswers,
            moris::moris_id              & aNodeId)
    {
        Cell<moris_index> tUnhandledRequestIndices;

        // iterate through received data
        for(moris::uint i = 0; i < aRequests.size(); i++)
        {
            uint tNumReceivedReqs = aRequests(i).n_cols();

            // avoid the dummy message
            if(aRequests(i)(0,0) != MORIS_INDEX_MAX)
            {
                // iterate through received requests
                for(moris::uint j = 0; j < tNumReceivedReqs; j++)
                {
                    moris_id        tParentId      = aRequests(i)(0,j);
                    enum EntityRank tParentRank    = (enum EntityRank) aRequests(i)(1,j);
                    moris_id        tSecondaryId   = aRequests(i)(2,j);
                    moris_index     tParentInd     = aBackgroundMesh->get_loc_entity_ind_from_entity_glb_id(tParentId,tParentRank);
                    bool            tRequestExists = false;
                    moris_index     tRequestIndex  = MORIS_INDEX_MAX;

                    if(aDecompData.mHasSecondaryIdentifier)
                    {
                        tRequestExists = aDecompData.request_exists(tParentInd,tSecondaryId,(EntityRank)tParentRank,tRequestIndex);
                    }
                    else
                    {
                        tRequestExists = aDecompData.request_exists(tParentInd,(EntityRank)tParentRank,tRequestIndex);
                    }

                    if(tRequestExists && aRequestAnswers(i)(j))
                    {
                        moris_id tNodeId =aRequestAnswers(i)(j);

                        // meaning the owning processor expected this and gave an answer
                        if(tNodeId < MORIS_ID_MAX && aDecompData.tNewNodeId(tRequestIndex) == MORIS_INDEX_MAX)
                        {
                            // set the new node id
                            aDecompData.tNewNodeId(tRequestIndex) = tNodeId;

                            aDecompData.mNumNewNodesWithIds++;
                        }
                        // The owner did not expect and did not return an answer
                        else
                        {   
                            // keep track of unhandled
                            tUnhandledRequestIndices.push_back(tRequestIndex);
                            // moris_index tNodeIndex = mBackgroundMesh.get_first_available_index(EntityRank::NODE);

                            // aDecompData.tNewNodeOwner(tRequestIndex) = par_rank();

                            // aDecompData.tNewNodeId(tRequestIndex) = tNodeId;
                            // aDecompData.tNewNodeIndex(tRequestIndex) = tNodeIndex;
                            // tNodeIndex++;

                            // // set the new node id
                            // aDecompData.tNewNodeId(tRequestIndex) = tNodeId;

                            // aDecompData.mNumNewNodesWithIds++;

                            // mBackgroundMesh.update_first_available_index(tNodeIndex, EntityRank::NODE);

                        }
                    }
                    else
                    {
                        MORIS_ASSERT(0,"Request does not exist.");
                    }
                }
            }
        }
    }
};

Integration_Mesh_Generator::Integration_Mesh_Generator( xtk::Model*                    aXTKModelPtr,
                                                        Cell<enum Subdivision_Method>  aMethods,
                                                        moris::Matrix<moris::IndexMat> aActiveGeometries):
    mXTKModel(aXTKModelPtr),
    mGeometryEngine(mXTKModel->get_geom_engine()),
    mActiveGeometries(aActiveGeometries)
{
    this->setup_subdivision_methods(aMethods);

    if(mActiveGeometries.numel() == 0)
    {
        mAllActiveGeometries = true;
    }
}

Integration_Mesh_Generator::~Integration_Mesh_Generator()
{
    
}
}
#endif