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
#include "cl_Tracer.hpp"

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
    moris::Cell<std::shared_ptr<moris::Matrix<moris::DDRMat>>>* mQueryEntityIndexedCoordinates;    

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
    set_coordinates_matrix(moris::Cell<std::shared_ptr<moris::Matrix<moris::DDRMat>>> * aCoordinates)
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

    moris::Cell<std::shared_ptr<moris::Matrix<moris::DDRMat>>> * get_query_indexed_coordinates() const
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
                                                moris::mtk::Mesh*                  aBackgroundMesh );

        bool
        allocate_child_meshes(Integration_Mesh_Generation_Data & aMeshGenerationData, 
                              Cut_Integration_Mesh* aCutIntegrationMesh,
                              moris::mtk::Mesh*     aBackgroundMesh  );

        bool
        determine_intersected_ig_edges_in_cell_group( 
            Cell<moris
            Integration_Mesh_Generation_Data & aMeshGenerationData,
            Cut_Integration_Mesh*              aCutIntegrationMesh,
            moris::mtk::Mesh*                  aBackgroundMesh )
        {
            // Generate edge to node connectivity for the entire cut integration mesh
            // could add direct substitution of connectivity from reg sub to speed up
            this->setup_child_mesh_edge_based_connectivity_of_ig_mesh(aMeshGenerationData, aCutIntegrationMesh);

            // Initialize geometric query
            Geometric_Query_XTK tGeometricQuery;
            tGeometricQuery.set_query_type(moris::ge::Query_Type::INTERSECTION_LOCATION);
            tGeometricQuery.set_coordinates_matrix(&aCutIntegrationMesh->mVertexCoordinates);

            

            return true;
        }

        void
        setup_child_mesh_edge_based_connectivity_of_ig_mesh(
            Integration_Mesh_Generation_Data & aMeshGenerationData,
            Cut_Integration_Mesh*              aCutIntegrationMesh)
        {
            Tracer tTracer( "XTK", "Integration", "setup_child_mesh_edge_based_connectivity_of_ig_mesh" );
            // number of child meshes
            moris::uint tNumChildMeshes = aCutIntegrationMesh->mIntegrationCellGroups.size();

            // allocate the edge based connectivity data in the integration mesh
            aCutIntegrationMesh->mIntegrationCellGroupsEdgeBasedConn.resize(tNumChildMeshes);

            for(moris::uint iCM = 0; iCM < tNumChildMeshes; iCM++ )
            {
                aCutIntegrationMesh->mIntegrationCellGroupsEdgeBasedConn(iCM) = std::make_shared<Edge_Based_Connectivity>();

                this->create_edges_from_element_to_node(
                    aCutIntegrationMesh->mIntegrationCellGroups(iCM)->mIgCellGroup,
                    aCutIntegrationMesh->mIntegrationCellGroupsEdgeBasedConn(iCM));
            }
            
            
        }

    moris_index 
    edge_exists(
        moris::Cell<moris::mtk::Vertex*> & aVerticesOnEdge,
        std::unordered_map<moris_index,moris_index> & aLocaLVertexMap,
        moris::Cell<moris::Cell<uint>> & aVertexToEdge,
        moris::Cell<moris::Cell<moris::mtk::Vertex*>> & aFullEdgeVertices)
     {
        moris_index tEdgeIndex = MORIS_INDEX_MAX;
        
        // get them in order based on id
        std::sort(aVerticesOnEdge.data().begin(), aVerticesOnEdge.data().end(), moris::comparePtrToVertexIdBased);

        // local vertex index
        auto tIter = aLocaLVertexMap.find(aVerticesOnEdge(0)->get_index());
        MORIS_ERROR(tIter != aLocaLVertexMap.end(),"Invalid vertex detected.");
        moris_index tLocalVertexIndex = tIter->second;


        // iterate through edges attached to the first vertex (depend on ascending order)
        for(moris::uint iEdge = 0; iEdge < aVertexToEdge(tLocalVertexIndex).size(); iEdge++)
        {
            moris_index tEdgeIndex = aVertexToEdge(tLocalVertexIndex)(iEdge);

            MORIS_ASSERT(aFullEdgeVertices(tEdgeIndex)(0)->get_index() == aVerticesOnEdge(0)->get_index(),"Numbering issues, edges should be in ascending order based on vertex id");

            // check the second vertex on the edge
            if(aFullEdgeVertices(tEdgeIndex)(1)->get_index() == aVerticesOnEdge(1)->get_index())
            {
                return tEdgeIndex;
            }
        }
        return tEdgeIndex;
     }

    void
    create_edges_from_element_to_node( moris::Cell<moris::mtk::Cell*> aCells,
                                       std::shared_ptr<Edge_Based_Connectivity> aEdgeConnectivity);


    void
    commit_new_ig_vertices_to_cut_mesh(
        Integration_Mesh_Generation_Data*  aMeshGenerationData,
        Decomposition_Data *               aDecompositionData,
        Cut_Integration_Mesh *             aCutIntegrationMesh,
        moris::mtk::Mesh *                 aBackgroundMesh,
        Decomposition_Algorithm*           aDecompositionAlgorithm);

    void
    link_new_vertices_to_geometry_engine(
        Decomposition_Data *               aDecompositionData,
        Decomposition_Algorithm*           aDecompAlg);

    void
    commit_new_ig_cells_to_cut_mesh(
        Integration_Mesh_Generation_Data*  aMeshGenerationData,
        Decomposition_Data *               aDecompositionData,
        Cut_Integration_Mesh *             aCutIntegrationMesh,
        moris::mtk::Mesh *                 aBackgroundMesh,
        Decomposition_Algorithm*           aDecompositionAlgorithm);


    private:
    void
    setup_subdivision_methods(Cell<enum Subdivision_Method> aMethods);
    // Parallel assignment of new integration vertex ids
    public:
    void
    assign_node_requests_identifiers( 
            Decomposition_Data &  aDecompData,
            Cut_Integration_Mesh* aCutIntegrationMesh,
            moris::mtk::Mesh*     aBackgroundMesh,
            moris::moris_index    aMPITag);

   void
   sort_new_node_requests_by_owned_and_not_owned(
            Decomposition_Data                    & tDecompData,
            moris::mtk::Mesh*                       aBackgroundMesh,
            Cell<uint>                            & aOwnedRequests,
            Cell<Cell<uint>>                      & aNotOwnedRequests,
            Cell<uint>                            & aProcRanks,
            std::unordered_map<moris_id,moris_id> & aProcRankToIndexInData);

    void
    assign_owned_request_id(
            Decomposition_Data & aDecompData,
            Cell<uint> const &   aOwnedRequest,
            moris::moris_id &    aNodeId);

    void
    setup_outward_requests(
            Decomposition_Data              const & aDecompData,
            moris::mtk::Mesh*                       aBackgroundMesh,
            Cell<Cell<uint>>                const & aNotOwnedRequests,
            Cell<uint>                      const & aProcRanks,
            std::unordered_map<moris_id,moris_id> & aProcRankToIndexInData,
            Cell<Matrix<IndexMat>>                & aOutwardRequests);

    void
    prepare_request_answers(
            Decomposition_Data           & aDecompData,
            moris::mtk::Mesh*              aBackgroundMesh,
            Cell<Matrix<IndexMat>> const & aReceiveData,
            Cell<Matrix<IndexMat>>       & aRequestAnswers);
    void
    handle_received_request_answers(
            Decomposition_Data           & aDecompData,
            moris::mtk::Mesh*              aBackgroundMesh,
            Cell<Matrix<IndexMat>> const & aRequests,
            Cell<Matrix<IndexMat>> const & aRequestAnswers,
            moris::moris_id              & aNodeId);
    
};


}


#endif