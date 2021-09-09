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
        std::unordered_map<moris_index, moris_index> mIntersectedBackgroundCellIndexToChildMeshIndex;
    };

    class Geometric_Query_XTK : public moris::ge::Geometric_Query_Interface
    {
    private:
        // tell the geometry engine what you are interested
        enum moris::ge::Query_Type mQueryType;

        // associated with a given child mesh
        enum EntityRank mQueryEntityRank;
        Matrix<IndexMat> mQueryEntityToVertices;
        Matrix<IndexMat> mQueryEntityParamCoords;

        // parent cell is the query object when checking for intersection
        moris::mtk::Cell *mQueryParentCell;

        // for edge based questions
        moris_index mCurrentEdgeIndex;
        std::shared_ptr<Edge_Based_Connectivity>       mEdgeConnectivity;
        moris::Cell<moris::mtk::Cell*>*                mAssociatedBgCellForEdge;
        moris::Cell<std::shared_ptr<IG_Vertex_Group>>* mAssociatedVertexGroup;
        Cut_Integration_Mesh*                          mCutIntegrationMesh;

        // parent info
        moris::mtk::Cell *mParentCell;

        // row-based and indexed based on the indexing in mQueryEntityToVertices.
        moris::Cell<std::shared_ptr<moris::Matrix<moris::DDRMat>>> *mQueryEntityIndexedCoordinates;

        // geometric index
        moris_index mGeometricIndex;

    public:
        Geometric_Query_XTK() :mQueryParentCell(nullptr), mParentCell(nullptr), mGeometricIndex(MORIS_INDEX_MAX)
        {
        }
        ~Geometric_Query_XTK() {}

        void
        set_query_type(enum moris::ge::Query_Type aQueryType)
        {
            mQueryType = aQueryType;
        }

        void
        set_query_entity_rank(enum moris::EntityRank aEntityRank)
        {
            mQueryEntityRank = aEntityRank;
        }

        void
        set_edge_connectivity(std::shared_ptr<Edge_Based_Connectivity> aEdgeConnectivity)
        {
            mEdgeConnectivity = aEdgeConnectivity;
        }

        void
        set_edge_associated_background_cell(moris::Cell<moris::mtk::Cell*>* aAssociatedBgCellsForEdge)
        {
            mAssociatedBgCellForEdge = aAssociatedBgCellsForEdge;
        }

        void
        set_associated_vertex_group(moris::Cell<std::shared_ptr<IG_Vertex_Group>>* aAssociatedVertexGroup)
        {
            mAssociatedVertexGroup= aAssociatedVertexGroup;
        }

        void
        set_cut_integration_mesh(Cut_Integration_Mesh* aCutIntegrationMesh)
        {
            mCutIntegrationMesh = aCutIntegrationMesh;
        }

        void
        set_current_edge_index(moris_index aCurrentIndex)
        {
            mCurrentEdgeIndex = aCurrentIndex;
            mQueryEntityToVertices.resize(1,2);
            mQueryEntityToVertices(0) = mEdgeConnectivity->mEdgeVertices(mCurrentEdgeIndex)(0)->get_index();
            mQueryEntityToVertices(1) = mEdgeConnectivity->mEdgeVertices(mCurrentEdgeIndex)(1)->get_index();
        
        }

        void
        set_parametric_coordinate_size(moris_index aParametricSize)
        {
            // hard coded for edges right now
            mQueryEntityParamCoords.resize(2,aParametricSize);
        }

        void
        set_parent_cell(moris::mtk::Cell *aParentCell)
        {
            mParentCell = aParentCell;
        }

        void
        set_query_cell(moris::mtk::Cell *aQueryParentCell)
        {
            mQueryParentCell = aQueryParentCell;
            mQueryEntityToVertices = aQueryParentCell->get_vertex_inds();
        }

        void
        set_query_edge_connectivity(std::shared_ptr<Edge_Based_Connectivity> aEdgeBasedConnectivity)
        {

        }

        void
        set_coordinates_matrix(moris::Cell<std::shared_ptr<moris::Matrix<moris::DDRMat>>> *aCoordinates)
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

        moris::Cell<std::shared_ptr<moris::Matrix<moris::DDRMat>>> *get_query_indexed_coordinates() const
        {
            return mQueryEntityIndexedCoordinates;
        }

        Matrix<DDRMat>
        get_vertex_local_coord_wrt_parent_entity(moris_index aVertexIndex) const
        {
            std::shared_ptr<IG_Vertex_Group> tCurrentVertexGroup = (*mAssociatedVertexGroup)(mCurrentEdgeIndex);
            return *tCurrentVertexGroup->get_vertex_local_coords(aVertexIndex);
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

        Matrix<DDRMat>
        get_query_parent_coordinates() const
        {
            return mParentCell->get_vertex_coords();
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
        Model *mXTKModel;
        moris::ge::Geometry_Engine *mGeometryEngine;
        moris::Matrix<moris::IndexMat> mActiveGeometries;
        enum Subdivision_Method mRegularSubdivision;   // Required to have a inheritance structure
        enum Subdivision_Method mConformalSubdivision; // Recursive no inheritance

    public:
        Integration_Mesh_Generator(
            xtk::Model *aXTKModelPtr,
            Cell<enum Subdivision_Method> aMethods,
            moris::Matrix<moris::IndexMat> aActiveGeometries);
        ~Integration_Mesh_Generator();

        bool
        perform();

        moris::Matrix<moris::IndexMat> const *
        get_active_geometries();

        moris::ge::Geometry_Engine *
        get_geom_engine();

        bool
        determine_intersected_background_cells(
            Integration_Mesh_Generation_Data &aMeshGenerationData,
            Cut_Integration_Mesh *aCutIntegrationMesh,
            moris::mtk::Mesh *aBackgroundMesh);


        bool
        allocate_child_meshes(Integration_Mesh_Generation_Data &aMeshGenerationData,
                              Cut_Integration_Mesh *aCutIntegrationMesh,
                              moris::mtk::Mesh *aBackgroundMesh);

        void
        collect_vertex_groups_for_background_cells(
           Integration_Mesh_Generation_Data*              aMeshGenerationData,
           Cut_Integration_Mesh*                          aCutIntegrationMesh,
           moris::Cell<moris::mtk::Cell*>*                aBackgroundCells,
           moris::Cell<std::shared_ptr<IG_Vertex_Group>>* aVertexGroups);

        moris_index
        edge_exists(
            moris::Cell<moris::mtk::Vertex *> &aVerticesOnEdge,
            std::unordered_map<moris_index, moris_index> &aLocaLVertexMap,
            moris::Cell<moris::Cell<uint>> &aVertexToEdge,
            moris::Cell<moris::Cell<moris::mtk::Vertex *>> &aFullEdgeVertices);

        void
        create_edges_from_element_to_node(
            moris::Cell<moris::mtk::Cell *> aCells,
            std::shared_ptr<Edge_Based_Connectivity> aEdgeConnectivity);

        void
        deduce_edge_ancestry(
            Cut_Integration_Mesh *                   aCutIntegrationMesh, 
            moris::mtk::Mesh*                        aBackgroundMesh,
            std::shared_ptr<Edge_Based_Connectivity> aIgCellGroupEdgeConnectivity, 
            moris::Cell<moris::mtk::Cell*> const &   aParentCellForDeduction, 
            std::shared_ptr<Edge_Based_Ancestry>     aIgEdgeAncestry);

        void
        commit_new_ig_vertices_to_cut_mesh(
            Integration_Mesh_Generation_Data *aMeshGenerationData,
            Decomposition_Data               *aDecompositionData,
            Cut_Integration_Mesh             *aCutIntegrationMesh,
            moris::mtk::Mesh          *aBackgroundMesh,
            Decomposition_Algorithm *aDecompositionAlgorithm);

        void
        link_new_vertices_to_geometry_engine(
            Decomposition_Data *aDecompositionData,
            Decomposition_Algorithm *aDecompAlg);

        void
        commit_new_ig_cells_to_cut_mesh(
            Integration_Mesh_Generation_Data *aMeshGenerationData,
            Decomposition_Data *aDecompositionData,
            Cut_Integration_Mesh *aCutIntegrationMesh,
            moris::mtk::Mesh *aBackgroundMesh,
            Decomposition_Algorithm *aDecompositionAlgorithm);

        void
        extract_cells_from_cell_groups(
            moris::Cell<std::shared_ptr<IG_Cell_Group>> const & aCellGroups,
            moris::Cell<moris::mtk::Cell *> &aCellsInGroups);

    private:
        void
        setup_subdivision_methods(Cell<enum Subdivision_Method> aMethods);
        // Parallel assignment of new integration vertex ids
    public:
        void
        assign_node_requests_identifiers(
            Decomposition_Data &aDecompData,
            Cut_Integration_Mesh *aCutIntegrationMesh,
            moris::mtk::Mesh *aBackgroundMesh,
            moris::moris_index aMPITag);

        void
        sort_new_node_requests_by_owned_and_not_owned(
            Decomposition_Data &tDecompData,
            moris::mtk::Mesh *aBackgroundMesh,
            Cell<uint> &aOwnedRequests,
            Cell<Cell<uint>> &aNotOwnedRequests,
            Cell<uint> &aProcRanks,
            std::unordered_map<moris_id, moris_id> &aProcRankToIndexInData);

        void
        assign_owned_request_id(
            Decomposition_Data &aDecompData,
            Cell<uint> const &aOwnedRequest,
            moris::moris_id &aNodeId);

        void
        setup_outward_requests(
            Decomposition_Data const &aDecompData,
            moris::mtk::Mesh *aBackgroundMesh,
            Cell<Cell<uint>> const &aNotOwnedRequests,
            Cell<uint> const &aProcRanks,
            std::unordered_map<moris_id, moris_id> &aProcRankToIndexInData,
            Cell<Matrix<IndexMat>> &aOutwardRequests);

        void
        prepare_request_answers(
            Decomposition_Data &aDecompData,
            moris::mtk::Mesh *aBackgroundMesh,
            Cell<Matrix<IndexMat>> const &aReceiveData,
            Cell<Matrix<IndexMat>> &aRequestAnswers);
        void
        handle_received_request_answers(
            Decomposition_Data &aDecompData,
            moris::mtk::Mesh *aBackgroundMesh,
            Cell<Matrix<IndexMat>> const &aRequests,
            Cell<Matrix<IndexMat>> const &aRequestAnswers,
            moris::moris_id &aNodeId);
    };

}

#endif