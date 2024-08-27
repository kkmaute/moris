/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Node_Hierarchy_Interface.hpp
 *
 */

#ifndef MORIS_CL_XTK_NODE_HIERARCHY_INTERFACE_HPP_
#define MORIS_CL_XTK_NODE_HIERARCHY_INTERFACE_HPP_

#include "cl_XTK_Decomposition_Algorithm.hpp"

namespace moris
{
    namespace ge
    {
        class Geometry_Engine;
    }
    namespace mtk
    {
        class Mesh;
    }
}    // namespace moris

namespace moris::xtk
{
    class Integration_Mesh_Generator;
    class Cut_Integration_Mesh;

    struct Integration_Mesh_Generation_Data;
    struct Decomposition_Data;
    struct IG_Cell_Group;

    struct Node_Hierarchy_Template
    {
        mtk::CellTopology  mCellTopology;
        moris_index        mNumCells;
        Matrix< IndexMat > mCellToNodeOrdinal;
    };

    class Node_Hierarchy_Template_Library
    {
      public:
        Node_Hierarchy_Template_Library()
        {
        }

        void
        load_template(
                moris_index              aSpatialDim,
                moris_index              aTemplateId,
                Node_Hierarchy_Template *aNodeHierTemplate );

        void
        load_2d_template(
                moris_index              aTemplateId,
                Node_Hierarchy_Template *aNodeHierTemplate );

        void
        load_3d_template(
                moris_index              aTemplateId,
                Node_Hierarchy_Template *aNodeHierTemplate );
    };

    class Node_Hierarchy_Interface : public Decomposition_Algorithm
    {
      private:
        moris::gen::Geometry_Engine      *mGeometryEngine;
        Integration_Mesh_Generation_Data *mMeshGenerationData;
        Decomposition_Data               *mDecompositionData;
        Cut_Integration_Mesh             *mCutIntegrationMesh;
        moris::mtk::Mesh                 *mBackgroundMesh;
        Integration_Mesh_Generator       *mGenerator;

      public:
        Node_Hierarchy_Interface( Parameter_List &aParameterList ) {}
        ~Node_Hierarchy_Interface() {}

        bool has_geometric_independent_vertices() const;

        /**
         * @brief perform the node hierarchy decomposition. We overwrite the existing default implementation because
         * this method iterates over geometries and constructs vertices with each one
         *
         * @param aMeshGenerationData Mesh generation data
         * @param aDecompositionData  Decomposition data
         * @param aCutIntegrationMesh Cut integration mesh
         * @param aBackgroundMesh     Background mesh
         * @param aMeshGenerator      Mesh generator
         */

        void perform(
                Integration_Mesh_Generation_Data *aMeshGenerationData,
                Decomposition_Data               *aDecompositionData,
                Cut_Integration_Mesh             *aCutIntegrationMesh,
                moris::mtk::Mesh                 *aBackgroundMesh,
                Integration_Mesh_Generator       *aMeshGenerator );
        void
        perform_impl_vertex_requests(
                Integration_Mesh_Generation_Data *aMeshGenerationData,
                Decomposition_Data               *aDecompositionData,
                Cut_Integration_Mesh             *aCutIntegrationMesh,
                moris::mtk::Mesh                 *aBackgroundMesh,
                Integration_Mesh_Generator       *aMeshGenerator );

        void
        perform_impl_generate_mesh(
                Integration_Mesh_Generation_Data *aMeshGenerationData,
                Decomposition_Data               *aDecompositionData,
                Cut_Integration_Mesh             *aCutIntegrationMesh,
                moris::mtk::Mesh                 *aBackgroundMesh,
                Integration_Mesh_Generator       *aMeshGenerator );

        enum Decomposition_Algorithm_Type
        get_algorithm_type() const;

        moris_index get_signature() const;

      private:
        /**
         *  @brief Pick out the integration cell groups from the cut integration mesh
         *
         * @param aIgCellGroups
         */
        void
        select_ig_cell_groups(
                Vector< std::shared_ptr< IG_Cell_Group > > &aIgCellGroups );

        bool
        determine_intersected_edges_and_make_requests(
                std::shared_ptr< Edge_Based_Connectivity >    aEdgeConnectivity,
                std::shared_ptr< Edge_Based_Ancestry >        aIgEdgeAncestry,
                Vector< moris::mtk::Cell * >                 *aBackgroundCellForEdge,
                Vector< std::shared_ptr< IG_Vertex_Group > > *aVertexGroups,
                Vector< moris_index >                        &aIntersectedEdges,
                Vector< moris::real >                        &aEdgeLocalCoordinate );

        /**
         * @brief Creates unique id for an edge based on the IDs of its two end vertices
         * using the Cantor pairing function
         *
         * @param aEdgeVertices list of mtk::vertices on edge
         * @return moris_index unique id for edge
         */

        moris_index
        hash_edge( Vector< moris::mtk::Vertex * > const &aEdgeVertices );

        bool
        associate_new_vertices_with_cell_groups(
                std::shared_ptr< Edge_Based_Connectivity >    aEdgeConnectivity,
                std::shared_ptr< Edge_Based_Ancestry >        aIgEdgeAncestry,
                Vector< moris::mtk::Cell * >                 *aBackgroundCellForEdge,
                Vector< std::shared_ptr< IG_Vertex_Group > > *aVertexGroups,
                Vector< moris_index >                        *aIntersectedEdges,
                Vector< moris::real >                        *aEdgeLocalCoordinate );

        void
        create_node_hierarchy_integration_cells(
                std::shared_ptr< Edge_Based_Connectivity > aEdgeConnectivity,
                std::shared_ptr< Edge_Based_Ancestry >     aIgEdgeAncestry,
                Vector< moris_index >                     *aIntersectedEdges );

        void
        determine_intersected_cell_information(
                std::shared_ptr< Edge_Based_Connectivity >                   aEdgeConnectivity,
                Vector< moris_index >                                       *aIntersectedEdges,
                Vector< std::shared_ptr< Vector< moris_index > > >          *aCellIndexIntersectedEdgeOrdinals,
                Vector< std::shared_ptr< Vector< moris::mtk::Vertex * > > > *aCellIndexIntersectedEdgeVertex );

        moris_index
        select_node_hier_2d_template(
                Vector< std::shared_ptr< Vector< moris_index > > >          *aCellIndexIntersectedEdgeOrdinals,
                Vector< std::shared_ptr< Vector< moris::mtk::Vertex * > > > *aCellIndexIntersectedEdgeVertex,
                Vector< std::shared_ptr< Vector< moris::mtk::Vertex * > > > *aNodesForTemplates,
                Vector< std::shared_ptr< Node_Hierarchy_Template > >        *aNHTemplate );

        moris_index
        select_node_hier_3d_template(
                Vector< std::shared_ptr< Vector< moris_index > > >          *aCellIndexIntersectedEdgeOrdinals,
                Vector< std::shared_ptr< Vector< moris::mtk::Vertex * > > > *aCellIndexIntersectedEdgeVertex,
                Vector< std::shared_ptr< Vector< moris::mtk::Vertex * > > > *aNodesForTemplates,
                Vector< std::shared_ptr< Node_Hierarchy_Template > >        *aNHTemplate );

        void
        sort_nodes_2d(
                moris::mtk::Cell const                           *aIgCell,
                Matrix< IndexMat >                               *aEdgeToVertexOrdinalMap,
                std::shared_ptr< Vector< moris_index > >          aCellIndexIntersectedEdgeOrdinals,
                std::shared_ptr< Vector< moris::mtk::Vertex * > > aCellIndexIntersectedEdgeVertex,
                moris_index                                      &aPermutation,
                std::shared_ptr< Vector< moris::mtk::Vertex * > > aSortedNodeInds );

        void
        sort_nodes_3d(
                moris::mtk::Cell const                           *aIgCell,
                Matrix< IndexMat >                               *aEdgeToVertexOrdinalMap,
                std::shared_ptr< Vector< moris_index > >          aCellIndexIntersectedEdgeOrdinals,
                std::shared_ptr< Vector< moris::mtk::Vertex * > > aCellIndexIntersectedEdgeVertex,
                moris_index                                      &aPermutation,
                std::shared_ptr< Vector< moris::mtk::Vertex * > > aSortedNodeInds );
    };

}    // namespace moris::xtk
#endif
