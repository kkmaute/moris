/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Elevate_Order_Interface.hpp
 *
 */

#ifndef SRC_cl_XTK_Elevate_Order_Interface
#define SRC_cl_XTK_Elevate_Order_Interface

#include "cl_XTK_Decomposition_Algorithm.hpp"

namespace moris::mtk
{
    class Mesh;
}

namespace moris::xtk
{

    class Integration_Mesh_Generator;
    class Cut_Integration_Mesh;

    struct Integration_Mesh_Generation_Data;
    struct Decomposition_Data;
    struct IG_Cell_Group;

    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------

    // pure virtual base class for 2D and 3D Regular Subdivision Data
    class Elevate_Order_Template
    {
      public:
        Elevate_Order_Template() {}

        virtual ~Elevate_Order_Template() {}

        /**
         * @brief Get number of nodes to be added per element
         *
         * @return moris_index number of nodes to be added per element
         */
        virtual uint
        get_num_new_nodes() const = 0;

        /**
         * @brief Get total number of IG vertices on new elements for Template
         *
         * @return moris_index total number of IG vertices on new elements
         */
        virtual uint
        get_total_ig_verts() const = 0;

        /**
         * @brief gives back whether new vertices will be created on the requested entity rank
         *
         * @param aEntityRank Entity Rank for which new vertices should be created
         * @return whether new vertices will be created or not
         */
        virtual bool
        has_new_vertices_on_entity( mtk::EntityRank aEntityRank ) const = 0;

        /**
         * @brief gives back the number of new vertices that will be created on the requested entity rank
         *
         * @param aEntityRank Entity Rank for which new vertices should be created
         * @return number of new vertices on that entity rank
         */
        virtual uint
        num_new_vertices_per_entity( mtk::EntityRank aEntityRank ) const = 0;

        /**
         * @brief Get parametric coords for new vertices wrt entity they are created on
         *
         * @param aEntityRank Entity Rank for which new vertices should be created
         * @return Vector< Matrix< DDRMat > > list of parametric coords for vertices to be created per entity
         */
        virtual Vector< Matrix< DDRMat > >
        get_new_vertex_parametric_coords_wrt_entity( mtk::EntityRank aEntityRank ) const = 0;

        /**
         * @brief Get the new ig cell's topology
         *
         * @return mtk::CellTopology new ig cell topology
         */
        virtual mtk::CellTopology
        get_ig_cell_topology() const = 0;

        /**
         * @brief Get number of spatial dimensions current template assumes
         *
         * @return uint number of spatial dimensions
         */
        virtual uint
        get_num_spatial_dims() const = 0;

        /**
         * @brief Get local vertex index wrt to the cell based on where the vertex sits
         *
         * @param aEntityRank entity rank where the vertex sits
         * @param aSignedLocalEntityIndex 1-based index of the entity the vertex sits on with a sign to indicate directionality
         * @param aVertexIndexOnEntity index position of the vertex on that entity
         * @return moris_index element local vertex index
         */
        virtual moris_index
        get_local_vertex_index( mtk::EntityRank aEntityRank, moris_index aSignedLocalEntityIndex, moris_index aVertexIndexOnEntity ) const = 0;

        /**
         * @brief Get the element local edge index and direction (sign) based on vertex indices attached to it
         * @note The return value is 1-BASED!
         *
         * @param aFirstLocalVertexIndex local index of the vertex where edge starts
         * @param aSecondLocalVertexIndex local index of the vertex where edge ends
         * @return moris_index 1-Based index of the edge. The sign indicates the direction.
         */
        virtual moris_index
        get_local_edge_index_based_on_vertex_indices( moris_index aFirstLocalVertexIndex, moris_index aSecondLocalVertexIndex ) const = 0;
    };

    // -------------------------------------------------------------------------

    class TRI3_to_TRI6 : public Elevate_Order_Template
    {
      public:
        TRI3_to_TRI6() {}

        uint
        get_num_new_nodes() const override
        {
            return 3;
        }

        uint
        get_total_ig_verts() const override
        {
            return 6;
        }

        bool
        has_new_vertices_on_entity( mtk::EntityRank aEntityRank ) const override
        {
            switch ( aEntityRank )
            {
                case mtk::EntityRank::NODE:
                {
                    MORIS_ERROR( false, "TRI3_to_TRI6::has_new_nodes_on_entity() - No new vertices on vertices themselves. This shouldn't be requested." );
                    return false;
                    break;
                }
                case mtk::EntityRank::EDGE:
                {
                    return true;
                    break;
                }
                case mtk::EntityRank::FACE:
                {
                    return false;
                    break;
                }
                case mtk::EntityRank::ELEMENT:
                {
                    return false;
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "TRI3_to_TRI6::has_new_nodes_on_entity() - Unknown Entity Rank." );
                    return false;
                    break;
                }
            }
        }

        uint
        num_new_vertices_per_entity( mtk::EntityRank aEntityRank ) const override
        {
            switch ( aEntityRank )
            {
                case mtk::EntityRank::NODE:
                {
                    MORIS_ERROR( false, "TRI3_to_TRI6::num_new_vertices_per_entity() - No new vertices on vertices themselves. This shouldn't be requested." );
                    return 0;
                    break;
                }
                case mtk::EntityRank::EDGE:
                {
                    return 1;
                    break;
                }
                case mtk::EntityRank::FACE:
                {
                    return 0;
                    break;
                }
                case mtk::EntityRank::ELEMENT:
                {
                    return 0;
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "TRI3_to_TRI6::num_new_vertices_per_entity() - Unknown Entity Rank." );
                    return 0;
                    break;
                }
            }
        }

        Vector< Matrix< DDRMat > >
        get_new_vertex_parametric_coords_wrt_entity( mtk::EntityRank aEntityRank ) const override
        {
            switch ( aEntityRank )
            {
                case mtk::EntityRank::NODE:
                {
                    MORIS_ERROR( false, "TRI3_to_TRI6::get_new_vertex_parametric_coords_wrt_entity() - No new vertices on vertices themselves. This shouldn't be requested." );
                    return Vector< Matrix< DDRMat > >( 0 );
                    break;
                }
                case mtk::EntityRank::EDGE:
                {
                    return { { { 0.0 } } };
                    break;
                }
                case mtk::EntityRank::FACE:
                {
                    MORIS_ERROR( false, "TRI3_to_TRI6::get_new_vertex_parametric_coords_wrt_entity() - No new vertices on faces. This shouldn't be requested." );
                    return Vector< Matrix< DDRMat > >( 0 );
                    break;
                }
                case mtk::EntityRank::ELEMENT:
                {
                    MORIS_ERROR( false, "TRI3_to_TRI6::get_new_vertex_parametric_coords_wrt_entity() - No new vertices inside cell. This shouldn't be requested." );
                    return Vector< Matrix< DDRMat > >( 0 );
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "TRI3_to_TRI6::get_new_vertex_parametric_coords_wrt_entity() - Unknown Entity Rank." );
                    return Vector< Matrix< DDRMat > >( 0 );
                    break;
                }
            }
        }

        mtk::CellTopology
        get_ig_cell_topology() const override
        {
            return mtk::CellTopology::TRI6;
        }

        uint
        get_num_spatial_dims() const override
        {
            return 2;
        }

        moris_index
        get_local_vertex_index( mtk::EntityRank aEntityRank, moris_index aSignedLocalEntityIndex, moris_index aVertexIndexOnEntity ) const override
        {
            switch ( aEntityRank )
            {
                case mtk::EntityRank::NODE:
                {
                    MORIS_ERROR( false, "TRI3_to_TRI6::get_cell_local_vertex_index() - No new vertices on vertices themselves. This shouldn't be requested." );
                    return -1;
                    break;
                }
                case mtk::EntityRank::EDGE:
                {
                    // check if requested vertex location makes sense
                    MORIS_ASSERT( aSignedLocalEntityIndex != 0,
                            "TRI3_to_TRI6::get_cell_local_vertex_index() - aSignedLocalEntityIndex is 1-based and must not be 0." );

                    MORIS_ASSERT( std::abs( aSignedLocalEntityIndex ) <= 3 && aVertexIndexOnEntity < 1,
                            "TRI3_to_TRI6::get_cell_local_vertex_index() - Index out of bounds. Only 3 edges with 1 new node per edge." );

                    // convert Signed 1-based index to unsigned 0-based index
                    uint tLocalEntityIndex = std::abs( aSignedLocalEntityIndex ) - 1;

                    // create look-up table for local vertex position
                    Matrix< IndexMat > tLookUp = {
                        { 3 },
                        { 4 },
                        { 5 }
                    };

                    return tLookUp( tLocalEntityIndex, aVertexIndexOnEntity );
                    break;
                }
                case mtk::EntityRank::FACE:
                {
                    MORIS_ERROR( false, "TRI3_to_TRI6::get_cell_local_vertex_index() - No new vertices on faces. This shouldn't be requested." );
                    return -1;
                    break;
                }
                case mtk::EntityRank::ELEMENT:
                {
                    MORIS_ERROR( false, "TRI3_to_TRI6::get_cell_local_vertex_index() - No new vertices inside cell. This shouldn't be requested." );
                    return -1;
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "TRI3_to_TRI6::get_cell_local_vertex_index() - Unknown Entity Rank." );
                    return -1;
                    break;
                }
            }
        }

        moris_index
        get_local_edge_index_based_on_vertex_indices( moris_index aFirstLocalVertexIndex, moris_index aSecondLocalVertexIndex ) const override
        {
            // check if input makes sense
            MORIS_ASSERT( aFirstLocalVertexIndex != aSecondLocalVertexIndex,
                    "TRI3_to_TRI6::get_local_edge_index_based_on_vertex_indices() - Vertices must have different indices" );

            // create LookUp
            Matrix< IndexMat > tLookUp = {
                { 0, 1, -3 },
                { -1, 0, 2 },
                { 3, -2, 0 }
            };

            // return edge index (1-Based) and direction as sign
            return tLookUp( aFirstLocalVertexIndex, aSecondLocalVertexIndex );
        }
    };

    // -------------------------------------------------------------------------

    class TET4_to_TET10 : public Elevate_Order_Template
    {
      public:
        TET4_to_TET10() {}

        uint
        get_num_new_nodes() const override
        {
            return 6;
        }

        uint
        get_total_ig_verts() const override
        {
            return 10;
        }

        bool
        has_new_vertices_on_entity( mtk::EntityRank aEntityRank ) const override
        {
            switch ( aEntityRank )
            {
                case mtk::EntityRank::NODE:
                {
                    MORIS_ERROR( false, "TET4_to_TET10::has_new_nodes_on_entity() - No new vertices on vertices themselves. This shouldn't be requested." );
                    return false;
                    break;
                }
                case mtk::EntityRank::EDGE:
                {
                    return true;
                    break;
                }
                case mtk::EntityRank::FACE:
                {
                    return false;
                    break;
                }
                case mtk::EntityRank::ELEMENT:
                {
                    return false;
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "TET4_to_TET10::has_new_nodes_on_entity() - Unknown Entity Rank." );
                    return false;
                    break;
                }
            }
        }

        uint
        num_new_vertices_per_entity( mtk::EntityRank aEntityRank ) const override
        {
            switch ( aEntityRank )
            {
                case mtk::EntityRank::NODE:
                {
                    MORIS_ERROR( false, "TET4_to_TET10::num_new_vertices_per_entity() - No new vertices on vertices themselves. This shouldn't be requested." );
                    return 0;
                    break;
                }
                case mtk::EntityRank::EDGE:
                {
                    return 1;
                    break;
                }
                case mtk::EntityRank::FACE:
                {
                    return 0;
                    break;
                }
                case mtk::EntityRank::ELEMENT:
                {
                    return 1;
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "TET4_to_TET10::num_new_vertices_per_entity() - Unknown Entity Rank." );
                    return 0;
                    break;
                }
            }
        }

        Vector< Matrix< DDRMat > >
        get_new_vertex_parametric_coords_wrt_entity( mtk::EntityRank aEntityRank ) const override
        {
            switch ( aEntityRank )
            {
                case mtk::EntityRank::NODE:
                {
                    MORIS_ERROR( false, "TET4_to_TET10::get_new_vertex_parametric_coords_wrt_entity() - No new vertices on vertices themselves. This shouldn't be requested." );
                    return Vector< Matrix< DDRMat > >( 0 );
                    break;
                }
                case mtk::EntityRank::EDGE:
                {
                    return { { { 0.0 } } };
                    break;
                }
                case mtk::EntityRank::FACE:
                {
                    MORIS_ERROR( false, "TET4_to_TET10::get_new_vertex_parametric_coords_wrt_entity() - No new vertices on faces. This shouldn't be requested." );
                    return Vector< Matrix< DDRMat > >( 0 );
                    break;
                }
                case mtk::EntityRank::ELEMENT:
                {
                    MORIS_ERROR( false, "TET4_to_TET10::get_new_vertex_parametric_coords_wrt_entity() - No new vertices inside cell. This shouldn't be requested." );
                    return Vector< Matrix< DDRMat > >( 0 );
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "TET4_to_TET10::get_new_vertex_parametric_coords_wrt_entity() - Unknown Entity Rank." );
                    return Vector< Matrix< DDRMat > >( 0 );
                    break;
                }
            }
        }

        mtk::CellTopology
        get_ig_cell_topology() const override
        {
            return mtk::CellTopology::TET10;
        }

        uint
        get_num_spatial_dims() const override
        {
            return 3;
        }

        moris_index
        get_local_vertex_index( mtk::EntityRank aEntityRank, moris_index aSignedLocalEntityIndex, moris_index aVertexIndexOnEntity ) const override
        {
            switch ( aEntityRank )
            {
                case mtk::EntityRank::NODE:
                {
                    MORIS_ERROR( false, "TET4_to_TET10::get_cell_local_vertex_index() - No new vertices on vertices themselves. This shouldn't be requested." );
                    return -1;
                    break;
                }
                case mtk::EntityRank::EDGE:
                {
                    // check if requested vertex location makes sense
                    MORIS_ASSERT( aSignedLocalEntityIndex != 0,
                            "TET4_to_TET10::get_cell_local_vertex_index() - aSignedLocalEntityIndex is 1-based and must not be 0." );

                    MORIS_ASSERT( std::abs( aSignedLocalEntityIndex ) <= 6 && aVertexIndexOnEntity < 1,
                            "TET4_to_TET10::get_cell_local_vertex_index() - Index out of bounds. Only 3 edges with 1 new node per edge." );

                    // convert Signed 1-based index to unsigned 0-based index
                    uint tLocalEntityIndex = std::abs( aSignedLocalEntityIndex ) - 1;

                    // create look-up table for local vertex position
                    Matrix< IndexMat > tLookUp = {
                        { 4 },
                        { 5 },
                        { 6 },
                        { 7 },
                        { 8 },
                        { 9 }
                    };

                    return tLookUp( tLocalEntityIndex, aVertexIndexOnEntity );
                    break;
                }
                case mtk::EntityRank::FACE:
                {
                    MORIS_ERROR( false, "TET4_to_TET10::get_cell_local_vertex_index() - No new vertices on faces. This shouldn't be requested." );
                    return -1;
                    break;
                }
                case mtk::EntityRank::ELEMENT:
                {
                    MORIS_ERROR( false, "TET4_to_TET10::get_cell_local_vertex_index() - No new vertices inside cell. This shouldn't be requested." );
                    return -1;
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "TET4_to_TET10::get_cell_local_vertex_index() - Unknown Entity Rank." );
                    return -1;
                    break;
                }
            }
        }

        moris_index
        get_local_edge_index_based_on_vertex_indices( moris_index aFirstLocalVertexIndex, moris_index aSecondLocalVertexIndex ) const override
        {
            // check if input makes sense
            MORIS_ASSERT( aFirstLocalVertexIndex != aSecondLocalVertexIndex,
                    "TET4_to_TET10::get_local_edge_index_based_on_vertex_indices() - Vertices must have different indices" );

            // create LookUp
            Matrix< IndexMat > tLookUp = {
                { 0, 1, -3, 4 },
                { -1, 0, 2, -5 },
                { 3, -2, 0, 6 },
                { -4, 5, -6, 0 }
            };

            // return edge index (1-Based) and direction as sign
            return tLookUp( aFirstLocalVertexIndex, aSecondLocalVertexIndex );
        }
    };

    // -------------------------------------------------------------------------

    class TRI3_to_TRI10 : public Elevate_Order_Template
    {
      public:
        TRI3_to_TRI10() {}

        uint
        get_num_new_nodes() const override
        {
            return 7;
        }

        uint
        get_total_ig_verts() const override
        {
            return 10;
        }

        bool
        has_new_vertices_on_entity( mtk::EntityRank aEntityRank ) const override
        {
            switch ( aEntityRank )
            {
                case mtk::EntityRank::NODE:
                {
                    MORIS_ERROR( false, "TRI3_to_TRI10::has_new_nodes_on_entity() - No new vertices on vertices themselves. This shouldn't be requested." );
                    return false;
                    break;
                }
                case mtk::EntityRank::EDGE:
                {
                    return true;
                    break;
                }
                case mtk::EntityRank::FACE:
                {
                    return false;
                    break;
                }
                case mtk::EntityRank::ELEMENT:
                {
                    return false;
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "TRI3_to_TRI10::has_new_nodes_on_entity() - Unknown Entity Rank." );
                    return false;
                    break;
                }
            }
        }

        uint
        num_new_vertices_per_entity( mtk::EntityRank aEntityRank ) const override
        {
            switch ( aEntityRank )
            {
                case mtk::EntityRank::NODE:
                {
                    MORIS_ERROR( false, "TRI3_to_TRI10::num_new_vertices_per_entity() - No new vertices on vertices themselves. This shouldn't be requested." );
                    return 0;
                    break;
                }
                case mtk::EntityRank::EDGE:
                {
                    return 2;
                    break;
                }
                case mtk::EntityRank::FACE:
                {
                    return 0;
                    break;
                }
                case mtk::EntityRank::ELEMENT:
                {
                    return 1;
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "TRI3_to_TRI10::num_new_vertices_per_entity() - Unknown Entity Rank." );
                    return 0;
                    break;
                }
            }
        }

        Vector< Matrix< DDRMat > >
        get_new_vertex_parametric_coords_wrt_entity( mtk::EntityRank aEntityRank ) const override
        {
            switch ( aEntityRank )
            {
                case mtk::EntityRank::NODE:
                {
                    MORIS_ERROR( false, "TRI3_to_TRI10::get_new_vertex_parametric_coords_wrt_entity() - No new vertices on vertices themselves. This shouldn't be requested." );
                    return Vector< Matrix< DDRMat > >( 0 );
                    break;
                }
                case mtk::EntityRank::EDGE:
                {
                    return {
                        { { -1.0 / 3.0 } },
                        { { 1.0 / 3.0 } }
                    };
                    break;
                }
                case mtk::EntityRank::FACE:
                {
                    MORIS_ERROR( false, "TRI3_to_TRI10::get_new_vertex_parametric_coords_wrt_entity() - No new vertices on faces. This shouldn't be requested." );
                    return Vector< Matrix< DDRMat > >( 0 );
                    break;
                }
                case mtk::EntityRank::ELEMENT:
                {
                    real tOneThird = 1.0 / 3.0;
                    return { { { tOneThird, tOneThird, tOneThird } } };
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "TRI3_to_TRI10::get_new_vertex_parametric_coords_wrt_entity() - Unknown Entity Rank." );
                    return Vector< Matrix< DDRMat > >( 0 );
                    break;
                }
            }
        }

        mtk::CellTopology
        get_ig_cell_topology() const override
        {
            return mtk::CellTopology::TRI10;
        }

        uint
        get_num_spatial_dims() const override
        {
            return 2;
        }
    };

    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------

    class Elevate_Order_Interface : public Decomposition_Algorithm
    {

        // -------------------------------------------------------------------------

      private:
        std::shared_ptr< Elevate_Order_Template > mElevateOrderTemplate;

        Integration_Mesh_Generation_Data* mMeshGenerationData = nullptr;
        Decomposition_Data*               mDecompositionData  = nullptr;
        Cut_Integration_Mesh*             mCutIntegrationMesh = nullptr;
        moris::mtk::Mesh*                 mBackgroundMesh     = nullptr;
        Integration_Mesh_Generator*       mGenerator          = nullptr;

        // -------------------------------------------------------------------------

      public:
        Elevate_Order_Interface( Parameter_List& aParameterList, enum Subdivision_Method aSubdivisionMethod );

        ~Elevate_Order_Interface() override {}

        moris_index get_signature() const override;

        enum Decomposition_Algorithm_Type get_algorithm_type() const override;

        bool has_geometric_dependent_vertices() const override;

        Vector< moris_index > get_decomposed_cell_indices() override;

        void perform(
                Integration_Mesh_Generation_Data* aMeshGenerationData,
                Decomposition_Data*               aDecompositionData,
                Cut_Integration_Mesh*             aCutIntegrationMesh,
                moris::mtk::Mesh*                 aBackgroundMesh,
                Integration_Mesh_Generator*       aMeshGenerator ) override;

        void
        perform_impl_vertex_requests(
                Integration_Mesh_Generation_Data* aMeshGenerationData,
                Decomposition_Data*               aDecompositionData,
                Cut_Integration_Mesh*             aCutIntegrationMesh,
                moris::mtk::Mesh*                 aBackgroundMesh,
                Integration_Mesh_Generator*       aMeshGenerator ) override
        {
            MORIS_ERROR( false, "Elevate_Order_Interface::perform_impl_vertex_requests() - This virtual function is not implemented in this Child Class." );
        }

        void
        perform_impl_generate_mesh(
                Integration_Mesh_Generation_Data* aMeshGenerationData,
                Decomposition_Data*               aDecompositionData,
                Cut_Integration_Mesh*             aCutIntegrationMesh,
                moris::mtk::Mesh*                 aBackgroundMesh,
                Integration_Mesh_Generator*       aMeshGenerator ) override
        {
            MORIS_ERROR( false, "Elevate_Order_Interface::perform_impl_generate_mesh() - This virtual function is not implemented in this Child Class." );
        }

        // -------------------------------------------------------------------------

      private:
        bool
        make_vertex_requests(
                const std::shared_ptr< Edge_Based_Connectivity >& aEdgeConnectivity,
                const std::shared_ptr< Edge_Based_Ancestry >&     aIgEdgeAncestry,
                Vector< moris::mtk::Cell* >*                      aIgCells,
                Vector< Vector< moris_index > >*                  aCellToNewLocalVertexIndices );

        bool
        associate_new_vertices_with_cell_groups(
                const std::shared_ptr< Edge_Based_Connectivity >& aEdgeConnectivity,
                const std::shared_ptr< Edge_Based_Ancestry >&     aIgEdgeAncestry,
                Vector< moris::mtk::Cell* >*                      aBackgroundCellForEdge,
                Vector< moris::mtk::Cell* >*                      aIgCells );

        void
        create_higher_order_integration_cells(
                const std::shared_ptr< Edge_Based_Connectivity >& aEdgeConnectivity,
                const std::shared_ptr< Edge_Based_Ancestry >&     aIgEdgeAncestry,
                Vector< moris::mtk::Cell* >*                      aIgCells,
                Vector< Vector< moris_index > >*                  aCellToNewLocalVertexIndices );

        // -------------------------------------------------------------------------

        /**
         * @brief Creates unique id for an edge based on the IDs of its two end vertices
         * using the Cantor pairing function
         *
         * @param aEdgeVertices list of mtk::vertices on edge
         * @return moris_index unique id for edge
         */
        moris_index
        hash_edge( Vector< moris::mtk::Vertex* > const & aEdgeVertices );

        /**
         * @brief Creates unique id for a face based on the IDs of its three end corner vertices
         * by recursively applying the Cantor pairing function
         *
         * @param aFaceVertices list of mtk::vertices at face corners
         * @return moris_index unique id for face
         */
        moris_index
        hash_face( Vector< moris::mtk::Vertex* > const & aFaceVertices );

        /**
         * @brief swaps the values of two inidices
         *
         * @param aInd1 first index, value to be swaped with ...
         * @param aInd2 the second index
         */
        void
        swap_indices( moris_index& aInd1, moris_index& aInd2 );

        // -------------------------------------------------------------------------

        Matrix< DDRMat >
        compute_edge_vertex_global_coordinates( Vector< moris::mtk::Vertex* > const & aEdgeVertices, Matrix< DDRMat > aEdgeCoord );

        Matrix< DDRMat >
        compute_tri_vertex_global_coordinates( Vector< moris::mtk::Vertex* > const & aTriVertices, const Matrix< DDRMat >& aTriCoords );

        Matrix< DDRMat >
        compute_tet_vertex_global_coordinates( Vector< moris::mtk::Vertex* > const & aTetVertices, const Matrix< DDRMat >& aTetCoords );
    };

}    // namespace moris::xtk
#endif /* cl_XTK_Elevate_Order_Interface.hpp */
