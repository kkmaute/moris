/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Integration_Mesh_Generator.hpp
 *
 */

#ifndef MORIS_CL_XTK_Integration_Mesh_Generator_HPP_
#define MORIS_CL_XTK_Integration_Mesh_Generator_HPP_

#include <utility>

#include "cl_XTK_Model.hpp"
#include "cl_GEN_Geometry_Engine.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Cell_Info.hpp"
#include "cl_MTK_Mesh_Core.hpp"
#include "cl_XTK_Background_Mesh.hpp"
#include "cl_XTK_Decomposition_Data.hpp"
#include "cl_XTK_Cut_Integration_Mesh.hpp"
#include "cl_MTK_Cell_Info_Factory.hpp"
#include "cl_XTK_Cell_No_CM.hpp"
#include "cl_XTK_Subphase_Group.hpp"
#include "cl_Tracer.hpp"

namespace moris::xtk
{
    // ----------------------------------------------------------------------------------

    // typedef std::unordered_map< moris::moris_index, moris::moris_index > IndexMap;

    // ----------------------------------------------------------------------------------

    struct Integration_Mesh_Generation_Data
    {
        // number of child meshes
        moris_index mNumChildMeshes = 0;

        // outer vector: geometry index (if the geometry is inactive the cell is empty)
        // inner vector: active child mesh index by cut mesh index
        Vector< Vector< moris_index > > mIntersectedBackgroundCellIndex;

        // All background cell indices that will be delaunay triangulated (have surface points as defined by GEN geometries)
        Vector< moris_index > mDelaunayBgCellInds;

        // All background cells that will use regular subdivision
        Vector< moris_index > mRegularSubdivisionBgCellInds;

        // All intersected background cells (uniques removed from the concatenated version of mIntersectedBackgroundCellIndex)
        Vector< moris_index > mAllIntersectedBgCellInds;

        // all non-intersected background cells needed to triangulate all BG-cells in post
        Vector< moris_index > mAllNonIntersectedBgCellInds;

        // // this maps from the background cell index to the child mesh index
        // std::unordered_map< moris_index, moris_index > mIntersectedBackgroundCellIndexToChildMeshIndex;

        // Stores the surface points from every geometry for each background cell that has surface points
        // Outer vector: background cell index
        // Matrix: each column is a surface point, each row is a spatial dimension
        Vector< Matrix< DDRMat > > mDelaunayPoints;

        // Stores the geometry index for each surface point for each background cell
        // Outer vector: background cell index
        // Inner vector: Geometry index. Has as many entries as columns in mDelaunayPoints
        Vector< Vector< moris_index > > mDelaunayGeometryIndices;

    };    // struct Integration_Mesh_Generation_Data

    // ----------------------------------------------------------------------------------

    class Geometric_Query
    {
      private:
        // associated with a given child mesh
        mtk::EntityRank mQueryEntityRank = mtk::EntityRank::INVALID;

        Matrix< IndexMat > mQueryEntityToVertices;
        Matrix< IndexMat > mQueryEntityParamCoords;

        // parent cell is the query object when checking for intersection
        moris::mtk::Cell* mQueryParentCell;

        // for edge based questions
        moris_index                                mCurrentEdgeIndex = MORIS_INDEX_MAX;
        std::shared_ptr< Edge_Based_Connectivity > mEdgeConnectivity;

        Vector< moris::mtk::Cell* >*                  mAssociatedBgCellForEdge = nullptr;
        Vector< std::shared_ptr< IG_Vertex_Group > >* mAssociatedVertexGroup   = nullptr;
        Cut_Integration_Mesh*                         mCutIntegrationMesh      = nullptr;

        // parent info
        moris::mtk::Cell* mParentCell;

        // row-based and indexed based on the indexing in mQueryEntityToVertices.
        Vector< std::shared_ptr< Matrix< DDRMat > > >* mQueryEntityIndexedCoordinates = nullptr;

        // geometric index
        moris_index mGeometricIndex;

      public:
        Geometric_Query()
                : mQueryParentCell( nullptr )
                , mParentCell( nullptr )
                , mGeometricIndex( MORIS_INDEX_MAX )
        {
        }
        ~Geometric_Query() {}

        void
        set_query_entity_rank( mtk::EntityRank aEntityRank )
        {
            mQueryEntityRank = aEntityRank;
        }

        void
        set_edge_connectivity( std::shared_ptr< Edge_Based_Connectivity > aEdgeConnectivity )
        {
            mEdgeConnectivity = std::move( aEdgeConnectivity );
        }

        void
        set_edge_associated_background_cell( Vector< moris::mtk::Cell* >* aAssociatedBgCellsForEdge )
        {
            mAssociatedBgCellForEdge = aAssociatedBgCellsForEdge;
        }

        void
        set_associated_vertex_group( Vector< std::shared_ptr< IG_Vertex_Group > >* aAssociatedVertexGroup )
        {
            mAssociatedVertexGroup = aAssociatedVertexGroup;
        }

        void
        set_cut_integration_mesh( Cut_Integration_Mesh* aCutIntegrationMesh )
        {
            mCutIntegrationMesh = aCutIntegrationMesh;
        }

        void
        set_current_edge_index( moris_index aCurrentIndex )
        {
            mCurrentEdgeIndex = aCurrentIndex;
            mQueryEntityToVertices.resize( 1, 2 );
            mQueryEntityToVertices( 0 ) = mEdgeConnectivity->mEdgeVertices( mCurrentEdgeIndex )( 0 )->get_index();
            mQueryEntityToVertices( 1 ) = mEdgeConnectivity->mEdgeVertices( mCurrentEdgeIndex )( 1 )->get_index();
        }

        void
        set_parametric_coordinate_size( moris_index aParametricSize )
        {
            // hard coded for edges right now
            mQueryEntityParamCoords.resize( 2, aParametricSize );
        }

        void
        set_parent_cell( moris::mtk::Cell* aParentCell )
        {
            mParentCell = aParentCell;
        }

        void
        set_query_cell( moris::mtk::Cell* aQueryParentCell )
        {
            mQueryParentCell       = aQueryParentCell;
            mQueryEntityToVertices = aQueryParentCell->get_vertex_inds();
        }

        void
        set_query_edge_connectivity( const std::shared_ptr< Edge_Based_Connectivity >& aEdgeBasedConnectivity )
        {
        }

        void
        set_coordinates_matrix( Vector< std::shared_ptr< Matrix< DDRMat > > >* aCoordinates )
        {
            mQueryEntityIndexedCoordinates = aCoordinates;
        }

        void
        set_geometric_index( moris_index aGeometricIndex )
        {
            mGeometricIndex = aGeometricIndex;
        }

        moris_index
        get_geometric_index() const
        {
            return mGeometricIndex;
        }

        mtk::EntityRank
        get_query_entity_rank() const
        {
            return mQueryEntityRank;
        }

        Matrix< IndexMat > const &
        get_query_entity_to_vertex_connectivity() const
        {
            return mQueryEntityToVertices;
        }

        Vector< std::shared_ptr< Matrix< DDRMat > > >*
        get_query_indexed_coordinates() const
        {
            return mQueryEntityIndexedCoordinates;
        }

        Matrix< DDRMat >
        get_vertex_local_coord_wrt_parent_entity( moris_index aVertexIndex ) const
        {
            std::shared_ptr< IG_Vertex_Group > tCurrentVertexGroup = ( *mAssociatedVertexGroup )( mCurrentEdgeIndex );
            return *tCurrentVertexGroup->get_vertex_local_coords( aVertexIndex );
        }

        mtk::EntityRank
        get_query_parent_entity_rank() const
        {
            return mtk::EntityRank::ELEMENT;
        }

        Matrix< IndexMat >
        get_query_parent_entity_connectivity() const
        {
            return mParentCell->get_vertex_inds();
        }

        Matrix< DDRMat >
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

        mtk::Geometry_Type get_geometry_type()
        {
            return mParentCell->get_geometry_type();
        }

        mtk::Interpolation_Order get_interpolation_order()
        {
            return mParentCell->get_interpolation_order();
        }

    };    // struct Geometric_Query

    // ----------------------------------------------------------------------------------

    class Decomposition_Algorithm;
    class Integration_Mesh_Generator
    {
        // ----------------------------------------------------------------------------------

      private:
        Model*                            mXTKModel       = nullptr;
        moris::gen::Geometry_Engine*      mGeometryEngine = nullptr;
        Vector< enum Subdivision_Method > mSubdivisionMethods;

        bool mOutputCutIgMesh = false;

      public:
        // ----------------------------------------------------------------------------------

        Integration_Mesh_Generator();

        // ----------------------------------------------------------------------------------

        Integration_Mesh_Generator(
                xtk::Model*                              aXTKModelPtr,
                const Vector< enum Subdivision_Method >& aMethods );

        // ----------------------------------------------------------------------------------

        ~Integration_Mesh_Generator();

        // ----------------------------------------------------------------------------------

        std::shared_ptr< Cut_Integration_Mesh >
        perform();

        // ----------------------------------------------------------------------------------

        moris::gen::Geometry_Engine*
        get_geom_engine();

        // ----------------------------------------------------------------------------------

        uint
        get_spatial_dim();

        // ----------------------------------------------------------------------------------

        uint
        get_ig_mesh_order();

        // ----------------------------------------------------------------------------------

        enum Subdivision_Method
        determine_reg_subdivision_template();

        // ----------------------------------------------------------------------------------

        enum Subdivision_Method
        determine_order_elevation_template();

        // ----------------------------------------------------------------------------------

        void
        set_cut_IG_mesh_output( bool aOutputCutIgMesh )
        {
            mOutputCutIgMesh = aOutputCutIgMesh;
        }

        // ----------------------------------------------------------------------------------

        bool
        get_cut_IG_mesh_output()
        {
            return mOutputCutIgMesh;
        }

        // ----------------------------------------------------------------------------------

        /**
         * @brief checks whether all intersected background cells are on the same level. The resultant bool is populated
         * in the cut integration mesh mChildMeshSameLevel. ultimately This triggers a remeshing of HMR in the background mesh
         *
         * @param aMeshGenerationData
         * @param aBackgroundMesh
         */

        void
        check_intersected_background_cell_levels(
                Integration_Mesh_Generation_Data& aMeshGenerationData,
                Cut_Integration_Mesh*             aCutIntegrationMesh,
                moris::mtk::Mesh*                 aBackgroundMesh );

        // ----------------------------------------------------------------------------------

        bool
        determine_intersected_background_cells(
                Integration_Mesh_Generation_Data& aMeshGenerationData,
                Cut_Integration_Mesh*             aCutIntegrationMesh,
                moris::mtk::Mesh*                 aBackgroundMesh );

        // ----------------------------------------------------------------------------------

        bool
        determine_non_intersected_background_cells(
                Integration_Mesh_Generation_Data& aMeshGenerationData,
                Cut_Integration_Mesh*             aCutIntegrationMesh,
                moris::mtk::Mesh*                 aBackgroundMesh );

        // ----------------------------------------------------------------------------------

        void
        deduce_interfaces(
                Cut_Integration_Mesh*                              aCutIntegrationMesh,
                const std::shared_ptr< Facet_Based_Connectivity >& aFacetConnectivity,
                Vector< moris_index >&                             aInterfaces );

        // ----------------------------------------------------------------------------------

        void
        construct_bulk_phase_to_bulk_phase_interface(
                Cut_Integration_Mesh*                                      aCutIntegrationMesh,
                Vector< moris_index >&                                     aInterfaces,
                Vector< Vector< std::shared_ptr< IG_Cell_Side_Group > > >& aInterfaceBulkPhaseToBulk );

        // ----------------------------------------------------------------------------------

        void
        construct_bulk_phase_to_bulk_phase_dbl_side_interface(
                Cut_Integration_Mesh*                                             aCutIntegrationMesh,
                Vector< moris_index >&                                            aInterfaces,
                Vector< Vector< std::shared_ptr< IG_Cell_Double_Side_Group > > >& aDblSideInterfaceBulkPhaseToBulk );

        // ----------------------------------------------------------------------------------

        std::string
        get_interface_side_set_name(
                moris_index aGeomIndex,
                moris_index aBulkPhaseIndex0,
                moris_index aBulkPhaseIndex1 );

        // ----------------------------------------------------------------------------------

        void
        construct_interface_sets(
                Cut_Integration_Mesh*                                      aCutIntegrationMesh,
                Vector< Vector< std::shared_ptr< IG_Cell_Side_Group > > >& aInterfaceBulkPhaseToBulk );

        // ----------------------------------------------------------------------------------

        void
        construct_bulk_phase_cell_groups(
                Cut_Integration_Mesh*                       aCutIntegrationMesh,
                Vector< std::shared_ptr< IG_Cell_Group > >& aBulkPhaseCellGroups );

        // ----------------------------------------------------------------------------------

        void
        construct_bulk_phase_blocks(
                Cut_Integration_Mesh*                       aCutIntegrationMesh,
                Vector< std::shared_ptr< IG_Cell_Group > >& aBulkPhaseCellGroups );

        // ----------------------------------------------------------------------------------

        bool
        allocate_child_meshes( Integration_Mesh_Generation_Data& aMeshGenerationData,
                Cut_Integration_Mesh*                            aCutIntegrationMesh,
                moris::mtk::Mesh*                                aBackgroundMesh );

        // ----------------------------------------------------------------------------------

        void
        collect_vertex_groups_for_background_cells(
                Integration_Mesh_Generation_Data*             aMeshGenerationData,
                Cut_Integration_Mesh*                         aCutIntegrationMesh,
                Vector< moris::mtk::Cell* >*                  aBackgroundCells,
                Vector< std::shared_ptr< IG_Vertex_Group > >* aVertexGroups );

        // ----------------------------------------------------------------------------------

        /**
         * this function checks if first vertex in first map (first input)
         * exists in last two inputs and returns MORIS_MAX_INDEX if not
         */
        moris_index
        edge_exists(
                Vector< moris::mtk::Vertex* >&                  aVerticesOnEdge,
                std::unordered_map< moris_index, moris_index >& aLocaLVertexMap,
                Vector< Vector< uint > >&                       aVertexToEdge,
                Vector< Vector< moris::mtk::Vertex* > >&        aFullEdgeVertices );

        // ----------------------------------------------------------------------------------

        moris_index
        facet_exists(
                Vector< moris::mtk::Vertex* >&                  aVerticesOnFacet,
                std::unordered_map< moris_index, moris_index >& aLocaLVertexMap,
                Vector< Vector< uint > >&                       aVertexToFacet,
                Vector< Vector< moris::mtk::Vertex* > >&        aFullFacetVertices );

        // ----------------------------------------------------------------------------------

        void
        construct_ig_cell_groups_facet_conn(
                Vector< std::shared_ptr< IG_Cell_Group > >&            aActiveIgCellGroups,
                Vector< std::shared_ptr< Facet_Based_Connectivity > >& aGroupFaceConn );

        // ----------------------------------------------------------------------------------

        void
        generate_cell_neighborhood(
                Vector< moris::mtk::Cell* >&                             aCells,
                const std::shared_ptr< Facet_Based_Connectivity >&       aFaceConnectivity,
                const std::shared_ptr< Cell_Neighborhood_Connectivity >& aNeighborhood );

        // ----------------------------------------------------------------------------------

        void
        create_facet_from_element_to_node(
                Vector< moris::mtk::Cell* >&                       aCells,
                const std::shared_ptr< Facet_Based_Connectivity >& aFaceConnectivity,
                const Cut_Integration_Mesh*                        aCutIntegrationMesh );

        [[deprecated]]
        void
        create_facet_from_element_to_node(
                Vector< moris::mtk::Cell* >&                       aCells,
                const std::shared_ptr< Facet_Based_Connectivity >& aFaceConnectivity );

        // ----------------------------------------------------------------------------------

        void
        create_vertex_facet_connectivity(
                const std::shared_ptr< Facet_Based_Connectivity >& aFaceConnectivity,
                Cut_Integration_Mesh*                              aCutIntegrationMesh );

        // ----------------------------------------------------------------------------------

        void
        create_edges_from_element_to_node(
                Vector< moris::mtk::Cell* >                       aCells,
                const std::shared_ptr< Edge_Based_Connectivity >& aEdgeConnectivity );

        // ----------------------------------------------------------------------------------

        void
        select_background_cell_for_edge(
                const std::shared_ptr< Edge_Based_Connectivity >& aEdgeBasedConnectivity,
                Cut_Integration_Mesh*                             aCutIntegrationMesh,
                Vector< moris::mtk::Cell* >&                      aBackgroundCellForEdge );

        // ----------------------------------------------------------------------------------

        void
        select_background_cell_for_facet(
                const std::shared_ptr< Facet_Based_Connectivity >& aEdgeBasedConnectivity,
                Cut_Integration_Mesh*                              aCutIntegrationMesh,
                Vector< moris::mtk::Cell* >&                       aBackgroundCellForEdge );

        // ----------------------------------------------------------------------------------

        void
        deduce_facet_ancestry(
                Cut_Integration_Mesh*                              aCutIntegrationMesh,
                moris::mtk::Mesh*                                  aBackgroundMesh,
                const std::shared_ptr< Facet_Based_Connectivity >& aIgCellGroupFacetConnectivity,
                Vector< moris::mtk::Cell* > const &                aParentCellForDeduction,
                const std::shared_ptr< Facet_Based_Ancestry >&     aIgFacetAncestry );

        // ----------------------------------------------------------------------------------

        void
        compute_bg_facet_to_child_facet_connectivity(
                Cut_Integration_Mesh*                                      aCutIntegrationMesh,
                moris::mtk::Mesh*                                          aBackgroundMesh,
                const std::shared_ptr< Facet_Based_Connectivity >&         aIgCellGroupFacetConnectivity,
                const std::shared_ptr< Facet_Based_Ancestry >&             aIgFacetAncestry,
                Vector< std::shared_ptr< Vector< moris::moris_index > > >& aBgFacetToIgFacet );

        // ----------------------------------------------------------------------------------

        void
        deduce_edge_ancestry(
                Cut_Integration_Mesh*                             aCutIntegrationMesh,
                moris::mtk::Mesh*                                 aBackgroundMesh,
                const std::shared_ptr< Edge_Based_Connectivity >& aIgCellGroupEdgeConnectivity,
                Vector< moris::mtk::Cell* > const &               aParentCellForDeduction,
                const std::shared_ptr< Edge_Based_Ancestry >&     aIgEdgeAncestry );

        // ----------------------------------------------------------------------------------

        void
        commit_new_ig_vertices_to_cut_mesh(
                Integration_Mesh_Generation_Data* aMeshGenerationData,
                Decomposition_Data*               aDecompositionData,
                Cut_Integration_Mesh*             aCutIntegrationMesh,
                moris::mtk::Mesh*                 aBackgroundMesh,
                Decomposition_Algorithm*          aDecompositionAlgorithm );

        // ----------------------------------------------------------------------------------

        void
        link_new_vertices_to_geometry_engine(
                Decomposition_Data*      aDecompositionData,
                Decomposition_Algorithm* aDecompAlg );

        // ----------------------------------------------------------------------------------

        void
        commit_new_ig_cells_to_cut_mesh(
                Integration_Mesh_Generation_Data* aMeshGenerationData,
                Decomposition_Data*               aDecompositionData,
                Cut_Integration_Mesh*             aCutIntegrationMesh,
                moris::mtk::Mesh*                 aBackgroundMesh,
                Decomposition_Algorithm*          aDecompositionAlgorithm );

        // ----------------------------------------------------------------------------------

        void
        identify_and_construct_subphases(
                Integration_Mesh_Generation_Data*                        aMeshGenerationData,
                Cut_Integration_Mesh*                                    aCutIntegrationMesh,
                moris::mtk::Mesh*                                        aBackgroundMesh,
                const std::shared_ptr< Cell_Neighborhood_Connectivity >& aCutNeighborhood );

        // ----------------------------------------------------------------------------------

        void
        construct_subphases_on_triangulated_non_cut_cells(
                Integration_Mesh_Generation_Data* aMeshGenerationData,
                Cut_Integration_Mesh*             aCutIntegrationMesh,
                moris::mtk::Mesh*                 aBackgroundMesh );

        // ----------------------------------------------------------------------------------

        void
        construct_subphase_neighborhood(
                Cut_Integration_Mesh*                                        aCutIntegrationMesh,
                moris::mtk::Mesh*                                            aBackgroundMesh,
                const std::shared_ptr< Facet_Based_Connectivity >&           aFacetConnectivity,
                Vector< std::shared_ptr< Vector< moris::moris_index > > >*   aBgFacetToChildFacet,
                const std::shared_ptr< Subphase_Neighborhood_Connectivity >& aSubphaseNeighborhood );

        // ----------------------------------------------------------------------------------

        void
        collect_subphases_attached_to_facet_on_cell(
                Cut_Integration_Mesh*                                  aCutIntegrationMesh,
                mtk::Cell const *                                      aBGCell,
                moris_index                                            aFacetOrdinal,
                const std::shared_ptr< Facet_Based_Connectivity >&     aFacetConnectivity,
                const std::shared_ptr< Vector< moris::moris_index > >& aBgFacetToChildrenFacets,
                Vector< moris_index >&                                 aSubphaseIndices,
                Vector< moris_index >&                                 aRepresentativeIgCells,
                Vector< moris_index >&                                 aRepresentativeIgCellsOrdinal );

        // ----------------------------------------------------------------------------------

        /**
         * @brief Gets the subphase information of all IG cells in aBGCell that are overlapping with the given facet from another BG cell
         */
        void
        collect_subphases_of_overlapping_facets_on_cell(
                Cut_Integration_Mesh*                                  aCutIntegrationMesh,
                moris_index                                            aMyIGFacetIndex,
                mtk::Cell const *                                      aBGCell,
                moris_index                                            aFacetOrdinal,
                const std::shared_ptr< Facet_Based_Connectivity >&     aFacetConnectivity,
                const std::shared_ptr< Vector< moris::moris_index > >& aBgFacetToChildrenFacets,
                Vector< moris_index >&                                 aSubphaseIndices,
                Vector< moris_index >&                                 aRepresentativeIgCells,
                Vector< moris_index >&                                 aRepresentativeIgCellsOrdinal );

        // ----------------------------------------------------------------------------------


        // ----------------------------------------------------------------------------------

        /**
         * @brief Gets the 2D coordinates in plane coordinates of a facet on a specific side ordinal
         *
         * @param aMyIGFacetIndex index of the facet in the cut integration mesh
         * @param aFacetConnectivity facet connectivity of the cut integration mesh, which holds information about the vertices of the facet
         * @param aFacetOrdinal background cell side ordinal that defines the plane for the facet
         */
        static Matrix< DDRMat >
        get_in_plane_facet_coordinates(
                moris_index                                        aMyIGFacetIndex,
                const std::shared_ptr< Facet_Based_Connectivity >& aFacetConnectivity,
                moris_index                                        aFacetOrdinal );

        // ----------------------------------------------------------------------------------

        void collect_ig_cells_and_side_ords_on_bg_facet(
                Cut_Integration_Mesh*        aCutIntegrationMesh,
                moris::moris_index           aBackgroundFacetIndex,
                Vector< moris::mtk::Cell* >& aIgCell,
                Vector< moris_index >&       aIgCellSideOrds );

        // ----------------------------------------------------------------------------------

        void append_subphase_information_on_facet(
                moris_index                                                     aChildCellFacetIndex,
                Cut_Integration_Mesh*                                           aCutIntegrationMesh,
                mtk::Cell const *                                               aBGCell,
                const std::shared_ptr< Facet_Based_Connectivity >&              aFacetConnectivity,
                Vector< moris_index >&                                          aSubphaseIndices,
                Vector< moris_index >&                                          aRepresentativeIgCells,
                Vector< moris_index >&                                          aRepresentativeIgCellsOrdinal,
                std::map< std::pair< moris_index, moris_index >, moris_index >& aSubphaseMap );
        // ----------------------------------------------------------------------------------

        /**
         * @brief finds the range of valid Subphase IDs on each processor and assigns them to the SPs
         *
         * @param aCutIntegrationMesh pointer to the Cut_Integration_Mesh
         * @param aBackgroundMesh pointer to the HMR Lagrange Background mesh
         */
        void
        assign_subphase_glob_ids(
                Cut_Integration_Mesh* aCutIntegrationMesh,
                moris::mtk::Mesh*     aBackgroundMesh );

        // ----------------------------------------------------------------------------------

        /**
         * @brief assign IDs to the subphases owned by the executing processor
         *
         * @param aCutIntegrationMesh pointer to the Cut_Integration_Mesh
         * @param tFirstID first subphase ID to be used by the executing processor
         */
        void
        assign_IDs_to_owned_subphases(
                Cut_Integration_Mesh* aCutIntegrationMesh,
                moris_id              aFirstID );

        // ----------------------------------------------------------------------------------

        /**
         * @brief Prepares lists of Subphases the executing proc has access to but doesn't own for each processor that actually owns these subphases.
         * Additionally, the parent IP cell IDs, IG cell group IDs, and the number of IG cells in each SP will be collected.
         *
         * @param aCutIntegrationMesh pointer to the Cut_Integration_Mesh
         * @param aNotOwnedSubphasesToProcs outer cell: owner proc index in XTK comm-table || inner cell: list of non-owned SPs whose IDs need to be requested from that owning proc
         * @param aParentCellIds outer cell: owner proc index in XTK comm-table || inner cell: parent Cell IDs corresponding to the SPs in the lists above
         * @param aChildCellIds outer cell: owner proc index in XTK comm-table || inner cell: IDs of the first IG cell in the IG cell groups corresponding to the SPs in the lists above
         * @param aNumChildCellsInSubphases outer cell: owner proc index in XTK comm-table || inner cell: number of IG cells in each of the SPs in the lists above
         */
        void
        prepare_requests_for_not_owned_subphase_IDs(
                Cut_Integration_Mesh*            aCutIntegrationMesh,
                Vector< Vector< moris_index > >& aNotOwnedSubphasesToProcs,
                Vector< Matrix< IdMat > >&       aParentCellIds,
                Vector< Matrix< IdMat > >&       aFirstChildIgCellIds,
                Vector< Matrix< IdMat > >&       aNumChildCellsInSubphases );

        // ----------------------------------------------------------------------------------

        /**
         * @brief find and answer with the IDs of the subphases owned by the executing processor and requested from other processors
         *
         * @param aCutIntegrationMesh pointer to the Cut_Integration_Mesh
         * @param aSubphaseIds outer cell: owner proc index in XTK comm-table || inner cell: IDs of the Subphase corresponding to the SPs specified by the above lists
         * @param aReceivedParentCellIds outer cell: owner proc index in XTK comm-table || inner cell: parent Cell IDs corresponding to the SPs whose IDs are requested
         * @param aReceivedFirstChildIgCellIds outer cell: owner proc index in XTK comm-table || inner cell: IDs of the first IG cell in the IG cell groups corresponding to the SPs whose IDs are requested
         * @param aReceivedNumChildCellsInSubphases outer cell: owner proc index in XTK comm-table || inner cell: number of IG cells in each of the SPs whose IDs are requested
         */
        void
        prepare_answers_for_owned_subphase_IDs(
                Cut_Integration_Mesh*             aCutIntegrationMesh,
                Vector< Matrix< IdMat > >&        aSubphaseIds,
                Vector< Matrix< IdMat > > const & aReceivedParentCellIds,
                Vector< Matrix< IdMat > > const & aReceivedFirstChildIgCellIds,
                Vector< Matrix< IdMat > > const & aReceivedNumChildCellsInSubphases );

        // ----------------------------------------------------------------------------------

        /**
         * @brief Assign the received subphase IDs to on the executing processor
         *
         * @param aCutIntegrationMesh pointer to the Cut_Integration_Mesh
         * @param tNotOwnedSubphasesToProcs outer cell: owner proc index in XTK comm-table || inner cell: list of non-owned SPs whose IDs need to be requested from that owning proc
         * @param tReceivedSubphaseIds outer cell: owner proc index in XTK comm-table || inner cell: IDs of the Subphase corresponding to the SPs specified by the above lists
         */
        void
        handle_requested_subphase_ID_answers(
                Cut_Integration_Mesh*                   aCutIntegrationMesh,
                Vector< Vector< moris_index > > const & tNotOwnedSubphasesToProcs,
                Vector< Matrix< IdMat > > const &       tReceivedSubphaseIds );

        // ----------------------------------------------------------------------------------
        // ----------------------------------------------------------------------------------

        /**
         * @brief finds the range of valid Subphase group IDs on each processor and assigns them to the SPGs
         *
         * @param aCutIntegrationMesh pointer to the Cut_Integration_Mesh
         * @param aBackgroundMesh pointer to the HMR Lagrange Background mesh
         * @param aBsplineMeshInfo pointer to the B-spline mesh info for the B-spline mesh wrt. which the SPG IDs should be assigned
         */
        void
        assign_subphase_group_glob_ids(
                Cut_Integration_Mesh* aCutIntegrationMesh,
                Bspline_Mesh_Info*    aBsplineMeshInfo );

        // ----------------------------------------------------------------------------------

        /**
         * @brief Prepares lists of subphase groups the executing proc has access to but doesn't own for each processor that actually owns these SPGs.
         *
         * @param aCutIntegrationMesh pointer to the Cut_Integration_Mesh
         * @param aBsplineMeshInfo pointer to the B-spline mesh info to be treated
         * @param aNotOwnedSpgsToProcs outer cell: owner proc index in XTK comm-table || inner cell: list of non-owned SPGs whose IDs need to be requested from that owning proc
         * @param aSubphaseIds outer cell: owner proc index in XTK comm-table || inner cell: first SP IDs in SPGs as listed in array above
         * @param aNumSpsInSpg outer cell: owner proc index in XTK comm-table || inner cell: number of SPs in the respective SPGs listed in array above
         */
        void
        prepare_requests_for_not_owned_subphase_group_IDs(
                Cut_Integration_Mesh*            aCutIntegrationMesh,
                Bspline_Mesh_Info*               aBsplineMeshInfo,
                Vector< Vector< moris_index > >& aNotOwnedSpgsToProcs,
                Vector< Matrix< IdMat > >&       aSubphaseIds,
                Vector< Matrix< IdMat > >&       aNumSpsInSpg );

        // ----------------------------------------------------------------------------------

        /**
         * @brief Find and answer with the IDs of the subphase groups owned by the executing processor and requested from other processors
         *
         * @param aCutIntegrationMesh pointer to the Cut_Integration_Mesh
         * @param aBsplineMeshInfo pointer to the B-spline mesh info to be treated
         * @param aSubphaseGroupIds output; outer cell: owner proc index in XTK comm-table || inner cell: list of SPG ID answers found to the requests below
         * @param aReceivedSubphaseIds input; outer cell: owner proc index in XTK comm-table || inner cell: first SP IDs in SPGs for which IDs were requested
         * @param aReceivedNumSpsInSpg input; outer cell: owner proc index in XTK comm-table || inner cell: number of SPs in the respective SPGs listed in array above
         */
        void
        prepare_answers_for_owned_subphase_group_IDs(
                Cut_Integration_Mesh*             aCutIntegrationMesh,
                Bspline_Mesh_Info*                aBsplineMeshInfo,
                Vector< Matrix< IdMat > >&        aSubphaseGroupIds,
                Vector< Matrix< IdMat > > const & aReceivedSubphaseIds,
                Vector< Matrix< IdMat > > const & aReceivedNumSpsInSpg );

        // ----------------------------------------------------------------------------------

        /**
         * @brief Assign the received subphase group IDs on the executing processor
         *
         * @param aCutIntegrationMesh pointer to the Cut_Integration_Mesh
         * @param aBsplineMeshInfo pointer to the B-spline mesh info to be treated
         * @param aNotOwnedSpgsToProcs outer cell: owner proc index in XTK comm-table || inner cell: list of non-owned SPGs whose IDs need to be requested from that owning proc
         * @param aReceivedSubphaseGroupIds output; outer cell: owner proc index in XTK comm-table || inner cell: list of SPG ID answers received from other procs
         */
        void
        handle_requested_subphase_group_ID_answers(
                Cut_Integration_Mesh*                   aCutIntegrationMesh,
                Bspline_Mesh_Info*                      aBsplineMeshInfo,
                Vector< Vector< moris_index > > const & aNotOwnedSpgsToProcs,
                Vector< Matrix< IdMat > > const &       aReceivedSubphaseGroupIds );

        // ----------------------------------------------------------------------------------
        // ----------------------------------------------------------------------------------

        /**
         * @brief perform flood fill on IG cell group and assign a subphase index to every IG cell
         *
         * @param aCutIntegrationMesh pointer to XTK Cut_Integration_Mesh
         * @param aCutNeighborhood Cell_Neighborhood_Connectivity of proc local mesh
         * @param aIgCellGroup pointer to IG_Cell_Group to perform flood fill on
         * @param aMaxValueAssigned (return value) number of subphases assigned during flood-fill (= max bg-cell local index of subphases assigned)
         * @return Matrix< IndexMat > List (vector) of subphase indices. Index in list corresponds to IG cell group local index of IG cell
         */
        Matrix< IndexMat >
        flood_fill_ig_cell_group(
                Cut_Integration_Mesh*                                    aCutIntegrationMesh,
                const std::shared_ptr< Cell_Neighborhood_Connectivity >& aCutNeighborhood,
                const std::shared_ptr< IG_Cell_Group >&                  aIgCellGroup,
                moris_index&                                             aMaxValueAssigned );

        // ----------------------------------------------------------------------------------

        void
        extract_cells_from_cell_groups(
                Vector< std::shared_ptr< IG_Cell_Group > > const & aCellGroups,
                Vector< moris::mtk::Cell* >&                       aCellsInGroups );

        // ----------------------------------------------------------------------------------

        void
        compute_ig_cell_bulk_phase( Cut_Integration_Mesh* aCutIntegrationMesh );

        moris_index
        deduce_ig_cell_bulk_phase_from_facets(
                moris::mtk::Cell const * aCell,
                Cut_Integration_Mesh*    aCutIntegrationMesh );

        // ----------------------------------------------------------------------------------

        moris_index
        get_max_index( Vector< moris::mtk::Cell* >& aCells );

        // ----------------------------------------------------------------------------------

        void
        assign_node_requests_identifiers(
                Decomposition_Data&   aDecompData,
                Cut_Integration_Mesh* aCutIntegrationMesh,
                moris::mtk::Mesh*     aBackgroundMesh );

        // ----------------------------------------------------------------------------------

        void
        sort_new_node_requests_by_owned_and_not_owned(
                Decomposition_Data&                       tDecompData,
                Cut_Integration_Mesh*                     aCutIntegrationMesh,
                moris::mtk::Mesh*                         aBackgroundMesh,
                Vector< uint >&                           aOwnedRequests,
                Vector< Vector< uint > >&                 aNotOwnedRequests,
                Vector< uint >&                           aProcRanks,
                std::unordered_map< moris_id, moris_id >& aProcRankToIndexInData );

        // ----------------------------------------------------------------------------------

        void
        assign_owned_request_id(
                Decomposition_Data&    aDecompData,
                Vector< uint > const & aOwnedRequest,
                moris::moris_id&       aNodeId );

        // ----------------------------------------------------------------------------------

        void
        setup_outward_requests(
                Decomposition_Data const &                aDecompData,
                moris::mtk::Mesh*                         aBackgroundMesh,
                Vector< Vector< uint > > const &          aNotOwnedRequests,
                Vector< uint > const &                    aProcRanks,
                std::unordered_map< moris_id, moris_id >& aProcRankToIndexInData,
                Vector< Matrix< IndexMat > >&             aOutwardRequests );

        // ----------------------------------------------------------------------------------

        void
        prepare_request_answers(
                Decomposition_Data&                  aDecompData,
                moris::mtk::Mesh*                    aBackgroundMesh,
                Vector< Matrix< IndexMat > > const & aReceiveData,
                Vector< Matrix< IndexMat > >&        aRequestAnswers );

        // ----------------------------------------------------------------------------------

        void
        handle_received_request_answers(
                Decomposition_Data&                  aDecompData,
                moris::mtk::Mesh*                    aBackgroundMesh,
                Vector< Matrix< IndexMat > > const & aRequests,
                Vector< Matrix< IndexMat > > const & aRequestAnswers,
                moris::moris_id&                     aNodeId );

        // ----------------------------------------------------------------------------------

        moris::uint
        verbosity_level();

        /**
         * @brief Removes subphases and deletes contained ig cells
         *
         * @param aSubphasesToRemove subphase indexes to delete
         */

        void
        remove_subphases_from_cut_mesh( Vector< moris_index > const & aSubphasesToRemove );

        // ----------------------------------------------------------------------------------
        // Functions for Constructing Subphase Groups and their neighborhood
        // ----------------------------------------------------------------------------------

        /**
         * @brief check whether it is necessary to construct SPGs and their connectivity
         */
        bool
        check_construct_subphase_groups();

        // ----------------------------------------------------------------------------------

        /**
         * @brief establish relationship to/ information about B-spline background mesh
         *
         * @param aCutIntegrationMesh pointer to the cut integration mesh which carries the information
         * @param aLagrangeMesh pointer to the Lagrange mesh for which the information needs to be established
         * @param aBsplineMeshInfo pointer to the B-spline mesh info to be used for construction
         * @param aMeshIndex discretization mesh index (NOT List index!, since this index is used to get B-spline mesh info from HMR)
         */
        void
        establish_bspline_mesh_info(
                Cut_Integration_Mesh* aCutIntegrationMesh,
                moris::mtk::Mesh*     aLagrangeMesh,
                Bspline_Mesh_Info*    aBsplineMeshInfo,
                const moris_index     aMeshIndex );

        // ----------------------------------------------------------------------------------

        /**
         * @brief find groups of subphases that belong together and the side ordinals through which it is connected to neighboring SPGs
         *
         * @param aCutIntegrationMesh pointer to cut integration mesh holding B-spline mesh info and SPGs
         * @param aBsplineMeshInfo pointer to the B-spline mesh info to be used for construction
         * @param aMeshListIndex discretization mesh list index (i.e. index of the B-spline mesh to be treated in the list of mesh indices to be enriched)
         */
        void
        construct_subphase_groups(
                Cut_Integration_Mesh* aCutIntegrationMesh,
                Bspline_Mesh_Info*    aBsplineMeshInfo,
                const moris_index     aMeshListIndex );

        // ----------------------------------------------------------------------------------

        /**
         * @brief finds all owned and non-owned SPGs on the current proc
         *
         * @param aCutIntegrationMesh pointer to cut integration mesh holding B-spline mesh info and SPGs
         * @param aLagrangeMesh pointer to Lagrange mesh
         * @param aBsplineMeshInfo pointer to the B-spline mesh info to be used for construction
         * @param aBsplineMeshListIndex index of the B-spline mesh within the list of mesh-indices to enrich
         */
        void
        communicate_subphase_groups(
                Cut_Integration_Mesh* aCutIntegrationMesh,
                Bspline_Mesh_Info*    aBsplineMeshInfo,
                const moris_index     aBsplineMeshListIndex );

        // ----------------------------------------------------------------------------------

        /**
         * @brief find the subphase indices within a given B-spline element on a given B-spline mesh
         *
         * @param aCutIntegrationMesh pointer to cut integration mesh
         * @param aMeshListIndex index of the B-spline mesh to be treated within the list of discretization mesh indices to be enriched
         * @param aCurrentBspCellIndex index of the B-spline element
         * @param aSubPhaseIndices output: list of subphase indices to be filled
         */
        void
        get_subphase_indices_in_bspline_cell(
                Cut_Integration_Mesh* aCutIntegrationMesh,
                moris_index           aMeshListIndex,
                uint                  aCurrentBspCellIndex,
                Matrix< IndexMat >&   aSubPhaseIndices );

        // ----------------------------------------------------------------------------------

        /**
         * @brief convert list of subphase indices in B-spline to unordered searchable map
         *
         * @param aSubphaseIndicesInBsplineCell list of the subphase indices inside the given B-spline element
         * @param aSubphaseIndexToBsplineCellIndex output: unordered index to index map holding the subphase indices
         */
        void
        construct_subphase_in_bspline_cell_map(
                Matrix< IndexMat > const & aSubphaseIndicesInBsplineCell,
                IndexMap&                  aSubphaseIndexToBsplineCellIndex );

        // ----------------------------------------------------------------------------------

        /**
         * @brief create SP to SP connectivity graph local to the B-spline element
         *
         * @param aCutIntegrationMesh pointer to cut integration mesh
         * @param aSubphasesInBsplineCell list of the subphase indices inside the given B-spline element
         * @param aSubphaseIndicesToBspline searchable unordered index to index map holding the subphase indices
         * @param aPrunedBsplineSubphaseToSubphase output: SP to SP connectivity graph local to the B-spline element
         */
        void
        generate_pruned_subphase_graph_in_bspline_cell(
                Cut_Integration_Mesh*      aCutIntegrationMesh,
                Matrix< IndexMat > const & aSubphasesInBsplineCell,
                IndexMap&                  aSubphaseIndicesToBspline,
                Matrix< IndexMat >&        aPrunedBsplineSubphaseToSubphase );

        // ----------------------------------------------------------------------------------

        /**
         * @brief perform flood fill on the SP to SP connectivity graph local to the B-spline element
         *
         * @param aSubphasesInBspCell list of the subphase indices inside the given B-spline element
         * @param aPrunedSubPhaseToSubphase SP to SP connectivity graph local to the B-spline element
         * @param aSubPhaseBinEnrichmentVals list of indices of disconnected group of subphases assigned to each subphase
         * @param aMaxEnrichmentLevel output: largest index of a disconnected group of subphases
         */
        void
        find_subphase_bin_enrichment_levels_in_bspline_cell(
                Matrix< IndexMat > const & aSubphasesInBspCell,
                Matrix< IndexMat > const & aPrunedSubPhaseToSubphase,
                Matrix< IndexMat >&        aSubPhaseBinEnrichmentVals,
                moris_index&               aMaxEnrichmentLevel );

        // ----------------------------------------------------------------------------------

        /**
         * @brief sort SPs into groups using flood fill indices
         *
         * @param aSubphaseBin list of indices of disconnected group of subphases assigned to each subphase
         * @param aSubphaseIndicesInBsplineCell list of the subphase indices inside the given B-spline element
         * @param aNumSPGs number of subphase groups to be created
         * @return Vector< Vector< moris_index > > groups of SP indices that are connected
         */
        Vector< Vector< moris_index > >
        split_flood_fill_bin(
                Matrix< IndexMat >& aSubphaseBin,
                Matrix< IndexMat >& aSubphaseIndicesInBsplineCell,
                uint                aNumSPGs );

        // ----------------------------------------------------------------------------------

        /**
         * @brief find the side ordinals through which an SPG is connected to neighboring SPGs
         *
         * @param aCutIntegrationMesh pointer to cut integration mesh
         * @param aSPsInGroup list of SP indices in current SPG
         * @param aSubphaseIndicesToBspline searchable unordered index to index map holding the subphase indices
         * @return Vector< bool > binary "punch-card" of the side ordinals having an SPG connection (in 2D length of 4, in 3D a length of 6)
         */
        Vector< bool >
        collect_subphase_group_ligament_side_ordinals(
                Cut_Integration_Mesh*  aCutIntegrationMesh,
                Vector< moris_index >& aSPsInGroup,
                IndexMap&              aSubphaseIndicesToBspline );

        // ----------------------------------------------------------------------------------

        /**
         * @brief collects the indices of the IG cells in all SPs given to the function and
         * compiles them in one linear list
         *
         * @param aCutIntegrationMesh pointer to cut integration mesh (needed to access IG_cell_groups)
         * @param aSPsInSPG list of SPs in the SPG
         * @param aIgCellIndicesInSPG output/to fill: list of IG cell indices in SPG
         */
        void
        collect_ig_cell_indices_in_SPG(
                Cut_Integration_Mesh*         aCutIntegrationMesh,
                Vector< moris_index > const & aSPsInSPG,
                Vector< moris_index >&        aIgCellIndicesInSPG );

        // ----------------------------------------------------------------------------------

        /**
         * @brief constructs the SPG Neighborhood Connectivity graph after the SPGs themselves have been constructed
         *
         * @param aCutIntegrationMesh pointer to cut integration mesh
         * @param aBsplineMeshInfo pointer to the B-spline mesh info to be used for construction
         * @param aMeshIndex discretization mesh index
         * @return const std::shared_ptr< Subphase_Neighborhood_Connectivity > pointer to the SPG Neighborhood Connectivity
         */
        std::shared_ptr< Subphase_Neighborhood_Connectivity >
        construct_subphase_group_neighborhood(
                Cut_Integration_Mesh* aCutIntegrationMesh,
                Bspline_Mesh_Info*    aBsplineMeshInfo,
                const moris_index     aMeshIndex );

        // ----------------------------------------------------------------------------------

        /**
         * @brief setup the material connectivity of the subphase groups wrt. to the coarsest
         * B-spline element related to a given base IP cell. creates lists of material sub-
         * domains wrt which basis extensions will be needed in the enrichment
         *
         * @param aCutIntegrationMesh pointer to the cut integration mesh
         */
        void
        construct_SPG_material_connectivity_information( Cut_Integration_Mesh* aCutIntegrationMesh );

        // ----------------------------------------------------------------------------------

        /**
         * @brief
         *
         * @param aCutIntegrationMesh
         * @param aLagrangeElementIndex
         * @param aBsplineMeshListIndex
         * @param aMaterialSpgIndices
         * @param aVoidSpgIndices
         */
        void
        find_material_SPGs_on_Lagrange_Cell(
                Cut_Integration_Mesh*  aCutIntegrationMesh,
                const moris_index      aLagrangeElementIndex,
                const moris_index      aBsplineMeshListIndex,
                Vector< moris_index >& aMaterialSpgIndices,
                Vector< moris_index >& aVoidSpgIndices );

        // ----------------------------------------------------------------------------------

    };    // class Integration_Mesh_Generator

    // ----------------------------------------------------------------------------------

}    // namespace moris::xtk

#endif
