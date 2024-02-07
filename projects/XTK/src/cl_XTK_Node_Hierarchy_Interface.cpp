/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Node_Hierarchy_Interface.cpp
 *
 */

#include "cl_XTK_Node_Hierarchy_Interface.hpp"
#include "cl_XTK_Decomposition_Algorithm.hpp"
#include "cl_XTK_Integration_Mesh_Generator.hpp"
#include "cl_XTK_Cut_Integration_Mesh.hpp"
#include "cl_MTK_Cell_Info_Factory.hpp"
#include "cl_MTK_Cell_Info.hpp"
#include "fn_Pairing.hpp"
#include "cl_GEN_Geometry_Engine.hpp"
#include <algorithm>    // std::sort, std::stable_sort
#include <numeric>
#include "cl_Tracer.hpp"
#include <chrono>
namespace xtk
{

    // ----------------------------------------------------------------------------------

    enum Decomposition_Algorithm_Type
    Node_Hierarchy_Interface::get_algorithm_type() const
    {
        return Decomposition_Algorithm_Type::NODE_HEIRARCHY;    // NODE_HIERARCHY
    }

    // ----------------------------------------------------------------------------------

    moris_index
    Node_Hierarchy_Interface::get_signature() const
    {
        return 102;
    }

    // ----------------------------------------------------------------------------------

    bool
    Node_Hierarchy_Interface::has_geometric_independent_vertices() const
    {
        return false;
    }

    // ----------------------------------------------------------------------------------

    void
    Node_Hierarchy_Interface::perform(
            Integration_Mesh_Generation_Data* aMeshGenerationData,
            Decomposition_Data*               aDecompositionData,
            Cut_Integration_Mesh*             aCutIntegrationMesh,
            mtk::Mesh*                        aBackgroundMesh,
            Integration_Mesh_Generator*       aMeshGenerator )
    {
        Tracer tTracer( "XTK", "Decomposition_Algorithm", "Node_Hierarchy_Interface", aMeshGenerator->verbosity_level(), 0 );

        // keep track of some useful classes (avoid passing to every function)
        mGeometryEngine     = aMeshGenerator->get_geom_engine();
        mMeshGenerationData = aMeshGenerationData;
        mDecompositionData  = aDecompositionData;
        mCutIntegrationMesh = aCutIntegrationMesh;
        mBackgroundMesh     = aBackgroundMesh;
        mGenerator          = aMeshGenerator;

        // iterate through geometries
        for ( moris::uint iGeom = 0; iGeom < mGeometryEngine->get_number_of_geometries(); iGeom++ )
        {
            // set a new decomposition data
            *aDecompositionData = Decomposition_Data();
            aDecompositionData->mDecompId = 10000 * iGeom + this->get_signature();

            // cell groups relevant to this geometry pick them out (i.e. the ones intersected by the current geometry)
            Vector< std::shared_ptr< IG_Cell_Group > > tIgCellGroups;
            this->select_ig_cell_groups( tIgCellGroups );

            // extract the cells from groups into a continuous list
            Vector< mtk::Cell* > tIgCellsInGroups;
            aMeshGenerator->extract_cells_from_cell_groups( tIgCellGroups, tIgCellsInGroups );

            // construct the edge connectivity for these ig cell groups
            std::shared_ptr< Edge_Based_Connectivity > tIgCellGroupEdgeConnectivity = std::make_shared< Edge_Based_Connectivity >();
            aMeshGenerator->create_edges_from_element_to_node( tIgCellsInGroups, tIgCellGroupEdgeConnectivity );

            // collect a representative background cell for each edges
            Vector< mtk::Cell* > tBackgroundCellForEdge;    // input: edge index || output: pointer to BG cell it belongs to
            aMeshGenerator->select_background_cell_for_edge( tIgCellGroupEdgeConnectivity, aCutIntegrationMesh, tBackgroundCellForEdge );

            // collect the vertex group related to the representative background cell
            Vector< std::shared_ptr< IG_Vertex_Group > > tVertexGroups;    // input: BG cell index || output: IG vertex group
            aMeshGenerator->collect_vertex_groups_for_background_cells( aMeshGenerationData, aCutIntegrationMesh, &tBackgroundCellForEdge, &tVertexGroups );

            // deduce the edge parent entity index and rank
            std::shared_ptr< Edge_Based_Ancestry > tIgEdgeAncestry = std::make_shared< Edge_Based_Ancestry >();
            aMeshGenerator->deduce_edge_ancestry( aCutIntegrationMesh, aBackgroundMesh, tIgCellGroupEdgeConnectivity, tBackgroundCellForEdge, tIgEdgeAncestry );

            // figure out which edges are intersected
            Vector< moris_index > tIntersectedEdgeIndices;
            Vector< real >        tIntersectedEdgeLocCoords;

            this->determine_intersected_edges_and_make_requests(
                    tIgCellGroupEdgeConnectivity,
                    tIgEdgeAncestry,
                    &tBackgroundCellForEdge,
                    &tVertexGroups,
                    tIntersectedEdgeIndices,
                    tIntersectedEdgeLocCoords );

            // give all these nodes IDs
            aMeshGenerator->assign_node_requests_identifiers( *aDecompositionData, aCutIntegrationMesh, aBackgroundMesh );

            // create associations between child meshes and the vertices we need this to commit the data to the integration mesh
            this->associate_new_vertices_with_cell_groups(
                    tIgCellGroupEdgeConnectivity,
                    tIgEdgeAncestry,
                    &tBackgroundCellForEdge,
                    &tVertexGroups,
                    &tIntersectedEdgeIndices,
                    &tIntersectedEdgeLocCoords );

            // update underlying IDs and owners of interpolation nodes in GE
            for ( uint iNodeRequest = 0; iNodeRequest < aDecompositionData->tNewNodeIndex.size(); iNodeRequest++ )
            {
                moris_index tNodeIndex = aDecompositionData->tNewNodeIndex( iNodeRequest );
                moris_id    tNodeId    = aDecompositionData->tNewNodeId( iNodeRequest );
                moris_index tNodeOwner = aDecompositionData->tNewNodeOwner( iNodeRequest );

                mGeometryEngine->update_intersection_node( tNodeIndex, tNodeId, tNodeOwner );
            }

            // commit vertices to the mesh
            aMeshGenerator->commit_new_ig_vertices_to_cut_mesh( aMeshGenerationData, aDecompositionData, aCutIntegrationMesh, aBackgroundMesh, this );

            // we are ready to create the new integration cells
            this->create_node_hierarchy_integration_cells( tIgCellGroupEdgeConnectivity, tIgEdgeAncestry, &tIntersectedEdgeIndices );

            // commit the cells to the mesh
            aMeshGenerator->commit_new_ig_cells_to_cut_mesh( aMeshGenerationData, aDecompositionData, aCutIntegrationMesh, aBackgroundMesh, this );

            // Advance geometry index
            if ( iGeom != mGeometryEngine->get_number_of_geometries() )
            {
                mGeometryEngine->advance_geometry_index();
            }
        }

        // trim data of CutIntegration mesh
        aCutIntegrationMesh->trim_data();

    }    // end function: Node_Hierarchy_Interface::perform()

    // ----------------------------------------------------------------------------------

    bool
    Node_Hierarchy_Interface::determine_intersected_edges_and_make_requests(
            std::shared_ptr< Edge_Based_Connectivity >         aEdgeConnectivity,
            std::shared_ptr< Edge_Based_Ancestry >             aIgEdgeAncestry,
            Vector< mtk::Cell* >*                         aBackgroundCellForEdge,
            Vector< std::shared_ptr< IG_Vertex_Group > >* aVertexGroups,
            Vector< moris_index >&                        aIntersectedEdges,
            Vector< real >&                               aEdgeLocalCoordinate )
    {
        Tracer tTracer( "XTK", "Decomposition_Algorithm", "Determine Intersected Edges", mGenerator->verbosity_level(), 1 );

        // get first unused index for nodes for numbering new nodes
        moris_index tNewNodeIndex = mCutIntegrationMesh->get_first_available_index( mtk::EntityRank::NODE );

        // get the number of edges to be treated
        uint tNumEdges = aEdgeConnectivity->mEdgeVertices.size();

        // initialize intersection information
        aIntersectedEdges.clear();
        aIntersectedEdges.reserve( tNumEdges );
        aEdgeLocalCoordinate.clear();
        aEdgeLocalCoordinate.reserve( tNumEdges );

        // the query interface
        // Initialize geometric query
        // (object which holds/collects all information needed to determine an intersection
        // (e.g. data of an edge and its relation to its Child Mesh and BG Cell))
        Geometric_Query tGeometricQuery;

        // setup the query data (fixed parts for this function)
        tGeometricQuery.set_coordinates_matrix( mCutIntegrationMesh->get_all_vertex_coordinates_loc_inds() );
        tGeometricQuery.set_cut_integration_mesh( mCutIntegrationMesh );
        tGeometricQuery.set_query_entity_rank( mtk::EntityRank::EDGE );
        tGeometricQuery.set_edge_connectivity( aEdgeConnectivity );
        tGeometricQuery.set_edge_associated_background_cell( aBackgroundCellForEdge );
        tGeometricQuery.set_associated_vertex_group( aVertexGroups );

        mDecompositionData->mHasSecondaryIdentifier = true;

        // iterate through the edges in aEdgeConnectivity ask the geometry engine if we are intersected
        for ( uint iEdge = 0; iEdge < tNumEdges; iEdge++ )
        {
            // initialize matrix storing intersection point on edge in parametric coords of BG cell
            if ( iEdge == 0 )
            {
                tGeometricQuery.set_parametric_coordinate_size( ( *aBackgroundCellForEdge )( iEdge )->get_cell_info()->get_loc_coord_dim() );
            }

            // update the current edge index treated
            tGeometricQuery.set_current_edge_index( iEdge );

            // update parent cell to BG cell of treated edge
            tGeometricQuery.set_parent_cell( ( *aBackgroundCellForEdge )( iEdge ) );

            // Get query info
            const Matrix< IndexMat >& tEdgeToVertex = tGeometricQuery.get_query_entity_to_vertex_connectivity();
            Matrix< IndexMat > tParentEntityIndices = tGeometricQuery.get_query_parent_entity_connectivity();

            // annoying copy until this is converted to using a vector
            Matrix< DDUMat >         tParentEntityIndiceUINT( tParentEntityIndices.numel() );
            for ( moris::uint i = 0; i < tParentEntityIndices.numel(); i++ )
            {
                tParentEntityIndiceUINT( i ) = (uint)tParentEntityIndices( i );
            }

            // see if the edge is intersected using the geometry engine
            bool tIsIntersected = mGeometryEngine->queue_intersection( tEdgeToVertex( 0 ),
                    tEdgeToVertex( 1 ),
                    tGeometricQuery.get_vertex_local_coord_wrt_parent_entity( tEdgeToVertex( 0 ) ),
                    tGeometricQuery.get_vertex_local_coord_wrt_parent_entity( tEdgeToVertex( 1 ) ),
                    tParentEntityIndiceUINT,
                    tGeometricQuery.get_geometry_type(),
                    tGeometricQuery.get_interpolation_order() );

            // for intersected edges, invoke intersection procedure
            if ( tIsIntersected )
            {
                bool tBothVerticesNotOnInterface =
                        !mGeometryEngine->queued_intersection_first_parent_on_interface()
                        && !mGeometryEngine->queued_intersection_second_parent_on_interface();

                // if one of the end vertices is on the interface, skip this general intersection procedure and use specific one instead
                if ( tBothVerticesNotOnInterface )
                {
                    // add index and intersection position to list of intersected edges
                    aIntersectedEdges.push_back( (moris_index)iEdge );
                    aEdgeLocalCoordinate.push_back( mGeometryEngine->get_queued_intersection_local_coordinate() );

                    // get edge parent entity index and rank
                    moris_index tParentIndex = aIgEdgeAncestry->mEdgeParentEntityIndex( iEdge );
                    moris_index tParentRank  = aIgEdgeAncestry->mEdgeParentEntityRank( iEdge );

                    // get unique edge id based on two end vertices of edge
                    moris_index tSecondaryId = this->hash_edge( aEdgeConnectivity->mEdgeVertices( iEdge ) );

                    // initialize variable holding possible new node index
                    moris_index tNewNodeIndexInSubdivision = MORIS_INDEX_MAX;

                    // check if new node for current edge has already been requested ...
                    bool tRequestExist = mDecompositionData->request_exists(
                            tParentIndex,
                            tSecondaryId,
                            (mtk::EntityRank)tParentRank,
                            tNewNodeIndexInSubdivision );

                    // ... if not request it
                    if ( !tRequestExist )
                    {
                        // find out which processor owns parent entity of currently treated edge
                        moris_index tOwningProc = mBackgroundMesh->get_entity_owner( tParentIndex, (mtk::EntityRank)tParentRank );

                        // Register new node request
                        tNewNodeIndexInSubdivision = mDecompositionData->register_new_request(
                                tParentIndex,
                                tSecondaryId,
                                tOwningProc,
                                (mtk::EntityRank)tParentRank,
                                mGeometryEngine->get_queued_intersection_global_coordinates() );

                        // create new node in GEN
                        mGeometryEngine->admit_queued_intersection();

                        // count number of new nodes created
                        tNewNodeIndex++;
                    }
                }
            }
        }

        // return success if finished
        return true;

    }    // end function: Node_Hierarchy_Interface::determine_intersected_edges_and_make_requests()

    // ----------------------------------------------------------------------------------

    bool
    Node_Hierarchy_Interface::associate_new_vertices_with_cell_groups(
            std::shared_ptr< Edge_Based_Connectivity >         aEdgeConnectivity,
            std::shared_ptr< Edge_Based_Ancestry >             aIgEdgeAncestry,
            Vector< mtk::Cell* >*                         aBackgroundCellForEdge,
            Vector< std::shared_ptr< IG_Vertex_Group > >* aVertexGroups,
            Vector< moris_index >*                        aIntersectedEdges,
            Vector< real >*                               aEdgeLocalCoordinate )
    {
        // trace this function
        Tracer tTracer( "XTK", "Decomposition_Algorithm", "Vertex Associations", mGenerator->verbosity_level(), 1 );

        // number of child meshes
        moris_index tNumChildMeshes = mCutIntegrationMesh->get_num_ig_cell_groups();

        // estimate required space for the data in decomposition data
        mDecompositionData->tCMNewNodeLoc        = Vector< Vector< moris_index > >( tNumChildMeshes );
        mDecompositionData->tCMNewNodeParamCoord = Vector< Vector< Matrix< DDRMat > > >( tNumChildMeshes );

        // get dimensionality of parent coords in BG cell
        moris_index tDimParamCoords = 0;
        if ( aBackgroundCellForEdge->size() > 0 )
        {
            tDimParamCoords = ( *aBackgroundCellForEdge )( 0 )->get_cell_info()->get_loc_coord_dim();
        }

        // list of node coords
        Matrix< DDRMat > tEdgeNodeParamCoordinates( 2, tDimParamCoords );

        Vector< uint > tCellGroups;
        tCellGroups.reserve( 20 );

        // iterate through the edges
        for ( uint iEdge = 0; iEdge < aIntersectedEdges->size(); iEdge++ )
        {
            moris_index tEdgeIndex = ( *aIntersectedEdges )( iEdge );

            // iterate through elements attached to this edge and collect the cell groups
            tCellGroups.clear();
            for ( uint iCell = 0; iCell < aEdgeConnectivity->mEdgeToCell( tEdgeIndex ).size(); iCell++ )
            {
                // integration cell just grab the first
                mtk::Cell* tCell = aEdgeConnectivity->mEdgeToCell( tEdgeIndex )( iCell );

                // fill list of Cell groups the edge is related to
                tCellGroups.push_back( mCutIntegrationMesh->get_ig_cell_group_memberships( tCell->get_index() )( 0 ) );
            }

            // remove duplicates from list
            moris::unique( tCellGroups );

            // iterate through the unique cell groups
            for ( uint iUCG = 0; iUCG < tCellGroups.size(); iUCG++ )
            {
                // get the background cell group associated with this cell
                MORIS_ERROR( aEdgeConnectivity->mEdgeToCell( tEdgeIndex ).size() > 0, "Edge not connected to any cells..." );

                // cell group membership
                moris_index tCellGroupMembershipIndex = tCellGroups( iUCG );

                // get the vertex group
                std::shared_ptr< IG_Vertex_Group > tVertGroup = mCutIntegrationMesh->get_vertex_group( tCellGroupMembershipIndex );

                // get the parametric coord wrt this parent cell
                std::shared_ptr< Matrix< DDRMat > > tVertex0LocalCoords = tVertGroup->get_vertex_local_coords( aEdgeConnectivity->mEdgeVertices( tEdgeIndex )( 0 )->get_index() );
                std::shared_ptr< Matrix< DDRMat > > tVertex1LocalCoords = tVertGroup->get_vertex_local_coords( aEdgeConnectivity->mEdgeVertices( tEdgeIndex )( 1 )->get_index() );

                // initialize variables to store coords of new vertex
                tEdgeNodeParamCoordinates.set_row( 0, *tVertex0LocalCoords );
                tEdgeNodeParamCoordinates.set_row( 1, *tVertex1LocalCoords );

                // get parametric coords wrt BG element where new vertex sits
                Matrix< DDRMat > tParametricCoordsRelativeToParentElem =
                        Interpolation::linear_interpolation_location( tEdgeNodeParamCoordinates, { { ( *aEdgeLocalCoordinate )( iEdge ) } } );

                // list of edges and locations to create new nodes on added to decomp. data
                mDecompositionData->tCMNewNodeLoc( tCellGroupMembershipIndex ).push_back( iEdge );
                mDecompositionData->tCMNewNodeParamCoord( tCellGroupMembershipIndex ).push_back( tParametricCoordsRelativeToParentElem );
            }
        }

        // return success if function is done
        return true;
    }

    moris_index
    Node_Hierarchy_Interface::hash_edge( Vector< mtk::Vertex* > const & aEdgeVertices )
    {
        MORIS_ERROR( aEdgeVertices.size() == 2, "Edge is expected to have two vertices" );
        moris_index tMinIdIndex = 0;
        moris_index tMaxIdIndex = 1;
        if ( aEdgeVertices( 1 )->get_id() < aEdgeVertices( 0 )->get_id() )
        {
            tMinIdIndex = 1;
            tMaxIdIndex = 0;
        }
        return xtk::cantor_pairing( aEdgeVertices( tMinIdIndex )->get_id(), aEdgeVertices( tMaxIdIndex )->get_id() );
    }

    void
    Node_Hierarchy_Interface::perform_impl_vertex_requests(
            Integration_Mesh_Generation_Data* aMeshGenerationData,
            Decomposition_Data*               aDecompositionData,
            Cut_Integration_Mesh*             aCutIntegrationMesh,
            mtk::Mesh*                        aBackgroundMesh,
            Integration_Mesh_Generator*       aMeshGenerator )
    {
    }

    void
    Node_Hierarchy_Interface::perform_impl_generate_mesh(
            Integration_Mesh_Generation_Data* aMeshGenerationData,
            Decomposition_Data*               aDecompositionData,
            Cut_Integration_Mesh*             aCutIntegrationMesh,
            mtk::Mesh*                        aBackgroundMesh,
            Integration_Mesh_Generator*       aMeshGenerator )
    {
    }

    void
    Node_Hierarchy_Interface::select_ig_cell_groups(
        Vector< std::shared_ptr< IG_Cell_Group > >& aIgCellGroups )
    {
        Tracer tTracer( "XTK", "Decomposition_Algorithm", "Select Ig Cell Groups", mGenerator->verbosity_level(), 1 );
        // intersected background cells for the current geometry
        Vector< moris_index > const& tIntersectedBackground = mMeshGenerationData->mIntersectedBackgroundCellIndex(
                mGeometryEngine->get_active_geometry_index() );

        // number of intersected background cells
        uint tNumIntersectedBackground = tIntersectedBackground.size();

        // allocate input data
        aIgCellGroups = Vector< std::shared_ptr< IG_Cell_Group > >( tNumIntersectedBackground, nullptr );

        // iterate through them and get the cell group index
        for ( uint iBGCell = 0; iBGCell < tNumIntersectedBackground; iBGCell++ )
        {
            moris_index tCurrentBgCellIndex = tIntersectedBackground( iBGCell );

            aIgCellGroups( iBGCell ) = mCutIntegrationMesh->get_ig_cell_group( tCurrentBgCellIndex );
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Node_Hierarchy_Interface::create_node_hierarchy_integration_cells(
            std::shared_ptr< Edge_Based_Connectivity > aEdgeConnectivity,
            std::shared_ptr< Edge_Based_Ancestry >     aIgEdgeAncestry,
            Vector< moris_index >*                aIntersectedEdges )
    {
        // time/log function
        Tracer tTracer( "XTK", "Node_Hierarchy_Interface", "Create NH IG Cells", mGenerator->verbosity_level(), 1 );

        // starting with a clean slate here
        mNumNewCells = 0;

        // get the number of elements currently in the IG mesh before applying the templates
        uint tNumElemsCurrentlyInCutIgMesh = mCutIntegrationMesh->get_num_entities( mtk::EntityRank::ELEMENT, 0 );

        // initialize map, input: IG Cell index || output: list of intersected edge ordinals (?)
        Vector< std::shared_ptr< Vector< moris_index > > >
                tCellIndexIntersectedEdgeOrdinals( tNumElemsCurrentlyInCutIgMesh, nullptr );

        // initialize map, input: IG Cell index || output: list of intersected vertex ordinals (?)
        Vector< std::shared_ptr< Vector< mtk::Vertex* > > >
                tCellIndexIntersectedEdgeVertex( tNumElemsCurrentlyInCutIgMesh, nullptr );

        // fill necessary information in maps initialized above
        this->determine_intersected_cell_information(
                aEdgeConnectivity,
                aIntersectedEdges,
                &tCellIndexIntersectedEdgeOrdinals,
                &tCellIndexIntersectedEdgeVertex );

        Vector< std::shared_ptr< Vector< mtk::Vertex* > > > tNodesForTemplates( tNumElemsCurrentlyInCutIgMesh, nullptr );
        Vector< std::shared_ptr< Node_Hierarchy_Template > >     tNHTemplate( tNumElemsCurrentlyInCutIgMesh, nullptr );

        // for each cell we need to select a node hierarchical template
        // number of cells we are creating
        moris_index                       tNumNewCells  = 0;
        moris_index                       tNodesPerCell = 0;
        std::shared_ptr< mtk::Cell_Info > tCellInfo     = nullptr;

        // Case: 3D mesh
        if ( mBackgroundMesh->get_spatial_dim() == 3 )
        {
            tNumNewCells = this->select_node_hier_3d_template(
                    &tCellIndexIntersectedEdgeOrdinals,
                    &tCellIndexIntersectedEdgeVertex,
                    &tNodesForTemplates,
                    &tNHTemplate );

            // constant parameters for 3d case
            tNodesPerCell = 4;
            mtk::Cell_Info_Factory tFactory;
            tCellInfo = tFactory.create_cell_info_sp( mtk::CellTopology::TET4 );
        }

        // Case: 2D mesh
        else if ( mBackgroundMesh->get_spatial_dim() == 2 )
        {
            tNumNewCells = this->select_node_hier_2d_template(
                    &tCellIndexIntersectedEdgeOrdinals,
                    &tCellIndexIntersectedEdgeVertex,
                    &tNodesForTemplates,
                    &tNHTemplate );

            // constant parameters for 2d case
            tNodesPerCell = 3;
            mtk::Cell_Info_Factory tFactory;
            tCellInfo = tFactory.create_cell_info_sp( mtk::CellTopology::TRI3 );
        }

        // Error for other number of spatial dims
        else
        {
            MORIS_ERROR( false, "Node_Hierarchy_Interface::create_node_hierarchy_integration_cells() - Only implemented in 2D and 3D." );
        }

        // allocate the algorithm data
        mNewCellToVertexConnectivity = Vector< Vector< moris_index > >( tNumNewCells, Vector< moris_index >( tNodesPerCell ) );
        mNewCellChildMeshIndex       = Vector< moris_index >( tNumNewCells );
        mNewCellCellIndexToReplace   = Vector< moris_index >( tNumNewCells );
        mNewCellCellInfo             = Vector< std::shared_ptr< mtk::Cell_Info > >( tNumNewCells, tCellInfo );

        moris_index tCurrentCellIndex = 0;
        for ( uint iCell = 0; iCell < tNumElemsCurrentlyInCutIgMesh; iCell++ )
        {
            if ( tNodesForTemplates( iCell ) != nullptr )
            {
                // cell group membership
                moris_index tCellGroupMembershipIndex = mCutIntegrationMesh->get_ig_cell_group_memberships( (moris_index)iCell )( 0 );

                // get the template
                std::shared_ptr< Node_Hierarchy_Template > tTemplate = tNHTemplate( iCell );

                // vertices sorted for this template
                Vector< mtk::Vertex* >& tSortedVertices = *(tNodesForTemplates)( iCell );

                // Matrix< IndexMat > mCellToNodeOrdinal;
                // num cells in template
                for ( moris_index iTemplateCell = 0; iTemplateCell < tTemplate->mNumCells; iTemplateCell++ )
                {
                    mNewCellChildMeshIndex( tCurrentCellIndex ) = tCellGroupMembershipIndex;

                    if ( iTemplateCell == 0 )
                    {
                        mNewCellCellIndexToReplace( tCurrentCellIndex ) = iCell;
                    }
                    else
                    {
                        mNewCellCellIndexToReplace( tCurrentCellIndex ) = MORIS_INDEX_MAX;
                    }

                    // iterate through vertices and put them in the algorithm data
                    for ( uint iV = 0; iV < tTemplate->mCellToNodeOrdinal.n_cols(); iV++ )
                    {
                        //
                        moris_index tVertexSortedOrd                            = tTemplate->mCellToNodeOrdinal( iTemplateCell, iV );
                        mNewCellToVertexConnectivity( tCurrentCellIndex )( iV ) = tSortedVertices( tVertexSortedOrd )->get_index();
                    }
                    tCurrentCellIndex++;
                }
            }

        }    // end for: each IG cell currently in the cut mesh

    }        // end function: Node_Hierarchy_Interface::create_node_hierarchy_integration_cells()

    // ----------------------------------------------------------------------------------

    void
    Node_Hierarchy_Interface::determine_intersected_cell_information(
            std::shared_ptr< Edge_Based_Connectivity >                     aEdgeConnectivity,
            Vector< moris_index >*                                    aIntersectedEdges,
            Vector< std::shared_ptr< Vector< moris_index > > >*  aCellIndexIntersectedEdgeOrdinals,
            Vector< std::shared_ptr< Vector< mtk::Vertex* > > >* aCellIndexIntersectedEdgeVertex )
    {
        // max new cell per intersected cell ~6
        uint tNumEdgesPerCell = 6;    // just use 6's since we operate on tets or TRIs. more conservative with 6 (slightly more memory)

        // iterate through intersected edges and figure out which side ordinals of an ig cell is intersected
        for ( uint iEdge = 0; iEdge < aIntersectedEdges->size(); iEdge++ )
        {
            // get index of currently treated edge
            moris_index tEdgeIndex = ( *aIntersectedEdges )( iEdge );

            // new vertex on the edge
            mtk::Vertex* tVertexOnEdge = mCutIntegrationMesh->get_mtk_vertex_pointer( mDecompositionData->tNewNodeIndex( iEdge ) );

            // go through cells attached to current edge
            for ( uint iCell = 0; iCell < aEdgeConnectivity->mEdgeToCell( tEdgeIndex ).size(); iCell++ )
            {
                // integration cell with this edge
                mtk::Cell* tCell = aEdgeConnectivity->mEdgeToCell( tEdgeIndex )( iCell );

                // integration cell edge ordinal
                moris_index tSideOrdinal = aEdgeConnectivity->mEdgeToCellEdgeOrdinal( tEdgeIndex )( iCell );

                moris_index tCellIndex = tCell->get_index();

                if ( ( *aCellIndexIntersectedEdgeOrdinals )( tCellIndex ) == nullptr )
                {
                    ( *aCellIndexIntersectedEdgeOrdinals )( tCellIndex ) = std::make_shared< Vector< moris_index > >( 0 );
                    ( *aCellIndexIntersectedEdgeOrdinals )( tCellIndex )->reserve( tNumEdgesPerCell );
                    ( *aCellIndexIntersectedEdgeVertex )( tCellIndex ) = std::make_shared< Vector< mtk::Vertex* > >( 0 );
                    ( *aCellIndexIntersectedEdgeVertex )( tCellIndex )->reserve( tNumEdgesPerCell );
                }

                ( *aCellIndexIntersectedEdgeOrdinals )( tCellIndex )->push_back( tSideOrdinal );
                ( *aCellIndexIntersectedEdgeVertex )( tCellIndex )->push_back( tVertexOnEdge );
            }
        }
    }

    // ----------------------------------------------------------------------------------

    moris_index
    Node_Hierarchy_Interface::select_node_hier_2d_template(
            Vector< std::shared_ptr< Vector< moris_index > > >*  aCellIndexIntersectedEdgeOrdinals,
            Vector< std::shared_ptr< Vector< mtk::Vertex* > > >* aCellIndexIntersectedEdgeVertex,
            Vector< std::shared_ptr< Vector< mtk::Vertex* > > >* aNodesForTemplates,
            Vector< std::shared_ptr< Node_Hierarchy_Template > >*     aNHTemplate )
    {
        // tally up the total number of cells we are going to construct
        moris_index tNumNewIgCells = 0;

        mNumNewCells = 0;

        std::unordered_map< moris_index, std::shared_ptr< Node_Hierarchy_Template > > tLoadedTemplates;
        std::shared_ptr< Matrix< IdMat > >                                            tEdgeToVertexOrdinalMap = nullptr;

        Node_Hierarchy_Template_Library tLibrary;

        // iterate through cells
        for ( uint iCell = 0; iCell < aCellIndexIntersectedEdgeOrdinals->size(); iCell++ )
        {
            if ( ( *aCellIndexIntersectedEdgeOrdinals )( iCell ) != nullptr )
            {
                mtk::Cell const * tIgCell = &mCutIntegrationMesh->get_mtk_cell( iCell );

                if ( tEdgeToVertexOrdinalMap == nullptr )
                {
                    tEdgeToVertexOrdinalMap = std::make_shared< Matrix< IdMat > >( tIgCell->get_cell_info()->get_node_to_edge_map() );
                }
                moris_index tPermutationId = MORIS_INDEX_MAX;

                ( *aNodesForTemplates )( iCell ) = std::make_shared< Vector< mtk::Vertex* > >();
                ( *aNodesForTemplates )( iCell )->reserve( 5 );

                this->sort_nodes_2d(
                        tIgCell,
                        tEdgeToVertexOrdinalMap.get(),
                        ( *aCellIndexIntersectedEdgeOrdinals )( iCell ),
                        ( *aCellIndexIntersectedEdgeVertex )( iCell ),
                        tPermutationId,
                        ( *aNodesForTemplates )( iCell ) );

                // if we haven't used this template yet, load it up
                if ( tLoadedTemplates.find( tPermutationId ) == tLoadedTemplates.end() )
                {
                    std::shared_ptr< Node_Hierarchy_Template > tNewTemplate = std::make_shared< Node_Hierarchy_Template >();
                    tLibrary.load_template( 2, tPermutationId, tNewTemplate.get() );

                    tLoadedTemplates[ tPermutationId ] = tNewTemplate;
                }
                ( *aNHTemplate )( iCell ) = tLoadedTemplates[ tPermutationId ];

                tNumNewIgCells = tNumNewIgCells + ( *aNHTemplate )( iCell )->mNumCells;
                mNumNewCells   = mNumNewCells + ( *aNHTemplate )( iCell )->mNumCells - 1;
            }
        }
        return tNumNewIgCells;
    }

    // ----------------------------------------------------------------------------------

    moris_index
    Node_Hierarchy_Interface::select_node_hier_3d_template(
            Vector< std::shared_ptr< Vector< moris_index > > >*  aCellIndexIntersectedEdgeOrdinals,
            Vector< std::shared_ptr< Vector< mtk::Vertex* > > >* aCellIndexIntersectedEdgeVertex,
            Vector< std::shared_ptr< Vector< mtk::Vertex* > > >* aNodesForTemplates,
            Vector< std::shared_ptr< Node_Hierarchy_Template > >*     aNHTemplate )
    {
        // tally up the total number of cells we are going to construct
        moris_index tNumNewIgCells = 0;

        mNumNewCells = 0;

        std::unordered_map< moris_index, std::shared_ptr< Node_Hierarchy_Template > > tLoadedTemplates;
        std::shared_ptr< Matrix< IdMat > >                                            tEdgeToVertexOrdinalMap = nullptr;

        Node_Hierarchy_Template_Library tLibrary;

        // iterate through cells
        for ( uint iCell = 0; iCell < aCellIndexIntersectedEdgeOrdinals->size(); iCell++ )
        {
            if ( ( *aCellIndexIntersectedEdgeOrdinals )( iCell ) != nullptr )
            {
                mtk::Cell const * tIgCell = &mCutIntegrationMesh->get_mtk_cell( iCell );

                if ( tEdgeToVertexOrdinalMap == nullptr )
                {
                    tEdgeToVertexOrdinalMap = std::make_shared< Matrix< IdMat > >( tIgCell->get_cell_info()->get_node_to_edge_map() );
                }
                moris_index tPermutationId = MORIS_INDEX_MAX;

                ( *aNodesForTemplates )( iCell ) = std::make_shared< Vector< mtk::Vertex* > >();
                ( *aNodesForTemplates )( iCell )->reserve( 8 );

                this->sort_nodes_3d(
                        tIgCell,
                        tEdgeToVertexOrdinalMap.get(),
                        ( *aCellIndexIntersectedEdgeOrdinals )( iCell ),
                        ( *aCellIndexIntersectedEdgeVertex )( iCell ),
                        tPermutationId,
                        ( *aNodesForTemplates )( iCell ) );

                // std::cout<<"tPermutationId = "<<tPermutationId<<std::endl;
                // moris::print(*(*aNodesForTemplates)(iCell),"(*aNodesForTemplates)(iCell)");

                // if we haven't used this template yet, load it up
                if ( tLoadedTemplates.find( tPermutationId ) == tLoadedTemplates.end() )
                {
                    std::shared_ptr< Node_Hierarchy_Template > tNewTemplate = std::make_shared< Node_Hierarchy_Template >();
                    tLibrary.load_template( 3, tPermutationId, tNewTemplate.get() );

                    tLoadedTemplates[ tPermutationId ] = tNewTemplate;
                }
                ( *aNHTemplate )( iCell ) = tLoadedTemplates[ tPermutationId ];

                tNumNewIgCells = tNumNewIgCells + ( *aNHTemplate )( iCell )->mNumCells;
                mNumNewCells   = mNumNewCells + ( *aNHTemplate )( iCell )->mNumCells - 1;
            }
        }
        return tNumNewIgCells;
    }

    // ----------------------------------------------------------------------------------

    void
    Node_Hierarchy_Interface::sort_nodes_2d(
            mtk::Cell const *                              aIgCell,
            Matrix< IndexMat >*                            aEdgeToVertexOrdinalMap,
            std::shared_ptr< Vector< moris_index > >  aCellIndexIntersectedEdgeOrdinals,
            std::shared_ptr< Vector< mtk::Vertex* > > aCellIndexIntersectedEdgeVertex,
            moris_index&                                   aPermutation,
            std::shared_ptr< Vector< mtk::Vertex* > > aSortedNodeInds )
    {
        // hier TRI3
        MORIS_ERROR( mBackgroundMesh->get_spatial_dim() == 2, "Node_Hierarchy_Interface::sort_nodes_2d() - number of spatial dimensions is not 2." );

        // initialize original index locations
        Vector< moris_index > tIndices( aCellIndexIntersectedEdgeOrdinals->size() );
        std::iota( tIndices.data().begin(), tIndices.data().end(), 0 );

        // get ascending order vertices
        std::stable_sort(
                tIndices.data().begin(),
                tIndices.data().end(),
                [ & ]( std::size_t i1, std::size_t i2 ) {
                    return ( *aCellIndexIntersectedEdgeOrdinals )( i1 ) < ( *aCellIndexIntersectedEdgeOrdinals )( i2 );
                } );

        Vector< mtk::Vertex* > tVertices = aIgCell->get_vertex_pointers();

        // intersection goes through two edges
        if ( aCellIndexIntersectedEdgeOrdinals->size() == 2 )
        {
            // get unique permutation id, identifying which template to use
            moris_index tLowVertexIdEdgeOrd  = ( *aCellIndexIntersectedEdgeOrdinals )( tIndices( 0 ) );
            moris_index tHighVertexIdEdgeOrd = ( *aCellIndexIntersectedEdgeOrdinals )( tIndices( 1 ) );
            aPermutation                     = tLowVertexIdEdgeOrd + tHighVertexIdEdgeOrd;

            aSortedNodeInds->resize( 5 );
            ( *aSortedNodeInds )( 0 ) = tVertices( 0 );
            ( *aSortedNodeInds )( 1 ) = tVertices( 1 );
            ( *aSortedNodeInds )( 2 ) = tVertices( 2 );
            ( *aSortedNodeInds )( 3 ) = ( *aCellIndexIntersectedEdgeVertex )( tIndices( 0 ) );
            ( *aSortedNodeInds )( 4 ) = ( *aCellIndexIntersectedEdgeVertex )( tIndices( 1 ) );
        }

        // intersection goes through one of the vertices
        else if ( aCellIndexIntersectedEdgeOrdinals->size() == 1 )
        {
            moris_index tVertexIdEdgeOrd = ( *aCellIndexIntersectedEdgeOrdinals )( tIndices( 0 ) );
            aPermutation                 = tVertexIdEdgeOrd + 10;

            aSortedNodeInds->resize( 4 );
            ( *aSortedNodeInds )( 0 ) = tVertices( 0 );
            ( *aSortedNodeInds )( 1 ) = tVertices( 1 );
            ( *aSortedNodeInds )( 2 ) = tVertices( 2 );
            ( *aSortedNodeInds )( 3 ) = ( *aCellIndexIntersectedEdgeVertex )( tIndices( 0 ) );
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Node_Hierarchy_Interface::sort_nodes_3d(
            mtk::Cell const *                              aIgCell,
            Matrix< IndexMat >*                            aEdgeToVertexOrdinalMap,
            std::shared_ptr< Vector< moris_index > >  aCellIndexIntersectedEdgeOrdinals,
            std::shared_ptr< Vector< mtk::Vertex* > > aCellIndexIntersectedEdgeVertex,
            moris_index&                                   aPermutation,
            std::shared_ptr< Vector< mtk::Vertex* > > aSortedNodeInds )
    {
        // hier tet 4
        MORIS_ERROR( mBackgroundMesh->get_spatial_dim() == 3, "Node_Hierarchy_Interface::sort_nodes_3d() - number of spatial dimensions is not 3." );

        // initialize original index locations
        Vector< moris_index > tIndices( aCellIndexIntersectedEdgeOrdinals->size() );
        std::iota( tIndices.data().begin(), tIndices.data().end(), 0 );

        // get ascending order vertices
        std::stable_sort(
                tIndices.data().begin(),
                tIndices.data().end(),
                [ & ]( std::size_t i1, std::size_t i2 ) {
                    return ( *aCellIndexIntersectedEdgeVertex )( i1 )->get_id() < ( *aCellIndexIntersectedEdgeVertex )( i2 )->get_id();
                } );

        Vector< mtk::Vertex* > tVertices = aIgCell->get_vertex_pointers();

        if ( aCellIndexIntersectedEdgeOrdinals->size() == 3 )
        {
            moris_index tLowIntersectVertexIdEdgeOrd  = ( *aCellIndexIntersectedEdgeOrdinals )( tIndices( 0 ) );
            moris_index tMidIntersectVertexIdEdgeOrd  = ( *aCellIndexIntersectedEdgeOrdinals )( tIndices( 1 ) );
            moris_index tHighIntersectVertexIdEdgeOrd = ( *aCellIndexIntersectedEdgeOrdinals )( tIndices( 2 ) );

            // the node shared by the edge with the low vertex id and the mid vertex id
            moris_index tN0Ordinal = aIgCell->get_cell_info()->get_shared_vertex_ordinal_between_edges( tLowIntersectVertexIdEdgeOrd, tMidIntersectVertexIdEdgeOrd );

            // moris_index tN0 = tVertices(tN0NodeOrdinal)->get_index();

            // the node on the other side of the edge with the low vertex id
            moris_index tN1Ordinal = tN0Ordinal == ( *aEdgeToVertexOrdinalMap )( tLowIntersectVertexIdEdgeOrd, 0 ) ? ( *aEdgeToVertexOrdinalMap )( tLowIntersectVertexIdEdgeOrd, 1 ) : ( *aEdgeToVertexOrdinalMap )( tLowIntersectVertexIdEdgeOrd, 0 );

            moris_index tN2Ordinal = tN0Ordinal == ( *aEdgeToVertexOrdinalMap )( tMidIntersectVertexIdEdgeOrd, 0 ) ? ( *aEdgeToVertexOrdinalMap )( tMidIntersectVertexIdEdgeOrd, 1 ) : ( *aEdgeToVertexOrdinalMap )( tMidIntersectVertexIdEdgeOrd, 0 );

            moris_index tN3Ordinal = tN0Ordinal == ( *aEdgeToVertexOrdinalMap )( tHighIntersectVertexIdEdgeOrd, 0 ) ? ( *aEdgeToVertexOrdinalMap )( tHighIntersectVertexIdEdgeOrd, 1 ) : ( *aEdgeToVertexOrdinalMap )( tHighIntersectVertexIdEdgeOrd, 0 );

            // Determine permutation
            // Rule:  1   * edge ordinal containing the lowest node ID
            //      + 10  * edge ordinal containing the middle node ID
            //      + 100 * edge ordinal containing the highest node ID
            aPermutation = tLowIntersectVertexIdEdgeOrd + 10 * tMidIntersectVertexIdEdgeOrd + 100 * tHighIntersectVertexIdEdgeOrd;

            aSortedNodeInds->resize( 7 );
            ( *aSortedNodeInds )( 0 ) = tVertices( tN0Ordinal );
            ( *aSortedNodeInds )( 1 ) = tVertices( tN1Ordinal );
            ( *aSortedNodeInds )( 2 ) = tVertices( tN2Ordinal );
            ( *aSortedNodeInds )( 3 ) = tVertices( tN3Ordinal );
            ( *aSortedNodeInds )( 4 ) = ( *aCellIndexIntersectedEdgeVertex )( tIndices( 0 ) );
            ( *aSortedNodeInds )( 5 ) = ( *aCellIndexIntersectedEdgeVertex )( tIndices( 1 ) );
            ( *aSortedNodeInds )( 6 ) = ( *aCellIndexIntersectedEdgeVertex )( tIndices( 2 ) );
        }
        else if ( aCellIndexIntersectedEdgeOrdinals->size() == 4 )
        {
            // figure out which interface nodes are away from each other in this case
            moris_index tLowIntersectVertexIdEdgeOrd     = ( *aCellIndexIntersectedEdgeOrdinals )( tIndices( 0 ) );
            moris_index tMidLowIntersectVertexIdEdgeOrd  = ( *aCellIndexIntersectedEdgeOrdinals )( tIndices( 1 ) );
            moris_index tMidHighIntersectVertexIdEdgeOrd = ( *aCellIndexIntersectedEdgeOrdinals )( tIndices( 2 ) );
            moris_index tHighIntersectVertexIdEdgeOrd    = ( *aCellIndexIntersectedEdgeOrdinals )( tIndices( 3 ) );

            // std::cout<<"tLowIntersectVertexIdEdgeOrd = "<<tLowIntersectVertexIdEdgeOrd<<std::endl;
            // std::cout<<"tMidLowIntersectVertexIdEdgeOrd = "<<tMidLowIntersectVertexIdEdgeOrd<<std::endl;
            // std::cout<<"tMidHighIntersectVertexIdEdgeOrd = "<<tMidHighIntersectVertexIdEdgeOrd<<std::endl;
            // std::cout<<"tHighIntersectVertexIdEdgeOrd = "<<tHighIntersectVertexIdEdgeOrd<<std::endl;

            // Determine permutation
            // Rule:  1    * edge ordinal containing the lowest node ID
            //      + 10   * edge ordinal containing the middle lowest node ID
            //      + 100  * edge ordinal containing the middle highest node ID
            //      + 1000 * edge ordinal containing the highest node ID
            aPermutation = tLowIntersectVertexIdEdgeOrd + 10 * tMidLowIntersectVertexIdEdgeOrd + 100 * tMidHighIntersectVertexIdEdgeOrd + 1000 * tHighIntersectVertexIdEdgeOrd;

            // determine which intersection node ids are across from eachother (meaning they share only an edge through the starting integration cell)
            moris_index tHLOppVertOrd   = aIgCell->get_cell_info()->get_shared_vertex_ordinal_between_edges( tLowIntersectVertexIdEdgeOrd, tHighIntersectVertexIdEdgeOrd );
            moris_index tHMHOppVertOrd  = aIgCell->get_cell_info()->get_shared_vertex_ordinal_between_edges( tMidHighIntersectVertexIdEdgeOrd, tHighIntersectVertexIdEdgeOrd );
            moris_index tHMLOppVertOrd  = aIgCell->get_cell_info()->get_shared_vertex_ordinal_between_edges( tMidLowIntersectVertexIdEdgeOrd, tHighIntersectVertexIdEdgeOrd );
            moris_index tLMLOppVertOrd  = aIgCell->get_cell_info()->get_shared_vertex_ordinal_between_edges( tLowIntersectVertexIdEdgeOrd, tMidLowIntersectVertexIdEdgeOrd );
            moris_index tLMHOppVertOrd  = aIgCell->get_cell_info()->get_shared_vertex_ordinal_between_edges( tLowIntersectVertexIdEdgeOrd, tMidHighIntersectVertexIdEdgeOrd );
            moris_index tMLMHOppVertOrd = aIgCell->get_cell_info()->get_shared_vertex_ordinal_between_edges( tMidLowIntersectVertexIdEdgeOrd, tMidHighIntersectVertexIdEdgeOrd );

            // std::cout<< "tHLOppVertOrd  = "<<tHLOppVertOrd <<std::endl;
            // std::cout<< "tHMHOppVertOrd = "<<tHMHOppVertOrd<<std::endl;
            // std::cout<< "tHMLOppVertOrd = "<<tHMLOppVertOrd<<std::endl;
            // std::cout<< "tLMLOppVertOrd = "<<tLMLOppVertOrd<<std::endl;
            // std::cout<< "tLMHOppVertOrd = "<<tLMHOppVertOrd<<std::endl;
            // std::cout<< "tMLMHOppVertOrd= "<<tMLMHOppVertOrd<<std::endl;

            moris_index tN0Ordinal = MORIS_INDEX_MAX;
            moris_index tN1Ordinal = MORIS_INDEX_MAX;
            moris_index tN2Ordinal = MORIS_INDEX_MAX;
            moris_index tN3Ordinal = MORIS_INDEX_MAX;
            if ( tHLOppVertOrd == MORIS_INDEX_MAX )
            {
                // orient it so that the first vertex is shared by the H and MH edges
                // N0 - shared by the H and MH edges
                // N1 - shared by the H and ML edges
                // N2 - shared by the L and MH edges
                // N3 - shared by the L and ML edges
                // shared vertex ordinal between H and MH edge
                tN0Ordinal = tHMHOppVertOrd;
                tN1Ordinal = tHMLOppVertOrd;
                tN2Ordinal = tLMHOppVertOrd;
                tN3Ordinal = tLMLOppVertOrd;
            }

            else if ( tHMHOppVertOrd == MORIS_INDEX_MAX )
            {
                // orient it so that the first vertex is shared by the H and MH edges
                // N0 - shared by the L  and H edges
                // N1 - shared by the ML and H edges
                // N2 - shared by the L  and MH edges
                // N3 - shared by the ML  and MH edges
                // shared vertex ordinal between H and MH edge
                tN0Ordinal = tHLOppVertOrd;
                tN1Ordinal = tHMLOppVertOrd;
                tN2Ordinal = tLMHOppVertOrd;
                tN3Ordinal = tMLMHOppVertOrd;
            }
            else if ( tHMLOppVertOrd == MORIS_INDEX_MAX )
            {
                // orient it so that the first vertex is shared by the H and MH edges
                // N0 - shared by the L  and H edges
                // N1 - shared by the MH and H edges
                // N2 - shared by the L  and ML edges
                // N3 - shared by the ML  and MH edges
                // shared vertex ordinal between H and MH edge
                tN0Ordinal = tHLOppVertOrd;
                tN1Ordinal = tHMHOppVertOrd;
                tN2Ordinal = tLMLOppVertOrd;
                tN3Ordinal = tMLMHOppVertOrd;
            }

            aSortedNodeInds->resize( 8 );
            ( *aSortedNodeInds )( 0 ) = tVertices( tN0Ordinal );
            ( *aSortedNodeInds )( 1 ) = tVertices( tN1Ordinal );
            ( *aSortedNodeInds )( 2 ) = tVertices( tN2Ordinal );
            ( *aSortedNodeInds )( 3 ) = tVertices( tN3Ordinal );
            ( *aSortedNodeInds )( 4 ) = ( *aCellIndexIntersectedEdgeVertex )( tIndices( 0 ) );
            ( *aSortedNodeInds )( 5 ) = ( *aCellIndexIntersectedEdgeVertex )( tIndices( 1 ) );
            ( *aSortedNodeInds )( 6 ) = ( *aCellIndexIntersectedEdgeVertex )( tIndices( 2 ) );
            ( *aSortedNodeInds )( 7 ) = ( *aCellIndexIntersectedEdgeVertex )( tIndices( 3 ) );
        }
        else if ( aCellIndexIntersectedEdgeOrdinals->size() == 2 )
        {
            moris_index tLowVertexIdEdgeOrd  = ( *aCellIndexIntersectedEdgeOrdinals )( tIndices( 0 ) );
            moris_index tHighVertexIdEdgeOrd = ( *aCellIndexIntersectedEdgeOrdinals )( tIndices( 1 ) );

            aSortedNodeInds->resize( 6 );
            ( *aSortedNodeInds )( 0 ) = tVertices( 0 );
            ( *aSortedNodeInds )( 1 ) = tVertices( 1 );
            ( *aSortedNodeInds )( 2 ) = tVertices( 2 );
            ( *aSortedNodeInds )( 3 ) = tVertices( 3 );
            ( *aSortedNodeInds )( 4 ) = ( *aCellIndexIntersectedEdgeVertex )( tIndices( 0 ) );
            ( *aSortedNodeInds )( 5 ) = ( *aCellIndexIntersectedEdgeVertex )( tIndices( 1 ) );
            aPermutation              = 10000 + 100 * ( tLowVertexIdEdgeOrd + 1 ) + 10 * ( tHighVertexIdEdgeOrd + 1 );
        }
        else if ( aCellIndexIntersectedEdgeOrdinals->size() == 1 )
        {
            aSortedNodeInds->resize( 5 );
            ( *aSortedNodeInds )( 0 ) = tVertices( 0 );
            ( *aSortedNodeInds )( 1 ) = tVertices( 1 );
            ( *aSortedNodeInds )( 2 ) = tVertices( 2 );
            ( *aSortedNodeInds )( 3 ) = tVertices( 3 );
            ( *aSortedNodeInds )( 4 ) = ( *aCellIndexIntersectedEdgeVertex )( tIndices( 0 ) );
            aPermutation              = 10000 + ( *aCellIndexIntersectedEdgeOrdinals )( tIndices( 0 ) );
        }
    }
    // ----------------------------------------------------------------------------------

}    // namespace xtk
