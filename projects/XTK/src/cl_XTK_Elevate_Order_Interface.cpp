/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Elevate_Order_Interface.cpp
 *
 */

#include "cl_XTK_Decomposition_Algorithm.hpp"
#include "cl_XTK_Elevate_Order_Interface.hpp"
#include "cl_XTK_Integration_Mesh_Generator.hpp"
#include "cl_XTK_Cut_Integration_Mesh.hpp"
#include "cl_MTK_Cell_Info_Factory.hpp"
#include "cl_MTK_Cell_Info.hpp"

namespace moris::xtk
{
    // ----------------------------------------------------------------------------------

    Elevate_Order_Interface::Elevate_Order_Interface(
            Parameter_List&         aParameterList,
            enum Subdivision_Method aSubdivisionMethod )
    {
        // get number of spatial dimensions and decide on subdivision template
        if ( aSubdivisionMethod == Subdivision_Method::P_ELEVATE_ORDER_TRI3_TRI6 )
        {
            mElevateOrderTemplate = std::make_shared< TRI3_to_TRI6 >();
        }
        // else if ( aSubdivisionMethod == Subdivision_Method::P_ELEVATE_ORDER_TRI3_TRI10 )
        // {
        //     mElevateOrderTemplate = std::make_shared< TRI3_to_TRI10 >();
        // }
        else if ( aSubdivisionMethod == Subdivision_Method::P_ELEVATE_ORDER_TET4_TET10 )
        {
            mElevateOrderTemplate = std::make_shared< TET4_to_TET10 >();
        }
        // else if ( aSubdivisionMethod == Subdivision_Method::P_ELEVATE_ORDER_TET4_TET20 )
        // {
        //     mElevateOrderTemplate = std::make_shared< TET4_to_TET20 >();
        // }
        else
        {
            MORIS_ERROR( false, "Elevate_Order_Interface::Elevate_Order_Interface() - order elevation template is unknown" );
        }
    }

    // ----------------------------------------------------------------------------------

    enum Decomposition_Algorithm_Type
    Elevate_Order_Interface::get_algorithm_type() const
    {
        return Decomposition_Algorithm_Type::ELEVATE_ORDER;
    }

    // ----------------------------------------------------------------------------------

    moris_index
    Elevate_Order_Interface::get_signature() const
    {
        return 1337;
    }

    // ----------------------------------------------------------------------------------

    bool
    Elevate_Order_Interface::has_geometric_independent_vertices() const
    {
        return false;
    }

    // ----------------------------------------------------------------------------------

    bool 
    Elevate_Order_Interface::is_eligible( std::pair< mtk::Cell*, Vector< Decomposition_Algorithm_Type > >& aElementContext,
                Cut_Integration_Mesh*                                           aCutIntegrationMesh,
                Integration_Mesh_Generator*                                     aMeshGenerator ) const
    {
        return true;
    }

    // ----------------------------------------------------------------------------------

    Vector< moris_index >
    Elevate_Order_Interface::get_decomposed_cell_indices()
    {
        return mMeshGenerationData->mRegularSubdivisionBgCellInds;
    }

    // ----------------------------------------------------------------------------------

    void
    Elevate_Order_Interface::perform(
            Integration_Mesh_Generation_Data* aMeshGenerationData,
            Decomposition_Data*               aDecompositionData,
            Cut_Integration_Mesh*             aCutIntegrationMesh,
            moris::mtk::Mesh*                 aBackgroundMesh,
            Integration_Mesh_Generator*       aMeshGenerator )
    {
        // log/trace this function
        Tracer tTracer( "XTK", "Decomposition_Algorithm", "Elevate_Order" );

        // keep track of some useful classes (avoid passing to every function)
        mMeshGenerationData = aMeshGenerationData;
        mDecompositionData  = aDecompositionData;
        mCutIntegrationMesh = aCutIntegrationMesh;
        mBackgroundMesh     = aBackgroundMesh;
        mGenerator          = aMeshGenerator;

        // set a new decomposition data
        *aDecompositionData = Decomposition_Data();

        // get list of all cell groups in mesh (which will be treated later)
        Vector< std::shared_ptr< IG_Cell_Group > > tIgCellGroups = mCutIntegrationMesh->get_all_cell_groups();

        // extract the cells from groups into a continuous list
        Vector< moris::mtk::Cell* > tIgCellsInGroups;
        aMeshGenerator->extract_cells_from_cell_groups( tIgCellGroups, tIgCellsInGroups );

        // construct the edge connectivity for these ig cell groups
        std::shared_ptr< Edge_Based_Connectivity > tIgCellGroupEdgeConnectivity = std::make_shared< Edge_Based_Connectivity >();
        aMeshGenerator->create_edges_from_element_to_node( tIgCellsInGroups, tIgCellGroupEdgeConnectivity );

        // collect a representative background cell for each edges
        Vector< moris::mtk::Cell* > tBackgroundCellForEdge;    // input: edge index || output: BG cell it belongs to
        aMeshGenerator->select_background_cell_for_edge( tIgCellGroupEdgeConnectivity, aCutIntegrationMesh, tBackgroundCellForEdge );

        // collect the vertex group related to the representative background cell
        Vector< std::shared_ptr< IG_Vertex_Group > > tVertexGroups;    // input: BG cell index || output: IG vertex group
        aMeshGenerator->collect_vertex_groups_for_background_cells( aMeshGenerationData, aCutIntegrationMesh, &tBackgroundCellForEdge, &tVertexGroups );

        // deduce the edge parent entity index and rank
        std::shared_ptr< Edge_Based_Ancestry > tIgEdgeAncestry = std::make_shared< Edge_Based_Ancestry >();
        aMeshGenerator->deduce_edge_ancestry( aCutIntegrationMesh, aBackgroundMesh, tIgCellGroupEdgeConnectivity, tBackgroundCellForEdge, tIgEdgeAncestry );

        // request vertices and make sure vertices on edges don't get requested twice, additionally build the local
        Vector< moris_index >           tEmptyCellVertexList( mElevateOrderTemplate->get_total_ig_verts(), -1 );
        Vector< Vector< moris_index > > tCellToLocalVertices( tIgCellsInGroups.size(), tEmptyCellVertexList );
        this->make_vertex_requests(
                tIgCellGroupEdgeConnectivity,
                tIgEdgeAncestry,
                &tIgCellsInGroups,
                &tCellToLocalVertices );

        // give all these nodes ids
        aMeshGenerator->assign_node_requests_identifiers( *aDecompositionData, aCutIntegrationMesh, aBackgroundMesh );

        // create associations between child meshes and the vertices we need this to commit the data to the integration mesh
        this->associate_new_vertices_with_cell_groups(
                tIgCellGroupEdgeConnectivity,
                tIgEdgeAncestry,
                &tBackgroundCellForEdge,
                &tIgCellsInGroups );

        // commit vertices to the mesh
        aMeshGenerator->commit_new_ig_vertices_to_cut_mesh( aMeshGenerationData, aDecompositionData, aCutIntegrationMesh, aBackgroundMesh, this );

        // we are ready to create the new integration cells
        this->create_higher_order_integration_cells(
                tIgCellGroupEdgeConnectivity,
                tIgEdgeAncestry,
                &tIgCellsInGroups,
                &tCellToLocalVertices );

        // commit the cells to the mesh
        aMeshGenerator->commit_new_ig_cells_to_cut_mesh( aMeshGenerationData, aDecompositionData, aCutIntegrationMesh, aBackgroundMesh, this );
    }

    // ----------------------------------------------------------------------------------

    bool
    Elevate_Order_Interface::make_vertex_requests(
            const std::shared_ptr< Edge_Based_Connectivity >& aEdgeConnectivity,
            const std::shared_ptr< Edge_Based_Ancestry >&     aIgEdgeAncestry,
            Vector< moris::mtk::Cell* >*                      aIgCells,
            Vector< Vector< moris_index > >*                  aCellToNewLocalVertexIndices )
    {
        Tracer tTracer( "XTK", "Elevate_Order_Interface", "make vertex requests" );

        // get first unused index for nodes for numbering new nodes
        mDecompositionData->mHasSecondaryIdentifier = true;

        // range of indices
        uint tNumElemsInCIM = mCutIntegrationMesh->get_num_entities( mtk::EntityRank::ELEMENT, 0 );

        // initialize proc-global to element-local map of vertex indices for each cell/element
        Vector< std::map< moris_index, uint > > tVertIndicesOnCell( tNumElemsInCIM );

        // initialize iCell-map, provides a relationship between the location of a cell in the provided list of cells and their indices
        std::unordered_map< moris_index, uint > tCellToIndexMap;

        // iterate through the elements get local cell vertex information and, if requested, create vertices on it
        for ( moris::uint iCell = 0; iCell < aIgCells->size(); iCell++ )
        {
            // get index of current cell
            moris_index tCurrentCellIndex = ( *aIgCells )( iCell )->get_index();

            // fill Cell to index map
            tCellToIndexMap[ tCurrentCellIndex ] = iCell;

            // get list of vertex indices for current element / fill list
            for ( uint iVert = 0; iVert < mElevateOrderTemplate->get_num_spatial_dims() + 1; iVert++ )
            {
                // fill map
                tVertIndicesOnCell( tCurrentCellIndex )[ ( *aIgCells )( iCell )->get_vertex_inds()( iVert ) ] = iVert;
            }

            // check if new vertices on faces get created (if true, assume 3D)
            if ( mElevateOrderTemplate->has_new_vertices_on_entity( mtk::EntityRank::ELEMENT ) )
            {
                // todo: assuming only one vertex gets created inside of elements
                // todo: (which is fine for standard quadratic and cubic elements, but maybe not for others)
                MORIS_ERROR( mElevateOrderTemplate->num_new_vertices_per_entity( mtk::EntityRank::ELEMENT ) == 1,
                        "Elevate_Order_Interface::make_vertex_requests() - currently only supports elements with one new vertex inside of each element." );

                // get Cell parent (BG cell) index
                moris_index tCellGroupMembershipIndex = mCutIntegrationMesh->get_ig_cell_group_memberships( tCurrentCellIndex )( 0 );
                moris_index tParentIndex              = mCutIntegrationMesh->get_ig_cell_group_parent_cell( tCellGroupMembershipIndex )->get_index();

                // iterate through vertices to be created inside each
                for ( uint iVert = 0; iVert < mElevateOrderTemplate->num_new_vertices_per_entity( mtk::EntityRank::ELEMENT ); iCell++ )
                {
                    // initialize variable holding possible new node index
                    moris_index tNewNodeIndexInDecompData = MORIS_INDEX_MAX;

                    // check if new node for current edge has already been requested ...
                    bool tRequestExist = mDecompositionData->request_exists(
                            tParentIndex,
                            iVert,    // use local new vertex ordinal as secondary ID, since inside of cell is not shared with other cells
                            mtk::EntityRank::ELEMENT,
                            tNewNodeIndexInDecompData );

                    // ... if not, request it
                    if ( !tRequestExist )
                    {
                        // find out which processor owns parent entity of currently treated edge
                        moris::moris_index tOwningProc = mBackgroundMesh->get_entity_owner( tParentIndex, mtk::EntityRank::ELEMENT );

                        // compute global coordinates for new node
                        Matrix< DDRMat > tNewVertexCoords;
                        if ( mElevateOrderTemplate->get_num_spatial_dims() == 2 )
                        {
                            tNewVertexCoords = this->compute_tri_vertex_global_coordinates(
                                    ( *aIgCells )( iCell )->get_vertex_pointers(),
                                    mElevateOrderTemplate->get_new_vertex_parametric_coords_wrt_entity( mtk::EntityRank::ELEMENT )( iVert ) );
                        }
                        else if ( mElevateOrderTemplate->get_num_spatial_dims() == 3 )
                        {
                            tNewVertexCoords = this->compute_tet_vertex_global_coordinates(
                                    ( *aIgCells )( iCell )->get_vertex_pointers(),
                                    mElevateOrderTemplate->get_new_vertex_parametric_coords_wrt_entity( mtk::EntityRank::ELEMENT )( iVert ) );
                        }
                        else
                        {
                            MORIS_ERROR( false, "Elevate_Order_Interface::make_vertex_requests() - num spatial dimensions returned by template must be 2 or 3." );
                        }

                        // Register new node request
                        tNewNodeIndexInDecompData = mDecompositionData->register_new_request(
                                tParentIndex,
                                iVert,    // use local new vertex ordinal as secondary ID, since inside of cell is not shared with other cells
                                tOwningProc,
                                mtk::EntityRank::ELEMENT,
                                tNewVertexCoords );

                        // count number of new nodes created
                        // tNewNodeIndex++;
                    }
                }    // end: loop over new vertices inside each cell
            }        // end: new vertices inside element
        }            // end: loop inside element

        // --------------------------------

        // check if new vertices on edges get created
        if ( mElevateOrderTemplate->has_new_vertices_on_entity( mtk::EntityRank::EDGE ) )
        {
            // get number of new vertices to be created on edge
            uint tNumNewVerticesPerEdge = mElevateOrderTemplate->num_new_vertices_per_entity( mtk::EntityRank::EDGE );

            // iterate through the edges in aEdgeConnectivity ask the geometry engine if we are intersected
            for ( moris::uint iEdge = 0; iEdge < aEdgeConnectivity->mEdgeVertices.size(); iEdge++ )
            {
                // get edge parent entity index and rank
                moris_index tParentIndex = aIgEdgeAncestry->mEdgeParentEntityIndex( iEdge );
                moris_index tParentRank  = aIgEdgeAncestry->mEdgeParentEntityRank( iEdge );

                // get unique edge id based on two end vertices of edge
                moris_index tSecondaryId = this->hash_edge( aEdgeConnectivity->mEdgeVertices( iEdge ) );

                // initialize variable holding possible new node index
                moris_index tNewNodeIndexInDecompData = MORIS_INDEX_MAX;

                // check if new node for current edge has already been requested ...
                bool tRequestExist = mDecompositionData->request_exists(
                        tParentIndex,
                        tSecondaryId,
                        (mtk::EntityRank)tParentRank,
                        tNewNodeIndexInDecompData );

                // ... if not, request it
                if ( !tRequestExist )
                {
                    // find out which processor owns parent entity of currently treated edge
                    moris::moris_index tOwningProc = mBackgroundMesh->get_entity_owner( tParentIndex, (mtk::EntityRank)tParentRank );

                    // collect new vertices' indices in decomp data
                    Matrix< IndexMat > tNewEdgeVertexIndicesInDecompData( tNumNewVerticesPerEdge, 1, -1 );

                    // iterate through vertices to be created on each edge
                    for ( uint iVert = 0; iVert < tNumNewVerticesPerEdge; iVert++ )
                    {
                        // compute global coordinates for new node
                        Matrix< DDRMat > tNewVertexCoords =
                                this->compute_edge_vertex_global_coordinates(
                                        aEdgeConnectivity->mEdgeVertices( iEdge ),
                                        mElevateOrderTemplate->get_new_vertex_parametric_coords_wrt_entity( mtk::EntityRank::EDGE )( iVert ) );

                        // fixme: can two nodes be requested with the same secondary ID?
                        // Register new node request
                        tNewNodeIndexInDecompData = mDecompositionData->register_new_request(
                                tParentIndex,
                                tSecondaryId,
                                tOwningProc,
                                (mtk::EntityRank)tParentRank,
                                tNewVertexCoords );

                        // record new vertex index in decomp data
                        tNewEdgeVertexIndicesInDecompData( iVert ) = tNewNodeIndexInDecompData;

                        // count number of new nodes created
                        // tNewNodeIndex++;

                    }    // end: loop over all new vertices on edge

                    // get index of edge vertices on proc
                    moris_index tFirstVertexIndex  = aEdgeConnectivity->mEdgeVertices( iEdge )( 0 )->get_index();
                    moris_index tSecondVertexIndex = aEdgeConnectivity->mEdgeVertices( iEdge )( 1 )->get_index();

                    // get number of cells attached to edge, to register new vertices to these cells
                    uint tNumCellsAttachedToEdge = aEdgeConnectivity->mEdgeToCell( iEdge ).size();

                    // go over cells attached to edge and register decomp data vertex indices for local vertices
                    for ( uint iCellAttachedToEdge = 0; iCellAttachedToEdge < tNumCellsAttachedToEdge; iCellAttachedToEdge++ )
                    {
                        // get Cell's index on proc
                        moris_index tCellIndex = aEdgeConnectivity->mEdgeToCell( iEdge )( iCellAttachedToEdge )->get_index();

                        // get the element local indices of the edge start and end vertices
                        uint tFirstVertLocalIndex  = tVertIndicesOnCell( tCellIndex ).find( tFirstVertexIndex )->second;
                        uint tSecondVertLocalIndex = tVertIndicesOnCell( tCellIndex ).find( tSecondVertexIndex )->second;

                        // get 1-based edge index and direction
                        moris_index tSignedEdgeIndex =
                                mElevateOrderTemplate->get_local_edge_index_based_on_vertex_indices( tFirstVertLocalIndex, tSecondVertLocalIndex );

                        // iterate through vertices to be created on each edge, and save decomp data index for new edges on each element
                        for ( uint iVert = 0; iVert < tNumNewVerticesPerEdge; iVert++ )
                        {
                            // get vertex' element local index
                            moris_index tLocalVertexIndex = mElevateOrderTemplate->get_local_vertex_index( mtk::EntityRank::EDGE, tSignedEdgeIndex, iVert );

                            // get integration cell index
                            uint iCell = tCellToIndexMap.find( tCellIndex )->second;

                            // store to decomp data index
                            ( *aCellToNewLocalVertexIndices )( iCell )( tLocalVertexIndex ) = tNewEdgeVertexIndicesInDecompData( iVert );
                        }
                    }
                }    // end: check for new request
            }        // end: loop over all edges
        }            // end: new vertices on edges

        // --------------------------------

        // check if new vertices on faces get created (if true, assume 3D)
        if ( mElevateOrderTemplate->has_new_vertices_on_entity( mtk::EntityRank::FACE ) )
        {
            // get face connectivity from Cut Integration Mesh
            std::shared_ptr< Facet_Based_Connectivity > tFaceConnectivity = mCutIntegrationMesh->get_face_connectivity();
            std::shared_ptr< Facet_Based_Ancestry >     tIgFaceAncestry   = mCutIntegrationMesh->get_face_ancestry();

            // todo: assuming only one vertex gets created per face (which is fine for standard quadratic and cubic elements, but maybe not for others)
            MORIS_ERROR( mElevateOrderTemplate->num_new_vertices_per_entity( mtk::EntityRank::FACE ) == 1,
                    "Elevate_Order_Interface::make_vertex_requests() - currently only supports elements with one new vertex per face max." );

            // iterate through the faces and create vertices on it
            for ( moris::uint iFace = 0; iFace < tFaceConnectivity->mFacetVertices.size(); iFace++ )
            {
                // get edge parent entity index and rank
                moris_index tParentIndex = tIgFaceAncestry->mFacetParentEntityIndex( iFace );
                moris_index tParentRank  = tIgFaceAncestry->mFacetParentEntityRank( iFace );

                // get unique face id based on three corner vertices of edge
                moris_index tSecondaryId = this->hash_face( tFaceConnectivity->mFacetVertices( iFace ) );

                // initialize variable holding possible new node index
                moris_index tNewNodeIndexInDecompData = MORIS_INDEX_MAX;

                // check if new node for current edge has already been requested ...
                bool tRequestExist = mDecompositionData->request_exists(
                        tParentIndex,
                        tSecondaryId,
                        (mtk::EntityRank)tParentRank,
                        tNewNodeIndexInDecompData );

                // ... if not, request it
                if ( !tRequestExist )
                {
                    // find out which processor owns parent entity of currently treated edge
                    moris::moris_index tOwningProc = mBackgroundMesh->get_entity_owner( tParentIndex, (mtk::EntityRank)tParentRank );

                    // compute global coordinates for new node
                    Matrix< DDRMat > tNewVertexCoords =
                            this->compute_tri_vertex_global_coordinates(
                                    tFaceConnectivity->mFacetVertices( iFace ),
                                    mElevateOrderTemplate->get_new_vertex_parametric_coords_wrt_entity( mtk::EntityRank::FACE )( 0 ) );

                    // Register new node request
                    tNewNodeIndexInDecompData = mDecompositionData->register_new_request(
                            tParentIndex,
                            tSecondaryId,
                            tOwningProc,
                            (mtk::EntityRank)tParentRank,
                            tNewVertexCoords );

                    // count number of new nodes created
                    // tNewNodeIndex++;
                }
            }    // end: loop over faces
        }        // end: new vertices on faces

        // --------------------------------

        // return success if finished
        return true;
    }

    // ----------------------------------------------------------------------------------

    bool
    Elevate_Order_Interface::associate_new_vertices_with_cell_groups(
            const std::shared_ptr< Edge_Based_Connectivity >& aEdgeConnectivity,
            const std::shared_ptr< Edge_Based_Ancestry >&     aIgEdgeAncestry,
            Vector< moris::mtk::Cell* >*                      aBackgroundCellForEdge,
            Vector< moris::mtk::Cell* >*                      aIgCells )
    {
        // trace this function
        Tracer tTracer( "XTK", "Decomposition_Algorithm", "Vertex Associations" );

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

        Vector< moris::uint > tCellGroups;
        tCellGroups.reserve( 20 );

        // iterate through the edges
        for ( moris::uint iEdge = 0; iEdge < aEdgeConnectivity->mEdgeVertices.size(); iEdge++ )
        {
            moris_index tEdgeIndex = (moris_index)iEdge;

            // iterate through elements attached to this edge and collect the cell groups
            tCellGroups.clear();
            for ( moris::uint iCell = 0; iCell < aEdgeConnectivity->mEdgeToCell( tEdgeIndex ).size(); iCell++ )
            {
                // integration cell just grab the first
                moris::mtk::Cell* tCell = aEdgeConnectivity->mEdgeToCell( tEdgeIndex )( iCell );

                // fill list of Cell groups the edge is related to
                tCellGroups.push_back( mCutIntegrationMesh->get_ig_cell_group_memberships( tCell->get_index() )( 0 ) );
            }

            // remove duplicates from list
            moris::unique( tCellGroups );

            // iterate through the unique cell groups
            for ( moris::uint iUCG = 0; iUCG < tCellGroups.size(); iUCG++ )
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
                        Interpolation::linear_interpolation_location(
                                tEdgeNodeParamCoordinates,
                                mElevateOrderTemplate->get_new_vertex_parametric_coords_wrt_entity( mtk::EntityRank::EDGE )( 0 ) );

                // list of edges and locations to create new nodes on added to decomp. data
                mDecompositionData->tCMNewNodeLoc( tCellGroupMembershipIndex ).push_back( iEdge );
                mDecompositionData->tCMNewNodeParamCoord( tCellGroupMembershipIndex ).push_back( tParametricCoordsRelativeToParentElem );
            }
        }

        // --------------------------------

        // return success if finished
        return true;
    }

    // ----------------------------------------------------------------------------------

    void
    Elevate_Order_Interface::create_higher_order_integration_cells(
            const std::shared_ptr< Edge_Based_Connectivity >& aEdgeConnectivity,
            const std::shared_ptr< Edge_Based_Ancestry >&     aIgEdgeAncestry,
            Vector< moris::mtk::Cell* >*                      aIgCells,
            Vector< Vector< moris_index > >*                  aCellToNewLocalVertexIndices )
    {
        // time/log function
        Tracer tTracer( "XTK", "Elevate_Order_Interface", "Create Higher Order Integration Cells" );

        // starting with a clean slate here
        mNumNewCells = 0;

        // gather information needed for new cells from order elevation template
        uint        tNumIgCellsInMesh = aIgCells->size();
        moris_index tNodesPerCell     = mElevateOrderTemplate->get_total_ig_verts();

        // create cell info for new cells
        moris::mtk::Cell_Info_Factory            tFactory;
        std::shared_ptr< moris::mtk::Cell_Info > tCellInfo = tFactory.create_cell_info_sp( mElevateOrderTemplate->get_ig_cell_topology() );

        // initialize the algorithm data
        mNewCellToVertexConnectivity = Vector< Vector< moris::moris_index > >( tNumIgCellsInMesh, Vector< moris::moris_index >( tNodesPerCell ) );
        mNewCellChildMeshIndex       = Vector< moris::moris_index >( tNumIgCellsInMesh );
        mNewCellCellIndexToReplace   = Vector< moris::moris_index >( tNumIgCellsInMesh );
        mNewCellCellInfo             = Vector< std::shared_ptr< moris::mtk::Cell_Info > >( tNumIgCellsInMesh, tCellInfo );

        // construct the Cell-Vertex-Connectivity for all new cells
        for ( moris::uint iCell = 0; iCell < tNumIgCellsInMesh; iCell++ )
        {
            // new cell replaces current cell, hence it has the same index
            mNewCellCellIndexToReplace( iCell ) = ( *aIgCells )( iCell )->get_index();

            // get cell group membership of current cell and assign it to new cell
            moris_index tCellGroupMembershipIndex = mCutIntegrationMesh->get_ig_cell_group_memberships( ( *aIgCells )( iCell )->get_index() )( 0 );
            mNewCellChildMeshIndex( iCell )       = tCellGroupMembershipIndex;

            // get list of vertices belonging to new cell, sorted in
            Vector< moris_index > tSortedVerticesForCell( mElevateOrderTemplate->get_total_ig_verts() );

            // iterate through old corner vertices on cell, get their indices, and construct the Cell-Vertex-Connectivity map
            for ( uint iVert = 0; iVert < mElevateOrderTemplate->get_num_spatial_dims() + 1; iVert++ )
            {
                // get vertex index of corner vertices from old element
                moris_index tOldVertexIndex = ( *aIgCells )( iCell )->get_vertex_inds()( iVert );

                // copy value onto IEN
                mNewCellToVertexConnectivity( iCell )( iVert ) = tOldVertexIndex;
            }

            // iterate through new vertices on cell, get their indices, and construct the Cell-Vertex-Connectivity map
            for ( uint iVert = mElevateOrderTemplate->get_num_spatial_dims() + 1; iVert < mElevateOrderTemplate->get_total_ig_verts(); iVert++ )
            {
                // get new vertex index from decomp data
                moris_index tNewVertexIndex = mDecompositionData->tNewNodeIndex( ( *aCellToNewLocalVertexIndices )(iCell)( iVert ) );

                // copy value onto IEN
                mNewCellToVertexConnectivity( iCell )( iVert ) = tNewVertexIndex;
            }
        }
    }

    // ----------------------------------------------------------------------------------
    // ----------------------------------------------------------------------------------

    moris_index
    Elevate_Order_Interface::hash_edge( Vector< moris::mtk::Vertex* > const & aEdgeVertices )
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

    // ----------------------------------------------------------------------------------

    moris_index
    Elevate_Order_Interface::hash_face( Vector< moris::mtk::Vertex* > const & aFaceVertices )
    {
        // check input
        MORIS_ERROR( aFaceVertices.size() == 3, "Face is expected to have three vertices" );

        const moris_id tVertId0 = aFaceVertices( 0 )->get_id();
        const moris_id tVertId1 = aFaceVertices( 1 )->get_id();
        const moris_id tVertId2 = aFaceVertices( 2 )->get_id();

        moris_index tMinIdIndex = 0;
        moris_index tMidIdIndex = 1;
        moris_index tMaxIdIndex = 2;

        // sort indices of vertices
        if ( tVertId0 > tVertId2 )
        {
            swap_indices( tMinIdIndex, tMaxIdIndex );
        }
        if ( tVertId0 > tVertId1 )
        {
            swap_indices( tMinIdIndex, tMidIdIndex );
        }
        if ( tVertId1 > tVertId2 )
        {
            swap_indices( tMidIdIndex, tMaxIdIndex );
        }

        // compute recursive cantor 3-tupel
        moris_index tFirstPairVal  = xtk::cantor_pairing( aFaceVertices( tMinIdIndex )->get_id(), aFaceVertices( tMidIdIndex )->get_id() );
        moris_index tSecondPairVal = xtk::cantor_pairing( tFirstPairVal, aFaceVertices( tMaxIdIndex )->get_id() );

        // check for likely overflow
        MORIS_ASSERT( tSecondPairVal > MORIS_INDEX_MAX / 10, "Elevate_Order_Interface::hash_face() - function is likely to lead to an overflow." );

        // return unique id for face
        return tSecondPairVal;
    }

    // ----------------------------------------------------------------------------------

    void
    Elevate_Order_Interface::swap_indices( moris_index& aInd1, moris_index& aInd2 )
    {
        moris_index tStore = aInd1;
        aInd1              = aInd2;
        aInd2              = tStore;
    }

    // -------------------------------------------------------------------------

    Matrix< DDRMat >
    Elevate_Order_Interface::compute_edge_vertex_global_coordinates(
            Vector< moris::mtk::Vertex* > const & aEdgeVertices,
            Matrix< DDRMat >                      aEdgeCoordinate )
    {
        // check inputs
        MORIS_ASSERT( aEdgeVertices.size() == 2,
                "Elevate_Order_Interface::compute_edge_vertex_global_coordinates() - two corner vertices must be provided." );
        MORIS_ASSERT( ( aEdgeCoordinate( 0 ) > -1.0 - 10.0 * MORIS_REAL_EPS ) && ( aEdgeCoordinate( 0 ) < 1.0 + 10.0 * MORIS_REAL_EPS ),
                "Elevate_Order_Interface::compute_edge_vertex_global_coordinates() - provided edge local coord out of bounds (-1,1)." );

        // compute interpolation
        real             tXi = aEdgeCoordinate( 0 );
        Matrix< DDRMat > tIntersectionCoords =
                0.5 * ( ( 1.0 - tXi ) * aEdgeVertices( 0 )->get_coords() + ( 1.0 + tXi ) * aEdgeVertices( 1 )->get_coords() );

        // return vertex coords
        return tIntersectionCoords;
    }

    // -------------------------------------------------------------------------

    Matrix< DDRMat >
    Elevate_Order_Interface::compute_tri_vertex_global_coordinates(
            Vector< moris::mtk::Vertex* > const & aTriVertices,
            const Matrix< DDRMat >&               aTriCoords )
    {
        MORIS_ERROR( false, "Elevate_Order_Interface::compute_tri_vertex_global_coordinates() - function not implemented yet." );
        return { { 0.0 } };
    }

    // -------------------------------------------------------------------------

    Matrix< DDRMat >
    Elevate_Order_Interface::compute_tet_vertex_global_coordinates(
            Vector< moris::mtk::Vertex* > const & aTriVertices,
            const Matrix< DDRMat >&               aTriCoords )
    {
        MORIS_ERROR( false, "Elevate_Order_Interface::compute_tet_vertex_global_coordinates() - function not implemented yet." );
        return { { 0.0 } };
    }

    // ----------------------------------------------------------------------------------

}    // namespace moris::xtk
