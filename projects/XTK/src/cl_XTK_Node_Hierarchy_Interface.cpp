#include "cl_XTK_Node_Hierarchy_Interface.hpp"
#include "cl_XTK_Decomposition_Algorithm.hpp"
#include "cl_XTK_Integration_Mesh_Generator.hpp"
#include "cl_XTK_Cut_Integration_Mesh.hpp"
#include "cl_MTK_Cell_Info_Factory.hpp"
#include "cl_MTK_Cell_Info.hpp"
#include "fn_Pairing.hpp"
#include "cl_GEN_Geometry_Engine.hpp"
#include <algorithm>// std::sort, std::stable_sort
#include <numeric>
#include "cl_Tracer.hpp"
#include <chrono>
namespace xtk
{
enum Decomposition_Algorithm_Type
Node_Hierarchy_Interface::get_algorithm_type() const
{
    return Decomposition_Algorithm_Type::NODE_HEIRARCHY;
}
moris_index
Node_Hierarchy_Interface::get_signature() const
{
    return 102;
}

bool
Node_Hierarchy_Interface::has_geometric_independent_vertices() const
{
    return false;
}

void
Node_Hierarchy_Interface::perform(
    Integration_Mesh_Generation_Data* aMeshGenerationData,
    Decomposition_Data*               aDecompositionData,
    Cut_Integration_Mesh*             aCutIntegrationMesh,
    moris::mtk::Mesh*                 aBackgroundMesh,
    Integration_Mesh_Generator*       aMeshGenerator )
{
    Tracer tTracer( "XTK", "Decomposition_Algorithm", "Node_Hierarchy_Interface" );
    // keep track of some useful classes (avoid passing to every function)
    mGeometryEngine     = aMeshGenerator->get_geom_engine();
    mMeshGenerationData = aMeshGenerationData;
    mDecompositionData  = aDecompositionData;
    mCutIntegrationMesh = aCutIntegrationMesh;
    mBackgroundMesh     = aBackgroundMesh;
    mGenerator          = aMeshGenerator;

    // moris::print(aBackgroundMesh->get_mtk_cell(0).get_vertex_coords(),"BG cell vertex coords");
    // moris::print(aBackgroundMesh->get_mtk_cell(0).get_vertex_inds(),  "BG cell inds");

    // // print the geometry data
    // mCurrentGeomIndex = 0;
    // for(moris::uint iV = 0; iV < aCutIntegrationMesh->get_num_entities(EntityRank::NODE,0); iV++)
    // {
    //     moris::mtk::Vertex* tVertex = aCutIntegrationMesh->get_mtk_vertex_pointer(iV);

    //     std::cout<<"tVertex = "<<tVertex->get_index()<<" | Geom Val = "<< mGeometryEngine->get_field_value(mCurrentGeomIndex, (uint)tVertex->get_index(), tVertex->get_coords())<<std::endl;
    // }

    // active geometries
    moris::Matrix< moris::IndexMat > const* tActiveGeometries = aMeshGenerator->get_active_geometries();

    // iterate through geometries
    for ( moris::uint iGeom = 0; iGeom < tActiveGeometries->numel(); iGeom++ )
    {
        // set a new decomposition data
        *aDecompositionData = Decomposition_Data();

        // set the current geometry index
        mCurrentGeomIndex = ( *tActiveGeometries )( iGeom );

        aDecompositionData->tDecompId = 10000 * iGeom + this->get_signature();

        // cell groups relevant to this geometry pick them out
        moris::Cell< std::shared_ptr< IG_Cell_Group > > tIgCellGroups;
        this->select_ig_cell_groups( tIgCellGroups );

        // extract the cells from groups into a continuos list
        moris::Cell< moris::mtk::Cell* > tIgCellsInGroups;
        aMeshGenerator->extract_cells_from_cell_groups( tIgCellGroups, tIgCellsInGroups );

        // construct the edge connectivity for these ig cell groups
        std::shared_ptr< Edge_Based_Connectivity > tIgCellGroupEdgeConnectivity = std::make_shared< Edge_Based_Connectivity >();
        aMeshGenerator->create_edges_from_element_to_node( tIgCellsInGroups, tIgCellGroupEdgeConnectivity );

        // collect a representative background cell for each edges
        moris::Cell< moris::mtk::Cell* > tBackgroundCellForEdge;
        aMeshGenerator->select_background_cell_for_edge( tIgCellGroupEdgeConnectivity, aCutIntegrationMesh, tBackgroundCellForEdge );

        // collect the vertex group related to the representative background cell
        moris::Cell< std::shared_ptr< IG_Vertex_Group > > tVertexGroups;
        aMeshGenerator->collect_vertex_groups_for_background_cells( aMeshGenerationData, aCutIntegrationMesh, &tBackgroundCellForEdge, &tVertexGroups );

        // deduce the edge parent entity index and rank
        std::shared_ptr< Edge_Based_Ancestry > tIgEdgeAncestry = std::make_shared< Edge_Based_Ancestry >();
        aMeshGenerator->deduce_edge_ancestry( aCutIntegrationMesh, aBackgroundMesh, tIgCellGroupEdgeConnectivity, tBackgroundCellForEdge, tIgEdgeAncestry );

        // figure out which edges are intersected
        moris::Cell< moris_index > tIntersectedEdgeIndices;
        moris::Cell< moris::real > tIntersectedEdgeLocCoords;
        this->determine_intersected_edges_and_make_requests(
            tIgCellGroupEdgeConnectivity,
            tIgEdgeAncestry,
            &tBackgroundCellForEdge,
            &tVertexGroups,
            tIntersectedEdgeIndices,
            tIntersectedEdgeLocCoords );

        // give all these nodes ids
        aMeshGenerator->assign_node_requests_identifiers( *aDecompositionData, aCutIntegrationMesh, aBackgroundMesh, this->get_signature() );

        // create associatios between child meshes and the vertices we need this to commit the data to the integration mesh
        this->associate_new_vertices_with_cell_groups(
            tIgCellGroupEdgeConnectivity,
            tIgEdgeAncestry,
            &tBackgroundCellForEdge,
            &tVertexGroups,
            &tIntersectedEdgeIndices,
            &tIntersectedEdgeLocCoords );

        // commit vertices to the mesh
        aMeshGenerator->commit_new_ig_vertices_to_cut_mesh( aMeshGenerationData, aDecompositionData, aCutIntegrationMesh, aBackgroundMesh, this );

        // we are ready to create the new integration cells
        this->create_node_hierarchy_integration_cells( tIgCellGroupEdgeConnectivity, tIgEdgeAncestry, &tIntersectedEdgeIndices );

        // commit the cells to the mesh
        aMeshGenerator->commit_new_ig_cells_to_cut_mesh( aMeshGenerationData, aDecompositionData, aCutIntegrationMesh, aBackgroundMesh, this );
    }
}

bool
Node_Hierarchy_Interface::determine_intersected_edges_and_make_requests(
    std::shared_ptr< Edge_Based_Connectivity >         aEdgeConnectivity,
    std::shared_ptr< Edge_Based_Ancestry >             aIgEdgeAncestry,
    moris::Cell< moris::mtk::Cell* >*                  aBackgroundCellForEdge,
    moris::Cell< std::shared_ptr< IG_Vertex_Group > >* aVertexGroups,
    moris::Cell< moris_index >&                        aIntersectedEdges,
    moris::Cell< moris::real >&                        aEdgeLocalCoordinate )
{
    Tracer tTracer( "XTK", "Decomposition_Algorithm", "Determine Intersected Edges" );

    moris::moris_index tNewNodeIndex = mCutIntegrationMesh->get_first_available_index( EntityRank::NODE );

    aIntersectedEdges.clear();
    aEdgeLocalCoordinate.clear();

    aIntersectedEdges.reserve( aEdgeConnectivity->mEdgeVertices.size() );
    aEdgeLocalCoordinate.reserve( aEdgeConnectivity->mEdgeVertices.size() );

    // the query interface
    // Initialize geometric query
    Geometric_Query_XTK tGeometricQuery;

    // setup the query data (fixed parts for this function)

    tGeometricQuery.set_query_type( moris::ge::Query_Type::INTERSECTION_LOCATION );
    tGeometricQuery.set_coordinates_matrix( mCutIntegrationMesh->get_all_vertex_coordinates_loc_inds() );
    tGeometricQuery.set_cut_integration_mesh( mCutIntegrationMesh );
    tGeometricQuery.set_query_entity_rank( EntityRank::EDGE );
    tGeometricQuery.set_edge_connectivity( aEdgeConnectivity );
    tGeometricQuery.set_edge_associated_background_cell( aBackgroundCellForEdge );
    tGeometricQuery.set_associated_vertex_group( aVertexGroups );
    tGeometricQuery.set_geometric_index( mCurrentGeomIndex );

    mDecompositionData->mHasSecondaryIdentifier = true;

    //iterate through the edges in aEdgeConnectivity ask the geometry engine if we are intersected
    for ( moris::uint iEdge = 0; iEdge < aEdgeConnectivity->mEdgeVertices.size(); iEdge++ )
    {
        if ( iEdge == 0 )
        {
            tGeometricQuery.set_parametric_coordinate_size( ( *aBackgroundCellForEdge )( iEdge )->get_cell_info()->get_loc_coord_dim() );
        }
        // update the current edge index
        tGeometricQuery.set_current_edge_index( iEdge );

        // change out the parent cell
        tGeometricQuery.set_parent_cell( ( *aBackgroundCellForEdge )( iEdge ) );

        // see if the edge is intersected using the geometry engine
        bool tIsIntersected = mGeometryEngine->geometric_query( &tGeometricQuery );

        if ( tIsIntersected )
        {
            if ( !mGeometryEngine->queued_intersection_first_parent_on_interface() && !mGeometryEngine->queued_intersection_second_parent_on_interface() )
            {

                aIntersectedEdges.push_back( (moris_index)iEdge );
                aEdgeLocalCoordinate.push_back( mGeometryEngine->get_queued_intersection_local_coordinate() );

                moris_index tParentIndex               = aIgEdgeAncestry->mEdgeParentEntityIndex( iEdge );
                moris_index tParentRank                = aIgEdgeAncestry->mEdgeParentEntityRank( iEdge );
                moris_index tSecondaryId               = this->hash_edge( aEdgeConnectivity->mEdgeVertices( iEdge ) );
                moris_index tNewNodeIndexInSubdivision = MORIS_INDEX_MAX;
                bool        tRequestExist              = mDecompositionData->request_exists( tParentIndex, tSecondaryId, (enum EntityRank)tParentRank, tNewNodeIndexInSubdivision );

                if ( !tRequestExist )
                {
                    moris::moris_index tOwningProc = mBackgroundMesh->get_entity_owner( tParentIndex, (enum EntityRank)tParentRank );


                    // Register new request
                    tNewNodeIndexInSubdivision = mDecompositionData->register_new_request(
                        tParentIndex,
                        tSecondaryId,
                        tOwningProc,
                        (enum EntityRank)tParentRank,
                        mGeometryEngine->get_queued_intersection_global_coordinates() );

                    mGeometryEngine->admit_queued_intersection( tNewNodeIndex );
                    tNewNodeIndex++;
                }
            }
        }
    }

    return true;
}

bool
Node_Hierarchy_Interface::associate_new_vertices_with_cell_groups(
    std::shared_ptr< Edge_Based_Connectivity >         aEdgeConnectivity,
    std::shared_ptr< Edge_Based_Ancestry >             aIgEdgeAncestry,
    moris::Cell< moris::mtk::Cell* >*                  aBackgroundCellForEdge,
    moris::Cell< std::shared_ptr< IG_Vertex_Group > >* aVertexGroups,
    moris::Cell< moris_index >*                        aIntersectedEdges,
    moris::Cell< moris::real >*                        aEdgeLocalCoordinate )
{
    Tracer tTracer( "XTK", "Decomposition_Algorithm", "Vertex Associations" );
    // number of child meshes
    moris_index tNumChildMeshes = mMeshGenerationData->mIntersectedBackgroundCellIndexToChildMeshIndex.size();

    // estimate required space for the data in decomposition data
    mDecompositionData->tCMNewNodeLoc        = Cell< Cell< moris_index > >( tNumChildMeshes );
    mDecompositionData->tCMNewNodeParamCoord = Cell< Cell< Matrix< DDRMat > > >( tNumChildMeshes );

    moris_index tDimParamCoords = 0;
    if ( aBackgroundCellForEdge->size() > 0 )
    {
        tDimParamCoords = ( *aBackgroundCellForEdge )( 0 )->get_cell_info()->get_loc_coord_dim();
    }

    moris::Matrix< moris::DDRMat > tEdgeNodeParamCoordinates( 2, tDimParamCoords );

    Cell< moris::uint > tCellGroups;
    tCellGroups.reserve( 20 );

    // iterate through the edges
    for ( moris::uint iEdge = 0; iEdge < aIntersectedEdges->size(); iEdge++ )
    {
        moris_index tEdgeIndex = ( *aIntersectedEdges )( iEdge );

        // iterate through elements attached to this edge and collect the cell groups
        tCellGroups.clear();
        for ( moris::uint iCell = 0; iCell < aEdgeConnectivity->mEdgeToCell( tEdgeIndex ).size(); iCell++ )
        {
            // integration cell just grab the first
            moris::mtk::Cell* tCell = aEdgeConnectivity->mEdgeToCell( tEdgeIndex )( iCell );

            tCellGroups.push_back( mCutIntegrationMesh->get_ig_cell_group_memberships( tCell->get_index() )( 0 ) );
        }

        moris::unique( tCellGroups );

        //iterate through the unique cell groups
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

            tEdgeNodeParamCoordinates.set_row( 0, *tVertex0LocalCoords );
            tEdgeNodeParamCoordinates.set_row( 1, *tVertex1LocalCoords );

            moris::Matrix< moris::DDRMat > tParametricCoordsRelativeToParentElem =
                Interpolation::linear_interpolation_location( tEdgeNodeParamCoordinates, { { ( *aEdgeLocalCoordinate )( iEdge ) } } );

            mDecompositionData->tCMNewNodeLoc( tCellGroupMembershipIndex ).push_back( iEdge );
            mDecompositionData->tCMNewNodeParamCoord( tCellGroupMembershipIndex ).push_back( tParametricCoordsRelativeToParentElem );
        }
    }

    return true;
}


moris_index
Node_Hierarchy_Interface::hash_edge( moris::Cell< moris::mtk::Vertex* > const& aEdgeVertices )
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
    moris::mtk::Mesh*                 aBackgroundMesh,
    Integration_Mesh_Generator*       aMeshGenerator )
{
}

void
Node_Hierarchy_Interface::perform_impl_generate_mesh(
    Integration_Mesh_Generation_Data* aMeshGenerationData,
    Decomposition_Data*               aDecompositionData,
    Cut_Integration_Mesh*             aCutIntegrationMesh,
    moris::mtk::Mesh*                 aBackgroundMesh,
    Integration_Mesh_Generator*       aMeshGenerator )
{
}

void
Node_Hierarchy_Interface::select_ig_cell_groups(
    moris::Cell< std::shared_ptr< IG_Cell_Group > >& aIgCellGroups )
{
    Tracer tTracer( "XTK", "Decomposition_Algorithm", "Select Ig Cell Groups" );
    // intersected background cells for the current geometry
    moris::Cell< moris_index > const& tIntersectedBackground = mMeshGenerationData->mIntersectedBackgroundCellIndex( mCurrentGeomIndex );

    // number of intersected background cells
    moris::uint tNumIntersectedBackground = tIntersectedBackground.size();

    // allocate input data
    aIgCellGroups = moris::Cell< std::shared_ptr< IG_Cell_Group > >( tNumIntersectedBackground, nullptr );

    // iterate through them and get the cell group index
    for ( moris::uint iBGCell = 0; iBGCell < tNumIntersectedBackground; iBGCell++ )
    {
        moris_index tCurrentBgCellIndex = tIntersectedBackground( iBGCell );

        auto tIter = mMeshGenerationData->mIntersectedBackgroundCellIndexToChildMeshIndex.find( tCurrentBgCellIndex );

        MORIS_ERROR( tIter != mMeshGenerationData->mIntersectedBackgroundCellIndexToChildMeshIndex.end(), "Issue with background cell data and allocated child mesh data" );

        moris_index tCMIndex = tIter->second;

        aIgCellGroups( iBGCell ) = mCutIntegrationMesh->get_ig_cell_group( tCMIndex );
    }
}

void
Node_Hierarchy_Interface::create_node_hierarchy_integration_cells(
    std::shared_ptr< Edge_Based_Connectivity > aEdgeConnectivity,
    std::shared_ptr< Edge_Based_Ancestry >     aIgEdgeAncestry,
    moris::Cell< moris_index >*                aIntersectedEdges )
{
    Tracer tTracer( "XTK", "Node_Hierarchy_Interface", "Create NH IG Cells" );
    // starting with a clean slate here
    mNumNewCells = 0;

    moris::Cell< std::shared_ptr< moris::Cell< moris_index > > >         tCellIndexIntersectedEdgeOrdinals( mCutIntegrationMesh->get_num_entities( EntityRank::ELEMENT, 0 ), nullptr );
    moris::Cell< std::shared_ptr< moris::Cell< moris::mtk::Vertex* > > > tCellIndexIntersectedEdgeVertex( mCutIntegrationMesh->get_num_entities( EntityRank::ELEMENT, 0 ), nullptr );
    this->determine_intersected_cell_information( aEdgeConnectivity, aIntersectedEdges, &tCellIndexIntersectedEdgeOrdinals, &tCellIndexIntersectedEdgeVertex );

    moris::Cell< std::shared_ptr< moris::Cell< moris::mtk::Vertex* > > > tNodesForTemplates( mCutIntegrationMesh->get_num_entities( EntityRank::ELEMENT, 0 ), nullptr );
    moris::Cell< std::shared_ptr< Node_Hierarchy_Template > >            tNHTemplate( mCutIntegrationMesh->get_num_entities( EntityRank::ELEMENT, 0 ), nullptr );

    // for each cell we need to select a node hier template
    // number of cells we are creating
    moris_index                              tNumNewCells  = 0;
    moris_index                              tNodesPerCell = 0;
    std::shared_ptr< moris::mtk::Cell_Info > tCellInfo     = nullptr;
    if ( mBackgroundMesh->get_spatial_dim() == 3 )
    {
        tNumNewCells = this->select_node_hier_3d_template(
            &tCellIndexIntersectedEdgeOrdinals,
            &tCellIndexIntersectedEdgeVertex,
            &tNodesForTemplates,
            &tNHTemplate );

        // constant parameters for 3d case
        tNodesPerCell = 4;
        moris::mtk::Cell_Info_Factory tFactory;
        tCellInfo = tFactory.create_cell_info_sp( CellTopology::TET4 );
    }

    else if ( mBackgroundMesh->get_spatial_dim() == 2 )
    {
        MORIS_ERROR( 0, "TODO" );
    }

    else
    {
        MORIS_ERROR( 0, "Only implemented in 2d and 3d." );
    }


    // allocate the algorithm data
    mNewCellToVertexConnectivity = moris::Cell< moris::Cell< moris::moris_index > >( tNumNewCells, moris::Cell< moris::moris_index >( tNodesPerCell ) );
    mNewCellChildMeshIndex       = moris::Cell< moris::moris_index >( tNumNewCells );
    mNewCellCellIndexToReplace   = moris::Cell< moris::moris_index >( tNumNewCells );
    mNewCellCellInfo             = moris::Cell< std::shared_ptr< moris::mtk::Cell_Info > >( tNumNewCells, tCellInfo );

    // std::cout<<"tNodesForTemplates.size() = "<<tNodesForTemplates.size()<<Std:

    moris_index tCurrentCellIndex = 0;
    for ( moris::uint iCell = 0; iCell < tNodesForTemplates.size(); iCell++ )
    {
        if ( tNodesForTemplates( iCell ) != nullptr )
        {
            // cell group membership
            moris_index tCellGroupMembershipIndex = mCutIntegrationMesh->get_ig_cell_group_memberships( (moris_index)iCell )( 0 );

            // get the template
            std::shared_ptr< Node_Hierarchy_Template > tTemplate = tNHTemplate( iCell );

            // vertices sorted for this template
            moris::Cell< moris::mtk::Vertex* >& tSortedVertices = *(tNodesForTemplates)( iCell );

            // moris::Matrix< moris::IndexMat > mCellToNodeOrdinal;
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
                for ( moris::uint iV = 0; iV < tTemplate->mCellToNodeOrdinal.n_cols(); iV++ )
                {
                    //
                    moris_index tVertexSortedOrd                            = tTemplate->mCellToNodeOrdinal( iTemplateCell, iV );
                    mNewCellToVertexConnectivity( tCurrentCellIndex )( iV ) = tSortedVertices( tVertexSortedOrd )->get_index();
                }
                tCurrentCellIndex++;
            }
        }
    }
}

void
Node_Hierarchy_Interface::determine_intersected_cell_information(
    std::shared_ptr< Edge_Based_Connectivity >                            aEdgeConnectivity,
    moris::Cell< moris_index >*                                           aIntersectedEdges,
    moris::Cell< std::shared_ptr< moris::Cell< moris_index > > >*         aCellIndexIntersectedEdgeOrdinals,
    moris::Cell< std::shared_ptr< moris::Cell< moris::mtk::Vertex* > > >* aCellIndexIntersectedEdgeVertex )
{
    // max new cell per intersected cell ~6
    moris::uint tNumEdgesPerCell = 6;// just use 6's since we operate on tets or tris. more conservative with 6 (slightly more memory)

    // iterate through intersected edges and figure out which side ordinals of an ig cell is intersected
    for ( moris::uint iEdge = 0; iEdge < aIntersectedEdges->size(); iEdge++ )
    {
        moris_index tEdgeIndex = ( *aIntersectedEdges )( iEdge );

        // vertex on the edge
        moris::mtk::Vertex* tVertexOnEdge = mCutIntegrationMesh->get_mtk_vertex_pointer( mDecompositionData->tNewNodeIndex( iEdge ) );

        for ( moris::uint iCell = 0; iCell < aEdgeConnectivity->mEdgeToCell( tEdgeIndex ).size(); iCell++ )
        {
            // integration cell with this edge
            moris::mtk::Cell* tCell = aEdgeConnectivity->mEdgeToCell( tEdgeIndex )( iCell );

            // integration cell edge ordinal
            moris_index tSideOrdinal = aEdgeConnectivity->mEdgeToCellEdgeOrdinal( tEdgeIndex )( iCell );

            if ( ( *aCellIndexIntersectedEdgeOrdinals )( tCell->get_index() ) == nullptr )
            {
                ( *aCellIndexIntersectedEdgeOrdinals )( tCell->get_index() ) = std::make_shared< moris::Cell< moris_index > >( 0 );
                ( *aCellIndexIntersectedEdgeOrdinals )( tCell->get_index() )->reserve( tNumEdgesPerCell );
                ( *aCellIndexIntersectedEdgeVertex )( tCell->get_index() ) = std::make_shared< moris::Cell< moris::mtk::Vertex* > >( 0 );
                ( *aCellIndexIntersectedEdgeVertex )( tCell->get_index() )->reserve( tNumEdgesPerCell );
            }

            ( *aCellIndexIntersectedEdgeOrdinals )( tCell->get_index() )->push_back( tSideOrdinal );
            ( *aCellIndexIntersectedEdgeVertex )( tCell->get_index() )->push_back( tVertexOnEdge );
        }
    }
}

moris_index
Node_Hierarchy_Interface::select_node_hier_3d_template(
    moris::Cell< std::shared_ptr< moris::Cell< moris_index > > >*         aCellIndexIntersectedEdgeOrdinals,
    moris::Cell< std::shared_ptr< moris::Cell< moris::mtk::Vertex* > > >* aCellIndexIntersectedEdgeVertex,
    moris::Cell< std::shared_ptr< moris::Cell< moris::mtk::Vertex* > > >* aNodesForTemplates,
    moris::Cell< std::shared_ptr< Node_Hierarchy_Template > >*            aNHTemplate )
{
    // tally up the total number of cells we are going to construct
    moris_index tNumNewIgCells = 0;

    mNumNewCells = 0;

    std::unordered_map< moris_index, std::shared_ptr< Node_Hierarchy_Template > > tLoadedTemplates;
    std::shared_ptr< Matrix< IdMat > >                                            tEdgeToVertexOrdinalMap = nullptr;

    Node_Hierachy_Template_Library tLibrary;

    // iterate through cells
    for ( moris::uint iCell = 0; iCell < aCellIndexIntersectedEdgeOrdinals->size(); iCell++ )
    {
        if ( ( *aCellIndexIntersectedEdgeOrdinals )( iCell ) != nullptr )
        {
            moris::mtk::Cell const* tIgCell = &mCutIntegrationMesh->get_mtk_cell( iCell );

            if ( tEdgeToVertexOrdinalMap == nullptr )
            {
                tEdgeToVertexOrdinalMap = std::make_shared< Matrix< IdMat > >( tIgCell->get_cell_info()->get_node_to_edge_map() );
            }
            moris_index tPermutationId = MORIS_INDEX_MAX;

            ( *aNodesForTemplates )( iCell ) = std::make_shared< moris::Cell< moris::mtk::Vertex* > >();
            ( *aNodesForTemplates )( iCell )->reserve( 8 );

            this->sort_nodes(
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

                tLoadedTemplates[tPermutationId] = tNewTemplate;
            }
            ( *aNHTemplate )( iCell ) = tLoadedTemplates[tPermutationId];


            tNumNewIgCells = tNumNewIgCells + ( *aNHTemplate )( iCell )->mNumCells;
            mNumNewCells   = mNumNewCells + ( *aNHTemplate )( iCell )->mNumCells - 1;
        }
    }
    return tNumNewIgCells;
}

void
Node_Hierarchy_Interface::sort_nodes(
    moris::mtk::Cell const*                               aIgCell,
    Matrix< IndexMat >*                                   aEdgeToVertexOrdinalMap,
    std::shared_ptr< moris::Cell< moris_index > >         aCellIndexIntersectedEdgeOrdinals,
    std::shared_ptr< moris::Cell< moris::mtk::Vertex* > > aCellIndexIntersectedEdgeVertex,
    moris_index&                                          aPermutation,
    std::shared_ptr< moris::Cell< moris::mtk::Vertex* > > aSortedNodeInds )
{
    // hier tet 4
    if ( mBackgroundMesh->get_spatial_dim() == 3 )
    {
        // initialize original index locations
        moris::Cell< moris_index > tIndices( aCellIndexIntersectedEdgeOrdinals->size() );
        std::iota( tIndices.data().begin(), tIndices.data().end(), 0 );

        // get ascending order vertices
        std::stable_sort( tIndices.data().begin(), tIndices.data().end(), [&]( std::size_t i1, std::size_t i2 ) { return ( *aCellIndexIntersectedEdgeVertex )( i1 )->get_id() < ( *aCellIndexIntersectedEdgeVertex )( i2 )->get_id(); } );


        Cell< moris::mtk::Vertex* > tVertices = aIgCell->get_vertex_pointers();


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

        else
        {
            MORIS_ERROR( 0, "Unhandled template 3d node hierarchy" );
        }
    }
}

void
Node_Hierachy_Template_Library::load_template(
    moris_index              aSpatialDim,
    moris_index              aTemplateId,
    Node_Hierarchy_Template* aNodeHierTemplate )
{
    if ( aSpatialDim == 3 )
    {

        this->load_3d_template( aTemplateId, aNodeHierTemplate );
    }
}

void
Node_Hierachy_Template_Library::load_3d_template(
    moris_index              aTemplateId,
    Node_Hierarchy_Template* aNodeHierTemplate )
{
    aNodeHierTemplate->mCellTopology = CellTopology::TET4;
    switch ( aTemplateId )
    {
    case ( 320 ):
    {
        // Permutation 320
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 5, 6 }, { 4, 1, 5, 6 }, { 1, 2, 5, 6 }, { 1, 3, 2, 6 } } );
        aNodeHierTemplate->mNumCells          = 4;
        break;
    }

    case ( 32 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 5, 6 }, { 4, 1, 5, 6 }, { 1, 2, 5, 6 }, { 1, 3, 2, 6 } } );
        aNodeHierTemplate->mNumCells          = 4;
        break;
    }


    case ( 203 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 5, 6 }, { 4, 1, 5, 6 }, { 1, 2, 5, 6 }, { 1, 3, 2, 6 } } );
        aNodeHierTemplate->mNumCells          = 4;
        break;
    }
    case ( 251 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 5, 6 }, { 4, 1, 5, 6 }, { 1, 2, 5, 6 }, { 1, 3, 2, 6 } } );
        aNodeHierTemplate->mNumCells          = 4;
        break;
    }

    case ( 512 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 5, 6 }, { 4, 1, 5, 6 }, { 1, 2, 5, 6 }, { 1, 3, 2, 6 } } );
        aNodeHierTemplate->mNumCells          = 4;
        break;
    }
    case ( 125 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 5, 6 }, { 4, 1, 5, 6 }, { 1, 2, 5, 6 }, { 1, 3, 2, 6 } } );
        aNodeHierTemplate->mNumCells          = 4;
        break;
    }
    case ( 140 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 5, 6 }, { 4, 1, 5, 6 }, { 1, 2, 5, 6 }, { 1, 3, 2, 6 } } );
        aNodeHierTemplate->mNumCells          = 4;
        break;
    }

    case ( 401 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 5, 6 }, { 4, 1, 5, 6 }, { 1, 2, 5, 6 }, { 1, 3, 2, 6 } } );
        aNodeHierTemplate->mNumCells          = 4;
        break;
    }
    case ( 14 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 5, 6 }, { 4, 1, 5, 6 }, { 1, 2, 5, 6 }, { 1, 3, 2, 6 } } );
        aNodeHierTemplate->mNumCells          = 4;
        break;
    }
    case ( 453 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 5, 6 }, { 4, 1, 5, 6 }, { 1, 2, 5, 6 }, { 1, 3, 2, 6 } } );
        aNodeHierTemplate->mNumCells          = 4;
        break;
    }
    case ( 534 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 5, 6 }, { 4, 1, 5, 6 }, { 1, 2, 5, 6 }, { 1, 3, 2, 6 } } );
        aNodeHierTemplate->mNumCells          = 4;
        break;
    }
    case ( 345 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 5, 6 }, { 4, 1, 5, 6 }, { 1, 2, 5, 6 }, { 1, 3, 2, 6 } } );
        aNodeHierTemplate->mNumCells          = 4;
        break;
    }

    case ( 230 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 5, 6 }, { 1, 4, 5, 6 }, { 2, 1, 5, 6 }, { 3, 1, 2, 6 } } );
        aNodeHierTemplate->mNumCells          = 4;
        break;
    }
    case ( 302 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 5, 6 }, { 1, 4, 5, 6 }, { 2, 1, 5, 6 }, { 3, 1, 2, 6 } } );
        aNodeHierTemplate->mNumCells          = 4;
        break;
    }
    case ( 23 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 5, 6 }, { 1, 4, 5, 6 }, { 2, 1, 5, 6 }, { 3, 1, 2, 6 } } );
        aNodeHierTemplate->mNumCells          = 4;
        break;
    }
    case ( 521 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 5, 6 }, { 1, 4, 5, 6 }, { 2, 1, 5, 6 }, { 3, 1, 2, 6 } } );
        aNodeHierTemplate->mNumCells          = 4;
        break;
    }
    case ( 152 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 5, 6 }, { 1, 4, 5, 6 }, { 2, 1, 5, 6 }, { 3, 1, 2, 6 } } );
        aNodeHierTemplate->mNumCells          = 4;
        break;
    }
    case ( 215 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 5, 6 }, { 1, 4, 5, 6 }, { 2, 1, 5, 6 }, { 3, 1, 2, 6 } } );
        aNodeHierTemplate->mNumCells          = 4;
        break;
    }
    case ( 410 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 5, 6 }, { 1, 4, 5, 6 }, { 2, 1, 5, 6 }, { 3, 1, 2, 6 } } );
        aNodeHierTemplate->mNumCells          = 4;
        break;
    }
    case ( 41 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 5, 6 }, { 1, 4, 5, 6 }, { 2, 1, 5, 6 }, { 3, 1, 2, 6 } } );
        aNodeHierTemplate->mNumCells          = 4;
        break;
    }
    case ( 104 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 5, 6 }, { 1, 4, 5, 6 }, { 2, 1, 5, 6 }, { 3, 1, 2, 6 } } );
        aNodeHierTemplate->mNumCells          = 4;
        break;
    }
    case ( 543 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 5, 6 }, { 1, 4, 5, 6 }, { 2, 1, 5, 6 }, { 3, 1, 2, 6 } } );
        aNodeHierTemplate->mNumCells          = 4;
        break;
    }
    case ( 354 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 5, 6 }, { 1, 4, 5, 6 }, { 2, 1, 5, 6 }, { 3, 1, 2, 6 } } );
        aNodeHierTemplate->mNumCells          = 4;
        break;
    }
    case ( 435 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 5, 6 }, { 1, 4, 5, 6 }, { 2, 1, 5, 6 }, { 3, 1, 2, 6 } } );
        aNodeHierTemplate->mNumCells          = 4;
        break;
    }
    case ( 5420 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 3, 7, 6 }, { 4, 3, 6, 7 }, { 4, 3, 7, 5 }, { 1, 2, 7, 5 }, { 2, 4, 7, 5 }, { 4, 2, 7, 6 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 5240 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 3, 7, 6 }, { 4, 3, 6, 7 }, { 4, 3, 7, 5 }, { 1, 2, 7, 5 }, { 2, 4, 7, 5 }, { 4, 2, 7, 6 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 425 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 3, 7, 6 }, { 4, 3, 6, 7 }, { 4, 3, 7, 5 }, { 1, 2, 7, 5 }, { 2, 4, 7, 5 }, { 4, 2, 7, 6 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 3501 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 3, 7, 6 }, { 4, 3, 6, 7 }, { 4, 3, 7, 5 }, { 1, 2, 7, 5 }, { 2, 4, 7, 5 }, { 4, 2, 7, 6 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 1503 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 3, 7, 6 }, { 4, 3, 6, 7 }, { 4, 3, 7, 5 }, { 1, 2, 7, 5 }, { 2, 4, 7, 5 }, { 4, 2, 7, 6 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 3051 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 3, 7, 6 }, { 4, 3, 6, 7 }, { 4, 3, 7, 5 }, { 1, 2, 7, 5 }, { 2, 4, 7, 5 }, { 4, 2, 7, 6 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 1053 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 3, 7, 6 }, { 4, 3, 6, 7 }, { 4, 3, 7, 5 }, { 1, 2, 7, 5 }, { 2, 4, 7, 5 }, { 4, 2, 7, 6 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 4312 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 3, 7, 6 }, { 4, 3, 6, 7 }, { 4, 3, 7, 5 }, { 1, 2, 7, 5 }, { 2, 4, 7, 5 }, { 4, 2, 7, 6 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 2314 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 3, 7, 6 }, { 4, 3, 6, 7 }, { 4, 3, 7, 5 }, { 1, 2, 7, 5 }, { 2, 4, 7, 5 }, { 4, 2, 7, 6 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 4132 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 3, 7, 6 }, { 4, 3, 6, 7 }, { 4, 3, 7, 5 }, { 1, 2, 7, 5 }, { 2, 4, 7, 5 }, { 4, 2, 7, 6 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 2134 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 3, 7, 6 }, { 4, 3, 6, 7 }, { 4, 3, 7, 5 }, { 1, 2, 7, 5 }, { 2, 4, 7, 5 }, { 4, 2, 7, 6 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 245 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 3, 7, 6 }, { 4, 3, 6, 7 }, { 4, 3, 7, 5 }, { 1, 2, 7, 5 }, { 2, 4, 7, 5 }, { 4, 2, 7, 6 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 4502 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 3, 0, 7, 6 }, { 3, 4, 6, 7 }, { 3, 4, 7, 5 }, { 2, 1, 7, 5 }, { 4, 2, 7, 5 }, { 2, 4, 7, 6 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 4052 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 3, 0, 7, 6 }, { 3, 4, 6, 7 }, { 3, 4, 7, 5 }, { 2, 1, 7, 5 }, { 4, 2, 7, 5 }, { 2, 4, 7, 6 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 1243 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 3, 0, 7, 6 }, { 3, 4, 6, 7 }, { 3, 4, 7, 5 }, { 2, 1, 7, 5 }, { 4, 2, 7, 5 }, { 2, 4, 7, 6 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 2504 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 3, 0, 7, 6 }, { 3, 4, 6, 7 }, { 3, 4, 7, 5 }, { 2, 1, 7, 5 }, { 4, 2, 7, 5 }, { 2, 4, 7, 6 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 2054 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 3, 0, 7, 6 }, { 3, 4, 6, 7 }, { 3, 4, 7, 5 }, { 2, 1, 7, 5 }, { 4, 2, 7, 5 }, { 2, 4, 7, 6 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 5310 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 3, 0, 7, 6 }, { 3, 4, 6, 7 }, { 3, 4, 7, 5 }, { 2, 1, 7, 5 }, { 4, 2, 7, 5 }, { 2, 4, 7, 6 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 5130 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 3, 0, 7, 6 }, { 3, 4, 6, 7 }, { 3, 4, 7, 5 }, { 2, 1, 7, 5 }, { 4, 2, 7, 5 }, { 2, 4, 7, 6 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 135 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 3, 0, 7, 6 }, { 3, 4, 6, 7 }, { 3, 4, 7, 5 }, { 2, 1, 7, 5 }, { 4, 2, 7, 5 }, { 2, 4, 7, 6 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 315 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 3, 0, 7, 6 }, { 3, 4, 6, 7 }, { 3, 4, 7, 5 }, { 2, 1, 7, 5 }, { 4, 2, 7, 5 }, { 2, 4, 7, 6 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 3421 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 3, 0, 7, 6 }, { 3, 4, 6, 7 }, { 3, 4, 7, 5 }, { 2, 1, 7, 5 }, { 4, 2, 7, 5 }, { 2, 4, 7, 6 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 3241 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 3, 0, 7, 6 }, { 3, 4, 6, 7 }, { 3, 4, 7, 5 }, { 2, 1, 7, 5 }, { 4, 2, 7, 5 }, { 2, 4, 7, 6 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 1423 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 3, 0, 7, 6 }, { 3, 4, 6, 7 }, { 3, 4, 7, 5 }, { 2, 1, 7, 5 }, { 4, 2, 7, 5 }, { 2, 4, 7, 6 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 4250 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 6, 7 }, { 0, 3, 6, 7 }, { 3, 5, 6, 7 }, { 5, 1, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 2450 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 6, 7 }, { 0, 3, 6, 7 }, { 3, 5, 6, 7 }, { 5, 1, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 4205 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 6, 7 }, { 0, 3, 6, 7 }, { 3, 5, 6, 7 }, { 5, 1, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 2405 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 6, 7 }, { 0, 3, 6, 7 }, { 3, 5, 6, 7 }, { 5, 1, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 5031 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 6, 7 }, { 0, 3, 6, 7 }, { 3, 5, 6, 7 }, { 5, 1, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 5013 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 6, 7 }, { 0, 3, 6, 7 }, { 3, 5, 6, 7 }, { 5, 1, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 531 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 6, 7 }, { 0, 3, 6, 7 }, { 3, 5, 6, 7 }, { 5, 1, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 513 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 6, 7 }, { 0, 3, 6, 7 }, { 3, 5, 6, 7 }, { 5, 1, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 3124 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 6, 7 }, { 0, 3, 6, 7 }, { 3, 5, 6, 7 }, { 5, 1, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 1342 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 6, 7 }, { 0, 3, 6, 7 }, { 3, 5, 6, 7 }, { 5, 1, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 1324 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 6, 7 }, { 0, 3, 6, 7 }, { 3, 5, 6, 7 }, { 5, 1, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 3142 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 6, 7 }, { 0, 3, 6, 7 }, { 3, 5, 6, 7 }, { 5, 1, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 5024 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 6, 7 }, { 3, 0, 6, 7 }, { 5, 3, 6, 7 }, { 1, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 524 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 6, 7 }, { 3, 0, 6, 7 }, { 5, 3, 6, 7 }, { 1, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 5042 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 6, 7 }, { 3, 0, 6, 7 }, { 5, 3, 6, 7 }, { 1, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 542 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 6, 7 }, { 3, 0, 6, 7 }, { 5, 3, 6, 7 }, { 1, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 3150 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 6, 7 }, { 3, 0, 6, 7 }, { 5, 3, 6, 7 }, { 1, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 1350 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 6, 7 }, { 3, 0, 6, 7 }, { 5, 3, 6, 7 }, { 1, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 3105 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 6, 7 }, { 3, 0, 6, 7 }, { 5, 3, 6, 7 }, { 1, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 1305 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 6, 7 }, { 3, 0, 6, 7 }, { 5, 3, 6, 7 }, { 1, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 4231 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 6, 7 }, { 3, 0, 6, 7 }, { 5, 3, 6, 7 }, { 1, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 2431 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 6, 7 }, { 3, 0, 6, 7 }, { 5, 3, 6, 7 }, { 1, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 4213 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 6, 7 }, { 3, 0, 6, 7 }, { 5, 3, 6, 7 }, { 1, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 2413 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 6, 7 }, { 3, 0, 6, 7 }, { 5, 3, 6, 7 }, { 1, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 4520 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 5, 7 }, { 0, 3, 5, 7 }, { 5, 3, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 5, 7 }, { 2, 5, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 2540 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 5, 7 }, { 0, 3, 5, 7 }, { 5, 3, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 5, 7 }, { 2, 5, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 4025 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 5, 7 }, { 0, 3, 5, 7 }, { 5, 3, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 5, 7 }, { 2, 5, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 2045 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 5, 7 }, { 0, 3, 5, 7 }, { 5, 3, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 5, 7 }, { 2, 5, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 5301 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 5, 7 }, { 0, 3, 5, 7 }, { 5, 3, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 5, 7 }, { 2, 5, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 5103 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 5, 7 }, { 0, 3, 5, 7 }, { 5, 3, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 5, 7 }, { 2, 5, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 351 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 5, 7 }, { 0, 3, 5, 7 }, { 5, 3, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 5, 7 }, { 2, 5, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 153 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 5, 7 }, { 0, 3, 5, 7 }, { 5, 3, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 5, 7 }, { 2, 5, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 3412 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 5, 7 }, { 0, 3, 5, 7 }, { 5, 3, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 5, 7 }, { 2, 5, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 3214 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 5, 7 }, { 0, 3, 5, 7 }, { 5, 3, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 5, 7 }, { 2, 5, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 1432 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 5, 7 }, { 0, 3, 5, 7 }, { 5, 3, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 5, 7 }, { 2, 5, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 1234 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 4, 0, 5, 7 }, { 0, 3, 5, 7 }, { 5, 3, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 5, 7 }, { 2, 5, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 5402 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 5, 7 }, { 3, 0, 5, 7 }, { 3, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 5, 7 }, { 5, 2, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 452 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 5, 7 }, { 3, 0, 5, 7 }, { 3, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 5, 7 }, { 5, 2, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 5204 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 5, 7 }, { 3, 0, 5, 7 }, { 3, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 5, 7 }, { 5, 2, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 254 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 5, 7 }, { 3, 0, 5, 7 }, { 3, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 5, 7 }, { 5, 2, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 3510 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 5, 7 }, { 3, 0, 5, 7 }, { 3, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 5, 7 }, { 5, 2, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 1530 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 5, 7 }, { 3, 0, 5, 7 }, { 3, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 5, 7 }, { 5, 2, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 3015 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 5, 7 }, { 3, 0, 5, 7 }, { 3, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 5, 7 }, { 5, 2, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 1035 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 5, 7 }, { 3, 0, 5, 7 }, { 3, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 5, 7 }, { 5, 2, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 4321 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 5, 7 }, { 3, 0, 5, 7 }, { 3, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 5, 7 }, { 5, 2, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 2341 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 5, 7 }, { 3, 0, 5, 7 }, { 3, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 5, 7 }, { 5, 2, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 4123 ):
    {

        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 5, 7 }, { 3, 0, 5, 7 }, { 3, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 5, 7 }, { 5, 2, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 2143 ):
    {
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 5, 7 }, { 3, 0, 5, 7 }, { 3, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 5, 7 }, { 5, 2, 6, 7 } } );
        aNodeHierTemplate->mNumCells          = 6;
        break;
    }
    case ( 10000 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 2, 3 }, { 4, 1, 2, 3 } } );
        aNodeHierTemplate->mNumCells          = 2;
        break;
    case ( 10001 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 1, 4, 3 }, { 0, 4, 2, 3 } } );
        aNodeHierTemplate->mNumCells          = 2;
        break;
    case ( 10002 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 1, 4, 3 }, { 4, 1, 2, 3 } } );
        aNodeHierTemplate->mNumCells          = 2;
        break;
    case ( 10003 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 1, 2, 4 }, { 4, 1, 2, 3 } } );
        aNodeHierTemplate->mNumCells          = 2;
        break;
    case ( 10004 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 2, 3 }, { 0, 1, 2, 4 } } );
        aNodeHierTemplate->mNumCells          = 2;
        break;
    case ( 10005 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 1, 4, 3 }, { 0, 1, 2, 4 } } );
        aNodeHierTemplate->mNumCells          = 2;
        break;


    case ( 10250 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 1, 4, 5 }, { 0, 4, 2, 5 }, { 0, 5, 2, 3 } } );
        aNodeHierTemplate->mNumCells          = 3;
        break;
    case ( 10520 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 1, 5, 4 }, { 0, 4, 5, 3 }, { 0, 5, 2, 3 } } );
        aNodeHierTemplate->mNumCells          = 3;
        break;
    case ( 10260 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 1, 4, 5 }, { 0, 1, 5, 3 }, { 0, 4, 2, 5 } } );
        aNodeHierTemplate->mNumCells          = 3;
        break;
    case ( 10620 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 1, 5, 3 }, { 0, 5, 4, 3 }, { 0, 5, 2, 4 } } );
        aNodeHierTemplate->mNumCells          = 3;
        break;
    case ( 10560 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 1, 5, 4 }, { 0, 1, 2, 5 }, { 0, 4, 5, 3 } } );
        aNodeHierTemplate->mNumCells          = 3;
        break;
    case ( 10650 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 1, 2, 5 }, { 0, 5, 2, 4 }, { 0, 5, 4, 3 } } );
        aNodeHierTemplate->mNumCells          = 3;
        break;
    case ( 10640 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 1, 2, 5 }, { 5, 1, 4, 3 }, { 1, 2, 5, 4 } } );
        aNodeHierTemplate->mNumCells          = 3;
        break;
    case ( 10460 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 1, 2, 5 }, { 0, 1, 5, 4 }, { 4, 1, 5, 3 } } );
        aNodeHierTemplate->mNumCells          = 3;
        break;
    case ( 10230 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 5, 4, 2, 3 }, { 1, 4, 5, 3 }, { 0, 1, 5, 3 } } );
        aNodeHierTemplate->mNumCells          = 3;
        break;
    case ( 10320 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 5, 2, 4, 3 }, { 0, 5, 4, 3 }, { 0, 1, 5, 3 } } );
        aNodeHierTemplate->mNumCells          = 3;
        break;
    case ( 10130 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 5, 1, 2, 3 }, { 5, 4, 1, 3 }, { 0, 4, 5, 3 } } );
        aNodeHierTemplate->mNumCells          = 3;
        break;
    case ( 10310 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 5, 1, 2, 3 }, { 4, 5, 2, 3 }, { 0, 5, 4, 3 } } );
        aNodeHierTemplate->mNumCells          = 3;
        break;
    case ( 10120 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 5, 2, 3 }, { 0, 4, 5, 3 }, { 4, 1, 5, 3 } } );
        aNodeHierTemplate->mNumCells          = 3;
        break;
    case ( 10210 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 5, 2, 3 }, { 5, 4, 2, 3 }, { 5, 1, 4, 3 } } );
        aNodeHierTemplate->mNumCells          = 3;
        break;
    case ( 10140 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 2, 5 }, { 4, 1, 2, 5 }, { 5, 1, 2, 3 } } );
        aNodeHierTemplate->mNumCells          = 3;
        break;
    case ( 10410 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 5, 2, 4 }, { 4, 5, 2, 3 }, { 5, 1, 2, 3 } } );
        aNodeHierTemplate->mNumCells          = 3;
        break;
    case ( 10150 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 4, 2, 5 }, { 4, 1, 2, 5 }, { 0, 5, 2, 3 } } );
        aNodeHierTemplate->mNumCells          = 3;
        break;
    case ( 10510 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 5, 1, 2, 4 }, { 5, 4, 2, 3 }, { 0, 5, 2, 3 } } );
        aNodeHierTemplate->mNumCells          = 3;
        break;
    case ( 10450 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 5, 2, 4 }, { 0, 1, 2, 5 }, { 4, 5, 2, 3 } } );
        aNodeHierTemplate->mNumCells          = 3;
        break;
    case ( 10540 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 1, 2, 5 }, { 5, 1, 2, 4 }, { 5, 4, 2, 3 } } );
        aNodeHierTemplate->mNumCells          = 3;
        break;
    case ( 10360 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 1, 5, 3 }, { 0, 1, 4, 5 }, { 1, 2, 4, 5 } } );
        aNodeHierTemplate->mNumCells          = 3;
        break;
    case ( 10630 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 1, 5, 3 }, { 1, 4, 5, 3 }, { 1, 2, 5, 4 } } );
        aNodeHierTemplate->mNumCells          = 3;
        break;
    case ( 10340 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 1, 4, 5 }, { 4, 1, 2, 5 }, { 1, 2, 5, 3 } } );
        aNodeHierTemplate->mNumCells          = 3;
        break;
    case ( 10430 ):
        aNodeHierTemplate->mCellToNodeOrdinal = moris::Matrix< moris::IndexMat >( { { 0, 1, 5, 4 }, { 4, 1, 5, 3 }, { 1, 2, 5, 3 } } );
        aNodeHierTemplate->mNumCells          = 3;
        break;
    default:
    {
        std::cout << "WARNING INVALID TEMPLATE ID: " << aTemplateId << std::endl;
        break;
    }
    }
}
}// namespace xtk
