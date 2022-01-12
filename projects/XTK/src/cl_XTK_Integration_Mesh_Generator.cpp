#include "cl_XTK_Integration_Mesh_Generator.hpp"
#include "cl_XTK_Decomposition_Algorithm_Factory.hpp"
#include "cl_XTK_Decomposition_Algorithm.hpp"
#include "fn_determine_cell_topology.hpp"

using namespace moris;
namespace xtk
{
Integration_Mesh_Generator::Integration_Mesh_Generator()
{
}

Integration_Mesh_Generator::Integration_Mesh_Generator(
    xtk::Model*                      aXTKModelPtr,
    Cell< enum Subdivision_Method >  aMethods,
    moris::Matrix< moris::IndexMat > aActiveGeometries ) :
    mXTKModel( aXTKModelPtr ),
    mGeometryEngine( mXTKModel->get_geom_engine() ), mActiveGeometries( aActiveGeometries ), mSubdivisionMethods( aMethods )
{
}
// ----------------------------------------------------------------------------------

Integration_Mesh_Generator::~Integration_Mesh_Generator()
{
}
// ----------------------------------------------------------------------------------

std::shared_ptr< Cut_Integration_Mesh >
Integration_Mesh_Generator::perform()
{
    Tracer tTracer( "XTK", "Integration_Mesh_Generator", "perform", mXTKModel->mVerboseLevel, 0 );

    // data structure to pass things around
    Integration_Mesh_Generation_Data tGenerationData;

    // pointer to the background mesh
    moris::mtk::Mesh* tBackgroundMesh = &mXTKModel->get_background_mesh();

    // the cut integration mesh
    std::shared_ptr< Cut_Integration_Mesh > tCutIntegrationMesh = std::make_shared< Cut_Integration_Mesh >( tBackgroundMesh, mXTKModel );

    mXTKModel->set_cut_ig_mesh(tCutIntegrationMesh);

    // Allocate a child mesh for each background cell
    this->allocate_child_meshes( tGenerationData, tCutIntegrationMesh.get(), tBackgroundMesh );

    // figure out which background cells are intersected and by which geometry they are intersected
    this->determine_intersected_background_cells( tGenerationData, tCutIntegrationMesh.get(), tBackgroundMesh );

    // verify levels of intersected background cells
    this->check_intersected_background_cell_levels(tGenerationData,tCutIntegrationMesh.get(),tBackgroundMesh);

    // make mSameLevelChildMeshes flag global consistent

    if(!tCutIntegrationMesh->mSameLevelChildMeshes)
    {
        return tCutIntegrationMesh;
    }

    // iterate through the subdivision methods
    for ( moris::uint iSubMethod = 0; iSubMethod < mSubdivisionMethods.size(); iSubMethod++ )
    {
        // create the subdivision routine with a factory
        std::shared_ptr< Decomposition_Algorithm > tDecompAlg =
                create_decomposition_algorithm( mSubdivisionMethods( iSubMethod ), mXTKModel->get_parameter_list() );

        // perform the decomposition
        Decomposition_Data tDecompositionData;

        tDecompAlg->perform( &tGenerationData, &tDecompositionData, tCutIntegrationMesh.get(), tBackgroundMesh, this );
    }

    tCutIntegrationMesh->finalize_cut_mesh_construction();

    // set the bulk phase of each cell
    this->compute_ig_cell_bulk_phase( tCutIntegrationMesh.get() );

    // pick out the cell groups
    moris::Cell< std::shared_ptr< IG_Cell_Group > >& tActiveIgCellGroups = tCutIntegrationMesh->get_all_cell_groups();

    moris::Cell< moris::mtk::Cell* > tActiveIgCells;
    extract_cells_from_cell_groups( tActiveIgCellGroups, tActiveIgCells );

    // create facet connectivity in the mesh
    std::shared_ptr< Facet_Based_Connectivity > tFaceConnectivity = std::make_shared< Facet_Based_Connectivity >();
    this->create_facet_from_element_to_node( tActiveIgCells, tFaceConnectivity );
    tCutIntegrationMesh->set_face_connectivity( tFaceConnectivity );

    // create facet ancestry
    moris::Cell< moris::mtk::Cell* > tBGCellForFacet;
    this->select_background_cell_for_facet( tFaceConnectivity, tCutIntegrationMesh.get(), tBGCellForFacet );

    std::shared_ptr< Facet_Based_Ancestry > tFacetAncestry = std::make_shared< Facet_Based_Ancestry >();
    this->deduce_facet_ancestry( tCutIntegrationMesh.get(), tBackgroundMesh, tFaceConnectivity, tBGCellForFacet, tFacetAncestry );
    tCutIntegrationMesh->set_face_ancestry( tFacetAncestry );

    // compute the facets attached to a given bg facet. (useful to deduce side sets later and to construct the subphase neighborhood
    // for the enrichment strategy)
    moris::Cell< std::shared_ptr< moris::Cell< moris::moris_index > > > tBgFacetToChildFacet;
    this->compute_bg_facet_to_child_facet_connectivity( tCutIntegrationMesh.get(), tBackgroundMesh, tFaceConnectivity, tFacetAncestry, tBgFacetToChildFacet );
    tCutIntegrationMesh->set_background_facet_to_child_facet_connectivity( tBgFacetToChildFacet );

    // create element to element connectivity
    std::shared_ptr< Cell_Neighborhood_Connectivity > tNeighborhood = std::make_shared< Cell_Neighborhood_Connectivity >();
    this->generate_cell_neighborhood( tActiveIgCells, tFaceConnectivity, tNeighborhood );

    // figure out the interface facets
    moris::Cell< moris_index > tInterfaces;
    this->deduce_interfaces( tCutIntegrationMesh.get(), tFaceConnectivity, tInterfaces );
    tCutIntegrationMesh->set_interface_facets( tInterfaces );

    // determine which bulk phase these interfaces are between
    // size: num_bulk_phase x num_bulk_phase
    moris::Cell< moris::Cell< std::shared_ptr< IG_Cell_Side_Group > > > tInterfaceByBulkPhase;
    this->construct_bulk_phase_to_bulk_phase_interface( tCutIntegrationMesh.get(), tInterfaces, tInterfaceByBulkPhase );

    moris::Cell< moris::Cell< std::shared_ptr< IG_Cell_Double_Side_Group > > > tDoubleSidedInterface;
    this->construct_bulk_phase_to_bulk_phase_dbl_side_interface( tCutIntegrationMesh.get(), tInterfaces, tDoubleSidedInterface );
    tCutIntegrationMesh->set_bulk_phase_to_bulk_phase_dbl_side_interface( tDoubleSidedInterface );

    // construct interface_sets
    this->construct_interface_sets( tCutIntegrationMesh.get(), tInterfaceByBulkPhase );

    // print diagonstic information if requested
    if ( mXTKModel->mDiagnostics )
    {
        std::string tCellDiagFile    = mXTKModel->get_diagnostic_file_name( std::string( "Cells" ) );
        std::string tVertDiagFile    = mXTKModel->get_diagnostic_file_name( std::string( "Vertex" ) );
        std::string tGroupDiagFile   = mXTKModel->get_diagnostic_file_name( std::string( "Groups" ) );
        std::string tGENVertDiagFile = mXTKModel->get_diagnostic_file_name( std::string( "Vertex_GEN" ) );
        tCutIntegrationMesh->print_cells( false, tCellDiagFile );
        tCutIntegrationMesh->print_vertices( false, tVertDiagFile );
        tCutIntegrationMesh->print_groupings( tGroupDiagFile );
        mXTKModel->get_geom_engine()->print_gen_vertices( tGENVertDiagFile, tCutIntegrationMesh.get() );
    }

    // identify and construct the subphase groups
    this->identify_and_construct_subphases( &tGenerationData, tCutIntegrationMesh.get(), tBackgroundMesh, tNeighborhood );

    // construct subphase neighborhood
    std::shared_ptr< Subphase_Neighborhood_Connectivity > tSubphaseNeighborhood = std::make_shared< Subphase_Neighborhood_Connectivity >();
    this->construct_subphase_neighborhood( tCutIntegrationMesh.get(), tBackgroundMesh, tFaceConnectivity, &tBgFacetToChildFacet, tSubphaseNeighborhood );
    tCutIntegrationMesh->set_subphase_neighborhood( tSubphaseNeighborhood );

    // construct the bulk phase blocks
    moris::Cell< std::shared_ptr< IG_Cell_Group > > tBulkPhaseCellGroups;
    this->construct_bulk_phase_cell_groups( tCutIntegrationMesh.get(), tBulkPhaseCellGroups );

    // check if order elevation has been requested
    if ( this->get_ig_mesh_order() > 1 )
    {
        // get the order elevation template
        enum Subdivision_Method tOrderElevationMethod = this->determine_order_elevation_template();

        // create the subdivision routine with a factory
        std::shared_ptr< Decomposition_Algorithm > tElevateOrderAlg = 
            create_decomposition_algorithm( tOrderElevationMethod, mXTKModel->get_parameter_list() );

        // perform the decomposition
        Decomposition_Data tDecompositionData;
        tElevateOrderAlg->perform( &tGenerationData, &tDecompositionData, tCutIntegrationMesh.get(), tBackgroundMesh, this );
    }

    this->construct_bulk_phase_blocks( tCutIntegrationMesh.get(), tBulkPhaseCellGroups );

    // output cut IG mesh for debugging
    if ( mOutputCutIgMesh )
    {
        tCutIntegrationMesh->write_mesh("./","cut_ig_mesh.exo");
    }
    
    return tCutIntegrationMesh;
}

// ----------------------------------------------------------------------------------

moris::Matrix< moris::IndexMat > const*
Integration_Mesh_Generator::get_active_geometries()
{
    return &mActiveGeometries;
}

moris::ge::Geometry_Engine*
Integration_Mesh_Generator::get_geom_engine()
{
    return mGeometryEngine;
}

// ----------------------------------------------------------------------------------

uint 
Integration_Mesh_Generator::get_spatial_dim()
{
    return mXTKModel->get_spatial_dim();
}

// ----------------------------------------------------------------------------------

uint 
Integration_Mesh_Generator::get_ig_mesh_order()
{
    return this->mXTKModel->ig_element_order();
}

// ----------------------------------------------------------------------------------

enum Subdivision_Method
Integration_Mesh_Generator::determine_order_elevation_template()
{
    if( this->get_spatial_dim() == 2 )
    {
        switch ( this->get_ig_mesh_order() )
        {
            case 2:
            {
                return Subdivision_Method::P_ELEVATE_ORDER_TRI3_TRI6;
                break;
            }

            case 3:
            {
                MORIS_ERROR( false, "Integration_Mesh_Generator::determine_order_elevation_template() - elevate order to TRI10 not implemented yet" );
                return Subdivision_Method::P_ELEVATE_ORDER_TRI3_TRI10;
                break;
            }
            
            default:
            {
                MORIS_ERROR( false, "Integration_Mesh_Generator::determine_order_elevation_template() - integration element order not suported. max is 3" );
                return Subdivision_Method::NO_METHOD;
                break;
            }
        }
    }

    else if( this->get_spatial_dim() == 3 )
    {
        switch ( this->get_ig_mesh_order() )
        {
            case 2:
            {
                return Subdivision_Method::P_ELEVATE_ORDER_TET4_TET10;
                break;
            }

            case 3:
            {
                MORIS_ERROR( false, "Integration_Mesh_Generator::determine_order_elevation_template() - elevate order to TET20 not implemented yet" );
                return Subdivision_Method::P_ELEVATE_ORDER_TET4_TET20;
                break;
            }
            
            default:
            {
                MORIS_ERROR( false, "Integration_Mesh_Generator::determine_order_elevation_template() - integration element order not suported. max is 3" );
                return Subdivision_Method::NO_METHOD;
                break;
            }
        }
    }
    
    else
    {
        MORIS_ERROR( false, "Integration_Mesh_Generator::determine_order_elevation_template() - spatial dim must be 2 or 3." );
        return Subdivision_Method::NO_METHOD;
    }
}

// ----------------------------------------------------------------------------------

void
Integration_Mesh_Generator::check_intersected_background_cell_levels(
    Integration_Mesh_Generation_Data& aMeshGenerationData,
    Cut_Integration_Mesh*             aCutIntegrationMesh,
    moris::mtk::Mesh*                 aBackgroundMesh)
{
    Tracer tTracer( "XTK", "Integration_Mesh_Generator", "check_intersected_background_cell_levels " ,mXTKModel->mVerboseLevel, 1  );

    bool tFlag = true;

    moris_index tReferenceLevel = 0;
    if(aMeshGenerationData.mAllIntersectedBgCellInds.size() > 0)
    {
        tReferenceLevel = aBackgroundMesh->get_mtk_cell(aMeshGenerationData.mAllIntersectedBgCellInds(0)).get_level();
    }
    
    for (size_t iBgCellIndex = 1; iBgCellIndex < aMeshGenerationData.mAllIntersectedBgCellInds.size(); iBgCellIndex++)
    {
        moris_index tLevel = aBackgroundMesh->get_mtk_cell(aMeshGenerationData.mAllIntersectedBgCellInds(iBgCellIndex)).get_level();
        if(tReferenceLevel != tLevel)
        {
            tFlag = false;
        }
    }

    aCutIntegrationMesh->mSameLevelChildMeshes = all_land(tFlag);
}

// ----------------------------------------------------------------------------------

bool
Integration_Mesh_Generator::determine_intersected_background_cells(
    Integration_Mesh_Generation_Data& aMeshGenerationData,
    Cut_Integration_Mesh*             aCutIntegrationMesh,
    moris::mtk::Mesh*                 aBackgroundMesh )
{
    Tracer tTracer( "XTK", "Integration_Mesh_Generator", "Determine intersected background cells" ,mXTKModel->mVerboseLevel, 1  );
    uint   tNumGeometries = mActiveGeometries.numel();

    uint tNumCells = aBackgroundMesh->get_num_elems();

    aMeshGenerationData.mIntersectedBackgroundCellIndex.resize( tNumGeometries );

    aMeshGenerationData.mAllIntersectedBgCellInds.reserve( tNumCells );

    // Initialize geometric query
    Geometric_Query_XTK tGeometricQuery;

    // say I am just interested in a yes or no answer
    tGeometricQuery.set_query_type( moris::ge::Query_Type::INTERSECTION_NO_LOCATION );

    // large coord matrix that I want to keep in scope for a long time avoid copying coordinate all the time.
    tGeometricQuery.set_coordinates_matrix( &aCutIntegrationMesh->mVertexCoordinates );

    tGeometricQuery.set_query_entity_rank( EntityRank::ELEMENT );

    // reserve memory for list of indices of intersected background cells
    // size estimate: tNumCells / tNumGeometries
    for ( moris::size_t iGeom = 0; iGeom < tNumGeometries; iGeom++ )
    {
        aMeshGenerationData.mIntersectedBackgroundCellIndex( iGeom ).reserve ( tNumCells / tNumGeometries );
    }

    // iterate through all cells
    for ( moris::uint iCell = 0; iCell < tNumCells; iCell++ )
    {
        // setup geometric query with this current cell information
        tGeometricQuery.set_parent_cell( &aBackgroundMesh->get_mtk_cell( (moris_index)iCell ) );
        tGeometricQuery.set_query_cell( &aBackgroundMesh->get_mtk_cell( (moris_index)iCell ) );

        // iterate through all geometries for current cell
        for ( moris::size_t iGeom = 0; iGeom < tNumGeometries; iGeom++ )
        {
            // current index for this geometry
            moris_index tGeometryIndex = mActiveGeometries( iGeom );

            // tell the query which geometric index we are working on
            tGeometricQuery.set_geometric_index( tGeometryIndex );

            if ( mXTKModel->get_geom_engine()->geometric_query( &tGeometricQuery ) )
            {
                // add background cell to the list for iGEOM
                aMeshGenerationData.mIntersectedBackgroundCellIndex( iGeom ).push_back( iCell );

                aMeshGenerationData.mNumChildMeshes++;

                aMeshGenerationData.mAllIntersectedBgCellInds.push_back( iCell );
            }
            // add to the global list of intersected cells,
            else if( mXTKModel->triangulate_all() )
            {
                aMeshGenerationData.mAllIntersectedBgCellInds.push_back( iCell );
            }
        }
    }

    // remove the excess space
    shrink_to_fit_all( aMeshGenerationData.mIntersectedBackgroundCellIndex );

    unique( aMeshGenerationData.mAllIntersectedBgCellInds );

    return true;
}

// ----------------------------------------------------------------------------------

void
Integration_Mesh_Generator::commit_new_ig_cells_to_cut_mesh(
    Integration_Mesh_Generation_Data* aMeshGenerationData,
    Decomposition_Data*               aDecompositionData,
    Cut_Integration_Mesh*             aCutIntegrationMesh,
    moris::mtk::Mesh*                 aBackgroundMesh,
    Decomposition_Algorithm*          aDecompositionAlgorithm )
{
    Tracer tTracer( "XTK", "Integration_Mesh_Generator", "Commit IG Cells To Mesh" ,mXTKModel->mVerboseLevel, 1  );
    
    // iterate through cells that the decomposition constructed
    moris::uint tNumNewCells = aDecompositionAlgorithm->mNumNewCells;

    MORIS_ERROR(
            aDecompositionAlgorithm->mNewCellToVertexConnectivity.size() == aDecompositionAlgorithm->mNewCellChildMeshIndex.size() &&
            aDecompositionAlgorithm->mNewCellChildMeshIndex.size() == aDecompositionAlgorithm->mNewCellCellIndexToReplace.size(),
            "Inconsistent size from decomposition algorithm" );

    // add space to the mesh
    moris::uint tNumStartingCellsControlled = aCutIntegrationMesh->mControlledIgCells.size();
    moris::uint tNumStartingTotalIgCells    = aCutIntegrationMesh->mIntegrationCells.size();

    aCutIntegrationMesh->mControlledIgCells.resize( tNumStartingCellsControlled + tNumNewCells );
    aCutIntegrationMesh->mIntegrationCells.resize( tNumStartingTotalIgCells + tNumNewCells );
    aCutIntegrationMesh->mIntegrationCellToCellGroupIndex.resize( tNumStartingTotalIgCells + tNumNewCells, 0 );
    aCutIntegrationMesh->mIntegrationCellBulkPhase.resize( tNumStartingTotalIgCells + tNumNewCells, MORIS_INDEX_MAX );

    // current index
    moris_index tCellIndex = tNumStartingTotalIgCells;

    // iterate through new and add to the mesh
    for ( moris::uint iCell = 0; iCell < aDecompositionAlgorithm->mNewCellToVertexConnectivity.size(); iCell++ )
    {
        moris_index tCellGroupIndex = aDecompositionAlgorithm->mNewCellChildMeshIndex( iCell );

        std::shared_ptr< IG_Cell_Group > tCellGroup = aCutIntegrationMesh->get_ig_cell_group( tCellGroupIndex );

        // parent cell owner
        moris_index tOwner = aCutIntegrationMesh->get_ig_cell_group_parent_cell( tCellGroupIndex )->get_owner();

        // collect the vertex pointers for the cell
        moris::Cell< moris::mtk::Vertex* > tVertexPointers( aDecompositionAlgorithm->mNewCellToVertexConnectivity( iCell ).size() );

        for ( moris::uint iV = 0; iV < aDecompositionAlgorithm->mNewCellToVertexConnectivity( iCell ).size(); iV++ )
        {
            tVertexPointers( iV ) = aCutIntegrationMesh->get_mtk_vertex_pointer( aDecompositionAlgorithm->mNewCellToVertexConnectivity( iCell )( iV ) );
        }

        bool tReplaceExistingCell = aDecompositionAlgorithm->mNewCellCellIndexToReplace( iCell ) != MORIS_INDEX_MAX;

        // cell index (if I replace one that is the index of this cell)
        moris_index tNewCellIndex = tReplaceExistingCell ? aDecompositionAlgorithm->mNewCellCellIndexToReplace( iCell ) : tCellIndex++;

        std::shared_ptr< xtk::Cell_XTK_No_CM > tNewCell = nullptr;

        // replace the cell, we should only replace cells that are in the same group
        if ( tReplaceExistingCell )
        {
            aCutIntegrationMesh->replace_controlled_ig_cell(
                    tNewCellIndex,
                    aCutIntegrationMesh->get_mtk_cell( tNewCellIndex ).get_id(),
                    aDecompositionAlgorithm->mNewCellCellInfo( iCell ),
                    tVertexPointers );
        }
        // create the new cell no id
        else
        {
            tNewCell = std::make_shared< xtk::Cell_XTK_No_CM >(
                    tNewCellIndex + 1,
                    tNewCellIndex,
                    tOwner,
                    aDecompositionAlgorithm->mNewCellCellInfo( iCell ),
                    tVertexPointers );

            // add the cell to the mesh
            aCutIntegrationMesh->set_integration_cell( tNewCellIndex, tNewCell );

            // add the cell to a child mesh group only if we aren't
            aCutIntegrationMesh->add_cell_to_cell_group( tNewCellIndex, tCellGroupIndex );
        }
    }
}

// ----------------------------------------------------------------------------------

void
Integration_Mesh_Generator::extract_cells_from_cell_groups(
    moris::Cell< std::shared_ptr< IG_Cell_Group > > const& aCellGroups,
    moris::Cell< moris::mtk::Cell* >&                      aCellsInGroups )
{
    Tracer tTracer( "XTK", "Integration_Mesh_Generator", "Extract cells from groups" ,mXTKModel->mVerboseLevel, 1  );
    // count the number of total cells
    moris::uint tNumCells = 0;
    for ( moris::uint i = 0; i < aCellGroups.size(); i++ )
    {
        tNumCells = tNumCells + aCellGroups( i )->mIgCellGroup.size();
    }

    // allocate output
    aCellsInGroups = moris::Cell< moris::mtk::Cell* >();
    aCellsInGroups.reserve( tNumCells );

    for ( moris::uint i = 0; i < aCellGroups.size(); i++ )
    {
        aCellsInGroups.append( aCellGroups( i )->mIgCellGroup );
    }
}
// ----------------------------------------------------------------------------------
void
Integration_Mesh_Generator::compute_ig_cell_bulk_phase(
    Cut_Integration_Mesh* aCutIntegrationMesh )
{
    for ( moris::size_t iCell = 0; iCell < aCutIntegrationMesh->get_num_entities( EntityRank::ELEMENT, 0 ); iCell++ )
    {
        moris_index tBulkPhaseIndex = this->deduce_ig_cell_bulk_phase_index( &aCutIntegrationMesh->get_mtk_cell( iCell ) );

        aCutIntegrationMesh->mIntegrationCellBulkPhase( iCell ) = tBulkPhaseIndex;
    }
}

// ----------------------------------------------------------------------------------

moris_index
Integration_Mesh_Generator::get_max_index( moris::Cell< moris::mtk::Cell* >& aCells )
{
    moris_index tMax = 0;

    for ( moris_index i = 0; i < (moris_index)aCells.size(); i++ )
    {
        if ( aCells( i )->get_index() > tMax )
        {
            tMax = aCells( i )->get_index();
        }
    }

    return tMax;
}

// ----------------------------------------------------------------------------------

moris_index
Integration_Mesh_Generator::deduce_ig_cell_bulk_phase_index( moris::mtk::Cell const* aCell )
{
    // cell vertices
    moris::Cell< moris::mtk::Vertex* > tVertices = aCell->get_vertex_pointers();
    moris::size_t                      tNumGeom  = mGeometryEngine->get_num_geometries();

    // allocate phase on or off value (either 0 or 1)
    Matrix< IndexMat > tPhaseVotes( 1, 2 );
    tPhaseVotes.fill( 0 );

    uint                             tMaxRow = 0;
    uint                             tMaxCol = 0;
    moris::Matrix< moris::IndexMat > tNodalPhaseVals( 1, tNumGeom, MORIS_INDEX_MAX );

    for ( moris::uint i = 0; i < tNumGeom; i++ )
    {
        bool tFoundNonInterfaceNode = false;

        for ( moris::uint iV = 0; iV < tVertices.size(); iV++ )
        {
            if ( !mGeometryEngine->is_interface_vertex( tVertices( iV )->get_index(), i ) )
            {
                moris_index tPhaseIndex = mGeometryEngine->get_node_phase_index_wrt_a_geometry( tVertices( iV )->get_index(), i );
                tFoundNonInterfaceNode  = true;
                tPhaseVotes( tPhaseIndex )++;
            }
        }

        // take the phase with the maximum number of votes
        tPhaseVotes.max( tMaxRow, tMaxCol );
        tNodalPhaseVals( 0, i ) = tMaxCol;
        tPhaseVotes.fill( 0 );
        if ( !tFoundNonInterfaceNode )
        {
            return mGeometryEngine->get_num_phases();
            std::cout << "WARNING Did not find a non-interface node for this element, set to dummy:  " << mGeometryEngine->get_num_phases() << std::endl;
        }
        // MORIS_ERROR( tFoundNonInterfaceNode, "Did not find a non-interface node for this element" );
    }

    return mGeometryEngine->get_elem_phase_index( tNodalPhaseVals );
}

// ----------------------------------------------------------------------------------

void
Integration_Mesh_Generator::deduce_interfaces(
    Cut_Integration_Mesh*                       aCutIntegrationMesh,
    std::shared_ptr< Facet_Based_Connectivity > aFacetConnectivity,
    moris::Cell< moris_index >&                 aInterfaces )
{
    Tracer tTracer( "XTK", "Integration_Mesh_Generator", "Deduce Interface" ,mXTKModel->mVerboseLevel, 1  );

    // all cell groups at this point should be such that an integration cell does not appear twice
    aInterfaces.clear();
    aInterfaces.reserve( aFacetConnectivity->mFacetVertices.size() );


    for ( moris::uint iFacet = 0; iFacet < aFacetConnectivity->mFacetVertices.size(); iFacet++ )
    {
        for ( moris::uint iG = 0; iG < mActiveGeometries.numel(); iG++ )
        {
            moris_index tGeomIndex        = mActiveGeometries( iG );
            bool        tIsInterfaceFacet = true;

            for ( moris::uint iV = 0; iV < aFacetConnectivity->mFacetVertices( iFacet ).size(); iV++ )
            {
                moris::mtk::Vertex* tVertex = aFacetConnectivity->mFacetVertices( iFacet )( iV );
                if ( !mGeometryEngine->is_interface_vertex( tVertex->get_index(), tGeomIndex ) )
                {
                    tIsInterfaceFacet = false;
                    break;
                }
            }

            if ( tIsInterfaceFacet )
            {
                aInterfaces.push_back( iFacet );
            }
        }
    }
}

// ----------------------------------------------------------------------------------

void
Integration_Mesh_Generator::identify_and_construct_subphases(
    Integration_Mesh_Generation_Data*                 aMeshGenerationData,
    Cut_Integration_Mesh*                             aCutIntegrationMesh,
    moris::mtk::Mesh*                                 aBackgroundMesh,
    std::shared_ptr< Cell_Neighborhood_Connectivity > aCutNeighborhood )
{
    Tracer tTracer( "XTK", "Integration_Mesh_Generator", "Identity subphases" ,mXTKModel->mVerboseLevel, 1  );

    // get the number of children meshes
    moris::uint tNumChildMeshes = aCutIntegrationMesh->get_num_ig_cell_groups();

    moris::uint tNumIgCells = aCutIntegrationMesh->get_num_entities( EntityRank::ELEMENT, 0 );

    moris::uint tTotalNumBgCells = aBackgroundMesh->get_num_entities( EntityRank::ELEMENT, 0 );

    // first subphase index (this is incremented by child meshes on call to set_elemental_subphase)
    // proc local subphase index
    moris::moris_index tSubPhaseIndex = tTotalNumBgCells;

    moris_index tMaxSubPhase = 0;

    // subphase groupings
    moris::uint tReserveSize = 2 * tNumChildMeshes + tTotalNumBgCells;
    aCutIntegrationMesh->mSubPhaseCellGroups.reserve( tReserveSize );
    aCutIntegrationMesh->mSubPhaseCellGroups.resize( tTotalNumBgCells, nullptr );
    aCutIntegrationMesh->mSubPhaseBulkPhase.reserve( tReserveSize );
    aCutIntegrationMesh->mSubPhaseBulkPhase.resize( tTotalNumBgCells );
    aCutIntegrationMesh->mSubPhaseParentCell.reserve( tReserveSize );
    aCutIntegrationMesh->mSubPhaseParentCell.resize( tTotalNumBgCells );
    aCutIntegrationMesh->mSubPhaseIds.reserve( tReserveSize );
    aCutIntegrationMesh->mSubPhaseIds.resize( tTotalNumBgCells );
    aCutIntegrationMesh->mParentCellToSubphase.resize( tTotalNumBgCells );
    aCutIntegrationMesh->mParentCellHasChildren.resize( tTotalNumBgCells );

    aCutIntegrationMesh->mIntegrationCellToSubphaseIndex = moris::Cell< moris_index >( tNumIgCells, MORIS_INDEX_MAX );

    moris::Cell< moris_index > tSubphaseIndices;
    tSubphaseIndices.reserve( 10 );

    // iterate over children meshes and perform local flood-fill
    for ( auto& iCell : aMeshGenerationData->mAllIntersectedBgCellInds )
    {
        std::shared_ptr< IG_Cell_Group > tIgCellGroup = aCutIntegrationMesh->get_ig_cell_group( iCell );

        moris::mtk::Cell* tParentCell = aCutIntegrationMesh->get_ig_cell_group_parent_cell( iCell );

        MORIS_ASSERT( tParentCell->get_index() == iCell, "Index mismatch parent index should allign wiht parent cell index" );

        // flood fill this group using the bulk phases
        moris::Matrix< moris::IndexMat > tLocalFloodFill = this->flood_fill_ig_cell_group( aCutIntegrationMesh, aCutNeighborhood, tIgCellGroup, tMaxSubPhase );

        // put first subphase in spot of parent cell
        aCutIntegrationMesh->mSubPhaseCellGroups( tParentCell->get_index() ) = std::make_shared< IG_Cell_Group >();
        aCutIntegrationMesh->mSubPhaseBulkPhase( tParentCell->get_index() )  = MORIS_INDEX_MAX;
        aCutIntegrationMesh->mSubPhaseParentCell( tParentCell->get_index() ) = tParentCell;
        aCutIntegrationMesh->mSubPhaseIds( tParentCell->get_index() )        = tParentCell->get_id();
        aCutIntegrationMesh->mParentCellToSubphase( tParentCell->get_index() ).push_back( tParentCell->get_index() );
        aCutIntegrationMesh->mParentCellHasChildren( tParentCell->get_index() ) = ( moris_index ) true;
        tSubphaseIndices.push_back( tParentCell->get_index() );

        // allocate the other subphases
        for ( moris::uint iSP = 1; iSP < (uint)tMaxSubPhase + 1; iSP++ )
        {
            aCutIntegrationMesh->mSubPhaseCellGroups.push_back( std::make_shared< IG_Cell_Group >() );
            aCutIntegrationMesh->mSubPhaseBulkPhase.push_back( MORIS_INDEX_MAX );
            aCutIntegrationMesh->mSubPhaseParentCell.push_back( tParentCell );
            aCutIntegrationMesh->mSubPhaseIds.push_back( MORIS_INDEX_MAX );
            aCutIntegrationMesh->mParentCellToSubphase( tParentCell->get_index() ).push_back( tSubPhaseIndex );
            tSubphaseIndices.push_back( tSubPhaseIndex++ );
        }

        // iterate through the local flood fill results
        for ( moris::uint iFF = 0; iFF < tLocalFloodFill.numel(); iFF++ )
        {
            moris_index tFFVal   = tLocalFloodFill( iFF );
            moris_index tSPIndex = tSubphaseIndices( tFFVal );

            // put the cell into the subphase group
            aCutIntegrationMesh->mSubPhaseCellGroups( tSPIndex )->mIgCellGroup.push_back( tIgCellGroup->mIgCellGroup( iFF ) );

            aCutIntegrationMesh->mIntegrationCellToSubphaseIndex( tIgCellGroup->mIgCellGroup( iFF )->get_index() ) = tSPIndex;

            // add bulk phase if I haven't yet
            if ( aCutIntegrationMesh->mSubPhaseBulkPhase( tSPIndex ) == MORIS_INDEX_MAX )
            {
                aCutIntegrationMesh->mSubPhaseBulkPhase( tSPIndex ) = aCutIntegrationMesh->get_cell_bulk_phase( tIgCellGroup->mIgCellGroup( iFF )->get_index() );
            }
        }

        aCutIntegrationMesh->set_child_mesh_subphase( iCell, tSubphaseIndices );

        tSubphaseIndices.clear();
    }

    // iterate over background cells and make the subphase group contain only them
    for ( moris::size_t i = 0; i < aBackgroundMesh->get_num_elems(); i++ )
    {
        if ( aCutIntegrationMesh->mSubPhaseCellGroups( i ) == nullptr )
        {
            moris::mtk::Cell* tCell                                          = &aCutIntegrationMesh->get_mtk_cell( (moris_index)i );
            aCutIntegrationMesh->mSubPhaseCellGroups( i )                    = std::make_shared< IG_Cell_Group >( 1 );
            aCutIntegrationMesh->mSubPhaseCellGroups( i )->mIgCellGroup( 0 ) = tCell;

            aCutIntegrationMesh->mSubPhaseBulkPhase( i )  = aCutIntegrationMesh->get_cell_bulk_phase( i );
            aCutIntegrationMesh->mSubPhaseParentCell( i ) = tCell;
            aCutIntegrationMesh->mSubPhaseIds( i )        = tCell->get_id();
            aCutIntegrationMesh->mParentCellToSubphase( i ).push_back( i );
            aCutIntegrationMesh->mIntegrationCellToSubphaseIndex( i ) = i;
            aCutIntegrationMesh->mParentCellHasChildren( i )          = ( moris_index ) false;
        }
    }


    aCutIntegrationMesh->mSubPhaseCellGroups.shrink_to_fit();
    aCutIntegrationMesh->mSubPhaseBulkPhase.shrink_to_fit();
    aCutIntegrationMesh->mSubPhaseParentCell.shrink_to_fit();
    aCutIntegrationMesh->mSubPhaseIds.shrink_to_fit();

    aCutIntegrationMesh->mOwnedSubphaseGroupsInds.reserve( aCutIntegrationMesh->mSubPhaseCellGroups.size() );
    aCutIntegrationMesh->mNotOwnedSubphaseGroupsInds.reserve( aCutIntegrationMesh->mSubPhaseCellGroups.size() );

    moris_index tParRank = moris::par_rank();
    // sort into ownership groups
    for ( moris::uint i = 0; i < aCutIntegrationMesh->mSubPhaseCellGroups.size(); i++ )
    {
        aCutIntegrationMesh->mSubPhaseParentCell( i )->get_owner() == tParRank ? aCutIntegrationMesh->mOwnedSubphaseGroupsInds.push_back( (moris_index)i ) : aCutIntegrationMesh->mNotOwnedSubphaseGroupsInds.push_back( (moris_index)i );
    }

    aCutIntegrationMesh->mOwnedSubphaseGroupsInds.shrink_to_fit();
    aCutIntegrationMesh->mNotOwnedSubphaseGroupsInds.shrink_to_fit();

    this->assign_subphase_glob_ids( aCutIntegrationMesh, aBackgroundMesh );
}

// ----------------------------------------------------------------------------------

void
Integration_Mesh_Generator::construct_subphase_neighborhood(
    Cut_Integration_Mesh*                                                aCutIntegrationMesh,
    moris::mtk::Mesh*                                                    aBackgroundMesh,
    std::shared_ptr< Facet_Based_Connectivity >                          aFacetConnectivity,
    moris::Cell< std::shared_ptr< moris::Cell< moris::moris_index > > >* aBgFacetToChildFacet,
    std::shared_ptr< Subphase_Neighborhood_Connectivity >                aSubphaseNeighborhood )
{
    Tracer tTracer( "XTK", "Integration_Mesh_Generator", "Subphase Neighborhood", mXTKModel->mVerboseLevel, 1 );
    aSubphaseNeighborhood->mSubphaseToSubPhase.resize( aCutIntegrationMesh->get_num_subphases() );
    aSubphaseNeighborhood->mSubphaseToSubPhaseMySideOrds.resize( aCutIntegrationMesh->get_num_subphases() );
    aSubphaseNeighborhood->mSubphaseToSubPhaseNeighborSideOrds.resize( aCutIntegrationMesh->get_num_subphases() );
    aSubphaseNeighborhood->mTransitionNeighborCellLocation.resize( aCutIntegrationMesh->get_num_subphases() );

    // allocate the pointers in subphase neighborhood
    for ( moris::uint i = 0; i < aCutIntegrationMesh->get_num_subphases(); i++ )
    {
        aSubphaseNeighborhood->mSubphaseToSubPhase( i )                 = std::make_shared< moris::Cell< moris_index > >();
        aSubphaseNeighborhood->mSubphaseToSubPhaseMySideOrds( i )       = std::make_shared< moris::Cell< moris_index > >();
        aSubphaseNeighborhood->mSubphaseToSubPhaseNeighborSideOrds( i ) = std::make_shared< moris::Cell< moris_index > >();
        aSubphaseNeighborhood->mTransitionNeighborCellLocation( i )     = std::make_shared< moris::Cell< moris_index > >();

        aSubphaseNeighborhood->mSubphaseToSubPhase( i )->reserve( 4 );
        aSubphaseNeighborhood->mSubphaseToSubPhaseMySideOrds( i )->reserve( 4 );
        aSubphaseNeighborhood->mSubphaseToSubPhaseNeighborSideOrds( i )->reserve( 4 );
        aSubphaseNeighborhood->mTransitionNeighborCellLocation( i )->reserve( 4 );
    }

    // I loop throuhg cells then side ordinals to capture the transition between adaptively refined meshes.
    // alternatively you could try looping through aBgFacetToChildFacet and then handle the subphase neigborhood for the transition locations seperately
    //iterate through background cells
    moris::uint tBGCells = aBackgroundMesh->get_num_entities( EntityRank::ELEMENT );
    for ( moris::moris_index iC = 0; iC < (moris::moris_index)tBGCells; iC++ )
    {
        mtk::Cell const* tCurrentCell = &aBackgroundMesh->get_mtk_cell( iC );

        // get the cells attached to the facet
        Matrix< IndexMat > tCellToCellSideIndex = aBackgroundMesh->get_elements_connected_to_element_and_face_ind_loc_inds( iC );
        Matrix< IndexMat > tCellToCellSideOrd   = aBackgroundMesh->get_elements_connected_to_element_and_face_ord_loc_inds( iC );

        Cell< mtk::Cell const* > tCells( tCellToCellSideOrd.numel() );
        aBackgroundMesh->get_mtk_cells( tCellToCellSideOrd.get_row( 0 ), tCells );

        // iterate through neighbor
        for ( moris::uint iN = 0; iN < tCellToCellSideOrd.n_cols(); iN++ )
        {
            // neighbor cell
            mtk::Cell const* tOtherCell = &aBackgroundMesh->get_mtk_cell( tCellToCellSideIndex( 0, iN ) );

            // facet ordinal shared for current neighbors
            moris_index tFacetIndex             = tCellToCellSideIndex( 1, iN );
            moris_index tMyOrdinal              = tCellToCellSideOrd( 1, iN );
            moris_index tNeighborOrdinal        = tCellToCellSideOrd( 2, iN );
            moris_index tTransitionCellLocation = tCellToCellSideOrd( 3, iN );

            Cell< moris::moris_index > tMyCellSubphaseIndices( 0 );
            Cell< moris::moris_index > tRepresentativeIgCells( 0 );
            Cell< moris::moris_index > tRepresentativeIgCellsOrdinal( 0 );
            this->collect_subphases_attached_to_facet_on_cell( aCutIntegrationMesh, tCurrentCell, tMyOrdinal, tFacetIndex, aFacetConnectivity, ( *aBgFacetToChildFacet )( tFacetIndex ), tMyCellSubphaseIndices, tRepresentativeIgCells, tRepresentativeIgCellsOrdinal );
            
            // transitioning between mesh levels
            if ( tTransitionCellLocation != MORIS_INDEX_MAX )
            {
                Matrix<IndexMat> tNeighborCellToFacetIndex = aBackgroundMesh->get_entity_connected_to_entity_loc_inds(tOtherCell->get_index(),EntityRank::ELEMENT, aBackgroundMesh->get_facet_rank());
                Cell< moris::moris_index > tNeighborSubphaseIndices( 0 );
                Cell< moris::moris_index > tNeighborRepresentativeIgCells( 0 );
                Cell< moris::moris_index > tNeighborRepresentativeIgCellsOrdinal( 0 );
                this->collect_subphases_attached_to_facet_on_cell( aCutIntegrationMesh, tOtherCell, tNeighborOrdinal, tNeighborCellToFacetIndex(tNeighborOrdinal), aFacetConnectivity, ( *aBgFacetToChildFacet )( tNeighborCellToFacetIndex(tNeighborOrdinal) ), tNeighborSubphaseIndices, tNeighborRepresentativeIgCells, tNeighborRepresentativeIgCellsOrdinal );

                for ( moris::uint i = 0; i < tMyCellSubphaseIndices.size(); i++ )
                {
                    moris_index tMySubphaseIndex  = tMyCellSubphaseIndices( i );
                    moris_index tMyBulkIndex = aCutIntegrationMesh->get_subphase_bulk_phase( tMySubphaseIndex );
                    for ( moris::uint j = 0; j < tNeighborSubphaseIndices.size(); j++ )
                    {
                        moris_index tNeighborBulkIndex     = aCutIntegrationMesh->get_subphase_bulk_phase( tNeighborSubphaseIndices(j) );
                        moris_index tNeighborSubphaseIndex = tNeighborSubphaseIndices( j );

                        if ( tMyBulkIndex == tNeighborBulkIndex )
                        {
                            aSubphaseNeighborhood->mSubphaseToSubPhase( tMySubphaseIndex )->push_back( tNeighborSubphaseIndex );
                            aSubphaseNeighborhood->mSubphaseToSubPhaseMySideOrds( tMySubphaseIndex )->push_back( tMyOrdinal );
                            aSubphaseNeighborhood->mSubphaseToSubPhaseNeighborSideOrds( tMySubphaseIndex )->push_back( tNeighborOrdinal );
                            aSubphaseNeighborhood->mTransitionNeighborCellLocation( tMySubphaseIndex )->push_back( tTransitionCellLocation );
                        }
                    }
                }
            }

            else
            {
                // iterate over subphases and add to neighborhood
                for ( moris::uint i = 0; i < tMyCellSubphaseIndices.size(); i++ )
                {
                    moris_index tMySubphaseIndex = tMyCellSubphaseIndices( i );
                    moris_index tMyIgCellIndex   = tRepresentativeIgCells( i );
                    moris_index tMyIgCellSideOrd = tRepresentativeIgCellsOrdinal( i );


                    // handle the case where we transition between background cell and triangulated cells
                    if ( !aCutIntegrationMesh->parent_cell_has_children( tCurrentCell->get_index() ) || !aCutIntegrationMesh->parent_cell_has_children( tOtherCell->get_index() ) )
                    {
                        Cell< moris::moris_index > tNeighborSubphaseIndices( 0 );
                        Cell< moris::moris_index > tNeighborRepresentativeIgCells( 0 );
                        Cell< moris::moris_index > tNeighborRepresentativeIgCellsOrdinal( 0 );
                        this->collect_subphases_attached_to_facet_on_cell( aCutIntegrationMesh, tOtherCell, tNeighborOrdinal, tFacetIndex, aFacetConnectivity, ( *aBgFacetToChildFacet )( tFacetIndex ), tNeighborSubphaseIndices, tNeighborRepresentativeIgCells, tNeighborRepresentativeIgCellsOrdinal );

                        // iterate through neighbors
                        for ( const auto& iNeighSp : tNeighborSubphaseIndices )
                        {
                            aSubphaseNeighborhood->mSubphaseToSubPhase( tMySubphaseIndex )->push_back( iNeighSp );
                            aSubphaseNeighborhood->mSubphaseToSubPhaseMySideOrds( tMySubphaseIndex )->push_back( tMyOrdinal );
                            aSubphaseNeighborhood->mSubphaseToSubPhaseNeighborSideOrds( tMySubphaseIndex )->push_back( tNeighborOrdinal );
                            aSubphaseNeighborhood->mTransitionNeighborCellLocation( tMySubphaseIndex )->push_back( tTransitionCellLocation );
                        }
                    }

                    else
                    {
                        // figure out the neighbor subphase index
                        const moris_index& tMyIgCellOrdInFacetConn = aFacetConnectivity->get_cell_ordinal( tMyIgCellIndex );
                        moris_index        tIgFacetIndex           = aFacetConnectivity->mCellToFacet( tMyIgCellOrdInFacetConn )( tMyIgCellSideOrd );
                        moris_index        tNeighborSubphaseIndex  = MORIS_INDEX_MAX;
                        // iterate through cells on facet and get the one that is not my cell
                        for ( const auto& iCell : aFacetConnectivity->mFacetToCell( tIgFacetIndex ) )
                        {
                            if ( iCell->get_index() != tMyIgCellIndex )
                            {
                                tNeighborSubphaseIndex = aCutIntegrationMesh->get_ig_cell_subphase_index( iCell->get_index() );
                            }
                        }
                        // MORIS_ASSERT( aCutIntegrationMesh->get_subphase_bulk_phase( tNeighborSubphaseIndex ) == aCutIntegrationMesh->get_subphase_bulk_phase( tMySubphaseIndex ), "Subphase bulk phase mismatch" );

                        if ( aCutIntegrationMesh->get_subphase_bulk_phase( tNeighborSubphaseIndex ) == aCutIntegrationMesh->get_subphase_bulk_phase( tMySubphaseIndex ) )
                        {
                            aSubphaseNeighborhood->mSubphaseToSubPhase( tMySubphaseIndex )->push_back( tNeighborSubphaseIndex );
                            aSubphaseNeighborhood->mSubphaseToSubPhaseMySideOrds( tMySubphaseIndex )->push_back( tMyOrdinal );
                            aSubphaseNeighborhood->mSubphaseToSubPhaseNeighborSideOrds( tMySubphaseIndex )->push_back( tNeighborOrdinal );
                            aSubphaseNeighborhood->mTransitionNeighborCellLocation( tMySubphaseIndex )->push_back( tTransitionCellLocation );
                        }
                    }
                }
            }
        }
    }
}

// ----------------------------------------------------------------------------------

void
Integration_Mesh_Generator::collect_subphases_attached_to_facet_on_cell(
    Cut_Integration_Mesh*                                aCutIntegrationMesh,
    moris::mtk::Cell const*                              aBGCell,
    moris::moris_index                                   aFacetOrdinal,
    moris::moris_index                                   aSharedFacetIndex,
    std::shared_ptr< Facet_Based_Connectivity >          aFacetConnectivity,
    std::shared_ptr< moris::Cell< moris::moris_index > > aBgFacetToChildrenFacets,
    Cell< moris::moris_index >&                          aSubphaseIndices,
    Cell< moris::moris_index >&                          aRepresentativeIgCells,
    Cell< moris::moris_index >&                          aRepresentativeIgCellsOrdinal )
{
    aSubphaseIndices.clear();

    bool tParentHasChildren = aCutIntegrationMesh->parent_cell_has_children( aBGCell->get_index() );

    if ( tParentHasChildren )
    {
        std::unordered_map< moris_index, moris_index > tSubphaseMap;
        MORIS_ASSERT( aBgFacetToChildrenFacets != nullptr, "Null ptr on facet that should have children facets" );

        // iterate through child facets attached to bg facet
        for ( moris::uint iChildFacet = 0; iChildFacet < aBgFacetToChildrenFacets->size(); iChildFacet++ )
        {
            moris_index tChildCellFacetIndex = ( *aBgFacetToChildrenFacets )( iChildFacet );

            for ( moris::uint iCell = 0; iCell < aFacetConnectivity->mFacetToCell( tChildCellFacetIndex ).size(); iCell++ )
            {
                moris_index tCellIndex       = aFacetConnectivity->mFacetToCell( tChildCellFacetIndex )( iCell )->get_index();
                moris_index tSubphaseIndex   = aCutIntegrationMesh->get_ig_cell_subphase_index( tCellIndex );
                moris_index tParentCellIndex = aCutIntegrationMesh->get_subphase_parent_cell( tSubphaseIndex )->get_index();
                if ( tParentCellIndex == aBGCell->get_index() )
                {
                    if ( tSubphaseMap.find( tSubphaseIndex ) == tSubphaseMap.end() )
                    {
                        aSubphaseIndices.push_back( tSubphaseIndex );
                        aRepresentativeIgCells.push_back( tCellIndex );
                        aRepresentativeIgCellsOrdinal.push_back( aFacetConnectivity->mFacetToCellEdgeOrdinal( tChildCellFacetIndex )( iCell ) );
                        tSubphaseMap[tSubphaseIndex] = 1;
                    }
                }
            }
        }
    }
    else
    {
        MORIS_ASSERT( aCutIntegrationMesh->get_parent_cell_subphases( aBGCell->get_index() ).size() == 1, "one subphase needs to be present in this case" );
        aSubphaseIndices = { { aBGCell->get_index() } };
        aRepresentativeIgCells.push_back( aBGCell->get_index() );
        aRepresentativeIgCellsOrdinal.push_back( aFacetOrdinal );
    }
}

// ----------------------------------------------------------------------------------

void
Integration_Mesh_Generator::collect_ig_cells_and_side_ords_on_bg_facet(
    Cut_Integration_Mesh*             aCutIntegrationMesh,
    moris::moris_index                aBackgroundFacetIndex,
    moris::Cell< moris::mtk::Cell* >& aIgCell,
    moris::Cell< moris_index >&       aIgCellSideOrds )
{
    aIgCell.clear();
    aIgCellSideOrds.clear();

    moris::Cell< std::shared_ptr< moris::Cell< moris::moris_index > > > const& tBGFacetToChildFacet = aCutIntegrationMesh->get_background_facet_to_child_facet_connectivity();


    if ( tBGFacetToChildFacet( aBackgroundFacetIndex ) == nullptr )
    {
        return;
    }

    else
    {
        // facet connectivity
        std::shared_ptr< Facet_Based_Connectivity > tFacetConn = aCutIntegrationMesh->get_face_connectivity();

        moris::Cell< moris::moris_index > const& tBGFacetToIgFacet = *tBGFacetToChildFacet( aBackgroundFacetIndex );

        // iterate through child facets attached to bg facet
        for ( moris::uint iChildFacet = 0; iChildFacet < tBGFacetToIgFacet.size(); iChildFacet++ )
        {
            moris_index tChildCellFacetIndex = tBGFacetToIgFacet( iChildFacet );

            for ( moris::uint iCell = 0; iCell < tFacetConn->mFacetToCell( tChildCellFacetIndex ).size(); iCell++ )
            {
                aIgCell.push_back( tFacetConn->mFacetToCell( tChildCellFacetIndex )( iCell ) );
                aIgCellSideOrds.push_back( tFacetConn->mFacetToCellEdgeOrdinal( tChildCellFacetIndex )( iCell ) );
            }
        }
    }
}

// ----------------------------------------------------------------------------------

void
Integration_Mesh_Generator::assign_subphase_glob_ids(
    Cut_Integration_Mesh* aCutIntegrationMesh,
    moris::mtk::Mesh*     aBackgroundMesh )
{
    // Get the number of subphases
    moris_id tNumSubphases = (moris_id)aCutIntegrationMesh->get_num_subphases();

    // Allocate global element ids starting at the maximum id in the background mesh (these need to be give to the children meshes)
    moris::moris_id tSubphaseIdOffset = aCutIntegrationMesh->allocate_subphase_ids( tNumSubphases );

    // // set subphase ids in the children meshes which I own
    moris::Cell< moris_index >& tOwnedSubphases = aCutIntegrationMesh->get_owned_subphase_indices();
    for ( moris::size_t i = 0; i < tOwnedSubphases.size(); i++ )
    {
        moris_index tSubphaseIndex = tOwnedSubphases( i );
        // moris::mtk::Cell* tParentCell = aCutIntegrationMesh->get_subphase_parent_cell(tSubphaseIndex);

        // moris_id tCellId = tParentCell->get_id();

        if ( aCutIntegrationMesh->mSubPhaseIds( tSubphaseIndex ) == MORIS_INDEX_MAX )
        {
            aCutIntegrationMesh->mSubPhaseIds( tSubphaseIndex ) = tSubphaseIdOffset++;
        }
    }

    // prepare outward requests
    Cell< Cell< moris_id > >                 tNotOwnedSubphasesToProcs;
    Cell< moris::Matrix< IdMat > >           tParentCellIds;
    Cell< moris::Matrix< IdMat > >           tChildCellIds;
    Cell< moris::Matrix< IdMat > >           tNumChildCellsInSubphase;
    Cell< uint >                             tProcRanks;
    std::unordered_map< moris_id, moris_id > tProcRankToDataIndex;
    this->prepare_subphase_identifier_requests(
        aCutIntegrationMesh,
        tNotOwnedSubphasesToProcs,
        tParentCellIds,
        tChildCellIds,
        tNumChildCellsInSubphase,
        tProcRanks,
        tProcRankToDataIndex );

    // send requests
    moris::uint tMPITag = 221;
    mXTKModel->send_outward_requests( tMPITag, tProcRanks, tParentCellIds );
    mXTKModel->send_outward_requests( tMPITag + 1, tProcRanks, tChildCellIds );
    mXTKModel->send_outward_requests( tMPITag + 2, tProcRanks, tNumChildCellsInSubphase );

    barrier();

    // receive requests
    moris::Cell< Matrix< IndexMat > > tReceivedParentCellIds;
    moris::Cell< Matrix< IndexMat > > tFirstChildCellIds;
    moris::Cell< Matrix< IndexMat > > tReceivedNumChildCellsInSubphase;
    moris::Cell< uint >               tProcsReceivedFrom1;
    moris::Cell< uint >               tProcsReceivedFrom2;
    mXTKModel->inward_receive_requests( tMPITag, 1, tReceivedParentCellIds, tProcsReceivedFrom1 );
    mXTKModel->inward_receive_requests( tMPITag + 1, 1, tFirstChildCellIds, tProcsReceivedFrom2 );
    mXTKModel->inward_receive_requests( tMPITag + 2, 1, tReceivedNumChildCellsInSubphase, tProcsReceivedFrom2 );
    MORIS_ASSERT( tProcsReceivedFrom1.size() == tProcsReceivedFrom2.size(), "Size mismatch between procs received from child cell ids and number of child cells" );

    // prepare answers
    Cell< Matrix< IndexMat > > tSubphaseIds;
    this->prepare_subphase_id_answers( aCutIntegrationMesh, aBackgroundMesh, tReceivedParentCellIds, tFirstChildCellIds, tReceivedNumChildCellsInSubphase, tSubphaseIds );

    // return information
    mXTKModel->return_request_answers( tMPITag + 2, tSubphaseIds, tProcsReceivedFrom1 );

    barrier();

    // receive the answers
    Cell< Matrix< IndexMat > > tReceivedSubphaseIds;
    mXTKModel->inward_receive_request_answers( tMPITag + 2, 1, tProcRanks, tReceivedSubphaseIds );

    // add child cell ids to not owned child meshes
    this->handle_received_subphase_id_request_answers( aCutIntegrationMesh, tNotOwnedSubphasesToProcs, tReceivedSubphaseIds );

    aCutIntegrationMesh->setup_glob_to_loc_subphase_map();

    barrier();
}

// ----------------------------------------------------------------------------------

void
Integration_Mesh_Generator::prepare_subphase_identifier_requests(
    Cut_Integration_Mesh*                     aCutIntegrationMesh,
    Cell< Cell< moris_id > >&                 aNotOwnedSubphasesToProcs,
    Cell< moris::Matrix< IdMat > >&           aParentCellIds,
    Cell< moris::Matrix< IdMat > >&           aChildCellIds,
    Cell< moris::Matrix< IdMat > >&           aNumChildCellsInSubphase,
    Cell< uint >&                             aProcRanks,
    std::unordered_map< moris_id, moris_id >& aProcRankToDataIndex )
{
    // access the communication table
    Matrix< IdMat > tCommTable = aCutIntegrationMesh->get_communication_table();

    // resize proc ranks and setup map to comm table
    aProcRanks.resize( tCommTable.numel() );
    for ( moris::uint i = 0; i < tCommTable.numel(); i++ )
    {
        aProcRankToDataIndex[tCommTable( i )] = i;
        aProcRanks( i )                       = ( tCommTable( i ) );
        aNotOwnedSubphasesToProcs.push_back( Cell< moris_id >( 0 ) );
    }

    // ask owning processor about child element ids
    moris::Cell< moris_index >& tNotOwnedSubphases = aCutIntegrationMesh->get_not_owned_subphase_indices();

    // not owned subphases
    for ( moris::size_t i = 0; i < tNotOwnedSubphases.size(); i++ )
    {
        moris_index tOwnerProc = aCutIntegrationMesh->get_subphase_parent_cell( tNotOwnedSubphases( i ) )->get_owner();

        moris_index tProcDataIndex = aProcRankToDataIndex[tOwnerProc];
        aNotOwnedSubphasesToProcs( tProcDataIndex ).push_back( tNotOwnedSubphases( i ) );
    }

    aParentCellIds.resize( aNotOwnedSubphasesToProcs.size() );
    aChildCellIds.resize( aNotOwnedSubphasesToProcs.size() );
    aNumChildCellsInSubphase.resize( aNotOwnedSubphasesToProcs.size() );

    // iterate through procs and child meshes shared with that processor
    for ( moris::size_t i = 0; i < aNotOwnedSubphasesToProcs.size(); i++ )
    {
        // number of child meshes shared with this processor
        moris::uint tNumSubphases = aNotOwnedSubphasesToProcs( i ).size();

        // allocate matrix
        aParentCellIds( i ).resize( 1, tNumSubphases );
        aChildCellIds( i ).resize( 1, tNumSubphases );
        aNumChildCellsInSubphase( i ).resize( 1, tNumSubphases );

        if ( tNumSubphases == 0 )
        {
            aParentCellIds( i ).resize( 1, 1 );
            aChildCellIds( i ).resize( 1, 1 );
            aNumChildCellsInSubphase( i ).resize( 1, 1 );
            aParentCellIds( i )( 0 )           = MORIS_INDEX_MAX;
            aChildCellIds( i )( 0 )            = MORIS_INDEX_MAX;
            aNumChildCellsInSubphase( i )( 0 ) = MORIS_INDEX_MAX;
        }

        for ( moris::uint j = 0; j < tNumSubphases; j++ )
        {
            moris_index tSubphaseIndex = aNotOwnedSubphasesToProcs( i )( j );

            // get the subphase cell group
            std::shared_ptr< IG_Cell_Group > tSubphase = aCutIntegrationMesh->get_subphase_ig_cells( tSubphaseIndex );

            moris::mtk::Cell* tParentCell = aCutIntegrationMesh->get_subphase_parent_cell( tSubphaseIndex );

            aParentCellIds( i )( j )           = tParentCell->get_id();
            aChildCellIds( i )( j )            = tSubphase->mIgCellGroup( 0 )->get_id();
            aNumChildCellsInSubphase( i )( j ) = tSubphase->mIgCellGroup.size();
        }
    }
}

// ----------------------------------------------------------------------------------

void
Integration_Mesh_Generator::prepare_subphase_id_answers(
    Cut_Integration_Mesh*       aCutIntegrationMesh,
    moris::mtk::Mesh*           aBackgroundMesh,
    Cell< Matrix< IndexMat > >& aReceivedParentCellIds,
    Cell< Matrix< IndexMat > >& aFirstChildCellIds,
    Cell< Matrix< IndexMat > >& aReceivedNumChildCellsInSubphase,
    Cell< Matrix< IndexMat > >& aSubphaseIds )
{
    MORIS_ASSERT( aReceivedParentCellIds.size() == aFirstChildCellIds.size(),
        "Mismatch in received parent cell ids and received parent cell number of children" );

    // allocate answer size
    aSubphaseIds.resize( aReceivedParentCellIds.size() );

    // iterate through received data
    for ( moris::uint i = 0; i < aReceivedParentCellIds.size(); i++ )
    {
        uint tNumReceivedReqs = aReceivedParentCellIds( i ).n_cols();

        aSubphaseIds( i ).resize( 1, tNumReceivedReqs );

        if ( aReceivedParentCellIds( i )( 0, 0 ) != MORIS_INDEX_MAX )
        {
            // iterate through received requests
            for ( moris::uint j = 0; j < tNumReceivedReqs; j++ )
            {
                // Child cell information
                moris_id tChildCellId = aFirstChildCellIds( i )( j );

                // subphase index on this processor
                moris_index tChildCellIndex = aCutIntegrationMesh->get_loc_entity_ind_from_entity_glb_id( tChildCellId, EntityRank::ELEMENT );

                // subphase index
                moris_index tSubphaseIndex = aCutIntegrationMesh->get_ig_cell_subphase_index( tChildCellIndex );

                moris_index tSubphaseId = aCutIntegrationMesh->get_subphase_id( tSubphaseIndex );

                MORIS_ERROR( aCutIntegrationMesh->get_subphase_parent_cell( tSubphaseIndex )->get_id() == aReceivedParentCellIds( i )( 0, j ), "Subphase parent cell id discrepency" );
                MORIS_ERROR( (moris_index)aCutIntegrationMesh->get_subphase_ig_cells( tSubphaseIndex )->mIgCellGroup.size() == aReceivedNumChildCellsInSubphase( i )( j ), "Number of cells in subphase discrepency" );

                MORIS_ERROR( tSubphaseId != MORIS_ID_MAX, "Child cell id not found in child mesh" );

                // place in return data
                aSubphaseIds( i )( j ) = tSubphaseId;
            }
        }
    }
}

// ----------------------------------------------------------------------------------

void
Integration_Mesh_Generator::handle_received_subphase_id_request_answers(
    Cut_Integration_Mesh*              aCutIntegrationMesh,
    Cell< Cell< moris_index > > const& aSubphaseIndices,
    Cell< Matrix< IndexMat > > const&  aReceivedSubphaseIds )
{
    // iterate through received data
    for ( moris::uint i = 0; i < aSubphaseIndices.size(); i++ )
    {
        uint tNumReceivedReqs = aSubphaseIndices( i ).size();

        // iterate through received requests
        for ( moris::uint j = 0; j < tNumReceivedReqs; j++ )
        {
            moris_index tSubphaseIndex                          = aSubphaseIndices( i )( j );
            aCutIntegrationMesh->mSubPhaseIds( tSubphaseIndex ) = aReceivedSubphaseIds( i )( j );
        }
    }
}

// ----------------------------------------------------------------------------------

moris::Matrix< moris::IndexMat >
Integration_Mesh_Generator::flood_fill_ig_cell_group(
    Cut_Integration_Mesh*                             aCutIntegrationMesh,
    std::shared_ptr< Cell_Neighborhood_Connectivity > aCutNeighborhood,
    std::shared_ptr< IG_Cell_Group >                  aIgCellGroup,
    moris_index&                                      aMaxValueAssigned )
{

    // Active phase index
    moris::size_t tPhaseIndex = 0;

    // Number of elements in the flood fill
    moris::size_t tNumElements = aIgCellGroup->mIgCellGroup.size();

    // Number of Elements with Set Phases (This allows for early termination of code if every element has been set)
    moris::size_t tNumPhasesSet = 0;

    // Current Element Index
    moris::moris_index tElementIndex = 0;

    // Neighbor Element Phase
    moris::size_t tNeighborPhase = 0;

    // Neighbor Element local Index in the list aElementsToInclude
    moris::moris_index tNeighborIndex = 0;
    moris::moris_index tNeighborOrd   = 0;

    // Current Subphase
    moris::size_t tCurrentSubphase = 0;

    // Track which elements have their phase set
    moris::Matrix< moris::DDBMat > tPhaseSet( 1, tNumElements, 0 );

    // Initialize element sub-phases
    moris::Matrix< moris::IndexMat > tElementSubphase( 1, tNumElements, MORIS_INDEX_MAX );

    // Initialize Active Front
    moris::size_t                    tActiveFrontCount   = 0;
    moris::size_t                    tActiveFrontElement = 0;
    moris::Matrix< moris::IndexMat > tActiveFront( 1, 10 * tNumElements + 1, 0 );

    // Map between the active element indexes provided and their corresponding iE (Only needed if all elements are not included)
    // key   - Element Index
    // value - flood fill local index
    std::unordered_map< moris::moris_index, moris::moris_index > tElementToLocalIndex;

    for ( moris::size_t iE = 0; iE < tNumElements; iE++ )
    {
        tElementToLocalIndex[aIgCellGroup->mIgCellGroup( iE )->get_index()] = iE;
    }

    // Loop over all elements
    for ( moris::size_t iE = 0; iE < tNumElements; iE++ )
    {
        tElementIndex = aIgCellGroup->mIgCellGroup( iE )->get_index();

        // If this element phase has not been set
        if ( !tPhaseSet( 0, iE ) )
        {
            // Phase Index of the element
            tPhaseIndex = aCutIntegrationMesh->get_cell_bulk_phase( tElementIndex );

            // Set the elements subphase value
            tElementSubphase( 0, iE ) = tCurrentSubphase;

            // Mark this element as set
            tPhaseSet( 0, iE ) = 1;

            // Update active front
            // iterate through the cell index neighbors
            for ( moris::size_t iN = 0; iN < aCutNeighborhood->mNeighborCells( tElementIndex )->size(); iN++ )
            {
                tNeighborIndex = (*aCutNeighborhood->mNeighborCells( tElementIndex ))( iN )->get_index();

                tNeighborPhase = aCutIntegrationMesh->get_cell_bulk_phase( tNeighborIndex );

                auto tIter = tElementToLocalIndex.find( tNeighborIndex );

                // if this isn't a neighbor in the current flood fill set move on
                if ( tIter == tElementToLocalIndex.end() )
                {
                    continue;
                }
                else
                {
                    tNeighborOrd = tIter->second;
                }

                // If this is an neighbor element to include in the subdomain, has not already
                // been set and its phase matches the current elements phase then
                // add it to the active front and increment the count
                if ( tPhaseSet( 0, tNeighborOrd ) != 1 && tNeighborPhase == tPhaseIndex )
                {
                    tActiveFront( 0, tActiveFrontCount ) = tNeighborIndex;
                    tActiveFrontCount++;
                }
            }

            // Iterate through active front until there are no more elements in the active front
            // We start at the end of the front and work backwards
            while ( tActiveFrontCount != 0 )
            {
                // Current Element Index in the Active Front
                tActiveFrontElement = tActiveFront( 0, tActiveFrontCount - 1 );

                // Get Neighbor index from map if we're not considering the full domain
                tNeighborOrd = tElementToLocalIndex[tActiveFrontElement];

                // Get the neighbors phase
                tNeighborPhase = aCutIntegrationMesh->get_cell_bulk_phase( tActiveFrontElement );

                // If the neighbor phase matches our phase, then we add it's neighbor to the active front
                // Unless it has already been set
                if ( tNeighborPhase == tPhaseIndex && tPhaseSet( 0, tNeighborOrd ) != 1 )
                {
                    // Set the neighbor elements subphase value
                    tElementSubphase( 0, tNeighborOrd ) = tCurrentSubphase;

                    // Mark element as set
                    tPhaseSet( 0, tNeighborOrd ) = 1;

                    // Increase the number of phases set
                    tNumPhasesSet++;


                    // Add the elements other neighbors to the active front
                    bool tReplaced = false;
                    for ( moris::size_t i = 0; i < aCutNeighborhood->mNeighborCells( tActiveFrontElement )->size(); i++ )
                    {
                        tElementIndex = (*aCutNeighborhood->mNeighborCells( tActiveFrontElement ))( i )->get_index();

                        auto tIter = tElementToLocalIndex.find( tElementIndex );

                        // if this isn't a neighbor in the current flood fill set move on
                        if ( tIter == tElementToLocalIndex.end() )
                        {
                            continue;
                        }
                        else
                        {
                            tNeighborOrd = tIter->second;
                        }

                        // If this element is active and its phase hasn't been set
                        if ( tPhaseSet( 0, tNeighborOrd ) != 1 )
                        {
                            // and the previous element hasn't been replaced, then replace it
                            // don't add to active front count
                            if ( !tReplaced )
                            {
                                tActiveFront( 0, tActiveFrontCount - 1 ) = tElementIndex;
                                tReplaced                                = true;
                            }

                            // Else add to end of active front and add to the count
                            else
                            {

                                MORIS_ASSERT( tActiveFrontCount < tActiveFront.numel(), " Active front in flood fill not big enough" );
                                tActiveFront( 0, tActiveFrontCount ) = tElementIndex;
                                tActiveFrontCount++;
                            }
                        }
                    }
                }

                // Else if the phase doesn't match we remove that element from the active front
                else
                {
                    tActiveFront( 0, tActiveFrontCount - 1 ) = MORIS_INDEX_MAX;
                    tActiveFrontCount--;
                }
            }

            tCurrentSubphase++;
        }
    }

    aMaxValueAssigned = tCurrentSubphase - 1;

    return tElementSubphase;
}

// ----------------------------------------------------------------------------------

void
Integration_Mesh_Generator::construct_bulk_phase_to_bulk_phase_interface(
    Cut_Integration_Mesh*                                                aCutIntegrationMesh,
    moris::Cell< moris_index >&                                          aInterfaces,
    moris::Cell< moris::Cell< std::shared_ptr< IG_Cell_Side_Group > > >& aInterfaceBulkPhaseToBulk )
{
    Tracer tTracer( "XTK", "Integration_Mesh_Generator", "Bulk phase to bulk phase interface" ,mXTKModel->mVerboseLevel, 1  );
    // number of bulk phases
    moris::uint tNumBulkPhase = mGeometryEngine->get_num_bulk_phase();

    // allocate the output
    aInterfaceBulkPhaseToBulk.resize( tNumBulkPhase, moris::Cell< std::shared_ptr< IG_Cell_Side_Group > >( tNumBulkPhase, nullptr ) );

    std::shared_ptr< Facet_Based_Connectivity > tFaceConn = aCutIntegrationMesh->get_face_connectivity();

    for ( moris::uint iIF = 0; iIF < aInterfaces.size(); iIF++ )
    {
        moris_index tFacetIndex = aInterfaces( iIF );

        if ( tFaceConn->mFacetToCell( tFacetIndex ).size() == 2 )
        {

            moris_index tMyCellBulkPhase   = aCutIntegrationMesh->get_cell_bulk_phase( tFaceConn->mFacetToCell( tFacetIndex )( 0 )->get_index() );
            moris_index tYourCellBulkPhase = aCutIntegrationMesh->get_cell_bulk_phase( tFaceConn->mFacetToCell( tFacetIndex )( 1 )->get_index() );

            if ( tMyCellBulkPhase != (moris_index)tNumBulkPhase && tYourCellBulkPhase != (moris_index)tNumBulkPhase )
            {

                if ( aInterfaceBulkPhaseToBulk( tMyCellBulkPhase )( tYourCellBulkPhase ) == nullptr )
                {
                    aInterfaceBulkPhaseToBulk( tMyCellBulkPhase )( tYourCellBulkPhase ) = std::make_shared< IG_Cell_Side_Group >( aInterfaces.size() );
                }
                if ( aInterfaceBulkPhaseToBulk( tYourCellBulkPhase )( tMyCellBulkPhase ) == nullptr )
                {
                    aInterfaceBulkPhaseToBulk( tYourCellBulkPhase )( tMyCellBulkPhase ) = std::make_shared< IG_Cell_Side_Group >( aInterfaces.size() );
                }

                aInterfaceBulkPhaseToBulk( tMyCellBulkPhase )( tYourCellBulkPhase )->mIgCells.push_back( tFaceConn->mFacetToCell( tFacetIndex )( 0 ) );
                aInterfaceBulkPhaseToBulk( tMyCellBulkPhase )( tYourCellBulkPhase )->mIgCellSideOrdinals.push_back( tFaceConn->mFacetToCellEdgeOrdinal( tFacetIndex )( 0 ) );

                aInterfaceBulkPhaseToBulk( tYourCellBulkPhase )( tMyCellBulkPhase )->mIgCells.push_back( tFaceConn->mFacetToCell( tFacetIndex )( 1 ) );
                aInterfaceBulkPhaseToBulk( tYourCellBulkPhase )( tMyCellBulkPhase )->mIgCellSideOrdinals.push_back( tFaceConn->mFacetToCellEdgeOrdinal( tFacetIndex )( 1 ) );
            }
        }
        else
        {
            // This else case happens when an interface is coincident with the boundary of the domain
            std::cout << "Warning: interface case not handled " << std::endl;
        }
    }
}

// ----------------------------------------------------------------------------------

void
Integration_Mesh_Generator::construct_bulk_phase_to_bulk_phase_dbl_side_interface(
    Cut_Integration_Mesh*                                                       aCutIntegrationMesh,
    moris::Cell< moris_index >&                                                 aInterfaces,
    moris::Cell< moris::Cell< std::shared_ptr< IG_Cell_Double_Side_Group > > >& aDblSideInterfaceBulkPhaseToBulk )
{
    Tracer tTracer( "XTK", "Integration_Mesh_Generator", "Bulk phase to bulk phase double sided interface",mXTKModel->mVerboseLevel, 1  );

    // number of bulk phases
    moris::uint tNumBulkPhase = mGeometryEngine->get_num_bulk_phase();

    // allocate the output
    aDblSideInterfaceBulkPhaseToBulk.resize( tNumBulkPhase, moris::Cell< std::shared_ptr< IG_Cell_Double_Side_Group > >( tNumBulkPhase, nullptr ) );

    std::shared_ptr< Facet_Based_Connectivity > tFaceConn = aCutIntegrationMesh->get_face_connectivity();

    // make sure a Facet Connectivity has been assigned to Cut integration mesh
    MORIS_ERROR( tFaceConn != nullptr, 
        "Integration_Mesh_Generator::construct_bulk_phase_to_bulk_phase_dbl_side_interface() - No Face Connectivity on Cut Integration Mesh." );

    for ( moris::uint iIF = 0; iIF < aInterfaces.size(); iIF++ )
    {
        moris_index tFacetIndex = aInterfaces( iIF );

        if ( tFaceConn->mFacetToCell( tFacetIndex ).size() == 2 )
        {

            moris_index tMyCellBulkPhase   = aCutIntegrationMesh->get_cell_bulk_phase( tFaceConn->mFacetToCell( tFacetIndex )( 0 )->get_index() );
            moris_index tYourCellBulkPhase = aCutIntegrationMesh->get_cell_bulk_phase( tFaceConn->mFacetToCell( tFacetIndex )( 1 )->get_index() );

            if ( aDblSideInterfaceBulkPhaseToBulk( tMyCellBulkPhase )( tYourCellBulkPhase ) == nullptr )
            {
                aDblSideInterfaceBulkPhaseToBulk( tMyCellBulkPhase )( tYourCellBulkPhase ) = std::make_shared< IG_Cell_Double_Side_Group >( 1 );
            }
            if ( aDblSideInterfaceBulkPhaseToBulk( tYourCellBulkPhase )( tMyCellBulkPhase ) == nullptr )
            {
                aDblSideInterfaceBulkPhaseToBulk( tYourCellBulkPhase )( tMyCellBulkPhase ) = std::make_shared< IG_Cell_Double_Side_Group >( 1 );
            }

            aDblSideInterfaceBulkPhaseToBulk( tMyCellBulkPhase )( tYourCellBulkPhase )->mLeaderIgCells.push_back( tFaceConn->mFacetToCell( tFacetIndex )( 0 ) );
            aDblSideInterfaceBulkPhaseToBulk( tMyCellBulkPhase )( tYourCellBulkPhase )->mLeaderIgCellSideOrdinals.push_back( tFaceConn->mFacetToCellEdgeOrdinal( tFacetIndex )( 0 ) );
            aDblSideInterfaceBulkPhaseToBulk( tMyCellBulkPhase )( tYourCellBulkPhase )->mFollowerIgCells.push_back( tFaceConn->mFacetToCell( tFacetIndex )( 1 ) );
            aDblSideInterfaceBulkPhaseToBulk( tMyCellBulkPhase )( tYourCellBulkPhase )->mFollowerIgCellSideOrdinals.push_back( tFaceConn->mFacetToCellEdgeOrdinal( tFacetIndex )( 1 ) );

            aDblSideInterfaceBulkPhaseToBulk( tYourCellBulkPhase )( tMyCellBulkPhase )->mLeaderIgCells.push_back( tFaceConn->mFacetToCell( tFacetIndex )( 1 ) );
            aDblSideInterfaceBulkPhaseToBulk( tYourCellBulkPhase )( tMyCellBulkPhase )->mLeaderIgCellSideOrdinals.push_back( tFaceConn->mFacetToCellEdgeOrdinal( tFacetIndex )( 1 ) );
            aDblSideInterfaceBulkPhaseToBulk( tYourCellBulkPhase )( tMyCellBulkPhase )->mFollowerIgCells.push_back( tFaceConn->mFacetToCell( tFacetIndex )( 0 ) );
            aDblSideInterfaceBulkPhaseToBulk( tYourCellBulkPhase )( tMyCellBulkPhase )->mFollowerIgCellSideOrdinals.push_back( tFaceConn->mFacetToCellEdgeOrdinal( tFacetIndex )( 0 ) );
        }
        else
        {
            std::cout << "WARNING, need to handle the single sided interface case here" << std::endl;
        }
    }
}

// ----------------------------------------------------------------------------------

std::string
Integration_Mesh_Generator::get_interface_side_set_name(
    moris_index aGeomIndex,
    moris_index aBulkPhaseIndex0,
    moris_index aBulkPhaseIndex1 )
{
    MORIS_ASSERT( aGeomIndex < (moris_index)mGeometryEngine->get_num_geometries(), "Geometry index out of bounds" );
    MORIS_ASSERT( aBulkPhaseIndex0 < (moris_index)mGeometryEngine->get_num_bulk_phase(), "Bulk phase index 0 out of bounds" );
    MORIS_ASSERT( aBulkPhaseIndex1 < (moris_index)mGeometryEngine->get_num_bulk_phase(), "Bulk phase index 1 out of bounds" );

    return "iside_b0_" + std::to_string( aBulkPhaseIndex0 ) + "_b1_" + std::to_string( aBulkPhaseIndex1 );
}

// ----------------------------------------------------------------------------------

void
Integration_Mesh_Generator::construct_interface_sets(
    Cut_Integration_Mesh*                                                aCutIntegrationMesh,
    moris::Cell< moris::Cell< std::shared_ptr< IG_Cell_Side_Group > > >& aInterfaceBulkPhaseToBulk )
{
    Tracer tTracer( "XTK", "Integration_Mesh_Generator", "Construct interface sets" ,mXTKModel->mVerboseLevel, 1  );
    // determine the side set names
    uint tNumBulkPhases = mGeometryEngine->get_num_bulk_phase();

    Cell< std::string > tInterfaceSideNames;
    tInterfaceSideNames.reserve( tNumBulkPhases * tNumBulkPhases );

    Cell< std::shared_ptr< IG_Cell_Side_Group > > tCellSideGroups;
    tCellSideGroups.reserve( tNumBulkPhases * tNumBulkPhases );

    for ( moris::moris_index iP0 = 0; iP0 < (moris_index)tNumBulkPhases; iP0++ )
    {
        for ( moris::moris_index iP1 = 0; iP1 < (moris_index)tNumBulkPhases; iP1++ )
        {
            if ( iP1 != iP0 )
            {
                std::string tInterfaceSideSetName = this->get_interface_side_set_name( 0, iP0, iP1 );

                tInterfaceSideNames.push_back( tInterfaceSideSetName );

                tCellSideGroups.push_back( aInterfaceBulkPhaseToBulk( iP0 )( iP1 ) );
            }
        }
    }

    Cell< moris_index > tSideSetOrds = aCutIntegrationMesh->register_side_set_names( tInterfaceSideNames );

    // iterate and add the side set groups to cut mehs
    for ( moris::uint iSet = 0; iSet < tSideSetOrds.size(); iSet++ )
    {
        aCutIntegrationMesh->mSideSetCellSides( tSideSetOrds( iSet ) ) = tCellSideGroups( iSet );
    }
}

// ----------------------------------------------------------------------------------

void
Integration_Mesh_Generator::construct_bulk_phase_cell_groups(
    Cut_Integration_Mesh*                            aCutIntegrationMesh,
    moris::Cell< std::shared_ptr< IG_Cell_Group > >& aBulkPhaseCellGroups )
{
    Tracer tTracer( "XTK", "Integration_Mesh_Generator", "Construct bulk phase cell groups" ,mXTKModel->mVerboseLevel, 1  );
    uint   tNumBulkPhases = mGeometryEngine->get_num_bulk_phase() + 1;// 1 for the err block

    aBulkPhaseCellGroups.resize( tNumBulkPhases, nullptr );

    // iterate through the controlled integration cells to avoid inactive bg cells
    for ( moris::uint iCell = 0; iCell < aCutIntegrationMesh->mControlledIgCells.size(); iCell++ )
    {
        moris::mtk::Cell* tCell = aCutIntegrationMesh->mControlledIgCells( iCell ).get();
        // cell bulk phase
        moris_index tBulkPhase = aCutIntegrationMesh->get_cell_bulk_phase( tCell->get_index() );

        if ( tBulkPhase != MORIS_INDEX_MAX )
        {
            if ( aBulkPhaseCellGroups( tBulkPhase ) == nullptr )
            {
                aBulkPhaseCellGroups( tBulkPhase ) = std::make_shared< IG_Cell_Group >( 0 );
            }
            aBulkPhaseCellGroups( tBulkPhase )->mIgCellGroup.push_back( &aCutIntegrationMesh->get_mtk_cell( tCell->get_index() ) );
        }
    }
}

// ----------------------------------------------------------------------------------

void
Integration_Mesh_Generator::construct_bulk_phase_blocks(
    Cut_Integration_Mesh*                            aCutIntegrationMesh,
    moris::Cell< std::shared_ptr< IG_Cell_Group > >& aBulkPhaseCellGroups )
{

    Tracer      tTracer( "XTK", "Integration_Mesh_Generator", "Construct bulk phase blocks",mXTKModel->mVerboseLevel, 1  );
    std::string tBlockBaseName = "cutblock";

    // determine the side set names
    uint tNumBulkPhases = mGeometryEngine->get_num_bulk_phase();

    MORIS_ERROR( aBulkPhaseCellGroups.size() == tNumBulkPhases + 1, "We expect there to be a bulk phase cell group for each bulk phase and 1 for problematic cells" );

    Cell< std::string > tBlockSetNames( tNumBulkPhases + 1 );
    
    // decide on cell topology of integration elements based on number of spatial dimensions
    enum CellTopology tCellTopo = xtk::determine_cell_topology( 
        aCutIntegrationMesh->get_spatial_dim(), 
        aCutIntegrationMesh->mXTKModel->ig_element_order(), 
        CellShape::SIMPLEX );

    // iterate through and construct the names of the blocks
    for ( moris::uint iBP = 0; iBP < tNumBulkPhases; iBP++ )
    {
        tBlockSetNames( iBP ) = tBlockBaseName + "_p_" + std::to_string( iBP );
    }

    tBlockSetNames( tNumBulkPhases ) = "err_block";

    Cell< moris_index > tBlockSetOrds = aCutIntegrationMesh->register_block_set_names( tBlockSetNames, tCellTopo );

    // iterate and add the side set groups to cut mehs
    for ( moris::uint iSet = 0; iSet < tBlockSetOrds.size(); iSet++ )
    {
        aCutIntegrationMesh->mBlockSetCellGroup( tBlockSetOrds( iSet ) ) = aBulkPhaseCellGroups( iSet );
    }
}

// ----------------------------------------------------------------------------------

moris_index
Integration_Mesh_Generator::edge_exists(
    moris::Cell< moris::mtk::Vertex* >&                aVerticesOnEdge,
    std::unordered_map< moris_index, moris_index >&    aLocaLVertexMap,
    moris::Cell< moris::Cell< uint > >&                aVertexToEdge,
    moris::Cell< moris::Cell< moris::mtk::Vertex* > >& aFullEdgeVertices )
{
    moris_index tEdgeIndex = MORIS_INDEX_MAX;

    // get them in order based on id
    std::sort( aVerticesOnEdge.data().begin(), aVerticesOnEdge.data().end(), moris::comparePtrToVertexIdBased );

    // local vertex index
    auto tIter = aLocaLVertexMap.find( aVerticesOnEdge( 0 )->get_index() );

    MORIS_ERROR( tIter != aLocaLVertexMap.end(), "Invalid vertex detected." );

    moris_index tLocalVertexIndex = tIter->second;

    // iterate through edges attached to the first vertex (depend on ascending order)
    for ( moris::uint iEdge = 0; iEdge < aVertexToEdge( tLocalVertexIndex ).size(); iEdge++ )
    {
        moris_index tEdgeIndex = aVertexToEdge( tLocalVertexIndex )( iEdge );

        MORIS_ASSERT( aFullEdgeVertices( tEdgeIndex )( 0 )->get_index() == aVerticesOnEdge( 0 )->get_index(), 
                "Numbering issues, edges should be in ascending order based on vertex id" );

        // check the second vertex on the edge
        if ( aFullEdgeVertices( tEdgeIndex )( 1 )->get_index() == aVerticesOnEdge( 1 )->get_index() )
        {
            return tEdgeIndex;
        }
    }
    return tEdgeIndex;
}

// ----------------------------------------------------------------------------------

moris_index
Integration_Mesh_Generator::facet_exists(
    moris::Cell< moris::mtk::Vertex* >&                aVerticesOnFacet,
    std::unordered_map< moris_index, moris_index >&    aLocaLVertexMap,
    moris::Cell< moris::Cell< uint > >&                aVertexToFacet,
    moris::Cell< moris::Cell< moris::mtk::Vertex* > >& aFullFacetVertices )
{
    moris_index tFacetIndex = MORIS_INDEX_MAX;

    // number of vertices per facet
    // moris::uint tNumVerticesPerFacet = aVerticesOnFacet.size();

    // get them in order based on id
    std::sort( aVerticesOnFacet.data().begin(), aVerticesOnFacet.data().end(), moris::comparePtrToVertexIdBased );

    // local vertex index
    auto tIter = aLocaLVertexMap.find( aVerticesOnFacet( 0 )->get_index() );
    MORIS_ERROR( tIter != aLocaLVertexMap.end(), "Invalid vertex detected." );

    moris_index tLocalVertexIndex = tIter->second;

    // iterate through edges attached to the first vertex (depend on ascending order)
    for ( moris::uint iFacet = 0; iFacet < aVertexToFacet( tLocalVertexIndex ).size(); iFacet++ )
    {
        tFacetIndex = aVertexToFacet( tLocalVertexIndex )( iFacet );

        MORIS_ASSERT( aFullFacetVertices( tFacetIndex )( 0 )->get_index() == aVerticesOnFacet( 0 )->get_index(), "Numbering issue" );

        // figure out if this is the same vertex
        bool tFoundFacet = std::equal( aFullFacetVertices( tFacetIndex ).begin(), aFullFacetVertices( tFacetIndex ).end(), aVerticesOnFacet.begin() );

        // for(moris::uint iV =1; iV<tNumVerticesPerFacet; iV++)
        // {
        //     if(aFullFacetVertices(tFacetIndex)(iV)->get_index() != aVerticesOnFacet(iV)->get_index())
        //     {
        //         tFoundFacet = false;
        //     }
        // }

        if ( tFoundFacet )
        {
            return tFacetIndex;
            break;
        }
    }
    return MORIS_INDEX_MAX;
}

// ----------------------------------------------------------------------------------

void
Integration_Mesh_Generator::create_facet_from_element_to_node(
    moris::Cell< moris::mtk::Cell* >&           aCells,
    std::shared_ptr< Facet_Based_Connectivity > aFaceConnectivity )
{
    Tracer tTracer( "XTK", "Integration_Mesh_Generator", "Creating Facets",mXTKModel->mVerboseLevel, 1 );

    // Note: this function assumes that all cells are of same type, i.e. have the same cell info

    if ( aCells.size() > 0 )
    {
        // cell information
        moris::mtk::Cell_Info const* tCellInfo = aCells( 0 )->get_cell_info();

        // node to facet map
        moris::Matrix< moris::IdMat > tElementToFacetMap = tCellInfo->get_node_to_facet_map();

        // list of unique vertices
        moris::uint                                    tNumNodes = 0;
        moris::Cell< moris::mtk::Vertex* >             tVertices;
        std::unordered_map< moris_index, moris_index > tVertexIndexToLocalIndexMap;

        // reserve memory for list of vertex points
        tVertices.reserve( aCells.size() );

        for ( moris::uint i = 0; i < aCells.size(); i++ )
        {
            moris::Cell< moris::mtk::Vertex* > tCellVerts = aCells( i )->get_vertex_pointers();

            for ( moris::uint iV = 0; iV < tCellVerts.size(); iV++ )
            {
                // check vertex has already been assigned local vertex index
                if ( tVertexIndexToLocalIndexMap.find( tCellVerts( iV )->get_index() ) == tVertexIndexToLocalIndexMap.end() )
                {
                    // add local vertex index to map
                    tVertexIndexToLocalIndexMap[tCellVerts( iV )->get_index()] = (moris_index)tNumNodes;

                    // store vertex pointer
                    tVertices.push_back( tCellVerts( iV ) );

                    // increase node counter
                    tNumNodes++;
                }
            }
        }

        // Maximum faces per node
        moris::uint tNumFacesPerElem  = tCellInfo->get_num_facets();
        moris::uint tNumElements      = aCells.size();
        moris::uint tMaxNumFacets     = tNumElements * tNumFacesPerElem;
        moris::uint tNumNodesPerFacet = tElementToFacetMap.n_cols();

        // set size for cell to facet map
        aFaceConnectivity->mCellToFacet.resize( aCells.size() );

        // Define by which the capacity of edge-based cells is increased
         uint tIncNumFacets = tMaxNumFacets / 10;

        // reserve memory for facet-based cells
        aFaceConnectivity->mFacetToCell.reserve( tIncNumFacets );
        aFaceConnectivity->mFacetToCellEdgeOrdinal.reserve( tIncNumFacets );
        aFaceConnectivity->mFacetVertices.reserve( tIncNumFacets );

        // define local list
        moris::Cell< moris::Cell< uint > > tVertexToFacetIndex( tNumNodes );
        moris::Cell< moris::mtk::Vertex* > tVerticesOnFacet( tNumNodesPerFacet, nullptr );

        // iterate through cells
        for ( moris::uint i = 0; i < aCells.size(); i++ )
        {
            // check that cell exits only once
            MORIS_ERROR( aFaceConnectivity->mCellIndexToCellOrdinal.find( aCells( i )->get_index() ) == aFaceConnectivity->mCellIndexToCellOrdinal.end(),
                    "Duplicate cell in aCells provided" );

            // set local cell index in map
            aFaceConnectivity->mCellIndexToCellOrdinal[aCells( i )->get_index()] = i;

            // get list of vertex pointers of cell
            moris::Cell< moris::mtk::Vertex* > tCellVerts = aCells( i )->get_vertex_pointers();

            // reserve memory for storing facets on cell
            aFaceConnectivity->mCellToFacet( i ).reserve( tNumFacesPerElem );

            // iterate through edges of cell
            for ( moris::uint iEdge = 0; iEdge < tElementToFacetMap.n_rows(); iEdge++ )
            {
                // get the vertices on the edge
                for ( moris::uint iVOnE = 0; iVOnE < tElementToFacetMap.n_cols(); iVOnE++ )
                {
                    tVerticesOnFacet( iVOnE ) = tCellVerts( tElementToFacetMap( iEdge, iVOnE ) );
                }

                // figure out if the edge exists and if so where
                moris_index tFacetIndex = this->facet_exists(
                        tVerticesOnFacet,
                        tVertexIndexToLocalIndexMap,
                        tVertexToFacetIndex,
                        aFaceConnectivity->mFacetVertices );

                // add new facet if it doesn't exist already
                if ( tFacetIndex == MORIS_INDEX_MAX )
                {
                    // set edge index
                    tFacetIndex = aFaceConnectivity->mFacetVertices.size();

                    // store pointers of vertices on facet
                    aFaceConnectivity->mFacetVertices.push_back( tVerticesOnFacet, tIncNumFacets );

                    // for current facet initialize list to store cells and cell ordinals
                    aFaceConnectivity->mFacetToCell.push_back( moris::Cell< moris::mtk::Cell* >(), tIncNumFacets );
                    aFaceConnectivity->mFacetToCellEdgeOrdinal.push_back( moris::Cell< moris::moris_index >(), tIncNumFacets );

                    // reserve memory for list storing cells and their side ordinal connected to edge
                    // guess of cell size: number of facets per element
                    aFaceConnectivity->mFacetToCell.back().reserve( tNumFacesPerElem );
                    aFaceConnectivity->mFacetToCellEdgeOrdinal.back().reserve( tNumFacesPerElem );

                    // store facet index with vertex
                    auto tIter = tVertexIndexToLocalIndexMap.find( tVerticesOnFacet( 0 )->get_index() );
                    MORIS_ERROR( tIter != tVertexIndexToLocalIndexMap.end(), "Invalid vertex detected." );
                    moris_index tLocalVertexIndex = tIter->second;

                    tVertexToFacetIndex( tLocalVertexIndex ).push_back( tFacetIndex );
                }

                // store facet index on cell
                aFaceConnectivity->mCellToFacet( i ).push_back( tFacetIndex );

                // store cell and cell ordinal with facet
                // if needed increase cell capacity by increments of number of faces per element
                aFaceConnectivity->mFacetToCell( tFacetIndex ).push_back( aCells( i ), tNumFacesPerElem );
                aFaceConnectivity->mFacetToCellEdgeOrdinal( tFacetIndex ).push_back( iEdge, tNumFacesPerElem );
            }
        }

        // trim inner and outer cells
        shrink_to_fit_all(aFaceConnectivity->mFacetVertices);
        shrink_to_fit_all(aFaceConnectivity->mFacetToCell);
        shrink_to_fit_all(aFaceConnectivity->mFacetToCellEdgeOrdinal);
    }
}

// ----------------------------------------------------------------------------------

void
Integration_Mesh_Generator::generate_cell_neighborhood(
    moris::Cell< moris::mtk::Cell* >&                 aCells,
    std::shared_ptr< Facet_Based_Connectivity >       aFaceConnectivity,
    std::shared_ptr< Cell_Neighborhood_Connectivity > aNeighborhood )
{
    Tracer tTracer( "XTK", "Integration_Mesh_Generator", "Generate Neighborhood" ,mXTKModel->mVerboseLevel, 1  );

    // Initialize Sizes and Variables used in routine
    moris_index tMaxIndex            = this->get_max_index( aCells );
    moris_index tMaxNumElementToFace = 2;

    // Initialize Element to Element with size number of elements x number of faces per element filled with a dummy value.
    aNeighborhood->mNeighborCells.resize( tMaxIndex + 1 );
    aNeighborhood->mMySideOrdinal.resize( tMaxIndex + 1 );
    aNeighborhood->mNeighborSideOrdinal.resize( tMaxIndex + 1 );

    for ( moris::uint i = 0; i < (moris::uint)tMaxIndex + 1; i++ )
    {
        aNeighborhood->mNeighborCells( i )       = std::make_shared< moris::Cell< moris::mtk::Cell* > >();
        aNeighborhood->mMySideOrdinal( i )       = std::make_shared< moris::Cell< moris_index > >();
        aNeighborhood->mNeighborSideOrdinal( i ) = std::make_shared< moris::Cell< moris_index > >();

        aNeighborhood->mNeighborCells( i )->reserve( tMaxNumElementToFace );
        aNeighborhood->mMySideOrdinal( i )->reserve( tMaxNumElementToFace );
        aNeighborhood->mNeighborSideOrdinal( i )->reserve( tMaxNumElementToFace );
    }

    // iterate through facet to cell connectivity. any cell that shares a facet is a neighbor
    for ( moris::uint iF = 0; iF < aFaceConnectivity->mFacetToCell.size(); iF++ )
    {
        // iterate through cells attached to this facet (either just 1 or 2)
        MORIS_ASSERT( aFaceConnectivity->mFacetToCell( iF ).size() == 1 || aFaceConnectivity->mFacetToCell( iF ).size() == 2,
                "Facet should either connect to no cell or one other cell" );

        if ( aFaceConnectivity->mFacetToCell( iF ).size() == 2 )
        {
            moris_index tCellIndex0 = aFaceConnectivity->mFacetToCell( iF )( 0 )->get_index();
            moris_index tCellIndex1 = aFaceConnectivity->mFacetToCell( iF )( 1 )->get_index();

            aNeighborhood->mNeighborCells( tCellIndex0 )->push_back( aFaceConnectivity->mFacetToCell( iF )( 1 ) );
            aNeighborhood->mNeighborCells( tCellIndex1 )->push_back( aFaceConnectivity->mFacetToCell( iF )( 0 ) );

            aNeighborhood->mMySideOrdinal( tCellIndex0 )->push_back( aFaceConnectivity->mFacetToCellEdgeOrdinal( iF )( 0 ) );
            aNeighborhood->mMySideOrdinal( tCellIndex1 )->push_back( aFaceConnectivity->mFacetToCellEdgeOrdinal( iF )( 1 ) );

            aNeighborhood->mNeighborSideOrdinal( tCellIndex0 )->push_back( aFaceConnectivity->mFacetToCellEdgeOrdinal( iF )( 1 ) );
            aNeighborhood->mNeighborSideOrdinal( tCellIndex0 )->push_back( aFaceConnectivity->mFacetToCellEdgeOrdinal( iF )( 0 ) );
        }
    }
}

// ----------------------------------------------------------------------------------

void
Integration_Mesh_Generator::create_edges_from_element_to_node(
    moris::Cell< moris::mtk::Cell* >           aCells,
    std::shared_ptr< Edge_Based_Connectivity > aEdgeConnectivity )
{
    Tracer tTracer( "XTK", "Integration_Mesh_Generator", "Creating Edges" ,mXTKModel->mVerboseLevel, 1  );

    // Note: this function assumes that all cells are of same type, i.e. have the same cell info

    if ( aCells.size() > 0 )
    {
        // cell information
        moris::mtk::Cell_Info const* tCellInfo = aCells( 0 )->get_cell_info();

        // number of edges per element
        moris::uint tNumEdgePerElem = tCellInfo->get_num_edges();

        // nodes to edge map
        moris::Matrix< moris::IdMat > tElementEdgeToNodeMap = tCellInfo->get_node_to_edge_map();

        // number of cells
        const moris::uint tNumElements = aCells.size();

        // number of nodes per edge
        const moris::uint tNumNodesPerEdge = 2;

        MORIS_ERROR( tNumNodesPerEdge == 2,
                "Only works on two node edges at the moment" );

        MORIS_ERROR( tNumNodesPerEdge == tElementEdgeToNodeMap.n_cols(),
                "Mismatch in number of nodes per edge (only operating on two node edges)" );

        // initialize counter for number of unique vertices/nodes
        moris::uint                                    tNumNodes = 0;
        moris::Cell< moris::mtk::Vertex* >             tVertices;

        // reserve memory for cell of vertex pointers
        tVertices.reserve( aCells.size() );

        // initialize map, pairing global (?) index of vertex to new numbering going through unique nodes
        std::unordered_map< moris_index, moris_index > tVertexIndexToLocalIndexMap;

        // loop over all cells/elements and vertices on them
        for ( moris::uint i = 0; i < aCells.size(); i++ )
        {
            moris::Cell< moris::mtk::Vertex* > tCellVerts = aCells( i )->get_vertex_pointers();

            for ( moris::uint iV = 0; iV < tCellVerts.size(); iV++ )
            {
                // check if vertex has already been given a new index
                if ( tVertexIndexToLocalIndexMap.find( tCellVerts( iV )->get_index() ) == tVertexIndexToLocalIndexMap.end() )
                {
                    // give new unique vertex new index
                    tVertexIndexToLocalIndexMap[tCellVerts( iV )->get_index()] = (moris_index)tNumNodes;

                    // store vertex pointer
                    tVertices.push_back( tCellVerts( iV ) );

                    // increment counter of unique nodes
                    tNumNodes++;
                }
            }
        }

        // Set size for cell to edge map
        aEdgeConnectivity->mCellToEdge.resize( aCells.size() );

        // Maximum number of edges
        uint tMaxNumEdges = tNumElements * tNumEdgePerElem;

        // Define by which the capacity of edge-based cells is increased
        uint tIncNumEdges = tMaxNumEdges * tNumNodesPerEdge / 10;

        // Reserve memory of edge-based cells
        aEdgeConnectivity->mEdgeToCell.reserve( tIncNumEdges );
        aEdgeConnectivity->mEdgeToCellEdgeOrdinal.reserve( tIncNumEdges );
        aEdgeConnectivity->mEdgeVertices.reserve( tIncNumEdges );

        // Set size of cell storing edges connected to node
        moris::Cell< moris::Cell< uint > > tVertexToEdgeIndex( tNumNodes );

        // Set size of auxiliary cell with vertex pointers
        moris::Cell< moris::mtk::Vertex* > tVerticesOnEdge( tNumNodesPerEdge, nullptr );

        // loop over all cells and their edges
        for ( moris::uint i = 0; i < aCells.size(); i++ )
        {
            // get pointers to vertices of cell
            moris::Cell< moris::mtk::Vertex* > tCellVerts = aCells( i )->get_vertex_pointers();

            // reserve memory for storing edges on cell
            aEdgeConnectivity->mCellToEdge( i ).reserve( tNumEdgePerElem );

            // iterate through facets of cell
            for ( moris::uint iEdge = 0; iEdge < tElementEdgeToNodeMap.n_rows(); iEdge++ )
            {
                // get the vertices on the edge
                for ( moris::uint iVOnE = 0; iVOnE < tElementEdgeToNodeMap.n_cols(); iVOnE++ )
                {
                    tVerticesOnEdge( iVOnE ) = tCellVerts( tElementEdgeToNodeMap( iEdge, iVOnE ) );
                }

                // figure out if the edge exists and if so where (if it doesn't exist it returns MORIS_INDEX_MAX as an index )
                moris_index tEdgeIndex = this->edge_exists(
                        tVerticesOnEdge,
                        tVertexIndexToLocalIndexMap,
                        tVertexToEdgeIndex,
                        aEdgeConnectivity->mEdgeVertices );

                // add new edge new if it doesn't exist already
                if ( tEdgeIndex == MORIS_INDEX_MAX )
                {
                    // set edge index
                    tEdgeIndex = aEdgeConnectivity->mEdgeVertices.size();

                    // store pointers of vertices on edge
                    aEdgeConnectivity->mEdgeVertices.push_back( tVerticesOnEdge, tIncNumEdges );

                    // for current edge initialize list to store connected cell and side ordinal
                    aEdgeConnectivity->mEdgeToCell.push_back( moris::Cell< moris::mtk::Cell* >( ), tIncNumEdges );
                    aEdgeConnectivity->mEdgeToCellEdgeOrdinal.push_back( moris::Cell< moris::moris_index >( ), tIncNumEdges );

                    // reserve memory for list storing cells and their side ordinal connected to edge
                    // guess of cell size: number of edges per element
                    aEdgeConnectivity->mEdgeToCell.back().reserve( tNumEdgePerElem );
                    aEdgeConnectivity->mEdgeToCellEdgeOrdinal.back().reserve( tNumEdgePerElem );

                    // store edge with vertex
                    auto tIter = tVertexIndexToLocalIndexMap.find( tVerticesOnEdge( 0 )->get_index() );

                    MORIS_ERROR( tIter != tVertexIndexToLocalIndexMap.end(), "Invalid vertex detected." );

                    moris_index tLocalVertexIndex = tIter->second;

                    // store edge index for vertex; if needed increase cell capacity by increments of 10
                    tVertexToEdgeIndex( tLocalVertexIndex ).push_back( tEdgeIndex, 10 );
                }

                // add edge index to cell
                aEdgeConnectivity->mCellToEdge( i ).push_back( tEdgeIndex );

                // add cell pointer and side ordinal to edge
                // if needed increase cell capacity by increments of number of edges per element
                aEdgeConnectivity->mEdgeToCell( tEdgeIndex ).push_back( aCells( i ), tNumEdgePerElem );
                aEdgeConnectivity->mEdgeToCellEdgeOrdinal( tEdgeIndex ).push_back( iEdge, tNumEdgePerElem );
            }
        }

        // trim outer and inner cells
        shrink_to_fit_all( aEdgeConnectivity->mEdgeVertices );
        shrink_to_fit_all( aEdgeConnectivity->mEdgeToCell );
        shrink_to_fit_all( aEdgeConnectivity->mEdgeToCellEdgeOrdinal );
    }
}

// ----------------------------------------------------------------------------------

void
Integration_Mesh_Generator::deduce_facet_ancestry(
    Cut_Integration_Mesh*                       aCutIntegrationMesh,
    moris::mtk::Mesh*                           aBackgroundMesh,
    std::shared_ptr< Facet_Based_Connectivity > aIgCellGroupFacetConnectivity,
    moris::Cell< moris::mtk::Cell* > const&     aParentCellForDeduction,
    std::shared_ptr< Facet_Based_Ancestry >     aIgFacetAncestry )
{
    Tracer tTracer( "XTK", "Integration_Mesh_Generator", "Facet Ancestry" ,mXTKModel->mVerboseLevel, 1  );

    // make sure we are starting clean
    aIgFacetAncestry->mFacetParentEntityIndex.clear();
    aIgFacetAncestry->mFacetParentEntityRank.clear();
    aIgFacetAncestry->mFacetParentEntityOrdinalWrtBackgroundCell.clear();

    // number of edges in the edge connectivity
    moris::uint tNumFacets = aIgCellGroupFacetConnectivity->mFacetVertices.size();

    // allocate the data in the edge ancestry
    aIgFacetAncestry->mFacetParentEntityIndex.resize( tNumFacets );
    aIgFacetAncestry->mFacetParentEntityRank.resize( tNumFacets );
    aIgFacetAncestry->mFacetParentEntityOrdinalWrtBackgroundCell.resize( tNumFacets );

    Cell< std::shared_ptr< Matrix< IndexMat > > > tEntityConnectedToParent( 4 );
    for ( moris::uint iInit = 0; iInit < 4; iInit++ )
    {
        tEntityConnectedToParent( iInit ) = std::make_shared< Matrix< IndexMat > >( 0, 0 );
    }

    Matrix< IndexMat > tActiveConnectivity( 1, 4, MORIS_INDEX_MAX );

    moris::Cell< moris_index > tCandidateFacetOrds;

    moris::Cell< moris::moris_index > tFacetVertexParentInds;
    moris::Cell< moris::moris_index > tFacetVertexParentRanks;
    moris::Cell< moris::moris_index > tFacetVertexParentOrds;


    // iterate through edges in the edge connectivity
    for ( moris::uint iFacet = 0; iFacet < tNumFacets; iFacet++ )
    {

        tCandidateFacetOrds.clear();
        tActiveConnectivity.fill( MORIS_INDEX_MAX );

        moris::mtk::Cell_Info const* tParentCellInfo = aParentCellForDeduction( iFacet )->get_cell_info();

        // vertices of the facet
        moris::Cell< moris::mtk::Vertex* > const& tFacetVertices = aIgCellGroupFacetConnectivity->mFacetVertices( iFacet );

        // get the parent of these vertices from the mesh
        tFacetVertexParentInds.resize( tFacetVertices.size() );
        tFacetVertexParentRanks.resize( tFacetVertices.size() );
        tFacetVertexParentOrds.resize( tFacetVertices.size() );

        // all we need to figure out is if the facet is a subfacet or a interior to the cell in this case.

        // iterate through vertices
        for ( moris::uint iV = 0; iV < tFacetVertices.size(); iV++ )
        {
            tFacetVertexParentInds( iV )  = aCutIntegrationMesh->mIgVertexParentEntityIndex( tFacetVertices( iV )->get_index() );
            tFacetVertexParentRanks( iV ) = aCutIntegrationMesh->mIgVertexParentEntityRank( tFacetVertices( iV )->get_index() );

            // iterior background cell is the parent
            if ( tFacetVertexParentRanks( iV ) == 3 )
            {
                // mark the entity parent as the deduction cell
                aIgFacetAncestry->mFacetParentEntityIndex( iFacet )                    = aParentCellForDeduction( iFacet )->get_index();
                aIgFacetAncestry->mFacetParentEntityRank( iFacet )                     = 3;
                aIgFacetAncestry->mFacetParentEntityOrdinalWrtBackgroundCell( iFacet ) = 0;
                tCandidateFacetOrds.clear();
                break;
            }

            if ( tActiveConnectivity( tFacetVertexParentRanks( iV ) ) == MORIS_INDEX_MAX )
            {
                ( *tEntityConnectedToParent( tFacetVertexParentRanks( iV ) ) ) =
                    aBackgroundMesh->get_entity_connected_to_entity_loc_inds(
                        aParentCellForDeduction( iFacet )->get_index(),
                        EntityRank::ELEMENT,
                        moris::get_entity_rank_from_index( tFacetVertexParentRanks( iV ) ) );
            }

            // figure out the ordinal
            for ( moris::uint iE = 0; iE < tEntityConnectedToParent( tFacetVertexParentRanks( iV ) )->numel(); iE++ )
            {
                if ( (*tEntityConnectedToParent( tFacetVertexParentRanks( iV ) ))( iE ) == tFacetVertexParentInds( iV ) )
                {
                    tFacetVertexParentOrds( iV ) = iE;
                    break;
                }
            }

            // iterate through facets of cell
            if ( iV == 0 )
            {
                for ( moris::uint iCellFacets = 0; iCellFacets < tParentCellInfo->get_num_facets(); iCellFacets++ )
                {
                    if ( tParentCellInfo->is_entity_connected_to_facet( (moris_index)iCellFacets, tFacetVertexParentOrds( iV ), tFacetVertexParentRanks( iV ) ) )
                    {
                        tCandidateFacetOrds.push_back( iCellFacets );
                    }
                }
            }
            else
            {
                for ( auto i = tCandidateFacetOrds.data().rbegin(); i != tCandidateFacetOrds.data().rend(); ++i )
                {
                    if ( !tParentCellInfo->is_entity_connected_to_facet( *i, tFacetVertexParentOrds( iV ), tFacetVertexParentRanks( iV ) ) )
                    {
                        *i = MORIS_INDEX_MAX;
                    }
                }
            }

            tCandidateFacetOrds.data().erase( std::remove( tCandidateFacetOrds.data().begin(), tCandidateFacetOrds.data().end(), MORIS_INDEX_MAX ), tCandidateFacetOrds.data().end() );
        }

        if ( tCandidateFacetOrds.size() == 1 )
        {
            if ( tActiveConnectivity( (uint)aCutIntegrationMesh->get_facet_rank() ) == MORIS_INDEX_MAX )
            {
                ( *tEntityConnectedToParent( (uint)aCutIntegrationMesh->get_facet_rank() ) ) =
                    aBackgroundMesh->get_entity_connected_to_entity_loc_inds(
                        aParentCellForDeduction( iFacet )->get_index(),
                        EntityRank::ELEMENT,
                        aCutIntegrationMesh->get_facet_rank() );
            }
            // mark the entity parent as the deduction cell
            aIgFacetAncestry->mFacetParentEntityIndex( iFacet )                    = (*tEntityConnectedToParent( (uint)aCutIntegrationMesh->get_facet_rank() ))( tCandidateFacetOrds( 0 ) );
            aIgFacetAncestry->mFacetParentEntityRank( iFacet )                     = (moris_index)aCutIntegrationMesh->get_facet_rank();
            aIgFacetAncestry->mFacetParentEntityOrdinalWrtBackgroundCell( iFacet ) = tCandidateFacetOrds( 0 );
        }
        else
        {
            // mark the entity parent as the deduction cell
            aIgFacetAncestry->mFacetParentEntityIndex( iFacet )                    = aParentCellForDeduction( iFacet )->get_index();
            aIgFacetAncestry->mFacetParentEntityRank( iFacet )                     = 3;
            aIgFacetAncestry->mFacetParentEntityOrdinalWrtBackgroundCell( iFacet ) = 0;
        }
    }
}

// ----------------------------------------------------------------------------------

void
Integration_Mesh_Generator::compute_bg_facet_to_child_facet_connectivity(
    Cut_Integration_Mesh*                                                aCutIntegrationMesh,
    moris::mtk::Mesh*                                                    aBackgroundMesh,
    std::shared_ptr< Facet_Based_Connectivity >                          aIgCellGroupFacetConnectivity,
    std::shared_ptr< Facet_Based_Ancestry >                              aIgFacetAncestry,
    moris::Cell< std::shared_ptr< moris::Cell< moris::moris_index > > >& aBgFacetToIgFacet )
{
    // number of facets in the bg mesh
    moris_index tNumBgFacets = aBackgroundMesh->get_num_entities( aBackgroundMesh->get_facet_rank() );

    // resize the bg facet cell to the size of facets in the bg mesh
    aBgFacetToIgFacet.resize( tNumBgFacets );

    // std::cout<<"(moris_index)aCutIntegrationMesh->get_facet_rank() = "<<(moris_index)aCutIntegrationMesh->get_facet_rank()<<std::endl;

    // iterate through facets in the Facet_Based_Connectivity
    for ( moris::uint iFacet = 0; iFacet < aIgCellGroupFacetConnectivity->mFacetVertices.size(); iFacet++ )
    {
        // std::cout<<"iFacet= "<<iFacet<<std::endl;
        moris_index tParentRank = aIgFacetAncestry->mFacetParentEntityRank( iFacet );

        if ( tParentRank == (moris_index)aCutIntegrationMesh->get_facet_rank() )
        {
            moris_index tParentIndex = aIgFacetAncestry->mFacetParentEntityIndex( iFacet );

            if ( aBgFacetToIgFacet( tParentIndex ) == nullptr )
            {
                aBgFacetToIgFacet( tParentIndex ) = std::make_shared< moris::Cell< moris::moris_index > >();
            }

            // std::cout<<" tParentFacet = "<<tParentIndex<<" | Igfacet = "<<iFacet<<std::endl;
            aBgFacetToIgFacet( tParentIndex )->push_back( iFacet );
        }
    }
}

// ----------------------------------------------------------------------------------

void
Integration_Mesh_Generator::deduce_edge_ancestry(
    Cut_Integration_Mesh*                      aCutIntegrationMesh,
    moris::mtk::Mesh*                          aBackgroundMesh,
    std::shared_ptr< Edge_Based_Connectivity > aIgCellGroupEdgeConnectivity,
    moris::Cell< moris::mtk::Cell* > const&    aParentCellForDeduction,
    std::shared_ptr< Edge_Based_Ancestry >     aIgEdgeAncestry )
{
    Tracer tTracer( "XTK", "Integration_Mesh_Generator", "Edge Ancestry",mXTKModel->mVerboseLevel, 1  );
    // make sure we are starting clean
    aIgEdgeAncestry->mEdgeParentEntityIndex.clear();
    aIgEdgeAncestry->mEdgeParentEntityRank.clear();
    aIgEdgeAncestry->mEdgeParentEntityOrdinalWrtBackgroundCell.clear();

    // number of edges in the edge connectivity
    moris::uint tNumEdges = aIgCellGroupEdgeConnectivity->mEdgeVertices.size();

    MORIS_ERROR( aParentCellForDeduction.size() == tNumEdges, "One representative parent cell is needed for each edge, to ensure all edges parents can be deduced." );

    // allocate the data in the edge ancestry
    aIgEdgeAncestry->mEdgeParentEntityIndex.resize( tNumEdges );
    aIgEdgeAncestry->mEdgeParentEntityRank.resize( tNumEdges );
    aIgEdgeAncestry->mEdgeParentEntityOrdinalWrtBackgroundCell.resize( tNumEdges );

    enum EntityRank tFacetRank = aCutIntegrationMesh->get_facet_rank();

    // iterate through edges in the edge connectivity
    for ( moris::uint iEdge = 0; iEdge < tNumEdges; iEdge++ )
    {
        // vertices of the edge
        moris::Cell< moris::mtk::Vertex* > const& tEdgeVertices = aIgCellGroupEdgeConnectivity->mEdgeVertices( iEdge );

        // get the parent of these vertices from the mesh
        moris::Cell< moris::moris_index > tEdgeVertexParentInds( tEdgeVertices.size() );
        moris::Cell< moris::moris_index > tEdgeVertexParentRanks( tEdgeVertices.size() );

        MORIS_ERROR( tEdgeVertices.size() == 2, "Edge should have two vertices" );

        // loop over vertices
        for ( moris::uint iV = 0; iV < tEdgeVertices.size(); iV++ )
        {
            tEdgeVertexParentInds( iV )  = aCutIntegrationMesh->mIgVertexParentEntityIndex( tEdgeVertices( iV )->get_index() );
            tEdgeVertexParentRanks( iV ) = aCutIntegrationMesh->mIgVertexParentEntityRank( tEdgeVertices( iV )->get_index() );
        }

        // min and max parent ranks
        auto tMinIter = std::min_element( tEdgeVertexParentRanks.data().begin(), tEdgeVertexParentRanks.data().end() );
        auto tMaxIter = std::max_element( tEdgeVertexParentRanks.data().begin(), tEdgeVertexParentRanks.data().end() );

        // max index
        moris_index tMinIndex = std::distance( tEdgeVertexParentRanks.data().begin(), tMinIter );
        moris_index tMaxIndex = std::distance( tEdgeVertexParentRanks.data().begin(), tMaxIter );


        // std::cout<<"\nVertex 0 = "<<tEdgeVertices(0)->get_id()
        //       <<" | Vertex 1 = "<<tEdgeVertices(1)->get_id()
        //       <<" | Min Rank = "<<*tMinIter
        //       <<" | Max Rank = "<<*tMaxIter;

        // one of the vertices is on interior
        if ( *tMaxIter == 3 )
        {
            // mark the entity parent as the deduction cell
            aIgEdgeAncestry->mEdgeParentEntityIndex( iEdge )                    = aParentCellForDeduction( iEdge )->get_index();
            aIgEdgeAncestry->mEdgeParentEntityRank( iEdge )                     = 3;
            aIgEdgeAncestry->mEdgeParentEntityOrdinalWrtBackgroundCell( iEdge ) = 0;
        }
        else if ( *tMaxIter == (moris_index)tFacetRank && *tMinIter == (moris_index)tFacetRank )
        {
            tMinIndex = 0;
            tMaxIndex = 1;


            if ( tEdgeVertexParentInds( tMinIndex ) == tEdgeVertexParentInds( tMaxIndex ) )
            {

                // not as convenient to get the edge/facet ordinal
                Matrix< IndexMat > tEntitiesConnectedToBaseCell = aBackgroundMesh->get_entity_connected_to_entity_loc_inds(
                    aParentCellForDeduction( iEdge )->get_index(),
                    EntityRank::ELEMENT,
                    moris::get_entity_rank_from_index( *tMaxIter ) );

                moris_index tSecondEntityOrdinal = MORIS_INDEX_MAX;
                for ( moris::uint iEnt = 0; iEnt < tEntitiesConnectedToBaseCell.numel(); iEnt++ )
                {
                    if ( tEntitiesConnectedToBaseCell( iEnt ) == tEdgeVertexParentInds( tMaxIndex ) )
                    {
                        tSecondEntityOrdinal = iEnt;
                        break;
                    }
                }

                // mark the entity parent as the deduction cell
                aIgEdgeAncestry->mEdgeParentEntityIndex( iEdge )                    = tEdgeVertexParentInds( 0 );
                aIgEdgeAncestry->mEdgeParentEntityRank( iEdge )                     = (moris_index)tFacetRank;
                aIgEdgeAncestry->mEdgeParentEntityOrdinalWrtBackgroundCell( iEdge ) = tSecondEntityOrdinal;
            }

            else
            {
                // mark the entity parent as the deduction cell
                aIgEdgeAncestry->mEdgeParentEntityIndex( iEdge )                    = aParentCellForDeduction( iEdge )->get_index();
                aIgEdgeAncestry->mEdgeParentEntityRank( iEdge )                     = 3;
                aIgEdgeAncestry->mEdgeParentEntityOrdinalWrtBackgroundCell( iEdge ) = 0;
            }
        }
        else if ( *tMinIter > 0 )
        {
            if ( tMinIndex == tMaxIndex )
            {
                tMinIndex = 0;
                tMaxIndex = 1;
            }


            // edge to edge or edge to facet
            Matrix< IndexMat > tEntitiesConnectedToBaseCellMaxRank = aBackgroundMesh->get_entity_connected_to_entity_loc_inds(
                aParentCellForDeduction( iEdge )->get_index(),
                EntityRank::ELEMENT,
                moris::get_entity_rank_from_index( *tMaxIter ) );

            Matrix< IndexMat > tEntitiesConnectedToBaseCellMinRank = aBackgroundMesh->get_entity_connected_to_entity_loc_inds(
                aParentCellForDeduction( iEdge )->get_index(),
                EntityRank::ELEMENT,
                moris::get_entity_rank_from_index( *tMinIter ) );

            moris_index tMaxRankOrd = MORIS_INDEX_MAX;
            moris_index tMinRankOrd = MORIS_INDEX_MAX;

            for ( moris::uint iEnt = 0; iEnt < tEntitiesConnectedToBaseCellMaxRank.numel(); iEnt++ )
            {
                if ( tEntitiesConnectedToBaseCellMaxRank( iEnt ) == tEdgeVertexParentInds( tMaxIndex ) )
                {
                    tMaxRankOrd = iEnt;
                    break;
                }
            }


            for ( moris::uint iEnt = 0; iEnt < tEntitiesConnectedToBaseCellMinRank.numel(); iEnt++ )
            {
                if ( tEntitiesConnectedToBaseCellMinRank( iEnt ) == tEdgeVertexParentInds( tMinIndex ) )
                {
                    tMinRankOrd = iEnt;
                    break;
                }
            }


            MORIS_ASSERT( *tMinIter == 1, "Not an edge" );
            moris::Cell< moris::moris_index > tParentOrdinalAndRank = 
                    aParentCellForDeduction( iEdge )->get_cell_info()->get_edge_path_to_entity_rank_and_ordinal( tMinRankOrd, tMaxRankOrd, *tMaxIter );

            // std::cout<< " | tMinRankOrd = "<<tMinRankOrd<<" | tMaxRankOrd = "<<tMaxRankOrd<<" Parent Ord = "<<tParentOrdinalAndRank(0)<<" | ParentRank = "<<tParentOrdinalAndRank(1)<<std::endl;
            // need connectivity wrt the minimum path rank
            Matrix< IndexMat > tParentEntityConnectivity = aBackgroundMesh->get_entity_connected_to_entity_loc_inds(
                aParentCellForDeduction( iEdge )->get_index(),
                EntityRank::ELEMENT,
                moris::get_entity_rank_from_index( tParentOrdinalAndRank( 1 ) ) );

            // moris::print(tParentEntityConnectivity,"tParentEntityConnectivity");

            // mark the entity parent as the deduction cell
            aIgEdgeAncestry->mEdgeParentEntityIndex( iEdge )                    = tParentEntityConnectivity( tParentOrdinalAndRank( 0 ) );
            aIgEdgeAncestry->mEdgeParentEntityRank( iEdge )                     = tParentOrdinalAndRank( 1 );
            aIgEdgeAncestry->mEdgeParentEntityOrdinalWrtBackgroundCell( iEdge ) = tParentOrdinalAndRank( 0 );
        }
        else
        {
            if ( tMinIndex == tMaxIndex )
            {
                MORIS_ERROR( *tMinIter == 0, "Same min max but not a double parent node edge" );
                tMinIndex = 0;
                tMaxIndex = 1;
            }

            //minimum entity is always a vertex so i can ask the cell
            moris_index tVertexOrdinal = 
                    aParentCellForDeduction( iEdge )->get_vertex_ordinal_wrt_cell( tEdgeVertices( tMinIndex )->get_index() );

            // not as convenient to get the edge/facet ordinal
            Matrix< IndexMat > tEntitiesConnectedToBaseCell = aBackgroundMesh->get_entity_connected_to_entity_loc_inds(
                aParentCellForDeduction( iEdge )->get_index(),
                EntityRank::ELEMENT,
                moris::get_entity_rank_from_index( *tMaxIter ) );

            moris_index tSecondEntityOrdinal = MORIS_INDEX_MAX;
            for ( moris::uint iEnt = 0; iEnt < tEntitiesConnectedToBaseCell.numel(); iEnt++ )
            {
                if ( tEntitiesConnectedToBaseCell( iEnt ) == tEdgeVertexParentInds( tMaxIndex ) )
                {
                    tSecondEntityOrdinal = iEnt;
                    break;
                }
            }

            moris::Cell< moris::moris_index > tParentOrdinalAndRank = 
                    aParentCellForDeduction( iEdge )->get_cell_info()->get_vertex_path_to_entity_rank_and_ordinal( tVertexOrdinal, tSecondEntityOrdinal, *tMaxIter );

            // need connectivity wrt the minimum path rank
            tEntitiesConnectedToBaseCell = aBackgroundMesh->get_entity_connected_to_entity_loc_inds(
                aParentCellForDeduction( iEdge )->get_index(),
                EntityRank::ELEMENT,
                moris::get_entity_rank_from_index( tParentOrdinalAndRank( 1 ) ) );

            if ( tSecondEntityOrdinal == MORIS_INDEX_MAX )
            {
                std::cout << "tSecondEntityOrdinal = " << tSecondEntityOrdinal << std::endl;
                moris::print( tEntitiesConnectedToBaseCell, "tEntitiesConnectedToBaseCell" );

                moris::print( tEdgeVertexParentInds, "tEdgeVertexParentInds" );
                moris::print( tEdgeVertexParentRanks, "tEdgeVertexParentRanks" );
                moris::print( tEdgeVertices( 0 )->get_coords(), "Edge Vertex 0 coords" );
                moris::print( tEdgeVertices( 1 )->get_coords(), "Edge Vertex 1 coords" );
            }

            MORIS_ASSERT( tSecondEntityOrdinal != MORIS_INDEX_MAX, "Not found" );

            // mark the entity parent as the deduction cell
            aIgEdgeAncestry->mEdgeParentEntityIndex( iEdge )                    = tEntitiesConnectedToBaseCell( tParentOrdinalAndRank( 0 ) );
            aIgEdgeAncestry->mEdgeParentEntityRank( iEdge )                     = tParentOrdinalAndRank( 1 );
            aIgEdgeAncestry->mEdgeParentEntityOrdinalWrtBackgroundCell( iEdge ) = tSecondEntityOrdinal;
        }
    }
}

// ----------------------------------------------------------------------------------

void
Integration_Mesh_Generator::commit_new_ig_vertices_to_cut_mesh(
    Integration_Mesh_Generation_Data* aMeshGenerationData,
    Decomposition_Data*               aDecompositionData,
    Cut_Integration_Mesh*             aCutIntegrationMesh,
    moris::mtk::Mesh*                 aBackgroundMesh,
    Decomposition_Algorithm*          aDecompAlg )
{
    Tracer tTracer( "XTK", "Integration_Mesh_Generator", "Commit Ig Vertices to Cut Mesh",mXTKModel->mVerboseLevel, 1  );

    if ( mXTKModel->mDiagnostics )
    {
        std::string tRequestDiagFile = 
                mXTKModel->get_diagnostic_file_name( std::string( "Requests_rid_" + std::to_string( aDecompositionData->tDecompId ) ) );
        aDecompositionData->print( *aBackgroundMesh, tRequestDiagFile );
    }

    // current index
    moris_index tControlledVertexIndex = aCutIntegrationMesh->mControlledIgVerts.size();

    // allocate new vertices
    moris::uint tNumNewIgVertices = aDecompositionData->tNewNodeIndex.size();
    aCutIntegrationMesh->mControlledIgVerts.resize( aCutIntegrationMesh->mControlledIgVerts.size() + tNumNewIgVertices );
    aCutIntegrationMesh->mIntegrationVertices.resize( aCutIntegrationMesh->mIntegrationVertices.size() + tNumNewIgVertices );
    aCutIntegrationMesh->mVertexCoordinates.resize( aCutIntegrationMesh->mIntegrationVertices.size(), nullptr );
    aCutIntegrationMesh->mIgVertexParentEntityRank.resize( aCutIntegrationMesh->mIntegrationVertices.size(), MORIS_INDEX_MAX );
    aCutIntegrationMesh->mIgVertexParentEntityIndex.resize( aCutIntegrationMesh->mIntegrationVertices.size(), MORIS_INDEX_MAX );

    // iterate and create new vertices
    for ( moris::uint iV = 0; iV < aDecompositionData->tNewNodeId.size(); iV++ )
    {
        // construct coordinate matrix
        aCutIntegrationMesh->mVertexCoordinates( aDecompositionData->tNewNodeIndex( iV ) ) = 
                std::make_shared< Matrix< DDRMat > >( aDecompositionData->tNewNodeCoordinate( iV ) );

        // create a controlled vertex (meaning I need to manage memory of it)
        aCutIntegrationMesh->mControlledIgVerts( tControlledVertexIndex ) = std::make_shared< moris::mtk::Vertex_XTK >(
            aDecompositionData->tNewNodeId( iV ),
            aDecompositionData->tNewNodeIndex( iV ),
            aDecompositionData->tNewNodeOwner( iV ),
            aCutIntegrationMesh->mVertexCoordinates( aDecompositionData->tNewNodeIndex( iV ) ) );

        // add vertex coordinates to the mesh data
        aCutIntegrationMesh->mIntegrationVertices( aDecompositionData->tNewNodeIndex( iV ) ) = 
                aCutIntegrationMesh->mControlledIgVerts( tControlledVertexIndex ).get();
        tControlledVertexIndex++;

        // add to the map
        MORIS_ASSERT( aCutIntegrationMesh->mIntegrationVertexIdToIndexMap.find( 
                aDecompositionData->tNewNodeId( iV ) ) == aCutIntegrationMesh->mIntegrationVertexIdToIndexMap.end(), "Id already in the map" );
        aCutIntegrationMesh->mIntegrationVertexIdToIndexMap[aDecompositionData->tNewNodeId( iV )] = aDecompositionData->tNewNodeIndex( iV );

        // add the ancestry information to the mesh
        aCutIntegrationMesh->mIgVertexParentEntityRank( aDecompositionData->tNewNodeIndex( iV ) )  = 
                (moris_index)aDecompositionData->tNewNodeParentRank( iV );
        aCutIntegrationMesh->mIgVertexParentEntityIndex( aDecompositionData->tNewNodeIndex( iV ) ) = 
                (moris_index)aDecompositionData->tNewNodeParentIndex( iV );
    }


    // iterate through child meshes and commit the vertices to their respective vertex groups
    for ( auto& iCell : aMeshGenerationData->mAllIntersectedBgCellInds )
    {
        MORIS_ERROR( aDecompositionData->tCMNewNodeLoc.size() == (moris::uint)aCutIntegrationMesh->mChildMeshes.size(), 
                "Mismatch in child mesh sizes. All child meshes need to be present in the decomposition data" );

        // add the vertices to child mesh groups
        moris_index tNumNewVertices = (moris_index)aDecompositionData->tCMNewNodeLoc( iCell ).size();

        // resize the vertices in the group
        aCutIntegrationMesh->mIntegrationVertexGroups( iCell )->reserve( tNumNewVertices + aCutIntegrationMesh->mIntegrationVertexGroups( iCell )->size() );

        for ( moris::moris_index iCMVerts = 0; iCMVerts < tNumNewVertices; iCMVerts++ )
        {
            moris_index      tNewNodeLocInDecomp      = aDecompositionData->tCMNewNodeLoc( iCell )( iCMVerts );
            moris_index      tNewNodeIndex            = aDecompositionData->tNewNodeIndex( tNewNodeLocInDecomp );
            Matrix< DDRMat > tParamCoordWrtParentCell = aDecompositionData->tCMNewNodeParamCoord( iCell )( iCMVerts );

            aCutIntegrationMesh->mIntegrationVertexGroups( iCell )->add_vertex( aCutIntegrationMesh->mIntegrationVertices( tNewNodeIndex ), std::make_shared< Matrix< DDRMat > >( tParamCoordWrtParentCell ) );
        }
    }

    // construct a relationship for the geometry engine to have geomtric information for this vertex
    this->link_new_vertices_to_geometry_engine( aDecompositionData, aDecompAlg );
}

// ----------------------------------------------------------------------------------

void
Integration_Mesh_Generator::link_new_vertices_to_geometry_engine(
    Decomposition_Data*      aDecompositionData,
    Decomposition_Algorithm* aDecompAlg )

{
    if ( aDecompAlg->has_geometric_independent_vertices() )
    {
        // pass the data in decomposition data to the geometry engine so it can keep track of these newly constructed vertices
        mGeometryEngine->create_new_child_nodes(
            &aDecompositionData->tNewNodeIndex,
            &aDecompositionData->mNewNodeParentCells,
            &aDecompositionData->mNewVertexLocalCoordWRTParentCell,
            &aDecompositionData->tNewNodeCoordinate );
    }
}

// ----------------------------------------------------------------------------------

void
Integration_Mesh_Generator::select_background_cell_for_edge(
    std::shared_ptr< Edge_Based_Connectivity > aEdgeBasedConnectivity,
    Cut_Integration_Mesh*                      aCutIntegrationMesh,
    moris::Cell< moris::mtk::Cell* >&          aBackgroundCellForEdge )
{
    Tracer tTracer( "XTK", "Decomposition_Algorithm", "Select BG Cell for Edge" ,mXTKModel->mVerboseLevel, 1  );
    // number of edges
    moris::uint tNumEdges = aEdgeBasedConnectivity->mEdgeVertices.size();

    aBackgroundCellForEdge.resize( tNumEdges );

    // iterate through edges
    for ( moris::uint iEdge = 0; iEdge < tNumEdges; iEdge++ )
    {
        // get the first cell attached to the edge
        MORIS_ERROR( aEdgeBasedConnectivity->mEdgeToCell( iEdge ).size() > 0, "Edge not connected to any cells..." );

        // integration cell just grab the first
        moris::mtk::Cell* tCell = aEdgeBasedConnectivity->mEdgeToCell( iEdge )( 0 );

        // cell group membership
        moris_index tCellGroupMembershipIndex = aCutIntegrationMesh->get_ig_cell_group_memberships( tCell->get_index() )( 0 );

        aBackgroundCellForEdge( iEdge ) = aCutIntegrationMesh->get_ig_cell_group_parent_cell( tCellGroupMembershipIndex );
    }
}

// ----------------------------------------------------------------------------------

void
Integration_Mesh_Generator::select_background_cell_for_facet(
    std::shared_ptr< Facet_Based_Connectivity > aFacetBasedConnectivity,
    Cut_Integration_Mesh*                       aCutIntegrationMesh,
    moris::Cell< moris::mtk::Cell* >&           aBackgroundCellForEdge )
{
    Tracer tTracer( "XTK", "Decomposition_Algorithm", "Select BG Cell for Facet" ,mXTKModel->mVerboseLevel, 1  );
    // number of edges
    moris::uint tNumFacets = aFacetBasedConnectivity->mFacetVertices.size();

    aBackgroundCellForEdge.resize( tNumFacets );

    // iterate through edges
    for ( moris::uint iEdge = 0; iEdge < tNumFacets; iEdge++ )
    {
        // get the first cell attached to the edge
        MORIS_ERROR( aFacetBasedConnectivity->mFacetToCell( iEdge ).size() > 0, "Facet not connected to any cells..." );

        // integration cell just grab the first
        moris::mtk::Cell* tCell = aFacetBasedConnectivity->mFacetToCell( iEdge )( 0 );

        // cell group membership
        moris_index tCellGroupMembershipIndex = aCutIntegrationMesh->get_ig_cell_group_memberships( tCell->get_index() )( 0 );

        aBackgroundCellForEdge( iEdge ) = aCutIntegrationMesh->get_ig_cell_group_parent_cell( tCellGroupMembershipIndex );
    }
}

// ----------------------------------------------------------------------------------

void
Integration_Mesh_Generator::collect_vertex_groups_for_background_cells(
    Integration_Mesh_Generation_Data*                  aMeshGenerationData,
    Cut_Integration_Mesh*                              aCutIntegrationMesh,
    moris::Cell< moris::mtk::Cell* >*                  aBackgroundCells,
    moris::Cell< std::shared_ptr< IG_Vertex_Group > >* aVertexGroups )
{
    Tracer tTracer( "XTK", "Integration_Mesh_Generator", "Collect Vertex Groups for BG Cell" ,mXTKModel->mVerboseLevel, 1  );
    aVertexGroups->resize( aBackgroundCells->size() );

    // iterate through background cells
    for ( moris::uint i = 0; i < aBackgroundCells->size(); i++ )
    {
        ( *aVertexGroups )( i ) = aCutIntegrationMesh->get_vertex_group( ( *aBackgroundCells )( i )->get_index() );
    }
}

// ----------------------------------------------------------------------------------

bool
Integration_Mesh_Generator::allocate_child_meshes(
    Integration_Mesh_Generation_Data& aMeshGenerationData,
    Cut_Integration_Mesh*             aCutIntegrationMesh,
    moris::mtk::Mesh*                 aBackgroundMesh )
{
    Tracer tTracer( "XTK", "Integration_Mesh_Generator", "Allocate child meshes",mXTKModel->mVerboseLevel, 1  );

    // allocate these data structures one per background cell
    aCutIntegrationMesh->mChildMeshes.resize( aBackgroundMesh->get_num_elems() );
    aCutIntegrationMesh->mIntegrationCellGroups.resize( aBackgroundMesh->get_num_elems(), nullptr );
    aCutIntegrationMesh->mIntegrationVertexGroups.resize( aBackgroundMesh->get_num_elems(), nullptr );
    aCutIntegrationMesh->mIntegrationCellGroupsParentCell.resize( aBackgroundMesh->get_num_elems(), nullptr );

    // create the child meshes
    for ( moris::uint iCell = 0; iCell < aBackgroundMesh->get_num_elems(); iCell++ )
    {
        moris_index       tCMIndex                                        = (moris_index)iCell;
        moris::mtk::Cell* tParentCell                                     = &aBackgroundMesh->get_mtk_cell( iCell );
        aCutIntegrationMesh->mChildMeshes( tCMIndex )                     = std::make_shared< Child_Mesh_Experimental >();
        aCutIntegrationMesh->mIntegrationCellGroups( tCMIndex )           = std::make_shared< IG_Cell_Group >( 0 );
        aCutIntegrationMesh->mChildMeshes( tCMIndex )->mIgCells           = aCutIntegrationMesh->mIntegrationCellGroups( tCMIndex );
        aCutIntegrationMesh->mIntegrationCellGroupsParentCell( tCMIndex ) = tParentCell;
        aCutIntegrationMesh->mChildMeshes( tCMIndex )->mParentCell        = tParentCell;
        aCutIntegrationMesh->mChildMeshes( tCMIndex )->mChildMeshIndex    = tCMIndex;

        aCutIntegrationMesh->mParentCellCellGroupIndex( tParentCell->get_index() ) = tCMIndex;

        // fixme: std::cout << "Integration_Mesh_Generator::allocate_child_meshes() - WARNING: GENERAlIZE NEEDED FOR MULTIPLE TOPOS" << std::endl;

        // initialize cell topology
        moris_index tNumGeometricVertices;
        
        // fixme: is there a more elegant way to decide on the tNumGeometricVertices?
        // get number of spatial dimensions and decide on cell topology of integration elements
        if ( this->get_spatial_dim() == 2 )
        {
            tNumGeometricVertices = 4;
        }
        else if ( this->get_spatial_dim() == 3 ) 
        {
            tNumGeometricVertices = 8;
        }
        else
        {
            MORIS_ERROR( false, "Integration_Mesh_Generator::construct_bulk_phase_blocks() - spatial dimension not 2 or 3" );
        }

        moris::Cell< moris::mtk::Vertex* > tParentCellVerts      = tParentCell->get_vertex_pointers();

        Matrix< DDRMat > tParamCoords;
        tParentCell->get_cell_info()->get_loc_coords_of_cell( tParamCoords );

        aCutIntegrationMesh->mIntegrationVertexGroups( tCMIndex ) = std::make_shared< IG_Vertex_Group >( tNumGeometricVertices );
        aCutIntegrationMesh->mChildMeshes( tCMIndex )->mIgVerts   = aCutIntegrationMesh->mIntegrationVertexGroups( tCMIndex );
        //FIXME: GET GEOMETRIC VERTICES FROM MTK CELL HARDCODED TO HEXFAMILY
        for ( moris::moris_index i = 0; i < tNumGeometricVertices; i++ )
        {
            aCutIntegrationMesh->mIntegrationVertexGroups( tCMIndex )->add_vertex( tParentCellVerts( i ), std::make_shared< Matrix< DDRMat > >( tParamCoords.get_row( i ) ) );
        }
    }

    return true;
}

// ----------------------------------------------------------------------------------

void
Integration_Mesh_Generator::assign_node_requests_identifiers(
    Decomposition_Data&   aDecompData,
    Cut_Integration_Mesh* aCutIntegrationMesh,
    moris::mtk::Mesh*     aBackgroundMesh,
    moris::moris_index    aMPITag )
{
    Tracer      tTracer( "XTK", "Decomposition_Algorithm", "Assign node ids" ,mXTKModel->mVerboseLevel, 1  );
    moris_index tNodeIndex = aCutIntegrationMesh->get_first_available_index( EntityRank::NODE );

    for ( moris::uint i = 0; i < aDecompData.tNewNodeIndex.size(); i++ )
    {
        // set the new node index
        aDecompData.tNewNodeIndex( i ) = tNodeIndex;
        tNodeIndex++;
    }

    barrier();
    // asserts
    MORIS_ERROR( aDecompData.tNewNodeId.size() == aDecompData.tNewNodeIndex.size(),
        "Dimension mismatch in assign_node_requests_identifiers" );
    MORIS_ERROR( aDecompData.tNewNodeId.size() == aDecompData.tNewNodeParentRank.size(),
        "Dimension mismatch in assign_node_requests_identifiers" );
    MORIS_ERROR( aDecompData.tNewNodeId.size() == aDecompData.tNewNodeParentIndex.size(),
        "Dimension mismatch in assign_node_requests_identifiers" );

    // owned requests and shared requests sorted by owning proc
    Cell< uint >                             tOwnedRequest;
    Cell< Cell< uint > >                     tNotOwnedRequests;
    Cell< uint >                             tProcRanks;
    std::unordered_map< moris_id, moris_id > tProcRankToDataIndex;
    this->sort_new_node_requests_by_owned_and_not_owned(
        aDecompData,
        aCutIntegrationMesh,
        aBackgroundMesh,
        tOwnedRequest,
        tNotOwnedRequests,
        tProcRanks,
        tProcRankToDataIndex );

    // allocate ids for nodes I own
    moris::moris_id tNodeId = aCutIntegrationMesh->allocate_entity_ids( aDecompData.tNewNodeId.size(), EntityRank::NODE );

    // Assign owned request identifiers
    this->assign_owned_request_id( aDecompData, tOwnedRequest, tNodeId );

    // prepare node information request data
    Cell< Matrix< IndexMat > > tOutwardRequests;
    this->setup_outward_requests( aDecompData, aBackgroundMesh, tNotOwnedRequests, tProcRanks, tProcRankToDataIndex, tOutwardRequests );

    // send requests to owning processor
    mXTKModel->send_outward_requests( aMPITag, tProcRanks, tOutwardRequests );

    // hold on to make sure everyone has sent all their information
    barrier();

    // receive the requests
    Cell< Matrix< IndexMat > > tReceivedRequests;
    Cell< uint >               tProcsReceivedFrom;
    mXTKModel->inward_receive_requests( aMPITag, 3, tReceivedRequests, tProcsReceivedFrom );

    // Prepare request answers
    Cell< Matrix< IndexMat > > tRequestAnwers;
    this->prepare_request_answers( aDecompData, aBackgroundMesh, tReceivedRequests, tRequestAnwers );

    // send the answers back
    mXTKModel->return_request_answers( aMPITag + 1, tRequestAnwers, tProcsReceivedFrom );

    barrier();

    // receive the answers
    Cell< Matrix< IndexMat > > tReceivedRequestsAnswers;
    mXTKModel->inward_receive_request_answers( aMPITag + 1, 1, tProcRanks, tReceivedRequestsAnswers );

    // handle received information
    this->handle_received_request_answers( aDecompData, aBackgroundMesh, tOutwardRequests, tReceivedRequestsAnswers, tNodeId );

    MORIS_ERROR( mXTKModel->verify_successful_node_assignment( aDecompData ),
        "Unsuccesssful node assignment detected." );

    barrier();
}

// ----------------------------------------------------------------------------------

void
Integration_Mesh_Generator::sort_new_node_requests_by_owned_and_not_owned(
    Decomposition_Data&                       tDecompData,
    Cut_Integration_Mesh*                     aCutIntegrationMesh,
    moris::mtk::Mesh*                         aBackgroundMesh,
    Cell< uint >&                             aOwnedRequests,
    Cell< Cell< uint > >&                     aNotOwnedRequests,
    Cell< uint >&                             aProcRanks,
    std::unordered_map< moris_id, moris_id >& aProcRankToIndexInData )
{
    // access the communication
    Matrix< IdMat > tCommTable = aCutIntegrationMesh->get_communication_table();

    // number of new nodes
    moris::uint tNumNewNodes = tDecompData.tNewNodeParentIndex.size();

    // reserve memory
    aOwnedRequests.reserve( tNumNewNodes );
    aNotOwnedRequests.reserve( tCommTable.numel() );

    // Par rank
    moris::moris_index tParRank = par_rank();

    // resize proc ranks and setup map to comm table
    aProcRanks.resize( tCommTable.numel() );
    for ( moris::uint i = 0; i < tCommTable.numel(); i++ )
    {
        aProcRankToIndexInData[tCommTable( i )] = i;
        aProcRanks( i )                         = ( tCommTable( i ) );

        // initialize and reserve memory for not-owned nodes
        aNotOwnedRequests.push_back( Cell< uint >( 0 ) );
        aNotOwnedRequests.back().reserve( 2*tNumNewNodes / tCommTable.numel() );
    }

    // iterate through each node request and figure out the owner
    for ( moris::uint i = 0; i < tNumNewNodes; i++ )
    {
        // Parent Rank
        enum EntityRank    tParentRank  = tDecompData.tNewNodeParentRank( i );
        moris::moris_index tParentIndex = tDecompData.tNewNodeParentIndex( i );

        // get the owner processor
        moris::moris_index tOwnerProc = aBackgroundMesh->get_entity_owner( tParentIndex, tParentRank );

        // If i own the request keep track of the index
        if ( tOwnerProc == tParRank )
        {
            aOwnedRequests.push_back( i );
        }
        else
        {
            moris_index tIndex = aProcRankToIndexInData[tOwnerProc];

            aNotOwnedRequests( tIndex ).push_back( i );
        }
    }

    // trim cells
    shrink_to_fit_all( aOwnedRequests );
    shrink_to_fit_all( aNotOwnedRequests );
}

// ----------------------------------------------------------------------------------

void
Integration_Mesh_Generator::assign_owned_request_id(
    Decomposition_Data& aDecompData,
    Cell< uint > const& aOwnedRequest,
    moris::moris_id&    aNodeId )
{
    for ( moris::uint i = 0; i < aOwnedRequest.size(); i++ )
    {
        moris_index tRequestIndex = aOwnedRequest( i );

        // set the new node id
        aDecompData.tNewNodeId( tRequestIndex ) = aNodeId;
        aNodeId++;

        // increment number of new nodes with set ids (for assertion purposes)
        aDecompData.mNumNewNodesWithIds++;
    }
}

// ----------------------------------------------------------------------------------

void
Integration_Mesh_Generator::setup_outward_requests(
    Decomposition_Data const&                 aDecompData,
    moris::mtk::Mesh*                         aBackgroundMesh,
    Cell< Cell< uint > > const&               aNotOwnedRequests,
    Cell< uint > const&                       aProcRanks,
    std::unordered_map< moris_id, moris_id >& aProcRankToIndexInData,
    Cell< Matrix< IndexMat > >&               aOutwardRequests )
{
    // size data
    aOutwardRequests.resize( aProcRanks.size() );

    // iterate through the processors we need information from and package the matrix
    for ( moris::uint i = 0; i < aProcRanks.size(); i++ )
    {
        uint tProcRank = aProcRanks( i );

        MORIS_ASSERT( aProcRankToIndexInData.find( tProcRank ) != aProcRankToIndexInData.end(), "Proc rank not in map" );
        uint tIndexInData = aProcRankToIndexInData[tProcRank];

        uint tNumRequests = aNotOwnedRequests( tIndexInData ).size();

        // size the sending matrix
        // column - request
        //   r0 - parent entity id
        //   r1 - parent entity rank
        //   r2 - Secondary id
        if ( tNumRequests > 0 )
        {
            aOutwardRequests( i ) = moris::Matrix< IndexMat >( 3, tNumRequests );
        }

        else
        {
            aOutwardRequests( i ) = moris::Matrix< IndexMat >( 3, 1, MORIS_INDEX_MAX );
        }

        // populate matrix to send;
        for ( moris::uint j = 0; j < tNumRequests; j++ )
        {
            moris_index     tRequestIndex = aNotOwnedRequests( tIndexInData )( j );
            moris_index     tParentIndex  = aDecompData.tNewNodeParentIndex( tRequestIndex );
            moris_index     tSecondaryId  = aDecompData.tSecondaryIdentifiers( tRequestIndex );
            enum EntityRank tParentRank   = aDecompData.tNewNodeParentRank( tRequestIndex );

            // swap out for hmr if needed (hmr calls edges in 2d faces)
            if(aBackgroundMesh->get_mesh_type() == MeshType::HMR )
            {
                if(aBackgroundMesh->get_spatial_dim() == 2)
                {
                    if(tParentRank == EntityRank::EDGE)
                    {
                        tParentRank =  EntityRank::FACE;
                    }
                }
            }

            aOutwardRequests( i )( 1, j ) = (moris_index)tParentRank;
            aOutwardRequests( i )( 0, j ) = aBackgroundMesh->get_glb_entity_id_from_entity_loc_index( tParentIndex, tParentRank );
            aOutwardRequests( i )( 2, j ) = tSecondaryId;
        }
    }
}

// ----------------------------------------------------------------------------------

void
Integration_Mesh_Generator::prepare_request_answers(
    Decomposition_Data&               aDecompData,
    moris::mtk::Mesh*                 aBackgroundMesh,
    Cell< Matrix< IndexMat > > const& aReceiveData,
    Cell< Matrix< IndexMat > >&       aRequestAnswers )
{
    // allocate answer size
    aRequestAnswers.resize( aReceiveData.size() );

    // iterate through received data
    for ( moris::uint i = 0; i < aReceiveData.size(); i++ )
    {
        uint tNumReceivedReqs = aReceiveData( i ).n_cols();

        aRequestAnswers( i ).resize( 1, tNumReceivedReqs );

        aRequestAnswers( i )( 0 ) = MORIS_INDEX_MAX;

        // avoid the dummy message
        if ( aReceiveData( i )( 0, 0 ) != MORIS_INDEX_MAX )
        {
            // iterate through received requests
            for ( moris::uint j = 0; j < tNumReceivedReqs; j++ )
            {
                moris_id        tParentId      = aReceiveData( i )( 0, j );
                enum EntityRank tParentRank    = (enum EntityRank)aReceiveData( i )( 1, j );
                moris_id        tSecondaryId   = aReceiveData( i )( 2, j );
                moris_index     tParentInd     = aBackgroundMesh->get_loc_entity_ind_from_entity_glb_id( tParentId, tParentRank );
                bool            tRequestExists = false;
                moris_index     tRequestIndex  = MORIS_INDEX_MAX;

                // swap out for hmr if needed (hmr calls edges in 2d faces)
                if(aBackgroundMesh->get_mesh_type() == MeshType::HMR )
                {
                    if(aBackgroundMesh->get_spatial_dim() == 2)
                    {
                        if(tParentRank == EntityRank::FACE)
                        {
                            tParentRank =  EntityRank::EDGE;
                        }
                    }
                }

                if ( aDecompData.mHasSecondaryIdentifier )
                {
                    tRequestExists = aDecompData.request_exists(
                        tParentInd,
                        tSecondaryId,
                        (EntityRank)tParentRank,
                        tRequestIndex );
                }
                else
                {
                    tRequestExists = aDecompData.request_exists(
                        tParentInd,
                        (EntityRank)tParentRank,
                        tRequestIndex );
                }

                if ( tRequestExists )
                {
                    moris_id tNodeId = aDecompData.tNewNodeId( tRequestIndex );

                    aRequestAnswers( i )( j ) = tNodeId;

                    if ( tNodeId == MORIS_ID_MAX )
                    {
                        std::cout << "tParentId = " << tParentId << " | Rank " << (uint)tParentRank << std::endl;
                        //                    MORIS_ERROR(0,"Max node");
                    }
                }
                else
                {
                    aRequestAnswers( i )( j ) = MORIS_ID_MAX;
                }
            }
        }
    }
}

// ----------------------------------------------------------------------------------

void
Integration_Mesh_Generator::handle_received_request_answers(
    Decomposition_Data&               aDecompData,
    moris::mtk::Mesh*                 aBackgroundMesh,
    Cell< Matrix< IndexMat > > const& aRequests,
    Cell< Matrix< IndexMat > > const& aRequestAnswers,
    moris::moris_id&                  aNodeId )
{
    Cell< moris_index > tUnhandledRequestIndices;

    // iterate through received data
    for ( moris::uint i = 0; i < aRequests.size(); i++ )
    {
        uint tNumReceivedReqs = aRequests( i ).n_cols();

        // avoid the dummy message
        if ( aRequests( i )( 0, 0 ) != MORIS_INDEX_MAX )
        {
            // iterate through received requests
            for ( moris::uint j = 0; j < tNumReceivedReqs; j++ )
            {
                moris_id        tParentId      = aRequests( i )( 0, j );
                enum EntityRank tParentRank    = (enum EntityRank)aRequests( i )( 1, j );
                moris_id        tSecondaryId   = aRequests( i )( 2, j );
                moris_index     tParentInd     = aBackgroundMesh->get_loc_entity_ind_from_entity_glb_id( tParentId, tParentRank );
                bool            tRequestExists = false;
                moris_index     tRequestIndex  = MORIS_INDEX_MAX;

                
                // swap out for hmr if needed (hmr calls edges in 2d faces)
                if(aBackgroundMesh->get_mesh_type() == MeshType::HMR )
                {
                    if(aBackgroundMesh->get_spatial_dim() == 2)
                    {
                        if(tParentRank == EntityRank::FACE)
                        {
                            tParentRank =  EntityRank::EDGE;
                        }
                    }
                }


                if ( aDecompData.mHasSecondaryIdentifier )
                {
                    tRequestExists = aDecompData.request_exists( tParentInd, tSecondaryId, (EntityRank)tParentRank, tRequestIndex );
                }
                else
                {
                    tRequestExists = aDecompData.request_exists( tParentInd, (EntityRank)tParentRank, tRequestIndex );
                }

                if ( tRequestExists && aRequestAnswers( i )( j ) )
                {
                    moris_id tNodeId = aRequestAnswers( i )( j );

                    // meaning the owning processor expected this and gave an answer
                    if ( tNodeId < MORIS_ID_MAX && aDecompData.tNewNodeId( tRequestIndex ) == MORIS_INDEX_MAX )
                    {
                        // set the new node id
                        aDecompData.tNewNodeId( tRequestIndex ) = tNodeId;

                        aDecompData.mNumNewNodesWithIds++;
                    }
                    // The owner did not expect and did not return an answer
                    else
                    {
                        // keep track of unhandled
                        tUnhandledRequestIndices.push_back( tRequestIndex );
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
                    MORIS_ASSERT( 0, "Request does not exist." );
                }
            }
        }
    }
}

// ----------------------------------------------------------------------------------

moris::uint
Integration_Mesh_Generator::verbosity_level()
{
    return mXTKModel->mVerboseLevel;
}

// ----------------------------------------------------------------------------------

void
Integration_Mesh_Generator::remove_subphases_from_cut_mesh(moris::Cell<moris_index> const & aSubphasesToRemove)
{

}

// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------

}// namespace xtk
