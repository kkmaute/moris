/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_XTK_Integration_Mesh_Generator.cpp
 *
 */

#include <unordered_map>
#include "cl_XTK_Integration_Mesh_Generator.hpp"
#include "cl_XTK_Decomposition_Algorithm_Factory.hpp"
#include "cl_XTK_Decomposition_Algorithm.hpp"
#include "fn_determine_cell_topology.hpp"
#include "fn_mesh_flood_fill.hpp"
#include "fn_XTK_find_most_frequent_int_in_cell.hpp"
#include "fn_XTK_Multiset_Operations.hpp"
#include "fn_stringify_matrix.hpp"

using namespace moris;
namespace xtk
{
    // ----------------------------------------------------------------------------------

    Integration_Mesh_Generator::Integration_Mesh_Generator()
    {
    }

    // ----------------------------------------------------------------------------------

    Integration_Mesh_Generator::Integration_Mesh_Generator(
            xtk::Model*                      aXTKModelPtr,
            Cell< enum Subdivision_Method >  aMethods )
            : mXTKModel( aXTKModelPtr )
            , mGeometryEngine( mXTKModel->get_geom_engine() )
            , mSubdivisionMethods( aMethods )
    {
        if ( mXTKModel->mParameterList.get< bool >( "has_parameter_list" ) )
        {
            mOutputCutIgMesh = mXTKModel->mParameterList.get< bool >( "output_cut_ig_mesh" );
        }
    }
    // ----------------------------------------------------------------------------------

    Integration_Mesh_Generator::~Integration_Mesh_Generator()
    {
    }

    // ----------------------------------------------------------------------------------

    std::shared_ptr< Cut_Integration_Mesh >
    Integration_Mesh_Generator::perform()
    {
        // log/trace this operation
        Tracer tTracer( "XTK", "Integration_Mesh_Generator", "perform", mXTKModel->mVerboseLevel, 0 );

        // data structure containing which BG cells are intersected to later trigger decomposition of those elements
        Integration_Mesh_Generation_Data tGenerationData;

        // pointer to the background mesh
        moris::mtk::Mesh* tBackgroundMesh = &mXTKModel->get_background_mesh();

        // initialize the cut integration mesh object
        std::shared_ptr< Cut_Integration_Mesh > tCutIntegrationMesh = std::make_shared< Cut_Integration_Mesh >( tBackgroundMesh, mXTKModel );

        // hand pointer to cut integration mesh to XTK-Model (which will own this shared pointer after the IMG is done)
        mXTKModel->set_cut_ig_mesh( tCutIntegrationMesh );

        // Allocate a child mesh for each background cell
        this->allocate_child_meshes( tGenerationData, tCutIntegrationMesh.get(), tBackgroundMesh );

        // figure out which background cells are intersected and by which geometry they are intersected
        this->determine_intersected_background_cells( tGenerationData, tCutIntegrationMesh.get(), tBackgroundMesh );

        // verify levels of intersected background cells
        this->check_intersected_background_cell_levels( tGenerationData, tCutIntegrationMesh.get(), tBackgroundMesh );

        // make mSameLevelChildMeshes flag global consistent
        if ( !tCutIntegrationMesh->mSameLevelChildMeshes )
        {
            return tCutIntegrationMesh;
        }

        // iterate through the subdivision methods
        for ( moris::uint iSubMethod = 0; iSubMethod < mSubdivisionMethods.size(); iSubMethod++ )
        {
            // create the subdivision routine with a factory
            std::shared_ptr< Decomposition_Algorithm > tDecompositionAlgorithm =
                    create_decomposition_algorithm( mSubdivisionMethods( iSubMethod ), mXTKModel->get_parameter_list() );

            // perform the decomposition
            Decomposition_Data tDecompositionData;

            tDecompositionAlgorithm->perform( &tGenerationData, &tDecompositionData, tCutIntegrationMesh.get(), tBackgroundMesh, this );
        }

        tCutIntegrationMesh->finalize_cut_mesh_construction();

        // pick out the cell groups
        moris::Cell< std::shared_ptr< IG_Cell_Group > >& tActiveIgCellGroups = tCutIntegrationMesh->get_all_cell_groups();

        moris::Cell< moris::mtk::Cell* > tActiveIgCells;
        extract_cells_from_cell_groups( tActiveIgCellGroups, tActiveIgCells );

        // create facet connectivity in the mesh
        std::shared_ptr< Facet_Based_Connectivity > tFaceConnectivity = std::make_shared< Facet_Based_Connectivity >();
        this->create_facet_from_element_to_node( tActiveIgCells, tFaceConnectivity );

        // TODO: MESH-CLEANUP
        // this->create_vertex_facet_connectivity(tFaceConnectivity, tCutIntegrationMesh.get() );

        tCutIntegrationMesh->set_face_connectivity( tFaceConnectivity );

        // TODO: MESH-CLEANUP
        // initialize MESH-GEN map
        // moris::Cell<moris_index> tGenMeshMap( tCutIntegrationMesh->get_num_nodes() );
        // std::iota(tGenMeshMap.begin(), tGenMeshMap.end(), 0);
        // mGeometryEngine->get_pdv_host_manager()->set_GenMeshMap(tGenMeshMap);

        // TODO: MESH-CLEANUP
        // perform mesh cleanup
        // Integration_Mesh_Cleanup* MeshCleanup = new Integration_Mesh_Cleanup( tCutIntegrationMesh.get(), tFaceConnectivity.get() );
        // MeshCleanup->perform( tActiveIgCells, mGeometryEngine );


        // set the bulk phase of each cell
        this->compute_ig_cell_bulk_phase( tCutIntegrationMesh.get() );

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

        // print diagnostic information if requested
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

        // if triangulation of all cells is requested, also subdivide the non-cut cells
        if ( mXTKModel->mTriangulateAllInPost )
        {
            MORIS_LOG_WARNING( "Triangulation of non-intersected background cells in post-processing requested." );
            MORIS_LOG_WARNING( "Mesh may be missing some side sets at the domain boundary." );

            // find non-intersected cells to triangulate in post
            this->determine_non_intersected_background_cells( tGenerationData, tCutIntegrationMesh.get(), tBackgroundMesh );

            // create the regular subdivision routine with a factory
            std::shared_ptr< Decomposition_Algorithm > tRegSubAlg =
                    create_decomposition_algorithm( this->determine_reg_subdivision_template(), mXTKModel->get_parameter_list() );

            // swap the list of intersected cells out with list of non-intersected cells to trigger regular subdivision on the non-intersected cells
            Cell< moris_index > tAllIntersectedBgCells = tGenerationData.mAllIntersectedBgCellInds;
            tGenerationData.mAllIntersectedBgCellInds  = tGenerationData.mAllNonIntersectedBgCellInds;

            // perform the regular subdivision on the non-intersected cells
            Decomposition_Data tDecompositionData;
            tRegSubAlg->perform( &tGenerationData, &tDecompositionData, tCutIntegrationMesh.get(), tBackgroundMesh, this );

            // assign list of intersected BG cells back to Mesh Generation Data
            tGenerationData.mAllIntersectedBgCellInds = tAllIntersectedBgCells;

            // re-finalize mesh
            tCutIntegrationMesh->finalize_cut_mesh_construction();
            this->compute_ig_cell_bulk_phase( tCutIntegrationMesh.get() );

            // assign newly generated IG cells to subphases
            this->construct_subphases_on_triangulated_non_cut_cells( &tGenerationData, tCutIntegrationMesh.get(), tBackgroundMesh );
        }

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

        // check if subphase-groups need to be constructed
        if ( this->check_construct_subphase_groups() )
        {
            // get the discretization mesh indices
            Matrix< IndexMat > tDiscretizationMeshIndices = mXTKModel->get_Bspline_mesh_indices();

            // get number of B-spline discretization meshes
            uint tNumBspMeshes = tDiscretizationMeshIndices.numel();

            // reserve memory in cell of b-spline mesh infos with correct number of meshes
            tCutIntegrationMesh->mBsplineMeshInfos.resize( tNumBspMeshes, nullptr );

            // initialize the size of the SPG connectivity
            tCutIntegrationMesh->mSubphaseGroupNeighborhood.resize( tNumBspMeshes );

            // establish B-spline mesh information for every mesh index
            for ( moris::size_t iBspMesh = 0; iBspMesh < tNumBspMeshes; iBspMesh++ )
            {
                // get B-spline mesh index currently treated
                moris_index tMeshIndex = tDiscretizationMeshIndices( iBspMesh );

                // Note: this "new" has its "delete" in the destructor of the cut integration mesh
                tCutIntegrationMesh->mBsplineMeshInfos( iBspMesh ) = new Bspline_Mesh_Info();

                // get the pointer to the Bspline_Mesh_Info currently treated
                Bspline_Mesh_Info* tBsplineMeshInfo = tCutIntegrationMesh->mBsplineMeshInfos( iBspMesh );

                // establish B-spline - Lagrange mesh relation
                this->establish_bspline_mesh_info( tCutIntegrationMesh.get(), tBackgroundMesh, tBsplineMeshInfo, tMeshIndex );

                // construct Subphase groups
                this->construct_subphase_groups( tCutIntegrationMesh.get(), tBsplineMeshInfo, iBspMesh );

                // establish the SP to SPG map
                tBsplineMeshInfo->create_SP_to_SPG_map( tCutIntegrationMesh->get_num_subphases() );

                // collect owned and non-owned SPGs and assign IDs
                this->communicate_subphase_groups( tCutIntegrationMesh.get(), tBsplineMeshInfo, iBspMesh );

                // construct and set Subphase Group Neighborhood
                const std::shared_ptr< Subphase_Neighborhood_Connectivity > tSpgNeighborhoodConnectivity =
                        this->construct_subphase_group_neighborhood( tCutIntegrationMesh.get(), tBsplineMeshInfo, tMeshIndex );
                tCutIntegrationMesh->mSubphaseGroupNeighborhood( iBspMesh ) = tSpgNeighborhoodConnectivity;
            }

            // construct information about which SPGs are present on which base IP cells and their material connectivity wrt. to the various B-spline meshes
            // Note: this information is needed later in the enrichment
            this->construct_SPG_material_connectivity_information( tCutIntegrationMesh.get() );

        }    // end: if SPGs are constructed

        // output cut IG mesh for debugging
        if ( mOutputCutIgMesh )
        {
            tCutIntegrationMesh->write_mesh( "./", "cut_ig_mesh.exo" );
        }

        return tCutIntegrationMesh;
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
    Integration_Mesh_Generator::determine_reg_subdivision_template()
    {
        if ( this->get_spatial_dim() == 2 )
        {
            return Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4;
        }
        else if ( this->get_spatial_dim() == 3 )
        {
            return Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8;
        }
        else
        {
            MORIS_ERROR( false, "Integration_Mesh_Generator::determine_reg_subdivision_template() - Spatial dimension is not 2 or 3." );
            return Subdivision_Method::NO_METHOD;
        }
    }

    // ----------------------------------------------------------------------------------

    enum Subdivision_Method
    Integration_Mesh_Generator::determine_order_elevation_template()
    {
        if ( this->get_spatial_dim() == 2 )
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

        else if ( this->get_spatial_dim() == 3 )
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
            moris::mtk::Mesh*                 aBackgroundMesh )
    {
        Tracer tTracer( "XTK", "Integration_Mesh_Generator", "check_intersected_background_cell_levels ", mXTKModel->mVerboseLevel, 1 );

        bool tFlag = true;

        moris_index tReferenceLevel = 0;
        if ( aMeshGenerationData.mAllIntersectedBgCellInds.size() > 0 )
        {
            tReferenceLevel = aBackgroundMesh->get_mtk_cell( aMeshGenerationData.mAllIntersectedBgCellInds( 0 ) ).get_level();
        }

        for ( size_t iBgCellIndex = 1; iBgCellIndex < aMeshGenerationData.mAllIntersectedBgCellInds.size(); iBgCellIndex++ )
        {
            moris_index tLevel = aBackgroundMesh->get_mtk_cell( aMeshGenerationData.mAllIntersectedBgCellInds( iBgCellIndex ) ).get_level();
            if ( tReferenceLevel != tLevel )
            {
                tFlag = false;
            }
        }

        aCutIntegrationMesh->mSameLevelChildMeshes = all_land( tFlag );
    }

    // ----------------------------------------------------------------------------------

    bool
    Integration_Mesh_Generator::determine_intersected_background_cells(
            Integration_Mesh_Generation_Data& aMeshGenerationData,
            Cut_Integration_Mesh*             aCutIntegrationMesh,
            moris::mtk::Mesh*                 aBackgroundMesh )
    {
        Tracer tTracer( "XTK", "Integration_Mesh_Generator", "Determine intersected background cells", mXTKModel->mVerboseLevel, 1 );
        uint   tNumGeometries = this->get_geom_engine()->get_number_of_geometries();

        uint tNumCells = aBackgroundMesh->get_num_elems();

        aMeshGenerationData.mIntersectedBackgroundCellIndex.resize( tNumGeometries );

        aMeshGenerationData.mAllIntersectedBgCellInds.reserve( tNumCells );

        // Initialize geometric query
        Geometric_Query tGeometricQuery;

        // large coord matrix that I want to keep in scope for a long time avoid copying coordinate all the time.
        tGeometricQuery.set_coordinates_matrix( &aCutIntegrationMesh->mVertexCoordinates );

        tGeometricQuery.set_query_entity_rank( mtk::EntityRank::ELEMENT );

        // reserve memory for list of indices of intersected background cells
        // size estimate: tNumCells / tNumGeometries
        for ( moris::size_t iGeom = 0; iGeom < tNumGeometries; iGeom++ )
        {
            aMeshGenerationData.mIntersectedBackgroundCellIndex( iGeom ).reserve( tNumCells / tNumGeometries );
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
                if ( mXTKModel->get_geom_engine()->is_intersected( iGeom, tGeometricQuery.get_query_entity_to_vertex_connectivity(), tGeometricQuery.get_query_indexed_coordinates() ) )
                {
                    // add background cell to the list for iGEOM
                    aMeshGenerationData.mIntersectedBackgroundCellIndex( iGeom ).push_back( iCell );

                    aMeshGenerationData.mNumChildMeshes++;

                    aMeshGenerationData.mAllIntersectedBgCellInds.push_back( iCell );
                }
                // add to the global list of intersected cells,
                else if ( mXTKModel->triangulate_all() )
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

    bool
    Integration_Mesh_Generator::determine_non_intersected_background_cells(
            Integration_Mesh_Generation_Data& aMeshGenerationData,
            Cut_Integration_Mesh*             aCutIntegrationMesh,
            moris::mtk::Mesh*                 aBackgroundMesh )
    {
        Tracer tTracer( "XTK", "Integration_Mesh_Generator", "Determine intersected background cells", mXTKModel->mVerboseLevel, 1 );
        uint   tNumGeometries = this->get_geom_engine()->get_number_of_geometries();

        // get number of Lagrange elements
        uint tNumCells = aBackgroundMesh->get_num_elems();

        // initialize memory for the list of non-intersected cells
        aMeshGenerationData.mAllNonIntersectedBgCellInds.reserve( tNumCells );

        // Initialize geometric query
        Geometric_Query tGeometricQuery;

        // large coord matrix that I want to keep in scope for a long time avoid copying coordinate all the time.
        tGeometricQuery.set_coordinates_matrix( &aCutIntegrationMesh->mVertexCoordinates );

        tGeometricQuery.set_query_entity_rank( mtk::EntityRank::ELEMENT );

        // iterate through all cells
        for ( moris::uint iCell = 0; iCell < tNumCells; iCell++ )
        {
            // assume the cell is non-cut
            bool tCellIsCut = false;

            // setup geometric query with this current cell information
            tGeometricQuery.set_parent_cell( &aBackgroundMesh->get_mtk_cell( (moris_index)iCell ) );
            tGeometricQuery.set_query_cell( &aBackgroundMesh->get_mtk_cell( (moris_index)iCell ) );

            // iterate through all geometries for current cell and check if it gets cut by any of the geometries
            for ( moris::size_t iGeom = 0; iGeom < tNumGeometries; iGeom++ )
            {
                // set to true if cell is cut by current or any previous geometry
                tCellIsCut = tCellIsCut || mXTKModel->get_geom_engine()->is_intersected(
                                     iGeom,
                                     tGeometricQuery.get_query_entity_to_vertex_connectivity(),
                                     tGeometricQuery.get_query_indexed_coordinates() );
            }

            // if the cell is not cut by any geometry, store the cell's index
            if ( !tCellIsCut )
            {
                aMeshGenerationData.mAllNonIntersectedBgCellInds.push_back( iCell );
            }
        }

        // remove the excess space
        shrink_to_fit_all( aMeshGenerationData.mAllIntersectedBgCellInds );

        // return success
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
        Tracer tTracer( "XTK", "Integration_Mesh_Generator", "Commit IG Cells To Mesh", mXTKModel->mVerboseLevel, 1 );

        // iterate through cells that the decomposition constructed
        moris::uint tNumNewCells = aDecompositionAlgorithm->mNumNewCells;

        MORIS_ERROR(
                aDecompositionAlgorithm->mNewCellToVertexConnectivity.size() == aDecompositionAlgorithm->mNewCellChildMeshIndex.size() &&    //
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
            moris::Cell< std::shared_ptr< IG_Cell_Group > > const & aCellGroups,
            moris::Cell< moris::mtk::Cell* >&                       aCellsInGroups )
    {
        Tracer tTracer( "XTK", "Integration_Mesh_Generator", "Extract cells from groups", mXTKModel->mVerboseLevel, 1 );
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
        // log/trace this function
        Tracer tTracer( "XTK", "Integration_Mesh_Generator", "Compute IG cell bulk phases", mXTKModel->mVerboseLevel, 1 );

        // get an estimate how many IG cells will have issues with vertex-based phase assignment
        uint tNumGeometries                   = mGeometryEngine->get_number_of_geometries();
        uint tNumIgCells                      = aCutIntegrationMesh->get_num_entities( mtk::EntityRank::ELEMENT, 0 );
        uint tExpectedNumNotAssignableIgCells = tNumGeometries * ( tNumIgCells / 1000 ) + 1;

        // initialize list of IG cells whose bulk-phase cannot be determined solely by the nodal Level-Set values
        moris::Cell< moris_index > tIgCellsWithoutBulkPhase( 0 );
        tIgCellsWithoutBulkPhase.reserve( tExpectedNumNotAssignableIgCells );

        // go over all IG cells in the mesh and try to establish the bulk-phase index based on nodal level-set values
        for ( moris::size_t iCell = 0; iCell < tNumIgCells; iCell++ )
        {
            // get the Bulk-phase index for the IG cell
            // moris_index tBulkPhaseIndex = this->deduce_ig_cell_bulk_phase_index( &aCutIntegrationMesh->get_mtk_cell( iCell ) );
            moris_index tBulkPhaseIndex = this->deduce_ig_cell_bulk_phase_from_vertices( &aCutIntegrationMesh->get_mtk_cell( iCell ) );

            // store the bulk-phase of the IG cell
            aCutIntegrationMesh->mIntegrationCellBulkPhase( iCell ) = tBulkPhaseIndex;

            // check if the bulk-phase could actually be assigned, if not...
            if ( tBulkPhaseIndex == MORIS_INDEX_MAX )
            {
                // ... mark the IG cell as unassigned
                tIgCellsWithoutBulkPhase.push_back( iCell );
            }
        }

        // determine whether all IG cells have already been assigned a bulk phase
        uint tNumIgCellsWithoutBulkPhase = tIgCellsWithoutBulkPhase.size();
        bool tAllIgCellsHaveBulkPhase    = ( tNumIgCellsWithoutBulkPhase == 0 );

        // initialize a list of IG cells that still has unassigned bulk-phases after neighbor-based bulk-phase voting
        moris::Cell< moris_index > tIgCellsStillWithoutBulkPhase;

        // keep assigning bulk-phase indices base on neighbors of the IG cells with unassigned bulk-phase
        while ( !tAllIgCellsHaveBulkPhase )
        {
            // log this information
            MORIS_LOG_INFO( "Could not assign bulk-phase to all IG cells. Perform neighbor based assignment." );
            MORIS_LOG_SPEC( "Number of IG cells to assign bulk-phases to", tNumIgCellsWithoutBulkPhase );

            // empty list of Ig cells that need subsequent passes
            tIgCellsStillWithoutBulkPhase.reserve( tNumIgCellsWithoutBulkPhase );
            tIgCellsStillWithoutBulkPhase.resize( 0 );

            // go over the IG cells with unassigned bulk-phases and try to find their bulk-phases based on their neighbors
            for ( moris::size_t iUnassignedCell = 0; iUnassignedCell < tNumIgCellsWithoutBulkPhase; iUnassignedCell++ )
            {
                // get the index for the IG cell with the unassigned bulk-phase
                moris_index tIgCellIndex = tIgCellsWithoutBulkPhase( iUnassignedCell );

                // get the bulk-phase for the IG cell based on the neighbor bulk-phases
                moris_index tBulkPhaseIndex = this->deduce_ig_cell_bulk_phase_from_facets( &aCutIntegrationMesh->get_mtk_cell( tIgCellIndex ), aCutIntegrationMesh );

                // store the bulk-phase of the IG cell
                aCutIntegrationMesh->mIntegrationCellBulkPhase( tIgCellIndex ) = tBulkPhaseIndex;

                // check if the bulk-phase could still not be assigned
                if ( tBulkPhaseIndex == MORIS_INDEX_MAX )
                {
                    // ... mark the IG cell as unassigned
                    tIgCellsStillWithoutBulkPhase.push_back( tIgCellIndex );
                }
            }    // end: loop over unassigned IG cells

            // update the number of un-assigned IG cells
            uint tNumIgCellsAssignedInLoop = tNumIgCellsWithoutBulkPhase - tIgCellsStillWithoutBulkPhase.size();
            tNumIgCellsWithoutBulkPhase    = tNumIgCellsWithoutBulkPhase - tNumIgCellsAssignedInLoop;

            // check that number of IG cells with assigned bulk-phases actually increases and that the code is not stuck in an infinite loop
            MORIS_ERROR( tNumIgCellsAssignedInLoop > 0,
                    "Integration_Mesh_Generator::compute_ig_cell_bulk_phase() - IG cell phase assignment stuck in infinite loop." );

            // check if anything is still unassigned
            tAllIgCellsHaveBulkPhase = ( tNumIgCellsWithoutBulkPhase == 0 );

            // copy over unassigned IG-cells and start over if needed
            tIgCellsWithoutBulkPhase = tIgCellsStillWithoutBulkPhase;
        }

        // Step 1: deduce IG cells bulk-phase indices based on vertex voting system
        //          + collect cells that don't have clear bulk phases in list

        // Step 2: for IG cells that don't have a clear bulk-phase, assign their bulk phase based on facet based neighbor voting
        //          + repeat this until all IG cells have their bulk-phase (encapsuled cells may not have any neighbors that are able to vote on the bulk phase)
        //          + have a check in the while loop that if there are unassigned IG cells but no new IG cells get their bulk-phases assigned, something's wrong
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
    Integration_Mesh_Generator::deduce_ig_cell_bulk_phase_index( moris::mtk::Cell const * aCell )
    {
        // cell vertices
        moris::Cell< moris::mtk::Vertex* > tVertices = aCell->get_vertex_pointers();
        moris::size_t                      tNumGeom  = mGeometryEngine->get_number_of_geometries();

        // allocate phase on or off value (either 0 or 1)
        Matrix< IndexMat > tPhaseVotes( 1, 2 );
        tPhaseVotes.fill( 0 );

        uint                             tMaxRow = 0;
        uint                             tMaxCol = 0;
        moris::Matrix< moris::IndexMat > tNodalPhaseVals( 1, tNumGeom, MORIS_INDEX_MAX );

        for ( moris::uint iGeom = 0; iGeom < tNumGeom; iGeom++ )
        {
            bool tFoundNonInterfaceNode = false;

            for ( moris::uint iVert = 0; iVert < tVertices.size(); iVert++ )
            {
                ge::Geometric_Region tRegion = mGeometryEngine->get_geometric_region( iGeom, tVertices( iVert )->get_index(), tVertices( iVert )->get_coords() );
                if ( tRegion == ge::Geometric_Region::NEGATIVE )
                {
                    tFoundNonInterfaceNode  = true;
                    tPhaseVotes( 0 )++;
                }
                else if ( tRegion == ge::Geometric_Region::POSITIVE )
                {
                    tFoundNonInterfaceNode  = true;
                    tPhaseVotes( 1 )++;
                }
            }    // end: loop over all vertices on IG cell

            // take the phase with the maximum number of votes
            tPhaseVotes.max( tMaxRow, tMaxCol );
            tNodalPhaseVals( 0, iGeom ) = tMaxCol;
            tPhaseVotes.fill( 0 );

            //
            if ( !tFoundNonInterfaceNode )
            {
                MORIS_LOG_WARNING(
                        "IMG::deduce_ig_cell_bulk_phase_index() - Did not find a non-interface node for this element, set to dummy: %zu",
                        mGeometryEngine->get_num_phases() );
                return mGeometryEngine->get_num_phases();
            }

            // MORIS_ERROR( tFoundNonInterfaceNode, "Did not find a non-interface node for this element" );

        }    // end: loop over all level-sets

        return mGeometryEngine->get_elem_phase_index( tNodalPhaseVals );
    }

    // ----------------------------------------------------------------------------------

    moris_index
    Integration_Mesh_Generator::deduce_ig_cell_bulk_phase_from_vertices( moris::mtk::Cell const * aCell )
    {
        // cell vertices
        moris::Cell< moris::mtk::Vertex* > tVertices = aCell->get_vertex_pointers();
        moris::size_t                      tNumGeom  = mGeometryEngine->get_number_of_geometries();

        // allocate phase on or off value (either 0 or 1)
        Matrix< IndexMat > tPhaseVotes( 1, 2 );
        tPhaseVotes.fill( 0 );

        uint                             tMaxRow = 0;
        uint                             tMaxCol = 0;
        moris::Matrix< moris::IndexMat > tNodalPhaseVals( 1, tNumGeom, MORIS_INDEX_MAX );

        for ( moris::uint iGeom = 0; iGeom < tNumGeom; iGeom++ )
        {
            bool tFoundNonInterfaceNode = false;

            for ( moris::uint iVert = 0; iVert < tVertices.size(); iVert++ )
            {
                uint tVertexIndex = tVertices( iVert )->get_index();
                ge::Geometric_Region tRegion = mGeometryEngine->get_geometric_region( iGeom, tVertexIndex, tVertices( iVert )->get_coords() );
                if ( tRegion == ge::Geometric_Region::NEGATIVE )
                {
                    tFoundNonInterfaceNode  = true;
                    tPhaseVotes( 0 )++;
                }
                else if ( tRegion == ge::Geometric_Region::POSITIVE )
                {
                    tFoundNonInterfaceNode  = true;
                    tPhaseVotes( 1 )++;
                }
            }    // end: loop over all vertices on IG cell

            // take the phase with the maximum number of votes
            tPhaseVotes.max( tMaxRow, tMaxCol );
            tNodalPhaseVals( 0, iGeom ) = tMaxCol;
            tPhaseVotes.fill( 0 );

            // do not return bulk-phase if non could be assigned based on the vertex Level-Set values
            if ( !tFoundNonInterfaceNode )
            {
                // MORIS_LOG_SPEC( "Vertex-based phase assignment failed for IG cell index", aCell->get_index() );
                return MORIS_INDEX_MAX;
            }

        }    // end: loop over all level-sets

        // return the dominant phase index
        return mGeometryEngine->get_elem_phase_index( tNodalPhaseVals );
    }

    // ----------------------------------------------------------------------------------

    moris_index
    Integration_Mesh_Generator::deduce_ig_cell_bulk_phase_from_facets(
            moris::mtk::Cell const * aCell,
            Cut_Integration_Mesh*    aCutIntegrationMesh )
    {
        // access the facet connectivity
        std::shared_ptr< Facet_Based_Connectivity > tFaceConn = aCutIntegrationMesh->get_face_connectivity();

        // get the current cell's index
        moris_index tIgCellIndex = aCell->get_index();

        // ignore background cells fed into this function
        // NOTE: QUADs/HEXs should always get their phase from the nodal based phase assignment, as they're not intersected by definition
        // NOTE: Therefore, they have a clear/trivial phase assignment before.
        if ( (uint)aCell->get_geometry_type() == 4 || (uint)aCell->get_geometry_type() == 2 )
        {
            return MORIS_INDEX_MAX - 1;
        }

        // get the current cell's position in the facet connectivity
        auto tMapPos = tFaceConn->mCellIndexToCellOrdinal.find( tIgCellIndex );
        MORIS_ERROR( tMapPos != tFaceConn->mCellIndexToCellOrdinal.end(),
                "Integration_Mesh_Generator::deduce_ig_cell_bulk_phase_from_facets() - "
                "Trying to assign bulk phase (based on neighbors) to IG cell that is neither in the facet connectivity graph, nor a QUAD or HEX." );
        moris_index tIgCellIndexInFacetConn = tMapPos->second;

        // get the facet indices for the current IG cell
        moris::Cell< moris_index > tFacetsOnCell    = tFaceConn->mCellToFacet( tIgCellIndexInFacetConn );
        uint                       tNumFacetsOnCell = tFacetsOnCell.size();

        // ballot on which neighboring IG cells can vote for their bulk phase based on their volume
        Mini_Map< moris_index, real > tBulkPhaseVoteBallot;

        // get the valid IG cells connected to current IG cell
        for ( uint iFacet = 0; iFacet < tNumFacetsOnCell; iFacet++ )
        {
            // get the current facet's index
            moris_index tFacetIndex = tFacetsOnCell( iFacet );

            // get the number of IG cells connected to the current facet
            uint tNumCellsOnFacet = tFaceConn->mFacetToCell( tFacetIndex ).size();

            // if this is not a domain boundary facet: proceed, otherwise: ignore
            if ( tNumCellsOnFacet > 1 )
            {
                // get the indices of the cells connected to the current facet
                moris::mtk::Cell* tFirstCellConnectedToFacet  = tFaceConn->mFacetToCell( tFacetIndex )( 0 );
                moris::mtk::Cell* tSecondCellConnectedToFacet = tFaceConn->mFacetToCell( tFacetIndex )( 1 );

                // store the neighbor cell
                moris::mtk::Cell* tNeighborIgCell;

                // if this cell is not the current cell itself, use the this cell as neighbor
                if ( tFirstCellConnectedToFacet->get_index() != tIgCellIndex )
                {
                    tNeighborIgCell = tFirstCellConnectedToFacet;
                }

                // otherwise the second connected cell is the neighbor cell
                else
                {
                    tNeighborIgCell = tSecondCellConnectedToFacet;
                }

                // check that the neighbor cell has valid bulk-phase already before listing valid neighbor, otherwise ignore
                moris_index tNeighborBulkPhase = aCutIntegrationMesh->mIntegrationCellBulkPhase( tNeighborIgCell->get_index() );
                if ( tNeighborBulkPhase != MORIS_INDEX_MAX )
                {
                    // get the valid IG cell's volume
                    real tCellVolume = tNeighborIgCell->compute_cell_measure();

                    // look for the bulk-phase on the voting ballot
                    auto tIter = tBulkPhaseVoteBallot.find( tNeighborBulkPhase );

                    // if it is not on the ballot add it
                    if ( tIter == tBulkPhaseVoteBallot.end() )
                    {
                        tBulkPhaseVoteBallot[ tNeighborBulkPhase ] = tCellVolume;
                    }

                    // if it is already on the ballot, add weight to the vote
                    else
                    {
                        tIter->second += tCellVolume;
                    }
                }    // end: only consider neighbors that already have a bulk-phase assigned to them

            }    // end: only treat facets with two IG cells connected to it
        }        // end: loop over all facets attached to current IG cell

        // initialize vote counter
        moris_index tWinnerBulkPhase = MORIS_INDEX_MAX;
        real        tWinnerVote      = -1.0;

        // go over ballot and take strongest vote
        for ( auto const & tIter : tBulkPhaseVoteBallot )
        {
            if ( tIter.second > tWinnerVote )
            {
                tWinnerBulkPhase = tIter.first;
                tWinnerVote      = tIter.second;
            }
        }

        // return the bulk-phase
        return tWinnerBulkPhase;
    }

    // ----------------------------------------------------------------------------------

    void
    Integration_Mesh_Generator::deduce_interfaces(
            Cut_Integration_Mesh*                       aCutIntegrationMesh,
            std::shared_ptr< Facet_Based_Connectivity > aFacetConnectivity,
            moris::Cell< moris_index >&                 aInterfaces )
    {
        // log/trace this function call
        Tracer tTracer( "XTK", "Integration_Mesh_Generator", "Deduce Interface", mXTKModel->mVerboseLevel, 1 );

        // get the number of facets in the mesh
        uint tNumFacets = aFacetConnectivity->mFacetToCell.size();

        // all cell groups at this point should be such that an integration cell does not appear twice
        aInterfaces.clear();
        aInterfaces.reserve( tNumFacets );

        // loop over all facets in the whole integration mesh
        for ( moris::uint iFacet = 0; iFacet < tNumFacets; iFacet++ )
        {
            // get the number of elements attached to the current facet
            uint tNumElemsAttachedToFacet = aFacetConnectivity->mFacetToCell( iFacet ).size();

            // check that each facet has, as expected would be expected, only one or two elements attached to it
            MORIS_ASSERT( ( tNumElemsAttachedToFacet > 0 ) && ( tNumElemsAttachedToFacet < 3 ),
                    "Integration_Mesh_Generator::deduce_interfaces() - Facet is expected to have exactly 1 or 2 elements attached to it." );

            // check if facet has two elements attached to it
            if ( tNumElemsAttachedToFacet == 2 )
            {
                // get access to the two elements
                moris::mtk::Cell* tFirstCell  = aFacetConnectivity->mFacetToCell( iFacet )( 0 );
                moris::mtk::Cell* tSecondCell = aFacetConnectivity->mFacetToCell( iFacet )( 1 );

                // get the bulk phase to the two elements
                moris_index tFirstCellBulkPhase  = aCutIntegrationMesh->get_cell_bulk_phase( tFirstCell->get_index() );
                moris_index tSecondCellBulkPhase = aCutIntegrationMesh->get_cell_bulk_phase( tSecondCell->get_index() );

                // check whether the two elements adjacent to the facet have different bulk phases
                if ( tFirstCellBulkPhase != tSecondCellBulkPhase )
                {
                    aInterfaces.push_back( iFacet );
                }
            }

        }    // end: loop over all facets in the mesh

        // size out unused memory
        aInterfaces.shrink_to_fit();
    }


    // Note: this is an old version of the function above, will be removed if new way of deducting interfaces turns out more robust
    // void
    // Integration_Mesh_Generator::deduce_interfaces(
    //     Cut_Integration_Mesh*                       aCutIntegrationMesh,
    //     std::shared_ptr< Facet_Based_Connectivity > aFacetConnectivity,
    //     moris::Cell< moris_index >&                 aInterfaces )
    // {
    //     Tracer tTracer( "XTK", "Integration_Mesh_Generator", "Deduce Interface" ,mXTKModel->mVerboseLevel, 1  );
    //
    //     // all cell groups at this point should be such that an integration cell does not appear twice
    //     aInterfaces.clear();
    //     aInterfaces.reserve( aFacetConnectivity->mFacetVertices.size() );
    //
    //     // loop over all facets in the whole integration mesh
    //     for ( moris::uint iFacet = 0; iFacet < aFacetConnectivity->mFacetVertices.size(); iFacet++ )
    //     {
    //         // loop over all geometries / Level-Sets
    //         for ( moris::uint iGeom = 0; iGeom < mActiveGeometries.numel(); iGeom++ )
    //         {
    //             // get the geometry index of the current Level-Set
    //             moris_index tGeomIndex        = mActiveGeometries( iGeom );
    //             bool        tIsInterfaceFacet = true;
    //
    //             // go over the vertices on the current facet for the current geometry
    //             for ( moris::uint iVert = 0; iVert < aFacetConnectivity->mFacetVertices( iFacet ).size(); iVert++ )
    //             {
    //                 // get access to the vertex being treated
    //                 moris::mtk::Vertex* tVertex = aFacetConnectivity->mFacetVertices( iFacet )( iVert );
    //
    //                 // determine whether the Level-Set is zero at the vertex within snapping tolerance
    //                 if ( !mGeometryEngine->get_geometric_region( tVertex->get_index(), tGeomIndex ) )
    //                 {
    //                     // if any of the vertices on the facet is not on the interface, mark the whole facet to not be an interface and ...
    //                     tIsInterfaceFacet = false;
    //
    //                     // ...skip the rest the vertex for-loop
    //                     break;
    //                 }
    //             }
    //
    //             // list facet as an interface facet ( this is repeated for each geometry )
    //             if ( tIsInterfaceFacet )
    //             {
    //                 aInterfaces.push_back( iFacet );
    //             }
    //
    //         } // end: loop over all Level-Sets
    //     } // end: loop over all facets in the mesh
    // }

    // ----------------------------------------------------------------------------------

    void
    Integration_Mesh_Generator::identify_and_construct_subphases(
            Integration_Mesh_Generation_Data*                 aMeshGenerationData,
            Cut_Integration_Mesh*                             aCutIntegrationMesh,
            moris::mtk::Mesh*                                 aBackgroundMesh,
            std::shared_ptr< Cell_Neighborhood_Connectivity > aCutNeighborhood )
    {
        // trace/log this function
        Tracer tTracer( "XTK", "Integration_Mesh_Generator", "Identify subphases", mXTKModel->mVerboseLevel, 1 );

        // get the number of children meshes (= number of Bg-Cells?)
        moris::uint tNumChildMeshes = aCutIntegrationMesh->get_num_ig_cell_groups();

        moris::uint tNumIgCells = aCutIntegrationMesh->get_num_entities( mtk::EntityRank::ELEMENT, 0 );

        moris::uint tTotalNumBgCells = aBackgroundMesh->get_num_entities( mtk::EntityRank::ELEMENT, 0 );

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

        // initialize list of sub-phases (indices) present on current Bg-cell
        moris::Cell< moris_index > tSubphaseIndices;
        tSubphaseIndices.reserve( 10 );

        // iterate over cut IP cells and perform local flood-fill
        for ( auto& iCell : aMeshGenerationData->mAllIntersectedBgCellInds )
        {
            // get pointer to IG-Cell-Group to currently treated child mesh
            std::shared_ptr< IG_Cell_Group > tIgCellGroup = aCutIntegrationMesh->get_ig_cell_group( iCell );

            // get bg-cell (= parent-cell) associated with current child mesh
            moris::mtk::Cell* tParentCell = aCutIntegrationMesh->get_ig_cell_group_parent_cell( iCell );

            // make sure assumption that intersected Bg-Cell index is equal to child mesh index holds true
            MORIS_ASSERT( tParentCell->get_index() == iCell, "Index mismatch parent index should align with parent cell index" );

            // flood fill this group using the bulk phases
            moris::Matrix< moris::IndexMat > tLocalFloodFill =
                    this->flood_fill_ig_cell_group( aCutIntegrationMesh, aCutNeighborhood, tIgCellGroup, tMaxSubPhase );

            // put first subphase in spot of parent cell (indices for bg-cells are reserved for first subphase within bg-cells)
            aCutIntegrationMesh->mSubPhaseCellGroups( tParentCell->get_index() ) = std::make_shared< IG_Cell_Group >();
            aCutIntegrationMesh->mSubPhaseBulkPhase( tParentCell->get_index() )  = MORIS_INDEX_MAX;
            aCutIntegrationMesh->mSubPhaseParentCell( tParentCell->get_index() ) = tParentCell;
            aCutIntegrationMesh->mSubPhaseIds( tParentCell->get_index() )        = tParentCell->get_id();
            aCutIntegrationMesh->mParentCellToSubphase( tParentCell->get_index() ).push_back( tParentCell->get_index() );
            aCutIntegrationMesh->mParentCellHasChildren( tParentCell->get_index() ) = ( moris_index ) true;
            tSubphaseIndices.push_back( tParentCell->get_index() );

            // allocate the other subphases (additional subphases within cut elements get appended to list of all subphases)
            for ( moris::uint iSP = 1; iSP < (uint)tMaxSubPhase + 1; iSP++ )
            {
                aCutIntegrationMesh->mSubPhaseCellGroups.push_back( std::make_shared< IG_Cell_Group >() );
                aCutIntegrationMesh->mSubPhaseBulkPhase.push_back( MORIS_INDEX_MAX );
                aCutIntegrationMesh->mSubPhaseParentCell.push_back( tParentCell );
                aCutIntegrationMesh->mSubPhaseIds.push_back( MORIS_ID_MAX );
                aCutIntegrationMesh->mParentCellToSubphase( tParentCell->get_index() ).push_back( tSubPhaseIndex );
                tSubphaseIndices.push_back( tSubPhaseIndex++ );    // increment sub-phase index counter to keep track of total number of sub-phases
            }

            // iterate through IG-cells on current Bg-Cell
            for ( moris::uint iIgCell = 0; iIgCell < tLocalFloodFill.numel(); iIgCell++ )
            {
                // get bg-cell-local index of subphase current IG-cell is associated with from flood fill
                moris_index tFFVal = tLocalFloodFill( iIgCell );

                // get (global but proc-local) index of subphase being treated
                moris_index tSPIndex = tSubphaseIndices( tFFVal );

                // put the cell into the subphase group
                aCutIntegrationMesh->mSubPhaseCellGroups( tSPIndex )->mIgCellGroup.push_back( tIgCellGroup->mIgCellGroup( iIgCell ) );

                // store sub-phase index related to integration cell index
                aCutIntegrationMesh->mIntegrationCellToSubphaseIndex( tIgCellGroup->mIgCellGroup( iIgCell )->get_index() ) = tSPIndex;

                // find (based on current IG cell) and store bulk phase index sub-phase belongs to
                if ( aCutIntegrationMesh->mSubPhaseBulkPhase( tSPIndex ) == MORIS_INDEX_MAX )
                {
                    aCutIntegrationMesh->mSubPhaseBulkPhase( tSPIndex ) =
                            aCutIntegrationMesh->get_cell_bulk_phase( tIgCellGroup->mIgCellGroup( iIgCell )->get_index() );
                }
            }

            // store sub-phases present on current Bg-cell to child mesh
            aCutIntegrationMesh->set_child_mesh_subphase( iCell, tSubphaseIndices );

            // clear list of sub-phases (indices) present on current Bg-cell
            tSubphaseIndices.clear();
        }

        // iterate over background cells and make the subphase contain only them
        for ( moris::size_t iBgCell = 0; iBgCell < aBackgroundMesh->get_num_elems(); iBgCell++ )
        {
            // check whether current Bg Cell has already been treated
            // (note: previous loop only treated cut Bg cells, this one treats the remaining non-cut cells)
            if ( aCutIntegrationMesh->mSubPhaseCellGroups( iBgCell ) == nullptr )
            {
                // on non-cut cell, create trivial IG-cell group consisting of single element
                // and store this one to the sub-phase information
                moris::mtk::Cell* tCell                                                = &aCutIntegrationMesh->get_mtk_cell( (moris_index)iBgCell );
                aCutIntegrationMesh->mSubPhaseCellGroups( iBgCell )                    = std::make_shared< IG_Cell_Group >( 1 );
                aCutIntegrationMesh->mSubPhaseCellGroups( iBgCell )->mIgCellGroup( 0 ) = tCell;

                aCutIntegrationMesh->mSubPhaseBulkPhase( iBgCell )  = aCutIntegrationMesh->get_cell_bulk_phase( iBgCell );
                aCutIntegrationMesh->mSubPhaseParentCell( iBgCell ) = tCell;
                aCutIntegrationMesh->mSubPhaseIds( iBgCell )        = tCell->get_id();
                aCutIntegrationMesh->mParentCellToSubphase( iBgCell ).push_back( iBgCell );
                aCutIntegrationMesh->mIntegrationCellToSubphaseIndex( iBgCell ) = iBgCell;
                aCutIntegrationMesh->mParentCellHasChildren( iBgCell )          = ( moris_index ) false;
            }
        }

        // shrink to fit lists with sub-phase information (in case size was over-estimated on initialization)
        aCutIntegrationMesh->mSubPhaseCellGroups.shrink_to_fit();
        aCutIntegrationMesh->mSubPhaseBulkPhase.shrink_to_fit();
        aCutIntegrationMesh->mSubPhaseParentCell.shrink_to_fit();
        aCutIntegrationMesh->mSubPhaseIds.shrink_to_fit();

        // --- PARALLEL-CONSIDERATIONS ---
        // initialize list of owned and non-owned indices
        aCutIntegrationMesh->mOwnedSubphaseGroupsInds.reserve( aCutIntegrationMesh->mSubPhaseCellGroups.size() );
        aCutIntegrationMesh->mNotOwnedSubphaseGroupsInds.reserve( aCutIntegrationMesh->mSubPhaseCellGroups.size() );

        // get current proc rank
        moris_index tParRank = moris::par_rank();

        // go through sub-phases and sort them into ownership groups based on their parent cells
        for ( moris::uint i = 0; i < aCutIntegrationMesh->mSubPhaseCellGroups.size(); i++ )
        {
            if ( aCutIntegrationMesh->mSubPhaseParentCell( i )->get_owner() == tParRank )
            {
                aCutIntegrationMesh->mOwnedSubphaseGroupsInds.push_back( (moris_index)i );
            }
            else
            {
                aCutIntegrationMesh->mNotOwnedSubphaseGroupsInds.push_back( (moris_index)i );
            }
        }

        // shrink to fit sub-phase ownership lists (in case size was over-estimated on initialization)
        aCutIntegrationMesh->mOwnedSubphaseGroupsInds.shrink_to_fit();
        aCutIntegrationMesh->mNotOwnedSubphaseGroupsInds.shrink_to_fit();

        // give all sub-phases a global ID (across all procs)
        this->assign_subphase_glob_ids( aCutIntegrationMesh, aBackgroundMesh );
    }

    // ----------------------------------------------------------------------------------

    void
    Integration_Mesh_Generator::construct_subphases_on_triangulated_non_cut_cells(
            Integration_Mesh_Generation_Data* aMeshGenerationData,
            Cut_Integration_Mesh*             aCutIntegrationMesh,
            moris::mtk::Mesh*                 aBackgroundMesh )
    {
        MORIS_ERROR( mXTKModel->mTriangulateAllInPost,
                "Integration_Mesh_Generator::construct_subphases_on_triangulated_non_cut_cells() - "
                "This function should only be invoked if triangulation of all elements in post has been requested." );

        // get the non-intersected BG cells
        Cell< moris_index > const & tNonIntersectedBgCells = aMeshGenerationData->mAllNonIntersectedBgCellInds;

        // get the number of new cell groups
        uint tNumNewTriangulatedBgCells = tNonIntersectedBgCells.size();

        // mark these BG cells to now also have children
        for ( uint iNonIntersectBgCell = 0; iNonIntersectBgCell < tNumNewTriangulatedBgCells; iNonIntersectBgCell++ )
        {
            // get the index of the non-intersected BG cell
            moris_index tBgCellIndex = tNonIntersectedBgCells( iNonIntersectBgCell );

            // mark the BG cell to have children
            aCutIntegrationMesh->mParentCellHasChildren( tBgCellIndex ) = ( moris_index ) true;

            // get the IG-Cell-Group that has been constructed on the current BG element by the regular subdivision on the non-intersected BG cells
            std::shared_ptr< IG_Cell_Group > tIgCellGroup = aCutIntegrationMesh->get_ig_cell_group( tBgCellIndex );

            // assign the new IG-Cell-Group to the subphase
            aCutIntegrationMesh->mSubPhaseCellGroups( tBgCellIndex )->mIgCellGroup = tIgCellGroup->mIgCellGroup;
        }

        // in debug check that all BG cells have been marked to have children
#ifdef MORIS_HAVE_DEBUG
        for ( uint iBgCell = 0; iBgCell < aCutIntegrationMesh->mParentCellHasChildren.size(); iBgCell++ )
        {
            MORIS_ASSERT( aCutIntegrationMesh->mParentCellHasChildren( iBgCell ) == ( moris_index ) true,
                    "Integration_Mesh_Generator::construct_subphases_on_triangulated_non_cut_cells() - "
                    "BG cell with index %i not marked to have children.",
                    iBgCell );
        }
#endif
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
        // log/trace this function
        Tracer tTracer( "XTK", "Integration_Mesh_Generator", "Subphase Neighborhood", mXTKModel->mVerboseLevel, 1 );

        // get the number of SPs in the mesh
        uint tNumSPs = aCutIntegrationMesh->get_num_subphases();

        // initialize lists in Subphase_Neighborhood_Connectivity
        aSubphaseNeighborhood->mSubphaseToSubPhase.resize( tNumSPs );
        aSubphaseNeighborhood->mSubphaseToSubPhaseMySideOrds.resize( tNumSPs );
        aSubphaseNeighborhood->mSubphaseToSubPhaseNeighborSideOrds.resize( tNumSPs );
        aSubphaseNeighborhood->mTransitionNeighborCellLocation.resize( tNumSPs );

        // allocate the pointers in subphase neighborhood ( for each sub-phase index there's a list )
        for ( moris::uint iSp = 0; iSp < tNumSPs; iSp++ )
        {
            aSubphaseNeighborhood->mSubphaseToSubPhase( iSp )                 = std::make_shared< moris::Cell< moris_index > >();
            aSubphaseNeighborhood->mSubphaseToSubPhaseMySideOrds( iSp )       = std::make_shared< moris::Cell< moris_index > >();
            aSubphaseNeighborhood->mSubphaseToSubPhaseNeighborSideOrds( iSp ) = std::make_shared< moris::Cell< moris_index > >();
            aSubphaseNeighborhood->mTransitionNeighborCellLocation( iSp )     = std::make_shared< moris::Cell< moris_index > >();

            aSubphaseNeighborhood->mSubphaseToSubPhase( iSp )->reserve( 4 );
            aSubphaseNeighborhood->mSubphaseToSubPhaseMySideOrds( iSp )->reserve( 4 );
            aSubphaseNeighborhood->mSubphaseToSubPhaseNeighborSideOrds( iSp )->reserve( 4 );
            aSubphaseNeighborhood->mTransitionNeighborCellLocation( iSp )->reserve( 4 );
        }

        /* Note:
         * I loop through cells then side ordinals to capture the transition between adaptively refined meshes.
         * alternatively you could try looping through aBgFacetToChildFacet and then handle the subphase
         * neighborhood for the transition locations separately
         */

        // get number of Bg-cells
        moris::uint tNumIpCells = aBackgroundMesh->get_num_entities( mtk::EntityRank::ELEMENT );

        // iterate through background cells
        for ( moris_index iIpCell = 0; iIpCell < (moris_index)tNumIpCells; iIpCell++ )
        {
            // get pointer to current Bg-cell
            mtk::Cell const * tCurrentCell = &aBackgroundMesh->get_mtk_cell( iIpCell );

            // get the Bg-cells connected to current Bg-cell and the corresponding facets and side ordinals through which these connections run
            Matrix< IndexMat > tCellToCellSideIndex = aBackgroundMesh->get_elements_connected_to_element_and_face_ind_loc_inds( iIpCell );
            Matrix< IndexMat > tCellToCellSideOrd   = aBackgroundMesh->get_elements_connected_to_element_and_face_ord_loc_inds( iIpCell );

            // get the number of IP cell neighbors
            uint tNumNeighborIpCells = tCellToCellSideOrd.n_cols();

            // iterate through neighboring Bg-cells
            for ( uint iNeighborIpCell = 0; iNeighborIpCell < tNumNeighborIpCells; iNeighborIpCell++ )
            {
                // currently treated neighbor cell
                mtk::Cell const * tOtherCell = &aBackgroundMesh->get_mtk_cell( tCellToCellSideIndex( 0, iNeighborIpCell ) );

                // facet (ordinal) shared for current neighbors
                moris_index tFacetIndex             = tCellToCellSideIndex( 1, iNeighborIpCell );
                moris_index tMyOrdinal              = tCellToCellSideOrd( 1, iNeighborIpCell );
                moris_index tNeighborOrdinal        = tCellToCellSideOrd( 2, iNeighborIpCell );
                moris_index tTransitionCellLocation = tCellToCellSideOrd( 3, iNeighborIpCell );

                // find Ig-cells that are representative for connection from center-Bg-cell through currently treated side ordinal
                Cell< moris_index > tMyCellSubphaseIndices( 0 );    // list of Sub-phases (indices) present on current Lag-elem
                Cell< moris_index > tRepresentativeIgCells( 0 );    // index of single Ig-cell that is representative (what does representative mean?)
                Cell< moris_index > tRepresentativeIgCellsOrdinal( 0 );
                this->collect_subphases_attached_to_facet_on_cell(
                        aCutIntegrationMesh,                         // mesh info
                        tCurrentCell,                                // center cell 'connection from'
                        tMyOrdinal,                                  // side ordinal of center cell the connection is coming from
                        tFacetIndex,                                 // index of facet through which connection goes
                        aFacetConnectivity,                          // facet connectivity info
                        ( *aBgFacetToChildFacet )( tFacetIndex ),    // further facet connectivity info
                        tMyCellSubphaseIndices,                      // to be filled
                        tRepresentativeIgCells,                      // to be filled
                        tRepresentativeIgCellsOrdinal );             // to be filled

                // transitioning between mesh levels (i.e. there's a hierarchical refinement boundary at the currently treated facet)
                if ( tTransitionCellLocation != MORIS_INDEX_MAX )
                {
                    Matrix< IndexMat > tNeighborCellToFacetIndex =
                            aBackgroundMesh->get_entity_connected_to_entity_loc_inds( tOtherCell->get_index(),
                                    mtk::EntityRank::ELEMENT,
                                    aBackgroundMesh->get_facet_rank() );
                    Cell< moris::moris_index > tNeighborSubphaseIndices( 0 );
                    Cell< moris::moris_index > tNeighborRepresentativeIgCells( 0 );
                    Cell< moris::moris_index > tNeighborRepresentativeIgCellsOrdinal( 0 );
                    this->collect_subphases_attached_to_facet_on_cell(
                            aCutIntegrationMesh,
                            tOtherCell,
                            tNeighborOrdinal,
                            tNeighborCellToFacetIndex( tNeighborOrdinal ),
                            aFacetConnectivity,
                            ( *aBgFacetToChildFacet )( tNeighborCellToFacetIndex( tNeighborOrdinal ) ),
                            tNeighborSubphaseIndices,
                            tNeighborRepresentativeIgCells,
                            tNeighborRepresentativeIgCellsOrdinal );

                    for ( moris::uint i = 0; i < tMyCellSubphaseIndices.size(); i++ )
                    {
                        moris_index tMySubphaseIndex = tMyCellSubphaseIndices( i );
                        moris_index tMyBulkIndex     = aCutIntegrationMesh->get_subphase_bulk_phase( tMySubphaseIndex );
                        for ( moris::uint j = 0; j < tNeighborSubphaseIndices.size(); j++ )
                        {
                            moris_index tNeighborBulkIndex     = aCutIntegrationMesh->get_subphase_bulk_phase( tNeighborSubphaseIndices( j ) );
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

                // trivial case: no transition between refinement levels
                else
                {
                    // iterate over sub-phases of center-Bg-cell relevant to current facet
                    for ( moris::uint iMySP = 0; iMySP < tMyCellSubphaseIndices.size(); iMySP++ )
                    {
                        // temporarily store sub-phase info for convenience
                        moris_index tMySubphaseIndex  = tMyCellSubphaseIndices( iMySP );
                        moris_index tMyIgCellIndex    = tRepresentativeIgCells( iMySP );
                        moris_index tMyIgCellSideOrd  = tRepresentativeIgCellsOrdinal( iMySP );
                        moris_index tMyBulkPhaseIndex = aCutIntegrationMesh->get_subphase_bulk_phase( tMySubphaseIndex );

                        // case: transition between background cell and another background cell or triangulated cells
                        // i.e. at least one of the Bg-cells is NOT cut
                        if ( !aCutIntegrationMesh->parent_cell_has_children( tCurrentCell->get_index() )    //
                                || !aCutIntegrationMesh->parent_cell_has_children( tOtherCell->get_index() ) )
                        {
                            // find Ig-cells that are representative for the currently treated neighbor-Bg-cell and its connection
                            Cell< moris::moris_index > tNeighborSubphaseIndices( 0 );
                            Cell< moris::moris_index > tNeighborRepresentativeIgCells( 0 );
                            Cell< moris::moris_index > tNeighborRepresentativeIgCellsOrdinal( 0 );
                            this->collect_subphases_attached_to_facet_on_cell(
                                    aCutIntegrationMesh,    // mesh info
                                    tOtherCell,             // neighbor cell 'connection to'
                                    tNeighborOrdinal,       // side ordinal of neighbor cell the connection is going to
                                    tFacetIndex,            // index of facet through which connection goes
                                    aFacetConnectivity,     // facet connectivity info
                                    ( *aBgFacetToChildFacet )( tFacetIndex ),
                                    tNeighborSubphaseIndices,
                                    tNeighborRepresentativeIgCells,
                                    tNeighborRepresentativeIgCellsOrdinal );

                            // iterate through neighbor Bg-cell's sub-phases connected to current facet and store them in the connectivity
                            // note: in the case that the transition includes one uncut element this should always be only one "my sub-phase" and one neighbor sub-phase
                            for ( const auto& iNeighSp : tNeighborSubphaseIndices )
                            {
                                // get the neighbors subphase's bulk-phase
                                moris_index tNeighborBulkPhaseIndex = aCutIntegrationMesh->get_subphase_bulk_phase( iNeighSp );

                                // check that the neighboring IG cell is of the same bulk-phase
                                if ( tNeighborBulkPhaseIndex == tMyBulkPhaseIndex )
                                {
                                    aSubphaseNeighborhood->mSubphaseToSubPhase( tMySubphaseIndex )->push_back( iNeighSp );
                                    aSubphaseNeighborhood->mSubphaseToSubPhaseMySideOrds( tMySubphaseIndex )->push_back( tMyOrdinal );
                                    aSubphaseNeighborhood->mSubphaseToSubPhaseNeighborSideOrds( tMySubphaseIndex )->push_back( tNeighborOrdinal );
                                    aSubphaseNeighborhood->mTransitionNeighborCellLocation( tMySubphaseIndex )->push_back( tTransitionCellLocation );
                                }
                            }
                        }

                        // case: both Bg-cells are cut
                        else
                        {
                            // get the index of the treated transition facet between the two IG cells considered
                            const moris_index& tMyIgCellOrdInFacetConn = aFacetConnectivity->get_cell_ordinal( tMyIgCellIndex );
                            moris_index        tIgFacetIndex           = aFacetConnectivity->mCellToFacet( tMyIgCellOrdInFacetConn )( tMyIgCellSideOrd );

                            // initialize variable
                            moris_index tNeighborSubphaseIndex = MORIS_INDEX_MAX;

                            // iterate through cells on facet and get the one that is not my cell, take the sub-phase index of this one
                            for ( const auto& iCell : aFacetConnectivity->mFacetToCell( tIgFacetIndex ) )
                            {
                                if ( iCell->get_index() != tMyIgCellIndex )
                                {
                                    // get sub-phase index from neighbor IG-cell
                                    tNeighborSubphaseIndex = aCutIntegrationMesh->get_ig_cell_subphase_index( iCell->get_index() );
                                }
                            }

                            // check subphase's validity
                            MORIS_ASSERT( tNeighborSubphaseIndex != MORIS_INDEX_MAX,
                                    "Integration_Mesh_Generator::construct_subphase_neighborhood() - "
                                    "Adjacent IG Cell not found." );

                            //// MORIS_ASSERT( aCutIntegrationMesh->get_subphase_bulk_phase( tNeighborSubphaseIndex ) == aCutIntegrationMesh->get_subphase_bulk_phase( tMySubphaseIndex ), "Subphase bulk phase mismatch" );

                            // get the neighboring subphases' bulk-phase
                            moris_index tNeighborBulkPhaseIndex = aCutIntegrationMesh->get_subphase_bulk_phase( tNeighborSubphaseIndex );

                            // if the connected sub-phases also have the same bulk-phase index, then ...
                            if ( tNeighborBulkPhaseIndex == tMyBulkPhaseIndex )
                            {
                                // ... register the sub-phase connection in SP Neighborhood Connection
                                aSubphaseNeighborhood->mSubphaseToSubPhase( tMySubphaseIndex )->push_back( tNeighborSubphaseIndex );
                                aSubphaseNeighborhood->mSubphaseToSubPhaseMySideOrds( tMySubphaseIndex )->push_back( tMyOrdinal );
                                aSubphaseNeighborhood->mSubphaseToSubPhaseNeighborSideOrds( tMySubphaseIndex )->push_back( tNeighborOrdinal );
                                aSubphaseNeighborhood->mTransitionNeighborCellLocation( tMySubphaseIndex )->push_back( tTransitionCellLocation );
                            }
                        }    // end if: case that both IP cells are cut
                    }        // end for: loop over subphases on current base IP cell
                }            // end if: case that there's no transition between refinement levels
            }                // end for: loop over the base IP cell's respective neighbors
        }                    // end for: loop over all base IP cells
    }                        // end function

    // ----------------------------------------------------------------------------------

    void
    Integration_Mesh_Generator::collect_subphases_attached_to_facet_on_cell(
            Cut_Integration_Mesh*                                aCutIntegrationMesh,
            moris::mtk::Cell const *                             aBGCell,
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

        // case: cluster is non-trivial
        if ( tParentHasChildren )
        {
            // initialize map that ?...
            std::unordered_map< moris_index, moris_index > tSubphaseMap;

            // make sure that a cut Bg-element has information which IG-cell facets comprise the Bg-cell facet
            MORIS_ASSERT( aBgFacetToChildrenFacets != nullptr,
                    "Integration_Mesh_Generator::collect_subphases_attached_to_facet_on_cell() - Null ptr on facet that should have children facets" );

            // get number of IG cells connected to current IP cell's facet
            uint tNumConnectedIgCells = aBgFacetToChildrenFacets->size();

            // iterate through child facets attached to bg facet
            for ( moris::uint iChildFacet = 0; iChildFacet < tNumConnectedIgCells; iChildFacet++ )
            {
                // get currently treated Ig-cell facet
                moris_index tChildCellFacetIndex = ( *aBgFacetToChildrenFacets )( iChildFacet );

                // get number of IG cells connected to current IG cell facet
                uint tNumIgCellsOnFacet = aFacetConnectivity->mFacetToCell( tChildCellFacetIndex ).size();

                // go over the Ig-cells attached to the current Ig-cell facet
                for ( moris::uint iCell = 0; iCell < tNumIgCellsOnFacet; iCell++ )
                {
                    // get index of current Ig-cell, the sub-phase it belongs to, and its parent cell
                    moris_index tCellIndex       = aFacetConnectivity->mFacetToCell( tChildCellFacetIndex )( iCell )->get_index();
                    moris_index tSubphaseIndex   = aCutIntegrationMesh->get_ig_cell_subphase_index( tCellIndex );
                    moris_index tParentCellIndex = aCutIntegrationMesh->get_subphase_parent_cell( tSubphaseIndex )->get_index();

                    // make sure current Ig-cell's parent and the treated Bg-cell are the same cell (otherwise different Ig-cell must be representative for current Bg-cell)
                    if ( tParentCellIndex == aBGCell->get_index() )
                    {
                        // if this sub-phase has not been registered in the sub-phase map yet, then ...
                        if ( tSubphaseMap.find( tSubphaseIndex ) == tSubphaseMap.end() )
                        {
                            // ... register this sub-phase as being connected to currently treated facet
                            aSubphaseIndices.push_back( tSubphaseIndex );

                            // store current Ig-cell and its facet as representative
                            aRepresentativeIgCells.push_back( tCellIndex );
                            aRepresentativeIgCellsOrdinal.push_back( aFacetConnectivity->mFacetToCellEdgeOrdinal( tChildCellFacetIndex )( iCell ) );

                            // register sub-phase index as being connected to current facet
                            tSubphaseMap[ tSubphaseIndex ] = 1;
                        }
                    }
                }
            }
        }

        // case: cluster is trivial
        else
        {
            // check that exactly one subphase has been assigned to trivial non-cut cell
            MORIS_ASSERT( aCutIntegrationMesh->get_parent_cell_subphases( aBGCell->get_index() ).size() == 1,
                    "Integration_Mesh_Generator::collect_subphases_attached_to_facet_on_cell() - one subphase needs to be present in this case" );

            // simply use Bg-cell information for everything, since trivial
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

        moris::Cell< std::shared_ptr< moris::Cell< moris::moris_index > > > const & tBGFacetToChildFacet = aCutIntegrationMesh->get_background_facet_to_child_facet_connectivity();


        if ( tBGFacetToChildFacet( aBackgroundFacetIndex ) == nullptr )
        {
            return;
        }

        else
        {
            // facet connectivity
            std::shared_ptr< Facet_Based_Connectivity > tFacetConn = aCutIntegrationMesh->get_face_connectivity();

            moris::Cell< moris::moris_index > const & tBGFacetToIgFacet = *tBGFacetToChildFacet( aBackgroundFacetIndex );

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
    // ----------------------------------------------------------------------------------

    void
    Integration_Mesh_Generator::assign_subphase_glob_ids(
            Cut_Integration_Mesh* aCutIntegrationMesh,
            moris::mtk::Mesh*     aBackgroundMesh )
    {
        // log this function when verbose output is requested
        Tracer tTracer( "XTK", "Integration Mesh Generator", "assign subphase IDs", mXTKModel->mVerboseLevel, 1 );

        // get the communication table
        Matrix< IdMat > tCommTable = aCutIntegrationMesh->get_communication_table();

        /* ---------------------------------------------------------------------------------------- */
        /* Step 1: Let each proc decide how many entity IDs it needs & communicate ID ranges */

        // Subphase IDs up to the number of IP cells have already been used. Hence, the first free ID is:
        moris_index tFirstSubphaseId = aBackgroundMesh->get_max_entity_id( mtk::EntityRank::ELEMENT ) + 1;

        // Get the number of subphases (on the current proc)
        moris_id tNumSubphases = (moris_id)aCutIntegrationMesh->get_num_subphases();

        // reserve IDs for his proc
        moris_id tMyFirstId = get_processor_offset( tNumSubphases ) + tFirstSubphaseId;

        /* ---------------------------------------------------------------------------------------- */
        /* Step 2: Assign IDs to owned entities */

        // assign IDs to owned sub-phases
        this->assign_IDs_to_owned_subphases( aCutIntegrationMesh, tMyFirstId );

        /* ---------------------------------------------------------------------------------------- */
        /* The following steps are only necessary if code runs in parallel */

        if ( par_size() == 1 )    // serial
        {
            // check that all IG cells are owned in serial
            MORIS_ASSERT( aCutIntegrationMesh->get_not_owned_subphase_indices().size() == 0,
                    "Integration_Mesh_Generator::assign_subphase_glob_ids() - "
                    "Code running in serial, but not all sub-phases are owned by proc 0." );
        }
        else    // parallel
        {

            /* ---------------------------------------------------------------------------------------- */
            /* Step 3: Prepare requests for non-owned entities */

            // initialize lists of information that identifies sub-phases (on other procs)
            Cell< Cell< moris_index > > tNotOwnedSubphasesToProcs;    // sub-phase index (local to current proc, just used for construction of arrays)
            Cell< Matrix< IdMat > >     tParentCellIds;               // IDs of the sub-phases' parent cells
            Cell< Matrix< IdMat > >     tFirstChildIgCellIds;         // IDs of the first IG cell in the IG cell groups corresponding to the SPs
            Cell< Matrix< IdMat > >     tNumChildCellsInSubphases;    // Number of IG cells in each of the sub-phases

            // fill information
            this->prepare_requests_for_not_owned_subphase_IDs(
                    aCutIntegrationMesh,
                    tNotOwnedSubphasesToProcs,
                    tParentCellIds,
                    tFirstChildIgCellIds,
                    tNumChildCellsInSubphases );

            /* ---------------------------------------------------------------------------------------- */
            /* Step 4: Send and Receive requests about non-owned entities to and from other procs */

            // initialize arrays for receiving
            Cell< Matrix< IdMat > > tReceivedParentCellIds;
            Cell< Matrix< IdMat > > tReceivedFirstChildIgCellIds;
            Cell< Matrix< IdMat > > tReceivedNumChildCellsInSubphases;

            // communicate information
            moris::communicate_mats( tCommTable, tParentCellIds, tReceivedParentCellIds );
            moris::communicate_mats( tCommTable, tFirstChildIgCellIds, tReceivedFirstChildIgCellIds );
            moris::communicate_mats( tCommTable, tNumChildCellsInSubphases, tReceivedNumChildCellsInSubphases );

            // clear memory not needed anymore
            tParentCellIds.clear();
            tFirstChildIgCellIds.clear();
            tNumChildCellsInSubphases.clear();

            /* ---------------------------------------------------------------------------------------- */
            /* Step 5: Find answers to the requests */

            // initialize lists of ID answers to other procs
            Cell< Matrix< IdMat > > tSubphaseIds;

            // prepare answers to the ID requests
            this->prepare_answers_for_owned_subphase_IDs(
                    aCutIntegrationMesh,
                    tSubphaseIds,
                    tReceivedParentCellIds,
                    tReceivedFirstChildIgCellIds,
                    tReceivedNumChildCellsInSubphases );

            // clear memory from requests (the answers to which have been found)
            tReceivedParentCellIds.clear();
            tReceivedFirstChildIgCellIds.clear();
            tReceivedNumChildCellsInSubphases.clear();

            /* ---------------------------------------------------------------------------------------- */
            /* Step 6: Send and receive answers to and from other procs */

            // initialize arrays for receiving
            Cell< Matrix< IdMat > > tReceivedSubphaseIds;

            // communicate answers
            moris::communicate_mats( tCommTable, tSubphaseIds, tReceivedSubphaseIds );

            // clear unused memory
            tSubphaseIds.clear();

            /* ---------------------------------------------------------------------------------------- */
            /* Step 7: Use answers to assign IDs to non-owned entities */

            this->handle_requested_subphase_ID_answers(
                    aCutIntegrationMesh,
                    tNotOwnedSubphasesToProcs,
                    tReceivedSubphaseIds );

        }    // end if: parallel

        // setup map relating the sub-phase indices on this processor to their IDs
        aCutIntegrationMesh->setup_glob_to_loc_subphase_map();

    }    // end function: Cut_Integration_Mesh::assign_subphase_glob_ids()

    // ----------------------------------------------------------------------------------

    void
    Integration_Mesh_Generator::assign_IDs_to_owned_subphases(
            Cut_Integration_Mesh* aCutIntegrationMesh,
            moris_id              aFirstID )
    {
        // access the list of sub-phases owned by the current proc
        moris::Cell< moris_index >& tOwnedSubphases = aCutIntegrationMesh->get_owned_subphase_indices();

        // assign  IDs to these sub-phases
        for ( moris::size_t iSP = 0; iSP < tOwnedSubphases.size(); iSP++ )
        {
            // index wrt. the current proc of the treated sub-phase
            moris_index tSubphaseIndex = tOwnedSubphases( iSP );

            // if the Sub-phase has not already been assigned
            if ( aCutIntegrationMesh->mSubPhaseIds( tSubphaseIndex ) == MORIS_ID_MAX )
            {
                aCutIntegrationMesh->mSubPhaseIds( tSubphaseIndex ) = aFirstID++;
            }
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Integration_Mesh_Generator::prepare_requests_for_not_owned_subphase_IDs(
            Cut_Integration_Mesh*        aCutIntegrationMesh,
            Cell< Cell< moris_index > >& aNotOwnedSubphasesToProcs,
            Cell< Matrix< IdMat > >&     aParentCellIds,
            Cell< Matrix< IdMat > >&     aFirstChildIgCellIds,
            Cell< Matrix< IdMat > >&     aNumChildCellsInSubphases )
    {
        // get the communication table and map
        Matrix< IdMat >                   tCommTable              = aCutIntegrationMesh->get_communication_table();
        uint                              tCommTableSize          = tCommTable.numel();
        std::map< moris_id, moris_index > tProcIdToCommTableIndex = aCutIntegrationMesh->get_communication_map();

        // initialize lists of identifying information
        aNotOwnedSubphasesToProcs.resize( tCommTableSize );
        aParentCellIds.resize( tCommTableSize );
        aFirstChildIgCellIds.resize( tCommTableSize );
        aNumChildCellsInSubphases.resize( tCommTableSize );

        // get the non-owned Subphases on the executing processor
        moris::Cell< moris_index >& tNotOwnedSubphases = aCutIntegrationMesh->get_not_owned_subphase_indices();

        // go through SPs that executing proc knows about, but doesn't own, ...
        for ( uint iNotOwnedSP = 0; iNotOwnedSP < tNotOwnedSubphases.size(); iNotOwnedSP++ )
        {
            // ... get their respective owners and position in the comm-table ...
            moris_index tOwnerProc = aCutIntegrationMesh->get_subphase_parent_cell( tNotOwnedSubphases( iNotOwnedSP ) )->get_owner();
            auto        tIter      = tProcIdToCommTableIndex.find( tOwnerProc );
            MORIS_ASSERT(
                    tIter != tProcIdToCommTableIndex.end(),
                    "Integration_Mesh_Generator::assign_subphase_glob_ids() - "
                    "Sub-phase owner (Proc #%i) not found in communication table of current proc #%i which is: %s",
                    tOwnerProc,
                    par_rank(),
                    ios::stringify_log( tCommTable ).c_str() );
            moris_index tProcDataIndex = tIter->second;

            // ... and finally add the non-owned SP in the list of SPs to be requested from that owning proc
            aNotOwnedSubphasesToProcs( tProcDataIndex ).push_back( tNotOwnedSubphases( iNotOwnedSP ) );
        }

        // iterate through processors
        for ( uint iProc = 0; iProc < tCommTableSize; iProc++ )
        {
            // get the number of non-owned subphases to be sent to each processor
            uint tNumSubphases = aNotOwnedSubphasesToProcs( iProc ).size();

            // allocate matrix
            aParentCellIds( iProc ).resize( 1, tNumSubphases );
            aFirstChildIgCellIds( iProc ).resize( 1, tNumSubphases );
            aNumChildCellsInSubphases( iProc ).resize( 1, tNumSubphases );

            // go through the sub-phases for which IDs will need to be answered by the other processor
            for ( uint iSP = 0; iSP < tNumSubphases; iSP++ )
            {
                // get the index of the sub-phase on the executing proc
                moris_index tSubphaseIndex = aNotOwnedSubphasesToProcs( iProc )( iSP );

                // get the subphase's IG cell group
                std::shared_ptr< IG_Cell_Group > tIgCellGroup = aCutIntegrationMesh->get_subphase_ig_cells( tSubphaseIndex );

                // get the subphase's parent IP cell
                mtk::Cell* tParentCell = aCutIntegrationMesh->get_subphase_parent_cell( tSubphaseIndex );

                // store the parent cell and IG cell group IDs, and the number of IG cells corresponding to the SPs to be requested in the remaining arrays
                aParentCellIds( iProc )( iSP )            = tParentCell->get_id();
                aFirstChildIgCellIds( iProc )( iSP )      = tIgCellGroup->mIgCellGroup( 0 )->get_id();
                aNumChildCellsInSubphases( iProc )( iSP ) = tIgCellGroup->mIgCellGroup.size();
            }
        }

        // size out unused memory
        aNotOwnedSubphasesToProcs.shrink_to_fit();
        aParentCellIds.shrink_to_fit();
        aFirstChildIgCellIds.shrink_to_fit();
        aNumChildCellsInSubphases.shrink_to_fit();
    }

    // ----------------------------------------------------------------------------------

    void
    Integration_Mesh_Generator::prepare_answers_for_owned_subphase_IDs(
            Cut_Integration_Mesh*           aCutIntegrationMesh,
            Cell< Matrix< IdMat > >&        aSubphaseIds,
            Cell< Matrix< IdMat > > const & aReceivedParentCellIds,
            Cell< Matrix< IdMat > > const & aReceivedFirstChildIgCellIds,
            Cell< Matrix< IdMat > > const & aReceivedNumChildCellsInSubphases )
    {
        // get the communication table and map
        Matrix< IdMat >                   tCommTable              = aCutIntegrationMesh->get_communication_table();
        uint                              tCommTableSize          = tCommTable.numel();
        std::map< moris_id, moris_index > tProcIdToCommTableIndex = aCutIntegrationMesh->get_communication_map();

        // initialize lists of ID answers to other procs
        aSubphaseIds.resize( tCommTableSize );

        // check that the received data is complete
        MORIS_ASSERT(
                aReceivedParentCellIds.size() == tCommTableSize && aReceivedFirstChildIgCellIds.size() == tCommTableSize && aReceivedNumChildCellsInSubphases.size() == tCommTableSize,
                "Integration_Mesh_Generator::prepare_answers_for_owned_subphase_IDs() - Received information incomplete." );

        // go through the list of processors in the array of ID requests
        for ( uint iProc = 0; iProc < tCommTableSize; iProc++ )
        {
            // get the number of SP IDs requested from the current proc position
            uint tNumReceivedReqs = aReceivedParentCellIds( iProc ).n_cols();

            // size the list of answers / IDs accordingly
            aSubphaseIds( iProc ).resize( 1, tNumReceivedReqs );

            // iterate through SPs for which the IDs are requested
            for ( uint iSP = 0; iSP < tNumReceivedReqs; iSP++ )
            {
                // get the first IG cell's ID in requested SP
                moris_id tFirstChildCellId = aReceivedFirstChildIgCellIds( iProc )( iSP );

                // get the IG cell's index on this processor
                moris_index tFirstChildCellIndex = aCutIntegrationMesh->get_loc_entity_ind_from_entity_glb_id( tFirstChildCellId, mtk::EntityRank::ELEMENT );

                // find the subphase index on this processor
                moris_index tSubphaseIndex = aCutIntegrationMesh->get_ig_cell_subphase_index( tFirstChildCellIndex );

                // get the subphase's ID
                moris_index tSubphaseId = aCutIntegrationMesh->get_subphase_id( tSubphaseIndex );

                // check that the SP found on this proc has the same parent IP cell as the one from the requesting proc
                MORIS_ERROR( aCutIntegrationMesh->get_subphase_parent_cell( tSubphaseIndex )->get_id() == aReceivedParentCellIds( iProc )( 0, iSP ),
                        "Integration_Mesh_Generator::prepare_answers_for_owned_subphase_IDs() - Subphase parent cell id discrepancy" );

                // check that the SP found on this proc has the same number of IG cells in it as the one from the requesting proc
                MORIS_ERROR( (moris_index)aCutIntegrationMesh->get_subphase_ig_cells( tSubphaseIndex )->mIgCellGroup.size() == aReceivedNumChildCellsInSubphases( iProc )( iSP ),
                        "Integration_Mesh_Generator::prepare_answers_for_owned_subphase_IDs() - Number of cells in subphase discrepancy" );

                // check that a valid Subphase ID has already been assigned
                MORIS_ERROR( tSubphaseId != MORIS_ID_MAX, "Integration_Mesh_Generator::prepare_answers_for_owned_subphase_IDs() - Subphase ID not found in child mesh" );

                // place Subphase ID in return data
                aSubphaseIds( iProc )( iSP ) = tSubphaseId;

            }    // end for: communicated sub-phases from current proc

        }    // end for: communication list for each processor
    }

    // ----------------------------------------------------------------------------------

    void
    Integration_Mesh_Generator::handle_requested_subphase_ID_answers(
            Cut_Integration_Mesh*               aCutIntegrationMesh,
            Cell< Cell< moris_index > > const & aNotOwnedSubphasesToProcs,
            Cell< Matrix< IdMat > > const &     aReceivedSubphaseIds )
    {
        // process answers from each proc communicated with
        for ( uint iProc = 0; iProc < aReceivedSubphaseIds.size(); iProc++ )
        {
            // get the number of requests and answers from the current proc
            uint tNumReceivedReqs = aReceivedSubphaseIds( iProc ).numel();

            MORIS_ERROR( tNumReceivedReqs == aNotOwnedSubphasesToProcs( iProc ).size(),
                    "Integration_Mesh_Generator::handle_requested_subphase_ID_answers() - Number of received and requested subphase IDs don't match." );

            // assign IDs to each communicated Sub-Phase
            for ( uint iSP = 0; iSP < tNumReceivedReqs; iSP++ )
            {
                moris_index tSubphaseIndex                          = aNotOwnedSubphasesToProcs( iProc )( iSP );
                aCutIntegrationMesh->mSubPhaseIds( tSubphaseIndex ) = aReceivedSubphaseIds( iProc )( iSP );
            }
        }
    }

    // ----------------------------------------------------------------------------------
    // ----------------------------------------------------------------------------------

    void
    Integration_Mesh_Generator::assign_subphase_group_glob_ids(
            Cut_Integration_Mesh* aCutIntegrationMesh,
            Bspline_Mesh_Info*    aBsplineMeshInfo )
    {
        // log this function when verbose output is requested
        Tracer tTracer( "XTK", "Integration Mesh Generator", "assign subphase group IDs", mXTKModel->mVerboseLevel, 1 );

        // get the communication table and map
        Matrix< IdMat >                   tCommTable              = aCutIntegrationMesh->get_communication_table();
        std::map< moris_id, moris_index > tProcIdToCommTableIndex = aCutIntegrationMesh->get_communication_map();

        /* ---------------------------------------------------------------------------------------- */
        /* Step 1: Let each proc decide how many entity IDs it needs & communicate ID ranges */

        // Get the number of subphase groups (on the current proc)
        moris_id tNumSubphaseGroups = aBsplineMeshInfo->get_num_SPGs();

        // reserve IDs for his proc
        moris_id tMyFirstId = get_processor_offset( tNumSubphaseGroups );

        /* ---------------------------------------------------------------------------------------- */
        /* Step 2: Assign IDs to owned entities */

        // assign IDs to owned sub-phase groups (recall: this function is repeated for each B-spline mesh and therefor B-spline mesh info)
        aBsplineMeshInfo->assign_owned_subphase_group_ids( tMyFirstId );

        /* ---------------------------------------------------------------------------------------- */
        /* The following steps are only necessary if code runs in parallel */

        if ( par_size() == 1 )    // serial
        {
            // check that all entities are owned in serial
            MORIS_ASSERT( aBsplineMeshInfo->mNotOwnedSubphaseGroupIndices.size() == 0,
                    "Integration_Mesh_Generator::assign_subphase_group_glob_ids() - "
                    "Code running in serial, but not all sub-phases groups are owned by proc 0." );
        }
        else    // parallel
        {

            // check that NOT all entities are owned in parallel
            MORIS_ASSERT( aBsplineMeshInfo->mNotOwnedSubphaseGroupIndices.size() > 0,
                    "Integration_Mesh_Generator::assign_subphase_group_glob_ids() - "
                    "Code running in parallel, but all sub-phases groups are owned by current proc #%i.",
                    par_rank() );

            /* ---------------------------------------------------------------------------------------- */
            /* Step 3: Prepare requests for non-owned entities */

            // initialize lists of information that identifies entities (on other procs)
            Cell< Cell< moris_index > >    tNotOwnedSpgsToProcs;    // SPG index (local to current proc, just used for construction of arrays)
            Cell< moris::Matrix< IdMat > > tSubphaseIds;            // first SP IDs in SPGs in each of the SPGs
            Cell< moris::Matrix< IdMat > > tNumSpsInSpg;            // Number of sub-phases in each SPG

            // fill identifying information
            this->prepare_requests_for_not_owned_subphase_group_IDs(
                    aCutIntegrationMesh,
                    aBsplineMeshInfo,
                    tNotOwnedSpgsToProcs,
                    tSubphaseIds,
                    tNumSpsInSpg );

            /* ---------------------------------------------------------------------------------------- */
            /* Step 4: Send and Receive requests about non-owned entities to and from other procs */

            // initialize arrays for receiving
            Cell< Matrix< IdMat > > tReceivedSubphaseIds;
            Cell< Matrix< IdMat > > tReceivedNumSpsInSpg;

            // communicate information
            moris::communicate_mats( tCommTable, tSubphaseIds, tReceivedSubphaseIds );
            moris::communicate_mats( tCommTable, tNumSpsInSpg, tReceivedNumSpsInSpg );

            // clear memory not needed anymore
            tSubphaseIds.clear();
            tNumSpsInSpg.clear();

            /* ---------------------------------------------------------------------------------------- */
            /* Step 5: Find answers to the requests */

            // initialize lists of ID answers to other procs
            Cell< Matrix< IdMat > > tSubphaseGroupIds;

            // answer requests from other procs
            this->prepare_answers_for_owned_subphase_group_IDs(
                    aCutIntegrationMesh,
                    aBsplineMeshInfo,
                    tSubphaseGroupIds,
                    tReceivedSubphaseIds,
                    tReceivedNumSpsInSpg );

            // clear memory from requests (the answers to which have been found)
            tReceivedSubphaseIds.clear();
            tReceivedNumSpsInSpg.clear();

            /* ---------------------------------------------------------------------------------------- */
            /* Step 6: Send and receive answers to and from other procs */

            // initialize arrays for receiving
            Cell< Matrix< IdMat > > tReceivedSubphaseGroupIds;

            // communicate answers
            moris::communicate_mats( tCommTable, tSubphaseGroupIds, tReceivedSubphaseGroupIds );

            // clear unused memory
            tSubphaseGroupIds.clear();

            /* ---------------------------------------------------------------------------------------- */
            /* Step 7: Use answers to assign IDs to non-owned entities */

            this->handle_requested_subphase_group_ID_answers(
                    aCutIntegrationMesh,
                    aBsplineMeshInfo,
                    tNotOwnedSpgsToProcs,
                    tReceivedSubphaseGroupIds );

        }    // end if: parallel

    }    // end function: Cut_Integration_Mesh::assign_subphase_group_glob_ids()

    // ----------------------------------------------------------------------------------

    void
    Integration_Mesh_Generator::prepare_requests_for_not_owned_subphase_group_IDs(
            Cut_Integration_Mesh*        aCutIntegrationMesh,
            Bspline_Mesh_Info*           aBsplineMeshInfo,
            Cell< Cell< moris_index > >& aNotOwnedSpgsToProcs,
            Cell< Matrix< IdMat > >&     aSubphaseIds,
            Cell< Matrix< IdMat > >&     aNumSpsInSpg )
    {
        // get the communication table and map
        Matrix< IdMat >                   tCommTable              = aCutIntegrationMesh->get_communication_table();
        uint                              tCommTableSize          = tCommTable.numel();
        std::map< moris_id, moris_index > tProcIdToCommTableIndex = aCutIntegrationMesh->get_communication_map();

        // initialize lists of identifying information
        aNotOwnedSpgsToProcs.resize( tCommTableSize );
        aSubphaseIds.resize( tCommTableSize );
        aNumSpsInSpg.resize( tCommTableSize );

        // get the non-owned Subphases on the executing processor
        moris::Cell< moris_index >& tNotOwnedSPGs    = aBsplineMeshInfo->mNotOwnedSubphaseGroupIndices;
        uint                        tNumNotOwnedSPGs = tNotOwnedSPGs.size();

        // go through SPGs that executing proc knows about, but doesn't own, ...
        for ( uint iNotOwnedSPG = 0; iNotOwnedSPG < tNumNotOwnedSPGs; iNotOwnedSPG++ )
        {
            // ... get their index ...
            moris_index tSpgIndex = tNotOwnedSPGs( iNotOwnedSPG );

            // ... get their respective owners, and position in the comm table ...
            moris_index tOwnerProc = aBsplineMeshInfo->mSubphaseGroups( tSpgIndex )->get_owner();
            auto        tIter      = tProcIdToCommTableIndex.find( tOwnerProc );
            MORIS_ASSERT(
                    tIter != tProcIdToCommTableIndex.end(),
                    "Integration_Mesh_Generator::prepare_requests_for_not_owned_subphase_group_IDs() - "
                    "Entity owner (Proc #%i) not found in communication table of current proc #%i which is: %s",
                    tOwnerProc,
                    par_rank(),
                    ios::stringify_log( tCommTable ).c_str() );
            moris_index tProcDataIndex = tIter->second;

            // ... and finally add the non-owned SPGs in the list of SPs to be requested from that owning proc
            aNotOwnedSpgsToProcs( tProcDataIndex ).push_back( tNotOwnedSPGs( iNotOwnedSPG ) );
        }

        // assemble identifying information for every processor communicated with
        for ( uint iProc = 0; iProc < tCommTableSize; iProc++ )
        {
            // get the number of non-owned SPGs to be sent to each processor processor
            uint tNumNotOwnedSpgsOnProc = aNotOwnedSpgsToProcs( iProc ).size();

            // allocate matrix
            aSubphaseIds( iProc ).resize( 1, tNumNotOwnedSpgsOnProc );
            aNumSpsInSpg( iProc ).resize( 1, tNumNotOwnedSpgsOnProc );

            // go through the Subphase groups for which IDs will be requested by the other processor
            for ( uint iSPG = 0; iSPG < tNumNotOwnedSpgsOnProc; iSPG++ )
            {
                // get the index of the subphase group on the executing proc
                moris_index tSpgIndex = aNotOwnedSpgsToProcs( iProc )( iSPG );

                // get the SPs in the Group
                const moris::Cell< moris_index >& tSpsInGroup    = aBsplineMeshInfo->mSubphaseGroups( tSpgIndex )->get_SP_indices_in_group();
                moris_index                       tSubphaseIndex = tSpsInGroup( 0 );

                // store the identifying information of the Subphase group in the output arrays
                aSubphaseIds( iProc )( iSPG ) = aCutIntegrationMesh->get_subphase_id( tSubphaseIndex );
                aNumSpsInSpg( iProc )( iSPG ) = tSpsInGroup.size();
            }
        }

        // size out unused memory
        aNotOwnedSpgsToProcs.shrink_to_fit();
        aSubphaseIds.shrink_to_fit();
        aNumSpsInSpg.shrink_to_fit();
    }

    // ----------------------------------------------------------------------------------

    void
    Integration_Mesh_Generator::prepare_answers_for_owned_subphase_group_IDs(
            Cut_Integration_Mesh*           aCutIntegrationMesh,
            Bspline_Mesh_Info*              aBsplineMeshInfo,
            Cell< Matrix< IdMat > >&        aSubphaseGroupIds,
            Cell< Matrix< IdMat > > const & aReceivedSubphaseIds,
            Cell< Matrix< IdMat > > const & aReceivedNumSpsInSpg )
    {
        // get the communication table and map
        Matrix< IdMat >                   tCommTable              = aCutIntegrationMesh->get_communication_table();
        uint                              tCommTableSize          = tCommTable.numel();
        std::map< moris_id, moris_index > tProcIdToCommTableIndex = aCutIntegrationMesh->get_communication_map();

        // initialize array of answers
        aSubphaseGroupIds.resize( tCommTableSize );

        // check that the received data is complete
        MORIS_ASSERT(
                aReceivedSubphaseIds.size() == tCommTableSize && aReceivedNumSpsInSpg.size() == tCommTableSize,
                "Integration_Mesh_Generator::prepare_answers_for_owned_subphase_group_IDs() - Received information incomplete." );

        // go through the list of processors in the array of ID requests
        for ( uint iProc = 0; iProc < tCommTableSize; iProc++ )
        {
            // get the number of entity IDs requested from the current proc position
            uint tNumReceivedReqs = aReceivedSubphaseIds( iProc ).numel();

            // size the list of answers / IDs accordingly
            aSubphaseGroupIds( iProc ).resize( 1, tNumReceivedReqs );

            // iterate through SPs for which the IDs are requested
            for ( uint iSPG = 0; iSPG < tNumReceivedReqs; iSPG++ )
            {
                // get the ID and index of the received SP
                moris_id    tSubphaseId    = aReceivedSubphaseIds( iProc )( iSPG );
                moris_index tSubphaseIndex = aCutIntegrationMesh->get_subphase_index( tSubphaseId );

                // get the subphase group's index that the subphase belongs to
                moris_index tSpgIndex = aBsplineMeshInfo->mSpToSpgMap( tSubphaseIndex );

                // get access to the Subphase group
                Subphase_Group const * tSubphaseGroup = aBsplineMeshInfo->mSubphaseGroups( tSpgIndex );

                // get the subphase group's ID
                moris_index tSubphaseGroupId = tSubphaseGroup->get_id();

                // check validity of subphase group's ID
                MORIS_ASSERT( tSubphaseGroupId != MORIS_ID_MAX,
                        "Integration_Mesh_Generator::prepare_answers_for_owned_subphase_group_IDs() - "
                        "Trying to request index from a subphase group whose ID has not been assigned yet." );

                // check that the SP found on this proc has the same number of IG cells in it as the one from the requesting proc
                MORIS_ERROR( (moris_index)tSubphaseGroup->get_num_SPs_in_group() == aReceivedNumSpsInSpg( iProc )( iSPG ),
                        "Integration_Mesh_Generator::prepare_answers_for_owned_subphase_group_IDs() - Number of cells in subphase discrepancy" );

                // place Subphase ID in return data
                aSubphaseGroupIds( iProc )( iSPG ) = tSubphaseGroupId;

            }    // end for: communication for each entity with current processor

        }    // end for: communication list for each processor
    }

    // ----------------------------------------------------------------------------------

    void
    Integration_Mesh_Generator::handle_requested_subphase_group_ID_answers(
            Cut_Integration_Mesh*               aCutIntegrationMesh,
            Bspline_Mesh_Info*                  aBsplineMeshInfo,
            Cell< Cell< moris_index > > const & aNotOwnedSpgsToProcs,
            Cell< Matrix< IdMat > > const &     aReceivedSubphaseGroupIds )
    {
        // process answers from each proc communicated with
        for ( uint iProc = 0; iProc < aReceivedSubphaseGroupIds.size(); iProc++ )
        {
            // get the number of requests and answers from the current proc
            uint tNumReceivedSpgIds = aReceivedSubphaseGroupIds( iProc ).numel();

            // check that the
            MORIS_ASSERT( tNumReceivedSpgIds = aNotOwnedSpgsToProcs( iProc ).size(),
                    "Integration_Mesh_Generator::handle_requested_subphase_group_ID_answers() - "
                    "Number of SPG ID requests and answers are not the same." );

            // assign IDs to each communicated entity
            for ( uint iSPG = 0; iSPG < tNumReceivedSpgIds; iSPG++ )
            {
                // get the current SPG index and ID from the data provided
                moris_index tSubphaseGroupIndex = aNotOwnedSpgsToProcs( iProc )( iSPG );
                moris_id    tSubphaseGroupId    = aReceivedSubphaseGroupIds( iProc )( iSPG );

                // store the received SPG ID
                aBsplineMeshInfo->mSubphaseGroups( tSubphaseGroupIndex )->set_id( tSubphaseGroupId );
            }
        }
    }

    // ----------------------------------------------------------------------------------
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

        // Number of elements within Lag-Bg-Element for the flood fill
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
        IndexMap tElementToLocalIndex;

        for ( moris::size_t iE = 0; iE < tNumElements; iE++ )
        {
            tElementToLocalIndex[ aIgCellGroup->mIgCellGroup( iE )->get_index() ] = iE;
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
                    tNeighborIndex = ( *aCutNeighborhood->mNeighborCells( tElementIndex ) )( iN )->get_index();

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
                    tNeighborOrd = tElementToLocalIndex[ tActiveFrontElement ];

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
                            tElementIndex = ( *aCutNeighborhood->mNeighborCells( tActiveFrontElement ) )( i )->get_index();

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
        // log/ trace this function
        Tracer tTracer( "XTK", "Integration_Mesh_Generator", "Bulk phase to bulk phase interface", mXTKModel->mVerboseLevel, 1 );

        // get the total number of bulk phases
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
                // std::cout << "Warning: interface case not handled " << std::endl;
                MORIS_LOG_WARNING( "IMG::construct_bulk_phase_to_bulk_phase_interface() - Interface case not handled" );
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
        Tracer tTracer( "XTK", "Integration_Mesh_Generator", "Bulk phase to bulk phase double sided interface", mXTKModel->mVerboseLevel, 1 );

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
        MORIS_ASSERT( aGeomIndex < (moris_index)mGeometryEngine->get_number_of_geometries(), "Geometry index out of bounds" );
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
        Tracer tTracer( "XTK", "Integration_Mesh_Generator", "Construct interface sets", mXTKModel->mVerboseLevel, 1 );
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

        // iterate and add the side set groups to cut mesh
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
        // log/trace this function
        Tracer tTracer( "XTK", "Integration_Mesh_Generator", "Construct bulk phase cell groups", mXTKModel->mVerboseLevel, 1 );

        // get number of bulk phases +1 for the err block
        uint tNumBulkPhases = mGeometryEngine->get_num_bulk_phase() + 1;

        // initialize list of groups of IG cells in one bulk phase
        aBulkPhaseCellGroups.resize( tNumBulkPhases, nullptr );

        // iterate through the controlled integration cells to avoid inactive bg cells
        for ( moris::uint iCell = 0; iCell < aCutIntegrationMesh->mControlledIgCells.size(); iCell++ )
        {
            // get pointer to current mtk IG cell
            moris::mtk::Cell* tCell = aCutIntegrationMesh->mControlledIgCells( iCell ).get();

            // cell bulk phase
            moris_index tBulkPhase = aCutIntegrationMesh->get_cell_bulk_phase( tCell->get_index() );

            // check that bulk phase index is used
            if ( tBulkPhase != MORIS_INDEX_MAX )
            {
                // create new IG cell group object for new bulk phase
                if ( aBulkPhaseCellGroups( tBulkPhase ) == nullptr )
                {
                    aBulkPhaseCellGroups( tBulkPhase ) = std::make_shared< IG_Cell_Group >( 0 );
                }

                // add current IG cell to list of IG cells in bulk phase
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
        // log/trace this function
        Tracer tTracer( "XTK", "Integration_Mesh_Generator", "Construct bulk phase blocks", mXTKModel->mVerboseLevel, 1 );

        // determine the side set names
        uint tNumBulkPhases = mGeometryEngine->get_num_bulk_phase();

        MORIS_ERROR( aBulkPhaseCellGroups.size() == tNumBulkPhases + 1,
                "Integration_Mesh_Generator::construct_bulk_phase_blocks() - We expect there to be a bulk phase cell group for each bulk phase and 1 for problematic cells" );

        // decide on cell topology of integration elements based on number of spatial dimensions
        mtk::CellTopology tCellTopo = xtk::determine_cell_topology(
                aCutIntegrationMesh->get_spatial_dim(),
                aCutIntegrationMesh->mXTKModel->ig_element_order(),
                mtk::CellShape::SIMPLEX );

        // iterate through and construct the names of the blocks
        std::string         tBlockBaseName = "cutblock";
        Cell< std::string > tBlockSetNames( tNumBulkPhases + 1 );
        for ( moris::uint iBP = 0; iBP < tNumBulkPhases; iBP++ )
        {
            tBlockSetNames( iBP ) = tBlockBaseName + "_p_" + std::to_string( iBP );
        }

        tBlockSetNames( tNumBulkPhases ) = "err_block";

        Cell< moris_index > tBlockSetOrds = aCutIntegrationMesh->register_block_set_names( tBlockSetNames, tCellTopo );

        // iterate and add the side set groups to cut mesh
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
    Integration_Mesh_Generator::create_vertex_facet_connectivity(
            std::shared_ptr< Facet_Based_Connectivity > aFaceConnectivity,
            Cut_Integration_Mesh*                       aCutIntegrationMesh )
    {
        // vertex to facet connectivity
        aFaceConnectivity->mVertexFacets.resize( aCutIntegrationMesh->mIntegrationVertices.size() );

        // check that facet connectivity has been constructed (which it isn't in some unit test)
        if ( aFaceConnectivity->mFacetVertices.size() == 0 )
        {
            // if it isn't, don't construct this part of the connectivity and exit the function
            return;
        }

        moris::uint tNumNodesPerFacet = aFaceConnectivity->mFacetVertices( 0 ).size();

        for ( moris::uint iFacet = 0; iFacet < aFaceConnectivity->mFacetVertices.size(); iFacet++ )
        {
            for ( moris::uint iVert = 0; iVert < tNumNodesPerFacet; iVert++ )
            {
                moris_index tVertIndex = aFaceConnectivity->mFacetVertices( iFacet )( iVert )->get_index();
                aFaceConnectivity->mVertexFacets( tVertIndex ).push_back( iFacet );
            }
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Integration_Mesh_Generator::create_facet_from_element_to_node(
            moris::Cell< moris::mtk::Cell* >&           aCells,
            std::shared_ptr< Facet_Based_Connectivity > aFaceConnectivity )
    {
        Tracer tTracer( "XTK", "Integration_Mesh_Generator", "Creating Facets" );

        // Note: this function assumes that all cells are of same type, i.e. have the same cell info

        if ( aCells.size() > 0 )
        {
            // cell information
            moris::mtk::Cell_Info const * tCellInfo = aCells( 0 )->get_cell_info();

            // node to facet map
            moris::Matrix< moris::IdMat > tElementToFacetMap = tCellInfo->get_node_to_facet_map();

            // list of unique vertices
            moris::uint                                    tNumNodes = 0;
            moris::Cell< moris::mtk::Vertex* >             tVertices;
            std::unordered_map< moris_index, moris_index > tVertexIndexToLocalIndexMap;

            // reserve memory for list of vertex points
            tVertices.reserve( aCells.size() );

            // this stack of loops counts the number of unique vertices in the mesh
            // go over all IG cells in the mesh
            for ( moris::uint iCell = 0; iCell < aCells.size(); iCell++ )
            {
                // get access to the list of vertices on the current IG cell
                moris::Cell< moris::mtk::Vertex* > tCellVerts = aCells( iCell )->get_vertex_pointers();

                // go over all vertices on the current IG cell
                for ( moris::uint iVertex = 0; iVertex < tCellVerts.size(); iVertex++ )
                {
                    // check if vertex has already been assigned local vertex index
                    if ( tVertexIndexToLocalIndexMap.find( tCellVerts( iVertex )->get_index() ) == tVertexIndexToLocalIndexMap.end() )
                    {
                        // add local vertex index to map
                        tVertexIndexToLocalIndexMap[ tCellVerts( iVertex )->get_index() ] = (moris_index)tNumNodes;

                        // store vertex pointer
                        tVertices.push_back( tCellVerts( iVertex ) );

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

            // go through all IG cells on mesh
            for ( moris::uint iCell = 0; iCell < aCells.size(); iCell++ )
            {
                // check that cell exits only once
                MORIS_ERROR( aFaceConnectivity->mCellIndexToCellOrdinal.find( aCells( iCell )->get_index() ) == aFaceConnectivity->mCellIndexToCellOrdinal.end(),
                        "Integration_Mesh_Generator::create_facet_from_element_to_node() - Duplicate cell in aCells provided" );

                // set local cell index in map
                aFaceConnectivity->mCellIndexToCellOrdinal[ aCells( iCell )->get_index() ] = iCell;

                // get list of vertex pointers of cell
                moris::Cell< moris::mtk::Vertex* > tCellVerts = aCells( iCell )->get_vertex_pointers();

                // reserve memory for storing facets on cell
                aFaceConnectivity->mCellToFacet( iCell ).reserve( tNumFacesPerElem );

                // iterate through edges of cell
                for ( moris::uint iFacet = 0; iFacet < tElementToFacetMap.n_rows(); iFacet++ )
                {
                    // get the vertices on the facet
                    for ( moris::uint iVertOnFace = 0; iVertOnFace < tElementToFacetMap.n_cols(); iVertOnFace++ )
                    {
                        tVerticesOnFacet( iVertOnFace ) = tCellVerts( tElementToFacetMap( iFacet, iVertOnFace ) );
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

                        // reserve memory for list storing cells and their side ordinal connected to facet
                        // guess of cell size: number of facets per element
                        aFaceConnectivity->mFacetToCell.back().reserve( tNumFacesPerElem );
                        aFaceConnectivity->mFacetToCellEdgeOrdinal.back().reserve( tNumFacesPerElem );

                        // store facet index with vertex
                        auto tIter = tVertexIndexToLocalIndexMap.find( tVerticesOnFacet( 0 )->get_index() );
                        MORIS_ERROR( tIter != tVertexIndexToLocalIndexMap.end(),
                                "Integration_Mesh_Generator::create_facet_from_element_to_node() - "
                                "Invalid vertex detected." );
                        moris_index tLocalVertexIndex = tIter->second;

                        tVertexToFacetIndex( tLocalVertexIndex ).push_back( tFacetIndex );
                    }

                    // store facet index on cell
                    aFaceConnectivity->mCellToFacet( iCell ).push_back( tFacetIndex );

                    // store cell and cell ordinal with facet
                    // if needed increase cell capacity by increments of number of faces per element
                    aFaceConnectivity->mFacetToCell( tFacetIndex ).push_back( aCells( iCell ), tNumFacesPerElem );
                    aFaceConnectivity->mFacetToCellEdgeOrdinal( tFacetIndex ).push_back( iFacet, tNumFacesPerElem );
                }
            }

            // trim inner and outer cells
            shrink_to_fit_all( aFaceConnectivity->mFacetVertices );
            shrink_to_fit_all( aFaceConnectivity->mFacetToCell );
            shrink_to_fit_all( aFaceConnectivity->mFacetToCellEdgeOrdinal );
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Integration_Mesh_Generator::generate_cell_neighborhood(
            moris::Cell< moris::mtk::Cell* >&                 aCells,               // list of mtk::Cells (pointers) that are active on processor's mesh
            std::shared_ptr< Facet_Based_Connectivity >       aFaceConnectivity,    // connectivity of Facets on processor's mesh
            std::shared_ptr< Cell_Neighborhood_Connectivity > aNeighborhood )       // to be filled ...
    {
        Tracer tTracer( "XTK", "Integration_Mesh_Generator", "Generate Neighborhood", mXTKModel->mVerboseLevel, 1 );

        // Initialize Sizes and Variables used in routine
        moris_index tMaxIndex            = this->get_max_index( aCells );
        moris_index tMaxNumElementToFace = 2;

        // Initialize Element to Element with size number of elements x number of faces per element filled with a dummy value.
        aNeighborhood->mNeighborCells.resize( tMaxIndex + 1 );
        aNeighborhood->mMySideOrdinal.resize( tMaxIndex + 1 );
        aNeighborhood->mNeighborSideOrdinal.resize( tMaxIndex + 1 );

        // loop over list of mtk::Cells active on current proc
        for ( moris::uint i = 0; i < (moris::uint)tMaxIndex + 1; i++ )
        {
            // initialize size of Lists within Lists in the Cell Neighborhood Connectivity
            aNeighborhood->mNeighborCells( i )       = std::make_shared< moris::Cell< moris::mtk::Cell* > >();
            aNeighborhood->mMySideOrdinal( i )       = std::make_shared< moris::Cell< moris_index > >();
            aNeighborhood->mNeighborSideOrdinal( i ) = std::make_shared< moris::Cell< moris_index > >();

            // note: the reserving of two entries is only an estimate, more mtk::Cells may be attached to a given mtk::Cell
            aNeighborhood->mNeighborCells( i )->reserve( tMaxNumElementToFace );
            aNeighborhood->mMySideOrdinal( i )->reserve( tMaxNumElementToFace );
            aNeighborhood->mNeighborSideOrdinal( i )->reserve( tMaxNumElementToFace );
        }

        // iterate through facets in cell connectivity. any cell that shares a facet is a neighbor
        for ( moris::uint iF = 0; iF < aFaceConnectivity->mFacetToCell.size(); iF++ )
        {
            // iterate through cells attached to this facet (either just 1 or 2)
            MORIS_ASSERT( aFaceConnectivity->mFacetToCell( iF ).size() == 1 || aFaceConnectivity->mFacetToCell( iF ).size() == 2,
                    "Facet should either connect to no cell or one other cell" );

            // only do something if facet connects two mtk::Cells
            if ( aFaceConnectivity->mFacetToCell( iF ).size() == 2 )
            {
                // get the pair of indices of mtk::Cells connected through facet with index iF
                moris_index tCellIndex0 = aFaceConnectivity->mFacetToCell( iF )( 0 )->get_index();
                moris_index tCellIndex1 = aFaceConnectivity->mFacetToCell( iF )( 1 )->get_index();

                // for each of the two mtk::Cells ...
                // ... store the respective other mtk::Cell (as a mtk::Cell connected to it)
                aNeighborhood->mNeighborCells( tCellIndex0 )->push_back( aFaceConnectivity->mFacetToCell( iF )( 1 ) );
                aNeighborhood->mNeighborCells( tCellIndex1 )->push_back( aFaceConnectivity->mFacetToCell( iF )( 0 ) );

                // ... store the respective side ordinal through which it connects to the other mtk::Cell
                aNeighborhood->mMySideOrdinal( tCellIndex0 )->push_back( aFaceConnectivity->mFacetToCellEdgeOrdinal( iF )( 0 ) );
                aNeighborhood->mMySideOrdinal( tCellIndex1 )->push_back( aFaceConnectivity->mFacetToCellEdgeOrdinal( iF )( 1 ) );

                // fixme: this can be deleted ?!
                // ... store the side ordinal of the respective other mtk::Cell through which the other mtk::Cell connects to the first mtk::Cell
                aNeighborhood->mNeighborSideOrdinal( tCellIndex0 )->push_back( aFaceConnectivity->mFacetToCellEdgeOrdinal( iF )( 1 ) );
                aNeighborhood->mNeighborSideOrdinal( tCellIndex1 )->push_back( aFaceConnectivity->mFacetToCellEdgeOrdinal( iF )( 0 ) );
            }
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Integration_Mesh_Generator::create_edges_from_element_to_node(
            moris::Cell< moris::mtk::Cell* >           aCells,
            std::shared_ptr< Edge_Based_Connectivity > aEdgeConnectivity )
    {
        Tracer tTracer( "XTK", "Integration_Mesh_Generator", "Creating Edges", mXTKModel->mVerboseLevel, 1 );

        // Note: this function assumes that all cells are of same type, i.e. have the same cell info

        if ( aCells.size() > 0 )
        {
            // cell information
            moris::mtk::Cell_Info const * tCellInfo = aCells( 0 )->get_cell_info();

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
            moris::uint                        tNumNodes = 0;
            moris::Cell< moris::mtk::Vertex* > tVertices;

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
                        tVertexIndexToLocalIndexMap[ tCellVerts( iV )->get_index() ] = (moris_index)tNumNodes;

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
                        aEdgeConnectivity->mEdgeToCell.push_back( moris::Cell< moris::mtk::Cell* >(), tIncNumEdges );
                        aEdgeConnectivity->mEdgeToCellEdgeOrdinal.push_back( moris::Cell< moris::moris_index >(), tIncNumEdges );

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
            moris::Cell< moris::mtk::Cell* > const &    aParentCellForDeduction,
            std::shared_ptr< Facet_Based_Ancestry >     aIgFacetAncestry )
    {
        Tracer tTracer( "XTK", "Integration_Mesh_Generator", "Facet Ancestry", mXTKModel->mVerboseLevel, 1 );

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

            moris::mtk::Cell_Info const * tParentCellInfo = aParentCellForDeduction( iFacet )->get_cell_info();

            // vertices of the facet
            moris::Cell< moris::mtk::Vertex* > const & tFacetVertices = aIgCellGroupFacetConnectivity->mFacetVertices( iFacet );

            // get the parent of these vertices from the mesh
            tFacetVertexParentInds.resize( tFacetVertices.size() );
            tFacetVertexParentRanks.resize( tFacetVertices.size() );
            tFacetVertexParentOrds.resize( tFacetVertices.size() );

            // all we need to figure out is if the facet is a sub-facet or a interior to the cell in this case.

            // iterate through vertices
            for ( moris::uint iV = 0; iV < tFacetVertices.size(); iV++ )
            {
                tFacetVertexParentInds( iV )  = aCutIntegrationMesh->mIgVertexParentEntityIndex( tFacetVertices( iV )->get_index() );
                tFacetVertexParentRanks( iV ) = aCutIntegrationMesh->mIgVertexParentEntityRank( tFacetVertices( iV )->get_index() );

                // interior background cell is the parent
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
                                    mtk::EntityRank::ELEMENT,
                                    moris::mtk::get_entity_rank_from_index( tFacetVertexParentRanks( iV ) ) );
                }

                // figure out the ordinal
                for ( moris::uint iE = 0; iE < tEntityConnectedToParent( tFacetVertexParentRanks( iV ) )->numel(); iE++ )
                {
                    if ( ( *tEntityConnectedToParent( tFacetVertexParentRanks( iV ) ) )( iE ) == tFacetVertexParentInds( iV ) )
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
                                    mtk::EntityRank::ELEMENT,
                                    aCutIntegrationMesh->get_facet_rank() );
                }
                // mark the entity parent as the deduction cell
                aIgFacetAncestry->mFacetParentEntityIndex( iFacet )                    = ( *tEntityConnectedToParent( (uint)aCutIntegrationMesh->get_facet_rank() ) )( tCandidateFacetOrds( 0 ) );
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

                // std::cout<<" tParentFacet = "<<tParentIndex<<" | IgFacet = "<<iFacet<<std::endl;
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
            moris::Cell< moris::mtk::Cell* > const &   aParentCellForDeduction,
            std::shared_ptr< Edge_Based_Ancestry >     aIgEdgeAncestry )
    {
        Tracer tTracer( "XTK", "Integration_Mesh_Generator", "Edge Ancestry", mXTKModel->mVerboseLevel, 1 );
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

        mtk::EntityRank tFacetRank = aCutIntegrationMesh->get_facet_rank();

        // iterate through edges in the edge connectivity
        for ( moris::uint iEdge = 0; iEdge < tNumEdges; iEdge++ )
        {
            // vertices of the edge
            moris::Cell< moris::mtk::Vertex* > const & tEdgeVertices = aIgCellGroupEdgeConnectivity->mEdgeVertices( iEdge );

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
                            mtk::EntityRank::ELEMENT,
                            moris::mtk::get_entity_rank_from_index( *tMaxIter ) );

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
                        mtk::EntityRank::ELEMENT,
                        moris::mtk::get_entity_rank_from_index( *tMaxIter ) );

                Matrix< IndexMat > tEntitiesConnectedToBaseCellMinRank = aBackgroundMesh->get_entity_connected_to_entity_loc_inds(
                        aParentCellForDeduction( iEdge )->get_index(),
                        mtk::EntityRank::ELEMENT,
                        moris::mtk::get_entity_rank_from_index( *tMinIter ) );

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
                        mtk::EntityRank::ELEMENT,
                        moris::mtk::get_entity_rank_from_index( tParentOrdinalAndRank( 1 ) ) );

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

                // minimum entity is always a vertex so i can ask the cell
                moris_index tVertexOrdinal =
                        aParentCellForDeduction( iEdge )->get_vertex_ordinal_wrt_cell( tEdgeVertices( tMinIndex )->get_index() );

                // not as convenient to get the edge/facet ordinal
                Matrix< IndexMat > tEntitiesConnectedToBaseCell = aBackgroundMesh->get_entity_connected_to_entity_loc_inds(
                        aParentCellForDeduction( iEdge )->get_index(),
                        mtk::EntityRank::ELEMENT,
                        moris::mtk::get_entity_rank_from_index( *tMaxIter ) );

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
                        mtk::EntityRank::ELEMENT,
                        moris::mtk::get_entity_rank_from_index( tParentOrdinalAndRank( 1 ) ) );

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
        Tracer tTracer( "XTK", "Integration_Mesh_Generator", "Commit Ig Vertices to Cut Mesh", mXTKModel->mVerboseLevel, 1 );

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
                                  aDecompositionData->tNewNodeId( iV ) )
                                  == aCutIntegrationMesh->mIntegrationVertexIdToIndexMap.end(),
                    "Id already in the map" );
            aCutIntegrationMesh->mIntegrationVertexIdToIndexMap[ aDecompositionData->tNewNodeId( iV ) ] = aDecompositionData->tNewNodeIndex( iV );

            // add the ancestry information to the mesh
            aCutIntegrationMesh->mIgVertexParentEntityRank( aDecompositionData->tNewNodeIndex( iV ) ) =
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

        // construct a relationship for the geometry engine to have geometric information for this vertex
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
            mGeometryEngine->create_new_derived_nodes(
                    aDecompositionData->mNewNodeParentCells,
                    aDecompositionData->mNewVertexLocalCoordWRTParentCell );
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Integration_Mesh_Generator::select_background_cell_for_edge(
            std::shared_ptr< Edge_Based_Connectivity > aEdgeBasedConnectivity,
            Cut_Integration_Mesh*                      aCutIntegrationMesh,
            moris::Cell< moris::mtk::Cell* >&          aBackgroundCellForEdge )
    {
        Tracer tTracer( "XTK", "Decomposition_Algorithm", "Select BG Cell for Edge", mXTKModel->mVerboseLevel, 1 );
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
        Tracer tTracer( "XTK", "Decomposition_Algorithm", "Select BG Cell for Facet", mXTKModel->mVerboseLevel, 1 );
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
        Tracer tTracer( "XTK", "Integration_Mesh_Generator", "Collect Vertex Groups for BG Cell", mXTKModel->mVerboseLevel, 1 );
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
        // log/trace this operation if increased verbose level (i.e. output detail) has been requested
        Tracer tTracer( "XTK", "Integration_Mesh_Generator", "Allocate child meshes", mXTKModel->mVerboseLevel, 1 );

        // allocate these data structures one per background cell
        aCutIntegrationMesh->mChildMeshes.resize( aBackgroundMesh->get_num_elems() );
        aCutIntegrationMesh->mIntegrationCellGroups.resize( aBackgroundMesh->get_num_elems(), nullptr );
        aCutIntegrationMesh->mIntegrationVertexGroups.resize( aBackgroundMesh->get_num_elems(), nullptr );
        aCutIntegrationMesh->mIntegrationCellGroupsParentCell.resize( aBackgroundMesh->get_num_elems(), nullptr );

        // create the child meshes
        for ( moris::uint iCell = 0; iCell < aBackgroundMesh->get_num_elems(); iCell++ )
        {
            // initialize the Child Mesh with its group of IG cells and populate their information
            moris_index       tCMIndex                                        = (moris_index)iCell;
            moris::mtk::Cell* tParentCell                                     = &aBackgroundMesh->get_mtk_cell( iCell );
            aCutIntegrationMesh->mChildMeshes( tCMIndex )                     = std::make_shared< Child_Mesh_Experimental >();
            aCutIntegrationMesh->mIntegrationCellGroups( tCMIndex )           = std::make_shared< IG_Cell_Group >( 0 );
            aCutIntegrationMesh->mChildMeshes( tCMIndex )->mIgCells           = aCutIntegrationMesh->mIntegrationCellGroups( tCMIndex );
            aCutIntegrationMesh->mIntegrationCellGroupsParentCell( tCMIndex ) = tParentCell;
            aCutIntegrationMesh->mChildMeshes( tCMIndex )->mParentCell        = tParentCell;
            aCutIntegrationMesh->mChildMeshes( tCMIndex )->mChildMeshIndex    = tCMIndex;

            // populate map linking the parent cell index to the child mesh index living on it
            aCutIntegrationMesh->mParentCellCellGroupIndex( tParentCell->get_index() ) = tCMIndex;

            // fixme: std::cout << "Integration_Mesh_Generator::allocate_child_meshes() - WARNING: GENERAlIZE NEEDED FOR MULTIPLE TOPOS" << std::endl;

            // get number of spatial dimensions and decide on cell topology of integration elements
            moris_index tNumGeometricVertices = determine_num_nodes( this->get_spatial_dim(), mtk::Interpolation_Order::LINEAR, mtk::CellShape::RECTANGULAR );

            // get list of vertices in current parent IP cell on background mesh
            moris::Cell< moris::mtk::Vertex* > tParentCellVerts = tParentCell->get_vertex_pointers();

            // get the parametric coordinates of the vertices of the parent cell
            Matrix< DDRMat > tParamCoords;
            tParentCell->get_cell_info()->get_loc_coords_of_cell( tParamCoords );

            // initialize and create a vertex group from the currently still un-cut background element
            aCutIntegrationMesh->mIntegrationVertexGroups( tCMIndex ) = std::make_shared< IG_Vertex_Group >( tNumGeometricVertices );
            aCutIntegrationMesh->mChildMeshes( tCMIndex )->mIgVerts   = aCutIntegrationMesh->mIntegrationVertexGroups( tCMIndex );
            // FIXME: GET GEOMETRIC VERTICES FROM MTK CELL HARDCODED TO HEX-FAMILY
            for ( moris::moris_index i = 0; i < tNumGeometricVertices; i++ )
            {
                aCutIntegrationMesh->mIntegrationVertexGroups( tCMIndex )->add_vertex( tParentCellVerts( i ), std::make_shared< Matrix< DDRMat > >( tParamCoords.get_row( i ) ) );
            }
        }

        // return successful child mesh creation
        return true;
    }

    // ----------------------------------------------------------------------------------

    void
    Integration_Mesh_Generator::assign_node_requests_identifiers(
            Decomposition_Data&   aDecompData,
            Cut_Integration_Mesh* aCutIntegrationMesh,
            moris::mtk::Mesh*     aBackgroundMesh )
    {

        moris_index tNodeIndex = aCutIntegrationMesh->get_first_available_index( mtk::EntityRank::NODE );

        for ( moris::uint i = 0; i < aDecompData.tNewNodeIndex.size(); i++ )
        {
            // set the new node index
            aDecompData.tNewNodeIndex( i ) = tNodeIndex;
            tNodeIndex++;
        }

        // perform sizing checks
        MORIS_ERROR( aDecompData.tNewNodeId.size() == aDecompData.tNewNodeIndex.size(),
                "Integration_Mesh_Generator::assign_node_requests_identifiers() - "
                "Dimension mismatch in assign_node_requests_identifiers" );
        MORIS_ERROR( aDecompData.tNewNodeId.size() == aDecompData.tNewNodeParentRank.size(),
                "Integration_Mesh_Generator::assign_node_requests_identifiers() - "
                "Dimension mismatch in assign_node_requests_identifiers" );
        MORIS_ERROR( aDecompData.tNewNodeId.size() == aDecompData.tNewNodeParentIndex.size(),
                "Integration_Mesh_Generator::assign_node_requests_identifiers() - "
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

        // convert the comm table from cell of procs to matrix format
        Matrix< IdMat > tCommTable( 1, tProcRanks.size() );
        for ( uint iProc = 0; iProc < tProcRanks.size(); iProc++ )
        {
            tCommTable( iProc ) = (moris_id)tProcRanks( iProc );
        }

        // allocate ids for nodes I own
        moris_id tNodeId = aCutIntegrationMesh->allocate_entity_ids( aDecompData.tNewNodeId.size(), mtk::EntityRank::NODE );

        // Assign owned request identifiers
        this->assign_owned_request_id( aDecompData, tOwnedRequest, tNodeId );

        // prepare node information request data
        Cell< Matrix< IndexMat > > tOutwardRequests;
        this->setup_outward_requests( aDecompData, aBackgroundMesh, tNotOwnedRequests, tProcRanks, tProcRankToDataIndex, tOutwardRequests );

        // send and receive the requests
        Cell< Matrix< IndexMat > > tReceivedRequests;
        moris::communicate_mats( tCommTable, tOutwardRequests, tReceivedRequests );

        // Prepare request answers
        Cell< Matrix< IndexMat > > tRequestAnswers;
        this->prepare_request_answers( aDecompData, aBackgroundMesh, tReceivedRequests, tRequestAnswers );

        // send and receive the answers
        Cell< Matrix< IndexMat > > tReceivedRequestsAnswers;
        moris::communicate_mats( tCommTable, tRequestAnswers, tReceivedRequestsAnswers );

        // handle received information
        this->handle_received_request_answers( aDecompData, aBackgroundMesh, tOutwardRequests, tReceivedRequestsAnswers, tNodeId );

        // check that all nodes have been assigned an ID
        MORIS_ERROR( mXTKModel->verify_successful_node_assignment( aDecompData ),
                "Integration_Mesh_Generator::assign_node_requests_identifiers() - "
                "Unsuccessful node assignment detected." );
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
            aProcRankToIndexInData[ tCommTable( i ) ] = i;
            aProcRanks( i )                           = ( tCommTable( i ) );

            // initialize and reserve memory for not-owned nodes
            aNotOwnedRequests.push_back( Cell< uint >( 0 ) );
            aNotOwnedRequests.back().reserve( 2 * tNumNewNodes / tCommTable.numel() );
        }

        // iterate through each node request and figure out the owner
        for ( moris::uint i = 0; i < tNumNewNodes; i++ )
        {
            // Parent Rank
            mtk::EntityRank    tParentRank  = tDecompData.tNewNodeParentRank( i );
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
                moris_index tIndex = aProcRankToIndexInData[ tOwnerProc ];

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
            Decomposition_Data&  aDecompData,
            Cell< uint > const & aOwnedRequest,
            moris::moris_id&     aNodeId )
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
            Decomposition_Data const &                aDecompData,
            moris::mtk::Mesh*                         aBackgroundMesh,
            Cell< Cell< uint > > const &              aNotOwnedRequests,
            Cell< uint > const &                      aProcRanks,
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
            uint tIndexInData = aProcRankToIndexInData[ tProcRank ];

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
                mtk::EntityRank tParentRank   = aDecompData.tNewNodeParentRank( tRequestIndex );

                // swap out for hmr if needed (hmr calls edges in 2d faces)
                if ( aBackgroundMesh->get_mesh_type() == mtk::MeshType::HMR )
                {
                    if ( aBackgroundMesh->get_spatial_dim() == 2 )
                    {
                        if ( tParentRank == mtk::EntityRank::EDGE )
                        {
                            tParentRank = mtk::EntityRank::FACE;
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
            Decomposition_Data&                aDecompData,
            moris::mtk::Mesh*                  aBackgroundMesh,
            Cell< Matrix< IndexMat > > const & aReceiveData,
            Cell< Matrix< IndexMat > >&        aRequestAnswers )
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
                    mtk::EntityRank tParentRank    = (mtk::EntityRank)aReceiveData( i )( 1, j );
                    moris_id        tSecondaryId   = aReceiveData( i )( 2, j );
                    moris_index     tParentInd     = aBackgroundMesh->get_loc_entity_ind_from_entity_glb_id( tParentId, tParentRank );
                    bool            tRequestExists = false;
                    moris_index     tRequestIndex  = MORIS_INDEX_MAX;

                    // swap out for hmr if needed (hmr calls edges in 2d faces)
                    if ( aBackgroundMesh->get_mesh_type() == mtk::MeshType::HMR )
                    {
                        if ( aBackgroundMesh->get_spatial_dim() == 2 )
                        {
                            if ( tParentRank == mtk::EntityRank::FACE )
                            {
                                tParentRank = mtk::EntityRank::EDGE;
                            }
                        }
                    }

                    if ( aDecompData.mHasSecondaryIdentifier )
                    {
                        tRequestExists = aDecompData.request_exists(
                                tParentInd,
                                tSecondaryId,
                                (mtk::EntityRank)tParentRank,
                                tRequestIndex );
                    }
                    else
                    {
                        tRequestExists = aDecompData.request_exists(
                                tParentInd,
                                (mtk::EntityRank)tParentRank,
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
            Decomposition_Data&                aDecompData,
            moris::mtk::Mesh*                  aBackgroundMesh,
            Cell< Matrix< IndexMat > > const & aRequests,
            Cell< Matrix< IndexMat > > const & aRequestAnswers,
            moris::moris_id&                   aNodeId )
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
                    mtk::EntityRank tParentRank    = (mtk::EntityRank)aRequests( i )( 1, j );
                    moris_id        tSecondaryId   = aRequests( i )( 2, j );
                    moris_index     tParentInd     = aBackgroundMesh->get_loc_entity_ind_from_entity_glb_id( tParentId, tParentRank );
                    bool            tRequestExists = false;
                    moris_index     tRequestIndex  = MORIS_INDEX_MAX;


                    // swap out for hmr if needed (hmr calls edges in 2d faces)
                    if ( aBackgroundMesh->get_mesh_type() == mtk::MeshType::HMR )
                    {
                        if ( aBackgroundMesh->get_spatial_dim() == 2 )
                        {
                            if ( tParentRank == mtk::EntityRank::FACE )
                            {
                                tParentRank = mtk::EntityRank::EDGE;
                            }
                        }
                    }


                    if ( aDecompData.mHasSecondaryIdentifier )
                    {
                        tRequestExists = aDecompData.request_exists( tParentInd, tSecondaryId, (mtk::EntityRank)tParentRank, tRequestIndex );
                    }
                    else
                    {
                        tRequestExists = aDecompData.request_exists( tParentInd, (mtk::EntityRank)tParentRank, tRequestIndex );
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

                            // mBackgroundMesh.update_first_available_index(tNodeIndex, mtk::EntityRank::NODE);
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
    Integration_Mesh_Generator::remove_subphases_from_cut_mesh( moris::Cell< moris_index > const & aSubphasesToRemove )
    {
    }


    // ----------------------------------------------------------------------------------
    // Functions for Constructing Subphase Groups and their neighborhood
    // ----------------------------------------------------------------------------------

    bool
    Integration_Mesh_Generator::check_construct_subphase_groups()
    {
        // initialize output
        bool tConstructSPGs = false;

        // only construct Subphase-Groups and their Neighborhood if SPG based enrichment has been requested
        // first check if XTK-Model has parameter-list (i.e. skip this in Unit tests that don't require it)
        if ( mXTKModel->mParameterList.get< bool >( "has_parameter_list" ) )
        {
            if ( mXTKModel->mParameterList.get< bool >( "use_SPG_based_enrichment" ) )
            {
                // check that the mesh indices for enrichment are defined
                MORIS_ERROR( !mXTKModel->mParameterList.get< std::string >( "enrich_mesh_indices" ).empty(),
                        "Integration_Mesh_Generator::check_construct_subphase_groups() - No B-spline mesh indices provided for enrichment. Unable to construct Subphase-groups." );

                // if all checks have past, set flag for constructing SPGs
                tConstructSPGs = true;
            }
        }

        // return flag
        return tConstructSPGs;
    }

    // ----------------------------------------------------------------------------------

    void
    Integration_Mesh_Generator::establish_bspline_mesh_info(
            Cut_Integration_Mesh* aCutIntegrationMesh,
            moris::mtk::Mesh*     aLagrangeMesh,
            Bspline_Mesh_Info*    aBsplineMeshInfo,
            const moris_index     aMeshIndex )
    {
        // log/trace the enrichment for specific B-spline mesh
        Tracer tTracer( "XTK", "Integration_Mesh_Generator", "Establish B-spline mesh info for mesh index " + std::to_string( aMeshIndex ) );

        // get from HMR which Lagrange elements sit in which B-spline elements
        aLagrangeMesh->get_lagrange_elements_in_bspline_elements(
                aMeshIndex,
                aBsplineMeshInfo->mExtractionCellsInBsplineCells,
                aBsplineMeshInfo->mExtractionCellsIndicesInBsplineCells,
                aBsplineMeshInfo->mExtractionCellToBsplineCell,
                aBsplineMeshInfo->mBsplineCellLevels,
                aBsplineMeshInfo->mBsplineCells );
    }

    // ----------------------------------------------------------------------------------

    void
    Integration_Mesh_Generator::construct_subphase_groups(
            Cut_Integration_Mesh* aCutIntegrationMesh,
            Bspline_Mesh_Info*    aBsplineMeshInfo,
            const moris_index     aMeshListIndex )
    {
        // log/trace the enrichment for specific B-spline mesh
        Tracer tTracer( "XTK", "Integration_Mesh_Generator", "Construct SPGs for mesh index " + std::to_string( aMeshListIndex ) );

        // get the number of active B-spline elements
        uint tNumBspElems = aBsplineMeshInfo->mExtractionCellsIndicesInBsplineCells.size();

        // reserve memory for list of SPGs in Bspline mesh info using estimated size
        aBsplineMeshInfo->mSubphaseGroups.reserve( 4 * tNumBspElems );

        // initialize SPG in B-spline element map with correct size
        aBsplineMeshInfo->mSpgIndicesInBsplineCells.resize( tNumBspElems );

        // establish SPGs on every B-spline element
        for ( uint iBspElem = 0; iBspElem < tNumBspElems; iBspElem++ )
        {
            // initialize and get list of subphase indices present on current B-spline element
            Matrix< IndexMat > tSubphaseIndicesInBsplineCell;
            this->get_subphase_indices_in_bspline_cell( aCutIntegrationMesh, aMeshListIndex, iBspElem, tSubphaseIndicesInBsplineCell );

            // fill searchable map with list of subphase indices obtained
            IndexMap tSubphaseIndexToBsplineCell;
            this->construct_subphase_in_bspline_cell_map( tSubphaseIndicesInBsplineCell, tSubphaseIndexToBsplineCell );

            // cut down the subphase neighborhood to only the B-spline cells
            moris::Matrix< moris::IndexMat > tPrunedSubphaseToSubphase;
            this->generate_pruned_subphase_graph_in_bspline_cell(
                    aCutIntegrationMesh,
                    tSubphaseIndicesInBsplineCell,
                    tSubphaseIndexToBsplineCell,
                    tPrunedSubphaseToSubphase );

            // get the number of SPGs in extraction cell cluster by performing a flood fill
            moris_index                      tMaxSpgInd = 0;
            moris::Matrix< moris::IndexMat > tSubphaseBin;
            this->find_subphase_bin_enrichment_levels_in_bspline_cell(
                    tSubphaseIndicesInBsplineCell,
                    tPrunedSubphaseToSubphase,
                    tSubphaseBin,
                    tMaxSpgInd );

            // increment to get actual number of Subphase groups (instead of just index)
            uint tNumSPGs = (uint)tMaxSpgInd + 1;

            // split subphase bin up into SPGs
            moris::Cell< moris::Cell< moris_index > > tSPsInBin = this->split_flood_fill_bin( tSubphaseBin, tSubphaseIndicesInBsplineCell, tNumSPGs );

            // debug
            MORIS_ASSERT( tNumSPGs == tSPsInBin.size(), "Integration_Mesh_Generator::create_subphase_groups() - Something doesn't line up..." );

            // initialize memory for the list of SPGs on the current B-spline element with correct size
            aBsplineMeshInfo->mSpgIndicesInBsplineCells( iBspElem ).reserve( tNumSPGs );

            // for each disconnected set of subphases ...
            for ( moris::size_t iSPG = 0; iSPG < tNumSPGs; iSPG++ )
            {
                // create SPGs and add to mesh
                aBsplineMeshInfo->add_subphase_group_to_bspline_cell( tSPsInBin( iSPG ), (moris_index)iBspElem );

                // create a list of IG cells in the SPG
                moris::Cell< moris_index > tIgCellIndicesInSPG;

                // go through SPs and collect the IG cells in them in one list
                this->collect_ig_cell_indices_in_SPG( aCutIntegrationMesh, tSPsInBin( iSPG ), tIgCellIndicesInSPG );

                // set the list of IG cells to the b-spline mesh info
                aBsplineMeshInfo->add_ig_cell_indices_to_last_admitted_subphase_group( tIgCellIndicesInSPG );

                // get bulk phase of the subphase group and set it
                moris_index tBulkPhaseIndex = aCutIntegrationMesh->get_subphase_bulk_phase( tSPsInBin( iSPG )( 0 ) );
                aBsplineMeshInfo->set_bulk_phase_of_last_admitted_subphase_group( tBulkPhaseIndex );
            }
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Integration_Mesh_Generator::communicate_subphase_groups(
            Cut_Integration_Mesh* aCutIntegrationMesh,
            Bspline_Mesh_Info*    aBsplineMeshInfo,
            const moris_index     aBsplineMeshListIndex )
    {
        // get the number of Subphase groups on the current B-spline mesh
        uint tNumSPGs     = aBsplineMeshInfo->get_num_SPGs();
        uint tHalfNumSPGs = (uint)std::floor( (real)tNumSPGs / 2.0 );

        // reserve memory for lists of owned and non-owned indices
        aBsplineMeshInfo->mOwnedSubphaseGroupIndices.reserve( tNumSPGs );
        aBsplineMeshInfo->mNotOwnedSubphaseGroupIndices.reserve( tHalfNumSPGs );

        // get current proc rank
        moris_index tCurrentParRank = moris::par_rank();

        // go through sub-phases groups and sort them into ownership groups based on their parent cells
        for ( moris::uint iSPG = 0; iSPG < tNumSPGs; iSPG++ )
        {
            // NOTE: the following operation assumes that the Lagrange mesh is at least as refined as the B-spline mesh and ...
            // NOTE: ... that any Lagrange element and its underlying B-spline element are on the same proc

            // get a representative IP cell for the Subphase group
            moris_index              tBsplineCellIndex = aBsplineMeshInfo->mSubphaseGroups( iSPG )->get_bspline_cell_index();
            moris::mtk::Cell const * tIpCell           = aBsplineMeshInfo->mExtractionCellsInBsplineCells( tBsplineCellIndex )( 0 );

            // get the owning proc from the IP cell
            moris_index tOwner = tIpCell->get_owner();

            // set the owner of the SPG for quick retrieval later
            aBsplineMeshInfo->mSubphaseGroups( iSPG )->set_owning_proc( tOwner );

            // depending on whether the current proc is the owner or not, list the SPG as owned or non-owned
            if ( tOwner == tCurrentParRank )
            {
                aBsplineMeshInfo->mOwnedSubphaseGroupIndices.push_back( (moris_index)iSPG );
            }
            else
            {
                aBsplineMeshInfo->mNotOwnedSubphaseGroupIndices.push_back( (moris_index)iSPG );
            }
        }

        // shrink to fit subphase group ownership lists (in case size was over-estimated on initialization)
        aBsplineMeshInfo->mOwnedSubphaseGroupIndices.shrink_to_fit();
        aBsplineMeshInfo->mNotOwnedSubphaseGroupIndices.shrink_to_fit();

        // give all sub-phases a global ID (across all procs)
        this->assign_subphase_group_glob_ids( aCutIntegrationMesh, aBsplineMeshInfo );

        // construct reverse map from SPG IDs to indices on each proc
        aCutIntegrationMesh->construct_spg_id_to_index_map( aBsplineMeshListIndex );
    }

    // ----------------------------------------------------------------------------------

    void
    Integration_Mesh_Generator::get_subphase_indices_in_bspline_cell(
            Cut_Integration_Mesh* aCutIntegrationMesh,
            moris_index           aMeshListIndex,
            uint                  aCurrentBspCellIndex,
            Matrix< IndexMat >&   aSubPhaseIndices )
    {
        // count the number of subphase clusters in support
        moris::uint tCount = 0;

        // store IP cell indices in current B-spline cell
        moris::Cell< moris_index > tLagElemInds =
                aCutIntegrationMesh->mBsplineMeshInfos( aMeshListIndex )->mExtractionCellsIndicesInBsplineCells( aCurrentBspCellIndex );

        // get number of Lagrange elements in current B-spline element
        moris::size_t tNumLagElems = tLagElemInds.size();

        // loop over Lagrange elements in bspline cell to count subphases
        for ( moris::size_t iE = 0; iE < tNumLagElems; iE++ )
        {
            // get the subphase (indices) living in current element
            moris::Cell< moris_index > const & tSubphaseIndicesOnLagElem = aCutIntegrationMesh->get_parent_cell_subphases( tLagElemInds( iE ) );

            // count number of subphases in support
            tCount = tCount + tSubphaseIndicesOnLagElem.size();
        }

        // initialize list of subphase indices in bspline cell
        aSubPhaseIndices.resize( 1, tCount );
        tCount = 0;

        // loop over Lagrange elements in bspline cell to store subphase indices
        for ( moris::size_t iE = 0; iE < tNumLagElems; iE++ )
        {
            // get the subphase (indices) living in current Lagrange element
            moris::Cell< moris_index > const & tSubphaseIndicesOnLagElem = aCutIntegrationMesh->get_parent_cell_subphases( tLagElemInds( iE ) );

            // go over subphases in current element and ...
            for ( moris::uint iSP = 0; iSP < tSubphaseIndicesOnLagElem.size(); iSP++ )
            {
                // ... fill the list of subphase indices in basis support
                aSubPhaseIndices( tCount++ ) = tSubphaseIndicesOnLagElem( iSP );
            }
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Integration_Mesh_Generator::construct_subphase_in_bspline_cell_map(
            moris::Matrix< moris::IndexMat > const & aSubphaseIndicesInBsplineCell,
            IndexMap&                                aSubphaseIndexToBsplineCellIndex )
    {
        // loop over subphases (indices) in B-spline cell
        for ( moris::moris_index i = 0; i < (moris::moris_index)aSubphaseIndicesInBsplineCell.numel(); i++ )
        {
            aSubphaseIndexToBsplineCellIndex[ aSubphaseIndicesInBsplineCell( i ) ] = i;
        }
    }

    //-------------------------------------------------------------------------------------

    void
    Integration_Mesh_Generator::generate_pruned_subphase_graph_in_bspline_cell(
            Cut_Integration_Mesh*                    aCutIntegrationMesh,
            moris::Matrix< moris::IndexMat > const & aSubphasesInBsplineCell,
            IndexMap&                                aSubphaseIndicesToBspline,
            moris::Matrix< moris::IndexMat >&        aPrunedBsplineSubphaseToSubphase )
    {
        // get pointer to Subphase Neighborhood Connectivity of mesh
        std::shared_ptr< Subphase_Neighborhood_Connectivity > tSubphaseNeighborhood = aCutIntegrationMesh->get_subphase_neighborhood();

        // get the number of SPs in the current B-spline cell
        uint tNumSpsInBspCell = aSubphasesInBsplineCell.numel();

        // get subphase neighborhood information form Subphase_Neighborhood_Connectivity and store in variable for easy handling
        moris::Cell< std::shared_ptr< moris::Cell< moris_index > > > const & tSubphaseToSubphase = tSubphaseNeighborhood->mSubphaseToSubPhase;

        // determine the maximum number of neighbors needed in graph
        uint tMaxNumNeighbors = 0;
        for ( moris::size_t iSP = 0; iSP < tNumSpsInBspCell; iSP++ )
        {
            moris::Cell< moris_index > const & tSingleSubphaseNeighbors = *tSubphaseToSubphase( aSubphasesInBsplineCell( iSP ) );
            uint                               tNumNeighbors            = tSingleSubphaseNeighbors.size();
            tMaxNumNeighbors                                            = std::max( tMaxNumNeighbors, tNumNeighbors );
        }

        // initialize matrix storing pruned subphase connectivity
        aPrunedBsplineSubphaseToSubphase.resize( tNumSpsInBspCell, tMaxNumNeighbors );
        aPrunedBsplineSubphaseToSubphase.fill( MORIS_INDEX_MAX );

        // go over subphases in B-spline Cell
        for ( moris::size_t iSP = 0; iSP < tNumSpsInBspCell; iSP++ )
        {
            // get list of current subphase's neighbors
            moris::Cell< moris_index > const & tSingleSubphaseNeighbors = *tSubphaseToSubphase( aSubphasesInBsplineCell( iSP ) );

            // iterate through neighbors and check if neighbors and prune if not in B-spline cell
            moris::uint tCount = 0;
            for ( moris::size_t jSpNeighbor = 0; jSpNeighbor < tSingleSubphaseNeighbors.size(); jSpNeighbor++ )
            {
                // get current neighbor subphase's index
                moris_index tNeighborSubphaseIndex = tSingleSubphaseNeighbors( jSpNeighbor );

                // find subphase index in list of subphases inside B-spline cell, and ...
                auto tNeighborIter = aSubphaseIndicesToBspline.find( tNeighborSubphaseIndex );

                // ... only add neighboring subphase, if it is found
                if ( tNeighborIter != aSubphaseIndicesToBspline.end() )
                {
                    aPrunedBsplineSubphaseToSubphase( iSP, tCount ) = tNeighborIter->second;
                    tCount++;
                }
            }
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Integration_Mesh_Generator::find_subphase_bin_enrichment_levels_in_bspline_cell(
            moris::Matrix< moris::IndexMat > const & aSubphasesInBspCell,
            moris::Matrix< moris::IndexMat > const & aPrunedSubPhaseToSubphase,
            moris::Matrix< moris::IndexMat >&        aSubPhaseBinEnrichmentVals,
            moris_index&                             aMaxEnrichmentLevel )
    {
        // Variables needed for flood-fill, consider removing these.
        // Active bins to include in flood-fill (We include all bins)
        moris::Matrix< moris::IndexMat > tActiveBins( 1, aPrunedSubPhaseToSubphase.n_rows() );

        for ( moris::size_t i = 0; i < aPrunedSubPhaseToSubphase.n_rows(); i++ )
        {
            ( tActiveBins )( 0, i ) = i;
        }

        // Mark all as included
        moris::Matrix< moris::IndexMat > tIncludedBins( 1, aSubphasesInBspCell.numel(), 1 );

        // Flood fill metric value (since all the subphases do not connect to dissimilar phases)
        moris::Matrix< moris::IndexMat > tDummyPhase( 1, aSubphasesInBspCell.numel(), 1 );

        aSubPhaseBinEnrichmentVals = flood_fill(
                aPrunedSubPhaseToSubphase,
                tDummyPhase,
                tActiveBins,
                tIncludedBins,
                mGeometryEngine->get_num_phases(),
                MORIS_INDEX_MAX,
                aMaxEnrichmentLevel,
                true );
    }

    //-------------------------------------------------------------------------------------

    moris::Cell< moris::Cell< moris_index > >
    Integration_Mesh_Generator::split_flood_fill_bin(
            moris::Matrix< moris::IndexMat >& aSubphaseBin,
            moris::Matrix< moris::IndexMat >& aSubphaseIndicesInBsplineCell,
            uint                              aNumSPGs )
    {
        // get number of subphases
        moris::size_t tNumSubphases = aSubphaseIndicesInBsplineCell.length();

        // initialize array of counters tracking how many subphases are in each group
        moris::Cell< uint > tSpCounters( aNumSPGs, 0 );

        // loop over all subphases in bin and ...
        for ( moris::size_t iSP = 0; iSP < tNumSubphases; iSP++ )
        {
            // ... increment counter for SPG current SP belongs to
            uint tSpgIndex = aSubphaseBin( iSP );
            tSpCounters( tSpgIndex ) += 1;
        }

        // initialize output Cell
        moris::Cell< moris::Cell< moris_index > > tSPGsInBin( aNumSPGs );
        for ( moris::size_t iSPG = 0; iSPG < aNumSPGs; iSPG++ )
        {
            tSPGsInBin( iSPG ) = moris::Cell< moris_index >( tSpCounters( iSPG ) );
        }

        // loop over all subphases in bin and ... // ... put their indices into groups
        tSpCounters = moris::Cell< uint >( aNumSPGs, 0 );
        for ( moris::size_t iSP = 0; iSP < tNumSubphases; iSP++ )
        {
            // ... put their indices into groups
            uint tSpgIndex                                  = aSubphaseBin( iSP );
            uint tCurrentSpIndexInSPG                       = tSpCounters( tSpgIndex );
            tSPGsInBin( tSpgIndex )( tCurrentSpIndexInSPG ) = aSubphaseIndicesInBsplineCell( iSP );

            // increment counter for SPG current SP belongs to
            tSpCounters( tSpgIndex ) += 1;
        }

        // return sorted bins of subphases
        return tSPGsInBin;
    }

    // ----------------------------------------------------------------------------------

    void
    Integration_Mesh_Generator::collect_ig_cell_indices_in_SPG(
            Cut_Integration_Mesh*              aCutIntegrationMesh,
            moris::Cell< moris_index > const & aSPsInSPG,
            moris::Cell< moris_index >&        aIgCellIndicesInSPG )
    {
        // count up all IG cells within SPG
        uint tIgCellCounter = 0;

        // for each SP get the number of Ig cells and add to total
        for ( moris::size_t iSP = 0; iSP < aSPsInSPG.size(); iSP++ )
        {
            moris_index tSpIndex = aSPsInSPG( iSP );
            tIgCellCounter += aCutIntegrationMesh->get_subphase_ig_cells( tSpIndex )->mIgCellGroup.size();
        }

        // use number of IG cells to initialize list
        aIgCellIndicesInSPG.resize( tIgCellCounter );

        // reset counter
        tIgCellCounter = 0;

        // fill the list of IG cells in SPG by going over each of the individual subphases and their IG cells
        for ( moris::size_t iSP = 0; iSP < aSPsInSPG.size(); iSP++ )
        {
            // get the index of the current subphase
            moris_index tSpIndex = aSPsInSPG( iSP );

            // get the group of IG cells associated with current subphase
            const moris::Cell< moris::mtk::Cell* >* tIgCellGroup = &( aCutIntegrationMesh->get_subphase_ig_cells( tSpIndex )->mIgCellGroup );

            // loop over IG cells in SP and add to list of IG cells in SPG
            for ( moris::size_t iIgCell = 0; iIgCell < tIgCellGroup->size(); iIgCell++ )
            {
                aIgCellIndicesInSPG( tIgCellCounter ) = ( *tIgCellGroup )( iIgCell )->get_index();
                tIgCellCounter++;
            }
        }
    }

    // ----------------------------------------------------------------------------------

    moris::Cell< bool >
    Integration_Mesh_Generator::collect_subphase_group_ligament_side_ordinals(
            Cut_Integration_Mesh*       aCutIntegrationMesh,
            moris::Cell< moris_index >& aSPsInGroup,
            IndexMap&                   aSubphaseIndicesToBspline )
    {
        // initialize punch-card for which side ordinals are used
        moris::Cell< bool > tUsedSideOrdinals;
        if ( mXTKModel->get_spatial_dim() == 2 )
        {
            tUsedSideOrdinals.resize( 4, false );
        }
        else if ( mXTKModel->get_spatial_dim() == 3 )
        {
            tUsedSideOrdinals.resize( 6, false );
        }
        else
        {
            MORIS_ERROR( false, "Integration_Mesh_Generator::collect_subphase_group_ligament_side_ordinals() - Number of spatial dims must be 2 or 3" );
        }

        // get pointer to Subphase Neighborhood Connectivity of mesh
        const std::shared_ptr< Subphase_Neighborhood_Connectivity > tSubphaseNeighborhood = aCutIntegrationMesh->get_subphase_neighborhood();

        // get subphase neighborhood information form Subphase_Neighborhood_Connectivity and store in variable for easy handling
        moris::Cell< std::shared_ptr< moris::Cell< moris_index > > > const & tSubphaseToSubphase   = tSubphaseNeighborhood->mSubphaseToSubPhase;
        moris::Cell< std::shared_ptr< moris::Cell< moris_index > > > const & tSubphaseSideOrdinals = tSubphaseNeighborhood->mSubphaseToSubPhaseMySideOrds;

        // go over subphases in subphase group
        for ( moris::size_t iSP = 0; iSP < aSPsInGroup.size(); iSP++ )
        {
            moris_index tCurrentSubphaseIndex = aSPsInGroup( iSP );

            // get list of current subphase's neighbors
            moris::Cell< moris_index > const & tSingleSubphaseNeighbors            = *tSubphaseToSubphase( tCurrentSubphaseIndex );
            moris::Cell< moris_index > const & tSingleSubphaseNeighborSideOrdinals = *tSubphaseSideOrdinals( tCurrentSubphaseIndex );

            // iterate through neighbors and check if neighbors outside of B-spline cell
            for ( moris::size_t iNeighbor = 0; iNeighbor < tSingleSubphaseNeighbors.size(); iNeighbor++ )
            {
                // get current neighbor subphase's index
                moris_index tNeighborSubphaseIndex        = tSingleSubphaseNeighbors( iNeighbor );
                moris_index tNeighborSubphaseSideOrdinals = tSingleSubphaseNeighborSideOrdinals( iNeighbor );

                // find subphase index in list of subphases inside B-spline cell, and ...
                auto tNeighborIter = aSubphaseIndicesToBspline.find( tNeighborSubphaseIndex );

                // ... add ligament side ordinal, if neighbor SP is outside B-spline Cell (i.e. it is NOT found in the list of SP indices in BSp elem)
                if ( tNeighborIter == aSubphaseIndicesToBspline.end() )
                {
                    tUsedSideOrdinals( tNeighborSubphaseSideOrdinals ) = true;
                }
            }
        }

        // return punch-card of used side ordinals
        return tUsedSideOrdinals;
    }

    // ----------------------------------------------------------------------------------

    const std::shared_ptr< Subphase_Neighborhood_Connectivity >
    Integration_Mesh_Generator::construct_subphase_group_neighborhood(
            Cut_Integration_Mesh* aCutIntegrationMesh,
            Bspline_Mesh_Info*    aBsplineMeshInfo,
            const moris_index     aMeshIndex )
    {
        // log/trace the enrichment for specific B-spline mesh
        Tracer tTracer( "XTK", "Integration_Mesh_Generator", "Construct SPG connectivity for mesh index " + std::to_string( aMeshIndex ) );

        // get the number of SPGs on current mesh index
        uint tNumSPGs = aBsplineMeshInfo->get_num_SPGs();

        // initialize SPG neighborhood
        std::shared_ptr< Subphase_Neighborhood_Connectivity > tSpgNeighborhood = std::make_shared< Subphase_Neighborhood_Connectivity >();

        // initialize size of SPG connectivity graph
        tSpgNeighborhood->mSubphaseToSubPhase.resize( tNumSPGs );
        tSpgNeighborhood->mSubphaseToSubPhaseMySideOrds.resize( tNumSPGs );
        tSpgNeighborhood->mSubphaseToSubPhaseNeighborSideOrds.resize( tNumSPGs );
        tSpgNeighborhood->mTransitionNeighborCellLocation.resize( tNumSPGs );

        // get the subphase neighborhood connectivity
        const std::shared_ptr< Subphase_Neighborhood_Connectivity > tSpNeighborhood = aCutIntegrationMesh->get_subphase_neighborhood();

        // loop over subphase groups on current mesh index
        for ( moris::size_t iSPG = 0; iSPG < tNumSPGs; iSPG++ )
        {
            // initialize sub-lists in SPG neighborhood ( for each SPG index there's a list )
            tSpgNeighborhood->mSubphaseToSubPhase( iSPG )                 = std::make_shared< moris::Cell< moris_index > >();
            tSpgNeighborhood->mSubphaseToSubPhaseMySideOrds( iSPG )       = std::make_shared< moris::Cell< moris_index > >();
            tSpgNeighborhood->mSubphaseToSubPhaseNeighborSideOrds( iSPG ) = std::make_shared< moris::Cell< moris_index > >();

            // leave the transition location empty since this is not needed
            tSpgNeighborhood->mTransitionNeighborCellLocation( iSPG ) = std::make_shared< moris::Cell< moris_index > >( 0 );

            // reserve memory for sub-lists in SPG neighborhood according to estimate
            tSpgNeighborhood->mSubphaseToSubPhase( iSPG )->reserve( 4 );
            tSpgNeighborhood->mSubphaseToSubPhaseMySideOrds( iSPG )->reserve( 4 );
            tSpgNeighborhood->mSubphaseToSubPhaseNeighborSideOrds( iSPG )->reserve( 4 );

            // get the SP indices on the current SPG
            const moris::Cell< moris_index >& tSpIndicesInGroup = aBsplineMeshInfo->mSubphaseGroups( iSPG )->get_SP_indices_in_group();

            // get the SPG's index
            const moris_index tSpgIndex = aBsplineMeshInfo->mSubphaseGroups( iSPG )->get_index();

            // initialize map storing which neighbor SPGs have already been found to avoid duplicate ligaments
            std::unordered_set< moris_index > tNeighborSPGsFound;

            // loop over the SPs in current SPG
            for ( moris::size_t iSP = 0; iSP < tSpIndicesInGroup.size(); iSP++ )
            {
                // get current SP's index
                const uint tCurrentSpIndex = (uint)tSpIndicesInGroup( iSP );

                // get the neighbors of the current subphase
                const std::shared_ptr< moris::Cell< moris_index > > tSpNeighborSPs = tSpNeighborhood->mSubphaseToSubPhase( tCurrentSpIndex );

                // loop over the neighbors and check if they're in separate SPGs
                for ( moris::size_t iSpNeighbor = 0; iSpNeighbor < tSpNeighborSPs->size(); iSpNeighbor++ )
                {
                    // get Neighbor SP's index
                    const moris_index tNeighborSpIndex = ( *tSpNeighborSPs )( iSpNeighbor );

                    // get current neighbor SP's SPG index
                    const moris_index tNeighborSpSpgIndex = aBsplineMeshInfo->mSpToSpgMap( tNeighborSpIndex );

                    // check whether the two SPs are in the same SPG; if not, establish connection in SPG connectivity
                    if ( tNeighborSpSpgIndex != tSpgIndex )
                    {
                        // check whether the connection to the neighboring SPG has already been found
                        if ( tNeighborSPGsFound.find( tNeighborSpSpgIndex ) == tNeighborSPGsFound.end() )
                        {
                            // register that the neighbor SPG has been found
                            tNeighborSPGsFound.insert( tNeighborSpSpgIndex );

                            // save neighbor SPG index in SPG neighborhood connectivity
                            tSpgNeighborhood->mSubphaseToSubPhase( iSPG )->push_back( tNeighborSpSpgIndex );

                            // get connectivity information for SPGs from SP neighborhood connectivity
                            const moris_index tMySideOrdinal    = ( *tSpNeighborhood->mSubphaseToSubPhaseMySideOrds( tCurrentSpIndex ) )( iSpNeighbor );
                            const moris_index tOtherSideOrdinal = ( *tSpNeighborhood->mSubphaseToSubPhaseNeighborSideOrds( tCurrentSpIndex ) )( iSpNeighbor );

                            // put connectivity information into SPG neighborhood connectivity
                            tSpgNeighborhood->mSubphaseToSubPhaseMySideOrds( iSPG )->push_back( tMySideOrdinal );
                            tSpgNeighborhood->mSubphaseToSubPhaseNeighborSideOrds( iSPG )->push_back( tOtherSideOrdinal );
                        }
                    }
                }    // end: loop over neighbor SPs
            }        // end: loop over SPs in SPG
        }            // end: loop over SPGs

        // hand SPG neighborhood connectivity over to cut integration mesh
        // aCutIntegrationMesh->mSubphaseGroupNeighborhood( iBspMesh ) = tSpNeighborhood;
        return tSpgNeighborhood;
    }

    // ----------------------------------------------------------------------------------

    void
    Integration_Mesh_Generator::construct_SPG_material_connectivity_information( Cut_Integration_Mesh* aCutIntegrationMesh )
    {
        // log/trace the enrichment for specific B-spline mesh
        Tracer tTracer( "XTK", "Integration_Mesh_Generator", "Construct SPG material connectivity information" );

        // get the number of base IP cells on the mesh
        uint tNumBaseIpCells = aCutIntegrationMesh->mBsplineMeshInfos( 0 )->mExtractionCellToBsplineCell.size();

        // debug - temporary test to make sure the right function is called
        MORIS_ASSERT( tNumBaseIpCells == aCutIntegrationMesh->get_num_base_ip_cells(),
                "Integration_Mesh_Generator::construct_SPG_material_connectivity_information() - "
                "Number of Lagrange elements reported by the xtk background mesh and the B-spline mesh info do not match" );

        // initialize list needed later
        aCutIntegrationMesh->mUnionVoidMsdIndices.resize( tNumBaseIpCells );
        aCutIntegrationMesh->mUnionVoidMsdIndexBulkPhases.resize( tNumBaseIpCells );

        // number of B-spline meshes
        Matrix< IndexMat > tBspMeshIndices = mXTKModel->get_Bspline_mesh_indices();
        uint               tNumBspMeshes   = tBspMeshIndices.numel();

        // initialize storage for generated data
        for ( uint iBspMesh = 0; iBspMesh < tNumBspMeshes; iBspMesh++ )
        {
            // get the current B-spline mesh info
            Bspline_Mesh_Info* tBsplineMeshInfo = aCutIntegrationMesh->mBsplineMeshInfos( iBspMesh );

            // initialize lists holding which SPGs are material and void on the base IP cells/Lagrange elements
            tBsplineMeshInfo->mExtractionCellMaterialSpgs.resize( tNumBaseIpCells );
            tBsplineMeshInfo->mExtractionCellVoidSpgs.resize( tNumBaseIpCells );
            tBsplineMeshInfo->mExtractionCellMaterialMsdIndices.resize( tNumBaseIpCells );
            tBsplineMeshInfo->mExtractionCellVoidMsdIndices.resize( tNumBaseIpCells );
            tBsplineMeshInfo->mExtractionCellFreeVoidMsdIndices.resize( tNumBaseIpCells );

            // initialize map assigning MSD indices to the SPGs
            tBsplineMeshInfo->mSpgToMsdIndex.resize( tBsplineMeshInfo->get_num_SPGs(), -1 );
        }

        // initialize the list storing which B-spline mesh index is the coarsest wrt. to a given Lagrange element
        Cell< moris_index >& tCoarsestBsplineMesh = aCutIntegrationMesh->mCoarsestBsplineMesh;
        tCoarsestBsplineMesh.resize( tNumBaseIpCells );

        // initialize list ticking off for which Lagrange elements the material indices have been constructed
        // (since there's a loop of Lag elems in the loop of Lag elems and we want to avoid unnecessarily treating elements twice)
        Cell< bool > tLagElemHasNotBeenTreated( tNumBaseIpCells, true );

        // all subsequent operations are per Lagrange element (= base IP cell)
        for ( uint iLagElem = 0; iLagElem < tNumBaseIpCells; iLagElem++ )
        {
            // --------------------------------
            // STEP 1: find the coarsest B-spline elements associated with each Lagrange element

            // initialize the coarsest level
            uint tCoarsestRefineLevel = MORIS_UINT_MAX;

            // find the mesh index of the B-spline element that contains a given Lagrange element
            for ( uint iBspMesh = 0; iBspMesh < tNumBspMeshes; iBspMesh++ )
            {
                // get the current B-spline mesh's discretization mesh index
                // moris_index tBspMeshIndex = tBspMeshIndices( iBspMesh );

                // get the B-spline mesh info for access
                Bspline_Mesh_Info* tBsplineMeshInfo = aCutIntegrationMesh->mBsplineMeshInfos( iBspMesh );

                // get the B-spline element index corresponding to the current Lagrange element
                uint tRefineLevel = tBsplineMeshInfo->get_B_spline_refinement_level_for_extraction_cell( iLagElem );

                // check whether the current B-spline mesh has a coarser B-spline element associated with
                // the current Lagrange element than any of the previous meshes
                if ( tRefineLevel < tCoarsestRefineLevel )
                {
                    // if so, safe the mesh as being the coarsest related mesh
                    tCoarsestRefineLevel             = tRefineLevel;
                    tCoarsestBsplineMesh( iLagElem ) = iBspMesh;
                }
            }

            // --------------------------------
            // STEP 2: relate the SPs on every Lagrange element to a material sub-domain (MSD) index relative to the coarsest mesh
            // Note: The below procedures are only valid as HMR imposes so-called "strong conditions"
            // Note: on the nested hierarchical domain boundaries, i.e. refined domain boundaries
            // Note: run along mesh lines that already exist on coarser levels

            // this step is done in batches of Lagrange elements belonging to the same coarsest B-spline element to save on re-computation of things
            // hence, skip elements that have already been treated in another batch
            if ( tLagElemHasNotBeenTreated( iLagElem ) )
            {
                // get the coarsest B-spline mesh
                moris_index tCoarsestBsplineMeshIndex = tCoarsestBsplineMesh( iLagElem );

                // get pointer to the B-spline mesh info for the coarsest element
                Bspline_Mesh_Info* tCoarsestBsplineMeshInfo =
                        aCutIntegrationMesh->mBsplineMeshInfos( tCoarsestBsplineMeshIndex );

                // get this coarsest B-spline element index
                const moris_index tCoarsestBsplineElemIndex =
                        tCoarsestBsplineMeshInfo->get_bspline_cell_index_for_extraction_cell( iLagElem );

                // get the list of Lagrange elements in this coarsest B-spline element
                const Cell< moris_index >& tLagElemsInBspElem =
                        tCoarsestBsplineMeshInfo->get_extraction_cell_indices_in_Bspline_cell( tCoarsestBsplineElemIndex );

                // get the subphase groups for the corresponding B-spline element
                const Cell< const Subphase_Group* > tSpgsOnCoarsestBspElem =
                        tCoarsestBsplineMeshInfo->get_SPGs_in_Bspline_cell( tCoarsestBsplineElemIndex );

                // collect SPs inside those SPGs and build map associating the SPs with MSD indices
                // Note: the local indices of the SPGs on the coarsest element correspond to the MSD indices
                IndexMap            tSpToMsdIndex;
                Cell< moris_index > tSpsInBspElem( 0 );
                tSpsInBspElem.reserve( tSpgsOnCoarsestBspElem.size() * tSpgsOnCoarsestBspElem( 0 )->get_SP_indices_in_group().size() );

                for ( uint iSPG = 0; iSPG < tSpgsOnCoarsestBspElem.size(); iSPG++ )
                {
                    // get the list of SPs in the current SPG
                    const moris::Cell< moris_index >& tSpIndicesInGroup =
                            tSpgsOnCoarsestBspElem( iSPG )->get_SP_indices_in_group();

                    // associate all SPs in the current SPG with their MSD index
                    for ( uint iSpInSpg = 0; iSpInSpg < tSpIndicesInGroup.size(); iSpInSpg++ )
                    {
                        tSpToMsdIndex[ tSpIndicesInGroup( iSpInSpg ) ] = iSPG;
                        tSpsInBspElem.push_back( tSpIndicesInGroup( iSpInSpg ) );
                    }
                }

                // get the size of the map
                uint tNumSpsInCoarsestBspElem = tSpsInBspElem.size();

                // use the constructed information to associate the SPGs on each B-spline mesh with the MSD indices
                for ( uint iBspMesh = 0; iBspMesh < tNumBspMeshes; iBspMesh++ )
                {
                    // get the current B-spline mesh info
                    Bspline_Mesh_Info* tBsplineMeshInfo = aCutIntegrationMesh->mBsplineMeshInfos( iBspMesh );

                    // go over the previously constructed map/list
                    for ( uint iSP = 0; iSP < tNumSpsInCoarsestBspElem; iSP++ )
                    {
                        // get the SP's index
                        moris_index tSpIndex = tSpsInBspElem( iSP );

                        // get the SPG corresponding to the current SP on the current mesh
                        moris_index tSpgIndex = tBsplineMeshInfo->mSpToSpgMap( tSpIndex );

                        // get the associated MSD index
                        auto tIter = tSpToMsdIndex.find( tSpIndex );
                        MORIS_ASSERT( tIter != tSpToMsdIndex.end(),
                                "Integration_Mesh_Generator::construct_SPG_material_connectivity_information() -"
                                "Something went wrong in the SP to MSD index map construction." );
                        moris_index tMsdIndex = tIter->second;

                        // add information to the SPG to MSD index map
                        tBsplineMeshInfo->mSpgToMsdIndex( tSpgIndex ) = tMsdIndex;
                    }
                }

                // --------------------------------
                // STEP 3: split the SPGs from every B-spline mesh related to a given Lag elem into void and material SPGs

                // treat all Lagrange elements in the B-spline element at the same time
                for ( uint iLagElemInBspElem = 0; iLagElemInBspElem < tLagElemsInBspElem.size(); iLagElemInBspElem++ )
                {
                    // get the Lagrange element index to be treated
                    moris_index tLagElemIndex = tLagElemsInBspElem( iLagElemInBspElem );

                    // check that the Lag elems in the current B-spline elem have not been treated yet
                    MORIS_ERROR( tLagElemHasNotBeenTreated( tLagElemIndex ),
                            "Integration_Mesh_Generator::construct_SPG_material_connectivity_information() - "
                            "Lagrange element in B-spline element marked as treated even though B-spline element has not been treated. "
                            "There must be a bug." );

                    // mark current Lagrange element as being treated
                    tLagElemHasNotBeenTreated( tLagElemIndex ) = false;

                    // the following information must be constructed wrt. each B-spline mesh
                    for ( uint iBspMesh = 0; iBspMesh < tNumBspMeshes; iBspMesh++ )
                    {
                        // get the current B-spline mesh info
                        Bspline_Mesh_Info* tBsplineMeshInfo = aCutIntegrationMesh->mBsplineMeshInfos( iBspMesh );

                        // access the lists storing which SPGs are material and void wrt. to the respective Lagrange element
                        moris::Cell< moris_index >& tMaterialSpgIndices = tBsplineMeshInfo->mExtractionCellMaterialSpgs( tLagElemIndex );
                        moris::Cell< moris_index >& tVoidSpgIndices     = tBsplineMeshInfo->mExtractionCellVoidSpgs( tLagElemIndex );

                        // fill these lists
                        this->find_material_SPGs_on_Lagrange_Cell( aCutIntegrationMesh, tLagElemIndex, iBspMesh, tMaterialSpgIndices, tVoidSpgIndices );

                        // --------------------------------
                        // STEP 4.1: derive the material MSD indices wrt. to each B-spline mesh from the void SPGs

                        // get the number of material SPGs and MSD indices
                        uint tNumMaterialSpgs = tMaterialSpgIndices.size();

                        // access the lists storing material MSD indices for each Lagrange element
                        moris::Cell< moris_index >& tMaterialMsdIndices = tBsplineMeshInfo->mExtractionCellMaterialMsdIndices( tLagElemIndex );

                        // resize this list to correct size
                        tMaterialMsdIndices.resize( tNumMaterialSpgs );

                        // go over material SPGs and store the corresponding MSD Indices
                        for ( uint iMaterialSpg = 0; iMaterialSpg < tNumMaterialSpgs; iMaterialSpg++ )
                        {
                            // get the current SPG's index
                            moris_index tSpgIndex = tMaterialSpgIndices( iMaterialSpg );

                            // get the corresponding MSD index
                            moris_index tMsdIndex = tBsplineMeshInfo->mSpgToMsdIndex( tSpgIndex );
                            MORIS_ASSERT( tMsdIndex > -1,
                                    "Integration_Mesh_Generator::construct_SPG_material_connectivity_information() - "
                                    "Subphase index not found in the SPG to MSD index map. Something must have gone wrong." );

                            // store the void MSD index
                            tMaterialMsdIndices( iMaterialSpg ) = tMsdIndex;
                        }

                        // --------------------------------
                        // STEP 4.2: derive the void MSD indices wrt. to each B-spline mesh from the void SPGs

                        // get the number of void SPGs and MSD indices
                        uint tNumVoidSpgs = tVoidSpgIndices.size();

                        // access the lists storing void MSD indices for each Lagrange element
                        moris::Cell< moris_index >& tVoidMsdIndices = tBsplineMeshInfo->mExtractionCellVoidMsdIndices( tLagElemIndex );

                        // resize this list to correct size
                        tVoidMsdIndices.resize( tNumVoidSpgs );

                        // go over void SPGs and store the corresponding MSD Indices
                        for ( uint iVoidSpg = 0; iVoidSpg < tNumVoidSpgs; iVoidSpg++ )
                        {
                            // get the current SPG's index
                            moris_index tSpgIndex = tVoidSpgIndices( iVoidSpg );

                            // get the corresponding MSD index
                            moris_index tMsdIndex = tBsplineMeshInfo->mSpgToMsdIndex( tSpgIndex );
                            MORIS_ASSERT( tMsdIndex > -1,
                                    "Integration_Mesh_Generator::construct_SPG_material_connectivity_information() - "
                                    "Subphase index not found in the SPG to MSD index map. Something must have gone wrong." );

                            // store the void MSD index
                            tVoidMsdIndices( iVoidSpg ) = tMsdIndex;
                        }

                    }    // end: loop over B-spline meshes for the current Lagrange element

                    // --------------------------------
                    // STEP 5: form a non-unique universal set of void MSD indices across all B-spline meshes

                    // initialize union set
                    aCutIntegrationMesh->mUnionVoidMsdIndices( tLagElemIndex ) =
                            aCutIntegrationMesh->mBsplineMeshInfos( 0 )->mExtractionCellVoidMsdIndices( tLagElemIndex );

                    // get the address to where the union multiset of MSD indices will be stored
                    moris::Cell< moris_index >* tUnionVoidMsdIndices = &( aCutIntegrationMesh->mUnionVoidMsdIndices( tLagElemIndex ) );

                    // form the union multiset of the void MSD indices across all B-spline meshes
                    for ( uint iBspMesh = 1; iBspMesh < tNumBspMeshes; iBspMesh++ )
                    {
                        moris::Cell< moris_index > tPreviousUnionVoidMsdIndices = *tUnionVoidMsdIndices;

                        // get the current B-spline mesh info
                        Bspline_Mesh_Info* tBsplineMeshInfo = aCutIntegrationMesh->mBsplineMeshInfos( iBspMesh );

                        // get the void MSD indices for the current B-spline mesh
                        moris::Cell< moris_index > tCurrentVoidMsdIndices =
                                tBsplineMeshInfo->mExtractionCellVoidMsdIndices( tLagElemIndex );

                        // form the union
                        xtk::multiset_union( tPreviousUnionVoidMsdIndices, tCurrentVoidMsdIndices, *tUnionVoidMsdIndices );
                    }

                    // --------------------------------
                    // STEP 6: deduce the free void MSD indices wrt. each B-spline mesh

                    // form the union multiset of the void MSD indices across all B-spline meshes
                    for ( uint iBspMesh = 0; iBspMesh < tNumBspMeshes; iBspMesh++ )
                    {
                        // get the current B-spline mesh info
                        Bspline_Mesh_Info* tBsplineMeshInfo = aCutIntegrationMesh->mBsplineMeshInfos( iBspMesh );

                        // get access to the list of void MSD indices
                        moris::Cell< moris_index >& tVoidMsdIndices = tBsplineMeshInfo->mExtractionCellVoidMsdIndices( tLagElemIndex );

                        // get access to the list of free void MSD indices
                        moris::Cell< moris_index >& tFreeVoidMsdIndices = tBsplineMeshInfo->mExtractionCellFreeVoidMsdIndices( tLagElemIndex );

                        // get the difference between the union and the void MSD indices to form the free void MSD indices
                        xtk::multiset_difference( *tUnionVoidMsdIndices, tVoidMsdIndices, tFreeVoidMsdIndices );
                    }

                    // --------------------------------
                    // STEP 7: find the bulk phases for union void MSDIs

                    // get the list of union void MSDIs for the current Lagrange element
                    moris::Cell< moris_index > const & tUnionVoidMsdIndicesFinal = aCutIntegrationMesh->mUnionVoidMsdIndices( tLagElemIndex );
                    uint                               tNumUnionVoidMSDIs        = tUnionVoidMsdIndicesFinal.size();

                    // initialize the list bulk-phase indices associated with the union void MSDIs
                    aCutIntegrationMesh->mUnionVoidMsdIndexBulkPhases( tLagElemIndex ).resize( tNumUnionVoidMSDIs );

                    // access the list of sub-phase groups associated with this Lagrange element on the coarsest B-spline mesh
                    moris::Cell< moris_index > const & tCoarsestSpgsAssociatedWithLagrangeCell =
                            tCoarsestBsplineMeshInfo->get_SPG_indices_associated_with_extraction_cell( tLagElemIndex );

                    // get the bulk-phases associated with each of the MSDIs
                    for ( uint iMSDI = 0; iMSDI < tNumUnionVoidMSDIs; iMSDI++ )
                    {
                        // get the MSD index
                        moris_index tMSDIndex = tUnionVoidMsdIndicesFinal( iMSDI );

                        // get the SPG associated with the current MSDI on the coarsest B-spline mesh
                        moris_index tAssociatedCoarsestSpg = tCoarsestSpgsAssociatedWithLagrangeCell( tMSDIndex );

                        // get the bulk-phase associated with this SPG
                        moris_index tMsdBulkPhase = tCoarsestBsplineMeshInfo->get_bulk_phase_for_subphase_group( tAssociatedCoarsestSpg );

                        // save this bulk phase as being associated with this union void MSDI
                        aCutIntegrationMesh->mUnionVoidMsdIndexBulkPhases( tLagElemIndex )( iMSDI ) = tMsdBulkPhase;
                    }

                }    // end for: loop over the Lagrange elements inside the current coarsest B-spline element

            }    // end if: only treat elements whose void MSD Indices have not been found yet

        }    // end for: loop over all lagrange elements

    }    // end function: Integration_Mesh_Generator::construct_SPG_material_connectivity_information(...)

    // ----------------------------------------------------------------------------------

    void
    Integration_Mesh_Generator::find_material_SPGs_on_Lagrange_Cell(
            Cut_Integration_Mesh*       aCutIntegrationMesh,
            const moris_index           aLagrangeElementIndex,
            const moris_index           aBsplineMeshListIndex,
            moris::Cell< moris_index >& aMaterialSpgIndices,
            moris::Cell< moris_index >& aVoidSpgIndices )
    {
        // get the number of SPs on the current IP cell
        const moris::Cell< moris_index >& tSPsOnCell = aCutIntegrationMesh->get_parent_cell_subphases( aLagrangeElementIndex );
        const uint                        tNumSPs    = tSPsOnCell.size();

        // initialize the list of material SPGs
        aMaterialSpgIndices.resize( tNumSPs );

        // get the pointer to the current B-spline mesh info
        Bspline_Mesh_Info* tBsplineMeshInfo = aCutIntegrationMesh->mBsplineMeshInfos( aBsplineMeshListIndex );

        // get the SPGs that are associated with the current IP cell
        moris::Cell< moris_index > const & tSPGsOnCell =
                tBsplineMeshInfo->get_SPG_indices_associated_with_extraction_cell( aLagrangeElementIndex );
        uint tNumSPGsOnCell = tSPGsOnCell.size();

        // initialize punch-card logging which SPGs have material on the current IP cell
        moris::Cell< bool > tVoidSPGs( tNumSPGsOnCell, true );

        // over the subphases on the current IP cell and mark the corresponding SPGs to have material
        for ( uint iSP = 0; iSP < tNumSPs; iSP++ )
        {
            // get the index of the current subphase
            moris_index tSpIndex = tSPsOnCell( iSP );

            // get the index of SPG the currently treated SP belongs to
            moris_index tSpgIndex = tBsplineMeshInfo->mSpToSpgMap( tSpIndex );

            // store SPG index containing material
            aMaterialSpgIndices( iSP ) = tSpgIndex;

            // find where the SPG is in the list of SPGs on the respective B-spline or Lagrange element
            moris_index tLocalSpgIndex = tBsplineMeshInfo->mSubphaseGroups( tSpgIndex )->get_local_index();

            // mark the SPG as having material in the punch card
            tVoidSPGs( tLocalSpgIndex ) = false;
        }

        // count the number of void IP cells that need to be constructed
        uint tNumVoidClusters = 0;
        for ( uint iSPG = 0; iSPG < tNumSPGsOnCell; iSPG++ )
        {
            tNumVoidClusters += tVoidSPGs( iSPG );
        }

        // set size of list of SPGs without material associated with current IP cell for the current B-spline mesh
        aVoidSpgIndices.resize( tNumVoidClusters );

        // store SPG indices for void clusters
        uint tVoidSpgCounter = 0;
        for ( uint iSPG = 0; iSPG < tNumSPGsOnCell; iSPG++ )
        {
            // if SPG has noted to not have material in it
            if ( tVoidSPGs( iSPG ) )
            {
                aVoidSpgIndices( tVoidSpgCounter ) = tSPGsOnCell( iSPG );
                tVoidSpgCounter++;
            }
        }
    }

    // ----------------------------------------------------------------------------------

}    // namespace xtk
