/*
 * cl_XTK_Model.cpp
 *
 *  Created on: Feb 18, 2019
 *      Author: doble
 */

#include "cl_XTK_Model.hpp"
#include "cl_XTK_Background_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"
#include "cl_XTK_Ghost_Stabilization.hpp"
#include "cl_XTK_Integration_Mesh_Generator.hpp"
#include "cl_XTK_Mesh_Cleanup.hpp"
#include "cl_XTK_Diagnostics.hpp"
//#include "cl_XTK_Contact_Sandbox.hpp"
#include "cl_XTK_Multigrid.hpp"
#include "fn_all_true.hpp"
#include "fn_unique.hpp"
#include "op_equal_equal.hpp"
#include "fn_sort.hpp"
#include "fn_iscol.hpp"
#include "fn_trans.hpp"
#include "fn_equal_to.hpp"
#include "fn_generate_element_to_element.hpp"
#include "fn_create_faces_from_element_to_node.hpp"
#include "fn_create_edges_from_element_to_node.hpp"
#include "HDF5_Tools.hpp"
#include "cl_MTK_Visualization_STK.hpp"
#include "cl_MTK_Cell_Info_Factory.hpp"
#include "cl_MTK_Cell_Info.hpp"
#include "cl_MTK_Writer_Exodus.hpp"
#include "cl_MTK_Periodic_Boundary_Condition_Helper.hpp"
#include "fn_Parsing_Tools.hpp"
#include "cl_TOL_Memory_Map.hpp"
#include "cl_Tracer.hpp"
#include "fn_stringify_matrix.hpp"

#include "cl_MTK_Intersection_Mesh.hpp"

using namespace moris;

namespace xtk
{
// ----------------------------------------------------------------------------------
// Constructor/Deconstructor Source code
// ----------------------------------------------------------------------------------

Model::~Model()
{
    delete mEnrichment;
    mEnrichment = nullptr;

    delete mGhostStabilization;
    mGhostStabilization = nullptr;

    for ( auto tIt : mEnrichedInterpMesh )
    {
        delete tIt;
    }

    mEnrichedInterpMesh.clear();

    for ( auto tIt : mEnrichedIntegMesh )
    {
        delete tIt;
    }

    mEnrichedIntegMesh.clear();

    if ( mIntersectionDetect != nullptr )
    {
        delete mIntersectionDetect;
    }

    if ( mIntersectionDetect2D != nullptr )
    {
        delete mIntersectionDetect2D;
    }
}

// ----------------------------------------------------------------------------------

/*
     * using the general geometry engine
     */
Model::Model(
    uint                            aModelDimension,
    moris::mtk::Interpolation_Mesh *aMeshData,
    moris::ge::Geometry_Engine *    aGeometryEngine,
    bool                            aLinkGeometryOnConstruction ) :
    mModelDimension( aModelDimension ),
    mBackgroundMesh( aMeshData ),
    mCutMesh( this, mModelDimension ),
    mGeometryEngine( aGeometryEngine ),
    mEnrichment( nullptr ),
    mGhostStabilization( nullptr ),
    mEnrichedInterpMesh( 0, nullptr ),
    mEnrichedIntegMesh( 0, nullptr ),
    mConvertedToTet10s( false )
{
    // flag this as a non-parameter list based run
    mParameterList.insert( "has_parameter_list", false );

}

// ----------------------------------------------------------------------------------

Model::Model( moris::ParameterList const &aParameterList ) :
    mParameterList( aParameterList ), mModelDimension( UINT_MAX ), mEnrichment( nullptr ), mGhostStabilization( nullptr ), mEnrichedInterpMesh( 0, nullptr ), mEnrichedIntegMesh( 0, nullptr ), mConvertedToTet10s( false )
{
    // flag this as a paramter list based run
    mParameterList.insert( "has_parameter_list", true );

    this->setup_diagnostics(
        mParameterList.get< bool >( "diagnostics" ),
        mParameterList.get< std::string >( "diagnostics_path" ),
        mParameterList.get< std::string >( "diagnostics_id" ) );
}

// ----------------------------------------------------------------------------------

void
Model::set_geometry_engine( moris::ge::Geometry_Engine *aGeometryEngine )
{
    mGeometryEngine = aGeometryEngine;
}

// ----------------------------------------------------------------------------------

void
Model::set_mtk_background_mesh( moris::mtk::Interpolation_Mesh *aMesh )
{
    this->initialize( aMesh );

    mInitializeCalled = true;
}

// ----------------------------------------------------------------------------------

void
Model::set_input_performer( std::shared_ptr< mtk::Mesh_Manager > aMTKPerformer )
{
    mMTKInputPerformer = aMTKPerformer;
}
// ----------------------------------------------------------------------------------

void
Model::set_cut_ig_mesh( std::shared_ptr< Cut_Integration_Mesh > aCutIgMesh )
{
    mCutIntegrationMesh = aCutIgMesh;
    mDecomposed         = true;
}
// ----------------------------------------------------------------------------------

void
Model::set_output_performer( std::shared_ptr< mtk::Mesh_Manager > aMTKPerformer )
{
    mMTKOutputPerformer = aMTKPerformer;
}

// ----------------------------------------------------------------------------------

void
Model::initialize( moris::mtk::Interpolation_Mesh *aMesh )
{   
    mBackgroundMesh     = aMesh;
    mModelDimension     = aMesh->get_spatial_dim();
    mCutMesh            = Cut_Mesh( this, mModelDimension );
    mEnrichment         = nullptr;
    mGhostStabilization = nullptr;
    mEnrichedInterpMesh = Cell< Enriched_Interpolation_Mesh * >( 0, nullptr );
    mEnrichedIntegMesh  = Cell< Enriched_Integration_Mesh * >( 0, nullptr );
    mConvertedToTet10s  = false;
}

// ----------------------------------------------------------------------------------

bool
Model::perform()
{
    Tracer tTracer( "XTK", "Overall", "Run" );

    mVerbose = mParameterList.get< bool >( "verbose" );

    if ( !mInitializeCalled )
    {
        MORIS_ERROR( mMTKInputPerformer != nullptr, "xtk::Model::perform(), mMTKInputPerformer not set!" );

        //FIXME hardcodes to mesh pair index 0
        moris::mtk::Interpolation_Mesh *tMesh = mMTKInputPerformer->get_interpolation_mesh( 0 );

        this->initialize( tMesh );
    }

    MORIS_ASSERT( this->has_parameter_list(), "Perform can only be called on a parameter list based XTK" );
    MORIS_ERROR( this->valid_parameters(), "Invalid parameters detected in XTK." );

    if ( mParameterList.get< std::string >( "probe_bg_cells" ) != "" )
    {
        Matrix< IdMat > tBgCellIds;
        string_to_mat( mParameterList.get< std::string >( "probe_bg_cells" ), tBgCellIds );
        print( tBgCellIds, "tBgCellIds" );
        this->probe_bg_cell( tBgCellIds );
    }

    if ( mParameterList.get< bool >( "decompose" ) )
    {
        mTriangulateAll = mParameterList.get< bool >( "triangulate_all" );

        mIgElementOrder = mParameterList.get< moris::uint >( "ig_element_order" );

        if ( mParameterList.get< bool >( "cleanup_cut_mesh" ) )
        {
            mCleanupMesh = true;
        }

        Cell< enum Subdivision_Method > tSubdivisionMethods = this->get_subdivision_methods();
        bool                            tSuccess            = this->decompose( tSubdivisionMethods );

        // Return false if cells are on different refinement levels
        if ( !tSuccess )
        {
            return false;
        }
    }

    if ( mParameterList.get< bool >( "cleanup_cut_mesh" ) )
    {
        mCleanupMesh = true;

        // cleanup the mesh
        Mesh_Cleanup tMeshCleanup( this, &mParameterList );
        tMeshCleanup.perform();
    }

    if ( mParameterList.get< bool >( "enrich" ) )
    {
        enum EntityRank tBasisRank = get_entity_rank_from_str( mParameterList.get< std::string >( "basis_rank" ) );

        Matrix< IndexMat > tMeshIndexCell;
        moris::string_to_mat( mParameterList.get< std::string >( "enrich_mesh_indices" ), tMeshIndexCell );

        this->perform_basis_enrichment( tBasisRank, tMeshIndexCell );

        // if high to low double side sets need to be created
        if ( mParameterList.get< bool >( "high_to_low_dbl_side_sets" ) )
        {
            for ( moris::uint i = 0; i < mGeometryEngine->get_num_bulk_phase(); i++ )
            {
                for ( moris::uint j = 0; j < mGeometryEngine->get_num_bulk_phase(); j++ )
                {
                    if ( i > j )
                    {
                        mEnrichedIntegMesh( 0 )->create_dbl_sided_interface_set( i, j );
                    }
                }
            }
        }

    }

    if ( mParameterList.get< bool >( "identify_hanging_nodes" ) )
    {
        this->perform_hanging_node_identification();
    }

    if ( mParameterList.get< bool >( "ghost_stab" ) )
    {
        this->construct_face_oriented_ghost_penalization_cells();

        if ( mParameterList.get< bool >( "visualize_ghost" ) )
        {
            Tracer tTracer( "XTK", "GhostStabilization", "Visualize" );

            for ( moris::moris_index i = 0; i < (moris_index)mGeometryEngine->get_num_bulk_phase(); i++ )
            {
                mGhostStabilization->visualize_ghost_on_mesh( i );
            }
        }
    }

    std::string tUnionBlockStr = mParameterList.get< std::string >( "union_blocks" );
    if ( !tUnionBlockStr.empty() )
    {
        // get the blocks to unionize
        moris::Cell< moris::Cell< std::string > > tUnionBlockCells;
        moris::string_to_cell_of_cell( tUnionBlockStr, tUnionBlockCells );

        // Row based
        Matrix< IndexMat > tUnionBlockColors      = string_to_mat< IndexMat >( mParameterList.get< std::string >( "union_block_colors" ) );
        std::string        tUnionNewBlockNamesStr = mParameterList.get< std::string >( "union_block_names" );

        moris::Cell< moris::Cell< std::string > > tNewBlockNames;
        moris::string_to_cell_of_cell( tUnionNewBlockNamesStr, tNewBlockNames );

        MORIS_ERROR( tUnionBlockCells.size() == tNewBlockNames.size(), "Dimension Mismatch in number of union operations for block" );
        MORIS_ERROR( tUnionBlockCells.size() == tUnionBlockColors.n_rows(), "Dimension Mismatch in number of union operations for block" );

        for ( moris::uint iUnion = 0; iUnion < tUnionBlockCells.size(); iUnion++ )
        {
            this->get_enriched_integ_mesh( 0 ).create_union_block( tUnionBlockCells( iUnion ), tNewBlockNames( iUnion )( 0 ), tUnionBlockColors.get_row( iUnion ) );
        }
    }

    std::string tUnionSideSetStr = mParameterList.get< std::string >( "union_side_sets" );
    if ( !tUnionSideSetStr.empty() )
    {
        // get the blocks to unionize
        moris::Cell< moris::Cell< std::string > > tUnionSideSetCells;
        moris::string_to_cell_of_cell( tUnionSideSetStr, tUnionSideSetCells );

        // Row based
        Matrix< IndexMat > tUnionSideSetColors      = string_to_mat< IndexMat >( mParameterList.get< std::string >( "union_side_set_colors" ) );
        std::string        tUnionNewSideSetNamesStr = mParameterList.get< std::string >( "union_side_set_names" );

        moris::Cell< moris::Cell< std::string > > tNewSideSetNames;
        moris::string_to_cell_of_cell( tUnionNewSideSetNamesStr, tNewSideSetNames );

        MORIS_ERROR( tUnionSideSetCells.size() == tNewSideSetNames.size(), "Dimension Mismatch in number of union operations for side set" );
        MORIS_ERROR( tUnionSideSetCells.size() == tUnionSideSetColors.n_rows(), "Dimension Mismatch in number of union operations for side set" );

        for ( moris::uint iUnion = 0; iUnion < tUnionSideSetCells.size(); iUnion++ )
        {
            this->get_enriched_integ_mesh( 0 ).create_union_side_set( tUnionSideSetCells( iUnion ), tNewSideSetNames( iUnion )( 0 ), tUnionSideSetColors.get_row( iUnion ) );
        }
    }

    std::string tDeactiveBlockstr = mParameterList.get< std::string >( "deactivate_all_but_blocks" );
    if ( !tDeactiveBlockstr.empty() )
    {
        // get the blocks to unionize
        moris::Cell< moris::Cell< std::string > > tBlocksToKeepStr;
        moris::string_to_cell_of_cell( tDeactiveBlockstr, tBlocksToKeepStr );

        MORIS_ERROR( tBlocksToKeepStr.size() == 1, "deactivate_all_but_block issue: This operation can only be performed on time" );

        this->get_enriched_integ_mesh( 0 ).deactive_all_blocks_but_selected( tBlocksToKeepStr( 0 ) );
    }

    std::string tDeactiveSideSetstr = mParameterList.get< std::string >( "deactivate_all_but_side_sets" );
    if ( !tDeactiveSideSetstr.empty() )
    {
        // get the blocks to unionize
        moris::Cell< moris::Cell< std::string > > tSideSetsToKeepStr;
        moris::string_to_cell_of_cell( tDeactiveSideSetstr, tSideSetsToKeepStr );

        MORIS_ERROR( tSideSetsToKeepStr.size() == 1, "deactive_all_side_sets_but_selected issue: This operation can only be performed on time" );

        this->get_enriched_integ_mesh( 0 ).deactive_all_side_sets_but_selected( tSideSetsToKeepStr( 0 ) );
    }

    if ( mParameterList.get< bool >( "multigrid" ) )
    {
        this->construct_multigrid();
    }

    if ( mEnriched )
    {
        xtk::Enriched_Interpolation_Mesh &tEnrInterpMesh = this->get_enriched_interp_mesh();
        xtk::Enriched_Integration_Mesh &  tEnrIntegMesh  = this->get_enriched_integ_mesh();
        // get meshes

        std::string tXTKMeshName = "XTKMesh";

        // place the pair in mesh manager
        mMTKOutputPerformer->register_mesh_pair( &tEnrInterpMesh, &tEnrIntegMesh, false, tXTKMeshName );

        //Periodic Boundary condition environment
        if ( mParameterList.get< std::string >( "periodic_side_set_pair" ) != "" )
        {
            //initialize the time tracer
            Tracer tTracer( "MTK", "Double Sided Set", " Periodic Boundary Condition " );

            //Construct intersection and perform
            mIntersectionDetect2D = new mtk::Intersection_Detect_2D( mMTKOutputPerformer, 0, mParameterList, mGeometryEngine->get_num_bulk_phase() );
            mIntersectionDetect2D->perform();
        }
        else
        {
            {
                //initialize the time tracer
                Tracer tTracer( "MTK", "Double Sided Set", " Periodic Boundary Condition " );


                mIntersectionDetect = new mtk::Intersection_Detect( mMTKOutputPerformer, 0, mParameterList, mGeometryEngine->get_num_bulk_phase() );
                mIntersectionDetect->perform();
            }

            {

                Tracer tTracer( "MTK", "Output Clusters", "Writing Mesh" );

                //Construct the intersection mesh
                mtk::Intersection_Mesh *tIscMesh = new mtk::Intersection_Mesh( &tEnrIntegMesh, mIntersectionDetect );

                //Write the mesh
                moris::mtk::Writer_Exodus tWriter2( tIscMesh );
                tWriter2.write_mesh( "", "VIS_ISC.exo", "", "temp.exo" );
                tWriter2.close_file();
            }
        }

        


        // if( mParameterList.get<bool>("contact_sandbox") )
        // {
        //     std::string tInterfaceSideSetName1 = tEnrIntegMesh.get_interface_side_set_name(0, 0, 2);
        //     std::string tInterfaceSideSetName2 = tEnrIntegMesh.get_interface_side_set_name(0, 1, 0);

        //     xtk::Contact_Sandbox tSandbox(&tEnrIntegMesh,
        //                                   tInterfaceSideSetName1,
        //                                   tInterfaceSideSetName2,
        //                                   mParameterList.get<real>("bb_epsilon"));

        //     // generate vertex displacement fields
        //     moris::real tInitialDisp = 0.0;
        //     moris::real tPredictedDisplX = 0.03;
        //     moris::real tPredictedDisplY = -0.03;
        //     moris::real tPredictedDisplZ = -0.01;
        //     Matrix<DDRMat> tCurrentDispl(tEnrIntegMesh.get_num_nodes(),this->get_spatial_dim(),tInitialDisp);
        //     Matrix<DDRMat> tPredictedDispl = tCurrentDispl;

        //     // get the vertices in bulk phase 1 and displace them through the current time step
        //     moris::mtk::Set * tSetC = tEnrIntegMesh.get_set_by_name( "HMR_dummy_c_p1");
        //     moris::mtk::Set * tSetN = tEnrIntegMesh.get_set_by_name( "HMR_dummy_n_p1");

        //     moris::Matrix< DDSMat > tVertsInChildBlock   = tSetC->get_ig_vertices_inds_on_block( true );
        //     moris::Matrix< DDSMat > tVertsInNoChildBlock = tSetN->get_ig_vertices_inds_on_block( true );

        //     // iterate through child verts block
        //     for(moris::uint i = 0; i < tVertsInChildBlock.numel(); i++)
        //     {
        //         moris_index tIndex = (moris_index)tVertsInChildBlock(i);
        //         tPredictedDispl(tIndex,0) = tInitialDisp + tPredictedDisplX;
        //         tPredictedDispl(tIndex,1) = tInitialDisp + tPredictedDisplY;
        //         if(this->get_spatial_dim() == 3)
        //         {
        //             tPredictedDispl(tIndex,2) = tInitialDisp + tPredictedDisplZ;
        //         }
        //     }

        //     // iterate through child verts block
        //     for(moris::uint i = 0; i < tVertsInNoChildBlock.numel(); i++)
        //     {
        //         moris_index tIndex = (moris_index)tVertsInNoChildBlock(i);
        //         tPredictedDispl(tIndex,0) = tInitialDisp + tPredictedDisplX;
        //         tPredictedDispl(tIndex,1) = tInitialDisp + tPredictedDisplY;
        //         if(this->get_spatial_dim() == 3)
        //         {
        //             tPredictedDispl(tIndex,2) = tInitialDisp + tPredictedDisplZ;
        //         }
        //     }

        //     tSandbox.perform_global_contact_search(tCurrentDispl,tPredictedDispl);
        // }

        if(mDiagnostics)
        {
            tEnrInterpMesh.write_diagnostics();
        }

        if ( mParameterList.get< bool >( "print_enriched_ig_mesh" ) )
        {
            tEnrIntegMesh.print();
        }

        if ( mParameterList.get< bool >( "exodus_output_XTK_ig_mesh" ) )
        {
            Tracer tTracer( "XTK", "Overall", "Visualize" );
            tEnrIntegMesh.write_mesh( &mParameterList );
        }
        // print the memory usage of XTK
        if ( mParameterList.get< bool >( "print_memory" ) )
        {
            moris::Memory_Map tXTKMM = this->get_memory_usage();
            tXTKMM.par_print( "XTK Model" );
        }
        // print
        MORIS_LOG_SPEC( "All_IG_verts", sum_all( tEnrIntegMesh.get_num_entities( EntityRank::NODE ) ) );
        MORIS_LOG_SPEC( "All_IG_cells", sum_all( tEnrIntegMesh.get_num_entities( EntityRank::ELEMENT ) ) );
        MORIS_LOG_SPEC( "All_IP_verts", sum_all( tEnrInterpMesh.get_num_entities( EntityRank::NODE ) ) );
        MORIS_LOG_SPEC( "All_IP_cells", sum_all( tEnrInterpMesh.get_num_entities( EntityRank::ELEMENT ) ) );
        MORIS_LOG_SPEC( "My_IG_verts", tEnrIntegMesh.get_num_entities( EntityRank::NODE ) );
        MORIS_LOG_SPEC( "My_IG_cells", tEnrIntegMesh.get_num_entities( EntityRank::ELEMENT ) );
        MORIS_LOG_SPEC( "My_IP_verts", tEnrInterpMesh.get_num_entities( EntityRank::NODE ) );
        MORIS_LOG_SPEC( "My_IP_cells", tEnrInterpMesh.get_num_entities( EntityRank::ELEMENT ) );
    }
    return true;
}


// ----------------------------------------------------------------------------------

bool
Model::has_parameter_list()
{
    return mParameterList.get< bool >( "has_parameter_list" );
}

// ----------------------------------------------------------------------------------

bool
Model::valid_parameters()
{
    bool tDecompose = mParameterList.get< bool >( "decompose" );
    bool tEnrich    = mParameterList.get< bool >( "enrich" );
    bool tGhost     = mParameterList.get< bool >( "ghost_stab" );
    bool tMultigrid = mParameterList.get< bool >( "multigrid" );

    if ( tEnrich == true )
    {
        MORIS_ERROR( tDecompose, "To perform basis enrichment, decomposition is also required." );
    }

    if ( tGhost == true )
    {
        MORIS_ERROR( tDecompose && tEnrich, "To perform ghost stabilization, decomposition and enrichment are also required." );
    }

    if ( tMultigrid == true )
    {
        MORIS_ERROR( tDecompose && tEnrich, "To perform multigrid, decomposition and enrichment are also required." );
    }

    return true;
}

// ----------------------------------------------------------------------------------

Cell< enum Subdivision_Method >
Model::get_subdivision_methods()
{
    MORIS_ASSERT( this->has_parameter_list(), "Perform can only be called on a parameter list based XTK" );

    moris::Cell< enum Subdivision_Method > tSubdivisionMethods;

    moris::uint       tSpatialDimension = this->get_spatial_dim();
    enum CellTopology tBGCellTopo       = this->get_parent_cell_topology();
    std::string       tDecompStr        = mParameterList.get< std::string >( "decomposition_type" );
    moris::lint       tOctreeRefLevel   = std::stoi( mParameterList.get< std::string >( "octree_refinement_level" ) );

    if ( tDecompStr.compare( "octree_only" ) == 0 )
    {
        return { Subdivision_Method::NC_OCTREE };
    }

    if ( tOctreeRefLevel > -1 )
    {
        tSubdivisionMethods.push_back( Subdivision_Method::NC_OCTREE );
    }

    // determine if we are going conformal or not
    bool tConformal = true;
    if ( tDecompStr.compare( "conformal" ) == 0 )
    {
        tConformal = true;
    }
    else if ( tDecompStr.compare( "nonconformal" ) == 0 )
    {
        tConformal = false;
    }
    else
    {
        MORIS_ERROR( 0, "Invalid decomposition_type provided. Recognized Options: Conformal and Nonconformal" );
    }

    if ( tSpatialDimension == 2 )
    {
        if ( tBGCellTopo == CellTopology::QUAD4 && tConformal )
        {
            moris::Cell< enum Subdivision_Method > tMethods = { Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3 };
            tSubdivisionMethods.append( tMethods );
        }
        else if ( tBGCellTopo == CellTopology::QUAD4 && !tConformal )
        {
            moris::Cell< enum Subdivision_Method > tMethods = { Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4 };
            tSubdivisionMethods.append( tMethods );
        }
    }
    else if ( tSpatialDimension == 3 )
    {
        if ( tBGCellTopo == CellTopology::HEX8 && tConformal )
        {
            moris::Cell< enum Subdivision_Method > tMethods = { Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4 };
            tSubdivisionMethods.append( tMethods );
        }
        else if ( tBGCellTopo == CellTopology::HEX8 && !tConformal )
        {
            moris::Cell< enum Subdivision_Method > tMethods = { Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8 };
            tSubdivisionMethods.append( tMethods );
        }
        else if ( tBGCellTopo == CellTopology::TET4 && tConformal )
        {
            moris::Cell< enum Subdivision_Method > tMethods = { Subdivision_Method::C_HIERARCHY_TET4 };
            tSubdivisionMethods.append( tMethods );
        }
    }
    else
    {
        MORIS_ASSERT( 0, "Invalid spatial dimension" );
    }

    return tSubdivisionMethods;
}

// ----------------------------------------------------------------------------------
// Decomposition Source code
// ----------------------------------------------------------------------------------

bool
Model::decompose( Cell< enum Subdivision_Method > aMethods )
{
    Tracer tTracer( "XTK", "Decomposition", "Decompose" );

    moris::Matrix< moris::IndexMat > tActiveGeometries( 1, mGeometryEngine->get_num_geometries() );

    for ( moris::uint i = 0; i < mGeometryEngine->get_num_geometries(); i++ )
    {
        tActiveGeometries( i ) = (moris_index)i;
    }

    Integration_Mesh_Generator tIntegrationGenerator( this, aMethods, tActiveGeometries );

    mCutIntegrationMesh = tIntegrationGenerator.perform();

    mDecomposed = mCutIntegrationMesh->mSameLevelChildMeshes;
    if(!mDecomposed)
    {
        MORIS_LOG_INFO( "Decomposition Failed. Not all child meshes on the same level." );
    }
    if ( mDiagnostics )
    {
        if ( interpolated_coordinate_check( mCutIntegrationMesh.get() ) )
        {
            MORIS_LOG_INFO( "Interpolated Coordinate Check: Pass" );
        }
        else
        {
            MORIS_LOG_INFO( "Interpolated Coordinate Check: Fail" );
        }
    }

    return mDecomposed;
}

// ----------------------------------------------------------------------------------

void
Model::create_new_node_association_with_geometry( Decomposition_Data &tDecompData )
{
//FIXME: REMOVE
    // // create geometry objects for each node
    // mGeometryEngine->create_new_child_nodes(
    //     tDecompData.tNewNodeIndex,
    //     tDecompData.tNewNodeParentTopology,
    //     tDecompData.tParamCoordRelativeToParent,
    //     mBackgroundMesh.get_all_node_coordinates_loc_inds() );
}

// ----------------------------------------------------------------------------------

bool
Model::verify_successful_node_assignment( Decomposition_Data &aDecompData )
{
    uint tNumUnsuccessful = 0;
    for ( moris::uint i = 0; i < aDecompData.tNewNodeId.size(); i++ )
    {
        if ( aDecompData.tNewNodeId( i ) == MORIS_INDEX_MAX )
        {
            tNumUnsuccessful++;
        }
    }

    if ( tNumUnsuccessful > 0 )
    {
        std::cout << "There were " << tNumUnsuccessful << " bad nodes of " << aDecompData.tNewNodeId.size() << " total nodes" << std::endl;
        return false;
    }

    return true;
}

// ----------------------------------------------------------------------------------

void
Model::send_outward_requests(
    moris_index const &         aMPITag,
    Cell< uint > const &        aProcRanks,
    Cell< Matrix< IndexMat > > &aOutwardRequests )
{
    // Cell of requests
    Cell< MPI_Request > tRequests( aProcRanks.size() );

    // iterate through owned requests and send
    for ( moris::uint i = 0; i < aProcRanks.size(); i++ )
    {
        tRequests( i ) = nonblocking_send(
            aOutwardRequests( i ),
            aOutwardRequests( i ).n_rows(),
            aOutwardRequests( i ).n_cols(),
            aProcRanks( i ),
            aMPITag );
    }
}

// ----------------------------------------------------------------------------------

void
Model::inward_receive_requests(
    moris_index const &         aMPITag,
    moris::uint                 aNumRows,
    Cell< Matrix< IndexMat > > &aReceivedData,
    Cell< uint > &              aProcRanksReceivedFrom )
{
    // ensure the sizes are correct.
    aReceivedData.resize( 0 );
    aProcRanksReceivedFrom.resize( 0 );

    // access the communication table
    Matrix< IdMat > tCommTable = mCutIntegrationMesh->get_communication_table();
    moris::uint     tCount     = 0;
    for ( moris::uint i = 0; i < tCommTable.numel(); i++ )
    {
        aReceivedData.push_back( Matrix< IndexMat >( 1, 1 ) );
        aProcRanksReceivedFrom.push_back( tCommTable( i ) );
        receive( aReceivedData( tCount ), aNumRows, tCommTable( i ), aMPITag );
        tCount++;
    }
}


// ----------------------------------------------------------------------------------

void
Model::inward_receive_request_answers(
    moris_index const &         aMPITag,
    moris::uint const &         aNumRows,
    Cell< uint > const &        aProcRanks,
    Cell< Matrix< IndexMat > > &aReceivedRequestAnswers )
{
    MPI_Status tStatus;

    for ( moris::uint i = 0; i < aProcRanks.size(); i++ )
    {
        bool tFlag = sent_message_exists( aProcRanks( i ), aMPITag, tStatus );
        while ( tFlag == false )
        {
            tFlag = sent_message_exists( aProcRanks( i ), aMPITag, tStatus );
        }

        aReceivedRequestAnswers.push_back( Matrix< IndexMat >( 1, 1 ) );
        receive( aReceivedRequestAnswers( i ), aNumRows, aProcRanks( i ), aMPITag );
    }
}

// ----------------------------------------------------------------------------------

void
Model::handle_received_request_answers(
    Decomposition_Data &              aDecompData,
    Cell< Matrix< IndexMat > > const &aRequests,
    Cell< Matrix< IndexMat > > const &aRequestAnswers,
    moris::moris_id &                 aNodeId )
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
                moris_index     tParentInd     = mBackgroundMesh->get_loc_entity_ind_from_entity_glb_id( tParentId, tParentRank );
                bool            tRequestExists = false;
                moris_index     tRequestIndex  = MORIS_INDEX_MAX;

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

    // handle the unhandled requests wiht current proc being the owner
    moris::moris_id tNodeId = mCutIntegrationMesh->allocate_entity_ids( tUnhandledRequestIndices.size(), EntityRank::NODE );

    for ( moris::uint i = 0; i < tUnhandledRequestIndices.size(); i++ )
    {
        moris_index tRequestIndex                  = tUnhandledRequestIndices( i );
        aDecompData.tNewNodeOwner( tRequestIndex ) = par_rank();
        aDecompData.tNewNodeId( tRequestIndex )    = tNodeId;
        tNodeId++;
    }
}

// ----------------------------------------------------------------------------------

void
Model::send_outward_requests_reals(
    moris_index const &       aMPITag,
    Cell< uint > const &      aProcRanks,
    Cell< Matrix< DDRMat > > &aOutwardRequests )
{
    // iterate through owned requests and send
    for ( moris::uint i = 0; i < aProcRanks.size(); i++ )
    {
        nonblocking_send(
            aOutwardRequests( i ),
            aOutwardRequests( i ).n_rows(),
            aOutwardRequests( i ).n_cols(),
            aProcRanks( i ),
            aMPITag );
    }
}

// ----------------------------------------------------------------------------------

void
Model::inward_receive_requests_reals(
    moris_index const &       aMPITag,
    moris::uint               aNumRows,
    Cell< Matrix< DDRMat > > &aReceivedData,
    Cell< uint > &            aProcRanksReceivedFrom )
{
    moris::moris_index tParRank = par_rank();
    moris::uint        tCount   = 0;
    MPI_Status         tStatus;
    for ( moris::uint i = 0; i < (moris::uint)par_size(); i++ )
    {
        if ( (moris_index)i != tParRank )
        {
            // if there is a sent message from a processor go receive it
            if ( sent_message_exists( i, aMPITag, tStatus ) )
            {
                aReceivedData.push_back( Matrix< DDRMat >( 1, 1 ) );

                aProcRanksReceivedFrom.push_back( i );

                receive( aReceivedData( tCount ), aNumRows, i, aMPITag );

                tCount++;
            }
        }
    }
}

// ----------------------------------------------------------------------------------

void
Model::return_request_answers_reals(
    moris_index const &             aMPITag,
    Cell< Matrix< DDRMat > > const &aRequestAnswers,
    Cell< uint > const &            aProcRanks )
{
    // iterate through owned requests and send
    for ( moris::uint i = 0; i < aProcRanks.size(); i++ )
    {
        nonblocking_send( aRequestAnswers( i ), aRequestAnswers( i ).n_rows(), aRequestAnswers( i ).n_cols(), aProcRanks( i ), aMPITag );
    }
}

// ----------------------------------------------------------------------------------

void
Model::inward_receive_request_answers_reals(
    moris_index const &       aMPITag,
    moris::uint const &       aNumRows,
    Cell< uint > const &      aProcRanks,
    Cell< Matrix< DDRMat > > &aReceivedData )
{
    for ( moris::uint i = 0; i < aProcRanks.size(); i++ )
    {
        aReceivedData.push_back( Matrix< DDRMat >( 1, 1 ) );

        receive( aReceivedData( i ), aNumRows, aProcRanks( i ), aMPITag );
    }
}

// ----------------------------------------------------------------------------------

void
Model::prepare_request_answers(
    Decomposition_Data &              aDecompData,
    Cell< Matrix< IndexMat > > const &aReceiveData,
    Cell< Matrix< IndexMat > > &      aRequestAnswers )
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
                moris_index     tParentInd     = mBackgroundMesh->get_loc_entity_ind_from_entity_glb_id( tParentId, tParentRank );
                bool            tRequestExists = false;
                moris_index     tRequestIndex  = MORIS_INDEX_MAX;

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
Model::return_request_answers(
    moris_index const &               aMPITag,
    Cell< Matrix< IndexMat > > const &aRequestAnswers,
    Cell< uint > const &              aProcRanks )
{
    // access the communication table
    Matrix< IdMat > tCommTable = mCutIntegrationMesh->get_communication_table();

    // iterate through owned requests and send
    for ( moris::uint i = 0; i < tCommTable.numel(); i++ )
    {
        nonblocking_send( aRequestAnswers( i ), aRequestAnswers( i ).n_rows(), aRequestAnswers( i ).n_cols(), tCommTable( i ), aMPITag );
    }
}

// ----------------------------------------------------------------------------------

moris::Cell< std::string >
Model::check_for_and_remove_internal_seacas_side_sets( moris::Cell< std::string > &aSideSetNames )
{
    for ( std::vector< std::string >::iterator iSet = aSideSetNames.begin(); iSet != aSideSetNames.end(); ++iSet )
    {
        if ( iSet->compare( "surface_1_quad4" ) == 0 )
        {
            aSideSetNames.data().erase( iSet-- );
        }

        else if ( iSet->compare( "surface_2_quad4" ) == 0 )
        {
            aSideSetNames.data().erase( iSet-- );
        }
        else if ( iSet->compare( "surface_3_quad4" ) == 0 )
        {
            aSideSetNames.data().erase( iSet-- );
        }
        else if ( iSet->compare( "surface_4_quad4" ) == 0 )
        {
            aSideSetNames.data().erase( iSet-- );
        }
        else if ( iSet->compare( "surface_5_quad4" ) == 0 )
        {
            aSideSetNames.data().erase( iSet-- );
        }
        else if ( iSet->compare( "surface_6_quad4" ) == 0 )
        {
            aSideSetNames.data().erase( iSet-- );
        }
        else if ( iSet->compare( "surface_hex8_quad_1" ) == 0 )
        {
            aSideSetNames.data().erase( iSet-- );
        }
        else if ( iSet->compare( "surface_hex8_quad_2" ) == 0 )
        {
            aSideSetNames.data().erase( iSet-- );
        }
        else if ( iSet->compare( "surface_hex8_quad4_1" ) == 0 )
        {
            aSideSetNames.data().erase( iSet-- );
        }
        else if ( iSet->compare( "surface_hex8_quad4_2" ) == 0 )
        {
            aSideSetNames.data().erase( iSet-- );
        }
        else if ( iSet->compare( "surface_hex8_quad4_3" ) == 0 )
        {
            aSideSetNames.data().erase( iSet-- );
        }
        else if ( iSet->compare( "surface_hex8_quad4_4" ) == 0 )
        {
            aSideSetNames.data().erase( iSet-- );
        }
    }

    return aSideSetNames;
}

// ----------------------------------------------------------------------------------
// Enrichment Source code
// ----------------------------------------------------------------------------------

void
Model::perform_basis_enrichment(
    enum EntityRank const &aBasisRank,
    moris_index const &    aMeshIndex )
{
    Tracer tTracer( "XTK", "Enrichment", "Enrich" );

    MORIS_ERROR( mDecomposed, "Prior to computing basis enrichment, the decomposition process must be called" );

    // allocate some new enriched interpolation and integration meshes
    mEnrichedIntegMesh.resize( aMeshIndex + 1, nullptr );
    mEnrichedInterpMesh.resize( aMeshIndex + 1, nullptr );

    this->perform_basis_enrichment_internal( aBasisRank, { { aMeshIndex } } );

    mEnrichment->write_diagnostics();

    // Change the enrichment flag
    mEnriched = true;
}


// ----------------------------------------------------------------------------------

void
Model::perform_basis_enrichment(
    enum EntityRank const &   aBasisRank,
    Matrix< IndexMat > const &aMeshIndex )
{
    Tracer tTracer( "XTK", "Enrichment" );

    MORIS_ERROR( mDecomposed, "Prior to computing basis enrichment, the decomposition process must be called" );

    // allocate some new enriched interpolation and integration meshes
    mEnrichedIntegMesh.resize( aMeshIndex.numel() + 1, nullptr );
    mEnrichedInterpMesh.resize( aMeshIndex.numel() + 1, nullptr );

    this->perform_basis_enrichment_internal( aBasisRank, aMeshIndex );

    if(mDiagnostics)
    {
        mEnrichment->write_diagnostics();
    }

    // Change the enrichment flag
    mEnriched = true;
}

void
Model::perform_hanging_node_identification()
{
    Tracer tTracer( "XTK", "Identify hanging nodes" );
    MORIS_ERROR(0,"NEEDS IMPLEMENTING");
    // MORIS_ERROR( mDecomposed, "Mesh needs to be decomposed prior to identifying hanging nodes" );

    // MORIS_ERROR( mEnriched, "Mesh needs to be enriched prior to identifying hanging nodes" );

    // // iterate through child meshes
    // for ( moris::uint iCM = 0; iCM < mCutMesh.get_num_child_meshes(); iCM++ )
    // {
    //     // active child mesh
    //     Child_Mesh &tChildMesh = mCutMesh.get_child_mesh( iCM );

    //     // get the neighbors
    //     Matrix< IndexMat > tElementNeighors = mBackgroundMesh->get_elements_connected_to_element_and_face_ind_loc_inds( tChildMesh.get_parent_element_index() );

    //     moris::Cell< moris_index > tTransitionFacets;

    //     // iterate through neighbor
    //     for ( moris::uint iN = 0; iN < tElementNeighors.n_cols(); iN++ )
    //     {
    //         moris_index tNeighborCellIndex = tElementNeighors( 0, iN );

    //         moris_index tSharedFaceIndex = tElementNeighors( 1, iN );

    //         if ( !mBackgroundMesh.entity_has_children( tNeighborCellIndex, EntityRank::ELEMENT ) )
    //         {
    //             tTransitionFacets.push_back( tSharedFaceIndex );
    //         }
    //     }

    //     tChildMesh.identify_hanging_nodes( tTransitionFacets );
    // }
}


void
Model::probe_bg_cell( Matrix< IndexMat > const &tBGCellIds )
{
    Tracer tTracer( "XTK", "BG Cell Probe" );

    for ( moris::uint i = 0; i < tBGCellIds.numel(); i++ )
    {
        Tracer tTracer( "XTK", "BG Cell Probe", "Cell Id " + std::to_string( tBGCellIds( i ) ) );

        moris_index                  tIndex      = mBackgroundMesh->get_loc_entity_ind_from_entity_glb_id( tBGCellIds( i ), EntityRank::ELEMENT );
        mtk::Cell &                  tCell       = mBackgroundMesh->get_mtk_cell( tIndex );
        Matrix< IndexMat >           tVertexIds  = tCell.get_vertex_ids();
        moris::Cell< mtk::Vertex * > tVertexPtrs = tCell.get_vertex_pointers();
        Matrix< IndexMat >           tVertexOwner( 1, tVertexPtrs.size() );

        MORIS_LOG_SPEC( "Cell Id", tBGCellIds( i ) );
        MORIS_LOG_SPEC( "Cell Index", tIndex );
        MORIS_LOG_SPEC( "Cell Owner", tCell.get_owner() );
        MORIS_LOG_SPEC( "Vertex Ids", ios::stringify_log( tVertexIds ) );
        MORIS_LOG_SPEC( "Vertex Owners", ios::stringify_log( tVertexOwner ) );

        // collect geometric info
        uint                     tNumGeom = mGeometryEngine->get_num_geometries();
        Cell< Matrix< DDRMat > > tVertexGeomVals( tNumGeom, Matrix< DDRMat >( 1, tVertexPtrs.size() ) );
        for ( moris::uint iG = 0; iG < tNumGeom; iG++ )
        {
            for ( moris::uint iV = 0; iV < tVertexPtrs.size(); iV++ )
            {
                tVertexGeomVals( iG )( iV ) = mGeometryEngine->get_field_value( iG, (uint)tVertexPtrs( iV )->get_index(), tVertexPtrs( iV )->get_coords() );
            }
            MORIS_LOG_SPEC( "Geom Field " + std::to_string( iG ) + " Vals", ios::stringify_log( tVertexGeomVals( iG ) ) );
        }
    }
}

// ----------------------------------------------------------------------------------
Cut_Integration_Mesh *
Model::get_cut_integration_mesh()
{
    MORIS_ASSERT( mDecomposed,
        "Cannot get cut integration mesh prior to the decomposition strategy " );

    return mCutIntegrationMesh.get();
}

Enrichment const &
Model::get_basis_enrichment()
{
    MORIS_ASSERT( mEnriched,
        "Cannot get basis enrichment from an XTK model which has not called perform_basis_enrichment " );

    return *mEnrichment;
}

// ----------------------------------------------------------------------------------

Enriched_Interpolation_Mesh &
Model::get_enriched_interp_mesh( moris::moris_index aIndex )
{
    MORIS_ASSERT( mEnriched,
        "Cannot get enriched interpolation mesh from an XTK model which has not called perform_basis_enrichment " );

    return *( mEnrichedInterpMesh( aIndex ) );
}

// ----------------------------------------------------------------------------------

Enriched_Integration_Mesh &
Model::get_enriched_integ_mesh( moris::moris_index aIndex )
{
    MORIS_ASSERT( mEnriched,
        "Cannot get enriched integration mesh from an XTK model which has not called perform_basis_enrichment " );

    return *( mEnrichedIntegMesh( aIndex ) );
}

// ----------------------------------------------------------------------------------

void
Model::perform_basis_enrichment_internal(
    enum EntityRank const &   aBasisRank,
    Matrix< IndexMat > const &aMeshIndex )
{

    // initialize enrichment (ptr because of circular dependency)
    mEnrichment = new Enrichment(
        Enrichment_Method::USE_INTERPOLATION_CELL_BASIS,
        aBasisRank,
        aMeshIndex,
        mGeometryEngine->get_num_phases(),
        this,
        mBackgroundMesh );

    // Set verbose flag to match XTK.
    mEnrichment->mVerbose = mVerbose;

    // perform the enrichment
    mEnrichment->perform_enrichment();
}

// ----------------------------------------------------------------------------------

void
Model::construct_face_oriented_ghost_penalization_cells()
{
    Tracer tTracer( "XTK", "GhostStabilization", "Stabilize" );

    MORIS_ERROR( mDecomposed, "Mesh needs to be decomposed prior to calling ghost penalization" );

    MORIS_ERROR( !mGhost, "Ghost penalization has already been called" );

    mGhostStabilization = new Ghost_Stabilization( this );

    mGhostStabilization->setup_ghost_stabilization();

    mGhost = true;
}

// ----------------------------------------------------------------------------------

Ghost_Stabilization &
Model::get_ghost_stabilization( moris::moris_index aIndex )
{
    MORIS_ERROR( mGhost, "Ghost has not been constructed on this model." );
    return *mGhostStabilization;
}

// ----------------------------------------------------------------------------------

void
Model::construct_multigrid()
{
    Tracer tTracer( "XTK", "Multigrid", "Run" );

    mMultigrid = std::make_shared< xtk::Multigrid >( this );

    mMultigrid->build_enriched_coeff_to_background_coeff_map();

    mMultigrid->create_fine_to_coarse_relationship();

    mMultigrid->create_coarse_to_fine_relationship();

    mMultigrid->create_coarse_to_fine_weights();

    std::string tName = "Enriched_bspline_1.exo";
    mMultigrid->build_basis_exodus_information( tName );
}

// ----------------------------------------------------------------------------------

Cut_Mesh &
Model::get_cut_mesh()
{
    return mCutMesh;
}

// ----------------------------------------------------------------------------------

Cut_Mesh const &
Model::get_cut_mesh() const
{
    return mCutMesh;
}

// ----------------------------------------------------------------------------------


moris::mtk::Interpolation_Mesh &
Model::get_background_mesh()
{
    return *mBackgroundMesh;
}

// ----------------------------------------------------------------------------------

moris::mtk::Interpolation_Mesh const &
Model::get_background_mesh() const
{
    return *mBackgroundMesh;
}

// ----------------------------------------------------------------------------------

moris::ge::Geometry_Engine *
Model::get_geom_engine()
{
    return mGeometryEngine;
}


bool
Model::subphase_is_in_child_mesh( moris_index aSubphaseIndex )
{
    MORIS_ERROR( 0, "TOREMOVE" );
    return false;
}

// ----------------------------------------------------------------------------------

moris::ParameterList &
Model::get_parameter_list()
{
    return mParameterList;
}

//------------------------------------------------------------------------------
void
Model::setup_diagnostics(
    bool               aDiagnostics,
    std::string const &aDiagnosticPath,
    std::string const &aDiagnosticLabel )
{
    mDiagnostics = aDiagnostics;

    if ( mDiagnostics )
    {
        mDiagnosticPath = aDiagnosticPath;
        mDiagnosticId   = aDiagnosticLabel;

        MORIS_ERROR( !mDiagnosticPath.empty(), "If diagnostics are turned on, a diagnostics path must be specified" );
        if ( mDiagnosticId.empty() )
        {
            mDiagnosticId = "no_spec";
        }
    }
}
std::string
Model::get_diagnostic_file_name( std::string const &aLabel ) const
{
    MORIS_ASSERT( mDiagnostics, "Only callable with diagnostics on" );
    return mDiagnosticPath + "/id_" + mDiagnosticId + "_ps_" + std::to_string( moris::par_size() ) + "_pr_" + std::to_string( moris::par_rank() ) + "_" + aLabel + ".csv";
}

bool
Model::triangulate_all()
{
    return mTriangulateAll;
}

uint
Model::ig_element_order()
{
    return mIgElementOrder;
}

//------------------------------------------------------------------------------

moris::mtk::Integration_Mesh *
Model::get_output_mesh( Output_Options const &aOutputOptions )

{
    MORIS_ERROR( 0, "Deprecated. (Removal in progress)" );

    return nullptr;
}

//------------------------------------------------------------------------------

moris::uint
Model::get_spatial_dim() const
{
    return mModelDimension;
}

//------------------------------------------------------------------------------

moris_index
Model::get_cell_xtk_index( moris_id aCellId )
{
    auto tIter = mCellGlbToLocalMap.find( aCellId );
    MORIS_ASSERT( tIter != mCellGlbToLocalMap.end(), "Id not in map" );
    return tIter->second;
}

//------------------------------------------------------------------------------

moris::Cell< moris::Cell< moris_index > > const &
Model::get_subphase_to_subphase()
{
    MORIS_ERROR( 0, "Deprecated." );
    return mSubphaseToSubPhase;
}

//------------------------------------------------------------------------------

moris::Cell< moris::Cell< moris_index > > const &
Model::get_subphase_to_subphase_my_side_ords()
{
    MORIS_ERROR( 0, "Deprecated." );
    return mSubphaseToSubPhaseMySideOrds;
}

//------------------------------------------------------------------------------

moris::Cell< moris::Cell< moris_index > > const &
Model::get_subphase_to_subphase_transition_loc()
{
    MORIS_ERROR( 0, "Deprecated." );
    return mTransitionNeighborCellLocation;
}

//------------------------------------------------------------------------------

moris::Cell< moris::Cell< moris_index > > const &
Model::get_subphase_to_subphase_neighbor_side_ords()
{
    MORIS_ERROR( 0, "Deprecated." );
    return mSubphaseToSubPhaseNeighborSideOrds;
}

//------------------------------------------------------------------------------

std::shared_ptr< Multigrid >
Model::get_multigrid_ptr()
{
    return mMultigrid;
}

//------------------------------------------------------------------------------
moris::Matrix< moris::IndexMat >
Model::get_num_subphase_neighbors()
{
    moris::Cell< moris::Cell< moris_index > > const &tSubPhaseToSubphase = this->get_subphase_to_subphase();
    moris::Matrix< moris::IndexMat >                 tSubphaseNumNeighbors( 1, tSubPhaseToSubphase.size() );
    for ( size_t iSP = 0; iSP < tSubPhaseToSubphase.size(); iSP++ )
    {
        tSubphaseNumNeighbors( iSP ) = tSubPhaseToSubphase( iSP ).size();
    }
    return tSubphaseNumNeighbors;
}
//------------------------------------------------------------------------------

moris_id
Model::get_subphase_id( moris_id aSubphaseIndex )
{

    return mCutIntegrationMesh->get_subphase_id( aSubphaseIndex );
}

// ----------------------------------------------------------------------------------

moris_index
Model::get_subphase_index( moris_id aSubphaseId )
{
    return mCutIntegrationMesh->get_subphase_index( aSubphaseId );
}

//------------------------------------------------------------------------------

moris::Memory_Map
Model::get_memory_usage()
{
    // memory map for model
    moris::Memory_Map tXTKModelMM;

    // member data that have memory maps
    moris::Memory_Map tCutMeshMM;
    moris::Memory_Map tBGMeshMM;
    moris::Memory_Map tEnrichmentMM;
    moris::Memory_Map tGhostMM;
    moris::Memory_Map tIgMeshMM;
    moris::Memory_Map tIpMeshMM;

    if ( mDecomposed )
    {
        tCutMeshMM = mCutMesh.get_memory_usage();
    }

    if ( mEnriched )
    {
        tEnrichmentMM = mEnrichment->get_memory_usage();
        tIgMeshMM     = this->get_enriched_integ_mesh().get_memory_usage();
        tIpMeshMM     = this->get_enriched_interp_mesh().get_memory_usage();
    }

    if ( mGhost )
    {
        tGhostMM = mGhostStabilization->get_memory_usage();
    }

    // make the sum of the cut mesh memory map the cut mesh memory
    tXTKModelMM.mMemoryMapData["Cut Mesh"]                            = tCutMeshMM.sum();
    tXTKModelMM.mMemoryMapData["Enrichment"]                          = tEnrichmentMM.sum();
    tXTKModelMM.mMemoryMapData["Enriched Ig Mesh"]                    = tIgMeshMM.sum();
    tXTKModelMM.mMemoryMapData["Enriched Ip Mesh"]                    = tIpMeshMM.sum();
    tXTKModelMM.mMemoryMapData["Ghost"]                               = tGhostMM.sum();
    tXTKModelMM.mMemoryMapData["Background Mesh"]                     = tBGMeshMM.sum();
    tXTKModelMM.mMemoryMapData["mElementToElement ptrs"]              = moris::internal_capacity( mElementToElement );
    tXTKModelMM.mMemoryMapData["mElementToElement ptrs"]              = moris::internal_capacity( mElementToElement );
    tXTKModelMM.mMemoryMapData["mSubphaseToSubPhase"]                 = moris::internal_capacity( mSubphaseToSubPhase );
    tXTKModelMM.mMemoryMapData["mSubphaseToSubPhaseMySideOrds"]       = moris::internal_capacity( mSubphaseToSubPhaseMySideOrds );
    tXTKModelMM.mMemoryMapData["mSubphaseToSubPhaseNeighborSideOrds"] = moris::internal_capacity( mSubphaseToSubPhaseNeighborSideOrds );

    tIgMeshMM.par_print( "Ig Mesh" );
    return tXTKModelMM;
}

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

}// namespace xtk
