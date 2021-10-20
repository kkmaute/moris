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
#include "cl_XTK_Mesh_Cleanup.hpp"
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

        for(auto tIt:mEnrichedInterpMesh)
        {
            delete tIt;
        }

        mEnrichedInterpMesh.clear();

        for(auto tIt:mEnrichedIntegMesh)
        {
            delete tIt;
        }

        mEnrichedIntegMesh.clear();

        if( mIntersectionDetect != nullptr)
        {
            delete mIntersectionDetect;
        }

        if( mIntersectionDetect2D != nullptr)
        {
            delete mIntersectionDetect2D;
        }

    }

    // ----------------------------------------------------------------------------------

    /*
     * using the general geometry engine
     */
    Model::Model(
            uint aModelDimension,
            moris::mtk::Interpolation_Mesh* aMeshData,
            moris::ge::Geometry_Engine* aGeometryEngine,
            bool aLinkGeometryOnConstruction )
    : mSameMesh(false),
      mModelDimension(aModelDimension),
      mBackgroundMesh(aMeshData,aGeometryEngine),
      mCutMesh(this,mModelDimension),
      mGeometryEngine(aGeometryEngine),
      mEnrichment(nullptr),
      mGhostStabilization(nullptr),
      mEnrichedInterpMesh(0,nullptr),
      mEnrichedIntegMesh(0,nullptr),
      mConvertedToTet10s(false)
    {
        // flag this as a non-parameter list based run
        mParameterList.insert("has_parameter_list", false);

        mBackgroundMesh.initialize_interface_node_flags(
                mBackgroundMesh.get_num_entities(EntityRank::NODE),
                mGeometryEngine->get_num_geometries());


    }

    // ----------------------------------------------------------------------------------

    Model::Model( moris::ParameterList const & aParameterList )
    : mSameMesh(false),
      mParameterList(aParameterList),
      mModelDimension(UINT_MAX),
      mEnrichment(nullptr),
      mGhostStabilization(nullptr),
      mEnrichedInterpMesh(0,nullptr),
      mEnrichedIntegMesh(0,nullptr),
      mConvertedToTet10s(false)
    {
        // flag this as a paramter list based run
        mParameterList.insert("has_parameter_list", true);
    }

    // ----------------------------------------------------------------------------------

    void
    Model::set_geometry_engine(moris::ge::Geometry_Engine* aGeometryEngine)
    {
        mGeometryEngine = aGeometryEngine;
    }

    // ----------------------------------------------------------------------------------

    void
    Model::set_mtk_background_mesh(moris::mtk::Interpolation_Mesh* aMesh)
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
    Model::set_output_performer( std::shared_ptr< mtk::Mesh_Manager > aMTKPerformer )
    {
        mMTKOutputPerformer = aMTKPerformer;
    }

    // ----------------------------------------------------------------------------------

    void
    Model::initialize( moris::mtk::Interpolation_Mesh* aMesh )
    {
        mSameMesh           = false;
        mModelDimension     = aMesh->get_spatial_dim();
        mCutMesh            = Cut_Mesh(this,mModelDimension);
        mEnrichment         = nullptr;
        mGhostStabilization = nullptr;
        mEnrichedInterpMesh = Cell<Enriched_Interpolation_Mesh*>(0, nullptr);
        mEnrichedIntegMesh  = Cell<Enriched_Integration_Mesh*>(0, nullptr);
        mConvertedToTet10s  = false;
        mBackgroundMesh = Background_Mesh(aMesh,mGeometryEngine);
        mBackgroundMesh.initialize_interface_node_flags(
                mBackgroundMesh.get_num_entities(EntityRank::NODE),
                mGeometryEngine->get_num_geometries());
    }

    // ----------------------------------------------------------------------------------

    bool
    Model::perform()
    {
        Tracer tTracer( "XTK", "Overall", "Run" );

        mVerbose = mParameterList.get<bool>("verbose");


        if( !mInitializeCalled )
        {
            MORIS_ERROR( mMTKInputPerformer != nullptr ,"xtk::Model::perform(), mMTKInputPerformer not set!");

            //FIXME hardcodes to mesh pair index 0
            moris::mtk::Interpolation_Mesh* tMesh = mMTKInputPerformer->get_interpolation_mesh( 0 );

            this->initialize( tMesh );
        }

        MORIS_ASSERT(this->has_parameter_list(),"Perform can only be called on a parameter list based XTK");
        MORIS_ERROR(this->valid_parameters(),"Invalid parameters detected in XTK.");

        if(mParameterList.get<std::string>("probe_bg_cells") != "")
        {
            Matrix<IdMat> tBgCellIds;
            string_to_mat( mParameterList.get< std::string >("probe_bg_cells"), tBgCellIds );
            print(tBgCellIds,"tBgCellIds");
            this->probe_bg_cell(tBgCellIds);
        }

        if(mParameterList.get<bool>("decompose"))
        {
            mTriangulateAll = mParameterList.get<bool>("triangulate_all");

            if(mParameterList.get<bool>("cleanup_cut_mesh"))
            {
                mCleanupMesh = true;
            }

            Cell<enum Subdivision_Method> tSubdivisionMethods = this->get_subdivision_methods();
            bool tSuccess = this->decompose(tSubdivisionMethods);

            // Return false if cells are on different refinement levels
            if(!tSuccess)
            {
                return false;
            }
        }

        if(mParameterList.get<bool>("cleanup_cut_mesh"))
        {
            mCleanupMesh = true;

            // cleanup the mesh
            Mesh_Cleanup tMeshCleanup(this,&mParameterList);
            tMeshCleanup.perform();
        }

        // at this point the cut mesh and background mesh is not going to change anymore
        this->finalize_mesh_data();


        if(mParameterList.get<bool>("enrich"))
        {
            enum EntityRank tBasisRank = get_entity_rank_from_str(mParameterList.get<std::string>("basis_rank"));

            Matrix<IndexMat> tMeshIndexCell;
            moris::string_to_mat(mParameterList.get< std::string >( "enrich_mesh_indices" ), tMeshIndexCell);

            this->perform_basis_enrichment(tBasisRank,tMeshIndexCell);

            // if high to low double side sets need to be created
            if(mParameterList.get<bool>("high_to_low_dbl_side_sets"))
            {
                for(moris::uint i = 0; i < mGeometryEngine->get_num_bulk_phase(); i++)
                {
                    for(moris::uint j = 0; j < mGeometryEngine->get_num_bulk_phase(); j++)
                    {
                        if(i > j)
                        {
                            mEnrichedIntegMesh(0)->create_dbl_sided_interface_set( i, j );
                        }
                    }
                }
            }
        }

        if(mParameterList.get<bool>("identify_hanging_nodes"))
        {
            this->perform_hanging_node_identification();
        }

        if(mParameterList.get<bool>("ghost_stab"))
        {
            this->construct_face_oriented_ghost_penalization_cells();

            if( mParameterList.get<bool>("visualize_ghost") )
            {
                Tracer tTracer( "XTK", "GhostStabilization", "Visualize" );

                for(moris::moris_index i = 0; i < (moris_index)mGeometryEngine->get_num_bulk_phase(); i++)
                {
                    mGhostStabilization->visualize_ghost_on_mesh(i);
                }
            }
        }

        std::string tUnionBlockStr = mParameterList.get<std::string>("union_blocks");
        if( !tUnionBlockStr.empty())
        {
            // get the blocks to unionize
            moris::Cell< moris::Cell< std::string > > tUnionBlockCells;
            moris::string_to_cell_of_cell(tUnionBlockStr,tUnionBlockCells);

            // Row based
            Matrix<IndexMat> tUnionBlockColors = string_to_mat<IndexMat>(mParameterList.get<std::string>("union_block_colors"));
            std::string tUnionNewBlockNamesStr = mParameterList.get<std::string>("union_block_names");

            moris::Cell< moris::Cell< std::string > > tNewBlockNames;
            moris::string_to_cell_of_cell(tUnionNewBlockNamesStr,tNewBlockNames);

            MORIS_ERROR(tUnionBlockCells.size() == tNewBlockNames.size(),"Dimension Mismatch in number of union operations for block");
            MORIS_ERROR(tUnionBlockCells.size() == tUnionBlockColors.n_rows(),"Dimension Mismatch in number of union operations for block");

            for(moris::uint iUnion = 0; iUnion < tUnionBlockCells.size(); iUnion++)
            {
                this->get_enriched_integ_mesh(0).create_union_block(tUnionBlockCells(iUnion),tNewBlockNames(iUnion)(0),tUnionBlockColors.get_row(iUnion));
            }

        }

        std::string tUnionSideSetStr = mParameterList.get<std::string>("union_side_sets");
        if( !tUnionSideSetStr.empty())
        {
            // get the blocks to unionize
            moris::Cell< moris::Cell< std::string > > tUnionSideSetCells;
            moris::string_to_cell_of_cell(tUnionSideSetStr,tUnionSideSetCells);

            // Row based
            Matrix<IndexMat> tUnionSideSetColors = string_to_mat<IndexMat>(mParameterList.get<std::string>("union_side_set_colors"));
            std::string tUnionNewSideSetNamesStr = mParameterList.get<std::string>("union_side_set_names");

            moris::Cell< moris::Cell< std::string > > tNewSideSetNames;
            moris::string_to_cell_of_cell(tUnionNewSideSetNamesStr,tNewSideSetNames);

            MORIS_ERROR(tUnionSideSetCells.size() == tNewSideSetNames.size(),"Dimension Mismatch in number of union operations for side set");
            MORIS_ERROR(tUnionSideSetCells.size() == tUnionSideSetColors.n_rows(),"Dimension Mismatch in number of union operations for side set");

            for(moris::uint iUnion = 0; iUnion < tUnionSideSetCells.size(); iUnion++)
            {
                this->get_enriched_integ_mesh(0).create_union_side_set(tUnionSideSetCells(iUnion),tNewSideSetNames(iUnion)(0),tUnionSideSetColors.get_row(iUnion));
            }

        }

        std::string tDeactiveBlockstr = mParameterList.get<std::string>("deactivate_all_but_blocks");
        if( !tDeactiveBlockstr.empty())
        {
            // get the blocks to unionize
            moris::Cell< moris::Cell< std::string > > tBlocksToKeepStr;
            moris::string_to_cell_of_cell(tDeactiveBlockstr,tBlocksToKeepStr);

            MORIS_ERROR(tBlocksToKeepStr.size() == 1,"deactivate_all_but_block issue: This operation can only be performed on time");

            this->get_enriched_integ_mesh(0).deactive_all_blocks_but_selected(tBlocksToKeepStr(0));
        }

        std::string tDeactiveSideSetstr = mParameterList.get<std::string>("deactivate_all_but_side_sets");
        if( !tDeactiveSideSetstr.empty())
        {
            // get the blocks to unionize
            moris::Cell< moris::Cell< std::string > > tSideSetsToKeepStr;
            moris::string_to_cell_of_cell(tDeactiveSideSetstr,tSideSetsToKeepStr);

            MORIS_ERROR(tSideSetsToKeepStr.size() == 1,"deactive_all_side_sets_but_selected issue: This operation can only be performed on time");

            this->get_enriched_integ_mesh(0).deactive_all_side_sets_but_selected(tSideSetsToKeepStr(0));
        }


        if( mParameterList.get<bool>("multigrid") )
        {
            this->construct_multigrid();
        }

        if(mEnriched)
        {
            // get meshes
            xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = this->get_enriched_interp_mesh();
            xtk::Enriched_Integration_Mesh   & tEnrIntegMesh  = this->get_enriched_integ_mesh();

            std::string tXTKMeshName = "XTKMesh";

            // place the pair in mesh manager
            mMTKOutputPerformer->register_mesh_pair( &tEnrInterpMesh, &tEnrIntegMesh, false, tXTKMeshName );

            //Periodic Boundary condition environment
            if( mParameterList.get< std::string >( "periodic_side_set_pair" ) != "" )
            {


                if(tEnrInterpMesh.get_spatial_dim() == 2 )
                {
                    //initialize the time tracer
                    Tracer tTracer( "MTK", "Double Sided Set", " Periodic Boundary Condition ");

                    //Construct intersection and perform
                    mIntersectionDetect2D = new mtk::Intersection_Detect_2D( mMTKOutputPerformer, 0, mParameterList, mGeometryEngine->get_num_bulk_phase()) ;
                    mIntersectionDetect2D->perform();
                }
                else
                {
                    {
                        //initialize the time tracer
                        Tracer tTracer( "MTK", "Double Sided Set", " Periodic Boundary Condition ");

                        mIntersectionDetect = new mtk::Intersection_Detect( mMTKOutputPerformer, 0, mParameterList, mGeometryEngine->get_num_bulk_phase()) ;

                        mIntersectionDetect->perform();
                    }

                    {
                        Tracer tTracer( "MTK", "Output Clusters", "Writing Mesh");

                        //Construct the intersection mesh
                        mtk::Intersection_Mesh* tIscMesh = new mtk::Intersection_Mesh(&tEnrIntegMesh,mIntersectionDetect );

                        //Write the mesh
                        moris::mtk::Writer_Exodus tWriter2(tIscMesh);
                        tWriter2.write_mesh("", "VIS_ISC.exo", "", "temp.exo");
                        tWriter2.close_file();
                    }


                }
                //constrcut the object for periodic boundary condition
                //                mtk::Periodic_Boundary_Condition_Helper tPBCHelper(mMTKOutputPerformer,0, mParameterList);
                //
                //                //perform periodic boundary condition
                //                tPBCHelper.setup_periodic_boundary_conditions();
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

            if( mParameterList.get<bool>("print_enriched_ig_mesh") )
            {
                tEnrIntegMesh.print();
            }

            if( mParameterList.get<bool>("exodus_output_XTK_ig_mesh") )
            {
                Tracer tTracer( "XTK", "Overall", "Visualize" );
                tEnrIntegMesh.write_mesh(&mParameterList);
            }

            // print the memory usage of XTK
            if( mParameterList.get<bool>("print_memory") )
            {
                moris::Memory_Map tXTKMM = this->get_memory_usage();
                tXTKMM.par_print("XTK Model");
            }

            // print 
            MORIS_LOG_SPEC("All_IG_verts",sum_all(tEnrIntegMesh.get_num_entities(EntityRank::NODE)));
            MORIS_LOG_SPEC("All_IG_cells",sum_all(tEnrIntegMesh.get_num_entities(EntityRank::ELEMENT)));
            MORIS_LOG_SPEC("All_IP_verts",sum_all(tEnrInterpMesh.get_num_entities(EntityRank::NODE)));
            MORIS_LOG_SPEC("All_IP_cells",sum_all(tEnrInterpMesh.get_num_entities(EntityRank::ELEMENT)));
            MORIS_LOG_SPEC("My_IG_verts",tEnrIntegMesh.get_num_entities(EntityRank::NODE));
            MORIS_LOG_SPEC("My_IG_cells",tEnrIntegMesh.get_num_entities(EntityRank::ELEMENT));
            MORIS_LOG_SPEC("My_IP_verts",tEnrInterpMesh.get_num_entities(EntityRank::NODE));
            MORIS_LOG_SPEC("My_IP_cells",tEnrInterpMesh.get_num_entities(EntityRank::ELEMENT));
        }

        return true;
    }

    // ----------------------------------------------------------------------------------

    bool
    Model::has_parameter_list()
    {
        return mParameterList.get<bool>("has_parameter_list");
    }

    // ----------------------------------------------------------------------------------

    bool
    Model::valid_parameters()
    {
        bool tDecompose = mParameterList.get<bool>("decompose");
        bool tEnrich    = mParameterList.get<bool>("enrich");
        bool tGhost     = mParameterList.get<bool>("ghost_stab");
        bool tMultigrid = mParameterList.get<bool>("multigrid");

        if(tEnrich == true)
        {
            MORIS_ERROR(tDecompose, "To perform basis enrichment, decomposition is also required.");
        }

        if(tGhost == true)
        {
            MORIS_ERROR(tDecompose && tEnrich, "To perform ghost stabilization, decomposition and enrichment are also required.");
        }

        if(tMultigrid == true)
        {
            MORIS_ERROR(tDecompose && tEnrich, "To perform multigrid, decomposition and enrichment are also required.");
        }

        return true;
    }

    // ----------------------------------------------------------------------------------

    Cell<enum Subdivision_Method>
    Model::get_subdivision_methods()
    {
        MORIS_ASSERT(this->has_parameter_list(),"Perform can only be called on a parameter list based XTK");

        moris::uint       tSpatialDimension = this->get_spatial_dim();
        enum CellTopology tBGCellTopo       = mBackgroundMesh.get_parent_cell_topology();
        std::string       tDecompStr        = mParameterList.get<std::string>("decomposition_type");

        // determine if we are going conformal or not
        bool tConformal = true;
        if(tDecompStr.compare("conformal") == 0 )
        {
            tConformal = true;
        }
        else if(tDecompStr.compare("nonconformal") == 0 )
        {
            tConformal = true;
        }
        else
        {
            MORIS_ERROR(0,"Invalid decomposition_type provided. Recognized Options: Conformal and Nonconformal");
        }

        if(tSpatialDimension == 2 )
        {
            if(tBGCellTopo == CellTopology::QUAD4  && tConformal)
            {
                return {Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3};
            }
            else if(tBGCellTopo == CellTopology::QUAD4  && !tConformal)
            {
                return {Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4};
            }
        }
        else if ( tSpatialDimension == 3 )
        {
            if(tBGCellTopo == CellTopology::HEX8  && tConformal)
            {
                return {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4};
            }
            else if(tBGCellTopo == CellTopology::HEX8  && !tConformal)
            {
                return {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8};
            }
            else if(tBGCellTopo == CellTopology::TET4  && tConformal)
            {
                return {Subdivision_Method::C_HIERARCHY_TET4};
            }
        }
        else
        {
            MORIS_ASSERT(0,"Invalid spatial dimension");
        }

        MORIS_ERROR(0,"Failed determining subdivision methods");

        return Cell<enum Subdivision_Method>(0);
    }

    // ----------------------------------------------------------------------------------
    // Decomposition Source code
    // ----------------------------------------------------------------------------------

    bool
    Model::decompose(Cell<enum Subdivision_Method> aMethods)
    {
        Tracer tTracer( "XTK", "Decomposition", "Decompose" );

        // Process for a decomposition
        uint tNumDecompositions = aMethods.size();
        uint tNumGeometries     = mGeometryEngine->get_num_geometries();

        // Tell the subdivision to assign node Ids if it is the only subdivision method (critical for outputting)
        // This is usually only going to happen in test cases
        // Note: the Conformal subdivision methods dependent on node ids for subdivision routine, the node Ids are set regardless of the below boolean
        bool tNonConformingMeshFlag = false;

        if(aMethods.size() == 1)
        {
            tNonConformingMeshFlag = true;
        }

        // outer cell - geometry index, inner cell active child mesh indices for each geometries
        moris::Cell<moris::Matrix<moris::IndexMat>> tActiveChildMeshIndicesByGeom(tNumGeometries);

        // Loop over each geometry and have an active child mesh indices list for each
        for(moris::size_t iGeom = 0; iGeom<tNumGeometries; iGeom++)
        {
            bool tFirstSubdivisionFlag = true;
            moris::Matrix< moris::IndexMat > tActiveChildMeshIndices(1,1,0);

            for (moris::size_t iDecomp = 0; iDecomp < tNumDecompositions; iDecomp++)
            {
                // Perform subdivision
                this->decompose_internal(aMethods(iDecomp), iGeom, tActiveChildMeshIndices, tFirstSubdivisionFlag, tNonConformingMeshFlag);

                // Change the first subdivision flag as false
                tFirstSubdivisionFlag = false;

                if(iDecomp == 0)
                {
                    if(!this->all_child_meshes_on_same_level())
                    {
                        MORIS_LOG_INFO("Intersected cells are not on the same level. ");
                        return false;
                    }
                }

            }

            // If it's not the last geometry tell the geometry engine we're moving on
            if(iGeom!= tNumGeometries-1)
            {
                mGeometryEngine->advance_geometry_index();
            }
        }

        // Tell the xtk mesh to set all necessary information to finalize decomposition allowing
        // i.e set element ids, indices for children elements
        this->finalize_decomp();

        MORIS_LOG_SPEC("Num Intersected BG Cell",mCutMesh.get_num_child_meshes());

        return true;
    }

    // ----------------------------------------------------------------------------------

    void
    Model::decompose_internal(
            enum Subdivision_Method    const & aSubdivisionMethod,
            moris::uint                        aGeomIndex,
            moris::Matrix< moris::IndexMat > & aActiveChildMeshIndices,
            bool const &                       aFirstSubdivision,
            bool const &                       aSetIds)
    {
        // FIXME" Keenan the code below needs to be modularized; each subdevision should be its own routine
        switch (aSubdivisionMethod)
        {
            case Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8:
            {
                Tracer tTracer( "XTK", "Decomposition", "DecomposeRegularHex8" );

                this->decompose_internal_reg_sub_hex8(aGeomIndex, aActiveChildMeshIndices, aFirstSubdivision, aSetIds);
                break;
            }
            case Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4:
            {
                Tracer tTracer( "XTK", "Decomposition", "DecomposeRegularQuad4" );

                this->decompose_internal_reg_sub_quad4(aGeomIndex, aActiveChildMeshIndices, aFirstSubdivision, aSetIds);
                break;

            }
            case Subdivision_Method::C_HIERARCHY_TET4:
            {
                Tracer tTracer( "XTK", "Decomposition", "DecomposeHierarchyTet4" );

                // If it the first subdivision we need to find the intersected before placing the conformal nodes
                // Intersected elements are flagged via the Geometry_Engine
                if(aFirstSubdivision)
                {
                    moris::Matrix< moris::IndexMat > tNewPairBool;
                    this->run_first_cut_routine(aGeomIndex, aActiveChildMeshIndices,tNewPairBool);

                    for(moris::size_t i = 0; i<aActiveChildMeshIndices.n_cols(); i++)
                    {
                        Child_Mesh & tChildMesh = mCutMesh.get_child_mesh(aActiveChildMeshIndices(0,i));
                        tChildMesh.generate_connectivities(true,true,true);
                    }

                    // set the child cell topology as tet 4s
                    mCutMesh.set_child_element_topology(CellTopology::TET4);
                }

                // For hex background meshes we have a three dimension parametric coordinate
                moris::size_t tDimParamCoord = 3;

                // For tet background meshes we have a 4-d parametric coordinate
                if(aFirstSubdivision)
                {
                    tDimParamCoord = 4;
                }

                // initialize a struct of all the data we are keeping track of in this decomposition
                // intended to reduce the clutter of function inputs etc
                Decomposition_Data tDecompData;
                tDecompData.mSubdivisionMethod      = Subdivision_Method::C_HIERARCHY_TET4;
                tDecompData.mConformalDecomp        = true;
                tDecompData.mHasSecondaryIdentifier = true;
                tDecompData.mFirstSubdivision       = aFirstSubdivision;

                // Initialize topologies used in this method (all local coordinates are with respect to an edge)
                Edge_Topology tEdgeTopology;

                // initialize a couple of commonly used matrices in this method
                moris::Matrix< moris::DDRMat > tLocalCoordRelativeToEdge(1,1, 0); // ALong an edge
                moris::Matrix< moris::DDRMat > tGlobalCoord(1,3, 0); // ALong an edge
                moris::Matrix< moris::DDRMat > tEdgeNodeParamCoordinates(2,tDimParamCoord); // parametric coordinate of end nodes wrt parent element

                // Check type specified as conformal (could change this to enum)
                moris::Matrix< moris::DDRMat > tNodeCoords = mBackgroundMesh.get_all_node_coordinates_loc_inds();

                // get the underlying background mesh data
                moris::mtk::Mesh & tMeshData = mBackgroundMesh.get_mesh_data();

                // resize child mesh to new node information
                tDecompData.tCMNewNodeLoc.resize(aActiveChildMeshIndices.n_cols());
                tDecompData.tCMNewNodeParamCoord.resize(aActiveChildMeshIndices.n_cols());

                for (moris::size_t j = 0; j < aActiveChildMeshIndices.n_cols(); j++)
                {
                    // Get the child mesh that is active
                    Child_Mesh & tChildMesh = mCutMesh.get_child_mesh(aActiveChildMeshIndices(0,j));

                    // edge to node connectivity from child mesh
                    moris::Matrix< moris::IndexMat > const & tEdgeToNode = tChildMesh.get_edge_to_node();

                    // Initialize number of new nodes
                    uint tNumNewNodes = 0;

                    // get reference to child mesh edge parent information
                    moris::Matrix< moris::IndexMat > const & tEdgeParentIndices = tChildMesh.get_edge_parent_inds();
                    moris::Matrix< moris::DDSTMat >  const & tEdgeParentRanks   = tChildMesh.get_edge_parent_ranks();

                    Cell<mtk::Vertex*> tBackgroundNodes = mBackgroundMesh.get_mtk_cell(tChildMesh.get_parent_element_index()).get_vertex_pointers();

                    for (moris::size_t tEdgeInd = 0; tEdgeInd < tEdgeToNode.n_rows(); tEdgeInd++)
                    {
                        Matrix<DDUMat> tElementNodeIndices(tBackgroundNodes.size(), 1);
                        Cell<Matrix<DDRMat>> tElementNodeCoordinates(tBackgroundNodes.size());
                        for (uint tNode = 0; tNode < tBackgroundNodes.size(); tNode++)
                        {
                            tElementNodeIndices(tNode) = tBackgroundNodes(tNode)->get_index();
                            tElementNodeCoordinates(tNode) = tBackgroundNodes(tNode)->get_coords();
                        }

                        if (mGeometryEngine->queue_intersection(
                                tEdgeToNode(tEdgeInd, 0),
                                tEdgeToNode(tEdgeInd, 1),
                                tChildMesh.get_parametric_coordinates(tEdgeToNode(tEdgeInd, 0)),
                                tChildMesh.get_parametric_coordinates(tEdgeToNode(tEdgeInd, 1)),
                                tNodeCoords.get_row(tEdgeToNode(tEdgeInd, 0)),
                                tNodeCoords.get_row(tEdgeToNode(tEdgeInd, 1)),
                                tElementNodeIndices,
                                tElementNodeCoordinates))
                        {
                            // Determine which parent nodes, if any, are on the interface
                            bool tFirstParentOnInterface = mGeometryEngine->queued_intersection_first_parent_on_interface();
                            bool tSecondParentOnInterface = mGeometryEngine->queued_intersection_second_parent_on_interface();

                            if (tFirstParentOnInterface)
                            {
                                // Tell the xtk mesh that the first node is an interface node
                                mBackgroundMesh.mark_node_as_interface_node(tEdgeToNode(tEdgeInd, 0), mGeometryEngine->get_active_geometry_index());
                            }
                            if (tSecondParentOnInterface)
                            {
                                // Tell the xtk mesh that the second node is an interface node
                                mBackgroundMesh.mark_node_as_interface_node(tEdgeToNode(tEdgeInd, 1), mGeometryEngine->get_active_geometry_index());
                            }
                            if (tFirstParentOnInterface and tSecondParentOnInterface)
                            {
                                // Tell the child mesh this edge is actually on the interface already
                                tChildMesh.mark_edge_as_on_interface(tEdgeInd);
                            }
                            if (not (tFirstParentOnInterface or tSecondParentOnInterface))
                            {
                                // get a local coordinate along the intersected edge [-1,1]
                                tLocalCoordRelativeToEdge(0,0) = mGeometryEngine->get_queued_intersection_local_coordinate();

                                // get the interpolated global coordinate
                                tGlobalCoord = mGeometryEngine->get_queued_intersection_global_coordinates();

                                // Add edge to the entity intersection connectivity
                                mCutMesh.add_entity_to_intersect_connectivity(aActiveChildMeshIndices(0,j), tNumNewNodes, tEdgeInd, 0);

                                // Edge nodes
                                moris::Matrix<moris::IndexMat> tEdgeNodes = tEdgeToNode.get_row(tEdgeInd);

                                // Compute new node parametric coordinate with respect to the current parent element
                                tEdgeNodeParamCoordinates.set_row(0, tChildMesh.get_parametric_coordinates(tEdgeNodes(0)));
                                tEdgeNodeParamCoordinates.set_row(1, tChildMesh.get_parametric_coordinates(tEdgeNodes(1)));
                                moris::Matrix< moris::DDRMat > tParametricCoordsRelativeToParentElem =
                                        Interpolation::linear_interpolation_location(tEdgeNodeParamCoordinates,tLocalCoordRelativeToEdge);

                                // Parent edge information
                                moris::size_t      tParentRank  = tEdgeParentRanks(0, tEdgeInd);
                                moris::moris_index tParentIndex = tEdgeParentIndices(0, tEdgeInd);

                                // get the owning processor for an entity
                                moris::moris_index tOwningProc = tMeshData.get_entity_owner(tParentIndex, (enum EntityRank)tParentRank);

                                // Convert to global id using mesh
                                tEdgeNodes(0, 0) = mBackgroundMesh.get_glb_entity_id_from_entity_loc_index(tEdgeNodes(0, 0), EntityRank::NODE);
                                tEdgeNodes(0, 1) = mBackgroundMesh.get_glb_entity_id_from_entity_loc_index(tEdgeNodes(0, 1), EntityRank::NODE);

                                // Order the nodes in ascending order
                                if(tEdgeNodes(0, 1) < tEdgeNodes(0, 0))
                                {
                                    moris::size_t tSwap = tEdgeNodes(0, 0);
                                    tEdgeNodes(0, 0) = tEdgeNodes(0, 1);
                                    tEdgeNodes(0, 1) = tSwap;
                                }

                                // Intersected edge is an existing  edge
                                // Make request in edge requests
                                // This does not require a supplemental identifier
                                // TODO: ADD OVERFLOW CHECK IN CANTOR PAIRING
                                moris::moris_index tSecondaryId = xtk::cantor_pairing(tEdgeNodes(0, 0),tEdgeNodes(0, 1));
                                moris_index tNewNodeIndexInSubdivision = MORIS_INDEX_MAX;
                                bool tRequestExist = tDecompData.request_exists(tParentIndex,tSecondaryId,(enum EntityRank)tParentRank,tNewNodeIndexInSubdivision);

                                // location for this face in the map
                                if(!tRequestExist)
                                {
                                    // Register new request
                                    tNewNodeIndexInSubdivision = tDecompData.register_new_request(
                                            tParentIndex,
                                            tSecondaryId,
                                            tOwningProc,
                                            (enum EntityRank)tParentRank,
                                            tGlobalCoord,
                                            new Edge_Topology(tEdgeToNode.get_row(tEdgeInd)), /*this is deleted in the decomp data deconstructor*/
                                            tLocalCoordRelativeToEdge.get_row(0));

                                    uint tNewNodeIndex = mBackgroundMesh.get_first_available_index(EntityRank::NODE);
                                    tDecompData.tNewNodeIndex(tDecompData.tNewNodeIndex.size() - 1) = tNewNodeIndex;
                                    mBackgroundMesh.update_first_available_index(tNewNodeIndex + 1, EntityRank::NODE);

                                    // Admit queued node in geometry engine
                                    mGeometryEngine->admit_queued_intersection(tNewNodeIndex);
                                }

                                // add to pending node pointers for child mesh
                                tDecompData.tCMNewNodeLoc(j).push_back(tNewNodeIndexInSubdivision);

                                // add parametric coordinate to decomp data
                                tDecompData.tCMNewNodeParamCoord(j).push_back(tParametricCoordsRelativeToParentElem);

                                // Creating a new node add 1 to count
                                tNumNewNodes++;
                            }
                        }
                    }

                    tChildMesh.mark_interface_faces_from_interface_coincident_faces();
                } // XTK Mesh loop

                moris_index tMessageTag = 60001; /*arbitrary tag for regular subdivision*/
                assign_node_requests_identifiers(tDecompData,tMessageTag);

                // Allocate interface flag space in XTK mesh even though these are not interface nodes
                mBackgroundMesh.allocate_space_in_interface_node_flags(tDecompData.tNewNodeIndex.size(),mGeometryEngine->get_num_geometries());

                // add nodes to the background mesh
                mBackgroundMesh.batch_create_new_nodes(tDecompData.tNewNodeId,tDecompData.tNewNodeIndex,tDecompData.tNewNodeOwner,tDecompData.tNewNodeCoordinate);

                //update underlying ids and owners of interpolation nodes in GE
                for( uint Ik = 0; Ik < tDecompData.tNewNodeIndex.size(); Ik++)
                {
                    moris_index tNodeIndex = tDecompData.tNewNodeIndex( Ik );
                    moris_id tNodeId       = tDecompData.tNewNodeId( Ik );
                    moris_index tNodeOwner = tDecompData.tNewNodeOwner( Ik );

                    mGeometryEngine->update_queued_intersection(tNodeIndex, tNodeId, tNodeOwner);
                }

                // add nodes to child mesh
                this->decompose_internal_set_new_nodes_in_child_mesh_nh(aActiveChildMeshIndices,tDecompData);

                // mark nodes as interface nodes
                moris_index tGeomIndex = mGeometryEngine->get_active_geometry_index();
                for(moris::uint i = 0; i <tDecompData.tNewNodeId.size(); i++)
                {
                    mBackgroundMesh.mark_node_as_interface_node(tDecompData.tNewNodeIndex(i),tGeomIndex);

                    // figure out if the new node is on any other interface by looking at the parent nodes of the edge
                    Topology * tEdgeTopo = tDecompData.tNewNodeParentTopology(i);

                    // get the vertex indices on the edge
                    moris::Matrix< moris::IndexMat > const & tVerticesOnParentEdge = tEdgeTopo->get_node_indices();

                    // make sure it is an edge with 2 vertices
                    MORIS_ASSERT(tEdgeTopo->get_topology_type() == Topology_Type::EDGE,"Parent topology needs to be an edge.");
                    MORIS_ASSERT(tVerticesOnParentEdge.numel()  == 2,"Parent edge needs to have two vertices.");

                    // determine if this vertex is on other interfaces
                    for(moris::uint j = 0; j < mGeometryEngine->get_num_geometries(); j++)
                    {

                        if(mGeometryEngine->is_interface_vertex(tVerticesOnParentEdge(0),j) and mGeometryEngine->is_interface_vertex(tVerticesOnParentEdge(1),j) )
                        {
                            mBackgroundMesh.mark_node_as_interface_node(tDecompData.tNewNodeIndex(i),j);
                        }
                    }
                }

                // Set Node Ids and tell the child mesh to update
                for (moris::size_t j = 0; j < aActiveChildMeshIndices.n_cols(); j++)
                {
                    moris::Matrix< moris::IndexMat > const & tNodeIndices = mCutMesh.get_node_indices(aActiveChildMeshIndices(0,j));
                    moris::Matrix< moris::IdMat > tNodeIds = mBackgroundMesh.get_glb_entity_id_from_entity_loc_index_range(tNodeIndices, EntityRank::NODE);

                    mCutMesh.set_node_ids(aActiveChildMeshIndices(0,j), tNodeIds);
                    mCutMesh.modify_templated_mesh(aActiveChildMeshIndices(0,j), TemplateType::HIERARCHY_TET4);
                }

                break;
            }
            case Subdivision_Method::C_TRI3:
            {
                Tracer tTracer( "XTK", "Decomposition", "DecomposeHierarchyTri3" );

                // If it the first subdivision we need to find the intersected before placing the conformal nodes
                // Intersected elements are flagged via the Geometry_Engine
                if(aFirstSubdivision)
                {

                    moris::Matrix< moris::IndexMat > tNewPairBool;
                    run_first_cut_routine( aGeomIndex, aActiveChildMeshIndices,tNewPairBool);

                    for(moris::size_t i = 0; i<aActiveChildMeshIndices.n_cols(); i++)
                    {
                        Child_Mesh & tChildMesh = mCutMesh.get_child_mesh(aActiveChildMeshIndices(0,i));
                        tChildMesh.generate_connectivities(false,true,true);
                    }

                    // set the child cell topology as tet 4s
                    mCutMesh.set_child_element_topology(CellTopology::TRI3);
                }

                // For quad background meshes we have a 2 dimension parametric coordinate
                moris::size_t tDimParamCoord = 2;

                // For tri background meshes we have a 3-d parametric coordinate
                if(aFirstSubdivision)
                {
                    tDimParamCoord = 3;
                }

                // initialize a struct of all the data we are keeping track of in this decomposition
                // intended to reduce the clutter of function inputs etc
                Decomposition_Data tDecompData;
                tDecompData.mSubdivisionMethod = Subdivision_Method::C_TRI3;
                tDecompData.mConformalDecomp = true;
                tDecompData.mHasSecondaryIdentifier = true;
                tDecompData.mFirstSubdivision = aFirstSubdivision;

                // Initialize topologies used in this method (all local coordinates are with respect to an edge)
                Edge_Topology tEdgeTopology;

                // initialize a couple of commonly used matrices in this method
                moris::Matrix< moris::DDRMat > tLocalCoordRelativeToEdge(1,1, 0); // ALong an edge
                moris::Matrix< moris::DDRMat > tGlobalCoord(1,2, 0); // ALong an edge
                moris::Matrix< moris::DDRMat > tEdgeNodeParamCoordinates(2,tDimParamCoord); // parametric coordinate of end nodes wrt parent element

                // Check type specified as conformal (could change this to enum)
                moris::Matrix< moris::DDRMat > tNodeCoords = mBackgroundMesh.get_all_node_coordinates_loc_inds();

                // get the underlying background mesh data
                moris::mtk::Mesh & tMeshData = mBackgroundMesh.get_mesh_data();

                // resize child mesh to new node information
                tDecompData.tCMNewNodeLoc.resize(aActiveChildMeshIndices.n_cols());
                tDecompData.tCMNewNodeParamCoord.resize(aActiveChildMeshIndices.n_cols());

                for (moris::size_t j = 0; j < aActiveChildMeshIndices.n_cols(); j++)
                {
                    // Get the child mesh that is active
                    Child_Mesh & tChildMesh = mCutMesh.get_child_mesh(aActiveChildMeshIndices(0,j));

                    // edge to node connectivity from child mesh
                    moris::Matrix< moris::IndexMat > const & tEdgeToNode = tChildMesh.get_edge_to_node();

                    // Initialize number of new nodes
                    uint tNumNewNodes = 0;

                    // get reference to child mesh edge parent information
                    moris::Matrix< moris::IndexMat > const & tEdgeParentIndices = tChildMesh.get_edge_parent_inds();
                    moris::Matrix< moris::DDSTMat >  const & tEdgeParentRanks   = tChildMesh.get_edge_parent_ranks();

                    Cell<mtk::Vertex*> tBackgroundNodes = mBackgroundMesh.get_mtk_cell(tChildMesh.get_parent_element_index()).get_vertex_pointers();

                    for (moris::size_t tEdgeInd = 0; tEdgeInd < tEdgeToNode.n_rows(); tEdgeInd++)
                    {
                        Matrix<DDUMat> tElementNodeIndices(tBackgroundNodes.size(), 1);
                        Cell<Matrix<DDRMat>> tElementNodeCoordinates(tBackgroundNodes.size());
                        for (uint tNode = 0; tNode < tBackgroundNodes.size(); tNode++)
                        {
                            tElementNodeIndices(tNode) = tBackgroundNodes(tNode)->get_index();
                            tElementNodeCoordinates(tNode) = tBackgroundNodes(tNode)->get_coords();
                        }

                        if (mGeometryEngine->queue_intersection(
                                tEdgeToNode(tEdgeInd, 0),
                                tEdgeToNode(tEdgeInd, 1),
                                tChildMesh.get_parametric_coordinates(tEdgeToNode(tEdgeInd, 0)),
                                tChildMesh.get_parametric_coordinates(tEdgeToNode(tEdgeInd, 1)),
                                tNodeCoords.get_row(tEdgeToNode(tEdgeInd, 0)),
                                tNodeCoords.get_row(tEdgeToNode(tEdgeInd, 1)),
                                tElementNodeIndices,
                                tElementNodeCoordinates))
                        {
                            // Determine which parent nodes, if any, are on the interface
                            bool tFirstParentOnInterface = mGeometryEngine->queued_intersection_first_parent_on_interface();
                            bool tSecondParentOnInterface = mGeometryEngine->queued_intersection_second_parent_on_interface();

                            if (tFirstParentOnInterface)
                            {
                                // Tell the xtk mesh that the first node is an interface node
                                mBackgroundMesh.mark_node_as_interface_node(tEdgeToNode(tEdgeInd, 0), mGeometryEngine->get_active_geometry_index());
                            }
                            if (tSecondParentOnInterface)
                            {
                                // Tell the xtk mesh that the second node is an interface node
                                mBackgroundMesh.mark_node_as_interface_node(tEdgeToNode(tEdgeInd, 1), mGeometryEngine->get_active_geometry_index());
                            }
                            if (tFirstParentOnInterface and tSecondParentOnInterface)
                            {
                                // Tell the child mesh this edge is actually on the interface already
                                tChildMesh.mark_edge_as_on_interface(tEdgeInd);
                            }
                            if (not (tFirstParentOnInterface or tSecondParentOnInterface))
                            {
                                // get a local coordinate along the intersected edge [-1,1]
                                tLocalCoordRelativeToEdge(0,0) = mGeometryEngine->get_queued_intersection_local_coordinate();

                                // get the interpolated global coordinate
                                tGlobalCoord = mGeometryEngine->get_queued_intersection_global_coordinates();

                                // Add edge to the entity intersection connectivity
                                mCutMesh.add_entity_to_intersect_connectivity(aActiveChildMeshIndices(0,j), tNumNewNodes, tEdgeInd, 0);

                                // Edge nodes
                                moris::Matrix<moris::IndexMat> tEdgeNodes = tEdgeToNode.get_row(tEdgeInd);

                                // Compute new node parametric coordinate with respect to the current parent element
                                tEdgeNodeParamCoordinates.set_row(0, tChildMesh.get_parametric_coordinates(tEdgeNodes(0)));
                                tEdgeNodeParamCoordinates.set_row(1, tChildMesh.get_parametric_coordinates(tEdgeNodes(1)));

                                moris::Matrix< moris::DDRMat > tParametricCoordsRelativeToParentElem =
                                        Interpolation::linear_interpolation_location(tEdgeNodeParamCoordinates,tLocalCoordRelativeToEdge);

                                // Parent edge information
                                moris::size_t      tParentRank  = tEdgeParentRanks(0, tEdgeInd);
                                moris::moris_index tParentIndex = tEdgeParentIndices(0, tEdgeInd);

                                // get the owning processor for an entity
                                moris::moris_index tOwningProc = tMeshData.get_entity_owner(tParentIndex, (enum EntityRank)tParentRank);

                                // Convert to global id using mesh
                                tEdgeNodes(0, 0) = mBackgroundMesh.get_glb_entity_id_from_entity_loc_index(tEdgeNodes(0, 0), EntityRank::NODE);
                                tEdgeNodes(0, 1) = mBackgroundMesh.get_glb_entity_id_from_entity_loc_index(tEdgeNodes(0, 1), EntityRank::NODE);

                                // Order the nodes in ascending order
                                if(tEdgeNodes(0, 1) < tEdgeNodes(0, 0))
                                {
                                    moris::size_t tSwap = tEdgeNodes(0, 0);
                                    tEdgeNodes(0, 0) = tEdgeNodes(0, 1);
                                    tEdgeNodes(0, 1) = tSwap;
                                }

                                // Intersected edge is an existing  edge
                                // Make request in edge requests
                                // This does not require a supplemental identifier
                                // TODO: ADD OVERFLOW CHECK IN CANTOR PAIRING!!!!!!
                                moris::moris_index tSecondaryId = xtk::cantor_pairing(tEdgeNodes(0, 0),tEdgeNodes(0, 1));
                                moris_index tNewNodeIndexInSubdivision = MORIS_INDEX_MAX;

                                bool tRequestExist = tDecompData.request_exists(
                                        tParentIndex,
                                        tSecondaryId,
                                        (enum EntityRank)tParentRank,
                                        tNewNodeIndexInSubdivision);

                                // location for this face in the map
                                if(!tRequestExist)
                                {
                                    // Register new request
                                    tNewNodeIndexInSubdivision = tDecompData.register_new_request(
                                            tParentIndex,
                                            tSecondaryId,
                                            tOwningProc,
                                            (enum EntityRank)tParentRank,
                                            tGlobalCoord,
                                            new Edge_Topology(tEdgeToNode.get_row(tEdgeInd)), /*Note: this is deleted in the decomp data deconstructor*/
                                            tLocalCoordRelativeToEdge.get_row(0));

                                    uint tNewNodeIndex = mBackgroundMesh.get_first_available_index(EntityRank::NODE);
                                    tDecompData.tNewNodeIndex(tDecompData.tNewNodeIndex.size() - 1) = tNewNodeIndex;
                                    mBackgroundMesh.update_first_available_index(tNewNodeIndex + 1, EntityRank::NODE);

                                    // Admit queued node in geometry engine
                                    mGeometryEngine->admit_queued_intersection(tNewNodeIndex);
                                }

                                // add to pending node pointers for child mesh
                                tDecompData.tCMNewNodeLoc(j).push_back(tNewNodeIndexInSubdivision);

                                // add parametric coordinate to decomp data
                                tDecompData.tCMNewNodeParamCoord(j).push_back(tParametricCoordsRelativeToParentElem);

                                // Creating a new node add 1 to count
                                tNumNewNodes++;
                            }
                        }
                    }

                    tChildMesh.mark_interface_faces_from_interface_coincident_faces();

                } // XTK Mesh loop

                moris_index tMessageTag = 60001; /*arbitrary tag for regular subdivision*/
                assign_node_requests_identifiers(tDecompData,tMessageTag);

                // Allocate interface flag space in XTK mesh even though these are not interface nodes
                mBackgroundMesh.allocate_space_in_interface_node_flags(tDecompData.tNewNodeIndex.size(),mGeometryEngine->get_num_geometries());

                // add nodes to the background mesh
                mBackgroundMesh.batch_create_new_nodes(tDecompData.tNewNodeId,tDecompData.tNewNodeIndex,tDecompData.tNewNodeOwner,tDecompData.tNewNodeCoordinate);

                //update underlying ids and owners of interpolation nodes in GE
                for( uint Ik = 0; Ik < tDecompData.tNewNodeIndex.size(); Ik++)
                {
                    moris_index tNodeIndex = tDecompData.tNewNodeIndex( Ik );
                    moris_id tNodeId       = tDecompData.tNewNodeId( Ik );
                    moris_index tNodeOwner = tDecompData.tNewNodeOwner( Ik );

                    mGeometryEngine->update_queued_intersection(tNodeIndex, tNodeId, tNodeOwner);
                }

                // add nodes to child mesh
                this->decompose_internal_set_new_nodes_in_child_mesh_nh(aActiveChildMeshIndices,tDecompData);

                // mark nodes as interface nodes
                moris_index tGeomIndex = mGeometryEngine->get_active_geometry_index();
                for(moris::uint i = 0; i <tDecompData.tNewNodeId.size(); i++)
                {
                    // this node is always on the geometry interface of current so mark this
                    mBackgroundMesh.mark_node_as_interface_node(tDecompData.tNewNodeIndex(i),tGeomIndex);

                    // figure out if the new node is on any other interface by looking at the parent nodes of the edge
                    Topology * tEdgeTopo = tDecompData.tNewNodeParentTopology(i);

                    // get the vertex indices on the edge
                    moris::Matrix< moris::IndexMat > const & tVerticesOnParentEdge = tEdgeTopo->get_node_indices();

                    // make sure it is an edge with 2 vertices
                    MORIS_ASSERT(tEdgeTopo->get_topology_type() == Topology_Type::EDGE,"Parent topology needs to be an edge.");
                    MORIS_ASSERT(tVerticesOnParentEdge.numel()  == 2,"Parent edge needs to have two vertices.");

                    // determine if this vertex is on other interfaces
                    for(moris::uint j = 0; j < mGeometryEngine->get_active_geometry_index(); j++)
                    {
                        if(mBackgroundMesh.is_interface_node(tVerticesOnParentEdge(0),j) and mBackgroundMesh.is_interface_node(tVerticesOnParentEdge(1),j) )
                        {
                            mBackgroundMesh.mark_node_as_interface_node(tDecompData.tNewNodeIndex(i),j);
                        }
                    }
                }

                // Set Node Ids and tell the child mesh to update
                for (moris::size_t j = 0; j < aActiveChildMeshIndices.n_cols(); j++)
                {
                    moris::Matrix< moris::IndexMat > const & tNodeIndices =
                            mCutMesh.get_node_indices(aActiveChildMeshIndices(0,j));

                    moris::Matrix< moris::IdMat > tNodeIds =
                            mBackgroundMesh.get_glb_entity_id_from_entity_loc_index_range(tNodeIndices, EntityRank::NODE);

                    mCutMesh.set_node_ids(aActiveChildMeshIndices(0,j), tNodeIds);
                    mCutMesh.modify_templated_mesh(aActiveChildMeshIndices(0,j), TemplateType::CONFORMAL_TRI3);
                }

                break;
            }
            default:
            {
                moris::size_t breaker = 0;
                MORIS_ERROR(breaker != 0, "formulate_node_request should not enter the default case, check to see if your aCheckType is undefined.");
            }
        }

    }

    // ----------------------------------------------------------------------------------

    void
    Model::decompose_internal_reg_sub_hex8(
            moris::uint                        aGeomIndex,
            moris::Matrix< moris::IndexMat > & aActiveChildMeshIndices,
            bool                       const & aFirstSubdivision,
            bool                       const & aSetIds )
    {
        MORIS_ASSERT(aFirstSubdivision,"NC_REGULAR_SUBDIVISION_HEX8 needs to be the first subdivision routine for each geometry");
        MORIS_ASSERT(mModelDimension == 3,"NC_REGULAR_SUBDIVISION_HEX8 needs to be done on a 3D mesh");

        // Runs the first cut routine to get the new active child mesh indices and indicate which are new and need to be regularly subdivided and which ones dont
        moris::Matrix< moris::IndexMat > tNewPairBool;
        run_first_cut_routine(aGeomIndex, aActiveChildMeshIndices, tNewPairBool);

        // set the child cell topology as tet 3s
        mCutMesh.set_child_element_topology(CellTopology::TET4);

        // initialize a struct of all the data we are keeping track of in this decomposition
        // intended to reduce the clutter of function inputs etc
        Decomposition_Data tDecompData;
        tDecompData.mSubdivisionMethod = Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8;

        // number of intersected elements
        moris::uint tIntersectedCount = aActiveChildMeshIndices.n_cols();

        // make node requests for each intersected element
        this->decompose_internal_reg_sub_hex8_make_requests(aActiveChildMeshIndices,tNewPairBool,tDecompData);

        // specify a dummy secondary id (not really needed for this type of decomposition)
        tDecompData.tSecondaryIdentifiers = Cell<moris_index>(tDecompData.tNewNodeParentIndex.size(), MORIS_INDEX_MAX);

        moris_index tMessageTag = 60000; /*arbitrary tag for regular subdivision*/
        assign_node_requests_identifiers(tDecompData,tMessageTag);

        // Allocate interface flag space in XTK mesh even though these are not interface nodes
        mBackgroundMesh.allocate_space_in_interface_node_flags(tDecompData.tNewNodeIndex.size(),mGeometryEngine->get_num_geometries());

        // add nodes to the background mesh
        this->assign_index(tDecompData);
        mBackgroundMesh.batch_create_new_nodes(tDecompData.tNewNodeId,tDecompData.tNewNodeIndex, tDecompData.tNewNodeOwner,tDecompData.tNewNodeCoordinate);

        // add nodes to child mesh
        this->decompose_internal_set_new_nodes_in_child_mesh_reg_sub(aActiveChildMeshIndices,tNewPairBool, 3, tDecompData);

        // associate new nodes with geometry objects
        this->create_new_node_association_with_geometry(tDecompData);

        for(moris::size_t i = 0; i< tIntersectedCount; i++)
        {
            if(tNewPairBool(0,i) == 0)
            {
                mCutMesh.generate_templated_mesh(aActiveChildMeshIndices(i),TemplateType::REGULAR_SUBDIVISION_HEX8);
            }
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Model::decompose_internal_reg_sub_hex8_make_requests(
            moris::Matrix< moris::IndexMat > & aActiveChildMeshIndices,
            moris::Matrix< moris::IndexMat > & tNewPairBool,
            Decomposition_Data               & tDecompData)
    {
        // mesh data accessor
        moris::mtk::Mesh & tXTKMeshData = mBackgroundMesh.get_mesh_data();

        // number of intersected elements
        moris::uint tIntersectedCount = aActiveChildMeshIndices.n_cols();

        // allocate child mesh to new node location
        tDecompData.tCMNewNodeLoc.resize(tIntersectedCount,7);
        tDecompData.tCMNewNodeParamCoord.resize(tIntersectedCount);

        // parametric coordinates relative to hex where we put the nodes
        // Parametric coordinates for this subdivision routine
        const moris::Matrix< moris::DDRMat > tParamCoordsRelativeToElem({
            { 0.0, -1.0,  0.0},
            { 1.0,  0.0,  0.0},
            { 0.0,  1.0,  0.0},
            {-1.0,  0.0,  0.0},
            { 0.0,  0.0, -1.0},
            { 0.0,  0.0,  1.0},
            { 0.0,  0.0,  0.0}});

        const moris::Matrix< moris::DDRMat > tParamCoordsRelativeToFace({
            { 0.0,  0.0},
            { 0.0,  0.0},
            { 0.0,  0.0},
            { 0.0,  0.0},
            { 0.0,  0.0},
            { 0.0,  0.0}});

        // get the underlying background mesh data
        moris::mtk::Mesh & tMeshData = mBackgroundMesh.get_mesh_data();

        // setup Child mesh to new node location
        for (moris::size_t i = 0; i < tIntersectedCount; i++)
        {
            if(tNewPairBool(0,i) == 0)
            {
                // Get element index
                moris::moris_index tElemInd = mCutMesh.get_parent_element_index(aActiveChildMeshIndices(0,i));

                // Get local index of faces connected to element using local element index
                moris::Matrix<moris::IndexMat> tFaceIndices =
                        tXTKMeshData.get_entity_connected_to_entity_loc_inds(
                                tElemInd,
                                moris::EntityRank::ELEMENT,
                                moris::EntityRank::FACE);

                // Loop over faces (6 in a hex 8) and set a node request.
                // Request will return a pointer to where the created node index will be placed
                for (moris::size_t fi = 0; fi < 6; fi++)
                {
                    moris_index tRequestLoc = MORIS_INDEX_MAX;
                    bool tRequestExists = tDecompData.request_exists(tFaceIndices(fi),EntityRank::FACE,tRequestLoc);

                    // if we haven't created a node on this face then create one
                    if(!tRequestExists)
                    {
                        // node indices attached to face fi
                        moris::Matrix<moris::IndexMat> tFaceNodes = tXTKMeshData.get_entity_connected_to_entity_loc_inds(
                                tFaceIndices(fi),
                                moris::EntityRank::FACE,
                                moris::EntityRank::NODE);

                        // face owner
                        moris::moris_index tOwningProc = tMeshData.get_entity_owner(tFaceIndices(fi), EntityRank::FACE);

                        // coordinates of nodes attached to the nodes of this face
                        moris::Matrix<moris::DDRMat> tCoordinates =
                                mBackgroundMesh.get_selected_node_coordinates_loc_inds(tFaceNodes);

                        // bilinearly interpolate to the center of this face fi
                        moris::Matrix<moris::DDRMat> tNewNodeCoordinates;
                        xtk::Interpolation::bilinear_interpolation(
                                tCoordinates,
                                tParamCoordsRelativeToFace.get_row(fi),
                                tNewNodeCoordinates);

                        // location for this face in the map
                        moris_index tNewNodeIndexInSubdivision = tDecompData.register_new_request(
                                tFaceIndices(fi),
                                tOwningProc,
                                EntityRank::FACE,
                                tNewNodeCoordinates,
                                new Quad_4_Topology(tFaceNodes), /*Note: this is deleted in the decomp data deconstructor*/
                                tParamCoordsRelativeToFace.get_row(fi));

                        // add to pending node pointers for child mesh
                        tDecompData.tCMNewNodeLoc(i)(fi) = tNewNodeIndexInSubdivision;

                        // add parametric coordinate to decomp data
                        tDecompData.tCMNewNodeParamCoord(i).push_back(tParamCoordsRelativeToElem.get_row(fi));
                    }

                    // if debug check the coordinate will be the same
                    else
                    {
                        tDecompData.tCMNewNodeLoc(i)(fi) = tRequestLoc;
                        tDecompData.tCMNewNodeParamCoord(i).push_back(tParamCoordsRelativeToElem.get_row(fi));
#ifdef DEBUG
                        moris::uint tNewNodeIndexInSubdivision = tRequestLoc;

                        // node indices attached to face fi
                        moris::Matrix<moris::IndexMat> tFaceNodes =
                                tXTKMeshData.get_entity_connected_to_entity_loc_inds(
                                        tFaceIndices(fi),
                                        moris::EntityRank::FACE,
                                        moris::EntityRank::NODE);

                        // coordinates of nodes attached to the nodes of this face
                        moris::Matrix<moris::DDRMat> tCoordinates = mBackgroundMesh.get_selected_node_coordinates_loc_inds(tFaceNodes);

                        // bilinearly interpolate to the center of this face fi
                        moris::Matrix<moris::DDRMat> tNewNodeCoordinates;
                        xtk::Interpolation::bilinear_interpolation(tCoordinates, tParamCoordsRelativeToFace.get_row(fi),tNewNodeCoordinates);

                        // other coordinate
                        moris::Matrix<moris::DDRMat> tExistingNodeCoordinate = tDecompData.tNewNodeCoordinate(tNewNodeIndexInSubdivision);

                        MORIS_ASSERT(all_true(tNewNodeCoordinates == tExistingNodeCoordinate) ,"Node coordinates created on same face do not match");
#endif
                    }
                }

                // Place node at center of element
                // get the nodes attached to the element
                moris::Matrix<moris::IndexMat>tElementNodes =
                        tXTKMeshData.get_entity_connected_to_entity_loc_inds(
                                tElemInd,
                                moris::EntityRank::ELEMENT,
                                moris::EntityRank::NODE);

                // coordinates of nodes attached to element
                moris::Matrix<moris::DDRMat> tCoordinates = mBackgroundMesh.get_selected_node_coordinates_loc_inds(tElementNodes);

                // trilinearly interpolate to the center of the element
                moris::Matrix<moris::DDRMat> tNewNodeCoordinates;
                xtk::Interpolation::trilinear_interpolation(tCoordinates, tParamCoordsRelativeToElem.get_row(6), tNewNodeCoordinates);

                // add the new node at center of element to the map
                // location for this face in the map
                moris_index tNewNodeIndexInSubdivision = MORIS_INDEX_MAX;

                // owner of element
                moris::moris_index tOwningProc = tMeshData.get_entity_owner(tElemInd, EntityRank::ELEMENT);

                MORIS_ASSERT(!tDecompData.request_exists(tElemInd,EntityRank::ELEMENT,tNewNodeIndexInSubdivision),
                        "All element requests should be unique, therefore tNewRequest is expected to be true here");

                tNewNodeIndexInSubdivision = tDecompData.register_new_request(
                        tElemInd,
                        tOwningProc,
                        EntityRank::ELEMENT,
                        tNewNodeCoordinates,
                        new Hexahedron_8_Topology(tElementNodes), /*Note: this is deleted in the decomp data deconstructor*/
                        tParamCoordsRelativeToElem.get_row(6));

                // add child mesh new node location and parametric coordinate relative to element
                tDecompData.tCMNewNodeLoc(i)(6) = tNewNodeIndexInSubdivision;
                tDecompData.tCMNewNodeParamCoord(i).push_back(tParamCoordsRelativeToElem.get_row(6));
            }
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Model::decompose_internal_reg_sub_quad4(
            moris::uint                        aGeomIndex,
            moris::Matrix< moris::IndexMat > & aActiveChildMeshIndices,
            bool                       const & aFirstSubdivision,
            bool                       const & aSetIds )
    {
        MORIS_ASSERT(aFirstSubdivision, "NC_REGULAR_SUBDIVISION_QUAD4 needs to be the first subdivision routine for each geometry.");
        MORIS_ASSERT(mModelDimension == 2, "NC_REGULAR_SUBDIVISION_QUAD4 needs to be done on a 2D mesh.");


        // Runs the first cut routine to get the new active child mesh indices and indicate which are new and need to be regularly subdivided and which ones don't
        moris::Matrix< moris::IndexMat > tNewPairBool;
        run_first_cut_routine(aGeomIndex, aActiveChildMeshIndices, tNewPairBool);

        // mark child cells as tri 3s
        mCutMesh.set_child_element_topology(CellTopology::TRI3);

        // initialize a struct of all the data we are keeping track of in this decomposition
        // intended to reduce the clutter of function inputs etc
        Decomposition_Data tDecompData;
        tDecompData.mSubdivisionMethod = Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4;

        // number of intersected elements
        moris::uint tIntersectedCount = aActiveChildMeshIndices.n_cols();

        // make node requests for each intersected element
        this->decompose_internal_reg_sub_quad4_make_requests(aActiveChildMeshIndices, tNewPairBool, tDecompData);

        // specify a dummy secondary id (not really needed for this type of decomposition)
        tDecompData.tSecondaryIdentifiers = Cell<moris_index>(tDecompData.tNewNodeParentIndex.size(), MORIS_INDEX_MAX);

        moris_index tMessageTag = 60000; /*arbitrary tag for regular subdivision*/
        assign_node_requests_identifiers(tDecompData,tMessageTag);

        // Allocate interface flag space in XTK mesh even though these are not interface nodes
        mBackgroundMesh.allocate_space_in_interface_node_flags(tDecompData.tNewNodeIndex.size(),mGeometryEngine->get_num_geometries());

        // add nodes to the background mesh
        this->assign_index(tDecompData);
        mBackgroundMesh.batch_create_new_nodes(tDecompData.tNewNodeId,tDecompData.tNewNodeIndex, tDecompData.tNewNodeOwner,tDecompData.tNewNodeCoordinate);

        // crate nodes in child mesh
        this->decompose_internal_set_new_nodes_in_child_mesh_reg_sub(aActiveChildMeshIndices,tNewPairBool,2,tDecompData);

        // associate new nodes with geometry objects
        create_new_node_association_with_geometry(tDecompData);

        for(moris::size_t i = 0; i< tIntersectedCount; i++)
        {
            if(tNewPairBool(0,i) == 0)
            {
                mCutMesh.generate_templated_mesh(aActiveChildMeshIndices(i),TemplateType::REGULAR_SUBDIVISION_QUAD4);
            }
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Model::decompose_internal_reg_sub_quad4_make_requests(
            moris::Matrix< moris::IndexMat > & aActiveChildMeshIndices,
            moris::Matrix< moris::IndexMat > & tNewPairBool,
            Decomposition_Data               & tDecompData)
    {
        // Get access to mesh data
        moris::mtk::Mesh & tXTKMeshData = mBackgroundMesh.get_mesh_data();

        // Get number of intersected elements
        moris::uint tIntersectedCount = aActiveChildMeshIndices.n_cols();

        // Allocate child mesh to new node location
        tDecompData.tCMNewNodeLoc.resize(tIntersectedCount,1);
        tDecompData.tCMNewNodeParamCoord.resize(tIntersectedCount);

        // Parametric coordinates relative to quad where we put the nodes
        const moris::Matrix< moris::DDRMat > tParamCoordsRelativeToElem({{ 0.0,  0.0}});

        // get the underlying background mesh data
        moris::mtk::Mesh & tMeshData = mBackgroundMesh.get_mesh_data();

        // setup Child mesh to new node location
        for (moris::size_t i = 0; i < tIntersectedCount; i++)
        {
            if(tNewPairBool(0,i) == 0)
            {
                // Get element index
                moris::moris_index tElemInd = mCutMesh.get_parent_element_index(aActiveChildMeshIndices(0,i));

                // Place node at center of element
                // get the nodes attached to the element
                moris::Matrix<moris::IndexMat>tElementNodes =
                        tXTKMeshData.get_entity_connected_to_entity_loc_inds(
                                tElemInd,
                                moris::EntityRank::ELEMENT,
                                moris::EntityRank::NODE);

                // coordinates of nodes attached to element
                moris::Matrix<moris::DDRMat> tCoordinates =
                        mBackgroundMesh.get_selected_node_coordinates_loc_inds(tElementNodes);

                // trilinearly interpolate to the center of the element
                moris::Matrix<moris::DDRMat> tNewNodeCoordinates;
                xtk::Interpolation::bilinear_interpolation(tCoordinates, tParamCoordsRelativeToElem.get_row(0), tNewNodeCoordinates);

                // add the new node at center of element to the map
                // location for this face in the map
                moris_index tNewNodeIndexInSubdivision = MORIS_INDEX_MAX;

                // owner of element
                moris::moris_index tOwningProc = tMeshData.get_entity_owner(tElemInd, EntityRank::ELEMENT);

                MORIS_ASSERT(!tDecompData.request_exists(tElemInd,EntityRank::ELEMENT,tNewNodeIndexInSubdivision),
                        "All element requests should be unique, therefore tNewRequest is expected to be true here");

                tNewNodeIndexInSubdivision = tDecompData.register_new_request(
                        tElemInd,
                        tOwningProc,
                        EntityRank::ELEMENT,
                        tNewNodeCoordinates,
                        new Quad_4_Topology(tElementNodes), /*Note: this is deleted in the decomp data deconstructor*/
                        tParamCoordsRelativeToElem.get_row(0));

                // add child mesh new node location and parametric coordinate relative to element
                tDecompData.tCMNewNodeLoc(i)(0) = tNewNodeIndexInSubdivision;
                tDecompData.tCMNewNodeParamCoord(i).push_back(tParamCoordsRelativeToElem.get_row(0));
            }
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Model::decompose_internal_set_new_nodes_in_child_mesh_reg_sub(
            moris::Matrix< moris::IndexMat > & aActiveChildMeshIndices,
            moris::Matrix< moris::IndexMat > & tNewPairBool,
            moris::real                        tNumParamCoords,
            Decomposition_Data &               tDecompData)
    {
        // number of intersected elements
        moris::uint tIntersectedCount = aActiveChildMeshIndices.n_cols();

        // iterate through active child mesh indices
        for(moris::uint i = 0 ; i <tIntersectedCount; i++)
        {
            // only regularly subdivide if it hasnt already been regularly subdivided
            if(tNewPairBool(0,i) == 0)
            {
                // number of new nodes for child mesh i
                moris::uint tNumNewNodesForCM = tDecompData.tCMNewNodeLoc(i).size();

                // matrix of new node indices
                moris::Matrix<IndexMat> tCMNewNodeInds(1,tNumNewNodesForCM);

                // matrix of new node ids
                moris::Matrix<IdMat> tCMNewNodeIds(1,tNumNewNodesForCM);

                // iterate through new nodes for child mesh i to collect index and id
                for(moris::uint iN =0; iN< tNumNewNodesForCM; iN++)
                {
                    // location relative to the decomposition data
                    moris::moris_index tNodeIndexInRequestVect = tDecompData.tCMNewNodeLoc(i)(iN);

                    // retreive node index and id
                    tCMNewNodeInds(iN) = tDecompData.tNewNodeIndex(tNodeIndexInRequestVect);
                    tCMNewNodeIds(iN)  = tDecompData.tNewNodeId(tNodeIndexInRequestVect);
                }

                // retrieve child mesh
                Child_Mesh & tChildMesh = mCutMesh.get_child_mesh(aActiveChildMeshIndices(i));

                // add node indices, ids, and vertices to child mesh
                tChildMesh.add_node_indices(tCMNewNodeInds);
                tChildMesh.add_node_ids(tCMNewNodeIds);

                // allocate space for parametric coordinates
                tChildMesh.allocate_parametric_coordinates(tNumNewNodesForCM,tNumParamCoords);

                // iterate through nods and add parametric coordinate
                for(moris::uint iN =0; iN< tNumNewNodesForCM; iN++)
                {
                    tChildMesh.add_node_parametric_coordinate( tCMNewNodeInds(iN),tDecompData.tCMNewNodeParamCoord(i)(iN));
                }
            }
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Model::decompose_internal_set_new_nodes_in_child_mesh_nh(
            moris::Matrix< moris::IndexMat > & aActiveChildMeshIndices,
            Decomposition_Data &               tDecompData)
    {
        // number of intersected elements
        moris::uint tIntersectedCount = aActiveChildMeshIndices.n_cols();

        // iterate through active child mesh indices
        for(moris::uint i = 0 ; i <tIntersectedCount; i++)
        {
            // number of new nodes for child mesh i
            moris::uint tNumNewNodesForCM = tDecompData.tCMNewNodeLoc(i).size();

            // matrix of new node indices
            moris::Matrix<IndexMat> tCMNewNodeInds(1,tNumNewNodesForCM);

            // matrix of new node ids
            moris::Matrix<IdMat> tCMNewNodeIds(1,tNumNewNodesForCM);

            // iterate through new nodes for child mesh i to collect index and id
            for(moris::uint iN =0; iN< tNumNewNodesForCM; iN++)
            {
                // location relative to the decomposition data
                moris::moris_index tNodeIndexInRequestVect = tDecompData.tCMNewNodeLoc(i)(iN);

                // retreive node index and id
                tCMNewNodeInds(iN) = tDecompData.tNewNodeIndex(tNodeIndexInRequestVect);
                tCMNewNodeIds(iN)  = tDecompData.tNewNodeId(tNodeIndexInRequestVect);
            }

            // retrieve child mesh
            Child_Mesh & tChildMesh = mCutMesh.get_child_mesh(aActiveChildMeshIndices(i));

            // add node indices, ids, and vertices to child mesh
            tChildMesh.add_node_indices(tCMNewNodeInds);
            tChildMesh.add_node_ids(tCMNewNodeIds);

            // allocate space for parametric coordinate

            // For hex background meshes we have a three dimension parametric coordinate
            enum CellTopology tBackgroundTopo = mBackgroundMesh.get_parent_cell_topology();

            moris::size_t tDimParamCoord = 0;

            // FIXME: Keenan - should be converted into switch statement if possible

            if(tBackgroundTopo == CellTopology::HEX8)
            {
                tDimParamCoord =3;
            }
            else if(tBackgroundTopo == CellTopology::TET4)
            {
                tDimParamCoord =4;
            }
            else if(tBackgroundTopo == CellTopology::QUAD4)
            {
                tDimParamCoord = 2;
            }
            else if(tBackgroundTopo == CellTopology::TRI3)
            {
                tDimParamCoord = 3;
            }
            else
            {
                MORIS_ERROR(0,"Invalid background cell topo");
            }

            tChildMesh.allocate_parametric_coordinates(tNumNewNodesForCM,tDimParamCoord);

            // iterate through nods and add parametric coordinate
            for(moris::uint iN =0; iN< tNumNewNodesForCM; iN++)
            {
                tChildMesh.add_node_parametric_coordinate( tCMNewNodeInds(iN),tDecompData.tCMNewNodeParamCoord(i)(iN));
            }
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Model::create_new_node_association_with_geometry(Decomposition_Data & tDecompData)
    {
        // create geometry objects for each node
        mGeometryEngine->create_new_child_nodes(
                tDecompData.tNewNodeIndex,
                tDecompData.tNewNodeParentTopology,
                tDecompData.tParamCoordRelativeToParent,
                mBackgroundMesh.get_all_node_coordinates_loc_inds());
    }

    // ----------------------------------------------------------------------------------

    void
    Model::catch_all_unhandled_interfaces()
    {

        MORIS_ERROR(this->check_for_all_cell_vertices_on_interface(),"All vertices of a cell on the interface");

        // MORIS_ERROR(this->check_for_degenerated_cells(),"Degenerated Cells Detected");

        // iterate through child meshes
        for(moris::uint iCM = 0; iCM < mCutMesh.get_num_child_meshes(); iCM++)
        {
            // active child mesh
            Child_Mesh & tChildMesh = mCutMesh.get_child_mesh(iCM);

            // child mesh vertex indices
            moris::Matrix< moris::IndexMat > const & tVertexIndices = tChildMesh.get_node_indices();

            // Keep track of which geomtries have interface vertices
            moris::Matrix< moris::IndexMat > tGeometryInterfaceBool(1, mGeometryEngine->get_num_geometries(), (moris_index)false);

            // iterate through the vertices of the child mesh and figure out which are interface nodes
            // this loop is so I don't have to loop over all facets
            for(moris::uint iV = 0; iV < tVertexIndices.numel(); iV++)
            {
                // iterate through geometries
                for(moris::uint iG = 0; iG < mGeometryEngine->get_num_geometries(); iG++)
                {
                    if(mGeometryEngine->is_interface_vertex(tVertexIndices(iV),(moris_index)iG))
                    {
                        // mark as interface relative to the geometry
                        tGeometryInterfaceBool(iG) = (moris_index) true;
                        break;
                    }
                }
            }

            // get the facet to node
            moris::Matrix< moris::IndexMat > const & tFacetToNode = tChildMesh.get_facet_to_node();

            // iterate through geometries, check and flag facets that are on the interface
            for(moris::uint iG = 0; iG < mGeometryEngine->get_num_geometries(); iG++)
            {
                if(tGeometryInterfaceBool(iG) == (moris_index)true)
                {
                    // iterate through the facets
                    for(moris::uint iFacet = 0;  iFacet < tFacetToNode.n_rows(); iFacet++)
                    {
                        // if one is not on the interface, this will be false and the facet is not on the interface
                        bool tIsInterfaceFacet = true;

                        // iterate through nodes on facet
                        for(moris::uint iV = 0; iV < tFacetToNode.n_cols(); iV++)
                        {
                            // if we aren't on the interface, then flip the flag and move onto the next facet in child mesh
                            if(! mGeometryEngine->is_interface_vertex(tFacetToNode(iFacet,iV),(moris_index)iG) )
                            {
                                tIsInterfaceFacet = false;
                                break;
                            }
                        }

                        // if this is an interface facet, we need to update the interface data
                        if(tIsInterfaceFacet)
                        {
                            tChildMesh.mark_facet_as_on_interface(iFacet,iG);
                        }
                    }
                }
            }
        }
    }

    // ----------------------------------------------------------------------------------

    bool
    Model::check_for_all_cell_vertices_on_interface()
    {
        bool tPassCheck = true;

        for(moris::uint iCM = 0; iCM < mCutMesh.get_num_child_meshes(); iCM++)
        {
            // active child mesh
            Child_Mesh & tChildMesh = mCutMesh.get_child_mesh(iCM);

            // Element to node
            moris::Matrix< moris::IndexMat > const & tElementToNode = tChildMesh.get_element_to_node();

            // iterate through geometries, check and flag facets that are on the interface
            for(moris::uint iG = 0; iG < mGeometryEngine->get_num_geometries(); iG++)
            {
                // iterate through cells
                for(moris::uint iC = 0 ; iC < tElementToNode.n_rows(); iC++)
                {
                    bool tAllVertsOnInterface = true;

                    for(moris::uint iV = 0; iV < tElementToNode.n_cols(); iV++)
                    {
                        if(! mGeometryEngine->is_interface_vertex(tElementToNode(iC,iV),(moris_index)iG) )
                        {
                            tAllVertsOnInterface = false;
                            break;
                        }
                    }

                    if(tAllVertsOnInterface)
                    {
                        tPassCheck = false;
                        break;
                    }
                }
            }
        }

        return tPassCheck;
    }

    // ----------------------------------------------------------------------------------

    bool
    Model::check_for_degenerated_cells( )
    {
        Cell<moris_index> tDegenerateCells; 

        moris::real tDegenerateTol = MORIS_REAL_MIN;

        // iterate through child meshes
        for(moris::uint iCM = 0; iCM < mCutMesh.get_num_child_meshes(); iCM++)
        {
            // active child mesh
            Child_Mesh & tChildMesh = mCutMesh.get_child_mesh(iCM);

            moris::Matrix< moris::IndexMat > const & tCellInds = tChildMesh.get_element_inds();

            for(moris::uint iC = 0; iC < tCellInds.numel(); iC++)
            {
                moris::mtk::Cell & tCell = mBackgroundMesh.get_mtk_cell(tCellInds(iC));
                moris::real tBulkMeasure = tCell.compute_cell_measure();

                if(tBulkMeasure < tDegenerateTol)
                {
                    tDegenerateCells.push_back(tCellInds(iC));
                }
            }
        }

        if(tDegenerateCells.size() > 0)
        {
            MORIS_LOG_SPEC("Number of Degenerated Cells: " , tDegenerateCells.size());
            return false;
        }

        return true;
    }
    // ----------------------------------------------------------------------------------

    void
    Model::assign_node_requests_identifiers(
            Decomposition_Data & aDecompData,
            moris::moris_index   aMPITag)
    {
        barrier();
        // asserts
        MORIS_ASSERT(aDecompData.tNewNodeId.size() == aDecompData.tNewNodeIndex.size(),
                "Dimension mismatch in assign_node_requests_identifiers");
        MORIS_ASSERT(aDecompData.tNewNodeId.size() == aDecompData.tNewNodeParentRank.size(),
                "Dimension mismatch in assign_node_requests_identifiers");
        MORIS_ASSERT(aDecompData.tNewNodeId.size() == aDecompData.tNewNodeParentIndex.size(),
                "Dimension mismatch in assign_node_requests_identifiers");

        // owned requests and shared requests sorted by owning proc
        Cell<uint> tOwnedRequest;
        Cell<Cell<uint>> tNotOwnedRequests;
        Cell<uint> tProcRanks;
        std::unordered_map<moris_id,moris_id> tProcRankToDataIndex;
        this->sort_new_node_requests_by_owned_and_not_owned(
                aDecompData,
                tOwnedRequest,
                tNotOwnedRequests,
                tProcRanks,
                tProcRankToDataIndex);

        // allocate ids for nodes I own
        moris::moris_id tNodeId  = mBackgroundMesh.allocate_entity_ids(aDecompData.tNewNodeId.size(), EntityRank::NODE);

        // Assign owned request identifiers
        this->assign_owned_request_id(aDecompData, tOwnedRequest, tNodeId);

        // prepare node information request data
        Cell<Matrix<IndexMat>> tOutwardRequests;
        this->setup_outward_requests(aDecompData, tNotOwnedRequests, tProcRanks, tProcRankToDataIndex, tOutwardRequests);

        // send requests to owning processor
        this->send_outward_requests(aMPITag,tProcRanks,tOutwardRequests);

        // hold on to make sure everyone has sent all their information
        barrier();

        // receive the requests
        Cell<Matrix<IndexMat>> tReceivedRequests;
        Cell<uint> tProcsReceivedFrom;
        this->inward_receive_requests(aMPITag, 3, tReceivedRequests, tProcsReceivedFrom);

        // Prepare request answers
        Cell<Matrix<IndexMat>> tRequestAnwers;
        this->prepare_request_answers(aDecompData,tReceivedRequests,tRequestAnwers);

        // send the answers back
        this->return_request_answers(aMPITag+1, tRequestAnwers, tProcsReceivedFrom);

        barrier();

        // receive the answers
        Cell<Matrix<IndexMat>> tReceivedRequestsAnswers;
        this->inward_receive_request_answers(aMPITag+1,1,tProcRanks,tReceivedRequestsAnswers);

        // handle received information
        this->handle_received_request_answers(aDecompData,tOutwardRequests,tReceivedRequestsAnswers,tNodeId);

        MORIS_ERROR(this->verify_successful_node_assignment(aDecompData),
                "Unsuccesssful node assignment detected.");

        barrier();
    }

    // ----------------------------------------------------------------------------------

    void
    Model::sort_new_node_requests_by_owned_and_not_owned(
            Decomposition_Data                    & tDecompData,
            Cell<uint>                            & aOwnedRequests,
            Cell<Cell<uint>>                      & aNotOwnedRequests,
            Cell<uint>                            & aProcRanks,
            std::unordered_map<moris_id,moris_id> & aProcRankToIndexInData)
    {
        // access the mesh data behind the background mesh
        moris::mtk::Mesh & tMeshData = mBackgroundMesh.get_mesh_data();

        // access the communication
        Matrix<IdMat> tCommTable = mBackgroundMesh.get_communication_table();

        // number of new nodes
        moris::uint tNumNewNodes = tDecompData.tNewNodeParentIndex.size();

        // Par rank
        moris::moris_index tParRank = par_rank();

        // resize proc ranks and setup map to comm table
        aProcRanks.resize(tCommTable.numel());
        for(moris::uint i = 0; i <tCommTable.numel(); i++)
        {
            aProcRankToIndexInData[tCommTable(i)] = i;
            aProcRanks(i) = (tCommTable(i));
            aNotOwnedRequests.push_back(Cell<uint>(0));
        }

        // iterate through each node request and figure out the owner
        for(moris::uint i = 0; i <tNumNewNodes; i++)
        {
            // Parent Rank
            enum EntityRank    tParentRank  = tDecompData.tNewNodeParentRank(i);
            moris::moris_index tParentIndex = tDecompData.tNewNodeParentIndex(i);

            // get the owner processor
            moris::moris_index tOwnerProc = tMeshData.get_entity_owner(tParentIndex,tParentRank);

            // If i own the request keep track of the index
            if(tOwnerProc == tParRank)
            {
                aOwnedRequests.push_back(i);
            }
            else
            {
                moris_index tIndex = aProcRankToIndexInData[tOwnerProc];

                aNotOwnedRequests(tIndex).push_back(i);
            }
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Model::assign_owned_request_id(
            Decomposition_Data & aDecompData,
            Cell<uint> const &   aOwnedRequest,
            moris::moris_id &    aNodeId)
    {
        for(moris::uint i = 0; i < aOwnedRequest.size(); i++)
        {
            moris_index tRequestIndex = aOwnedRequest(i);

            // set the new node id
            aDecompData.tNewNodeId(tRequestIndex) = aNodeId;
            aNodeId++;

            // increment number of new nodes with set ids (for assertion purposes)
            aDecompData.mNumNewNodesWithIds++;
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Model::assign_index(Decomposition_Data & aDecompData)
    {
        moris_index tNodeIndex = mBackgroundMesh.get_first_available_index(EntityRank::NODE);

        for(moris::uint i = 0; i < aDecompData.tNewNodeIndex.size(); i++)
        {
            // set the new node index
            aDecompData.tNewNodeIndex(i) = tNodeIndex;
            tNodeIndex++;
        }

        mBackgroundMesh.update_first_available_index(tNodeIndex, EntityRank::NODE);
    }

    // ----------------------------------------------------------------------------------

    void
    Model::setup_outward_requests(
            Decomposition_Data              const & aDecompData,
            Cell<Cell<uint>>                const & aNotOwnedRequests,
            Cell<uint>                      const & aProcRanks,
            std::unordered_map<moris_id,moris_id> & aProcRankToIndexInData,
            Cell<Matrix<IndexMat>>                & aOutwardRequests)
    {
        // size data
        aOutwardRequests.resize(aProcRanks.size());

        // iterate through the processors we need information from and package the matrix
        for(moris::uint i = 0; i < aProcRanks.size(); i++)
        {
            uint tProcRank = aProcRanks(i);

            MORIS_ASSERT(aProcRankToIndexInData.find(tProcRank) != aProcRankToIndexInData.end(),"Proc rank not in map");
            uint tIndexInData = aProcRankToIndexInData[tProcRank];

            uint tNumRequests = aNotOwnedRequests(tIndexInData).size();

            // size the sending matrix
            // column - request
            //   r0 - parent entity id
            //   r1 - parent entity rank
            //   r2 - Secondary id
            if(tNumRequests > 0)
            {
                aOutwardRequests(i) = moris::Matrix<IndexMat>(3,tNumRequests);
            }

            else
            {
                aOutwardRequests(i) = moris::Matrix<IndexMat>(3,1,MORIS_INDEX_MAX);
            }

            // populate matrix to send;
            for(moris::uint j = 0; j < tNumRequests; j++)
            {
                moris_index     tRequestIndex = aNotOwnedRequests(tIndexInData)(j);
                moris_index     tParentIndex  = aDecompData.tNewNodeParentIndex(tRequestIndex);
                moris_index     tSecondaryId  = aDecompData.tSecondaryIdentifiers(tRequestIndex);
                enum EntityRank tParentRank   = aDecompData.tNewNodeParentRank(tRequestIndex);

                aOutwardRequests(i)(0,j) = mBackgroundMesh.get_glb_entity_id_from_entity_loc_index(tParentIndex,tParentRank);
                aOutwardRequests(i)(1,j) = (moris_index)tParentRank;
                aOutwardRequests(i)(2,j) = tSecondaryId;
            }
        }
    }

    // ----------------------------------------------------------------------------------

    bool
    Model::verify_successful_node_assignment(Decomposition_Data & aDecompData)
    {
        uint tNumUnsuccessful = 0;
        for(moris::uint i = 0; i < aDecompData.tNewNodeId.size(); i++)
        {
            if(aDecompData.tNewNodeId(i) == MORIS_INDEX_MAX)
            {
                tNumUnsuccessful++;
            }
        }

        if(tNumUnsuccessful > 0)
        {
            std::cout<<"There were "<<tNumUnsuccessful<<" bad nodes of "<<aDecompData.tNewNodeId.size()<<" total nodes"<<std::endl;
            return false;
        }

        return true;
    }

    // ----------------------------------------------------------------------------------

    void
    Model::send_outward_requests(
            moris_index            const & aMPITag,
            Cell<uint>             const & aProcRanks,
            Cell<Matrix<IndexMat>> & aOutwardRequests)
    {
        // Cell of requests
        Cell<MPI_Request> tRequests(aProcRanks.size());

        // iterate through owned requests and send
        for(moris::uint i = 0; i < aProcRanks.size(); i++)
        {
            tRequests(i) = nonblocking_send(
                    aOutwardRequests(i),
                    aOutwardRequests(i).n_rows(),
                    aOutwardRequests(i).n_cols(),
                    aProcRanks(i),
                    aMPITag);
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Model::inward_receive_requests(
            moris_index            const & aMPITag,
            moris::uint                    aNumRows,
            Cell<Matrix<IndexMat>> &       aReceivedData,
            Cell<uint>             &       aProcRanksReceivedFrom)
    {
        // ensure the sizes are correct.
        aReceivedData.resize(0);
        aProcRanksReceivedFrom.resize(0);

        // access the communication table
        Matrix<IdMat> tCommTable = mBackgroundMesh.get_communication_table();
        moris::uint tCount = 0;
        for(moris::uint i = 0; i<tCommTable.numel(); i++)
        {
            aReceivedData.push_back(Matrix<IndexMat>(1,1));
            aProcRanksReceivedFrom.push_back(tCommTable(i));
            receive(aReceivedData(tCount),aNumRows, tCommTable(i),aMPITag);
            tCount++;
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Model::inward_receive_request_answers(
            moris_index            const & aMPITag,
            moris::uint            const & aNumRows,
            Cell<uint>             const & aProcRanks,
            Cell<Matrix<IndexMat>> &       aReceivedRequestAnswers)
    {
        MPI_Status tStatus;

        for(moris::uint i = 0; i<aProcRanks.size(); i++)
        {
            bool tFlag = sent_message_exists(aProcRanks(i),aMPITag,tStatus);
            while(tFlag == false)
            {
                tFlag = sent_message_exists(aProcRanks(i),aMPITag,tStatus);
            }

            aReceivedRequestAnswers.push_back(Matrix<IndexMat>(1,1));
            receive(aReceivedRequestAnswers(i),aNumRows, aProcRanks(i),aMPITag);
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Model::handle_received_request_answers(
            Decomposition_Data           & aDecompData,
            Cell<Matrix<IndexMat>> const & aRequests,
            Cell<Matrix<IndexMat>> const & aRequestAnswers,
            moris::moris_id              & aNodeId)
    {
        Cell<moris_index> tUnhandledRequestIndices;

        // iterate through received data
        for(moris::uint i = 0; i < aRequests.size(); i++)
        {
            uint tNumReceivedReqs = aRequests(i).n_cols();

            // avoid the dummy message
            if(aRequests(i)(0,0) != MORIS_INDEX_MAX)
            {
                // iterate through received requests
                for(moris::uint j = 0; j < tNumReceivedReqs; j++)
                {
                    moris_id        tParentId      = aRequests(i)(0,j);
                    enum EntityRank tParentRank    = (enum EntityRank) aRequests(i)(1,j);
                    moris_id        tSecondaryId   = aRequests(i)(2,j);
                    moris_index     tParentInd     = mBackgroundMesh.get_mesh_data().get_loc_entity_ind_from_entity_glb_id(tParentId,tParentRank);
                    bool            tRequestExists = false;
                    moris_index     tRequestIndex  = MORIS_INDEX_MAX;

                    if(aDecompData.mHasSecondaryIdentifier)
                    {
                        tRequestExists = aDecompData.request_exists(tParentInd,tSecondaryId,(EntityRank)tParentRank,tRequestIndex);
                    }
                    else
                    {
                        tRequestExists = aDecompData.request_exists(tParentInd,(EntityRank)tParentRank,tRequestIndex);
                    }

                    if(tRequestExists && aRequestAnswers(i)(j))
                    {
                        moris_id tNodeId =aRequestAnswers(i)(j);

                        // meaning the owning processor expected this and gave an answer
                        if(tNodeId < MORIS_ID_MAX && aDecompData.tNewNodeId(tRequestIndex) == MORIS_INDEX_MAX)
                        {
                            // set the new node id
                            aDecompData.tNewNodeId(tRequestIndex) = tNodeId;

                            aDecompData.mNumNewNodesWithIds++;
                        }
                        // The owner did not expect and did not return an answer
                        else
                        {   
                            // keep track of unhandled
                            tUnhandledRequestIndices.push_back(tRequestIndex);
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
                        MORIS_ASSERT(0,"Request does not exist.");
                    }
                }
            }
        }

        // handle the unhandled requests wiht current proc being the owner
        moris::moris_id tNodeId  = mBackgroundMesh.allocate_entity_ids(tUnhandledRequestIndices.size(), EntityRank::NODE);

        for (moris::uint i = 0; i < tUnhandledRequestIndices.size(); i++)
        {
            moris_index tRequestIndex = tUnhandledRequestIndices(i);
            aDecompData.tNewNodeOwner(tRequestIndex) = par_rank();
            aDecompData.tNewNodeId(tRequestIndex) = tNodeId;
            tNodeId++;
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Model::send_outward_requests_reals(
            moris_index const    & aMPITag,
            Cell<uint>  const    & aProcRanks,
            Cell<Matrix<DDRMat>> & aOutwardRequests)
    {
        // iterate through owned requests and send
        for(moris::uint i = 0; i < aProcRanks.size(); i++)
        {
            nonblocking_send(
                    aOutwardRequests(i),
                    aOutwardRequests(i).n_rows(),
                    aOutwardRequests(i).n_cols(),
                    aProcRanks(i),
                    aMPITag);
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Model::inward_receive_requests_reals(
            moris_index const &    aMPITag,
            moris::uint            aNumRows,
            Cell<Matrix<DDRMat>> & aReceivedData,
            Cell<uint>           & aProcRanksReceivedFrom)
    {
        moris::moris_index tParRank = par_rank();
        moris::uint tCount = 0;
        MPI_Status tStatus;
        for(moris::uint i = 0; i<(moris::uint)par_size(); i++)
        {
            if((moris_index)i != tParRank)
            {
                // if there is a sent message from a processor go receive it
                if(sent_message_exists(i,aMPITag,tStatus))
                {
                    aReceivedData.push_back(Matrix<DDRMat>(1,1));

                    aProcRanksReceivedFrom.push_back(i);

                    receive(aReceivedData(tCount),aNumRows, i,aMPITag);

                    tCount++;
                }
            }
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Model::return_request_answers_reals(
            moris_index          const & aMPITag,
            Cell<Matrix<DDRMat>> const & aRequestAnswers,
            Cell<uint>           const & aProcRanks)
    {
        // iterate through owned requests and send
        for(moris::uint i = 0; i < aProcRanks.size(); i++)
        {
            nonblocking_send(aRequestAnswers(i),aRequestAnswers(i).n_rows(),aRequestAnswers(i).n_cols(),aProcRanks(i),aMPITag);
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Model::inward_receive_request_answers_reals(
            moris_index            const & aMPITag,
            moris::uint            const & aNumRows,
            Cell<uint>             const & aProcRanks,
            Cell<Matrix<DDRMat>>         & aReceivedData)
    {
        for(moris::uint i = 0; i<aProcRanks.size(); i++)
        {
            aReceivedData.push_back(Matrix<DDRMat>(1,1));

            receive(aReceivedData(i),aNumRows, aProcRanks(i),aMPITag);
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Model::prepare_request_answers(
            Decomposition_Data           & aDecompData,
            Cell<Matrix<IndexMat>> const & aReceiveData,
            Cell<Matrix<IndexMat>>       & aRequestAnswers)
    {
        // allocate answer size
        aRequestAnswers.resize(aReceiveData.size());

        // iterate through received data
        for(moris::uint i = 0; i < aReceiveData.size(); i++)
        {
            uint tNumReceivedReqs = aReceiveData(i).n_cols();

            aRequestAnswers(i).resize(1,tNumReceivedReqs);

            aRequestAnswers(i)(0) = MORIS_INDEX_MAX;

            // avoid the dummy message
            if(aReceiveData(i)(0,0) != MORIS_INDEX_MAX)
            {
                // iterate through received requests
                for(moris::uint j = 0; j < tNumReceivedReqs; j++)
                {
                    moris_id        tParentId      = aReceiveData(i)(0,j);
                    enum EntityRank tParentRank    = (enum EntityRank) aReceiveData(i)(1,j);
                    moris_id        tSecondaryId   = aReceiveData(i)(2,j);
                    moris_index     tParentInd     = mBackgroundMesh.get_mesh_data().get_loc_entity_ind_from_entity_glb_id(tParentId,tParentRank);
                    bool            tRequestExists = false;
                    moris_index     tRequestIndex  = MORIS_INDEX_MAX;

                    if(aDecompData.mHasSecondaryIdentifier)
                    {
                        tRequestExists = aDecompData.request_exists(
                                tParentInd,
                                tSecondaryId,
                                (EntityRank)tParentRank,
                                tRequestIndex);
                    }
                    else
                    {
                        tRequestExists = aDecompData.request_exists(
                                tParentInd,
                                (EntityRank)tParentRank,
                                tRequestIndex);
                    }

                    if(tRequestExists)
                    {
                        moris_id tNodeId =aDecompData.tNewNodeId(tRequestIndex);

                        aRequestAnswers(i)(j) = tNodeId;

                        if(tNodeId == MORIS_ID_MAX)
                        {
                            std::cout<<"tParentId = "<<tParentId<<" | Rank "<<(uint)tParentRank<<std::endl;
                            //                    MORIS_ERROR(0,"Max node");
                        }
                    }
                    else
                    {
                        aRequestAnswers(i)(j) = MORIS_ID_MAX;
                    }
                }
            }
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Model::return_request_answers(
            moris_index            const & aMPITag,
            Cell<Matrix<IndexMat>> const & aRequestAnswers,
            Cell<uint>             const & aProcRanks)
    {
        // access the communication table
        Matrix<IdMat> tCommTable = mBackgroundMesh.get_communication_table();

        // iterate through owned requests and send
        for(moris::uint i = 0; i < tCommTable.numel(); i++)
        {
            nonblocking_send(aRequestAnswers(i),aRequestAnswers(i).n_rows(),aRequestAnswers(i).n_cols(),tCommTable(i),aMPITag);
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Model::finalize_decomp()
    {
        // Change XTK model decomposition state flag
        mDecomposed = true;

        // Sort the children meshes into groups
        this->sort_children_meshes_into_groups();

        // assign child element indices ( seperated to facilitate mesh cleanup )
        this->assign_child_element_indices(true);

        // give each child cell its id (parallel consistent) and index (not parallel consistent)
        this->assign_child_element_ids();

        // creates mtk cells for all child elements (parent elements are assumed to have mtk cells in the mtk mesh)
        this->create_child_element_mtk_cells();

        // add vertices to child meshes
        this->add_vertices_to_child_meshes();

        // Compute the child element phase using the geometry engine
        // a case where the phase may not be set is when we only do a
        // non-conformal decomposition
        this->set_element_phases();

        // identify local subphases in child mesh
        this->identify_local_subphase_clusters_in_child_meshes();

        // constructs the subphase double side sets internal to a child mesh 
        // I do this here because it also figures out if the child mesh has inter child mesh interfaces
        // this flag is needed for the cleanup cut mesh call
        // this->construct_internal_double_sides_between_subphases();

        // // cleanup the mesh
        // if(mCleanupMesh)
        // {
        //     std::cout<<"Mesh Cleanup"<<std::endl;
        //     Mesh_Cleanup tMeshCleanup(this,&mParameterList);
        //     tMeshCleanup.perform();
        // }


        // // identify local subphases in child mesh
        // this->identify_local_subphase_clusters_in_child_meshes();

        // // add child element to local to global map
        // this->add_child_elements_to_local_to_global_map();

        // // Associate nodes created during decomposition to their child meshes
        // this->associate_nodes_created_during_decomp_to_child_meshes();

        // // set the glb to loc map for all cells
        // this->setup_cell_glb_to_local_map();

        // // assign subphase ids
        // this->assign_subphase_glob_ids();

        // // setup global to local subphase map
        // this->setup_glob_to_loc_subphase_map();
    }

    // ----------------------------------------------------------------------------------

    void
    Model::finalize_mesh_data()
    {
        //
        mBackgroundMesh.setup_local_to_global_maps();

        //
        this->sort_children_meshes_into_groups();

        // set the element phases
        this->set_element_phases();

        // identify local subphases in child mesh
        this->identify_local_subphase_clusters_in_child_meshes();

        // Associate nodes created during decomposition to their child meshes
        this->associate_nodes_created_during_decomp_to_child_meshes();

        // constructs the subphase double side sets internal to a child mesh 
        // I do this here because it also figures out if the child mesh has inter child mesh interfaces
        // this flag is needed for the cleanup cut mesh call
        this->construct_internal_double_sides_between_subphases();

        // set the glb to loc map for all cells
        this->setup_cell_glb_to_local_map();

        // assign subphase ids
        this->assign_subphase_glob_ids();

        // setup global to local subphase map
        this->setup_glob_to_loc_subphase_map();
        // this catches the missed interfaces due to coincidence
        this->catch_all_unhandled_interfaces(); 


        mMeshDataFinalized = true;
    }

    void
    Model::assign_child_element_indices( bool aUpdateAvailable )
    {
        // Allocate global element ids (these need to be give to the children meshes)
        moris_index tElementIndOffset = mBackgroundMesh.get_first_available_index(EntityRank::ELEMENT);

        // set child elements ids in the children meshes which I own and dont share
        Cell<Child_Mesh*> const & tOwnedChildMeshes = mCutMesh.get_owned_child_meshes();
        for(moris::size_t i = 0; i<tOwnedChildMeshes.size(); i++)
        {
            tOwnedChildMeshes(i)->set_child_element_inds(tElementIndOffset);
        }

        // set the indices of the child elements in not owned but not the ids
        // set child elements ids in the children meshes which I own and dont share
        Cell<Child_Mesh*> const & tNotOwnedChildMeshes     = mCutMesh.get_not_owned_child_meshes();
        for(moris::size_t i = 0; i<tNotOwnedChildMeshes.size(); i++)
        {
            tNotOwnedChildMeshes(i)->set_child_element_inds(tElementIndOffset);
        }

        if(aUpdateAvailable)
        {
            // tell the background mesh about the new first available index
            mBackgroundMesh.update_first_available_index(tElementIndOffset,EntityRank::ELEMENT);
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Model::assign_child_element_ids()
    {
        // Set child element ids and indices
        moris::size_t tNumElementsInCutMesh = mCutMesh.get_num_entities(EntityRank::ELEMENT);

        // Allocate global element ids (these need to be give to the children meshes)
        moris_id    tElementIdOffset = mBackgroundMesh.allocate_entity_ids(tNumElementsInCutMesh, moris::EntityRank::ELEMENT);

        // set child elements ids in the children meshes which I own and dont share
        Cell<Child_Mesh*> const & tOwnedChildMeshes = mCutMesh.get_owned_child_meshes();
        for(moris::size_t i = 0; i<tOwnedChildMeshes.size(); i++)
        {
            tOwnedChildMeshes(i)->set_child_element_ids(tElementIdOffset);
        }

        // prepare outward requests
        Cell<Cell<moris_id>>        tNotOwnedChildMeshesToProcs;
        Cell<moris::Matrix<IdMat>>  tOwnedParentCellId;
        Cell<moris::Matrix<IdMat>>  tNumOwnedCellIdsOffsets;
        Cell<uint>                  tProcRanks;
        std::unordered_map<moris_id,moris_id>  tProcRankToDataIndex;

        this->prepare_child_element_identifier_requests(tNotOwnedChildMeshesToProcs, tOwnedParentCellId, tNumOwnedCellIdsOffsets, tProcRanks,tProcRankToDataIndex);

        // send requests
        moris::uint tMPITag = 141;
        this->send_outward_requests(tMPITag, tProcRanks,tOwnedParentCellId);
        this->send_outward_requests(tMPITag+1, tProcRanks, tNumOwnedCellIdsOffsets);

        barrier();

        // receive requests
        Cell<Matrix<IndexMat>> tReceivedParentCellIds;
        Cell<Matrix<IndexMat>> tReceivedParentCellNumChildren;
        Cell<uint> tProcsReceivedFrom1;
        Cell<uint> tProcsReceivedFrom2;
        this->inward_receive_requests(tMPITag, 1, tReceivedParentCellIds, tProcsReceivedFrom1);
        this->inward_receive_requests(tMPITag+1,1, tReceivedParentCellNumChildren, tProcsReceivedFrom2);

        MORIS_ASSERT(tProcsReceivedFrom1.size() == tProcsReceivedFrom2.size(),
                "Size mismatch between procs received from child cell ids and number of child cells");

        Cell<Matrix<IndexMat>> tChildIdOffsets;
        this->prepare_child_cell_id_answers(tReceivedParentCellIds,tReceivedParentCellNumChildren,tChildIdOffsets);

        // return information
        this->return_request_answers(tMPITag+2, tChildIdOffsets, tProcsReceivedFrom1);

        // receive the information
        barrier();

        // receive the answers
        Cell<Matrix<IndexMat>> tReceivedChildIdOffsets;
        this->inward_receive_request_answers(tMPITag+2,1,tProcRanks,tReceivedChildIdOffsets);

        // add child cell ids to not owned child meshes
        this->handle_received_child_cell_id_request_answers(tNotOwnedChildMeshesToProcs,tReceivedChildIdOffsets);

        barrier();
    }

    // ----------------------------------------------------------------------------------

    void
    Model::prepare_child_element_identifier_requests(
            Cell<Cell<moris_id>>       & aNotOwnedChildMeshesToProcs,
            Cell<moris::Matrix<IdMat>> & aOwnedParentCellId,
            Cell<moris::Matrix<IdMat>> & aNumOwnedCellIdsOffsets,
            Cell<uint>                 & aProcRanks,
            std::unordered_map<moris_id,moris_id> & aProcRankToDataIndex)
    {
        // ask owning processor about child element ids
        Cell<Child_Mesh*> const & tNotOwnedChildMeshes     = mCutMesh.get_not_owned_child_meshes();
        Cell<moris_id>    const & tNotOwnedChildMeshOwners = mCutMesh.get_not_owned_child_owners();
        Cell<moris_id>            tCounts(0);
        moris_index tCurrentIndex = 0;

        // access the communication table
        Matrix<IdMat> tCommTable = mBackgroundMesh.get_communication_table();

        for(moris::size_t i = 0; i < tCommTable.numel(); i++)
        {
            aProcRankToDataIndex[tCommTable(i)] = tCurrentIndex;
            aProcRanks.push_back(tCommTable(i));
            aNotOwnedChildMeshesToProcs.push_back(Cell<moris_id>(0));
            tCounts.push_back(0);
            tCurrentIndex++;
        }

        // sort child meshes by owner
        for(moris::size_t i = 0; i < tNotOwnedChildMeshes.size(); i++)
        {
            moris_index tOwnerProc = tNotOwnedChildMeshOwners(i);
            moris_index tProcDataIndex = aProcRankToDataIndex[tOwnerProc];
            aNotOwnedChildMeshesToProcs(tProcDataIndex).push_back(i);
        }

        aOwnedParentCellId.resize(aNotOwnedChildMeshesToProcs.size());
        aNumOwnedCellIdsOffsets.resize(aNotOwnedChildMeshesToProcs.size());

        // iterate through procs and child meshes shared with that processor
        for(moris::size_t i = 0; i < aNotOwnedChildMeshesToProcs.size(); i++)
        {
            // number of child meshes shared with this processor
            moris::uint tNumCM = aNotOwnedChildMeshesToProcs(i).size();

            // allocate matrix
            aOwnedParentCellId(i).resize(1,tNumCM);
            aNumOwnedCellIdsOffsets(i).resize(1,tNumCM);

            for(moris::uint j = 0; j < tNumCM; j++)
            {
                Child_Mesh* tCM = tNotOwnedChildMeshes(aNotOwnedChildMeshesToProcs(i)(j));

                aOwnedParentCellId(i)(j)      = mBackgroundMesh.get_mesh_data().get_glb_entity_id_from_entity_loc_index(tCM->get_parent_element_index(),EntityRank::ELEMENT);
                aNumOwnedCellIdsOffsets(i)(j) = tCM->get_num_entities(EntityRank::ELEMENT);
            }
        }

        for(moris::size_t i = 0; i < tCommTable.numel(); i++)
        {
            if(aNotOwnedChildMeshesToProcs(i).size() == 0)
            {
                aOwnedParentCellId(i).resize(1,1);
                aNumOwnedCellIdsOffsets(i).resize(1,1);
                aOwnedParentCellId(i)(0,0) = MORIS_INDEX_MAX;
                aNumOwnedCellIdsOffsets(i)(0,0) = MORIS_INDEX_MAX;
            }
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Model::prepare_child_cell_id_answers(
            Cell<Matrix<IndexMat>> & aReceivedParentCellIds,
            Cell<Matrix<IndexMat>> & aReceivedParentCellNumChildren,
            Cell<Matrix<IndexMat>> & aChildCellIdOffset)
    {
        MORIS_ASSERT(aReceivedParentCellIds.size() == aReceivedParentCellNumChildren.size(),
                "Mismatch in received parent cell ids and received parent cell number of children");

        // allocate answer size
        aChildCellIdOffset.resize(aReceivedParentCellIds.size());

        // iterate through received data
        for(moris::uint i = 0; i < aReceivedParentCellIds.size(); i++)
        {
            uint tNumReceivedReqs = aReceivedParentCellIds(i).n_cols();

            aChildCellIdOffset(i).resize(1,tNumReceivedReqs);

            if(aReceivedParentCellIds(i)(0) != MORIS_INDEX_MAX)
            {
                // iterate through received requests
                for(moris::uint j = 0; j < tNumReceivedReqs; j++)
                {
                    // parent cell information
                    moris_id tParentId           = aReceivedParentCellIds(i)(0,j);
                    moris_index tParentCellIndex = mBackgroundMesh.get_mesh_data().get_loc_entity_ind_from_entity_glb_id(tParentId,EntityRank::ELEMENT);

                    if(!mBackgroundMesh.entity_has_children(tParentCellIndex,EntityRank::ELEMENT))
                    {
                        std::cout<<"tParentId = "<<tParentId<<" on "<<par_rank()<<std::endl;
                    }
                    // get child mesh
                    MORIS_ASSERT(mBackgroundMesh.entity_has_children(tParentCellIndex,EntityRank::ELEMENT),
                            "Request is made for child element ids on a parent cell not intersected");

                    moris_index tCMIndex = mBackgroundMesh.child_mesh_index(tParentCellIndex,EntityRank::ELEMENT);

                    Child_Mesh & tCM = mCutMesh.get_child_mesh(tCMIndex);

                    MORIS_ASSERT((uint) par_rank() == mBackgroundMesh.get_mesh_data().get_entity_owner(tParentCellIndex,EntityRank::ELEMENT),
                            "I dont own this entity that had info requestsed.");

                    // place in return data
                    MORIS_ASSERT(tCM.get_num_entities(EntityRank::ELEMENT) == (uint)aReceivedParentCellNumChildren(i)(j),
                            "Number of child cells in child mesh do not match number on other processor");

                    // since hmr ownership is not correct
                    if(tCM.get_element_ids().numel()>0)
                    {
                        aChildCellIdOffset(i)(j) = tCM.get_element_ids()(0);
                    }
                }
            }
            else
            {
                aChildCellIdOffset(i)(0) =   MORIS_INDEX_MAX;
            }
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Model::handle_received_child_cell_id_request_answers(
            Cell<Cell<moris_index>> const & aChildMeshesInInNotOwned,
            Cell<Matrix<IndexMat>>  const & aReceivedChildCellIdOffset)
    {
        Cell<Child_Mesh*> const & tNotOwnedChildMeshes = mCutMesh.get_not_owned_child_meshes();

        // iterate through received data
        for(moris::uint i = 0; i < aChildMeshesInInNotOwned.size(); i++)
        {
            uint tNumReceivedReqs = aChildMeshesInInNotOwned(i).size();

            // iterate through received requests
            for(moris::uint j = 0; j < tNumReceivedReqs; j++)
            {
                moris_id tChildMeshInNotOwned = aChildMeshesInInNotOwned(i)(j);
                Child_Mesh* tCM = tNotOwnedChildMeshes(tChildMeshInNotOwned);
                moris_id tChildCellFirstId = aReceivedChildCellIdOffset(i)(j);

                tCM->set_child_element_ids(tChildCellFirstId);
            }
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Model::add_child_elements_to_local_to_global_map()
    {
        // get all child element ids and indexes
        Matrix<IndexMat> tChildElementInds = mCutMesh.get_all_element_inds();
        Matrix<IndexMat> tChildElementIds  = mCutMesh.get_all_element_ids();

        mBackgroundMesh.add_cells_to_global_to_local_map(tChildElementInds,tChildElementIds);
    }

    // ----------------------------------------------------------------------------------

    void
    Model::sort_children_meshes_into_groups()
    {
        // background mesh data
        moris::mtk::Mesh & tMeshData = mBackgroundMesh.get_mesh_data();

        // my proc rank
        moris_index tProcRank = par_rank();

        // number of children meshes
        uint tNumChildrenMeshes = mCutMesh.get_num_child_meshes();

        // allocate data
        Cell<Child_Mesh*>   tOwnedChildrenMeshes(tNumChildrenMeshes);
        Cell<Child_Mesh*>   tNotOwnedChildrenMeshes(tNumChildrenMeshes);
        Cell<moris_id>      tNotOwnedOwningProc(tNumChildrenMeshes);

        // keep track of the number in each group
        uint tOwnedCount    = 0;
        uint tNotOwnedCount = 0;

        for(moris::size_t i = 0; i<mCutMesh.get_num_child_meshes(); i++)
        {
            Child_Mesh & tChildMesh = mCutMesh.get_child_mesh(i);

            moris_index tParentCellInd = tChildMesh.get_parent_element_index();

            // get owner of parent cell
            moris_index tOwnerProc = tMeshData.get_entity_owner(tParentCellInd,EntityRank::ELEMENT);

            // if this processor does not own the element add it to the not owned shared list
            if(tOwnerProc != tProcRank)
            {
                tNotOwnedChildrenMeshes(tNotOwnedCount) = &tChildMesh;
                tNotOwnedOwningProc(tNotOwnedCount) = tOwnerProc;
                tNotOwnedCount++;
            }
            else
            {
                tOwnedChildrenMeshes(tOwnedCount) = &tChildMesh;
                tOwnedCount++;
            }
        }

        // size out extra space
        tOwnedChildrenMeshes.resize(tOwnedCount);
        tNotOwnedChildrenMeshes.resize(tNotOwnedCount);
        tNotOwnedOwningProc.resize(tNotOwnedCount);

        // add to cut mesh
        mCutMesh.add_child_mesh_groups( tOwnedChildrenMeshes, tNotOwnedChildrenMeshes, tNotOwnedOwningProc);
    }

    // ----------------------------------------------------------------------------------

    void
    Model::associate_nodes_created_during_decomp_to_child_meshes()
    {
        // Initialize the data in the XTK mesh
        mBackgroundMesh.allocate_external_node_to_child_mesh_associations();

        // Number of children meshes
        size_t tNumCM = mCutMesh.get_num_child_meshes();
        for(size_t i = 0 ; i < tNumCM; i++)
        {
            // Get reference to the child mesh
            Child_Mesh const & tChildMesh = mCutMesh.get_child_mesh(i);

            // Get reference to the nods in the child mesh node indices
            moris::Matrix<moris::IndexMat> const & tNodeIndices = tChildMesh.get_node_indices();

            // Associate these node indices with their child mesh index
            mBackgroundMesh.associate_external_nodes_to_child_mesh(i,tNodeIndices);
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Model::set_element_phases()
    {
        moris::size_t tNumElem = mBackgroundMesh.get_mesh_data().get_num_entities(EntityRank::ELEMENT);

        // Set element phase indices
        mBackgroundMesh.initialize_element_phase_indices(tNumElem + mCutMesh.get_num_entities(EntityRank::ELEMENT));

        for(moris::size_t i = 0; i<tNumElem; i++)
        {
            if(mBackgroundMesh.entity_has_children(i,EntityRank::ELEMENT))
            {
                moris::size_t tChildMeshIndex = mBackgroundMesh.child_mesh_index(i,EntityRank::ELEMENT);

                Child_Mesh & tChildMesh = mCutMesh.get_child_mesh(tChildMeshIndex);

                moris::Matrix< moris::IndexMat > tElemToNode = tChildMesh.get_element_to_node();

                moris::Matrix< moris::IndexMat > const & tElemInds  = tChildMesh.get_element_inds();

                tChildMesh.initialize_element_phase_mat();

                moris::size_t tNumElem = tChildMesh.get_num_entities(EntityRank::ELEMENT);

                for( moris::size_t j = 0; j<tNumElem; j++)
                {
                    moris::size_t tElemPhaseIndex = this->determine_element_phase_index(j,tElemToNode);
                    mBackgroundMesh.set_element_phase_index(tElemInds(0,j),tElemPhaseIndex);
                    tChildMesh.set_element_phase_index(j,tElemPhaseIndex);
                }
            }
            else
            {
                moris::Matrix< moris::IndexMat > tElementNodes = mBackgroundMesh.get_mesh_data().get_entity_connected_to_entity_loc_inds(i,moris::EntityRank::ELEMENT, moris::EntityRank::NODE);

                if(iscol(tElementNodes))
                {
                    tElementNodes = moris::trans(tElementNodes);
                }

                moris::size_t tElemPhaseIndex = this->determine_element_phase_index(0,tElementNodes);

                mBackgroundMesh.set_element_phase_index(i,tElemPhaseIndex);
            }
        }
    }

    // ----------------------------------------------------------------------------------

    void Model::set_downward_inheritance()
    {
        moris::size_t tNumChildMesh = mCutMesh.get_num_child_meshes();
        Cell<std::pair<moris::moris_index,moris::moris_index>> tXTKElementToCutMeshPairs(tNumChildMesh);

        for(moris::size_t iMesh = 0; iMesh<tNumChildMesh; iMesh++)
        {
            tXTKElementToCutMeshPairs(iMesh) = std::pair<moris::moris_index,moris::moris_index> (mCutMesh.get_parent_element_index(iMesh),iMesh);
        }

        mBackgroundMesh.register_new_downward_inheritance(tXTKElementToCutMeshPairs);
    }

    // ----------------------------------------------------------------------------------

    void  Model::run_first_cut_routine(
            moris::uint                        aGeomIndex,
            moris::Matrix< moris::IndexMat > & aActiveChildMeshIndices,
            moris::Matrix< moris::IndexMat > & aNewPairBool)
    {
        // Note this method is independent of node ids for this reason Background_Mesh is not given the node Ids during this subdivision
        moris::mtk::Mesh & tXTKMeshData = mBackgroundMesh.get_mesh_data();

        // Ensure size is zero
        aActiveChildMeshIndices.set_size(0, 0);
        aNewPairBool.set_size(0, 0);

        // Number of elements
        moris::size_t tNumElements = mBackgroundMesh.get_num_entities(EntityRank::ELEMENT);
        MORIS_ERROR(tNumElements > 0,"Model::run_first_cut_routine - Empty mesh passed to XTK");

        // get the linear background cell indo
        mtk::Cell_Info_Factory tCellInfoFactory;
        moris::mtk::Cell_Info* tGeometricCellInfo = tCellInfoFactory.create_cell_info(mBackgroundMesh.get_parent_cell_topology());

        // collect a vector of node to element connectivity for the entire mesh
        moris::uint tNumNodesPerElement = tGeometricCellInfo->get_num_verts();

        // New child meshes
        moris::size_t tNumNewChildMeshes = 0;

        // New child mesh data
        moris::moris_index tNewIndex = 0;
        Matrix<IndexMat>   tIntersectedElementIndices(0, 0);
        Matrix<IndexMat>   tElementNodeIndices(tNumElements, tNumNodesPerElement);

        Cell<std::pair<moris::moris_index, moris::moris_index>> tNewChildElementPair(0);

        // set intersection counter
        uint tIsectCounter=0;

        // store inital size of vectors
        uint tIsecElementInitialSize = tIntersectedElementIndices.n_cols();
        uint tActiveChildInitialSize = aActiveChildMeshIndices.n_cols();
        uint tNewPairBoolInitialSize = aNewPairBool.n_cols();

        // check that initial sizes are zero
        MORIS_ERROR( tIsecElementInitialSize == 0 && tActiveChildInitialSize == 0 && tNewPairBoolInitialSize == 0,
                "Model::run_first_cut_routine - incorrect initial array sizes.\n");

        // resize arrays
        tIntersectedElementIndices.resize(1, tNumElements);
        aActiveChildMeshIndices.resize(1,tNumElements);
        aNewPairBool.resize(1,tNumElements);

        // Loop over elements to check for intersections
        for (size_t tParentElementIndex = 0; tParentElementIndex < tNumElements; tParentElementIndex++)
        {
            Matrix<IndexMat> tElementNodeIndicesTemp = tXTKMeshData.get_entity_connected_to_entity_loc_inds(
                    tParentElementIndex,
                    moris::EntityRank::ELEMENT,
                    moris::EntityRank::NODE);

            for (moris::uint j = 0; j < tNumNodesPerElement; j++)
            {
                tElementNodeIndices(tParentElementIndex, j) = tElementNodeIndicesTemp(j);
            }

            // is the cell intersected
            bool tIsIntersected = mGeometryEngine->is_intersected(
                    tElementNodeIndicesTemp,
                    mBackgroundMesh.get_selected_node_coordinates_loc_inds(tElementNodeIndicesTemp));

            // Intersected elements are flagged via the Geometry_Engine
            if(tIsIntersected || mTriangulateAll)
            {
                // Element index
                tIntersectedElementIndices(tIsectCounter) = tParentElementIndex;

                // Determine how many meshes need to be registered
                if (not mBackgroundMesh.entity_has_children(tParentElementIndex, EntityRank::ELEMENT))
                {
                    tNewIndex = tNumNewChildMeshes + mCutMesh.get_num_child_meshes();

                    tNewChildElementPair.push_back(
                            std::pair<moris::moris_index, moris::moris_index>(tParentElementIndex, tNewIndex));

                    aActiveChildMeshIndices(tIsectCounter) = tNewIndex;

                    aNewPairBool(tIsectCounter) = 0;

                    tNumNewChildMeshes++;
                }
                else
                {
                    aActiveChildMeshIndices(tIsectCounter) =
                            mBackgroundMesh.child_mesh_index(tParentElementIndex, EntityRank::ELEMENT);

                    aNewPairBool(tIsectCounter) = 1;
                }

                // increase counter for intersected elements
                tIsectCounter++;
            }
        }

        //shrink-to-fit arrays
        tIntersectedElementIndices.resize(1, tIsectCounter);
        aActiveChildMeshIndices.resize(1,tIsectCounter);
        aNewPairBool.resize(1,tIsectCounter);

        // Add the downward pair to the mesh for all the newly created element pairs
        mBackgroundMesh.register_new_downward_inheritance(tNewChildElementPair);

        // Allocate space for more simple meshes in XTK mesh
        mCutMesh.inititalize_new_child_meshes(tNumNewChildMeshes, mModelDimension);

        // figure out the template to use
        enum CellTopology tParentTopo = mBackgroundMesh.get_parent_cell_topology();
        enum TemplateType tParentTemp = TemplateType::INVALID_TEMPLATE_TYPE;

        if (tParentTopo == CellTopology::HEX8 or tParentTopo == CellTopology::HEX27 or tParentTopo == CellTopology::HEX64 )
        {
            tParentTemp = TemplateType::HEX_8;
        }

        else if (tParentTopo == CellTopology::TET4 or tParentTopo == CellTopology::TET10)
        {
            tParentTemp = TemplateType::TET_4;
        }

        else if (tParentTopo == CellTopology::QUAD4 or tParentTopo == CellTopology::QUAD9 or tParentTopo == CellTopology::QUAD16)
        {
            tParentTemp = TemplateType::QUAD_4;
        }
        else if (tParentTopo == CellTopology::TRI3)
        {
            tParentTemp = TemplateType::TRI_3;
        }
        else
        {
            MORIS_ERROR(false, "Invalid background cell topo.");
        }

        moris::Matrix< moris::IndexMat > tPlaceHolder(1, 1);
        for (moris::size_t j = 0; j < tIntersectedElementIndices.length(); j++)
        {
            if (aNewPairBool(j) == 0)
            {
                // Get information to provide ancestry
                // This could be replaced with a proper topology implementation that knows faces, edges based on parent element nodes
                Matrix< IndexMat > tNodetoElemConnVec = tElementNodeIndices.get_row(tIntersectedElementIndices(j));
                Matrix< IndexMat > tEdgetoElemConnInd = tXTKMeshData.get_entity_connected_to_entity_loc_inds(
                        tIntersectedElementIndices(j),
                        moris::EntityRank::ELEMENT,
                        moris::EntityRank::EDGE);
                Matrix< IndexMat > tFacetoElemConnInd = tXTKMeshData.get_entity_connected_to_entity_loc_inds(
                        tIntersectedElementIndices(j),
                        moris::EntityRank::ELEMENT,
                        moris::EntityRank::FACE);

                Matrix< IndexMat > tElementMat = {{tIntersectedElementIndices(j)}};

                Cell<moris::Matrix< moris::IndexMat >> tAncestorInformation = {tPlaceHolder, tEdgetoElemConnInd, tFacetoElemConnInd, tElementMat};
                mCutMesh.initialize_new_mesh_from_parent_element(aActiveChildMeshIndices(j), tParentTemp, tNodetoElemConnVec, tAncestorInformation);

                // add node ids
                mBackgroundMesh.convert_loc_entity_ind_to_glb_entity_ids(EntityRank::NODE, tNodetoElemConnVec);

                // get child mesh
                moris::moris_index tCMIndex = mBackgroundMesh.child_mesh_index(tIntersectedElementIndices(j), EntityRank::ELEMENT);
                Child_Mesh & tChildMesh = mCutMesh.get_child_mesh(tCMIndex);

                // add geometry interface
                tChildMesh.add_new_geometry_interface(aGeomIndex);

                // add node ids
                tChildMesh.add_node_ids(tNodetoElemConnVec);
            }
            else
            {
                // get child mesh
                moris::moris_index tCMIndex = mBackgroundMesh.child_mesh_index(tIntersectedElementIndices(j), EntityRank::ELEMENT);
                Child_Mesh & tChildMesh = mCutMesh.get_child_mesh(tCMIndex);

                // add geometry interface
                tChildMesh.add_new_geometry_interface(aGeomIndex);
            }
        }

        delete tGeometricCellInfo;
    }

    // ----------------------------------------------------------------------------------
    bool
    Model::all_child_meshes_on_same_level()
    {
        moris::uint tNumChildMeshes = mCutMesh.get_num_child_meshes();

        moris::uint tLevel = 0;

        bool tFlag = true;

        // iterate through
        for(moris::uint i = 0; i< tNumChildMeshes; i++)
        {
            // get the child mesh
            Child_Mesh & tChildMesh = mCutMesh.get_child_mesh(i);

            // get the parent cell index
            moris_index tParentCellIndex = tChildMesh.get_parent_element_index();

            moris::uint tCellLevel = mBackgroundMesh.get_mesh_data().get_mtk_cell(tParentCellIndex).get_level();

            if(i == 0)
            {
                tLevel = tCellLevel;
            }
            else if (tCellLevel != tLevel)
            {
                tFlag = false;
            }
        }

        // TODO the current implementation might be a problem for small meshes,
        // when only one single element is intersected on a proc

        return all_land(tFlag);
    }

    void
    Model::create_child_element_mtk_cells()
    {
        moris::uint tNumChildMeshes = mCutMesh.get_num_child_meshes();

        moris::mtk::Mesh const & tMeshData = mBackgroundMesh.get_mesh_data();

        for(moris::uint i=0; i<tNumChildMeshes; i++)
        {
            Child_Mesh & tChildMesh = mCutMesh.get_child_mesh(i);

            moris::moris_index tOwnerProc = tMeshData.get_entity_owner(tChildMesh.get_parent_element_index(),EntityRank::ELEMENT);

            moris::Matrix< moris::IdMat >    const & tElementIds  = tChildMesh.get_element_ids();
            moris::Matrix< moris::IndexMat > const & tElementInds = tChildMesh.get_element_inds();
            // Iterate over elements
            for(moris::uint j = 0; j<tElementIds.numel(); j++)
            {
                mBackgroundMesh.add_child_element_to_mtk_cells(tElementInds(j),tElementIds(j),tOwnerProc,j, &tChildMesh);
            }
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Model::add_vertices_to_child_meshes()
    {
        moris::uint tNumChildMeshes = mCutMesh.get_num_child_meshes();

        for(moris::uint i=0; i<tNumChildMeshes; i++)
        {
            Child_Mesh & tChildMesh = mCutMesh.get_child_mesh(i);
            Cell<moris::mtk::Vertex const *> tVertices = mBackgroundMesh.get_mtk_vertices(tChildMesh.get_node_indices());
            tChildMesh.add_vertices(tVertices);
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Model::setup_cell_glb_to_local_map()
    {
        mCellGlbToLocalMap.clear();
        for(moris::uint i = 0; i < this->get_num_elements_total(); i++)
        {
            moris_id tId = mBackgroundMesh.get_glb_entity_id_from_entity_loc_index((moris_index)i,EntityRank::ELEMENT);

            MORIS_ASSERT(mCellGlbToLocalMap.find(tId) == mCellGlbToLocalMap.end(),"Id already in map");
            mCellGlbToLocalMap[tId] = (moris_index) i;
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Model::identify_local_subphase_clusters_in_child_meshes()
    {
        // get the number of children meshes
        moris::size_t tNumChildMeshes =  mCutMesh.get_num_child_meshes();

        // first subphase index (this is incremented by child meshes on call to set_elemental_subphase)
        // (proc local subphase index, never global)
        moris::moris_index tSubPhaseIndex = mBackgroundMesh.get_mesh_data().get_num_entities(EntityRank::ELEMENT);

        // iterate over children meshes and perform local flood-fill
        for(moris::size_t i = 0; i<tNumChildMeshes; i++)
        {
            // Get child mesh index
            Child_Mesh & tChildMesh = mCutMesh.get_child_mesh(i);

            // Perform local flood-fill on child mesh to identify subphase
            moris::Matrix< moris::IndexMat > tLocalFloodFill = local_child_mesh_flood_fill(tChildMesh);

            // Set the local floodfill data as the elemental subphase values in the child mesh
            // The child mesh then sorts the elements into bins
            tChildMesh.set_elemental_subphase(tSubPhaseIndex,tLocalFloodFill);
        }

        // tell the cut mesh how many subphases there are
        mCutMesh.set_num_subphases(tSubPhaseIndex);


        // tell the cut mesh to setup subphase to child mesh connectivity
        mCutMesh.setup_subphase_to_child_mesh_connectivity();
    }
    // ----------------------------------------------------------------------------------
    void
    Model::construct_internal_double_sides_between_subphases()
    {
        // background mesh
        moris_index tMyProcRank = par_rank();
        moris::mtk::Interpolation_Mesh & tBGMesh = this->get_background_mesh().get_mesh_data();

        uint tNumChildMeshes = this->get_cut_mesh().get_num_child_meshes();

        // iterate through children meshes
        for(moris::uint iCM = 0; iCM < tNumChildMeshes; iCM++)
        {
            // get the child mesh
            Child_Mesh * tChildMesh = & this->get_cut_mesh().get_child_mesh((moris_index)iCM);

            if(tBGMesh.get_entity_owner(tChildMesh->get_parent_element_index(),EntityRank::ELEMENT) == (uint)tMyProcRank)
            {
                // tell the child mesh to construct its double side sets between subphases
                tChildMesh->construct_internal_double_sides_between_subphases();
            }
        }
    }
    // ----------------------------------------------------------------------------------

    void
    Model::assign_subphase_glob_ids()
    {
        // Get the number of subphases
        moris_id tNumSubphases = (moris_id)mCutMesh.get_num_subphases();

        // Allocate global element ids starting at the maximum id in the background mesh (these need to be give to the children meshes)
        moris::moris_id tSubphaseIdOffset = mBackgroundMesh.allocate_entity_ids(tNumSubphases, EntityRank::ELEMENT);

        // set subphase ids in the children meshes which I own
        Cell<Child_Mesh*> const & tOwnedChildMeshes = mCutMesh.get_owned_child_meshes();
        for(moris::size_t i = 0; i<tOwnedChildMeshes.size(); i++)
        {
            moris_id tCellId = mBackgroundMesh.get_mesh_data().get_glb_entity_id_from_entity_loc_index(tOwnedChildMeshes(i)->get_parent_element_index(),EntityRank::ELEMENT);
            // iterate through subphase ids
            tOwnedChildMeshes(i)->set_subphase_id(0,tCellId);
            for(moris::uint j = 1; j <tOwnedChildMeshes(i)->get_num_subphase_bins(); j++)
            {
                tOwnedChildMeshes(i)->set_subphase_id(j,tSubphaseIdOffset);
                tSubphaseIdOffset++;
            }
        }

        // prepare outward requests
        Cell<Cell<moris_id>>        tNotOwnedSubphasesToProcs;
        Cell<Cell<moris_id>>        tCMSubphaseIndices;
        Cell<moris::Matrix<IdMat>>  tParentCellIds;
        Cell<moris::Matrix<IdMat>>  tChildCellIds;
        Cell<moris::Matrix<IdMat>>  tNumChildCellsInSubphase;
        Cell<uint>                  tProcRanks;
        std::unordered_map<moris_id,moris_id>  tProcRankToDataIndex;
        this->prepare_subphase_identifier_requests(tNotOwnedSubphasesToProcs, tCMSubphaseIndices, tParentCellIds, tChildCellIds, tNumChildCellsInSubphase,tProcRanks,tProcRankToDataIndex);

        // send requests
        moris::uint tMPITag = 221;
        this->send_outward_requests(tMPITag, tProcRanks,tParentCellIds);
        this->send_outward_requests(tMPITag+1, tProcRanks, tChildCellIds);
        this->send_outward_requests(tMPITag+2, tProcRanks, tNumChildCellsInSubphase);

        barrier();

        // receive requests
        Cell<Matrix<IndexMat>> tReceivedParentCellIds;
        Cell<Matrix<IndexMat>> tFirstChildCellIds;
        Cell<Matrix<IndexMat>> tReceivedNumChildCellsInSubphase;
        Cell<uint> tProcsReceivedFrom1;
        Cell<uint> tProcsReceivedFrom2;
        this->inward_receive_requests(tMPITag, 1, tReceivedParentCellIds, tProcsReceivedFrom1);
        this->inward_receive_requests(tMPITag+1,1, tFirstChildCellIds, tProcsReceivedFrom2);
        this->inward_receive_requests(tMPITag+2,1, tReceivedNumChildCellsInSubphase, tProcsReceivedFrom2);
        MORIS_ASSERT(tProcsReceivedFrom1.size() == tProcsReceivedFrom2.size(),"Size mismatch between procs received from child cell ids and number of child cells");

        // prepare answers
        Cell<Matrix<IndexMat>> tSubphaseIds;
        this->prepare_subphase_id_answers(tReceivedParentCellIds,tFirstChildCellIds,tReceivedNumChildCellsInSubphase,tSubphaseIds);

        // return information
        this->return_request_answers(tMPITag+2, tSubphaseIds, tProcsReceivedFrom1);

        barrier();

        // receive the answers
        Cell<Matrix<IndexMat>> tReceivedSubphaseIds;
        this->inward_receive_request_answers(tMPITag+2,1,tProcRanks,tReceivedSubphaseIds);

        // add child cell ids to not owned child meshes
        this->handle_received_subphase_id_request_answers(tNotOwnedSubphasesToProcs,tCMSubphaseIndices,tReceivedSubphaseIds);

        barrier();
    }

    // ----------------------------------------------------------------------------------

    void
    Model::prepare_subphase_identifier_requests(
            Cell<Cell<moris_id>>       & aNotOwnedSubphasesToProcs,
            Cell<Cell<moris_id>>       & aSubphaseCMIndices,
            Cell<moris::Matrix<IdMat>> & aParentCellIds,
            Cell<moris::Matrix<IdMat>> & aChildCellIds,
            Cell<moris::Matrix<IdMat>> & aNumChildCellsInSubphase,
            Cell<uint>                 & aProcRanks,
            std::unordered_map<moris_id,moris_id> & aProcRankToDataIndex)
    {
        // access the communication table
        Matrix<IdMat> tCommTable = mBackgroundMesh.get_communication_table();

        // resize proc ranks and setup map to comm table
        aProcRanks.resize(tCommTable.numel());
        for(moris::uint i = 0; i <tCommTable.numel(); i++)
        {
            aProcRankToDataIndex[tCommTable(i)] = i;
            aProcRanks(i) = (tCommTable(i));
            aNotOwnedSubphasesToProcs.push_back(Cell<moris_id>(0));
            aSubphaseCMIndices.push_back(Cell<moris_id>(0));
        }

        // ask owning processor about child element ids
        Cell<Child_Mesh*> const & tNotOwnedChildMeshes     = mCutMesh.get_not_owned_child_meshes();
        Cell<moris_id>    const & tNotOwnedChildMeshOwners = mCutMesh.get_not_owned_child_owners();

        // sort child meshes by owner
        for(moris::size_t i = 0; i < tNotOwnedChildMeshes.size(); i++)
        {
            moris_index tOwnerProc = tNotOwnedChildMeshOwners(i);

            moris_index tProcDataIndex = aProcRankToDataIndex[tOwnerProc];

            // iterate through subphases
            Child_Mesh* tCM = tNotOwnedChildMeshes(i);

            for(moris::uint iS = 0; iS < tCM->get_num_subphase_bins();iS++)
            {
                aNotOwnedSubphasesToProcs(tProcDataIndex).push_back(i);
                aSubphaseCMIndices(tProcDataIndex).push_back(iS);
            }
        }

        aParentCellIds.resize(aNotOwnedSubphasesToProcs.size());
        aChildCellIds.resize(aNotOwnedSubphasesToProcs.size());
        aNumChildCellsInSubphase.resize(aNotOwnedSubphasesToProcs.size());

        // iterate through procs and child meshes shared with that processor
        for(moris::size_t i = 0; i < aNotOwnedSubphasesToProcs.size(); i++)
        {
            // number of child meshes shared with this processor
            moris::uint tNumCM = aNotOwnedSubphasesToProcs(i).size();

            // allocate matrix
            aParentCellIds(i).resize(1,tNumCM);
            aChildCellIds(i).resize(1,tNumCM);
            aNumChildCellsInSubphase(i).resize(1,tNumCM);

            if(tNumCM == 0)
            {
                aParentCellIds(i).resize(1,1);
                aChildCellIds(i).resize(1,1);
                aNumChildCellsInSubphase(i).resize(1,1);
                aParentCellIds(i)(0) = MORIS_INDEX_MAX;
                aChildCellIds(i)(0) = MORIS_INDEX_MAX;
                aNumChildCellsInSubphase(i)(0) = MORIS_INDEX_MAX;
            }

            for(moris::uint j = 0; j < tNumCM; j++)
            {
                Child_Mesh* tCM = tNotOwnedChildMeshes(aNotOwnedSubphasesToProcs(i)(j));

                // subphase groups
                Cell<moris::Matrix< moris::IndexMat >> const & tSubphaseClusters = tCM->get_subphase_groups();

                // element ids
                Matrix<IdMat> const & tCellIds  = tCM->get_element_ids();

                // subphase index
                moris_index tSPIndex = aSubphaseCMIndices(i)(j);

                // cell index
                moris_index tCMCellInd = tSubphaseClusters(tSPIndex)(0);

                aParentCellIds(i)(j)           = mBackgroundMesh.get_glb_entity_id_from_entity_loc_index(tCM->get_parent_element_index(),EntityRank::ELEMENT);
                aChildCellIds(i)(j)            = tCellIds(tCMCellInd);
                aNumChildCellsInSubphase(i)(j) = tSubphaseClusters(tSPIndex).numel();
            }
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Model::prepare_subphase_id_answers(
            Cell<Matrix<IndexMat>> & aReceivedParentCellIds,
            Cell<Matrix<IndexMat>> & aFirstChildCellIds,
            Cell<Matrix<IndexMat>> & aReceivedNumChildCellsInSubphase,
            Cell<Matrix<IndexMat>> & aSubphaseIds)
    {
        MORIS_ASSERT(aReceivedParentCellIds.size() == aFirstChildCellIds.size(),
                "Mismatch in received parent cell ids and received parent cell number of children");

        // allocate answer size
        aSubphaseIds.resize(aReceivedParentCellIds.size());

        // iterate through received data
        for(moris::uint i = 0; i < aReceivedParentCellIds.size(); i++)
        {
            uint tNumReceivedReqs = aReceivedParentCellIds(i).n_cols();

            aSubphaseIds(i).resize(1,tNumReceivedReqs);

            if(aReceivedParentCellIds(i)(0,0) != MORIS_INDEX_MAX)
            {
                // iterate through received requests
                for(moris::uint j = 0; j < tNumReceivedReqs; j++)
                {
                    // parent cell information
                    moris_id tParentId           = aReceivedParentCellIds(i)(0,j);
                    moris_index tParentCellIndex = mBackgroundMesh.get_mesh_data().get_loc_entity_ind_from_entity_glb_id(tParentId,EntityRank::ELEMENT);

                    // Child cell information
                    moris_id tChildCellId = aFirstChildCellIds(i)(j);

                    // number of child cells in subphase
                    moris_id tNumChildCells = aReceivedNumChildCellsInSubphase(i)(j);

                    // get child mesh
                    MORIS_ASSERT(mBackgroundMesh.entity_has_children(tParentCellIndex,EntityRank::ELEMENT),
                            "Request is made for child element ids on a parent cell not intersected");

                    moris_index tCMIndex = mBackgroundMesh.child_mesh_index(tParentCellIndex,EntityRank::ELEMENT);
                    Child_Mesh & tCM     = mCutMesh.get_child_mesh(tCMIndex);

                    // figure out which subphase this child cell belongs to in the given child mesh
                    Matrix<IdMat> const & tChildMeshCellIds = tCM.get_element_ids();

                    moris_id tSubphaseId         = MORIS_ID_MAX;
                    moris_index tCMSubphaseIndex = MORIS_ID_MAX;

                    for(moris::uint iCM = 0; iCM < tChildMeshCellIds.numel(); iCM++)
                    {
                        if(tChildMeshCellIds(iCM) == tChildCellId)
                        {
                            tSubphaseId      = tCM.get_element_subphase_id(iCM);
                            tCMSubphaseIndex = tCM.get_element_subphase_index(iCM);

                            Cell<moris::Matrix< moris::IndexMat >> const & tSubphaseClusters = tCM.get_subphase_groups();

                            MORIS_ERROR( (moris_id) tSubphaseClusters(tCMSubphaseIndex).numel() == tNumChildCells,
                                    "Number of cells in subphase mismatch");
                        }
                    }

                    MORIS_ERROR(tSubphaseId!= MORIS_ID_MAX,"Child cell id not found in child mesh");

                    // place in return data
                    aSubphaseIds(i)(j) = tSubphaseId;
                }
            }
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Model::handle_received_subphase_id_request_answers(
            Cell<Cell<moris_index>>    const & aChildMeshesInNotOwned,
            Cell<Cell<moris_index>>    const & aCMSubphaseIndices,
            Cell<Matrix<IndexMat>>     const & aReceivedSubphaseIds)
    {
        Cell<Child_Mesh*> const & tNotOwnedChildMeshes     = mCutMesh.get_not_owned_child_meshes();

        // iterate through received data
        for(moris::uint i = 0; i < aChildMeshesInNotOwned.size(); i++)
        {
            uint tNumReceivedReqs = aChildMeshesInNotOwned(i).size();

            // iterate through received requests
            for(moris::uint j = 0; j < tNumReceivedReqs; j++)
            {
                moris_id tChildMeshInNotOwned = aChildMeshesInNotOwned(i)(j);
                Child_Mesh* tCM = tNotOwnedChildMeshes(tChildMeshInNotOwned);
                moris_id tSubphaseIndex = aCMSubphaseIndices(i)(j);
                moris_id tSubphaseId = aReceivedSubphaseIds(i)(j);

                tCM->set_subphase_id(tSubphaseIndex,tSubphaseId);
            }
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Model::setup_glob_to_loc_subphase_map()
    {

        mGlobalToLocalSubphaseMap.clear();
        for(moris::uint i = 0; i < mCutMesh.get_num_subphases(); i++)
        {
            moris_id tSubphaseId = this->get_subphase_id((moris_id)i);
            MORIS_ASSERT(mGlobalToLocalSubphaseMap.find(tSubphaseId) == mGlobalToLocalSubphaseMap.end(),"Subphase id already in map");
            mGlobalToLocalSubphaseMap[tSubphaseId] = i;
        }
    }

    // ----------------------------------------------------------------------------------
    // Unzipping Child Mesh Source code
    // ----------------------------------------------------------------------------------

    void
    Model::unzip_child_mesh()
    {
        // start the clock
        std::clock_t start = std::clock();

        MORIS_ERROR(mDecomposed,"Prior to unzip_child_mesh, the decomposition process must be called");

        // unzip the interface
        unzip_child_mesh_internal();

        mUnzipped = true;
        if(moris::par_rank() == 0  && mVerbose)
        {
            std::cout<<"XTK: Child mesh unzipping completed in "<< (std::clock() - start) / (double)(CLOCKS_PER_SEC)<<" s."<<std::endl;
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Model::unzip_child_mesh_internal()
    {
        // get the number of children meshes
        moris::size_t tNumChildMeshes =  mCutMesh.get_num_child_meshes();

        for(moris::size_t i = 0; i<tNumChildMeshes; i++)
        {
            // Get child mesh index
            //        Child_Mesh & tChildMesh = mCutMesh.get_child_mesh(i);
        }
    }

    // ----------------------------------------------------------------------------------
    // Unzipping Interface Source code
    // ----------------------------------------------------------------------------------
    void
    Model::unzip_interface()
    {
        // start the clock
        std::clock_t start = std::clock();

        // unzip the interface
        unzip_interface_internal();

        mUnzipped = true;
        if(moris::par_rank() == 0  && mVerbose)
        {
            std::cout<<"XTK: Interface unzipping completed in "<< (std::clock() - start) / (double)(CLOCKS_PER_SEC)<<" s."<<std::endl;
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Model::unzip_interface_internal()
    {
        // Get the number of geometries (we need to unzip each interface)
        uint tNumGeoms = mGeometryEngine->get_num_geometries();

        // Get the interface nodes (wrt all geometries)
        Cell<moris::Matrix<moris::IndexMat>> tAllInterfaceNodeIds = mBackgroundMesh.get_interface_nodes_glob_ids();
        Cell<moris::Matrix<moris::IndexMat>> tAllInterfaceNodeInds = mBackgroundMesh.get_interface_nodes_loc_inds();

        // Keep count of which interface node index in matrices were in
        for(uint iG = 0; iG<tNumGeoms; iG++)
        {
            // Interface node wrt geometry iG
            moris::Matrix<moris::IndexMat> const & tInterfaceNodeIds = tAllInterfaceNodeIds(iG);
            moris::Matrix<moris::IndexMat> const & tInterfaceNodeInds = tAllInterfaceNodeInds(iG);

            // Number of interface nodes wrt geometry iG (and check the sizes of the interface information)
            uint tNumInterfaceNodes = tInterfaceNodeIds.numel();
            MORIS_ASSERT(tInterfaceNodeIds.numel() == tInterfaceNodeInds.numel(), "Interface Ids and Indices dimension mismatch");

            // Assign node ids and indices ( row - 0 Node ids, row 1 - New node ids)
            moris::Matrix<moris::IdMat> tNewUnzippedNodeIds((size_t)tNumInterfaceNodes);
            moris::Matrix<moris::IndexMat> tNewUnzippedNodeInds((size_t)tNumInterfaceNodes);
            this->unzip_interface_internal_assign_node_identifiers(tNumInterfaceNodes,tNewUnzippedNodeInds,tNewUnzippedNodeIds);

            // Add new nodes to the mesh (as a copy of the existing node)
            mBackgroundMesh.batch_create_new_nodes_as_copy_of_other_nodes(tInterfaceNodeInds,tNewUnzippedNodeIds,tNewUnzippedNodeInds);

            // Allocate space in background mesh interface node flags
            mBackgroundMesh.allocate_space_in_interface_node_flags(tNumInterfaceNodes, mGeometryEngine->get_num_geometries());

            // Mark the newly created nodes as interface nodes
            mBackgroundMesh.mark_nodes_as_interface_node_loc_inds(tNewUnzippedNodeInds,iG);

            // unzip_child_mesh_index
            this->unzip_interface_internal_modify_child_mesh(iG,tInterfaceNodeInds,tNewUnzippedNodeInds,tNewUnzippedNodeIds);
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Model::unzip_interface_internal_assign_node_identifiers(
            moris::uint                   aNumNodes,
            moris::Matrix<moris::IdMat> & aUnzippedNodeIndices,
            moris::Matrix<moris::IdMat> & aUnzippedNodeIds)
    {
        // Verify sizes
        MORIS_ASSERT(aUnzippedNodeIndices.numel() == aNumNodes,
                "Size mismatch between aNumNodes and aUnzippedNodeIndices. Please pre-allocate these matrices ");

        MORIS_ASSERT(aUnzippedNodeIds.numel() == aNumNodes,
                "Size mismatch between aNumNodes and aUnzippedNodeIds.  Please pre-allocate these matrices ");

        // Ask the mesh for new node ids
        moris::moris_index tNodeIndexOffset = mBackgroundMesh.get_first_available_index(EntityRank::NODE);
        moris::moris_id    tNodeIdOffset    = mBackgroundMesh.allocate_entity_ids(aNumNodes,EntityRank::NODE);

        // Iterate new nodes and assign new node ids
        for( uint iN = 0; iN<aNumNodes; iN++)
        {

            // TODO: ADD PARALLEL OWNERSHIP STUFF HERE TO ASSIGN CONSISTENT NODE IDS
            // Give node global ids
            aUnzippedNodeIds(iN) = tNodeIdOffset;

            // increase the id offset
            tNodeIdOffset++;

            // Give nodes processor indices
            aUnzippedNodeIndices(iN) = tNodeIndexOffset;

            // increase the node index offset
            tNodeIndexOffset++;
        }

        // update the first available node index in our background mesh
        mBackgroundMesh.update_first_available_index(tNodeIndexOffset,EntityRank::NODE);

    }

    // ----------------------------------------------------------------------------------

    void
    Model::unzip_interface_internal_modify_child_mesh(
            moris::uint                         aGeometryIndex,
            moris::Matrix<moris::IdMat> const & aInterfaceNodeIndices,
            moris::Matrix<moris::IdMat> const & aUnzippedNodeIndices,
            moris::Matrix<moris::IdMat> const & aUnzippedNodeIds)
    {

        // from interface node indices, figure out which interface nodes live in which interface
        moris::Cell<moris::Cell< moris::moris_index >> tChildMeshInterfaceNodes =
                unzip_interface_internal_collect_child_mesh_to_interface_node(
                        aInterfaceNodeIndices,
                        aUnzippedNodeIndices,
                        aUnzippedNodeIds);

        // Flag indicating there is an interface without an element pair
        bool tNoPairFlag = false;

        // Iterate through child meshes and add new unzipped nodes
        uint tCMIndex = 0;
        for(auto iCM = tChildMeshInterfaceNodes.begin(); iCM != tChildMeshInterfaceNodes.end(); ++iCM)
        {
            // number of interface nodes in this child mesh
            uint tNumCMInterfaceNodes = iCM->size();

            // Allocate matrices of interface node indices
            moris::Matrix< moris::IndexMat > tCMInterfaceNodeIndices(1,tNumCMInterfaceNodes);
            moris::Matrix< moris::IndexMat > tCMUnzippedInterfaceNodeIndices(1,tNumCMInterfaceNodes);
            moris::Matrix< moris::IdMat >    tCMUnzippedInterfaceNodeIds(1,tNumCMInterfaceNodes);

            // Collect information on the interface nodes on this child mesh
            for(moris::uint iN = 0; iN<tNumCMInterfaceNodes; iN++)
            {
                // node index local to the numbering scheme in interface nodes
                moris::moris_index tInterfaceLocInd = (*iCM)(iN);

                tCMInterfaceNodeIndices(iN)         = aInterfaceNodeIndices(tInterfaceLocInd);
                tCMUnzippedInterfaceNodeIndices(iN) = aUnzippedNodeIndices(tInterfaceLocInd);
                tCMUnzippedInterfaceNodeIds(iN)     = aUnzippedNodeIds(tInterfaceLocInd);
            }

            // Tell the child mesh to unzip it's interface
            Child_Mesh & tChildMesh = mCutMesh.get_child_mesh(tCMIndex);

            // initialize the unzipping (which basically allocated node to element connectivity
            // in the child mesh because it is not typically needed)
            tChildMesh.initialize_unzipping();

            // Ask the child mesh to construct interface element pairs
            moris::Matrix<moris::IndexMat> tInterfaceElementPairsCMIndex;
            moris::Matrix<moris::IndexMat> tInterfaceSideOrdinals;

            tChildMesh.unzip_child_mesh_interface_get_interface_element_pairs(aGeometryIndex,tNoPairFlag,tInterfaceElementPairsCMIndex,tInterfaceSideOrdinals);

            // Convert the pairs to processor local indices because we need to be able to access the element phase index
            moris::Matrix<moris::IndexMat> tInterfaceElementPairs = tChildMesh.convert_to_proc_local_elem_inds(tInterfaceElementPairsCMIndex);

            // TODO: Add method to resolve cross child mesh element pairs for when the interface coincides with a parent face
            // NOTE: By using the sign of the geometry value, it really shouldnt take a whole lot of work to accomodated
            MORIS_ERROR(!tNoPairFlag,
                    " in unzip_interface_internal_modify_child_mesh, interface detected on a child mesh boundary. Currently, no method is implemented to resolve this");

            // Take the child mesh pairs and determine who gets which id
            // This output is either a 0 or 1, meaning the first or second element of the pair gets the unzipped nodes
            moris::Matrix< moris::IndexMat > tElementWhichKeepsOriginalNodes =
                    this->unzip_interface_internal_assign_which_element_uses_unzipped_nodes(aGeometryIndex,tInterfaceElementPairs);

            // Get the elements on the boundary
            tChildMesh.unzip_child_mesh_interface(
                    aGeometryIndex,
                    tInterfaceElementPairsCMIndex,
                    tElementWhichKeepsOriginalNodes,
                    tCMInterfaceNodeIndices,
                    tCMUnzippedInterfaceNodeIndices,
                    tCMUnzippedInterfaceNodeIds);

            // Construct interface elements
            unzip_interface_construct_interface_elements(aGeometryIndex,
                    tInterfaceElementPairs,
                    tInterfaceSideOrdinals);

            tChildMesh.finalize_unzipping();

            tCMIndex++;
        }

        unzip_interface_assign_element_identifiers();
    }

    // ----------------------------------------------------------------------------------

    moris::Matrix< moris::IndexMat >
    Model::unzip_interface_internal_assign_which_element_uses_unzipped_nodes(
            moris::moris_index                       aGeometryIndex,
            moris::Matrix< moris::IndexMat > const & aInterfaceElementPairs )
    {
        // specify which geometry sign gets to keep the nodes
        moris::moris_index tValWhichUsesUnzipped = 0;

        // The rule used here is whichever phase gets a 1 from the phase table with respect to the current geometry index
        // gets to keep the original node. The other element changes it's nodes to the unzipped indices.
        moris::Matrix< moris::IndexMat > tElementWhichKeepsUsesUnzippedNodes(aInterfaceElementPairs.n_cols());

        // number of pairs
        moris::uint tNumPairs = aInterfaceElementPairs.n_cols();

        // allocate
        moris::moris_index tElement0           = MORIS_INDEX_MAX;
        moris::moris_index tElement1           = MORIS_INDEX_MAX;

        // iterate through pairs
        for(moris::uint iP = 0; iP<tNumPairs; iP++)
        {
            // set this back to moris index max so we dont run into issues if there is actually only one element in the pair
            moris::moris_index tElement0GeomSign   = MORIS_INDEX_MAX; // sign of phase of element wrt to aGeometryIndex (0 for negative, 1 for positive)
            moris::moris_index tElement1GeomSign   = MORIS_INDEX_MAX;// sign of phase of element wrt to aGeometryIndex (0 for negative, 1 for positive)

            // Element indices
            tElement0 = aInterfaceElementPairs(0,iP);
            tElement1 = aInterfaceElementPairs(1,iP);

            // Figure out which element in the pair gets to keep the original
            if(tElement0 != MORIS_INDEX_MAX)
            {
                tElement0GeomSign = 0;

                if(tElement0GeomSign == tValWhichUsesUnzipped)
                {
                    tElementWhichKeepsUsesUnzippedNodes(iP) = 0;
                }
            }

            if(tElement1 != MORIS_INDEX_MAX)
            {
                tElement1GeomSign = 0;
                if(tElement1GeomSign == tValWhichUsesUnzipped)
                {
                    tElementWhichKeepsUsesUnzippedNodes(iP) = 1;
                }
            }

            // Make sure both don't end up with the same sign
            MORIS_ASSERT(tElement1GeomSign != tElement0GeomSign,"Both elements in an interface pair returned the same phase sign");

            MORIS_ASSERT(tElement0GeomSign != MORIS_INDEX_MAX,"tElement0GeomSign no pair");
            MORIS_ASSERT(tElement1GeomSign != MORIS_INDEX_MAX,"tElement1GeomSign no pair");

        }

        return tElementWhichKeepsUsesUnzippedNodes;
    }

    // ----------------------------------------------------------------------------------

    moris::Cell<moris::Cell< moris::moris_index >>
    Model::unzip_interface_internal_collect_child_mesh_to_interface_node(
            moris::Matrix<moris::IdMat> const & aInterfaceNodeIndices,
            moris::Matrix<moris::IdMat> const & aUnzippedNodeIndices,
            moris::Matrix<moris::IdMat> const & aUnzippedNodeIds)
            {
        // Allocate cell to keep track of the node indices of each child mesh
        moris::uint tNumChildMeshes = mCutMesh.get_num_child_meshes();

        moris::Cell<moris::Cell< moris::moris_index >> tChildMeshInterfaceNodes(tNumChildMeshes);

        // Iterate through interface node indices
        for(moris::uint iN = 0; iN < aInterfaceNodeIndices.numel(); iN++)
        {
            // Get the child meshes which have the interface node
            moris::Matrix< moris::IndexMat > tNodeChildMeshIndices = mBackgroundMesh.get_node_child_mesh_assocation(aInterfaceNodeIndices(iN));

            // Iterate through child meshes and mark interface node indices in this child mesh
            for(moris::uint iCM = 0; iCM<tNodeChildMeshIndices.numel(); iCM++)
            {
                moris::moris_index tCMIndex = tNodeChildMeshIndices(iCM);
                tChildMeshInterfaceNodes(tCMIndex).push_back(iN);
            }
        }

        return tChildMeshInterfaceNodes;
            }

    // ----------------------------------------------------------------------------------

    void
    Model::unzip_interface_construct_interface_elements(moris::uint aGeometryIndex,
            moris::Matrix< moris::IndexMat > const & aElementPairs,
            moris::Matrix< moris::IndexMat > const & aSideOrdinalPairs)
    {
        moris::uint tNumPairs = aElementPairs.n_cols();

        moris::Cell<const moris::mtk::Cell*> tPairCells(2);

        for(moris::uint i = 0; i <tNumPairs; i++)
        {
            tPairCells(0) = mBackgroundMesh.get_child_element_mtk_cell_ptr(aElementPairs(0,i));
            tPairCells(1) = mBackgroundMesh.get_child_element_mtk_cell_ptr(aElementPairs(1,i));

            // give ownership information to left element
            moris_index tOwnerOfElem0 = tPairCells(0)->get_owner();

            // construct an interface element
            Interface_Element tInterfaceElement;
            tInterfaceElement.set_element_pair_and_side_ordinal(tPairCells,aSideOrdinalPairs.get_column(i));
            tInterfaceElement.set_element_owner(tOwnerOfElem0);

            // Add interface element to cut mesh
            mCutMesh.add_interface_element(tInterfaceElement);
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Model::unzip_interface_assign_element_identifiers()
    {
        moris::Cell<Interface_Element> & tInterfaceElements = mCutMesh.get_interface_elements();
        moris::uint tNumInterfaceElements = tInterfaceElements.size();

        // Allocate ids
        moris::moris_index tIdOffset    = mBackgroundMesh.allocate_entity_ids(tNumInterfaceElements,EntityRank::ELEMENT);
        moris::moris_index tIndexOffset = mBackgroundMesh.get_first_available_index(EntityRank::ELEMENT);

        for(moris::uint i = 0; i <tNumInterfaceElements; i++)
        {
            tInterfaceElements(i).set_element_identifiers(tIndexOffset,tIdOffset);
            tIndexOffset++;
            tIdOffset++;
        }

        mBackgroundMesh.update_first_available_index(tIndexOffset,EntityRank::ELEMENT);
    }

    // ----------------------------------------------------------------------------------
    // Enrichment Source code
    // ----------------------------------------------------------------------------------

    void
    Model::perform_basis_enrichment(
            enum EntityRank  const & aBasisRank,
            moris_index      const & aMeshIndex)
    {
        Tracer tTracer( "XTK", "Enrichment", "Enrich" );

        if(!mMeshDataFinalized)
        {
            this->finalize_mesh_data();
        }

        MORIS_ERROR(mDecomposed,"Prior to computing basis enrichment, the decomposition process must be called");

        // allocate some new enriched interpolation and integration meshes
        mEnrichedIntegMesh.resize(aMeshIndex+1,nullptr);
        mEnrichedInterpMesh.resize(aMeshIndex+1,nullptr);

        this->perform_basis_enrichment_internal(aBasisRank,{{aMeshIndex}});

        // Change the enrichment flag
        mEnriched = true;
    }

    // ----------------------------------------------------------------------------------

    void
    Model::perform_basis_enrichment(
            enum EntityRank  const & aBasisRank,
            Matrix<IndexMat> const & aMeshIndex)
    {
        Tracer tTracer( "XTK", "Enrichment");

        if(!mMeshDataFinalized)
        {
            this->finalize_mesh_data();
        }
        MORIS_ERROR(mDecomposed,"Prior to computing basis enrichment, the decomposition process must be called");

        // allocate some new enriched interpolation and integration meshes
        mEnrichedIntegMesh.resize(aMeshIndex.numel()+1,nullptr);
        mEnrichedInterpMesh.resize(aMeshIndex.numel()+1,nullptr);

        this->perform_basis_enrichment_internal(aBasisRank,aMeshIndex);


        // Change the enrichment flag
        mEnriched = true;
    }

    void Model::perform_hanging_node_identification()
    {
        Tracer tTracer( "XTK", "Identify hanging nodes" );

        MORIS_ERROR(mDecomposed,"Mesh needs to be decomposed prior to identifying hanging nodes");

        MORIS_ERROR(mEnriched,"Mesh needs to be enriched prior to identifying hanging nodes");

        // get the interpolation mesh
        moris::mtk::Interpolation_Mesh & tInterpMesh = mBackgroundMesh.get_mesh_data();

        // iterate through child meshes
        for(moris::uint iCM = 0; iCM < mCutMesh.get_num_child_meshes(); iCM++)
        {
            // active child mesh
            Child_Mesh & tChildMesh = mCutMesh.get_child_mesh(iCM);

            // get the neighbors
            Matrix<IndexMat> tElementNeighors = tInterpMesh.
                    get_elements_connected_to_element_and_face_ind_loc_inds(tChildMesh.get_parent_element_index() );

            moris::Cell< moris_index > tTransitionFacets;

            // iterate through neighbor
            for(moris::uint iN = 0; iN < tElementNeighors.n_cols(); iN++)
            {
                moris_index tNeighborCellIndex = tElementNeighors( 0, iN );

                moris_index tSharedFaceIndex = tElementNeighors( 1, iN );

                if( !mBackgroundMesh.entity_has_children( tNeighborCellIndex, EntityRank::ELEMENT )  )
                {
                    tTransitionFacets.push_back( tSharedFaceIndex );
                }
            }

            tChildMesh.identify_hanging_nodes( tTransitionFacets );

        }
    }

    void
    Model::probe_bg_cell(Matrix<IndexMat> const & tBGCellIds)
    {
        Tracer tTracer( "XTK", "BG Cell Probe" );

        for(moris::uint i = 0; i < tBGCellIds.numel(); i++)
        {
            Tracer tTracer("XTK", "BG Cell Probe", "Cell Id " + std::to_string(tBGCellIds(i)));

            moris_index tIndex = mBackgroundMesh.get_mesh_data().get_loc_entity_ind_from_entity_glb_id(tBGCellIds(i), EntityRank::ELEMENT);
            mtk::Cell &tCell = mBackgroundMesh.get_mesh_data().get_mtk_cell(tIndex);
            Matrix<IndexMat> tVertexIds = tCell.get_vertex_ids();
            moris::Cell<mtk::Vertex *> tVertexPtrs = tCell.get_vertex_pointers();
            Matrix<IndexMat> tVertexOwner(1, tVertexPtrs.size());

            MORIS_LOG_SPEC("Cell Id", tBGCellIds(i));
            MORIS_LOG_SPEC("Cell Index", tIndex);
            MORIS_LOG_SPEC("Cell Owner", tCell.get_owner());
            MORIS_LOG_SPEC("Vertex Ids", ios::stringify_log(tVertexIds));
            MORIS_LOG_SPEC("Vertex Owners", ios::stringify_log(tVertexOwner));

            // collect geometric info
            uint tNumGeom = mGeometryEngine->get_num_geometries();
            Cell<Matrix<DDRMat>> tVertexGeomVals(tNumGeom, Matrix<DDRMat>(1, tVertexPtrs.size()));
            for (moris::uint iG = 0; iG < tNumGeom; iG++)
            {
                for (moris::uint iV = 0; iV < tVertexPtrs.size(); iV++)
                {
                    tVertexGeomVals(iG)(iV) = mGeometryEngine->get_field_value(iG, (uint)tVertexPtrs(iV)->get_index(), tVertexPtrs(iV)->get_coords());
                }
                MORIS_LOG_SPEC("Geom Field " + std::to_string(iG) + " Vals", ios::stringify_log(tVertexGeomVals(iG)));
            }
        }

    }

    // ----------------------------------------------------------------------------------

    Enrichment const &
    Model::get_basis_enrichment()
    {
        MORIS_ASSERT(mEnriched,
                "Cannot get basis enrichment from an XTK model which has not called perform_basis_enrichment ");

        return *mEnrichment;
    }

    // ----------------------------------------------------------------------------------

    Enriched_Interpolation_Mesh &
    Model::get_enriched_interp_mesh(moris::moris_index aIndex)
    {
        MORIS_ASSERT(mEnriched,
                "Cannot get enriched interpolation mesh from an XTK model which has not called perform_basis_enrichment ");

        return *(mEnrichedInterpMesh(aIndex));
    }

    // ----------------------------------------------------------------------------------

    Enriched_Integration_Mesh &
    Model::get_enriched_integ_mesh(moris::moris_index aIndex)
    {
        MORIS_ASSERT(mEnriched,
                "Cannot get enriched integration mesh from an XTK model which has not called perform_basis_enrichment ");

        return *(mEnrichedIntegMesh(aIndex));
    }

    // ----------------------------------------------------------------------------------

    void
    Model::perform_basis_enrichment_internal(
            enum EntityRank  const & aBasisRank,
            Matrix<IndexMat> const & aMeshIndex)
    {

        // initialize enrichment (ptr because of circular dependency)
        mEnrichment = new Enrichment(
                Enrichment_Method::USE_INTERPOLATION_CELL_BASIS,
                aBasisRank,
                aMeshIndex,
                mGeometryEngine->get_num_phases(),
                this,
                &mCutMesh,
                &mBackgroundMesh);

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

        MORIS_ERROR(mDecomposed,"Mesh needs to be decomposed prior to calling ghost penalization");

        MORIS_ERROR(!mGhost,"Ghost penalization has already been called");

        mGhostStabilization = new Ghost_Stabilization(this);

        mGhostStabilization->setup_ghost_stabilization();

        mGhost = true;
    }

    // ----------------------------------------------------------------------------------

    Ghost_Stabilization &
    Model::get_ghost_stabilization(moris::moris_index  aIndex)
    {
        MORIS_ERROR(mGhost,"Ghost has not been constructed on this model.");
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
        mMultigrid->build_basis_exodus_information(tName);
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

    Background_Mesh &
    Model::get_background_mesh()
    {
        return mBackgroundMesh;
    }

    // ----------------------------------------------------------------------------------

    Background_Mesh const &
    Model::get_background_mesh() const
    {
        return mBackgroundMesh;
    }

    // ----------------------------------------------------------------------------------

    moris::ge::Geometry_Engine*
    Model::get_geom_engine()
    {
        return mGeometryEngine;
    }

    // ----------------------------------------------------------------------------------
    // Tet 10 conversion Source code
    // ----------------------------------------------------------------------------------

    void
    Model::convert_mesh_tet4_to_tet10()
    {
        MORIS_ASSERT(0,"not currently working");
        mConvertedToTet10s = true;

        // start timing on this decomposition
        std::clock_t start = std::clock();
        //        convert_mesh_tet4_to_tet10_internal();

        if(moris::par_rank() == 0)
        {
            std::cout<<"Tet4 to Tet10 conversion completed in "<< (std::clock() - start) / (double)(CLOCKS_PER_SEC)<<" s."<<std::endl;
        }
    }


    //------------------------------------------------------------------------------

    void
    Model::construct_neighborhood()
    {
        // size the member data
        mElementToElement.resize(this->get_num_elements_total());
        mElementToElementSideOrds.resize(this->get_num_elements_total());
        mElementToElementNeighborSideOrds.resize(this->get_num_elements_total());

        // add uncut neighborhood to connectivity
        // keep track of boundaries with an uncut mesh next to a cut mesh
        moris::Cell<moris::Cell<moris_index>> tCutToUncutFace;
        this->construct_uncut_neighborhood(tCutToUncutFace);

        // since there are hanging nodes in HMR we need to do a little extra
        // to get the neighborhood correct
        if(mBackgroundMesh.get_mesh_data().get_mesh_type() == MeshType::HMR)
        {
            // create the simple relationship neighborhood between child meshes
            this->construct_cut_mesh_simple_neighborhood();  

            // create neighborhood between the portion of the mesh intersected and not intersected. Creates
            // connections between tetrahedral cells and hexahedral cells
            this->construct_cut_mesh_to_uncut_mesh_neighborhood(tCutToUncutFace);

            // Currently not doing anything. Intended to construct the neighborhood between intersected cells
            // on multiple refinement levels
            this->construct_complex_neighborhood();
        }
        else
        {
            // create the simple relationship neighborhood between child meshes
            this->construct_cut_mesh_simple_neighborhood();

            // create the link between uncut background cells and their neighboring children cells
            this->construct_cut_mesh_to_uncut_mesh_neighborhood(tCutToUncutFace);
        }

        //    print_neighborhood();
    }

    // ----------------------------------------------------------------------------------

    void
    Model::delete_neighborhood()
    {
        mElementToElement.resize(0);
    }

    // ----------------------------------------------------------------------------------

    void
    Model::construct_subphase_neighborhood()
    // {
    //     // get the interpolation mesh
    //     moris::mtk::Interpolation_Mesh &tInterpMesh = mBackgroundMesh.get_mesh_data();

    //     // allocate subphase to subphase connectivity
    //     mSubphaseToSubPhase = moris::Cell<moris::Cell<moris::moris_index>>(mCutMesh.get_num_subphases());
    //     mSubphaseToSubPhaseMySideOrds = moris::Cell<moris::Cell<moris::moris_index>>(mCutMesh.get_num_subphases());
    //     mSubphaseToSubPhaseNeighborSideOrds = moris::Cell<moris::Cell<moris::moris_index>>(mCutMesh.get_num_subphases());
    //     mTransitionNeighborCellLocation = moris::Cell<moris::Cell<moris::moris_index>>(mCutMesh.get_num_subphases());

    //     // non unique temporary data
    //     moris::Cell<moris::Cell<moris::moris_index>> tNonUniqueSubphaseToSubphase(mCutMesh.get_num_subphases());
    //     moris::Cell<moris::Cell<moris::moris_index>> tNonUniqueSubphaseToSubPhaseMySideOrds(mCutMesh.get_num_subphases());
    //     moris::Cell<moris::Cell<moris::moris_index>> tNonUniqueSubphaseToSubPhaseNeighborSideOrds(mCutMesh.get_num_subphases());
    //     moris::Cell<moris::Cell<moris::moris_index>> tNonUniqueTransitionLocation(mCutMesh.get_num_subphases());

    //     //iterate through facets
    //     moris::uint tNumFacets = tInterpMesh.get_num_entities(EntityRank::ELEMENT);
    //     for (moris::moris_index iC = 0; iC < (moris::moris_index)tNumFacets; iC++)
    //     {
    //         // current cell
    //         mtk::Cell const *tCurrentCell = &tInterpMesh.get_mtk_cell(iC);

    //         // get the cells attached to the facet
    //         Matrix<IndexMat> tCellToCellSideIndex = tInterpMesh.get_elements_connected_to_element_and_face_ind_loc_inds(iC);
    //         Matrix<IndexMat> tCellToCellSideOrd = tInterpMesh.get_elements_connected_to_element_and_face_ord_loc_inds(iC);

    //         // get the neighboring cells
    //         Cell<mtk::Cell const *> tCells(tCellToCellSideOrd.numel());
    //         tInterpMesh.get_mtk_cells(tCellToCellSideOrd.get_row(0), tCells);

    //         // iterate through neighbor
    //         for (moris::uint iN = 0; iN < tCellToCellSideOrd.n_cols(); iN++)
    //         {
    //             // facet ordinal shared for current neighbors
    //             moris_index tMyOrdinal = tCellToCellSideOrd(1, iN);
    //             moris_index tNeighborOrdinal = tCellToCellSideOrd(2, iN);
    //             moris_index tTransitionCellLocation = tCellToCellSideOrd(3, iN);

    //             // neighbor cell
    //             mtk::Cell const *tOtherCell = &tInterpMesh.get_mtk_cell(tCellToCellSideIndex(0, iN));

    //             // get the subphase indices attached to the facet which is connected to the current cell
    //             Cell<moris::moris_index> tMyCellSubphaseIndices(0);
    //             Cell<moris::moris_index> tMyCellSubphaseBulkIndices(0);
    //             Cell<moris::moris_index> tMyRepresentativeCellIndex(0);
    //             Cell<moris::moris_index> tMyRepresentativeCellSideOrdinal(0);
    //             this->collect_subphases_attached_to_facet_on_cell(tCurrentCell->get_index(), tMyOrdinal, tMyCellSubphaseIndices, tMyCellSubphaseBulkIndices,tMyRepresentativeCellIndex,tMyRepresentativeCellSideOrdinal);

    //             // get the subphase indices attached to the facet which is connected to the other cell in the neighborhood
    //             Cell<moris::moris_index> tNeighborSubphaseIndices(0);
    //             Cell<moris::moris_index> tNeighborSubphaseBulkIndices(0);
    //             Cell<moris::moris_index> tNeighborRepresentativeCellIndex(0);
    //             Cell<moris::moris_index> tNeighborRepresentativeCellSideOrdinal(0);
    //             this->collect_subphases_attached_to_facet_on_cell(tOtherCell->get_index(), tNeighborOrdinal, tNeighborSubphaseIndices, tNeighborSubphaseBulkIndices,tNeighborRepresentativeCellIndex,tNeighborRepresentativeCellSideOrdinal);

    //             // iterate over subphases and add to neighborhood
    //             for (moris::uint i = 0; i < tMyCellSubphaseIndices.size(); i++)
    //             {
    //                 moris_index tMyBulkIndex = tMyCellSubphaseBulkIndices(i);
    //                 moris_index tMySubphaseIndex = tMyCellSubphaseIndices(i);

    //                 for (moris::uint j = 0; j < tNeighborSubphaseIndices.size(); j++)
    //                 {
    //                     moris_index tNeighborBulkIndex = tNeighborSubphaseBulkIndices(j);
    //                     moris_index tNeighborSubphaseIndex = tNeighborSubphaseIndices(j);

    //                     if (tMyBulkIndex == tNeighborBulkIndex)
    //                     {
    //                         mSubphaseToSubPhase(tMySubphaseIndex).push_back(tNeighborSubphaseIndex);
    //                         mSubphaseToSubPhaseMySideOrds(tMySubphaseIndex).push_back(tMyOrdinal);
    //                         mSubphaseToSubPhaseNeighborSideOrds(tMySubphaseIndex).push_back(tNeighborOrdinal);
    //                         mTransitionNeighborCellLocation(tMySubphaseIndex).push_back(tTransitionCellLocation);

    //                         //                        tNonUniqueSubphaseToSubphase(tNeighborSubphaseIndex).push_back(tMySubphaseIndex);
    //                         //                        tNonUniqueSubphaseToSubPhaseMySideOrds(tNeighborSubphaseIndex).push_back(tNeighborOrdinal);
    //                         //                        tNonUniqueSubphaseToSubPhaseNeighborSideOrds(tNeighborSubphaseIndex).push_back(tMyOrdinal);
    //                         //                        tNonUniqueTransitionLocation(tNeighborSubphaseIndex).push_back(MORIS_INDEX_MAX);
    //                     }
    //                 }
    //             }
    //         }
    //     }

    //     std::cout<<"PRINTING SUBPHASE"<<std::endl;
    //     this->print_subphase_neighborhood();
    // }
    {
        // get the interpolation mesh
        moris::mtk::Interpolation_Mesh & tInterpMesh = mBackgroundMesh.get_mesh_data();

        // construct element to subphase index
        moris::Matrix<moris::IndexMat> tCellToSubphase = this->get_element_to_subphase();

        // allocate subphase to subphase connectivity
        mSubphaseToSubPhase                 = moris::Cell<moris::Cell<moris::moris_index>>(mCutMesh.get_num_subphases());
        mSubphaseToSubPhaseMySideOrds       = moris::Cell<moris::Cell<moris::moris_index>>(mCutMesh.get_num_subphases());
        mSubphaseToSubPhaseNeighborSideOrds = moris::Cell<moris::Cell<moris::moris_index>>(mCutMesh.get_num_subphases());
        mTransitionNeighborCellLocation     = moris::Cell<moris::Cell<moris::moris_index>>(mCutMesh.get_num_subphases());

        // temporary map
        moris::Cell<std::unordered_map<moris_index,moris_index>> tSubphaseToSubphaseTracker(mCutMesh.get_num_subphases());

        // non unique temporary data
        moris::Cell<moris::Cell<moris::moris_index>> tNonUniqueSubphaseToSubphase(mCutMesh.get_num_subphases());
        moris::Cell<moris::Cell<moris::moris_index>> tNonUniqueSubphaseToSubPhaseMySideOrds(mCutMesh.get_num_subphases());
        moris::Cell<moris::Cell<moris::moris_index>> tNonUniqueSubphaseToSubPhaseNeighborSideOrds(mCutMesh.get_num_subphases());
        moris::Cell<moris::Cell<moris::moris_index>> tNonUniqueTransitionLocation(mCutMesh.get_num_subphases());

        //iterate through facets
        moris::uint tNumFacets = tInterpMesh.get_num_entities(EntityRank::ELEMENT);
        for(moris::moris_index iC = 0; iC < (moris::moris_index)tNumFacets; iC++)
        {
            // current cell
            mtk::Cell const * tCurrentCell = & tInterpMesh.get_mtk_cell(iC);

            // get the cells attached to the facet
            Matrix<IndexMat> tCellToCellSideIndex = tInterpMesh.get_elements_connected_to_element_and_face_ind_loc_inds(iC);
            Matrix<IndexMat> tCellToCellSideOrd  = tInterpMesh.get_elements_connected_to_element_and_face_ord_loc_inds(iC);

            // get the neighboring cells
            Cell<mtk::Cell const *> tCells(tCellToCellSideOrd.numel());
            tInterpMesh.get_mtk_cells(tCellToCellSideOrd.get_row(0),tCells);

            // iterate through neighbor
            for(moris::uint iN = 0; iN < tCellToCellSideOrd.n_cols(); iN++)
            {
                // facet ordinal shared for current neighbors
                moris_index tMyOrdinal  = tCellToCellSideOrd(1,iN);
                moris_index tNeighborOrdinal = tCellToCellSideOrd(2,iN);
                moris_index tTransitionCellLocation = tCellToCellSideOrd(3,iN);

                // neighbor cell
                mtk::Cell const * tOtherCell = & tInterpMesh.get_mtk_cell(tCellToCellSideIndex(0,iN));

                // get the subphase indices attached to the facet which is connected to the current cell
                Cell<moris::moris_index> tMyCellSubphaseIndices(0);
                Cell<moris::moris_index> tMyCellSubphaseBulkIndices(0);
                Cell<moris::moris_index> tMyRepresentativeCellIndex(0);
                Cell<moris::moris_index> tMyRepresentativeCellSideOrdinal(0);
                this->collect_subphases_attached_to_facet_on_cell( tCurrentCell->get_index(), tMyOrdinal, tMyCellSubphaseIndices, tMyCellSubphaseBulkIndices,tMyRepresentativeCellIndex,tMyRepresentativeCellSideOrdinal);

                // // get the subphase indices attached to the facet which is connected to the other cell in the neighborhood
                Cell<moris::moris_index> tNeighborSubphaseIndices(0);
                Cell<moris::moris_index> tNeighborSubphaseBulkIndices(0);
                Cell<moris::moris_index> tNeighborRepresentativeCellId(0);
                Cell<moris::moris_index> tNeighborRepresentativeCellSideOrdinal(0);
                this->collect_subphases_attached_to_facet_on_cell( tOtherCell->get_index(), tNeighborOrdinal, tNeighborSubphaseIndices, tNeighborSubphaseBulkIndices,tNeighborRepresentativeCellId,tNeighborRepresentativeCellSideOrdinal);

                std::unordered_map<moris_index,moris_index> tActiveNeighborSubphases;
                for(moris::uint iM = 0; iM < tNeighborSubphaseIndices.size(); iM++ )
                {
                    tActiveNeighborSubphases[tNeighborSubphaseIndices(iM)] =1;
                }

                moris::Cell<moris::Cell<moris::mtk::Cell*>> const & tIGCellNeighborhood = this->get_element_to_element();

                // iterate over subphases and add to neighborhood
                for(moris::uint i = 0; i < tMyCellSubphaseIndices.size(); i++)
                {
                    moris_index tMyBulkIndex     = tMyCellSubphaseBulkIndices(i);
                    moris_index tMySubphaseIndex = tMyCellSubphaseIndices(i);

                    // iterate through my representative cell neighborhood
                    moris_index tMyCellIndex = tMyRepresentativeCellIndex(i);

                    // iterate through neighborhood
                    for (moris::uint iNeighbor = 0; iNeighbor < tIGCellNeighborhood(tMyCellIndex).size(); iNeighbor++)
                    {
                        moris_index tNeighborIndex         = tIGCellNeighborhood(tMyCellIndex)(iNeighbor)->get_index();
                        moris_index tNeighborSubphaseIndex = tCellToSubphase(tNeighborIndex);
                        moris_index tNeighborBulkPhase     = mBackgroundMesh.get_element_phase_index(tNeighborIndex);

                        if(tNeighborBulkPhase == tMyBulkIndex && tNeighborSubphaseIndex != tMySubphaseIndex)
                        {
                            if(tActiveNeighborSubphases.find(tNeighborSubphaseIndex) !=tActiveNeighborSubphases.end())
                            {
                                 if(tSubphaseToSubphaseTracker(tMySubphaseIndex).find(tNeighborSubphaseIndex) == tSubphaseToSubphaseTracker(tMySubphaseIndex).end() )
                                {
                                    mSubphaseToSubPhase(tMySubphaseIndex).push_back(tNeighborSubphaseIndex);
                                    mSubphaseToSubPhaseMySideOrds(tMySubphaseIndex).push_back(tMyOrdinal);
                                    mSubphaseToSubPhaseNeighborSideOrds(tMySubphaseIndex).push_back(tNeighborOrdinal);
                                    mTransitionNeighborCellLocation(tMySubphaseIndex).push_back(tTransitionCellLocation);
                                    tSubphaseToSubphaseTracker(tMySubphaseIndex)[tNeighborSubphaseIndex] = 1;
                                }
                            }
                        }

                    }
                  }
            }
        }

    // std::cout<<"PRINTING SUBPHASE"<<std::endl;
    // this->print_subphase_neighborhood();
    }

    // ----------------------------------------------------------------------------------

    void
    Model::collect_subphases_attached_to_facet_on_cell(
            moris::moris_index         aCellIndex,
            moris::moris_index         aFacetOrdinal,
            Cell<moris::moris_index> & aCellSubphaseIndices,
            Cell<moris::moris_index> & aCellSubphaseBulkIndices,
            Cell<moris::moris_index> & aRepresentativeCellInd,
            Cell<moris::moris_index> & aRepresentativeCellSideOrdinal)
    {
        aRepresentativeCellInd.clear();
        aRepresentativeCellSideOrdinal.clear();

        Matrix<IndexMat> tCellFacets = mBackgroundMesh.get_mesh_data().get_entity_connected_to_entity_loc_inds(aCellIndex,EntityRank::ELEMENT, mBackgroundMesh.get_mesh_data().get_facet_rank());

        moris_index tFacetIndex = tCellFacets(aFacetOrdinal);

        // set pointer to child mesh if it has children
        if(mBackgroundMesh.entity_has_children(aCellIndex, EntityRank::ELEMENT))
        {
            // get the child mesh ptr
            moris::moris_index tCMIndex = mBackgroundMesh.child_mesh_index(aCellIndex,EntityRank::ELEMENT);
            Child_Mesh const * tCMCell = & mCutMesh.get_child_mesh(tCMIndex);


            // get subphases attached to facet
            Cell<moris::moris_index> tCellCMSubphaseIndices;
            tCMCell->get_subphases_attached_to_facet(tFacetIndex, tCellCMSubphaseIndices, aRepresentativeCellInd, aRepresentativeCellSideOrdinal);

            // reference to subphase indices and bulkphases of subphases
            Cell<moris::moris_index> const & tCMSubphaseBulkIndices = tCMCell->get_subphase_bin_bulk_phase();
            Cell<moris::moris_index> const & tCMSubphaseIndices     = tCMCell->get_subphase_indices();

            // put the information in processor local indices with bulk phase
            aCellSubphaseIndices.resize(tCellCMSubphaseIndices.size());
            aCellSubphaseBulkIndices.resize(tCellCMSubphaseIndices.size());
            

            for(moris::uint i = 0; i < tCellCMSubphaseIndices.size(); i++)
            {
                aCellSubphaseBulkIndices(i) = tCMSubphaseBulkIndices(tCellCMSubphaseIndices(i));
                aCellSubphaseIndices(i) = tCMSubphaseIndices(tCellCMSubphaseIndices(i));
            }
        }
        else
        {
            // get the cell subphase indices
            aCellSubphaseBulkIndices = {{mBackgroundMesh.get_element_phase_index(aCellIndex)}};
            aCellSubphaseIndices     = {{aCellIndex}};
            aRepresentativeCellInd    = {{aCellIndex}};
            aRepresentativeCellSideOrdinal.push_back(tFacetIndex);
        }
    }

    // ----------------------------------------------------------------------------------

    bool
    Model::subphase_is_in_child_mesh(moris_index aSubphaseIndex)
    {
        if(aSubphaseIndex >= (moris_index)mBackgroundMesh.get_mesh_data().get_num_entities(EntityRank::ELEMENT))
        {
            return true;
        }

        else if (mBackgroundMesh.entity_has_children(aSubphaseIndex,EntityRank::ELEMENT))
        {
            return true;
        }

        return false;
    }

    // ----------------------------------------------------------------------------------

    void
    Model::construct_cut_mesh_simple_neighborhood()
    {
        // Collect element to node of child mesh and create a full element to element graph for child mesh
        Matrix<IndexMat> tCMElementToNode = mCutMesh.get_full_element_to_node_loc_inds();

        // generate connectivities (face then element conn) we do not keep faces around though
        // face connectivity
        CellTopology tCellTopo = mCutMesh.get_child_element_topology();
        Matrix<IndexMat> tElementToFace;
        Matrix<IndexMat> tFaceToNode;
        Matrix<IndexMat> tNodeToFace;
        Matrix<IndexMat> tFaceToElement;

        moris::mtk::Cell_Info_Factory tFactory;
        moris::mtk::Cell_Info* tCellInfo = tFactory.create_cell_info(tCellTopo);
        xtk::create_faces_from_element_to_node(tCellInfo, mBackgroundMesh.get_num_entities(EntityRank::NODE), tCMElementToNode, tElementToFace, tFaceToNode, tNodeToFace, tFaceToElement);

        // remove this data right away since it is not needed
        tNodeToFace.set_size(0,0);

        // element connectivity
        moris::size_t tNumFacePerElem = tCellInfo->get_num_facets();

        // generate the element to element connectivity. this only captures the easy neighborhood relationships
        // where we consider elements neighbors if they share a full face
        moris::Matrix< moris::IndexMat > tElementToElementSharedFacet;
        Matrix<IndexMat> tElementToElementMat = generate_element_to_element(tFaceToElement, mCutMesh.get_num_entities(EntityRank::ELEMENT), tNumFacePerElem, MORIS_INDEX_MAX,tElementToElementSharedFacet);

        moris::Matrix< moris::IndexMat > tCMElementInds = mCutMesh.get_all_element_inds();

        moris::Cell<moris::mtk::Cell*> tCMCellPtrs(tCMElementInds.numel());
        for(moris::uint iC = 0; iC<tCMCellPtrs.size(); iC++ )
        {
            tCMCellPtrs(iC) = &mBackgroundMesh.get_mtk_cell(tCMElementInds(iC));
        }

        // add to member data
        for(moris::uint iC = 0; iC<tElementToElementMat.n_rows(); iC++ )
        {
            for(moris::uint iN = 0; iN< tElementToElementMat.n_cols(); iN++)
            {
                if(tElementToElementMat(iC,iN) == MORIS_INDEX_MAX)
                {
                    continue;
                }

                moris_index tFacetIndex         = tElementToElementSharedFacet(iC,iN);
                moris_index tNeighborCellIndex  = tElementToElementMat(iC,iN);


                // figure out my side ordinal
                moris_index tMysideOrd   = MORIS_INDEX_MAX;
                for(moris::uint iSO = 0; iSO < tElementToFace.n_cols(); iSO++)
                {
                    if(tElementToFace(iC,iSO) == (moris_index) tFacetIndex)
                    {
                        tMysideOrd = iSO;
                    }
                }

                // figure out my neighbors side ordinal
                moris_index tNeighborSideOrd = MORIS_INDEX_MAX;
                for(moris::uint iSO = 0; iSO < tElementToFace.n_cols(); iSO++)
                {
                    if(tElementToFace(tNeighborCellIndex,iSO) == (moris_index) tFacetIndex)
                    {
                        tNeighborSideOrd = iSO;
                    }
                }

                MORIS_ERROR(tNeighborSideOrd!= MORIS_INDEX_MAX,"Neighbor Side Ord Not Found");
                MORIS_ERROR(tMysideOrd!= MORIS_INDEX_MAX,"My Cell Side Ord Not Found");
                mElementToElement(tCMElementInds(iC)).push_back(tCMCellPtrs(tNeighborCellIndex));
                mElementToElementSideOrds(tCMElementInds(iC)).push_back(tMysideOrd);
                mElementToElementNeighborSideOrds(tCMElementInds(iC)).push_back(tNeighborSideOrd);

            }
        }

        delete tCellInfo;
    }

    // ----------------------------------------------------------------------------------

    void
    Model::construct_complex_neighborhood()
    {

    }

    // ----------------------------------------------------------------------------------

    void
    Model::construct_cut_mesh_to_uncut_mesh_neighborhood(moris::Cell<moris::Cell<moris_index>> & aCutToUncutFace)
    {
        // estimate maximum number of elements on face
        const uint tMaxElemOnFace = 100;

        // matrices used throughout routine
        moris::Matrix< moris::IdMat >    tChildElemsIdsOnFace(1,tMaxElemOnFace);
        moris::Matrix< moris::IndexMat > tChildElemsCMIndOnFace(1,tMaxElemOnFace);
        moris::Matrix< moris::IndexMat > tChildElemOnFaceOrdinal(1,tMaxElemOnFace);

        // iterate through the cut to uncut relationships
        for(moris::uint iR = 0; iR < aCutToUncutFace.size(); iR++)
        {
            // get data for readability from aCutToUncutFace
            moris_index tCellUnCutInd = aCutToUncutFace(iR)(0);
            moris_index tCellCutInd   = aCutToUncutFace(iR)(1);
            moris_index tFaceIndex    = aCutToUncutFace(iR)(2);

            // uncut cell
            moris::mtk::Cell* tCellUnCut = &mBackgroundMesh.get_mtk_cell(tCellUnCutInd);

            Child_Mesh &  tChildMesh =
                    mCutMesh.get_child_mesh(mBackgroundMesh.child_mesh_index(tCellCutInd,EntityRank::ELEMENT));

            Matrix<IndexMat> const & tChildCellInds = tChildMesh.get_element_inds();

            // define variable for actual number of child elements on face
            uint tNumberOfChildElemsOnFace;

            tChildMesh.get_child_elements_connected_to_parent_facet(
                    tFaceIndex,
                    tNumberOfChildElemsOnFace,
                    tChildElemsIdsOnFace,
                    tChildElemsCMIndOnFace,
                    tChildElemOnFaceOrdinal);

            // get child cells and add cut cell to neighborhood, cut cells to neighborhood of uncut
            for(moris::uint  i = 0; i < tNumberOfChildElemsOnFace; i++ )
            {
                mElementToElement(tChildCellInds(tChildElemsCMIndOnFace(i))).push_back(tCellUnCut);

                mElementToElement(tCellUnCutInd).
                        push_back(&mBackgroundMesh.get_mtk_cell(tChildCellInds(tChildElemsCMIndOnFace(i))));
            }
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Model::construct_uncut_neighborhood(moris::Cell<moris::Cell<moris_index>> & aCutToUncutFace)
    {
        moris::mtk::Interpolation_Mesh & tInterpMesh = mBackgroundMesh.get_mesh_data();
        moris::uint tNumCells                        = tInterpMesh.get_num_entities(EntityRank::ELEMENT);

        // iterate through background cells
        for(moris::uint iC = 0; iC < tNumCells; iC++)
        {
            // if current cell has no children, i.e. is uncut
            if(!mBackgroundMesh.entity_has_children((moris_index)iC,EntityRank::ELEMENT))
            {
                // get the neighbors
                Matrix<IndexMat> tElementNeighors = tInterpMesh.get_elements_connected_to_element_and_face_ind_loc_inds(iC);

                // iterate through neighbors
                for(moris::uint iN = 0; iN<tElementNeighors.n_cols(); iN++ )
                {
                    // if the neighbor has no children, i.e., is uncut
                    if(!mBackgroundMesh.entity_has_children(tElementNeighors(0,iN),EntityRank::ELEMENT))
                    {
                        // add neighbor cell to list of neighbors for current cell
                        mElementToElement(iC).push_back(&mBackgroundMesh.get_mtk_cell(tElementNeighors(0,iN)));
                    }

                    // mark as a cut to uncut boundary
                    else
                    {
                        // store
                        aCutToUncutFace.push_back( {(moris_index)iC , tElementNeighors(0,iN), tElementNeighors(1,iN)} );
                    }
                }
            }
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Model::print_neighborhood()
    {
        for(moris::uint iE = 0; iE < mElementToElement.size(); iE++)
        {
            moris::mtk::Cell & tCell = mBackgroundMesh.get_mtk_cell(iE);

            std::cout<<"Element Id: "<<std::setw(8)<<tCell.get_id()<<" | ";

            for(moris::uint iN = 0; iN < mElementToElement(iE).size(); iN++)
            {
                std::cout<<std::setw(8)<<mElementToElement(iE)(iN)->get_id();
            }
            std::cout<<std::endl;
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Model::print_cells()
    {
        for(moris::uint  i = 0; i < this->get_num_elements_total(); i++ )
        {
            mtk::Cell & tCell = mBackgroundMesh.get_mtk_cell((moris_index)i);

            moris::Cell<moris::mtk::Vertex *> tVertices = tCell.get_vertex_pointers();

            std::cout<<"Cell Id: "<<std::setw(8)<<tCell.get_id();
            std::cout<<" | Cell Index: "<<std::setw(8)<<tCell.get_index();
            std::cout<<" | Phase: "<<std::setw(8)<<mBackgroundMesh.get_element_phase_index(i) << " | Verts: ";
            for(moris::uint j = 0; j < tVertices.size(); j++)
            {
                std::cout<<std::setw(8)<<tVertices(j)->get_id();
            }

            std::cout<<std::endl;
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Model::print_interface_vertices()
    {
        mBackgroundMesh.print_interface_node_flags();
    }

    //------------------------------------------------------------------------------

    void
    Model::print_subphase_neighborhood()
    {

        std::cout<<"Subphases"<<std::endl;
        for(moris::uint iC = 0; iC<mSubphaseToSubPhase.size(); iC++ )
        {
            std::cout<<std::setw(6)<<iC<<" | ";

            for(moris::uint iN = 0; iN< mSubphaseToSubPhase(iC).size(); iN++)
            {
                std::cout<<std::setw(6)<<mSubphaseToSubPhase(iC)(iN);
            }
            std::cout<<std::endl;
        }

        std::cout<<"Subphases My Side Ordinals"<<std::endl;
        for(moris::uint iC = 0; iC<mSubphaseToSubPhaseMySideOrds.size(); iC++ )
        {
            std::cout<<std::setw(6)<<iC<<" | ";

            for(moris::uint iN = 0; iN< mSubphaseToSubPhaseMySideOrds(iC).size(); iN++)
            {
                std::cout<<std::setw(6)<<mSubphaseToSubPhaseMySideOrds(iC)(iN);
            }
            std::cout<<std::endl;
        }

        std::cout<<"Subphases Neighbor Side Ordinals"<<std::endl;
        for(moris::uint iC = 0; iC<mSubphaseToSubPhaseNeighborSideOrds.size(); iC++ )
        {
            std::cout<<std::setw(6)<<iC<<" | ";

            for(moris::uint iN = 0; iN< mSubphaseToSubPhaseNeighborSideOrds(iC).size(); iN++)
            {
                std::cout<<std::setw(6)<<mSubphaseToSubPhaseNeighborSideOrds(iC)(iN);
            }
            std::cout<<std::endl;
        }

        if(mBackgroundMesh.get_mesh_data().get_mesh_type() == MeshType::HMR)
        {
            std::cout<<"Transition Neighbor Locations"<<std::endl;
            for(moris::uint iC = 0; iC<mTransitionNeighborCellLocation.size(); iC++ )
            {
                std::cout<<std::setw(6)<<iC<<" | ";

                for(moris::uint iN = 0; iN< mTransitionNeighborCellLocation(iC).size(); iN++)
                {
                    std::cout<<std::setw(12)<<mTransitionNeighborCellLocation(iC)(iN);
                }
                std::cout<<std::endl;
            }
        }
    }

    //------------------------------------------------------------------------------

    moris::mtk::Integration_Mesh*
    Model::get_output_mesh(Output_Options const & aOutputOptions)

    {
        // start timing on this output
        std::clock_t start = std::clock();

        // create the output mtk mesh
        moris::mtk::Integration_Mesh* tOutputMesh = construct_output_mesh(aOutputOptions);

        if(moris::par_rank() == 0  && mVerbose)
        {
            std::cout<<"XTK: Mesh output completed in "<< (std::clock() - start) / (double)(CLOCKS_PER_SEC)<<" s."<<std::endl;
        }
        return tOutputMesh;
    }

    //------------------------------------------------------------------------------

    moris::uint
    Model::get_spatial_dim() const
    {
        return mModelDimension;
    }

    //------------------------------------------------------------------------------

    moris::uint
    Model::get_num_elements_total()
    {
        MORIS_ASSERT(mDecomposed,"Prior to using get_num_elements, the decomposition process must be finished");

        return mBackgroundMesh.get_num_entities(EntityRank::ELEMENT);
    }

    //------------------------------------------------------------------------------

    moris::uint
    Model::get_num_elements_unzipped()
    {
        return this->get_num_elements_total() - mCutMesh.get_num_child_meshes();
    }

    //------------------------------------------------------------------------------

    moris_index
    Model::get_cell_xtk_index(moris_id aCellId)
    {
        auto tIter = mCellGlbToLocalMap.find(aCellId);
        MORIS_ASSERT(tIter != mCellGlbToLocalMap.end(),"Id not in map");
        return tIter->second;
    }

    //------------------------------------------------------------------------------

    moris::Cell<moris::Cell<moris_index>>  const &
    Model::get_subphase_to_subphase()
    {
        return mSubphaseToSubPhase;
    }

    //------------------------------------------------------------------------------

    moris::Cell<moris::Cell<moris_index>>  const &
    Model::get_subphase_to_subphase_my_side_ords()
    {
        return mSubphaseToSubPhaseMySideOrds;
    }

    //------------------------------------------------------------------------------

    moris::Cell<moris::Cell<moris_index>>  const &
    Model::get_subphase_to_subphase_transition_loc()
    {
        return mTransitionNeighborCellLocation;
    }

    //------------------------------------------------------------------------------

    moris::Cell<moris::Cell<moris_index>>  const &
    Model::get_subphase_to_subphase_neighbor_side_ords()
    {
        return mSubphaseToSubPhaseNeighborSideOrds;
    }

    //------------------------------------------------------------------------------

    std::shared_ptr< Multigrid >
    Model::get_multigrid_ptr()
    {
        return mMultigrid;
    }

    //------------------------------------------------------------------------------

    moris::Matrix<moris::IndexMat>
    Model::get_element_to_subphase()
    {
        // child mesh subphases
        moris::uint tNumElem = this->get_num_elements_total();
        Matrix<IndexMat> tSubphase(1,tNumElem);
        mCutMesh.populate_subphase_vector(tSubphase);

        // populate the non intersected background cells
        for(moris::uint i = 0; i < mBackgroundMesh.get_mesh_data().get_num_entities(EntityRank::ELEMENT); i++)
        {
            tSubphase(i) = i;
        }

        return tSubphase;
    }

    moris::Matrix<moris::IndexMat>
    Model::get_num_subphase_neighbors()
    {
        moris::Cell<moris::Cell<moris_index>>  const & tSubPhaseToSubphase = this->get_subphase_to_subphase();
        moris::Matrix<moris::IndexMat> tSubphaseNumNeighbors( 1, tSubPhaseToSubphase.size());
        for (size_t iSP = 0; iSP < tSubPhaseToSubphase.size(); iSP++)
        {
            tSubphaseNumNeighbors(iSP) = tSubPhaseToSubphase(iSP).size();
        }
        return tSubphaseNumNeighbors;
    }

    //------------------------------------------------------------------------------

    moris_id
    Model::get_subphase_id(moris_id aSubphaseIndex)
    {
        if(this->subphase_is_in_child_mesh(aSubphaseIndex))
        {
            return mCutMesh.get_subphase_id(aSubphaseIndex);
        }
        else
        {
            return mBackgroundMesh.get_glb_entity_id_from_entity_loc_index(aSubphaseIndex,EntityRank::ELEMENT);
        }
    }

    // ----------------------------------------------------------------------------------

    moris_index
    Model::get_subphase_index(moris_id aSubphaseId)
    {
        auto tIter = mGlobalToLocalSubphaseMap.find(aSubphaseId);

        MORIS_ASSERT(tIter != mGlobalToLocalSubphaseMap.end(),"Subphase id not in map");

        return tIter->second;
    }

    //------------------------------------------------------------------------------


    moris::mtk::Integration_Mesh*
    Model::construct_output_mesh( Output_Options const & aOutputOptions )
    {

        // start timing on this decomposition
        std::clock_t start = std::clock();

        // Get mesh information ready for outputting
        // Package element to Node Connectivity
        moris::uint tSpatialDim = mBackgroundMesh.get_mesh_data().get_spatial_dim();

        // Children element nodes connected to elements
        moris::Cell<moris::Matrix<moris::IdMat>>  tElementToNodeChildrenByPhase = mCutMesh.get_full_element_to_node_by_phase_glob_ids(mGeometryEngine->get_num_bulk_phase(),mBackgroundMesh.get_mesh_data());

        // Child element ids
        moris::Cell<moris::Matrix<moris::IdMat>>  tChildElementsByPhase = mCutMesh.get_child_elements_by_phase(mGeometryEngine->get_num_bulk_phase(),mBackgroundMesh.get_mesh_data());

        // Parent elements without children
        Cell<moris::Matrix<moris::IdMat>>  tElementNoChildrenIdsByPhase = mBackgroundMesh.get_all_non_intersected_elements_by_phase(mGeometryEngine->get_num_bulk_phase());

        // Connectivity of parent elements without children
        Cell<moris::Matrix<moris::IdMat>>  tElementToNodeNoChildrenByPhase = mBackgroundMesh.get_non_intersected_element_to_node_by_phase(mGeometryEngine->get_num_bulk_phase());

        // Node map  of nodes in the phase we are interested in
        moris::Matrix<moris::IndexMat> tOutputtedNodeInds;
        moris::Matrix<moris::IdMat>  tLocalToGlobalNodeMap = this->get_node_map_restricted_to_output_phases(aOutputOptions,tOutputtedNodeInds);

        // All node coordinates
        moris::Matrix<moris::DDRMat> tNodeCoordinates = mBackgroundMesh.get_selected_node_coordinates_loc_inds(tOutputtedNodeInds);

        // Number of bulk phases
        uint tNumBulkPhases = mGeometryEngine->get_num_phases();

        // Get non-interescted parent elements by phase
        Cell<moris::Matrix<moris::IdMat>> tNoChildElementsByPhase = mBackgroundMesh.get_all_non_intersected_elements_by_phase(tNumBulkPhases);

        // combination of the elements by phase (if specified)
        Cell<moris::Matrix<moris::IdMat>> tCombinedElementsByPhase(tNoChildElementsByPhase.size());
        if(!aOutputOptions.mSeparateInterfaceBlock)
        {
            MORIS_ASSERT(mBackgroundMesh.get_parent_cell_topology() == CellTopology::TET4,
                    " Combining the interface block w/ non-interface block only valid on tet background mesh");

            tCombinedElementsByPhase = combine_interface_and_non_interface_blocks(tChildElementsByPhase,tNoChildElementsByPhase);
        }

        // Interface nodes
        Cell<moris::Matrix<moris::IndexMat>> tInterfaceNodes = mBackgroundMesh.get_interface_nodes_glob_ids();

        // number of phases being output
        moris::uint tNumPhasesOutput = get_num_phases_to_output(aOutputOptions);

        // Set up field data structure for MTK input
        moris::mtk::MtkFieldsInfo tFieldsInfo;

        // External Fields - real cell fields
        uint tNumExtRealCellScalarFields = aOutputOptions.mRealElementExternalFieldNames.size();
        moris::Cell<moris::mtk::Scalar_Field_Info<DDRMat>> tExternalRealCellScalarFields(tNumExtRealCellScalarFields);
        for(uint i = 0; i<tNumExtRealCellScalarFields; i++)
        {
            tExternalRealCellScalarFields(i).set_field_name(aOutputOptions.mRealElementExternalFieldNames(i));
            tExternalRealCellScalarFields(i).set_field_entity_rank(moris::EntityRank::ELEMENT);
            add_field_for_mesh_input(&tExternalRealCellScalarFields(i),tFieldsInfo);
        }

        // External Fields - real vertex fields
        uint tNumExtRealVertexScalarFields = aOutputOptions.mRealNodeExternalFieldNames.size();
        moris::Cell<moris::mtk::Scalar_Field_Info<DDRMat>> tExternalRealVertexScalarFields(tNumExtRealVertexScalarFields);
        for(uint i = 0; i<tNumExtRealVertexScalarFields; i++)
        {
            tExternalRealVertexScalarFields(i).set_field_name(aOutputOptions.mRealNodeExternalFieldNames(i));
            tExternalRealVertexScalarFields(i).set_field_entity_rank(moris::EntityRank::NODE);
            add_field_for_mesh_input(&tExternalRealVertexScalarFields(i),tFieldsInfo);
        }

        // sensitivity fields
        moris::Cell<moris::Matrix<DDRMat>> adxdpData;
        moris::Cell<std::string>           adxdpNames;
        moris::Cell<moris::Matrix<DDRMat>> aDesVars;
        moris::Cell<std::string>           aDesVarsName;
        moris::Matrix<moris::DDRMat>       aNumDesVars;
        std::string                        aNumDesVarsName;
        moris::Cell<moris::mtk::Scalar_Field_Info<DDRMat>> tdxdpDataFields;
        moris::Cell<moris::mtk::Scalar_Field_Info<DDRMat>> tDesVarFields;
        moris::Cell<moris::mtk::Scalar_Field_Info<DDRMat>> tNumDesVarsField;

        //TODO: implement node owner (currently set to owned by this proc)
        //    moris::Matrix<moris::IdMat> tNodeOwner(1,tOutputtedNodeInds.numel(),moris::par_rank());

        moris::Matrix<moris::IdMat> tNodeOwner = mBackgroundMesh.get_vertices_owner(tOutputtedNodeInds);

        // Set up mesh sets
        // Initialize Sets information structure
        moris::mtk::MtkSetsInfo tMtkMeshSets;

        //
        moris::uint tNumBlocksPerPhase = 2;
        if(!aOutputOptions.mSeparateInterfaceBlock)
        {
            MORIS_ASSERT(mBackgroundMesh.get_parent_cell_topology() == CellTopology::TET4,
                    " Combining the interface block w/ non-interface block only valid on tet background mesh");
            tNumBlocksPerPhase = 1;
        }

        // Setup block sets
        Cell<moris::mtk::MtkBlockSetInfo> tBlockSets(tNumPhasesOutput*tNumBlocksPerPhase);
        uint tCount= 0;

        for(uint i = 0; i <tNumBulkPhases; i++)
        {
            if(aOutputOptions.output_phase(i) && aOutputOptions.mSeparateInterfaceBlock)
            {
                // Children of material phase i
                tBlockSets(tCount).mCellIdsInSet = &tChildElementsByPhase(i);
                tBlockSets(tCount).mBlockSetName = "child_"+std::to_string(i);
                tBlockSets(tCount).mBlockSetTopo = mCutMesh.get_child_element_topology();

                tMtkMeshSets.add_block_set(&tBlockSets(tCount));
                tCount++;

                // Children of material phase i
                tBlockSets(tCount).mCellIdsInSet = &tNoChildElementsByPhase(i);
                tBlockSets(tCount).mBlockSetName = "parent_"+std::to_string(i);
                tBlockSets(tCount).mBlockSetTopo = mBackgroundMesh.get_parent_cell_topology();

                tMtkMeshSets.add_block_set(&tBlockSets(tCount));
                tCount++;
            }

            else if(aOutputOptions.output_phase(i) && !aOutputOptions.mSeparateInterfaceBlock)
            {
                // Children of material phase i
                tBlockSets(tCount).mCellIdsInSet = &tCombinedElementsByPhase(i);
                tBlockSets(tCount).mBlockSetName = "phase_"+std::to_string(i);
                tBlockSets(tCount).mBlockSetTopo = mBackgroundMesh.get_parent_cell_topology();

                tMtkMeshSets.add_block_set(&tBlockSets(tCount));
                tCount++;
            }
        }

        // Interface elements
        moris::Matrix<moris::IndexMat> tInterfaceElements(0,0);
        moris::Matrix<moris::IndexMat> tInterfaceElementIds(0,0);
        moris::mtk::MtkBlockSetInfo tUnzippedInterfaceBlockSet;
        if(mUnzipped && aOutputOptions.mHaveInterface)
        {
            // get the interface elements local node index element connectivity
            tInterfaceElements = mCutMesh.get_extracted_interface_elements_loc_inds();

            mBackgroundMesh.convert_loc_entity_ind_to_glb_entity_ids(EntityRank::NODE,tInterfaceElements);

            tInterfaceElementIds = mCutMesh.get_interface_element_ids();

            tUnzippedInterfaceBlockSet.mCellIdsInSet = &tInterfaceElementIds;
            tUnzippedInterfaceBlockSet.mBlockSetName = "interface";
            tUnzippedInterfaceBlockSet.mBlockSetTopo = CellTopology::PRISM6;
            tMtkMeshSets.add_block_set(&tUnzippedInterfaceBlockSet);
        }

        // propogate background mesh node sets
        moris::Cell<moris::Matrix<IndexMat>> tBackgroundNodeSetData;
        moris::Cell<moris::mtk::MtkNodeSetInfo> tBackgroundNodeSets;
        if(aOutputOptions.mAddNodeSets)
        {
            tBackgroundNodeSets = propogate_background_node_sets(tBackgroundNodeSetData,aOutputOptions);

            for(moris::uint i = 0; i<tBackgroundNodeSets.size(); i++)
            {
                tMtkMeshSets.add_node_set(&tBackgroundNodeSets(i));
            }
        }

        moris::Cell<moris::mtk::MtkNodeSetInfo> tInterfaceNodeSets(tInterfaceNodes.size());
        if(aOutputOptions.mHaveInterface)
        {

            for(uint i = 0; i<tInterfaceNodes.size(); i++)
            {
                tInterfaceNodeSets(i).mNodeIds     = &tInterfaceNodes(i);
                tInterfaceNodeSets(i).mNodeSetName = "inodes_" +std::to_string(i) ;
                tMtkMeshSets.add_node_set(&tInterfaceNodeSets(i));
            }
        }

        // Get the packaged interface side sets from the cut mesh
        Cell<moris::Matrix<moris::IdMat>> tInterfaceElemIdandSideOrd;
        Cell<std::string> tInterfaceNames;
        Cell<moris::mtk::MtkSideSetInfo> tInterfaceSideSets;
        if(aOutputOptions.mHaveInterface)
        {
            this->setup_interface_single_side_sets(aOutputOptions,tInterfaceElemIdandSideOrd,tInterfaceNames);

            tInterfaceSideSets.resize(tInterfaceElemIdandSideOrd.size());
            for(moris::uint i = 0; i < tInterfaceElemIdandSideOrd.size(); i++)
            {
                tInterfaceSideSets(i).mElemIdsAndSideOrds = &tInterfaceElemIdandSideOrd(i);
                tInterfaceSideSets(i).mSideSetName        = tInterfaceNames(i);
                tMtkMeshSets.add_side_set(&tInterfaceSideSets(i));
            }
        }

        // propogate side sets from background mesh
        moris::Cell<moris::Matrix<IndexMat>> tSideSetData;
        moris::Cell<moris::mtk::MtkSideSetInfo> tBackgroundSideSets;
        if(aOutputOptions.mAddSideSets)
        {
            // collect information about background side set
            tBackgroundSideSets = this->propogate_background_side_sets(tSideSetData,aOutputOptions);

            // add to mesh input structure
            for(moris::uint i = 0; i<tBackgroundSideSets.size(); i++)
            {
                tMtkMeshSets.add_side_set(&tBackgroundSideSets(i));
            }
        }

        // Mesh data input structure (with multi element type mesh)
        moris::mtk::MtkMeshData tMeshDataInput;

        moris::uint tNumElemTypes = tNumPhasesOutput*2;
        if(mUnzipped)
        {
            tNumElemTypes = tNumElemTypes + 1;
        }

        tMeshDataInput.ElemConn             = moris::Cell<moris::Matrix<IdMat>*>(tNumElemTypes);
        tMeshDataInput.LocaltoGlobalElemMap = moris::Cell<moris::Matrix<IdMat>*>(tNumElemTypes);
        tMeshDataInput.CellTopology         = moris::Cell<enum CellTopology>(tNumElemTypes,CellTopology::INVALID);

        tMeshDataInput.Verbose                 = mVerbose;
        tMeshDataInput.SpatialDim              = &tSpatialDim;

        tCount = 0;
        for(moris::uint  i = 0 ; i <mGeometryEngine->get_num_bulk_phase(); i++)
        {
            if(aOutputOptions.output_phase(i))
            {
                tMeshDataInput.ElemConn(tCount)             = &tElementToNodeChildrenByPhase(i);
                tMeshDataInput.LocaltoGlobalElemMap(tCount) = &tChildElementsByPhase(i);
                tMeshDataInput.CellTopology(tCount)         = mCutMesh.get_child_element_topology();
                tCount++;
                tMeshDataInput.ElemConn(tCount)             = &tElementToNodeNoChildrenByPhase(i);
                tMeshDataInput.LocaltoGlobalElemMap(tCount) = &tElementNoChildrenIdsByPhase(i);
                tMeshDataInput.CellTopology(tCount)         = mBackgroundMesh.get_parent_cell_topology();
                tCount++;
            }
        }

        tMeshDataInput.NodeCoords              = &tNodeCoordinates;
        tMeshDataInput.FieldsInfo              = &tFieldsInfo;
        tMeshDataInput.LocaltoGlobalNodeMap    = &tLocalToGlobalNodeMap;
        tMeshDataInput.SetsInfo                = &tMtkMeshSets;
        tMeshDataInput.MarkNoBlockForIO        = false;
        tMeshDataInput.CreateAllEdgesAndFaces  = true;
        tMeshDataInput.AutoAuraOptionInSTK     = false;

        //Add clustering information
        moris::mtk::Cell_Cluster_Input tCellClusterInput;                // Cell clusters
        moris::Cell<Matrix<IdMat>>     tClusterCellIds;                  // Cell cluster Ids
        moris::mtk::Side_Cluster_Input tInterfaceSideClusterInput;       // Side clusters
        moris::Cell<Matrix<IdMat>>     tInterfaceCellIdsandSideOrds;     // side cluster ids and side ordinals
        moris::Cell<Matrix<DDRMat>>    tInterfaceSideClusterParamCoords; // side cluster vertex parametric coordinates
        //
        if(aOutputOptions.mAddClusters)
        {
            // cell clustering
            this->setup_cell_clusters_for_output(tCellClusterInput,aOutputOptions,tClusterCellIds);

            tMeshDataInput.CellClusterInput = &tCellClusterInput;

            MORIS_ASSERT(mGeometryEngine->get_num_geometries() == 1,"This has not been setup for multi geometry problems.");

            //fixme: support changes to interface side
            this->setup_interface_side_cluster(tInterfaceNames(0), tInterfaceSideClusterInput, aOutputOptions, tInterfaceCellIdsandSideOrds, tInterfaceSideClusterParamCoords);

            tMeshDataInput.SideClusterInput = &tInterfaceSideClusterInput;
        }

        // Interface elements
        if(mUnzipped)
        {
            tMeshDataInput.ElemConn(tNumElemTypes-1)             = &tInterfaceElements;
            tMeshDataInput.LocaltoGlobalElemMap(tNumElemTypes-1) = &tInterfaceElementIds;
            tMeshDataInput.CellTopology(tNumElemTypes-1)         = CellTopology::PRISM6;
        }

        // add parallel information
        moris::mtk::Visualization_STK tVizTool;
        if(aOutputOptions.mAddParallelFields && par_size()>1)
        {
            moris::mtk::MtkFieldsInfo* tParFields = tVizTool.setup_parallel_cell_fields_for_declaration();
            tFieldsInfo.combine_fields_info(*tParFields);
        }

        if(moris::par_rank() == 0 && mVerbose)
        {
            std::cout<<"XTK: Mesh data setup completed in " <<(std::clock() - start) / (double)(CLOCKS_PER_SEC)<<" s."<<std::endl;
        }

        start = std::clock();

        // cast background mesh to an interpolation mesh and pass in
        moris::mtk::Integration_Mesh* tMeshData = nullptr;
        if(aOutputOptions.mAddClusters)
        {
            moris::mtk::Interpolation_Mesh* tInterpMesh = dynamic_cast<moris::mtk::Interpolation_Mesh*>(&mBackgroundMesh.get_mesh_data());

            tMeshData = moris::mtk::create_integration_mesh( MeshType::STK, tMeshDataInput, tInterpMesh );
        }
        else
        {
            tMeshData = moris::mtk::create_integration_mesh( MeshType::STK, tMeshDataInput);
        }

        if(aOutputOptions.mAddParallelFields && par_size()>1)
        {
            tVizTool.populate_parallel_cell_fields_on_mesh(tMeshData);
        }

        if(moris::par_rank() == 0 && mVerbose)
        {
            std::cout<<"XTK: Write to MTK completed in " <<(std::clock() - start) / (double)(CLOCKS_PER_SEC)<<" s."<<std::endl;
            tMeshDataInput.print_summary();
        }

        return tMeshData;
    }

    //------------------------------------------------------------------------------

    moris::Matrix<moris::IndexMat>
    Model::get_node_map_restricted_to_output_phases(
            Output_Options                 const & aOutputOptions,
            moris::Matrix<moris::IndexMat>       & aOutputtedNodeInds)
    {
        moris::Matrix<moris::IndexMat> tNodeMap = mBackgroundMesh.get_local_to_global_map(EntityRank::NODE);

        moris_index tMyProcRank = par_rank();

        aOutputtedNodeInds.resize(tNodeMap.n_rows(),tNodeMap.n_cols());
        moris::uint tNumNodes = tNodeMap.numel();
        // if we are returning all phases there is no reason to restrict the map
        if(aOutputOptions.output_all_phases())
        {
            for(moris::uint i = 0; i <tNumNodes; i++)
            {
                aOutputtedNodeInds(i) = mBackgroundMesh.get_loc_entity_ind_from_entity_glb_id(tNodeMap(i),EntityRank::NODE);
            }

            return tNodeMap;
        }
        else
        {
            moris::Matrix<moris::IndexMat> tRestrictedNodeMap(tNodeMap.n_rows(),tNodeMap.n_cols());

            moris::uint tCount = 0;
            for(moris::uint i = 0; i <tNumNodes; i++)
            {
                if(output_node(i,aOutputOptions) && mBackgroundMesh.get_vertex_owner(i) == tMyProcRank)
                {
                    moris_index   tVertexIndex = mBackgroundMesh.get_loc_entity_ind_from_entity_glb_id(tNodeMap(i),EntityRank::NODE);
                    aOutputtedNodeInds(tCount) = tVertexIndex;
                    tRestrictedNodeMap(tCount) = tNodeMap(i);

                    tCount++;
                }
            }

            tRestrictedNodeMap.resize(1,tCount);
            aOutputtedNodeInds.resize(1,tCount);

            return tRestrictedNodeMap;
        }
    }

    //------------------------------------------------------------------------------

    void
    Model::setup_interface_single_side_sets(
            Output_Options              const & aOutputOptions,
            Cell<moris::Matrix<moris::IdMat>> & aCellIdsAndSideOrds,
            Cell<std::string>                 & aInterfaceSetNames)
    {
        std::string tSetNameBase = "iside_";

        for(moris_index iG = 0; iG< (moris_index)mGeometryEngine->get_num_geometries(); iG++)
        {
            for(moris_index iP0 = 0; iP0 < (moris_index)mGeometryEngine->get_num_bulk_phase(); iP0++)
            {
                for(moris_index iP1 = iP0+1; iP1 < (moris_index)mGeometryEngine->get_num_bulk_phase(); iP1++)
                {
                    if(aOutputOptions.output_phase(iP0) && aOutputOptions.output_phase(iP1))
                    {
                        std::string tSetName = tSetNameBase+"g_"+ std::to_string(iG) + "_p0_"+std::to_string(iP0)+"_p1_"+std::to_string(iP1);

                        aInterfaceSetNames.push_back(tSetName);
                        aCellIdsAndSideOrds.push_back(mCutMesh.pack_interface_sides(iG,iP0,iP1));
                    }
                }
            }
        }
    }

    //------------------------------------------------------------------------------

    moris::Cell<moris::mtk::MtkNodeSetInfo>
    Model::propogate_background_node_sets(
            moris::Cell<moris::Matrix<IndexMat>>       & aNodeSetData,
            Output_Options                       const & aOutputOptions)
    {
        // access background mesh data
        moris::mtk::Mesh & tMeshData = mBackgroundMesh.get_mesh_data();

        // get all node set names in background mesh
        moris::Cell<std::string> tSetNames = tMeshData.get_set_names(EntityRank::NODE);

        // allocate output
        aNodeSetData = moris::Cell<moris::Matrix<IndexMat>>(tSetNames.size());
        moris::Cell<moris::mtk::MtkNodeSetInfo> tNodeSetInfo(tSetNames.size());
        for(moris::uint i = 0; i <tSetNames.size(); i++)
        {
            moris::uint tCount = 0;
            moris::Cell<moris::mtk::Vertex const *> tNodesInSetInds = tMeshData.get_vertices_in_vertex_set_no_aura(tSetNames(i));

            aNodeSetData(i) = moris::Matrix<moris::IndexMat>(tNodesInSetInds.size(),1);

            for(moris::uint iNode =0; iNode<tNodesInSetInds.size(); iNode++)
            {
                moris_index tVertexInd = tNodesInSetInds(iNode)->get_index();

                if(this->output_node(tVertexInd,aOutputOptions))
                {
                    aNodeSetData(i)(tCount) = tNodesInSetInds(iNode)->get_index();
                    tCount++;
                }
            }

            aNodeSetData(i).resize(tCount,1);

            mBackgroundMesh.convert_loc_entity_ind_to_glb_entity_ids(EntityRank::NODE,aNodeSetData(i));

            tNodeSetInfo(i).mNodeIds = &aNodeSetData(i);
            tNodeSetInfo(i).mNodeSetName = tSetNames(i);
        }

        return tNodeSetInfo;
    }

    //------------------------------------------------------------------------------

    moris::Cell<moris::mtk::MtkSideSetInfo>
    Model::propogate_background_side_sets(moris::Cell<moris::Matrix<IndexMat>> & aSideSetElemIdsAndOrds,
            Output_Options const & aOutputOptions)
    {
        // access background mesh data
        moris::mtk::Mesh & tMeshData = mBackgroundMesh.get_mesh_data();

        // get all side set names in background mesh
        moris::Cell<std::string> tSetNames = tMeshData.get_set_names(tMeshData.get_facet_rank());

        // remove internal side sets which show up with a generated string
        tSetNames = check_for_and_remove_internal_seacas_side_sets(tSetNames);

        // allocate output side sets
        moris::Cell<moris::mtk::MtkSideSetInfo> tSideSets(2*tSetNames.size());
        aSideSetElemIdsAndOrds = moris::Cell<moris::Matrix<IndexMat>>(2*tSetNames.size());

        for(moris::uint i = 0; i <tSetNames.size(); i++)
        {
            moris::uint tNoChildInd = 2*i;
            moris::uint tChildInd = 2*i+1;

            this->propogate_background_side_set(
                    tSetNames(i),
                    tNoChildInd,
                    tChildInd,
                    aSideSetElemIdsAndOrds,
                    tSideSets,
                    aOutputOptions,
                    false);
        }

        return tSideSets;
    }

    // ----------------------------------------------------------------------------------

    void
    Model::propogate_background_side_set( 
            std::string                         const & aSideSetName,
            moris::moris_index                          aNoChildIndex,
            moris::moris_index                          aChildIndex,
            moris::Cell<moris::Matrix<IndexMat>>      & aElementIdsAndSideOrd,
            moris::Cell<moris::mtk::MtkSideSetInfo>   & aSideSetData,
            Output_Options                      const & aOutputOptions,
            bool                                        aOutputIndices)
    {
        // access background mesh data
        moris::mtk::Mesh & tMeshData = mBackgroundMesh.get_mesh_data();

        // estimate maximum number of elements on face
        const uint tMaxElemOnFace = 100;

        // matrices used throughout routine
        moris::Matrix< moris::IndexMat > tElementsAttachedToFace(1,1);
        moris::Matrix< moris::IdMat >    tChildElemsIdsOnFace(1,tMaxElemOnFace);
        moris::Matrix< moris::IndexMat > tChildElemsCMIndOnFace(1,tMaxElemOnFace);
        moris::Matrix< moris::IndexMat > tChildElemOnFaceOrdinal(1,tMaxElemOnFace);

        moris::uint tElementIndex          = 0;
        moris::uint tPhaseIndex            = 0;
        moris::uint tFaceOrdinal           = 0;
        moris::moris_index tChildMeshIndex = 0;
        moris::moris_id    tElementId      = 0;
        moris::moris_index tMyProcRank     = par_rank();
        bool tHasChildren                  = false;

        // get cells and sides in side set
        moris::Cell< mtk::Cell const * > tCellsInSideSet;
        moris::Matrix< moris::IndexMat > tSideSetOrdinals;
        tMeshData.get_sideset_cells_and_ords(aSideSetName, tCellsInSideSet, tSideSetOrdinals );

        // side set data non-intersected
        aElementIdsAndSideOrd(aNoChildIndex)   = Matrix<IndexMat>(tCellsInSideSet.size()*2,2);

        // intersected data
        //TODO: FIGURE OUT MAXIMUM VALUE
        aElementIdsAndSideOrd(aChildIndex) = Matrix<IndexMat>(tCellsInSideSet.size()*2*10,2);

        // keep count
        moris::Cell<moris::uint> tCount(2,0);

        // iterate through sides in set i
        for(moris::uint iSide= 0; iSide<tCellsInSideSet.size(); iSide++)
        {
            tElementIndex = tCellsInSideSet(iSide)->get_index();

            // sides attached to cell
            moris::Matrix<moris::IdMat> tElementFaces = tMeshData.get_entity_connected_to_entity_loc_inds(tElementIndex,EntityRank::ELEMENT,EntityRank::FACE);

            moris_index tSideIndex = tElementFaces(tSideSetOrdinals(iSide));

            if(tMeshData.get_entity_owner(tElementIndex, EntityRank::ELEMENT) == (uint)tMyProcRank)
            {
                tHasChildren = mBackgroundMesh.entity_has_children(tElementIndex,EntityRank::ELEMENT);
                // get the faces from the child mesh
                if(tHasChildren)
                {
                    tChildMeshIndex = mBackgroundMesh.child_mesh_index(tElementIndex,EntityRank::ELEMENT);

                    Child_Mesh const & tChildMesh = mCutMesh.get_child_mesh(tChildMeshIndex);

                    // define variable for actual number of child elements on face
                    uint tNumberOfChildElemsOnFace;

                    tChildMesh.get_child_elements_connected_to_parent_facet(
                            tSideIndex,
                            tNumberOfChildElemsOnFace,
                            tChildElemsIdsOnFace,
                            tChildElemsCMIndOnFace,
                            tChildElemOnFaceOrdinal);

                    moris::Matrix< moris::IndexMat > const & tChildElementPhaseIndices = tChildMesh.get_element_phase_indices();
                    moris::Matrix< moris::IndexMat > const & tChildElementIndices = tChildMesh.get_element_inds();
                    moris::Matrix< moris::IndexMat > const & tElementIds = tChildMesh.get_element_ids();

                    for(moris::moris_index iCElem  = 0; iCElem < (moris::moris_index)tNumberOfChildElemsOnFace; iCElem++)
                    {
                        tPhaseIndex = tChildElementPhaseIndices(0,tChildElemsCMIndOnFace(0,iCElem));

                        if(aOutputOptions.output_phase(tPhaseIndex))
                        {
                            // Child Element Id
                            if(!aOutputIndices)
                            {
                                tElementId   = tElementIds(tChildElemsCMIndOnFace(iCElem));
                                tFaceOrdinal = tChildElemOnFaceOrdinal(iCElem);

                                aElementIdsAndSideOrd(aChildIndex)(tCount(1),0) = tElementId;
                                aElementIdsAndSideOrd(aChildIndex)(tCount(1),1) = tFaceOrdinal;

                                tCount(1)++;
                            }
                            else
                            {
                                tFaceOrdinal   = tChildElemOnFaceOrdinal(iCElem);

                                aElementIdsAndSideOrd(aChildIndex)(tCount(1),0) = tChildElementIndices(tChildElemsCMIndOnFace(iCElem));
                                aElementIdsAndSideOrd(aChildIndex)(tCount(1),1) = tFaceOrdinal;

                                tCount(1)++;
                            }
                        }
                    }
                }
                else
                {
                    tFaceOrdinal = tSideSetOrdinals(iSide);
                    tElementId   = tCellsInSideSet(iSide)->get_id();

                    if(aOutputOptions.output_phase(mBackgroundMesh.get_element_phase_index(tElementIndex)))
                    {
                        if(!aOutputIndices)
                        {
                            aElementIdsAndSideOrd(aNoChildIndex)(tCount(0),0) = tElementId;
                            aElementIdsAndSideOrd(aNoChildIndex)(tCount(0),1) = tFaceOrdinal;
                            tCount(0)++;
                        }
                        else
                        {
                            aElementIdsAndSideOrd(aNoChildIndex)(tCount(0),0) = tElementIndex;
                            aElementIdsAndSideOrd(aNoChildIndex)(tCount(0),1) = tFaceOrdinal;
                            tCount(0)++;
                        }
                    }
                }
            }
        }

        // resize data
        aElementIdsAndSideOrd(aChildIndex).resize(tCount(1),2);
        aElementIdsAndSideOrd(aNoChildIndex).resize(tCount(0),2);

        // Add data to side set info
        // no child
        aSideSetData(aNoChildIndex).mElemIdsAndSideOrds = &aElementIdsAndSideOrd(aNoChildIndex);
        aSideSetData(aNoChildIndex).mSideSetName        = aSideSetName;
        aSideSetData(aNoChildIndex).mSideTopology       = CellTopology::QUAD4;
        aSideSetData(aChildIndex).mElemIdsAndSideOrds   = &aElementIdsAndSideOrd(aChildIndex);
        aSideSetData(aChildIndex).mSideSetName          = aSideSetName + "_i";
        aSideSetData(aChildIndex).mSideTopology         = CellTopology::TRI3;
    }

    // ----------------------------------------------------------------------------------

    moris::Cell<std::string>
    Model::check_for_and_remove_internal_seacas_side_sets(moris::Cell<std::string> & aSideSetNames)
    {
        for(std::vector<std::string>::iterator iSet = aSideSetNames.begin(); iSet != aSideSetNames.end(); ++iSet)
        {
            if(iSet->compare("surface_1_quad4") == 0)
            {
                aSideSetNames.data().erase(iSet--);
            }

            else if(iSet->compare("surface_2_quad4") == 0)
            {
                aSideSetNames.data().erase(iSet--);
            }
            else if(iSet->compare("surface_3_quad4") == 0)
            {
                aSideSetNames.data().erase(iSet--);
            }
            else if(iSet->compare("surface_4_quad4") == 0)
            {
                aSideSetNames.data().erase(iSet--);
            }
            else if(iSet->compare("surface_5_quad4") == 0)
            {
                aSideSetNames.data().erase(iSet--);
            }
            else if(iSet->compare("surface_6_quad4") == 0)
            {
                aSideSetNames.data().erase(iSet--);
            }
            else if(iSet->compare("surface_hex8_quad_1") == 0)
            {
                aSideSetNames.data().erase(iSet--);
            }
            else if(iSet->compare("surface_hex8_quad_2") == 0)
            {
                aSideSetNames.data().erase(iSet--);
            }
            else if(iSet->compare("surface_hex8_quad4_1") == 0)
            {
                aSideSetNames.data().erase(iSet--);
            }
            else if(iSet->compare("surface_hex8_quad4_2") == 0)
            {
                aSideSetNames.data().erase(iSet--);
            }
            else if(iSet->compare("surface_hex8_quad4_3") == 0)
            {
                aSideSetNames.data().erase(iSet--);
            }
            else if(iSet->compare("surface_hex8_quad4_4") == 0)
            {
                aSideSetNames.data().erase(iSet--);
            }
        }

        return aSideSetNames;
    }

    //------------------------------------------------------------------------------

    Cell<moris::Matrix<moris::IdMat>>
    Model::combine_interface_and_non_interface_blocks(
            Cell<moris::Matrix<moris::IdMat>> & aChildElementsByPhase,
            Cell<moris::Matrix<moris::IdMat>> & aNoChildElementsByPhase)
            {
        moris::uint tNumPhase = aChildElementsByPhase.size();

        Cell<moris::Matrix<moris::IdMat>> tCombinedElementsByPhase(tNumPhase);

        for(moris::uint i =0; i<tNumPhase; i++)
        {
            moris::uint tNumChildElems   = aChildElementsByPhase(i).numel();
            moris::uint tNumNoChildElems = aNoChildElementsByPhase(i).numel();

            tCombinedElementsByPhase(i) = moris::Matrix<moris::IdMat>(1,tNumChildElems + tNumNoChildElems);

            tCombinedElementsByPhase(i)({0,0},{0,tNumChildElems-1}) = aChildElementsByPhase(i).get_row(0);
            tCombinedElementsByPhase(i)({0,0},{tNumChildElems,tNumChildElems + tNumNoChildElems -1}) = aNoChildElementsByPhase(i).get_row(0);
        }

        return tCombinedElementsByPhase;
            }

    //------------------------------------------------------------------------------

    uint
    Model::get_num_phases_to_output(Output_Options const & aOutputOptions)
    {
        uint tNumPhasesOutput = 0;
        if(aOutputOptions.output_all_phases())
        {
            tNumPhasesOutput = mGeometryEngine->get_num_bulk_phase();
        }
        else
        {
            tNumPhasesOutput = aOutputOptions.num_phases_to_output();
        }

        return tNumPhasesOutput;
    }

    //------------------------------------------------------------------------------

    void
    Model::setup_cell_clusters_for_output(
            moris::mtk::Cell_Cluster_Input       & aCellClusterInput,
            Output_Options                 const & aOutputOptions,
            moris::Cell<Matrix<IdMat>>           & aCellIds)
    {
        // iterate through child meshes and construct cells
        uint tNumChildMeshes = mCutMesh.get_num_child_meshes();

        for(moris::uint i = 0; i < tNumChildMeshes; i ++)
        {
            // Get child mesh
            Child_Mesh const & tChildMesh = mCutMesh.get_child_mesh(i);

            // pack the element ids into phase grouping
            Cell<moris::Matrix< moris::IdMat >> tElementIds;
            Cell<moris::Matrix< moris::IdMat >> tCMElementInds;
            tChildMesh.pack_child_mesh_by_phase(mGeometryEngine->get_num_bulk_phase(), tElementIds, tCMElementInds);

            // add them to cell to keep in scope
            aCellIds.push_back(tElementIds(0));
            aCellIds.push_back(tElementIds(1));
        }

        for(moris::uint i = 0; i < tNumChildMeshes; i ++)
        {
            // Get child mesh
            Child_Mesh const & tChildMesh = mCutMesh.get_child_mesh(i);

            // primary index
            moris_index tPrimaryCellIndex = 2*i;
            moris_index tVoidCellIndex    = 2*i+1;

            // parent index
            moris_index tParentCellIndex = tChildMesh.get_parent_element_index();

            // access the parent element from the background mesh
            moris::mtk::Cell* tInterpCell = &mBackgroundMesh.get_mesh_data().get_mtk_cell(tParentCellIndex);

            // add to cluster
            aCellClusterInput.add_cluster_data(
                    tInterpCell,
                    &aCellIds(tPrimaryCellIndex),
                    &aCellIds(tVoidCellIndex),
                    &tChildMesh.get_node_ids(),
                    &tChildMesh.get_parametric_coordinates());
        }
    }

    //------------------------------------------------------------------------------

    void
    Model::setup_interface_side_cluster(
            std::string                      aInterfaceSideLabelBase,
            moris::mtk::Side_Cluster_Input & aSideClusterInput,
            Output_Options           const & aOutputOptions,
            moris::Cell<Matrix<IdMat>>     & aCellIdsandSideOrds,
            moris::Cell<Matrix<DDRMat>>    & aParametricCoordinates)
    {
        moris::uint tNumPhases = mGeometryEngine->get_num_bulk_phase();

        moris::uint tNumChildMeshes = mCutMesh.get_num_child_meshes();

        for(moris::uint  iP = 0; iP<tNumPhases; iP++)
        {
            // if we are outputting this phase
            //        if(aOutputOptions.output_phase((size_t)iP))
            if(iP == 0)
            {
                // add side set to output
                std::string tSetName = aInterfaceSideLabelBase;

                //iterate through children meshes
                for(moris::uint iC = 0; iC < tNumChildMeshes; iC ++)
                {
                    // Get child mesh
                    Child_Mesh const & tChildMesh = mCutMesh.get_child_mesh(iC);

                    // package this child element by bulk phase
                    moris::Matrix< moris::IdMat > tInterfaceElementIdsAndSideOrd = tChildMesh.pack_interface_sides( 0, 0, 1 );

                    // add to data which will stay in scope
                    aCellIdsandSideOrds.push_back(tInterfaceElementIdsAndSideOrd);
                }
            }
        }

        uint tCount = 0;
        for(moris::uint  iP = 0; iP<tNumPhases; iP++)
        {
            // if we are outputting this phase
            //        if(aOutputOptions.output_phase((size_t)iP))
            if(iP == 0)
            {
                // add side set to output
                std::string tSetName = aInterfaceSideLabelBase;
                moris_index tSideSetOrd = aSideClusterInput.add_side_set_label(tSetName);

                //iterate through children meshes
                for(moris::uint iC = 0; iC < tNumChildMeshes; iC ++)
                {
                    // Get child mesh
                    Child_Mesh const & tChildMesh = mCutMesh.get_child_mesh(iC);

                    // parent cell index
                    moris_index tParentCellIndex = tChildMesh.get_parent_element_index();

                    // access the parent element from the background mesh
                    moris::mtk::Cell* tInterpCell = &mBackgroundMesh.get_mesh_data().get_mtk_cell(tParentCellIndex);

                    // add to cluster input data
                    //fixme: Add only vertex indices on the interface to cluster. Adding all.
                    aSideClusterInput.add_cluster_data(false,tSideSetOrd,tInterpCell,&aCellIdsandSideOrds(tCount),&tChildMesh.get_node_ids(),&tChildMesh.get_parametric_coordinates());

                    tCount++;
                }
            }
        }
    }

    //------------------------------------------------------------------------------

    bool
    Model::output_node(
            moris::moris_index     aNodeIndex,
            Output_Options const & aOutputOptions)
    {
        bool tIsInterface = mBackgroundMesh.is_interface_node(aNodeIndex,0);
        moris::size_t tPhaseIndex = 0;
        mGeometryEngine->get_phase_index(aNodeIndex,
                mBackgroundMesh.get_selected_node_coordinates_loc_inds({{aNodeIndex}}));

        if(aOutputOptions.output_phase(tPhaseIndex) && !tIsInterface)
        {
            return true;
        }
        else if(tIsInterface)
        {
            return true;
        }

        return false;
    }

    //------------------------------------------------------------------------------

    moris::size_t
    Model::determine_element_phase_index(
            moris::size_t                            aRowIndex,
            moris::Matrix< moris::IndexMat > const & aElementToNodeIndex)
    {
        moris::size_t tNumGeom = mGeometryEngine->get_num_geometries();
        moris::size_t tNumNodesPerElem = aElementToNodeIndex.n_cols();
        moris::Matrix< moris::IndexMat > tNodalPhaseVals(1,tNumGeom,MORIS_INDEX_MAX);

        // allocate phase on or off value (either 0 or 1)
        Matrix<IndexMat> tPhaseVotes(1,2);
        tPhaseVotes.fill(0);

        // maximum row index and column index
        uint tMaxRow = 0;
        uint tMaxCol = 0;
        for (moris::uint i = 0; i < tNumGeom; i++)
        {
            bool tFoundNonInterfaceNode = false;

            for( moris::size_t j = 0; j<tNumNodesPerElem; j++)
            {
                if(!mBackgroundMesh.is_interface_node(aElementToNodeIndex(aRowIndex,j),i))
                {
                    moris::uint tNodeIndex    = (moris::uint)aElementToNodeIndex(aRowIndex, j);
                    moris_index tPhaseIndex   = mGeometryEngine->get_node_phase_index_wrt_a_geometry(tNodeIndex,i);
                    tFoundNonInterfaceNode    = true;
                    tPhaseVotes(tPhaseIndex)++;
                }
            }

            // take the phase with the maximum number of votes
            tPhaseVotes.max(tMaxRow,tMaxCol);
            tNodalPhaseVals(0,i) = tMaxCol;

            tPhaseVotes.fill(0);

            MORIS_ERROR(tFoundNonInterfaceNode,"Did not find a non-interface node for this element");
        }

        moris::moris_index tElemPhaseVal = mGeometryEngine->get_elem_phase_index(tNodalPhaseVals);

        return tElemPhaseVal;
    }

    // ----------------------------------------------------------------------------------

    void
    Model::print_decompsition_preamble(Cell<enum Subdivision_Method> aMethods)
    {
        // Only process with rank 0 prints the preamble

        if(moris::par_rank() == 0 && mVerbose)
        {
            std::cout<<"XTK: Specified Decomposition Routines: ";

            for(moris::size_t i = 0 ; i<aMethods.size(); i++)
            {
                std::cout<<"["<<get_enum_str(aMethods(i))<<  "] ";
            }

            std::cout<<std::endl;
        }
    }

    //------------------------------------------------------------------------------

    moris::Cell<std::string>
    Model::assign_geometry_data_names()
    {
        uint tNumGeometries = mGeometryEngine->get_num_geometries();

        // base string of geometry data
        std::string tBaseName = "gd_";

        // Allocate output
        moris::Cell<std::string> tGeometryFieldName(tNumGeometries);

        //Iterate through geometries
        for(uint iG = 0; iG <tNumGeometries; iG++)
        {
            tGeometryFieldName(iG) = tBaseName+std::to_string(iG);
        }

        return tGeometryFieldName;
    }

    //------------------------------------------------------------------------------

    moris::Cell < enum moris::EntityRank >
    Model::assign_geometry_data_field_ranks()
    {
        uint tNumGeometries = mGeometryEngine->get_num_geometries();

        // base string of geometry data
        std::string tBaseName = "gd_";

        // Allocate output
        // Note: for now this is a nodal field always
        moris::Cell<enum moris::EntityRank> tGeometryFieldRank(tNumGeometries,moris::EntityRank::NODE);

        return tGeometryFieldRank;
    }

    //------------------------------------------------------------------------------

    Matrix<DDRMat>
    Model::get_timing_data() const
    {
        Matrix<DDRMat> tTimingMat(mTimingData.size(),1);
        for(moris::uint i = 0; i < mTimingData.size(); i++)
        {
            tTimingMat(i) = mTimingData(i);
        }

        return tTimingMat;
    }

    //------------------------------------------------------------------------------

    Cell<std::string>
    Model::get_timing_labels() const
    {
        return mTimingLabels;
    }

    //------------------------------------------------------------------------------

    void
    Model::print_timing_data() const
    {
        if(par_rank() == 0)
        {
            std::cout<<"XTK Timing Data:"<<std::endl;
            for(moris::uint i = 0; i < mTimingData.size(); i++)
            {
                std::cout<<std::setw(36)<<mTimingLabels(i)<<": "<<std::setw(8)<< std::scientific<<mTimingData(i)<<std::endl;
            }
        }
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

        if(mDecomposed)
        {
            tCutMeshMM    = mCutMesh.get_memory_usage();
            tBGMeshMM     = mBackgroundMesh.get_memory_usage();
        }

        if(mEnriched)
        {
            tEnrichmentMM = mEnrichment->get_memory_usage();
            tIgMeshMM     = this->get_enriched_integ_mesh().get_memory_usage();
            tIpMeshMM     = this->get_enriched_interp_mesh().get_memory_usage();
        }

        if(mGhost)
        {
            tGhostMM    = mGhostStabilization->get_memory_usage();
        }

        // make the sum of the cut mesh memory map the cut mesh memory
        tXTKModelMM.mMemoryMapData["Cut Mesh"] = tCutMeshMM.sum();
        tXTKModelMM.mMemoryMapData["Enrichment"] = tEnrichmentMM.sum();
        tXTKModelMM.mMemoryMapData["Enriched Ig Mesh"] = tIgMeshMM.sum();
        tXTKModelMM.mMemoryMapData["Enriched Ip Mesh"] = tIpMeshMM.sum();
        tXTKModelMM.mMemoryMapData["Ghost"] = tGhostMM.sum();
        tXTKModelMM.mMemoryMapData["Background Mesh"] = tBGMeshMM.sum();
        tXTKModelMM.mMemoryMapData["mElementToElement ptrs"] = moris::internal_capacity(mElementToElement);
        tXTKModelMM.mMemoryMapData["mElementToElement ptrs"] = moris::internal_capacity(mElementToElement);
        tXTKModelMM.mMemoryMapData["mSubphaseToSubPhase"] = moris::internal_capacity(mSubphaseToSubPhase);
        tXTKModelMM.mMemoryMapData["mSubphaseToSubPhaseMySideOrds"] = moris::internal_capacity(mSubphaseToSubPhaseMySideOrds);
        tXTKModelMM.mMemoryMapData["mSubphaseToSubPhaseNeighborSideOrds"] = moris::internal_capacity(mSubphaseToSubPhaseNeighborSideOrds);

        tIgMeshMM.par_print("Ig Mesh");
        return tXTKModelMM;
    }

    //------------------------------------------------------------------------------

    void
    Model::save_timing_to_hdf5( const std::string & aFilePath ) const
    {
        if(par_rank() == 0)
        {
            // timing data
            Matrix<DDRMat>    tTimingData   = this->get_timing_data();
            Cell<std::string> tTimingLabels = this->get_timing_labels();

            // Create a new file using default properties
            hid_t tFileID = H5Fcreate(
                    aFilePath.c_str(),
                    H5F_ACC_TRUNC,
                    H5P_DEFAULT,
                    H5P_DEFAULT);

            // error handler
            herr_t tStatus;

            // save some generic data about this run
            // number of procs
            save_scalar_to_hdf5_file(
                    tFileID,
                    "par_size",
                    par_size(),
                    tStatus );

            // iterate through timing labels and save to the file
            for(moris::uint iL = 0; iL < tTimingLabels.size(); iL++)
            {
                std::string tLabel = "Label_" + std::to_string(iL);

                // save label of this timing data
                save_string_to_hdf5_file(
                        tFileID,
                        tLabel,
                        tTimingLabels(iL),
                        tStatus );
            }

            // save matrix of timing data to file
            save_matrix_to_hdf5_file(
                    tFileID,
                    "Timing Data",
                    tTimingData,
                    tStatus );

            // Close file
            close_hdf5_file(tFileID);
        }
    }

    //------------------------------------------------------------------------------

    void
    Model::save_model_statistics_to_file( const std::string & aFilePath )
    {
        MORIS_ERROR(mDecomposed,"save_model_statistics_to_file() requires model to be at least decomposed");

        // collect global data

        // save the number of intersected background elements
        moris::uint tNumOwnedChildMeshes = this->get_cut_mesh().get_owned_child_meshes().size();
        moris::uint tGlobalChildMeshes = sum_all(tNumOwnedChildMeshes);

        if(par_rank() == 0)
        {
            // Create a new file using default properties
            hid_t tFileID = H5Fcreate(
                    aFilePath.c_str(),
                    H5F_ACC_TRUNC,
                    H5P_DEFAULT,
                    H5P_DEFAULT);

            // error handler
            herr_t tStatus;

            // save the number of background cells
            save_scalar_to_hdf5_file(
                    tFileID,
                    "num background cells",
                    mBackgroundMesh.get_num_entities_background(EntityRank::ELEMENT),
                    tStatus );

            save_scalar_to_hdf5_file(
                    tFileID,
                    "num intersect",
                    tGlobalChildMeshes,
                    tStatus );

            // save the number of sub-phases
            save_scalar_to_hdf5_file(
                    tFileID,
                    "num subphase",
                    this->get_cut_mesh().get_num_subphases(),
                    tStatus );

            // save the number of bulk-phases
            save_scalar_to_hdf5_file(
                    tFileID,
                    "num bulkphase",
                    this->get_geom_engine()->get_num_phases(),
                    tStatus );

            // save the number of geometries
            save_scalar_to_hdf5_file(
                    tFileID,
                    "num geometries",
                    this->get_geom_engine()->get_num_geometries(),
                    tStatus );

            if(mEnriched)
            {
                // add enrichment specific information
                // - number of enriched basis functions
                // - ratio between
            }

            if(mGhost)
            {
                // add ghost specific information
                // - number of ghost facets
                // - number of transition facets
            }

            // Close file
            close_hdf5_file(tFileID);
        }
    }

    //------------------------------------------------------------------------------

    void
    Model::add_timing_data(
            real        const & aTime,
            std::string const & aLabel,
            std::string const & aCategory)
    {
        mTimingData.push_back(aTime);
        mTimingLabels.push_back(aLabel);
    }
}
