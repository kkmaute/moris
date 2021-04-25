/*
 * ut_Snapping.cpp
 *
 *  Created on: Feb 2, 2021
 *      Author: doble
 */

#include "catch.hpp"
#include <string>
#include "cl_Cell.hpp"

// XTKL: Geometry  Include
#include "cl_Logger.hpp"
#include "cl_XTK_Model.hpp"
#include "cl_HMR.hpp"
#include "cl_GEN_Geometry_Engine.hpp"


#include "fn_PRM_FEM_Parameters.hpp"
#include "fn_PRM_HMR_Parameters.hpp"
#include "fn_PRM_GEN_Parameters.hpp"
#include "fn_PRM_XTK_Parameters.hpp"
#include "cl_Param_List.hpp"
#include "fn_stringify_matrix.hpp"

namespace xtk
{


TEST_CASE("Mesh Cleanup","[Mesh_Cleanup_2d]")
{
    if(par_size() == 1)
    {
 
        bool tIs3D     = false;

        // HMR parameters
        real tXDomainDim = 2.0;
        real tYDomainDim = 1.0;
        real tZDomainDim = 0.1;

        uint tNumElemsX = tIs3D ? 6 : 80;
        uint tNumElemsY = tIs3D ? 3 : 40;
        uint tNumElemsZ = 1;

        std::string tNumElemsPerDim     = tIs3D ? std::to_string(tNumElemsX)  + "," + std::to_string(tNumElemsY) + "," + std::to_string(tNumElemsZ) :    
                                                std::to_string(tNumElemsX)  + "," + std::to_string(tNumElemsY);
        std::string tDomainDims         = tIs3D ? std::to_string(tXDomainDim) + "," + std::to_string(tYDomainDim) + "," + std::to_string(tZDomainDim) :
                                                std::to_string(tXDomainDim) + "," + std::to_string(tYDomainDim);
        std::string tDomainOffset       = tIs3D ? "0.0, 0.0, 0.0" : "0.0, 0.0";
        std::string tDomainSidesets     = tIs3D ? "1,2,3,4,5,6" : "1,2,3,4";
        std::string tInterpolationOrder = "1";

        int tRefineBuffer      = 1;

        // element size
        moris::real hy = tYDomainDim/ (real)tNumElemsY;

        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // plane 1
        moris::real tPlaneXCenter1 = 0.501*tXDomainDim ;
        moris::real tPlaneYCenter1 = 0.501*tYDomainDim;
        moris::real tPlaneXNormal1 = 0.0;
        moris::real tPlaneYNormal1 = 1.0;

        // plane 2
        moris::real tPlaneXCenter2 = 0.25*tXDomainDim + 0.01;
        moris::real tPlaneYCenter2 = 0.50*tYDomainDim + 0.00;
        moris::real tPlaneXNormal2 = 1.0;
        moris::real tPlaneYNormal2 = 0.0;

        // plane 3
        moris::real tPlaneXCenter3 = 0.75*tXDomainDim + 0.01;
        moris::real tPlaneYCenter3 = 0.50*tYDomainDim + 0.00;
        moris::real tPlaneXNormal3 = 1.0;
        moris::real tPlaneYNormal3 = 0.0;

        // plane 4
        moris::real tPlaneXCenter4 = 0.50*tXDomainDim + 0.5*hy;
        moris::real tPlaneYCenter4 = 0.50*tYDomainDim + 0.5*hy;
        moris::real tPlaneXNormal4 = 1.0;
        moris::real tPlaneYNormal4 = 0.0;
        
        // XTK parameter list
        moris::Cell< moris::Cell< moris::ParameterList > > tXTKParams(1);
        tXTKParams( 0 ).resize( 1 );
        tXTKParams( 0 )( 0 ) = prm::create_xtk_parameter_list();
        tXTKParams( 0 )( 0 ).set( "decompose",                 true );
        tXTKParams( 0 )( 0 ).set( "decomposition_type",        "conformal" );
        tXTKParams( 0 )( 0 ).set( "enrich",                    true );
        tXTKParams( 0 )( 0 ).set( "basis_rank",                "bspline" );
        tXTKParams( 0 )( 0 ).set( "enrich_mesh_indices",       "0" );
        tXTKParams( 0 )( 0 ).set( "ghost_stab",                true );
        tXTKParams( 0 )( 0 ).set( "multigrid",                 false );
        tXTKParams( 0 )( 0 ).set( "verbose",                   false );
        tXTKParams( 0 )( 0 ).set( "print_enriched_ig_mesh",    false );
        tXTKParams( 0 )( 0 ).set( "exodus_output_XTK_ig_mesh", true );
        tXTKParams( 0 )( 0 ).set( "visualize_ghost", true );
        tXTKParams( 0 )( 0 ).set( "cleanup_cut_mesh", true );
        

        moris::Cell< moris::Cell< moris::ParameterList > > tGenParamList;
        tGenParamList.resize( 3 );
        tGenParamList( 0 ).resize( 1 );
        tGenParamList( 1 ).resize( 1 );


        tGenParamList( 0 )( 0 ) = prm::create_gen_parameter_list();
        tGenParamList( 0 )( 0 ).set("isocontour_threshold", 1e-12);      // Level-set isocontour level
        tGenParamList( 0 )( 0 ).set("isocontour_tolerance", 1e-12);      // Level-set isocontour level
        
        Matrix<DDUMat> tPhaseMap(16,1,7);
        tPhaseMap(8)  = 0;
        tPhaseMap(12) = 0;
        tPhaseMap(13) = 1;
        tPhaseMap(15) = 2;
        tPhaseMap(0)  = 3;
        tPhaseMap(4)  = 4;
        tPhaseMap(5)  = 5;
        tPhaseMap(7)  = 5;
        tGenParamList( 0 )( 0 ).set("phase_table",moris::ios::stringify(tPhaseMap));\

        tGenParamList( 0 )( 0 ).set("initial_advs"      , "0.4999");
        tGenParamList( 0 )( 0 ).set("lower_bounds"      , "0.0001");
        tGenParamList( 0 )( 0 ).set("upper_bounds"      , "0.4999");
        tGenParamList( 0 )( 0 ).set("IQI_types"         , "IQIBulkStrainEnergy");
        
        // Geometry parameter lists
        moris::uint tGeoCounter = 0;


        tGenParamList( 1 )( tGeoCounter ) = prm::create_geometry_parameter_list();
        tGenParamList( 1 )( tGeoCounter ).set( "type", "plane");
        std::string tPlaneStr1 = std::to_string(tPlaneXCenter1)  + "," + std::to_string(tPlaneYCenter1)  + "," + std::to_string(tPlaneXNormal1)  + "," + std::to_string(tPlaneYNormal1); 
        tGenParamList( 1 )( tGeoCounter ).set( "constant_parameters", tPlaneStr1);
        tGeoCounter++;
        
        tGenParamList( 1 ).push_back( prm::create_geometry_parameter_list() ); 
        tGenParamList( 1 )( tGeoCounter ).set( "type", "plane");
        tGenParamList( 1 )( tGeoCounter ).set( "constant_parameters", std::to_string(tPlaneXCenter2) + "," + std::to_string(tPlaneYCenter2)  + "," + std::to_string(tPlaneXNormal2) + "," + std::to_string(tPlaneYNormal2));
        tGeoCounter++;

        tGenParamList( 1 ).push_back( prm::create_geometry_parameter_list() ); 
        tGenParamList( 1 )( tGeoCounter ).set( "type", "plane");
        tGenParamList( 1 )( tGeoCounter ).set( "constant_parameters", std::to_string(tPlaneXCenter3) + "," + std::to_string(tPlaneYCenter3)  + "," + std::to_string(tPlaneXNormal3) + "," + std::to_string(tPlaneYNormal3));
        tGeoCounter++;

        tGenParamList( 1 ).push_back( prm::create_geometry_parameter_list() ); 
        tGenParamList( 1 )( tGeoCounter ).set( "type", "plane");
        tGenParamList( 1 )( tGeoCounter ).set( "constant_parameters", std::to_string(tPlaneXCenter4) + "," + std::to_string(tPlaneYCenter4)  + "," + std::to_string(tPlaneXNormal4) + "," + std::to_string(tPlaneYNormal4));
        tGeoCounter++;

        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // HMR parameters

        moris::Cell< moris::Cell< moris::ParameterList > > tHMRParams(1);
        tHMRParams.resize( 1 );
        tHMRParams( 0 ).resize( 1 );

        tHMRParams( 0 )( 0 ) = prm::create_hmr_parameter_list();
        tHMRParams( 0 )( 0 ).set( "number_of_elements_per_dimension", tNumElemsPerDim );
        tHMRParams( 0 )( 0 ).set( "domain_dimensions",                tDomainDims );
        tHMRParams( 0 )( 0 ).set( "domain_offset",                    tDomainOffset );
        tHMRParams( 0 )( 0 ).set( "domain_sidesets",                  tDomainSidesets);
        tHMRParams( 0 )( 0 ).set( "lagrange_output_meshes",           "0");
        tHMRParams( 0 )( 0 ).set( "lagrange_orders",  tInterpolationOrder );
        tHMRParams( 0 )( 0 ).set( "lagrange_pattern", std::string( "0" )  );
        tHMRParams( 0 )( 0 ).set( "bspline_orders",   tInterpolationOrder );
        tHMRParams( 0 )( 0 ).set( "bspline_pattern",  std::string( "0" )  );
        tHMRParams( 0 )( 0 ).set( "lagrange_to_bspline", "0" );
        tHMRParams( 0 )( 0 ).set( "truncate_bsplines",  1 );
        tHMRParams( 0 )( 0 ).set( "refinement_buffer",  tRefineBuffer );
        tHMRParams( 0 )( 0 ).set( "staircase_buffer",   tRefineBuffer );
        tHMRParams( 0 )( 0 ).set( "initial_refinement", "0" );
        tHMRParams( 0 )( 0 ).set( "initial_refinement_pattern", "0" );
        tHMRParams( 0 )( 0 ).set( "use_number_aura", 1);
        tHMRParams( 0 )( 0 ).set( "use_multigrid",  0 );
        tHMRParams( 0 )( 0 ).set( "severity_level", 0 );

        // PSEUDO Workflow
        // HMR initialize
        std::shared_ptr< hmr::HMR > pHMR = std::make_shared< hmr::HMR >( tHMRParams( 0 )( 0 ) );

        // Geometry engine initialize
        std::shared_ptr< ge::Geometry_Engine > pGEN = std::make_shared< ge::Geometry_Engine >( tGenParamList, nullptr );

        // Initialize  Background Mesh Mesh
        std::shared_ptr< mtk::Mesh_Manager > pBGMTK = std::make_shared< mtk::Mesh_Manager >();

        // Initialize output mesh from XTK
        std::shared_ptr< mtk::Mesh_Manager > pOutputMTK = std::make_shared< mtk::Mesh_Manager >();

        // XTK initialize
        std::shared_ptr< xtk::Model > pXTK = std::make_shared< xtk::Model >( tXTKParams( 0 )( 0 ) );
        
        // Set performer to HMR
        pHMR->set_performer( pBGMTK );
        
        // Set XTK Cooperations
        pXTK->set_geometry_engine( pGEN.get() );
        pXTK->set_input_performer( pBGMTK );
        pXTK->set_output_performer( pOutputMTK );

        // uniform initial refinement
        pHMR->perform_initial_refinement();

        // HMR finalize
        pHMR->perform();
        //
        pGEN->distribute_advs( pBGMTK->get_mesh_pair( 0 ),{} );

        // Output GEN fields, if requested
        pGEN->output_fields(pBGMTK->get_interpolation_mesh( 0 ));

        // XTK perform - decompose - enrich - ghost - multigrid
        pXTK->perform();


        // verify after mesh cleanup that the maps are consistent
        xtk::Background_Mesh* tBGMesh = & pXTK->get_background_mesh();

        for(moris::uint i = 0; i < tBGMesh->get_num_entities(EntityRank::ELEMENT); i++)
        {
            mtk::Cell & tCell = tBGMesh->get_mtk_cell((moris_index)i);
            CHECK(tCell.get_index() == (moris_index) i);

            if(tBGMesh->entity_has_children((moris_index)i, EntityRank::ELEMENT))
            {
                moris_index tCMIndex = pXTK->get_background_mesh().child_mesh_index((moris_index)i,EntityRank::ELEMENT);
                Child_Mesh & tCM = pXTK->get_cut_mesh().get_child_mesh(tCMIndex);
                CHECK(tCM.get_parent_element_index() == (moris_index) i);
            }
        }

        

    }
}

TEST_CASE("Plato Example","[Plato]")
{
    if(par_size() == 1)
    {
 
 // Get default XTK parameters
      moris::ParameterList tParameterlist = moris::prm::create_xtk_parameter_list();

        // HMR parameters
        moris::real tXDomainDim = 10.0;
        moris::real tYDomainDim = 10.0;
        moris::real tZDomainDim = 10.0;

        moris::uint tNumElemsX = 10;
        moris::uint tNumElemsY = 10;
        moris::uint tNumElemsZ = 10;

        std::string tNumElemsPerDim     = std::to_string(tNumElemsX)  + "," + std::to_string(tNumElemsY) + "," + std::to_string(tNumElemsZ);
        std::string tDomainDims         = std::to_string(tXDomainDim) + "," + std::to_string(tYDomainDim) + "," + std::to_string(tZDomainDim);
        std::string tDomainOffset       = "0.0, 0.0, 0.0";
        std::string tDomainSidesets     = "1,2,3,4,5,6";
        std::string tInterpolationOrder = "1";

        int tRefineBuffer      = 1;

        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // circle center
        moris::real tCenterCircleX1 = 0.500*tXDomainDim;
        moris::real tCenterCircleX2 = 0.500*tYDomainDim;
        moris::real tCenterCircleX3 = 0.500*tZDomainDim;
        moris::real tCenterCircleR  = 0.15011*tXDomainDim;

        // right side nbc 2
        moris::real tNBCCircleX1 = 1.000*tXDomainDim;
        moris::real tNBCCircleX2 = 0.500*tYDomainDim;
        moris::real tNBCCircleX3 = 0.500*tZDomainDim;
        moris::real tNBCCircleR  = 0.0511*tXDomainDim;

        // dbc circle 1
        moris::real tDBCCircle1X1 = 0.000*tXDomainDim;
        moris::real tDBCCircle1X2 = 0.000*tYDomainDim;
        moris::real tDBCCircle1X3 = 0.000*tZDomainDim;
        moris::real tDBCCircle1R  = 0.111*tXDomainDim;

        // dbc circle 2
        moris::real tDBCCircle2X1 = 0.000*tXDomainDim;
        moris::real tDBCCircle2X2 = 0.000*tYDomainDim;
        moris::real tDBCCircle2X3 = 1.000*tZDomainDim;
        moris::real tDBCCircle2R  = 0.111*tXDomainDim;

        // dbc circle 3
        moris::real tDBCCircle3X1 = 0.000*tXDomainDim;
        moris::real tDBCCircle3X2 = 1.000*tYDomainDim;
        moris::real tDBCCircle3X3 = 1.000*tZDomainDim;
        moris::real tDBCCircle3R  = 0.111*tXDomainDim;

        // dbc circle 4
        moris::real tDBCCircle4X1 = 0.000*tXDomainDim;
        moris::real tDBCCircle4X2 = 1.000*tYDomainDim;
        moris::real tDBCCircle4X3 = 0.000*tZDomainDim;
        moris::real tDBCCircle4R  = 0.111*tXDomainDim;




        

        
        // XTK parameter list
        moris::Cell< moris::Cell< moris::ParameterList > > tXTKParams(1);
        tXTKParams( 0 ).resize( 1 );
        tXTKParams( 0 )( 0 ) = moris::prm::create_xtk_parameter_list();
        tXTKParams( 0 )( 0 ).set( "decompose",                 true );
        tXTKParams( 0 )( 0 ).set( "decomposition_type",        "conformal" );
        tXTKParams( 0 )( 0 ).set( "enrich",                    true );
        tXTKParams( 0 )( 0 ).set( "basis_rank",                "bspline" );
        tXTKParams( 0 )( 0 ).set( "enrich_mesh_indices",       "0" );
        tXTKParams( 0 )( 0 ).set( "ghost_stab",                false );
        tXTKParams( 0 )( 0 ).set( "multigrid",                 false );
        tXTKParams( 0 )( 0 ).set( "verbose",                   false );
        tXTKParams( 0 )( 0 ).set( "print_enriched_ig_mesh",    false );
        tXTKParams( 0 )( 0 ).set( "exodus_output_XTK_ig_mesh", true );
        tXTKParams( 0 )( 0 ).set( "triangulate_all", true );
        tXTKParams( 0 )( 0 ).set( "output_file", "plato_output.exo" );
        tXTKParams( 0 )( 0 ).set( "union_blocks" ,"HMR_dummy_c_p47,HMR_dummy_c_p55,HMR_dummy_c_p59;HMR_dummy_c_p62,HMR_dummy_c_p61,HMR_dummy_c_p63");
        tXTKParams( 0 )( 0 ).set( "union_block_names" ,"Material;Material2");
        tXTKParams( 0 )( 0 ).set( "union_block_colors" ,"0;1");
        tXTKParams( 0 )( 0 ).set( "union_side_sets" ,"SideSet_4_c_p61,SideSet_4_c_p62;SideSet_4_c_p59,SideSet_4_c_p55");
        tXTKParams( 0 )( 0 ).set( "union_side_set_names" ,"SS1;SS2");
        tXTKParams( 0 )( 0 ).set( "union_side_set_colors" ,"0;1");
        tXTKParams( 0 )( 0 ).set( "deactivate_all_but_blocks","Material,Material2");
        tXTKParams( 0 )( 0 ).set( "deactivate_all_but_side_sets","SS1,SS2");

        moris::Cell< moris::Cell< moris::ParameterList > > tGenParamList;
        tGenParamList.resize( 3 );
        tGenParamList( 0 ).resize( 1 );
        tGenParamList( 1 ).resize( 6 );


        tGenParamList( 0 )( 0 ) = moris::prm::create_gen_parameter_list();
        tGenParamList( 0 )( 0 ).set("isocontour_threshold", 1e-12);      // Level-set isocontour level
        tGenParamList( 0 )( 0 ).set("isocontour_tolerance", 1e-12);      // Level-set isocontour level
        
        // moris::Matrix<moris::DDUMat> tPhaseMap(16,1,7);
        // tPhaseMap(8)  = 0;
        // tPhaseMap(12) = 0;
        // tPhaseMap(13) = 1;
        // tPhaseMap(15) = 2;
        // tPhaseMap(0)  = 3;
        // tPhaseMap(4)  = 4;
        // tPhaseMap(5)  = 5;
        // tPhaseMap(7)  = 5;
        // tGenParamList( 0 )( 0 ).set("phase_table",moris::ios::stringify(tPhaseMap));

        tGenParamList( 0 )( 0 ).set("initial_advs"      , "0.4999");
        tGenParamList( 0 )( 0 ).set("lower_bounds"      , "0.0001");
        tGenParamList( 0 )( 0 ).set("upper_bounds"      , "0.4999");
        
        // Geometry parameter lists
        moris::uint tGeoCounter = 0;
        tGenParamList( 1 )( tGeoCounter ) = moris::prm::create_geometry_parameter_list();
        tGenParamList( 1 )( tGeoCounter ).set( "type", "sphere");
        std::string tGeomStr = std::to_string(tCenterCircleX1)  + "," + std::to_string(tCenterCircleX2)  + "," + std::to_string(tCenterCircleX3)  + "," + std::to_string(tCenterCircleR); 
        tGenParamList( 1 )( tGeoCounter ).set( "constant_parameters", tGeomStr);
        tGeoCounter++;


        // right side nbc 
        tGenParamList( 1 )( tGeoCounter ) = moris::prm::create_geometry_parameter_list();
        tGenParamList( 1 )( tGeoCounter ).set( "type", "sphere");
        tGeomStr = std::to_string(tNBCCircleX1)  + "," + std::to_string(tNBCCircleX2)  + "," + std::to_string(tNBCCircleX3)  + "," + std::to_string(tNBCCircleR); 
        tGenParamList( 1 )( tGeoCounter ).set( "constant_parameters", tGeomStr);
        tGeoCounter++;

        tGenParamList( 1 )( tGeoCounter ) = moris::prm::create_geometry_parameter_list();
        tGenParamList( 1 )( tGeoCounter ).set( "type", "sphere");
        tGeomStr = std::to_string(tDBCCircle1X1)  + "," 
                  + std::to_string(tDBCCircle1X2)  + "," 
                  + std::to_string(tDBCCircle1X3)  + "," 
                  + std::to_string(tDBCCircle1R); 
        tGenParamList( 1 )( tGeoCounter ).set( "constant_parameters", tGeomStr);
        tGeoCounter++;

        tGenParamList( 1 )( tGeoCounter ) = moris::prm::create_geometry_parameter_list();
        tGenParamList( 1 )( tGeoCounter ).set( "type", "sphere");
        tGeomStr = std::to_string(tDBCCircle2X1)  + "," 
                 + std::to_string(tDBCCircle2X2)  + "," 
                 + std::to_string(tDBCCircle2X3)  + "," 
                 + std::to_string(tDBCCircle2R); 
        tGenParamList( 1 )( tGeoCounter ).set( "constant_parameters", tGeomStr);
        tGeoCounter++;

        tGenParamList( 1 )( tGeoCounter ) = moris::prm::create_geometry_parameter_list();
        tGenParamList( 1 )( tGeoCounter ).set( "type", "sphere");
        tGeomStr = std::to_string(tDBCCircle3X1)  + "," 
                             + std::to_string(tDBCCircle3X2)  + "," 
                             + std::to_string(tDBCCircle3X3)  + "," 
                             + std::to_string(tDBCCircle3R); 
        tGenParamList( 1 )( tGeoCounter ).set( "constant_parameters", tGeomStr);
        tGeoCounter++;

        tGenParamList( 1 )( tGeoCounter ) = moris::prm::create_geometry_parameter_list();
        tGenParamList( 1 )( tGeoCounter ).set( "type", "sphere");
        tGeomStr = std::to_string(tDBCCircle4X1)  + "," 
                 + std::to_string(tDBCCircle4X2)  + "," 
                 + std::to_string(tDBCCircle4X3)  + "," 
                 + std::to_string(tDBCCircle4R); 
        tGenParamList( 1 )( tGeoCounter ).set( "constant_parameters", tGeomStr);
        tGeoCounter++;

        
        
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // HMR parameters

        moris::Cell< moris::Cell< moris::ParameterList > > tHMRParams(1);
        tHMRParams.resize( 1 );
        tHMRParams( 0 ).resize( 1 );

        tHMRParams( 0 )( 0 ) = moris::prm::create_hmr_parameter_list();
        tHMRParams( 0 )( 0 ).set( "number_of_elements_per_dimension", tNumElemsPerDim );
        tHMRParams( 0 )( 0 ).set( "domain_dimensions",                tDomainDims );
        tHMRParams( 0 )( 0 ).set( "domain_offset",                    tDomainOffset );
        tHMRParams( 0 )( 0 ).set( "domain_sidesets",                  tDomainSidesets);
        tHMRParams( 0 )( 0 ).set( "lagrange_output_meshes",           "0");
        tHMRParams( 0 )( 0 ).set( "lagrange_orders",  tInterpolationOrder );
        tHMRParams( 0 )( 0 ).set( "lagrange_pattern", std::string( "0" )  );
        tHMRParams( 0 )( 0 ).set( "bspline_orders",   tInterpolationOrder );
        tHMRParams( 0 )( 0 ).set( "bspline_pattern",  std::string( "0" )  );
        tHMRParams( 0 )( 0 ).set( "lagrange_to_bspline", "0" );
        tHMRParams( 0 )( 0 ).set( "truncate_bsplines",  1 );
        tHMRParams( 0 )( 0 ).set( "refinement_buffer",  tRefineBuffer );
        tHMRParams( 0 )( 0 ).set( "staircase_buffer",   tRefineBuffer );
        tHMRParams( 0 )( 0 ).set( "initial_refinement", "0" );
        tHMRParams( 0 )( 0 ).set( "initial_refinement_pattern", "0" );
        tHMRParams( 0 )( 0 ).set( "use_number_aura", 1);
        tHMRParams( 0 )( 0 ).set( "use_multigrid",  0 );
        tHMRParams( 0 )( 0 ).set( "severity_level", 0 );

        // PSEUDO Workflow
        // HMR initialize
        std::shared_ptr< hmr::HMR > pHMR = std::make_shared< hmr::HMR >( tHMRParams( 0 )( 0 ) );

        // Geometry engine initialize
        std::shared_ptr< ge::Geometry_Engine > pGEN = std::make_shared< ge::Geometry_Engine >( tGenParamList, nullptr );

        // Initialize  Background Mesh Mesh
        std::shared_ptr< mtk::Mesh_Manager > pBGMTK = std::make_shared< mtk::Mesh_Manager >();

        // Initialize output mesh from XTK
        std::shared_ptr< mtk::Mesh_Manager > pOutputMTK = std::make_shared< mtk::Mesh_Manager >();

        // XTK initialize
        std::shared_ptr< xtk::Model > pXTK = std::make_shared< xtk::Model >( tXTKParams( 0 )( 0 ) );
        
        // Set performer to HMR
        pHMR->set_performer( pBGMTK );
        
        // Set XTK Cooperations
        pXTK->set_geometry_engine( pGEN.get() );
        pXTK->set_input_performer( pBGMTK );
        pXTK->set_output_performer( pOutputMTK );

        // uniform initial refinement
        pHMR->perform_initial_refinement();

        // HMR finalize
        pHMR->perform();
        //
        pGEN->distribute_advs( pBGMTK->get_mesh_pair( 0 ),{} );

        // Output GEN fields, if requested
        pGEN->output_fields(pBGMTK->get_interpolation_mesh( 0 ));

        // XTK perform - decompose - enrich - ghost - multigrid
        pXTK->perform();




        

    }
}

}
