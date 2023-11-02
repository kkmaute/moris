/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_WRK_Workflow.cpp
 *
 */

#include "catch.hpp"
#include "paths.hpp"

#include "op_times.hpp"
#include "op_equal_equal.hpp"

#include <fn_equal_to.hpp>
#include <fn_PRM_print.hpp>
#include <fn_PRM_HMR_Parameters.hpp>
#include <fn_PRM_GEN_Parameters.hpp>
#include <fn_PRM_XTK_Parameters.hpp>
#include <fn_PRM_STK_Parameters.hpp>

#include <cl_HMR.hpp>
#include "cl_MTK_Enums.hpp"
#include "cl_XTK_Model.hpp"
#include "cl_Library_Factory.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_GEN_Geometry_Engine.hpp"
#include "cl_Communication_Tools.hpp"
#include "cl_WRK_Workflow_HMR_XTK.hpp"
#include "cl_WRK_Performer_Manager.hpp"
#include "cl_MTK_Integration_Mesh_STK.hpp"
#include "cl_MTK_Interpolation_Mesh_STK.hpp"

namespace MorphTestWRK
{

void 
setWorkFlowParams(
  const std::string                                  & aWorkflowType,
  moris::Cell< moris::Cell< moris::ParameterList > > & aXTKParameters
)
{
  std::string tBasisRank = "";
  std::string tEnrMeshIndices = "";
  if (aWorkflowType == "HMR")
  {
    tBasisRank = "bspline";
    tEnrMeshIndices = "0,1";
  }
  else if (aWorkflowType == "STK")
  {
    tBasisRank = "node";
    tEnrMeshIndices = "0";
  }
  aXTKParameters(0)(0).set("basis_rank", tBasisRank);
  aXTKParameters(0)(0).set("enrich_mesh_indices", tEnrMeshIndices);
}

void 
setFixedBlockParams(
  moris::Cell< moris::Cell< moris::ParameterList > > & aXTKParameters
)
{
  // input keyword in the Morph xtkApp is "fixed_blocks"
  //
  auto tFixedBlocks = std::string("");
  if( !tFixedBlocks.empty() ) {
    aXTKParameters(0)(0).set("fixed_blocks",tFixedBlocks);
  }
  // input keyword in the Morph xtkApp is "fixed_block_phase"
  //
  auto tFixedBlocksPhase = std::string("");
  if( !tFixedBlocksPhase.empty() ){
    aXTKParameters(0)(0).set("fixed_block_phase",tFixedBlocksPhase);
  }
}

void 
setParamsXTK(
  const std::string                                  & aWorkflowType,
  moris::Cell< moris::Cell< moris::ParameterList > > & aXTKParameters,
  bool aVerbose = false
)
{
  aXTKParameters(0)(0).set("verbose", aVerbose);
  MorphTestWRK::setWorkFlowParams(aWorkflowType,aXTKParameters);

  aXTKParameters(0)(0).set("output_file", "ut_xtk_bolted_bracket_mod.exo");
  aXTKParameters(0)(0).set("output_path", "./");
  aXTKParameters(0)(0).set("union_blocks", "block_1_n_p1,block_1_c_p1,block_2_n_p1,block_2_c_p1;block_1_n_p0,block_1_c_p0,block_2_n_p0,block_2_c_p0");
  aXTKParameters(0)(0).set("union_block_names", "block_1;block_2");
  aXTKParameters(0)(0).set("union_block_colors", "1;1");
  aXTKParameters(0)(0).set("union_side_sets", "ss_1_n_p1,ss_1_c_p1;ss_2_n_p1,ss_2_c_p1");
  aXTKParameters(0)(0).set("union_side_set_names", "ss_1;ss_2");
  aXTKParameters(0)(0).set("union_side_set_colors", "1;1");
  aXTKParameters(0)(0).set("deactivate_all_but_blocks", "block_1,block_2");
  aXTKParameters(0)(0).set("deactivate_all_but_side_sets", "ss_1,ss_2");
  aXTKParameters(0)(0).set("decompose", true);
  aXTKParameters(0)(0).set("decomposition_type", "conformal");
  aXTKParameters(0)(0).set("enrich", true);
}

void 
setParamsHMR(
  moris::Cell< moris::Cell< moris::ParameterList > > & aHMRParameters
)
{
  // "number_of_elements_per_dimension" keyword
  std::string tNumElemsStr = "1;1;1";
  // "domain_dimensions" keyword
  std::string tDomSizeStr = "1;1;1";
  // "domain_offset" keyword
  std::string tDomOffsetStr = "1;1;1";
  // "domain_sidesets" keyword
  std::string tSideSetStr = "ss_1;ss_2";
  // "state_refinement_level" keyword
  std::string tStateRefinement = "1;1;1";
  // "geom_refinement_level" keyword
  std::string tGeomRefinementLev = "1;1;1";

  std::string tLagrangeOrder = "1";
  std::string tBsplineOrder = "2,1";
  std::string tLagrangePattern = "1";
  std::string tInitialRef = tGeomRefinementLev + "," + tStateRefinement;
  aHMRParameters(0)(0).set("number_of_elements_per_dimension", tNumElemsStr);
  aHMRParameters(0)(0).set("domain_dimensions", tDomSizeStr);
  aHMRParameters(0)(0).set("domain_offset", tDomOffsetStr);
  aHMRParameters(0)(0).set("domain_sidesets", tSideSetStr);
  aHMRParameters(0)(0).set("lagrange_output_meshes", std::string("0"));
  aHMRParameters(0)(0).set("lagrange_orders", tLagrangeOrder);
  aHMRParameters(0)(0).set("lagrange_pattern", tLagrangePattern);
  aHMRParameters(0)(0).set("bspline_orders", tBsplineOrder);
  aHMRParameters(0)(0).set("bspline_pattern", std::string("0,1"));
  aHMRParameters(0)(0).set("lagrange_to_bspline", "0,1");
  aHMRParameters(0)(0).set("truncate_bsplines", 1);
  aHMRParameters(0)(0).set("refinement_buffer", 1);
  aHMRParameters(0)(0).set("staircase_buffer", 1);
  aHMRParameters(0)(0).set("initial_refinement", tInitialRef);
  aHMRParameters(0)(0).set("initial_refinement_pattern", "0,1");
  aHMRParameters(0)(0).set("use_number_aura", 1);
  aHMRParameters(0)(0).set("use_multigrid", 0);
  aHMRParameters(0)(0).set("severity_level", 0);
  aHMRParameters(0)(0).set("write_lagrange_output_mesh", "./HMR_Output_Lagrange.e");
  aHMRParameters(0)(0).set("use_advanced_T_matrix_scheme", 1);
}

void 
setParamsSTK(
  moris::Cell< moris::Cell< moris::ParameterList > > & aSTKParameters
)
{
  aSTKParameters(0)(0).set("input_file", "ut_xtk_bolted_bracket.exo");
}

void 
setSwissCheese(
  moris::Cell< moris::Cell< moris::ParameterList > > & aGENParameters
)
{
  // initialize swiss cheese geometry
  //
  size_t tGeoCount = 0;
  aGENParameters(1).push_back(moris::prm::create_swiss_cheese_slice_parameter_list());
  aGENParameters(1)(tGeoCount).set("left_bound","0.5");
  aGENParameters(1)(tGeoCount).set("right_bound","2.0");
  aGENParameters(1)(tGeoCount).set("top_bound","1.5");
  aGENParameters(1)(tGeoCount).set("bottom_bound","-1.5");
  aGENParameters(1)(tGeoCount).set("hole_x_semidiameter","0.211");
  aGENParameters(1)(tGeoCount).set("hole_y_semidiameter","0.211");
  aGENParameters(1)(tGeoCount).set("number_of_x_holes","3");
  aGENParameters(1)(tGeoCount).set("number_of_y_holes","3");
  //
  aGENParameters(1)(tGeoCount).set("superellipse_exponent","4");         // Superellipse exponent
  aGENParameters(1)(tGeoCount).set("superellipse_scaling","0.25");       // Superellipse scaling
  aGENParameters(1)(tGeoCount).set("superellipse_regularization","0.0"); // Superellipse regularization
}

void 
setParamsDiscretization(
  const size_t                                             & aGeometryID,
        moris::Cell< moris::Cell< moris::ParameterList > > & aGENParameters
)
{
  // set dicretization parameter 
  //
  auto tDesignType = aGENParameters(0)(aGeometryID).get<bool>("design");
  if( tDesignType )
  {
    // Index of B-spline mesh to create level set field on (-1 = none)
    //
    aGENParameters(1)(aGeometryID).set("discretization_mesh_index", 0);  
    // Lower bound of level set field (if bspline_mesh_index >= 0)
    //
    aGENParameters(1)(aGeometryID).set("discretization_lower_bound", -3.0); 
    // Upper bound of level set field (if bspline_mesh_index >= 0)
    //
    aGENParameters(1)(aGeometryID).set("discretization_upper_bound", 3.0);  
  }
}

void 
setParamsGeomBase(
  moris::Cell< moris::Cell< moris::ParameterList > > & aGENParameters
)
{
  // simple use case: only one geometry defined and two phases (solid and void)
  //
  aGENParameters(0)(0).set("IQI_types", "Strain_Energy,Volume");
  aGENParameters(0)(0).set("isocontour_threshold", "0.0"); // Level-set isocontour level
  aGENParameters(0)(0).set("isocontour_tolerance", "1e-12"); // Level-set isocontour level
  aGENParameters(0)(0).set("output_mesh_file", "./gen_fields.exo");
  aGENParameters(0)(0).set("num_phases", "2");
  aGENParameters(0)(0).set("default_phase", "-1");
  // geometry parameters
  //
  aGENParameters(0)(0).set("design", "true");
  aGENParameters(0)(0).set("type", "swiss_cheese_slice");
}

void 
setParamsGEN(
  moris::Cell< moris::Cell< moris::ParameterList > > & aGENParameters,
  bool aVerbose = false
)
{
  // set fundamental geometry parameters
  //
  MorphTestWRK::setParamsGeomBase(aGENParameters);
  //
  size_t tGeoCount = 0;
  if (aVerbose) {
    moris::prm::print(aGENParameters(0)(tGeoCount));
  }
  // set initial guess parameters 
  //
  auto tGeomType = aGENParameters(0)(0).get<std::string>("type");
  if( tGeomType == "swiss_cheese_slice" ){
    MorphTestWRK::setSwissCheese(aGENParameters);
  }
  else {
    aGENParameters(1).push_back(moris::prm::create_geometry_parameter_list());
    aGENParameters(1)(tGeoCount).set("type",tGeomType);
    aGENParameters(1)(tGeoCount).set("constant_parameters","");
  }
  // set discretization parameters
  //
  MorphTestWRK::setParamsDiscretization(tGeoCount,aGENParameters);
  if (aVerbose) {
    moris::prm::print(aGENParameters(1)(tGeoCount));
  }
}

std::shared_ptr< moris::mtk::Mesh_Manager > 
registerBackGroundMeshPair(
   moris::Cell< moris::Cell< moris::ParameterList > > & aSTKParameters
)
{
  // create background mesh metadata manager
  //
  std::shared_ptr< moris::mtk::Mesh_Manager > tCurrentBGMTK = 
    std::make_shared< moris::mtk::Mesh_Manager >();
  // construct the stk
  std::string tMeshFile = aSTKParameters(0)(0).get<std::string>("input_file");
  auto tSTKIpMesh = std::make_shared<moris::mtk::Interpolation_Mesh_STK>( tMeshFile, nullptr );
  auto tSTKIgMesh = std::make_shared<moris::mtk::Integration_Mesh_STK> ( *tSTKIpMesh.get(), nullptr );
  tCurrentBGMTK->register_mesh_pair(tSTKIpMesh.get(),tSTKIgMesh.get());
  return tCurrentBGMTK;
}

// ##########
// unit tests
// ##########

TEST_CASE( "WRK_morph_test", "[moris],[WRK_morph_xtk_parse_params_test]" )
{
  moris::Cell< moris::Cell< moris::ParameterList > > tXTKParameters;
  tXTKParameters.resize(1);
  tXTKParameters(0).resize(1);
  tXTKParameters(0)(0) = moris::prm::create_xtk_parameter_list();
  // set xtk params
  //
  auto tVerbose = true;
  auto tWorkflow = std::string("STK");
  MorphTestWRK::setParamsXTK(tWorkflow,tXTKParameters,tVerbose);
  // print to console
  //
  if (tVerbose) {
    moris::prm::print(tXTKParameters(0)(0));
  }
}

TEST_CASE( "WRK_morph_test", "[moris],[WRK_morph_stk_parse_params_test]" )
{
  // parse workflow params
  //
  moris::Cell< moris::Cell< moris::ParameterList > > tSTKParameters;
  tSTKParameters.resize(1);
  tSTKParameters(0).resize(1);
  tSTKParameters(0)(0) = moris::prm::create_stk_parameter_list();
  MorphTestWRK::setParamsSTK(tSTKParameters);
}

TEST_CASE( "WRK_morph_test", "[moris],[WRK_morph_hmr_parse_params_test]" )
{
  // parse workflow params
  //
  moris::Cell< moris::Cell< moris::ParameterList > > tHMRParameters;
  tHMRParameters.resize(1);
  tHMRParameters(0).resize(1);
  tHMRParameters(0)(0) = moris::prm::create_hmr_parameter_list();
  MorphTestWRK::setParamsHMR(tHMRParameters);
  // print to console
  //
  auto tVerbose = true;
  if ( tVerbose ) {
    moris::prm::print(tHMRParameters(0)(0));
  }
}

TEST_CASE( "WRK_morph_test", "[moris],[WRK_morph_xtk_parse_gen_params_test]" )
{
  // parse geometry params
  //
  moris::Cell< moris::Cell< moris::ParameterList > > tGENParameters;
  tGENParameters.resize(3);
  tGENParameters(0).resize(1);
  tGENParameters(0)(0) = moris::prm::create_gen_parameter_list();
  MorphTestWRK::setParamsGEN(tGENParameters);
}

TEST_CASE( "WRK_morph_test", "[moris],[WRK_morph_xtk_initialize_hmr_params_test]" )
{
  // initialize the background mesh
  // this happens here so we can export the scalar field data map
  std::shared_ptr< moris::hmr::HMR > tCurrentHMR;
  auto aCurrentBGMTK = std::make_shared<moris::mtk::Mesh_Manager>();
  // Set performer to HMR
  tCurrentHMR->set_performer(aCurrentBGMTK);
  // uniform initial refinement
  tCurrentHMR->perform_initial_refinement();
  // HMR finalize
  tCurrentHMR->perform();
}

TEST_CASE( "WRK_morph_test", "[moris],[WRK_morph_xtk_initialize_stk_params_test]" )
{
  // parse parameters
  //
  moris::Cell< moris::Cell< moris::ParameterList > > tSTKParameters;
  tSTKParameters.resize(1);
  tSTKParameters(0).resize(1);
  tSTKParameters(0)(0) = moris::prm::create_stk_parameter_list();
  MorphTestWRK::setParamsSTK(tSTKParameters);
  // construct stk workflow
  //
  std::shared_ptr< moris::mtk::Mesh_Manager > tCurrentBGMTK = 
    MorphTestWRK::registerBackGroundMeshPair(tSTKParameters);
}

TEST_CASE( "WRK_morph_test", "[moris],[WRK_morph_xtk_generate_model_operation_stk_not_initialize]" )
{
  // Set STK parameters
  //
  moris::Cell< moris::Cell< moris::ParameterList > > tSTKParameters;
  tSTKParameters.resize(1);
  tSTKParameters(0).resize(1);
  tSTKParameters(0)(0) = moris::prm::create_stk_parameter_list();
  MorphTestWRK::setParamsSTK(tSTKParameters);
  // register background mesh interpolation-integration metadata pair
  std::shared_ptr< moris::mtk::Mesh_Manager > tCurrentBGMTK = 
    MorphTestWRK::registerBackGroundMeshPair(tSTKParameters);
  // set geometry params
  //
  moris::Cell< moris::Cell< moris::ParameterList > > tGENParameters;
  tGENParameters.resize(3);
  tGENParameters(0).resize(1);
  tGENParameters(0)(0) = moris::prm::create_gen_parameter_list();
  MorphTestWRK::setParamsGEN(tGENParameters);
  // initialize geometry engine
  //
  auto tCurrentGEN = std::make_shared< moris::ge::Geometry_Engine >( tGENParameters, nullptr );
  // distribute advs 
  //
  tCurrentGEN->distribute_advs( tCurrentBGMTK->get_mesh_pair( 0 ),{}, moris::mtk::EntityRank::NODE );
  // test advs metadata
  //
  moris::Matrix<moris::DDRMat> tADVs = tCurrentGEN->get_advs();
  moris::real tTol = 1e-6; 
  CHECK( tADVs.numel() == 1 );
  CHECK( tADVs.min()   == 0 );
  CHECK( tADVs.max()   == 1 );
  CHECK( moris::equal_to(moris::norm(tADVs),1,tTol) );
  // initialize xtk param list
  //
  moris::Cell< moris::Cell< moris::ParameterList > > tXTKParameters;
  tXTKParameters.resize(1);
  tXTKParameters(0).resize(1);
  tXTKParameters(0)(0) = moris::prm::create_xtk_parameter_list();
  // set xtk params
  //
  auto tWorkflow = std::string("STK");
  MorphTestWRK::setParamsXTK(tWorkflow,tXTKParameters);
  // initialize xtk app
  //
  auto tCurrentOutputMTK = std::make_shared< moris::mtk::Mesh_Manager >();
  auto tCurrentXTK = std::make_shared< xtk::Model >( tXTKParameters( 0 )( 0 ) );
  // set XTK Cooperations
  //
  tCurrentXTK->set_geometry_engine(  tCurrentGEN.get() );
  tCurrentXTK->set_input_performer(  tCurrentBGMTK );
  tCurrentXTK->set_output_performer( tCurrentOutputMTK ); 
  // set gen mesh metadata
  //
  tCurrentGEN->output_fields( tCurrentBGMTK->get_interpolation_mesh( 0 ) );
  tCurrentGEN->reset_mesh_information( tCurrentBGMTK->get_interpolation_mesh(0) );
  // XTK perform - decompose - enrich - ghost - multigrid
  //
  tCurrentXTK->perform();
  // Assign PDVs
  //
  tCurrentGEN->create_pdvs( tCurrentOutputMTK->get_mesh_pair(0) );
}

} // namespace MorphTestWRK