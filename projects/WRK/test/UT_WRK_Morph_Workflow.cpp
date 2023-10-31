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

#include "cl_Communication_Tools.hpp"
#include "paths.hpp"

#include "op_times.hpp"
#include "op_equal_equal.hpp"

#include <fn_PRM_HMR_Parameters.hpp>
#include <fn_PRM_GEN_Parameters.hpp>
#include <fn_PRM_XTK_Parameters.hpp>
#include <fn_PRM_STK_Parameters.hpp>
#include <fn_PRM_print.hpp>

#include "cl_WRK_Performer_Manager.hpp"
#include "cl_WRK_Workflow_HMR_XTK.hpp"
#include "cl_Library_Factory.hpp"
#include "cl_Communication_Tools.hpp"

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
setParamsGEN(
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

  auto tGeomType = aGENParameters(0)(0).get<std::string>("type");
  if( tGeomType == "swiss_cheese_slice" ){
    MorphTestWRK::setSwissCheese(aGENParameters);
  }
}

TEST_CASE( "WRK_morph_test", "[moris],[WRK_morph_parse_params_test]" )
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
  // parse workflow params
  //
  moris::Cell< moris::Cell< moris::ParameterList > > tSTKParameters;
  moris::Cell< moris::Cell< moris::ParameterList > > tHMRParameters;
  if( tWorkflow == "HMR"){
    tHMRParameters.resize(1);
    tHMRParameters(0).resize(1);
    tHMRParameters(0)(0) = moris::prm::create_hmr_parameter_list();
    MorphTestWRK::setParamsHMR(tHMRParameters);
    if (tVerbose) {
      moris::prm::print(tHMRParameters(0)(0));
    }
  }
  else if (tWorkflow == "STK") {
    tSTKParameters.resize(1);
    tSTKParameters(0).resize(1);
    tSTKParameters(0)(0) = moris::prm::create_stk_parameter_list();
    MorphTestWRK::setParamsSTK(tSTKParameters);
  }
  // parse geometry params
  //
  moris::Cell< moris::Cell< moris::ParameterList > > tGENParameters;
  MorphTestWRK::setParamsGEN(tGENParameters);
}

} // namespace MorphTestWRK