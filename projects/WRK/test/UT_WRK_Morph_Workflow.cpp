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
#include <cl_XTK_Enriched_Integration_Mesh.hpp>

namespace MorphTestWRK
{


struct InputMetaData
{
  bool mVerbose = false;
  std::string mWorkflow = "";
  size_t mSpatialDimensions = 3;
  
  moris::Cell< moris::ParameterList > mGeometryParameters;

  moris::hmr::Parameters mMorisHMRParameters;

  moris::Cell< moris::Cell< moris::ParameterList > > mXTKParameters;
  moris::Cell< moris::Cell< moris::ParameterList > > mSTKParameters;
  moris::Cell< moris::Cell< moris::ParameterList > > mHMRParameters;
  moris::Cell< moris::Cell< moris::ParameterList > > mGENParameters;
};

struct Performers
{
  std::shared_ptr< moris::hmr::HMR >            mCurrentHMR;
  std::shared_ptr< moris::ge::Geometry_Engine > mCurrentGEN;
  std::shared_ptr< mtk::Mesh_Manager >          mCurrentBGMTK;
  std::shared_ptr< mtk::Mesh_Manager >          mCurrentOutputMTK;
  std::shared_ptr< xtk::Model >                 mCurrentXTK;
  std::shared_ptr<mtk::Interpolation_Mesh_STK>  mSTKIpMesh;
  std::shared_ptr<mtk::Integration_Mesh_STK>    mSTKIgMesh;
};

struct MetaDataXTK
{
  // phase table member data
  int mNumPhases;
  int mDefaultPhase;

  // input parameters
  moris::uint mStep = 0;

  /// @brief output mesh index
  uint mOutputMeshIndex;

  // Cell Mesh Fields
  moris::Cell<moris::Matrix<DDRMat>> mCriteriaSensitivity;
  moris::Cell<std::string> mCriteriaSensitivityFieldName;

  // seeded field data
  std::string                  mOutputSeededMeshFile;
  bool                         mSeeded = false;
  std::string                  mSeededOpName;
  moris::Matrix<moris::DDRMat> mSeededField;

  // fixed design variables
  moris::Matrix<moris::IndexMat> mFixedBGNodesGD;

  // Fields on Cut Mesh
  std::vector<std::string> mNodeFieldsOnCutMesh;

  // Fields on BG Mesh
  std::vector<std::string> mNodeFieldsOnBGMesh;
};

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
  MorphTestWRK::InputMetaData & InputParameters
)
{
  InputParameters.mXTKParameters(0)(0).set("verbose", InputParameters.mVerbose);
  MorphTestWRK::setWorkFlowParams(InputParameters.mWorkflow,InputParameters.mXTKParameters);
  InputParameters.mXTKParameters(0)(0).set("output_file", "ut_xtk_bolted_bracket_mod.exo");
  InputParameters.mXTKParameters(0)(0).set("output_path", "./");
  InputParameters.mXTKParameters(0)(0).set("union_blocks", "block_1_n_p1,block_1_c_p1,block_2_n_p1,block_2_c_p1;block_1_n_p0,block_1_c_p0,block_2_n_p0,block_2_c_p0");
  InputParameters.mXTKParameters(0)(0).set("union_block_names", "block_1;block_2");
  InputParameters.mXTKParameters(0)(0).set("union_block_colors", "1;1");
  InputParameters.mXTKParameters(0)(0).set("union_side_sets", "ss_1_n_p1,ss_1_c_p1;ss_2_n_p1,ss_2_c_p1");
  InputParameters.mXTKParameters(0)(0).set("union_side_set_names", "ss_1;ss_2");
  InputParameters.mXTKParameters(0)(0).set("union_side_set_colors", "1;1");
  InputParameters.mXTKParameters(0)(0).set("deactivate_all_but_blocks", "block_1,block_2");
  InputParameters.mXTKParameters(0)(0).set("deactivate_all_but_side_sets", "ss_1,ss_2");
  InputParameters.mXTKParameters(0)(0).set("decompose", true);
  InputParameters.mXTKParameters(0)(0).set("decomposition_type", "conformal");
  InputParameters.mXTKParameters(0)(0).set("enrich", true);
}

void 
setParamsHMR(
  moris::Cell< moris::Cell< moris::ParameterList > > & aHMRParameters
)
{
  // "number_of_elements_per_dimension" keyword
  std::string tNumElemsStr = "1,1,1";
  // "domain_dimensions" keyword
  std::string tDomSizeStr = "1,1,1";
  // "domain_offset" keyword
  std::string tDomOffsetStr = "1,1,1";
  // "domain_sidesets" keyword
  std::string tSideSetStr = "ss_1,ss_2";
  // "state_refinement_level" keyword
  std::string tStateRefinement = "1,1,1";
  // "geom_refinement_level" keyword
  std::string tGeomRefinementLev = "1,1,1";
  //
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
setParamsHMR(
  MorphTestWRK::InputMetaData & aInputParameters
)
{
  // initialize moris data structures
  //
  moris::Matrix<DDLUMat> tElementsPerDimenson;
  moris::Matrix<DDRMat> tDimensons;
  moris::Matrix<DDRMat> tOffset;
  moris::Matrix<DDUMat> tSideSets;

  moris::Cell<enum Subdivision_Method> tDecompositionMethods(2);
  if ( aInputParameters.mSpatialDimensions == 2u )
  {
    tElementsPerDimenson = { {22}, {8} };
    tDimensons = { {6}, {2} };
    tOffset = { {-3.0}, {-1.0} };
    tSideSets = { {1},{2},{3},{4} };
    tDecompositionMethods( 0 ) = Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4;
    tDecompositionMethods( 1 ) = Subdivision_Method::C_TRI3;
  }
  else if ( aInputParameters.mSpatialDimensions == 3u )
  {
    tElementsPerDimenson = { {22}, {8}, {8} };
    tDimensons = { {6}, {2}, {2} };
    tOffset = { {-3.0}, {-1.0},{-1.0} };
    tSideSets = {{1},{2},{3},{4},{5},{6}};
    tDecompositionMethods( 0 ) = Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8;
    tDecompositionMethods( 1 ) = Subdivision_Method::C_HIERARCHY_TET4;
  }
  // set lagrange and bspline integration order
  //
  moris::uint tBsplineOrder  = 1;
  moris::uint tLagrangeOrder = tBsplineOrder;
  // set hmr parameters 
  //
  aInputParameters.mMorisHMRParameters.set_number_of_elements_per_dimension( tElementsPerDimenson );
  aInputParameters.mMorisHMRParameters.set_domain_dimensions( tDimensons );
  aInputParameters.mMorisHMRParameters.set_domain_offset( tOffset );
  aInputParameters.mMorisHMRParameters.set_bspline_truncation( true );
  aInputParameters.mMorisHMRParameters.set_output_meshes( {{ {0} }} );
  aInputParameters.mMorisHMRParameters.set_lagrange_orders  ( { {tLagrangeOrder} });
  aInputParameters.mMorisHMRParameters.set_lagrange_patterns({ {0} });
  aInputParameters.mMorisHMRParameters.set_bspline_orders   ( { {tBsplineOrder} } );
  aInputParameters.mMorisHMRParameters.set_bspline_patterns ( { {0} } );
  aInputParameters.mMorisHMRParameters.set_side_sets( tSideSets );
  aInputParameters.mMorisHMRParameters.set_refinement_buffer( 2 );
  aInputParameters.mMorisHMRParameters.set_staircase_buffer( 2 );
  aInputParameters.mMorisHMRParameters.set_lagrange_to_bspline_mesh( {{ {0} }});
}

void 
setParamsSTK(
  moris::Cell< moris::Cell< moris::ParameterList > > & aSTKParameters
)
{
  aSTKParameters(0)(0).set("input_file", "ut_xtk_bolted_bracket.exo");
  aSTKParameters(0)(0).set("input_path", "/home/maguilo/Morphorm/build/debug/moris/");
}

void 
setSwissCheeseParams(
  moris::ParameterList                               & aGEOParameters,
  moris::Cell< moris::Cell< moris::ParameterList > > & aGENParameters
)
{
  // initialize swiss cheese geometry
  //
  auto tGeometryID = aGEOParameters.get<int>("geometry_id");
  aGENParameters(1).push_back(moris::prm::create_swiss_cheese_slice_parameter_list());
  aGENParameters(1)(tGeometryID).set("left_bound",0.5);
  aGENParameters(1)(tGeometryID).set("right_bound",2.0);
  aGENParameters(1)(tGeometryID).set("top_bound",1.5);
  aGENParameters(1)(tGeometryID).set("bottom_bound",-1.5);
  aGENParameters(1)(tGeometryID).set("hole_x_semidiameter",0.211);
  aGENParameters(1)(tGeometryID).set("hole_y_semidiameter",0.211);
  aGENParameters(1)(tGeometryID).set("number_of_x_holes",3);
  aGENParameters(1)(tGeometryID).set("number_of_y_holes",3);
  aGENParameters(1)(tGeometryID).set("superellipse_exponent",4.0);         // Superellipse exponent
  aGENParameters(1)(tGeometryID).set("superellipse_scaling",0.25);       // Superellipse scaling
  aGENParameters(1)(tGeometryID).set("superellipse_regularization",0.0); // Superellipse regularization
}

void 
setParamsGEODiscretization(
  moris::ParameterList                               & aGEOParameters,
  moris::Cell< moris::Cell< moris::ParameterList > > & aGENParameters
)
{
  // set dicretization parameter 
  //
  auto tDesignType = aGEOParameters.get<bool>("design");
  if( tDesignType )
  {
    auto tGeometryID = aGEOParameters.get<int>("geometry_id");
    // Index of B-spline mesh to create level set field on (-1 = none)
    //
    aGENParameters(1)(tGeometryID).set("discretization_mesh_index", 0);  
    // Lower bound of level set field (if bspline_mesh_index >= 0)
    //
    aGENParameters(1)(tGeometryID).set("discretization_lower_bound", -3.0); 
    // Upper bound of level set field (if bspline_mesh_index >= 0)
    //
    aGENParameters(1)(tGeometryID).set("discretization_upper_bound", 3.0);  
  }
}

void 
setParamsGEO(
  moris::ParameterList                               & aGEOParameters,
  moris::Cell< moris::Cell< moris::ParameterList > > & aGENParameters
)
{
  // set geometry parameters 
  //
  auto tGeoID   = aGEOParameters.get<int>("geometry_id");
  auto tGeoType = aGEOParameters.get<std::string>("geometry_type");
  if( tGeoType == "swiss_cheese_slice" ){
    MorphTestWRK::setSwissCheeseParams(aGEOParameters,aGENParameters);
  }
  else {
    aGENParameters(1).push_back(moris::prm::create_geometry_parameter_list());
    aGENParameters(1)(tGeoID).set("constant_parameters","");
  }
}

void
createParamListXTK(
  moris::Cell< moris::Cell< moris::ParameterList > > & aXTKParameters
)
{
  aXTKParameters.resize(1);
  aXTKParameters(0).resize(1);
  aXTKParameters(0)(0) = moris::prm::create_xtk_parameter_list();  
}

void
createParamListSTK(
  moris::Cell< moris::Cell< moris::ParameterList > > & aSTKParameters
)
{
  aSTKParameters.resize(1);
  aSTKParameters(0).resize(1);
  aSTKParameters(0)(0) = moris::prm::create_stk_parameter_list();
  aSTKParameters(0)(0).insert("input_path","./");
}

void
createParamListHMR(
  moris::Cell< moris::Cell< moris::ParameterList > > & aHMRParameters
)
{
  aHMRParameters.resize(1);
  aHMRParameters(0).resize(1);
  aHMRParameters(0)(0) = moris::prm::create_hmr_parameter_list();
}

void
createParamListGEN(
  moris::Cell< moris::Cell< moris::ParameterList > > & aGENParameters
)
{
  // integrated quantities of interests and adv parameter definitions
  //
  aGENParameters.resize(3);
  aGENParameters(0).resize(1);
  aGENParameters(0)(0) = moris::prm::create_gen_parameter_list();
}

void 
createParamLists(
  MorphTestWRK::InputMetaData & aInputParameterLists
)
{
  MorphTestWRK::createParamListXTK(aInputParameterLists.mXTKParameters);
  MorphTestWRK::createParamListSTK(aInputParameterLists.mSTKParameters);
  MorphTestWRK::createParamListHMR(aInputParameterLists.mHMRParameters);
  MorphTestWRK::createParamListGEN(aInputParameterLists.mGENParameters);
}

void 
setParamLists(
  MorphTestWRK::InputMetaData & aInputParameters
)
{
  // set xtk parameters
  //
  MorphTestWRK::setParamsXTK(aInputParameters);
  // set hmr parameters
  //
  //MorphTestWRK::setParamsHMR(aInputParameters.mHMRParameters);
  MorphTestWRK::setParamsHMR(aInputParameters);
  // set stk parameters
  //
  MorphTestWRK::setParamsSTK(aInputParameters.mSTKParameters);
  // set gen parameters
  //
  auto tNumGeometries = aInputParameters.mGeometryParameters.size();
  for(decltype(tNumGeometries) tIndex = 0; tIndex < tNumGeometries; tIndex++) {
    MorphTestWRK::setParamsGEO(aInputParameters.mGeometryParameters(tIndex),aInputParameters.mGENParameters);
    MorphTestWRK::setParamsGEODiscretization(aInputParameters.mGeometryParameters(tIndex),aInputParameters.mGENParameters);
  }
}

void 
initBackgroundMeshHMR(
  MorphTestWRK::InputMetaData & aInputParameters,
  MorphTestWRK::Performers    & aPerformers
)
{
  // initialize hmr performer 
  aPerformers.mCurrentHMR = std::make_shared< moris::hmr::HMR >( aInputParameters.mMorisHMRParameters );
  // set performer to HMR
  aPerformers.mCurrentHMR->set_performer(aPerformers.mCurrentBGMTK);
  // uniform initial refinement
  aPerformers.mCurrentHMR->perform_initial_refinement();
  // hmr finalize
  aPerformers.mCurrentHMR->perform();  
}

void 
initBackgroundMeshSTK(
  MorphTestWRK::InputMetaData & aInputParameters,
  MorphTestWRK::Performers    & aPerformers
)
{
  auto tMeshPath = aInputParameters.mSTKParameters(0)(0).get<std::string>("input_path");
  tMeshPath = tMeshPath.back() == '/' ? tMeshPath : tMeshPath + "/";
  auto tMeshFile = aInputParameters.mSTKParameters(0)(0).get<std::string>("input_file");
  std::string tMeshFilePlusPath = tMeshPath + tMeshFile;
  aPerformers.mSTKIpMesh = std::make_shared<mtk::Interpolation_Mesh_STK>(tMeshFilePlusPath, nullptr);
  aPerformers.mSTKIgMesh = std::make_shared<mtk::Integration_Mesh_STK>(*aPerformers.mSTKIpMesh, nullptr);
  aPerformers.mCurrentBGMTK->register_mesh_pair(aPerformers.mSTKIpMesh.get(), aPerformers.mSTKIgMesh.get());
}

void 
initBackgroundMesh(
  MorphTestWRK::InputMetaData & aInputParameters,
  MorphTestWRK::Performers    & aPerformers
)
{
  if( aInputParameters.mWorkflow == "HMR" ) {
    MorphTestWRK::initBackgroundMeshHMR(aInputParameters,aPerformers);
  }
  else if( aInputParameters.mWorkflow == "STK" ) {
    MorphTestWRK::initBackgroundMeshSTK(aInputParameters,aPerformers);
  }
  else {
    MORIS_ERROR(0, "Invalid workflow type: STK or HMR");
  }
}

void 
generateModel(
  MorphTestWRK::InputMetaData & tInputParameters,
  MorphTestWRK::Performers    & aPerformers
)
{
  // initialize the background mesh
  //
  MorphTestWRK::initBackgroundMesh(tInputParameters, aPerformers);
  // initialize geometry engine
  //
  aPerformers.mCurrentGEN = 
    std::make_shared< moris::ge::Geometry_Engine >( tInputParameters.mGENParameters, nullptr );
  // distribute advs 
  //
  // fails in cl_MTK_Mesh_Core.cpp:843, get_coefficient_owners_of_node, the function is not implemented
  aPerformers.mCurrentGEN->distribute_advs( 
    aPerformers.mCurrentBGMTK->get_mesh_pair(0),{}, moris::mtk::EntityRank::NODE);
}

void 
regenerateModel(
  MorphTestWRK::MetaDataXTK & aMetaDataXTK,
  MorphTestWRK::Performers  & aPerformers
)
{
  // add criteria sensitivity fields to mesh (for visualization)
  //
  xtk::Enriched_Integration_Mesh & tEnrIGMesh = aPerformers.mCurrentXTK->get_enriched_integ_mesh(0);
  size_t tNumPDVs = aMetaDataXTK.mCriteriaSensitivity.size();
  for(size_t tIndex = 0; tIndex < tNumPDVs; tIndex++) {
    moris::moris_index tFieldIndex = 
      tEnrIGMesh.create_field(aMetaDataXTK.mCriteriaSensitivityFieldName(tIndex),moris::mtk::EntityRank::NODE);
    tEnrIGMesh.add_field_data(tFieldIndex,moris::mtk::EntityRank::NODE,aMetaDataXTK.mCriteriaSensitivity(tIndex));
  }
  // write the mesh 
  //
  moris::mtk::Writer_Exodus tWriter( aPerformers.mCurrentOutputMTK->get_integration_mesh(0) );
  std::string tOptIterStrg = std::to_string(aMetaDataXTK.mOutputMeshIndex);
  std::string tTempFile = "./xtk_temp2." + std::to_string( aMetaDataXTK.mOutputMeshIndex ) + ".exo";
  std::string tMeshFile = "./xtk_evolve.e-s." + std::string(4-tOptIterStrg.length(),'0') + tOptIterStrg;
  tWriter.write_mesh( /*path=*/"", tMeshFile, /*path=*/"", tTempFile);
  tWriter.set_time(aMetaDataXTK.mOutputMeshIndex);
  // write sensitivities to output mesh
  //
  tWriter.set_nodal_fields( aMetaDataXTK.mCriteriaSensitivityFieldName );
  for(size_t tIndex = 0; tIndex < tNumPDVs; tIndex++)
  {
    if(aMetaDataXTK.mCriteriaSensitivity(tIndex).numel() > 0)
    {
      tWriter.write_nodal_field(
        aMetaDataXTK.mCriteriaSensitivityFieldName(tIndex), aMetaDataXTK.mCriteriaSensitivity(tIndex)
      );
    }
  }
  aMetaDataXTK.mOutputMeshIndex++;
  tWriter.close_file();  
}

void 
initAppXTK(
  const MorphTestWRK::InputMetaData & aInputParameters,
        MorphTestWRK::Performers    & aPerformers
)
{
  // initialize xtk app
  //
  aPerformers.mCurrentOutputMTK = std::make_shared< moris::mtk::Mesh_Manager >();
  aPerformers.mCurrentXTK = std::make_shared< xtk::Model >( aInputParameters.mXTKParameters( 0 )( 0 ) );
  // set XTK Cooperations
  //
  aPerformers.mCurrentXTK->set_geometry_engine(  aPerformers.mCurrentGEN.get() );
  aPerformers.mCurrentXTK->set_input_performer(  aPerformers.mCurrentBGMTK );
  aPerformers.mCurrentXTK->set_output_performer( aPerformers.mCurrentOutputMTK ); 
  // set gen mesh metadata
  //
  aPerformers.mCurrentGEN->output_fields( aPerformers.mCurrentBGMTK->get_interpolation_mesh( 0 ) );
  aPerformers.mCurrentGEN->reset_mesh_information( aPerformers.mCurrentBGMTK->get_interpolation_mesh(0) );
  // XTK perform - decompose - enrich - ghost - multigrid
  //
  aPerformers.mCurrentXTK->perform(); // fails in cl_HMR_Mesh.cpp:1660 at get_bspline(). there is a FIXME comment?
  // Assign PDVs
  //
  aPerformers.mCurrentGEN->create_pdvs( aPerformers.mCurrentOutputMTK->get_mesh_pair(0) );
}

// ##########
// unit tests
// ##########

TEST_CASE( "WRK_morph_xtk_parse_params_test", "[WRK_morph_test]" )
{
  // create xtk parameter list
  //
  MorphTestWRK::InputMetaData tInputParameters;
  MorphTestWRK::createParamListXTK(tInputParameters.mXTKParameters);
  // set xtk params
  //
  tInputParameters.mVerbose = true;
  tInputParameters.mWorkflow = std::string("STK");
  MorphTestWRK::setParamsXTK(tInputParameters);
  // print to console
  //
  moris::prm::print(tInputParameters.mXTKParameters(0)(0));
}

TEST_CASE( "WRK_morph_stk_parse_params_test", "[WRK_morph_test]" )
{
  // parse workflow params
  //
  moris::Cell< moris::Cell< moris::ParameterList > > tSTKParameters;
  MorphTestWRK::createParamListSTK(tSTKParameters);
  MorphTestWRK::setParamsSTK(tSTKParameters);
  moris::prm::print(tSTKParameters(0)(0));
}

TEST_CASE( "WRK_morph_hmr_parse_params_test", "[WRK_morph_test]" )
{
  // parse workflow params
  //
  moris::Cell< moris::Cell< moris::ParameterList > > tHMRParameters;
  MorphTestWRK::createParamListHMR(tHMRParameters);
  MorphTestWRK::setParamsHMR(tHMRParameters);
  moris::prm::print(tHMRParameters(0)(0));
}

TEST_CASE( "WRK_morph_set_hmr_params_test", "[WRK_morph_test]" )
{
  // parse workflow params
  //
  MorphTestWRK::InputMetaData tInputParameters;
  MorphTestWRK::setParamsHMR(tInputParameters);
  tInputParameters.mMorisHMRParameters.print();
}

TEST_CASE( "WRK_morph_xtk_parse_gen_params_test", "[WRK_morph_test]" )
{
  // create core geometry engine params
  //
  moris::Cell< moris::Cell< moris::ParameterList > > tGENParameters;
  MorphTestWRK::createParamListGEN(tGENParameters);
  // set my geometry params
  //
  moris::ParameterList tGeometryParameters;
  tGeometryParameters.insert("design", true);
  tGeometryParameters.insert("geometry_id", 0);
  tGeometryParameters.insert("geometry_type", "swiss_cheese_slice");
  MorphTestWRK::setParamsGEO(tGeometryParameters,tGENParameters);
  // set discretization parameters
  //
  MorphTestWRK::setParamsGEODiscretization(tGeometryParameters,tGENParameters);
  // print to console
  //
  moris::prm::print(tGENParameters(0)(0));
  moris::prm::print(tGENParameters(1)(0));
}

TEST_CASE( "WRK_morph_xtk_initialize_hmr_params_test", "[WRK_morph_test]" )
{
  // create parameter lists
  //
  MorphTestWRK::InputMetaData tInputParameters;
  tInputParameters.mWorkflow = "HMR";
  MorphTestWRK::createParamLists(tInputParameters);
  // set parameters
  //
  moris::ParameterList tGeometryParameters;
  tGeometryParameters.insert("design", true);
  tGeometryParameters.insert("geometry_id", 0);
  tGeometryParameters.insert("geometry_type", "swiss_cheese_slice");
  tInputParameters.mGeometryParameters.push_back(tGeometryParameters);
  MorphTestWRK::setParamLists(tInputParameters);
  // initialize geometry engine
  //
  MorphTestWRK::Performers tPerformers;
  tPerformers.mCurrentGEN = std::make_shared< ge::Geometry_Engine >( tInputParameters.mGENParameters, nullptr );
  tPerformers.mCurrentBGMTK = std::make_shared< mtk::Mesh_Manager >();
  tPerformers.mCurrentOutputMTK = std::make_shared< mtk::Mesh_Manager >();
  // initialize the background mesh
  //
  MorphTestWRK::initBackgroundMesh(tInputParameters, tPerformers);
}

TEST_CASE( "WRK_morph_xtk_initialize_stk_params_test", "[WRK_morph_test]" )
{
  // create parameter lists
  //
  MorphTestWRK::InputMetaData tInputParameters;
  tInputParameters.mWorkflow = "STK";
  MorphTestWRK::createParamLists(tInputParameters);
  // set parameters
  //
  moris::ParameterList tGeometryParameters;
  tGeometryParameters.insert("design", true);
  tGeometryParameters.insert("geometry_id", 0);
  tGeometryParameters.insert("geometry_type", "swiss_cheese_slice");
  tInputParameters.mGeometryParameters.push_back(tGeometryParameters);
  MorphTestWRK::setParamLists(tInputParameters);
  // initialize geometry engine
  //
  MorphTestWRK::Performers tPerformers;
  tPerformers.mCurrentGEN = std::make_shared< ge::Geometry_Engine >( tInputParameters.mGENParameters, nullptr );
  tPerformers.mCurrentBGMTK = std::make_shared< mtk::Mesh_Manager >();
  tPerformers.mCurrentOutputMTK = std::make_shared< mtk::Mesh_Manager >();
  // initialize the background mesh
  //
  MorphTestWRK::initBackgroundMesh(tInputParameters, tPerformers);
}

TEST_CASE( "WRK_morph_xtk_generate_model_operation_hmr", "[WRK_morph_test]" )
{
  // create parameter lists
  //
  MorphTestWRK::InputMetaData tInputParameters;
  tInputParameters.mWorkflow = "HMR";
  MorphTestWRK::createParamLists(tInputParameters);
  // set parameters
  //
  moris::ParameterList tGeometryParameters;
  tGeometryParameters.insert("design", true);
  tGeometryParameters.insert("geometry_id", 0);
  tGeometryParameters.insert("geometry_type", "swiss_cheese_slice");
  tInputParameters.mGeometryParameters.push_back(tGeometryParameters);
  MorphTestWRK::setParamLists(tInputParameters);
  // initialize geometry engine
  //
  MorphTestWRK::Performers tPerformers;
  tPerformers.mCurrentGEN = std::make_shared< ge::Geometry_Engine >( tInputParameters.mGENParameters, nullptr );
  tPerformers.mCurrentBGMTK = std::make_shared< mtk::Mesh_Manager >();
  tPerformers.mCurrentOutputMTK = std::make_shared< mtk::Mesh_Manager >();
  // initialize background mesh 
  //
  MorphTestWRK::generateModel(tInputParameters,tPerformers);
  // initialize xtk app
  //
  MorphTestWRK::initAppXTK(tInputParameters,tPerformers);
  // test advs metadata
  //
  //moris::Matrix<moris::DDRMat> tADVs = tPerformers.mCurrentGEN->get_advs();
  //moris::real tTol = 1e-6; 
  //CHECK( tADVs.numel() == 1 );
  //CHECK( tADVs.min()   == 0 );
  //CHECK( tADVs.max()   == 1 );
  //CHECK( moris::equal_to(moris::norm(tADVs),1,tTol) );
}

TEST_CASE( "WRK_morph_xtk_generate_model_operation_stk", "[WRK_morph_test]" )
{
  // create parameter lists
  //
  MorphTestWRK::InputMetaData tInputParameters;
  tInputParameters.mWorkflow = "STK";
  MorphTestWRK::createParamLists(tInputParameters);
  // set parameters
  //
  moris::ParameterList tGeometryParameters;
  tGeometryParameters.insert("design", true);
  tGeometryParameters.insert("geometry_id", 0);
  tGeometryParameters.insert("geometry_type", "swiss_cheese_slice");
  tInputParameters.mGeometryParameters.push_back(tGeometryParameters);
  MorphTestWRK::setParamLists(tInputParameters);
  // initialize geometry engine
  //
  MorphTestWRK::Performers tPerformers;
  tPerformers.mCurrentGEN = std::make_shared< ge::Geometry_Engine >( tInputParameters.mGENParameters, nullptr );
  tPerformers.mCurrentBGMTK = std::make_shared< mtk::Mesh_Manager >();
  tPerformers.mCurrentOutputMTK = std::make_shared< mtk::Mesh_Manager >();
  // initialize background mesh 
  //
  MorphTestWRK::generateModel(tInputParameters,tPerformers);
  // initialize xtk app
  //
  //MorphTestWRK::initAppXTK(tInputParameters,tPerformers);
  // test advs metadata
  //
  //moris::Matrix<moris::DDRMat> tADVs = tPerformers.mCurrentGEN->get_advs();
  //moris::real tTol = 1e-6; 
  //CHECK( tADVs.numel() == 1 );
  //CHECK( tADVs.min()   == 0 );
  //CHECK( tADVs.max()   == 1 );
  //CHECK( moris::equal_to(moris::norm(tADVs),1,tTol) );
}

TEST_CASE( "WRK_morph_xtk_regenerate_model_operation_stk", "[WRK_morph_test]" )
{
  // set parameters
  //
  MorphTestWRK::InputMetaData tInputParameters;
  tInputParameters.mWorkflow = "STK";
  MorphTestWRK::createParamLists(tInputParameters);
  // initialize model
  //
  MorphTestWRK::Performers tPerformers;
  MorphTestWRK::generateModel(tInputParameters,tPerformers);
  // initialize xtk app
  //
  MorphTestWRK::initAppXTK(tInputParameters,tPerformers);
  // regenerate model
  //
  MorphTestWRK::MetaDataXTK tMetaDataXTK;
  MorphTestWRK::regenerateModel(tMetaDataXTK,tPerformers);
  // reinitialize model
  //
  MorphTestWRK::initAppXTK(tInputParameters,tPerformers);
  // test advs metadata
  //
  moris::Matrix<moris::DDRMat> tADVs = tPerformers.mCurrentGEN->get_advs();
  moris::real tTol = 1e-6; 
  CHECK( tADVs.numel() == 1 );
  CHECK( tADVs.min()   == 0 );
  CHECK( tADVs.max()   == 1 );
  CHECK( moris::equal_to(moris::norm(tADVs),1,tTol) );
}

TEST_CASE( "WRK_morph_xtk_sensitivity_operation", "[WRK_morph_test]" )
{
  // set parameters
  //
  MorphTestWRK::InputMetaData tInputParameters;
  tInputParameters.mWorkflow = "STK";
  MorphTestWRK::createParamLists(tInputParameters);
  // initialize model
  //
  MorphTestWRK::Performers tPerformers;
  MorphTestWRK::generateModel(tInputParameters,tPerformers);
  // initialize xtk app
  //
  MorphTestWRK::initAppXTK(tInputParameters,tPerformers);
}

} // namespace MorphTestWRK