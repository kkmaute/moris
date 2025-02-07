/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_XTK_HMR_Enrich_Inclusion.cpp
 *
 */

// #include "catch.hpp"

// #include "cl_Param_List.hpp"

// #include "cl_XTK_Model.hpp"
// #include "moris_typedefs.hpp"
// #include "cl_MTK_Mesh_Manager.hpp"

// #include "cl_MTK_Vertex.hpp"    //MTK
// #include "cl_MTK_Cell.hpp"
// #include "cl_MTK_Enums.hpp"
// #include "cl_MTK_Mesh.hpp"

// #include "cl_MTK_Mesh_Manager.hpp"
// #include "cl_MTK_Integration_Mesh_STK.hpp"
// #include "cl_MTK_Interpolation_Mesh.hpp"
// #include "cl_MTK_Integration_Mesh.hpp"
// #include "cl_MTK_Writer_Exodus.hpp"
// #include "cl_XTK_Enriched_Integration_Mesh.hpp"

// #include "cl_Matrix.hpp"        //LINALG
// #include "linalg_typedefs.hpp"
// #include "fn_equal_to.hpp" // ALG/src

// #include "cl_HMR_Mesh_Interpolation.hpp"
// #include "cl_HMR.hpp"
// #include "cl_HMR_Background_Mesh.hpp" //HMR/src
// #include "cl_HMR_BSpline_Mesh_Base.hpp" //HMR/src
// #include "cl_HMR_Element.hpp" //HMR/src
// #include "cl_HMR_Factory.hpp" //HMR/src
// #include "cl_HMR_Field.hpp"
// #include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src
// #include "fn_PRM_HMR_Parameters.hpp"

// #include "cl_GEN_Circle.hpp"
// #include "cl_GEN_Plane.hpp"
// #include "cl_GEN_User_Defined_Field.hpp"
// #include "fn_norm.hpp"

// namespace UT_Setup
// {
// ParameterList
// set_inclusion_hmr_params()
// {

//     // background mesh parameters
//     std::string tNumElementsPerDir = "4,4";
//     std::string tDimensions        = "4,4";
//     std::string tOffSet            = "-2.0,-2.0";
//     std::string tSideSets          = "1,2,3,4";

//     int tLevelsetOrder       = 2;
//     int tDispOrder           = 1;

//     int tLevelsetInitialRef  = 1;
//     int tDispInitialRef      = 3;
//     int tRefineBuffer        = 0;

//     // note: pattern 0 - Levelset field  pattern 1 - displacement field
//     std::string tLagrangeOrder   = std::to_string(std::max(tLevelsetOrder,tDispOrder));
//     std::string tBsplineOrder    = std::to_string(tLevelsetOrder)      + "," + std::to_string(tDispOrder);
//     std::string tInitialRef      = std::to_string(tLevelsetInitialRef) + "," + std::to_string(tDispInitialRef);
//     std::string tLagrangePattern = tLevelsetInitialRef > tDispInitialRef ? "0" : "1";

//     ParameterList tParameterlist = prm::create_hmr_parameter_list();
//     tParameterlist.set( "number_of_elements_per_dimension", tNumElementsPerDir);
//     tParameterlist.set( "domain_dimensions",                tDimensions    );
//     tParameterlist.set( "domain_offset",                    tOffSet);
//     tParameterlist.set( "lagrange_output_meshes",           "0");
//     tParameterlist.set( "lagrange_orders",  tLagrangeOrder   );
//     tParameterlist.set( "lagrange_pattern", tLagrangePattern );
//     tParameterlist.set( "bspline_orders",   tBsplineOrder );
//     tParameterlist.set( "bspline_pattern",  "0,1" );
//     tParameterlist.set( "initial_refinement",         tInitialRef );
//     tParameterlist.set( "initial_refinement_pattern", "0,1" );
//     tParameterlist.set( "lagrange_to_bspline", "0,1");
//     tParameterlist.set( "refinement_buffer",  tRefineBuffer );
//     tParameterlist.set( "staircase_buffer",   tRefineBuffer );
//     tParameterlist.set( "use_number_aura", 1);
//     tParameterlist.set( "use_multigrid",  0 );
//     tParameterlist.set( "severity_level", 0 );

//         return tParameterlist
// }
// }
// namespace moris::xtk
// {
//     // Tests:
//     // Enrichment of an inclusion with multiple mesh interpolations to enrich
//     TEST_CASE("Enrichment of Inclusion","[XTK_ENR_INCL]")
//     {
//         if(par_size()<=1)
//         {
//             // HMR Parameter List
//             ParameterList tHMRParamList = UT_Setup::set_inclusion_hmr_params();

//             // create the HMR parameter file
//             tHMRParamList = std::make_shared< hmr::HMR >( tHMRParamList );

//         }
//     }
// }
