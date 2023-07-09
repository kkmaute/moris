/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Vessel_2D.cpp
 *
 */

#include <string>
#include <iostream>
#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_MSI_Equation_Object.hpp"
#include "cl_TSA_Time_Solver.hpp"
#include "cl_DLA_Solver_Interface.hpp"
#include "cl_DLA_Linear_Solver_Aztec.hpp"
#include "fn_PRM_FEM_Parameters.hpp"
#include "fn_PRM_MSI_Parameters.hpp"
#include "fn_PRM_SOL_Parameters.hpp"
#include "fn_PRM_VIS_Parameters.hpp"
#include "fn_PRM_HMR_Parameters.hpp"
#include "fn_PRM_GEN_Parameters.hpp"
#include "fn_PRM_XTK_Parameters.hpp"
#include "fn_PRM_OPT_Parameters.hpp"
#include "fn_equal_to.hpp"

#include "AztecOO.h"

#ifdef  __cplusplus
extern "C"
{
#endif
//------------------------------------------------------------------------------
namespace moris
{

// Geometry Parameters
moris::real tCenterX         = 0.5;   /* y bottom plane (m) */
moris::real tCenterY         = 0.5;   /* y bottom plane (m) */
moris::real tOuterRad        = 0.4;   /* y top plane    (m) */
moris::real tInnerRad        = 0.385; /* y top plane    (m) */

    // Constant function for properties
    void Func_Const( moris::Matrix<
            moris::DDRMat >                                & aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
            moris::fem::Field_Interpolator_Manager         * aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    void Func_Initial_Condition(
            moris::Matrix< moris::DDRMat >                 & aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
            moris::fem::Field_Interpolator_Manager         * aFIManager )
    {
        aPropMatrix = {{ 489.0 }};
    }

    bool Output_Criterion( moris::tsa::Time_Solver * aTimeSolver )
    {
        return true;
    }

    moris::real Outer_Ring(
            const moris::Matrix< DDRMat >     & aCoordinates,
            const moris::Cell< moris::real* > & aGeometryParameters )
    {
        moris::real aReturnValue = tOuterRad - std::sqrt( std::pow(aCoordinates(0) - tCenterX,2.0) + std::pow(aCoordinates(1) - tCenterY,2.0));
        return std::abs(aReturnValue)<1e-8 ? 1e-8 : aReturnValue;
    }

    moris::real Inner_Ring(
            const moris::Matrix< DDRMat >     & aCoordinates,
            const moris::Cell< moris::real* > & aGeometryParameters )
    {
        moris::real aReturnValue = tInnerRad - std::sqrt( std::pow(aCoordinates(0) - tCenterX,2.0) + std::pow(aCoordinates(1) - tCenterY,2.0));
        return std::abs(aReturnValue)<1e-8 ? 1e-8 : aReturnValue;
    }

    moris::Matrix< DDRMat > Func_Sensitivity(
            const moris::Matrix< DDRMat >     & aCoordinates,
            const moris::Cell< moris::real* > & aGeometryParameters )
    {
        moris::Matrix< DDRMat > aReturnValue;
        return aReturnValue;
    }

    void OPTParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_opt_problem_parameter_list();

        tParameterlist( 0 )( 0 ).set( "is_optimization_problem", false);
    }

    void HMRParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_hmr_parameter_list();

    tParameterlist( 0 )( 0 ).set( "number_of_elements_per_dimension", "32,32");
        tParameterlist( 0 )( 0 ).set( "domain_dimensions",                "1,1");
        tParameterlist( 0 )( 0 ).set( "domain_offset",                    "0,0" );
        tParameterlist( 0 )( 0 ).set( "domain_sidesets",                  "1,2,3,4");
        tParameterlist( 0 )( 0 ).set( "lagrange_output_meshes",           "0");

        tParameterlist( 0 )( 0 ).set( "lagrange_orders",  std::string( "1" ));
        tParameterlist( 0 )( 0 ).set( "lagrange_pattern", std::string( "0" ));
        tParameterlist( 0 )( 0 ).set( "bspline_orders",   std::string( "1" ));
        tParameterlist( 0 )( 0 ).set( "bspline_pattern",  std::string( "0" ));

        tParameterlist( 0 )( 0 ).set( "lagrange_to_bspline", "0" );

        tParameterlist( 0 )( 0 ).set( "truncate_bsplines",  1 );
        tParameterlist( 0 )( 0 ).set( "refinement_buffer",  2 );
        tParameterlist( 0 )( 0 ).set( "staircase_buffer",   2 );
        tParameterlist( 0 )( 0 ).set( "initial_refinement", "1" );
        tParameterlist( 0 )( 0 ).set( "initial_refinement_pattern", "0" );

        tParameterlist( 0 )( 0 ).set( "use_number_aura", 1);

        tParameterlist( 0 )( 0 ).set( "use_multigrid",  0 );
        tParameterlist( 0 )( 0 ).set( "severity_level", 0 );

        tParameterlist( 0 )( 0 ).set( "adaptive_refinement_level", 4 );
    }

    void XTKParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_xtk_parameter_list();
        tParameterlist( 0 )( 0 ).set( "decompose",                 true );
        tParameterlist( 0 )( 0 ).set( "decomposition_type",        "conformal" );
        tParameterlist( 0 )( 0 ).set( "enrich",                    true );
        tParameterlist( 0 )( 0 ).set( "basis_rank",                "bspline" );
        tParameterlist( 0 )( 0 ).set( "enrich_mesh_indices",       "0" );
        tParameterlist( 0 )( 0 ).set( "ghost_stab",                true );
        tParameterlist( 0 )( 0 ).set( "multigrid",                 false );
        tParameterlist( 0 )( 0 ).set( "verbose",                   true );
        tParameterlist( 0 )( 0 ).set( "print_enriched_ig_mesh",    true );
        tParameterlist( 0 )( 0 ).set( "exodus_output_XTK_ig_mesh", false );
        tParameterlist( 0 )( 0 ).set( "high_to_low_dbl_side_sets", true );
    }

    void GENParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 3 );
        tParameterlist( 0 ).resize( 1 );

        // Main GEN parameter list
        tParameterlist( 0 )( 0 ) = prm::create_gen_parameter_list();
    tParameterlist( 0 )( 0 ).set( "HMR_refinements", 2 );

        // init geometry counter
        uint tGeoCounter = 0;

        // Geometry parameter lists
        tParameterlist( 1 ).push_back( prm::create_user_defined_geometry_parameter_list() );
        tParameterlist( 1 )( tGeoCounter ).set( "field_function_name", "Outer_Ring");
        tParameterlist( 1 )( tGeoCounter ).set( "sensitivity_function_name", "Func_Sensitivity");
        tParameterlist( 1 )( tGeoCounter ).set( "constant_parameters", "");
        tGeoCounter++;

        tParameterlist( 1 ).push_back( prm::create_user_defined_geometry_parameter_list() );
        tParameterlist( 1 )( tGeoCounter ).set( "field_function_name", "Inner_Ring");
        tParameterlist( 1 )( tGeoCounter ).set( "sensitivity_function_name", "Func_Sensitivity");
        tParameterlist( 1 )( tGeoCounter ).set( "constant_parameters", "");
        tGeoCounter++;
    }

    void FEMParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterList )
    {
        // create a cell of cell of parameter list for fem
        tParameterList.resize( 5 );

        //------------------------------------------------------------------------------
        // init property counter
        uint tPropCounter = 0;

        // create parameter list for outer ring

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ) = prm::create_property_parameter_list();
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropDensity_Outer" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      "1.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropHeatCapacity_Outer" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      "1.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropConductivity_Outer" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      "0.01" );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const" );
        tPropCounter++;

        // create parameter list for inner circle

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ) = prm::create_property_parameter_list();
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropDensity_Inner" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      "1.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropHeatCapacity_Inner" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      "1.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropConductivity_Inner" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      "0.01" );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropLatentHeat_Inner" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      "250.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropPCTemp_Inner" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      "500." );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropPhaseState_Inner" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      "2.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropPCconst_Inner" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      "20.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const" );
        tPropCounter++;

        // flux on outer ring

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropImposedFlux" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      "100.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const" );
        tPropCounter++;

        // algorithmic parameters

        // create parameter list for property 6
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropWeightCurrent" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      "100.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const" );
        tPropCounter++;

        // create parameter list for property 7
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropWeightPrevious" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      "100.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const" );
        tPropCounter++;

        // create parameter list for property 8
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropInitialCondition" );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Initial_Condition" );
        tPropCounter++;

        //------------------------------------------------------------------------------
        // init CM counter
        uint tCMCounter = 0;

        // create parameter list for constitutive model - outer ring
        tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_name", "CMDiffusion_Outer" );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::DIFF_LIN_ISO ) );
        tParameterList( 1 )( tCMCounter ).set( "dof_dependencies",  std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( 1 )( tCMCounter ).set( "properties",
                "PropConductivity_Outer , Conductivity;"
                "PropDensity_Outer      , Density;"
                "PropHeatCapacity_Outer , Heat_Capacity" );
        tCMCounter++;

        // create parameter list for constitutive model - inner circle
        tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_name", "CMDiffusion_Inner" );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::DIFF_LIN_ISO_PC ) );
        tParameterList( 1 )( tCMCounter ).set( "dof_dependencies",  std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( 1 )( tCMCounter ).set( "properties",
                "PropConductivity_Inner, Conductivity;"
                "PropDensity_Inner     , Density;"
                "PropHeatCapacity_Inner, Heat_Capacity;"
                "PropLatentHeat_Inner  , Latent_Heat;"
                "PropPCTemp_Inner      , PC_Temp;"
                "PropPhaseState_Inner  , Phase_State_Function;"
                "PropPCconst_Inner     , Phase_Change_Const"    );
        tCMCounter++;

        //------------------------------------------------------------------------------
        // init SP counter
        uint tSPCounter = 0;

        // create parameter list for stabilization parameter 2
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name",      "SPGPTemp_Outer" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type",      static_cast< uint >( fem::Stabilization_Type::GHOST_DISPL ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters",     "0.01" );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties",       "PropConductivity_Outer,Material" );
        tSPCounter++;

        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name",      "SPGPTemp_Inner" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type",      static_cast< uint >( fem::Stabilization_Type::GHOST_DISPL ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters",     "0.01" );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties",       "PropConductivity_Inner,Material" );
        tSPCounter++;

        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name",      "SPGGLSDiffusion_Inner" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type",      static_cast< uint >( fem::Stabilization_Type::GGLS_DIFFUSION_PC ) );
        tParameterList( 2 )( tSPCounter ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties",
                "PropConductivity_Inner , Conductivity;"
                "PropDensity_Inner      , Density;"
                "PropHeatCapacity_Inner , Heat_Capacity;"
                "PropLatentHeat_Inner   , Latent_Heat;"
                "PropPCTemp_Inner       , PC_Temp;"
                "PropPhaseState_Inner   , Phase_State_Function;"
                "PropPCconst_Inner      , Phase_Change_Const"    );
        tSPCounter++;

        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name",  "SPInterfaceNitsche" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type",  static_cast< uint >( fem::Stabilization_Type::NITSCHE_INTERFACE ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", "100.0" );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties",   "PropConductivity_Outer,Material" );
        tParameterList( 2 )( tSPCounter ).set( "follower_properties",    "PropConductivity_Inner,Material" );
        tSPCounter++;

        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name",  "SPInterfaceLeaderWeight" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type",  static_cast< uint >( fem::Stabilization_Type::LEADER_WEIGHT_INTERFACE ) );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties",   "PropConductivity_Outer,Material" );
        tParameterList( 2 )( tSPCounter ).set( "follower_properties",    "PropConductivity_Inner,Material" );
        tSPCounter++;

        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name",  "SPInterfaceFollowerWeight" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type",  static_cast< uint >( fem::Stabilization_Type::FOLLOWER_WEIGHT_INTERFACE ) );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties",   "PropConductivity_Outer,Material" );
        tParameterList( 2 )( tSPCounter ).set( "follower_properties",    "PropConductivity_Inner,Material" );
        tSPCounter++;

        //------------------------------------------------------------------------------
        // init IWG counter
        uint tIWGCounter = 0;

        // create parameter  for outer ring
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWGDiffusionBulk" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::SPATIALDIFF_BULK ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies",    "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMDiffusion_Outer,Diffusion" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             "HMR_dummy_n_p2,HMR_dummy_c_p2" );
        tIWGCounter++;

        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWGDiffusionBulk" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::SPATIALDIFF_BULK ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies",    "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMDiffusion_Inner,Diffusion" );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   "SPGGLSDiffusion_Inner,GGLS_Param" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             "HMR_dummy_n_p3,HMR_dummy_c_p3" );
        tIWGCounter++;

        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWGInletFlux" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::SPATIALDIFF_NEUMANN ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies",    "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties",          "PropImposedFlux,Neumann" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             "iside_b0_2_b1_0" );
        tIWGCounter++;

		tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
		tParameterList( 3 )( tIWGCounter ).set( "IWG_name", 				  "IWGGPTemp" );
		tParameterList( 3 )( tIWGCounter ).set( "IWG_type", 				  static_cast< uint >( fem::IWG_Type::SPATIALDIFF_GHOST ) );
		tParameterList( 3 )( tIWGCounter ).set( "dof_residual", 			  "TEMP" );
		tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies",    "TEMP" );
		tParameterList( 3 )( tIWGCounter ).set( "follower_dof_dependencies",	  "TEMP" );
		tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   "SPGPTemp_Outer,GhostDispl" );
		tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",			  "ghost_p2" );
		tIWGCounter++;

		tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
		tParameterList( 3 )( tIWGCounter ).set( "IWG_name", 				  "IWGGPTemp" );
		tParameterList( 3 )( tIWGCounter ).set( "IWG_type", 				  static_cast< uint >( fem::IWG_Type::SPATIALDIFF_GHOST ) );
		tParameterList( 3 )( tIWGCounter ).set( "dof_residual", 			  "TEMP" );
		tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies",    "TEMP" );
		tParameterList( 3 )( tIWGCounter ).set( "follower_dof_dependencies",	  "TEMP" );
		tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   "SPGPTemp_Inner,GhostDispl" );
		tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",			  "ghost_p3" );
		tIWGCounter++;

        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWGInterfaceTEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::SPATIALDIFF_INTERFACE ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies",    "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "follower_dof_dependencies",     "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMDiffusion_Outer,Diffusion" );
        tParameterList( 3 )( tIWGCounter ).set( "follower_constitutive_models",  "CMDiffusion_Inner,Diffusion" );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",
                "SPInterfaceNitsche     ,NitscheInterface;"
                "SPInterfaceLeaderWeight,LeaderWeightInterface;"
                "SPInterfaceFollowerWeight ,FollowerWeightInterface"   );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             "dbl_iside_p0_2_p1_3" );
        tIWGCounter++;

        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWGTimeContinuityTemp" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::TIME_CONTINUITY_DOF ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies",    "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties",
                "PropWeightCurrent,WeightCurrent;"
                "PropWeightPrevious,WeightPrevious;"
                "PropInitialCondition,InitialCondition" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             "HMR_dummy_n_p2,HMR_dummy_c_p2,HMR_dummy_n_p3,HMR_dummy_c_p3" );
        tParameterList( 3 )( tIWGCounter ).set( "time_continuity",            true );
        tIWGCounter++;

        //------------------------------------------------------------------------------
        // init IQI counter
        uint tIQICounter = 0;

        // create parameter list for IQI 4
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name",                   "IQIBulkTEMP" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::DOF ) );
        tParameterList( 4 )( tIQICounter ).set( "IQI_output_type",            static_cast< uint >( vis::Output_Type::TEMP ) );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies",    "TEMP" );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index",      0 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names",             "HMR_dummy_n_p2,HMR_dummy_c_p2,HMR_dummy_n_p3,HMR_dummy_c_p3" );
        tIQICounter++;
    }

    void SOLParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 7 );
        for( uint Ik = 0; Ik < 7; Ik ++)
        {
            tParameterlist( Ik ).resize( 1 );
        }

        tParameterlist( 0 )( 0 ) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL );

        tParameterlist( 1 )( 0 ) = moris::prm::create_linear_solver_parameter_list();

        tParameterlist( 2 )( 0 ) = moris::prm::create_nonlinear_algorithm_parameter_list();
        tParameterlist( 2 )( 0 ).set("NLA_rel_res_norm_drop",    1e-04 );
        tParameterlist( 2 )( 0 ).set("NLA_relaxation_parameter", 1.00  );
        tParameterlist( 2 )( 0 ).set("NLA_max_iter", 20 );

        tParameterlist( 3 )( 0 ) = moris::prm::create_nonlinear_solver_parameter_list();
        tParameterlist( 3 )( 0 ).set("NLA_DofTypes"      , "TEMP" );

        tParameterlist( 4 )( 0 ) = moris::prm::create_time_solver_algorithm_parameter_list();
        tParameterlist( 4 )( 0 ).set("TSA_Num_Time_Steps",     50 );
        tParameterlist( 4 )( 0 ).set("TSA_Time_Frame",         0.5 );

        tParameterlist( 5 )( 0 ) = moris::prm::create_time_solver_parameter_list();
        tParameterlist( 5 )( 0 ).set("TSA_DofTypes",           "TEMP" );
        tParameterlist( 5 )( 0 ).set("TSA_Initialize_Sol_Vec", "TEMP,0.0" );
        tParameterlist( 5 )( 0 ).set("TSA_Output_Indices",     "0" );
        tParameterlist( 5 )( 0 ).set("TSA_Output_Criteria",     "Output_Criterion" );
        tParameterlist( 5 )( 0 ).set("TSA_time_level_per_type","TEMP,2" );

        tParameterlist( 6 )( 0 ) = moris::prm::create_solver_warehouse_parameterlist();
    }

    void MSIParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_msi_parameter_list();
    }

    void VISParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_vis_parameter_list();
        tParameterlist( 0 )( 0 ).set( "File_Name"  , std::pair< std::string, std::string >( "./", "Vessel_2D.exo" ) );
        tParameterlist( 0 )( 0 ).set( "Mesh_Type"  , static_cast< uint >( vis::VIS_Mesh_Type::STANDARD ) );
        tParameterlist( 0 )( 0 ).set( "Set_Names"  , std::string( "HMR_dummy_n_p2,HMR_dummy_c_p2,HMR_dummy_n_p3,HMR_dummy_c_p3" ) );
        tParameterlist( 0 )( 0 ).set( "Field_Names", std::string( "TEMP" ) );
        tParameterlist( 0 )( 0 ).set( "Field_Type" , std::string( "NODAL" ) );
        tParameterlist( 0 )( 0 ).set( "Output_Type", std::string( "TEMP" ) );
    }

    //------------------------------------------------------------------------------
}

//------------------------------------------------------------------------------
#ifdef  __cplusplus
}
#endif

