/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Comsol.cpp
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
    moris::real tPlaneBottom     = 0.0; /* y bottom plane (m) */
    moris::real tPlaneTop        = 1.0; /* y top plane    (m) */
    moris::real tPlaneLeft       = 0.0; /* x left plane   (m) */
    moris::real tPlaneRight      = 1.0; /* x right plane  (m) */

    moris::sint tStep = 300;
    moris::real tTmax = 57600.0;
    moris::real tDirichletRampUp = 3.0;

    // Constant function for properties
    void Func_Const( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
            moris::fem::Field_Interpolator_Manager         * aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    void Func_Initial_Condition( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
            moris::fem::Field_Interpolator_Manager         * aFIManager )
    {
        aPropMatrix = {{ 313.0 }};
    }

    void Func_Wall_Condition( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
            moris::fem::Field_Interpolator_Manager 	* aFIManager )
    {
        real tT = aFIManager->get_IP_geometry_interpolator()->valt()( 0 );

        real value = std::min(350.0, 313.0+(350-313)*tT/(tDirichletRampUp*tTmax/tStep));

        aPropMatrix      = {{ 0.0 }};
        aPropMatrix(0,0) = value;
    }

    bool Output_Criterion( moris::tsa::Time_Solver * aTimeSolver )
    {
        return true;
    }

    moris::real Func_Bottom_Plane( const moris::Matrix< DDRMat >     & aCoordinates,
            const moris::Cell< moris::real* > & aGeometryParameters )
    {
        moris::real aReturnValue = aCoordinates( 1 ) - 10000; //tPlaneBottom - 0.01;
        return aReturnValue;
    }

    moris::real Func_Top_Plane( const moris::Matrix< DDRMat >        & aCoordinates,
            const moris::Cell< moris::real* > & aGeometryParameters )
    {
        moris::real aReturnValue = aCoordinates( 1 ) - tPlaneTop + 0.01;
        return aReturnValue;
    }

    moris::real Func_Left_Plane( const moris::Matrix< DDRMat >       & aCoordinates,
            const moris::Cell< moris::real* > & aGeometryParameters )
    {
        moris::real aReturnValue = aCoordinates( 0 ) - tPlaneLeft - 0.01;
        return aReturnValue;
    }

    moris::real Func_Right_Plane( const moris::Matrix< DDRMat >      & aCoordinates,
            const moris::Cell< moris::real* > & aGeometryParameters )
    {
        moris::real aReturnValue = aCoordinates( 0 ) - tPlaneRight + 0.01;
        return aReturnValue;
    }

    moris::Matrix< DDRMat > Func_Sensitivity( const moris::Matrix< DDRMat >     & aCoordinates,
            const moris::Cell< moris::real* > & aGeometryParameters )
        {
        moris::Matrix< DDRMat > aReturnValue;
        return aReturnValue;
        }

    void HMRParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_hmr_parameter_list();

        tParameterlist( 0 )( 0 ).set( "number_of_elements_per_dimension", "500,1");
        tParameterlist( 0 )( 0 ).set( "domain_dimensions",                "0.28,0.0005");
        tParameterlist( 0 )( 0 ).set( "domain_offset",                    "0.0,0.0" );
        tParameterlist( 0 )( 0 ).set( "domain_sidesets",                  "1,2,3,4");
        tParameterlist( 0 )( 0 ).set( "lagrange_output_meshes",           "0");

        tParameterlist( 0 )( 0 ).set( "lagrange_orders",  std::string( "1" ));
        tParameterlist( 0 )( 0 ).set( "lagrange_pattern", std::string( "0" ));
        tParameterlist( 0 )( 0 ).set( "bspline_orders",   std::string( "1" ));
        tParameterlist( 0 )( 0 ).set( "bspline_pattern",  std::string( "0" ));

        tParameterlist( 0 )( 0 ).set( "lagrange_to_bspline", "0" );

        tParameterlist( 0 )( 0 ).set( "truncate_bsplines",  1 );
        tParameterlist( 0 )( 0 ).set( "refinement_buffer",  3 );
        tParameterlist( 0 )( 0 ).set( "staircase_buffer",   3 );
        tParameterlist( 0 )( 0 ).set( "initial_refinement", "0" );
        tParameterlist( 0 )( 0 ).set( "initial_refinement_pattern", "0" );

        tParameterlist( 0 )( 0 ).set( "use_multigrid",  0 );
        tParameterlist( 0 )( 0 ).set( "severity_level", 0 );

        tParameterlist( 0 )( 0 ).set( "adaptive_refinement_level", 1 );
    }

    void OPTParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_opt_problem_parameter_list();

        tParameterlist( 0 )( 0 ).set( "is_optimization_problem", false);
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
        tParameterlist( 0 )( 0 ).set( "exodus_output_XTK_ig_mesh", true );
        tParameterlist( 0 )( 0 ).set( "high_to_low_dbl_side_sets", true );
    }

    void GENParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 3 );
        tParameterlist( 0 ).resize( 1 );

        // Main GEN parameter list
        tParameterlist( 0 )( 0 ) = prm::create_gen_parameter_list();
        tParameterlist( 0 )( 0 ).set( "HMR_refinements", 1 );

        // init geometry counter
        uint tGeoCounter = 0;

        // Geometry parameter lists
        tParameterlist( 1 ).push_back( prm::create_user_defined_geometry_parameter_list() );
        tParameterlist( 1 )( tGeoCounter ).set( "field_function_name", "Func_Bottom_Plane");
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

        // create parameter list for property 1
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ) = prm::create_property_parameter_list();
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropDensity" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      "0.75" );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const" );
        tPropCounter++;

        // create parameter list for property 2
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropConductivity" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      "2.1e-7" );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const" );
        tPropCounter++;

        // create parameter list for property 3
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropHeatCapacity" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      "2.4" );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const" );
        tPropCounter++;

        // create parameter list for property 3
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropLatentHeat" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      "175.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const" );
        tPropCounter++;

        // create parameter list for property 9
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropPCTemp" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      "314.5" );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropPhaseState" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      "2.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const" );
        tPropCounter++;

        // create parameter list for property 11
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropPCconst" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      "3.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const" );
        tPropCounter++;

        // create parameter list for property 5
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropImposedFlux" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      "0.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const" );
        tPropCounter++;

        // create parameter list for property 4
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropImposedTemp" );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Wall_Condition" );
        tPropCounter++;

        // create parameter list for property 6
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropWeightCurrent" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      "1000000.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const" );
        tPropCounter++;

        // create parameter list for property 7
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropWeightPrevious" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      "1000000.0" );
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

        // create parameter list for constitutive model 1
        tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_name", "CMDiffusionPhaseChange" );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::DIFF_LIN_ISO_PC ) );
        tParameterList( 1 )( tCMCounter ).set( "dof_dependencies",  std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( 1 )( tCMCounter ).set( "properties",        "PropConductivity,Conductivity;PropDensity,Density;PropHeatCapacity,Heat_Capacity;PropLatentHeat,Latent_Heat;PropPCTemp,PC_Temp;PropPhaseState,Phase_State_Function;PropPCconst,Phase_Change_Const" );
        tCMCounter++;

        //------------------------------------------------------------------------------
        // init SP counter
        uint tSPCounter = 0;

        // create parameter list for stabilization parameter 2
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name",      "SPGPTemp" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type",      static_cast< uint >( fem::Stabilization_Type::GHOST_DISPL ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters",     "0.01" );
        tParameterList( 2 )( tSPCounter ).set( "master_properties",       "PropConductivity,Material" );
        tSPCounter++;

        // create parameter list for stabilization parameter 3
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name",      "SPGGLSDiffusion" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type",      static_cast< uint >( fem::Stabilization_Type::GGLS_DIFFUSION_PC ) );
        tParameterList( 2 )( tSPCounter ).set( "master_properties",       "PropConductivity,Conductivity;PropDensity,Density;PropHeatCapacity,Heat_Capacity" );
        tParameterList( 2 )( tSPCounter ).set( "master_dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tSPCounter++;

        // create parameter list for stabilization parameter 1
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name",      "SPNitscheTemp" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type",      static_cast< uint >( fem::Stabilization_Type::DIRICHLET_NITSCHE ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters",     "1000.0" );
        tParameterList( 2 )( tSPCounter ).set( "master_properties",       "PropConductivity,Material" );
        tSPCounter++;

        //------------------------------------------------------------------------------
        // init IWG counter
        uint tIWGCounter = 0;

        // create parameter list for IWG 1
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWGDiffusionBulk" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::SPATIALDIFF_BULK ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", "CMDiffusionPhaseChange,Diffusion" );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   "SPGGLSDiffusion,GGLS_Param" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             "HMR_dummy_n_p0" );
        tIWGCounter++;

        // create parameter list for IWG 3
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWGInletFlux" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::SPATIALDIFF_NEUMANN ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "master_properties",          "PropImposedFlux,Neumann" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             "SideSet_2_n_p0" );
        tIWGCounter++;

        // create parameter list for IWG 2
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", 		 "IWGOutletTemp" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", 		 static_cast< uint >( fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual",		 "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",   "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "master_properties",	 "PropImposedTemp,Dirichlet" );
        tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models","CMDiffusionPhaseChange,Diffusion" );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",  "SPNitscheTemp,DirichletNitsche" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",		 "SideSet_4_n_p0" );
        tIWGCounter++;

        // create parameter list for IWG 5
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWGTimeContinuityTemp" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::TIME_CONTINUITY_DOF ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "master_properties",          "PropWeightCurrent,WeightCurrent;PropWeightPrevious,WeightPrevious;PropInitialCondition,InitialCondition" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             "HMR_dummy_n_p0" );
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
        tParameterList( 4 )( tIQICounter ).set( "master_dof_dependencies",    "TEMP" );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index",      0 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names",             "HMR_dummy_n_p0" );
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
        tParameterlist( 2 )( 0 ).set("NLA_rel_res_norm_drop", 1e-05 );
        tParameterlist( 2 )( 0 ).set("NLA_relaxation_parameter", 0.9 );
        tParameterlist( 2 )( 0 ).set("NLA_max_iter", 30 );

        tParameterlist( 3 )( 0 ) = moris::prm::create_nonlinear_solver_parameter_list();
        tParameterlist( 3 )( 0 ).set("NLA_DofTypes"      , "TEMP" );

        tParameterlist( 4 )( 0 ) = moris::prm::create_time_solver_algorithm_parameter_list();
        tParameterlist( 4 )( 0 ).set("TSA_Num_Time_Steps",     tStep );
        tParameterlist( 4 )( 0 ).set("TSA_Time_Frame",         tTmax );

        tParameterlist( 5 )( 0 ) = moris::prm::create_time_solver_parameter_list();
        tParameterlist( 5 )( 0 ).set("TSA_DofTypes",            "TEMP" );
        tParameterlist( 5 )( 0 ).set("TSA_Initialize_Sol_Vec",  "TEMP,0.0" );
        tParameterlist( 5 )( 0 ).set("TSA_Output_Indices",      "0" );
        tParameterlist( 5 )( 0 ).set("TSA_Output_Criteria",      "Output_Criterion" );
        tParameterlist( 5 )( 0 ).set("TSA_time_level_per_type", "TEMP,2" );

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
        tParameterlist( 0 )( 0 ).set( "File_Name"  , std::pair< std::string, std::string >( "./", "Comsol.exo" ) );
        tParameterlist( 0 )( 0 ).set( "Mesh_Type"  , static_cast< uint >( vis::VIS_Mesh_Type::STANDARD ) );
        tParameterlist( 0 )( 0 ).set( "Set_Names"  , std::string( "HMR_dummy_n_p0" ) );
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

