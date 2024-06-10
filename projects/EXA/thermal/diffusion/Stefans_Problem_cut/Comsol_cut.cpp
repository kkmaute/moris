/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Comsol_cut.cpp
 *
 */

#include <string>
#include <iostream>
#include "moris_typedefs.hpp"
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

#ifdef __cplusplus
extern "C" {

// global variables
extern uint gInterpolationOrder;

#endif
//------------------------------------------------------------------------------
namespace moris
{

    // Geometry Parameters
    moris::real tXlength = 0.028;
    moris::real tYlength = 0.0015;
    moris::real tXcenter = 0.5 * tXlength;
    moris::real tYcenter = 0.5 * tYlength;
    moris::real tEps     = 1.0e-4;

    // mesh
    // std::string tNumElemsPerDim = "520,6";
    // std::string tHMRDomainDimensions = "0.2912,0.0030";
    std::string tNumElemsPerDim      = "26,3";
    std::string tHMRDomainDimensions = "0.030,0.0030";

    // time solver parameters
    moris::sint tStep = 15;
    moris::real tTmax = 480.0;

    // ramp up of Dirichlet BC (number of time slabs to ramp up the value on the BC)
    moris::real tRampUp = 9.0;

    // Constant function for properties
    void
    Func_Const(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    void
    Func_Initial_Condition(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix = { { 313.0 } };
    }

    void
    Func_Wall_Condition(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        real tT = aFIManager->get_IP_geometry_interpolator()->valt()( 0 );

        real value = std::min( 350.0, 313.0 + ( 350 - 313 ) * tT / ( tRampUp * tTmax / tStep ) );

        aPropMatrix         = { { 0.0 } };
        aPropMatrix( 0, 0 ) = value;
    }

    bool
    Output_Criterion( moris::tsa::Time_Solver* aTimeSolver )
    {
        return true;
    }

    // Level set functions for Rectangle
    moris::real Top_Boundary(
            const moris::Matrix< DDRMat >& aCoordinates,
            const Vector< real >&     aGeometryParameters )
    {
        moris::real tLSval = 0.5 * tYlength - aCoordinates( 1 ) + tYcenter;

        // clean return value to return non-zero value
        return tLSval;
    }

    moris::real Bottom_Boundary(
            const moris::Matrix< DDRMat >& aCoordinates,
            const Vector< real >&     aGeometryParameters )
    {
        moris::real tLSval = 0.5 * tYlength + aCoordinates( 1 ) - tYcenter;

        // clean return value to return non-zero value
        return tLSval;
    }

    moris::real Left_Boundary(
            const moris::Matrix< DDRMat >& aCoordinates,
            const Vector< real >&     aGeometryParameters )
    {
        moris::real tLSval = 0.5 * tXlength + aCoordinates( 0 ) - tXcenter;

        // clean return value to return non-zero value
        return tLSval;
    }

    moris::real Right_Boundary(
            const moris::Matrix< DDRMat >& aCoordinates,
            const Vector< real >&     aGeometryParameters )
    {
        moris::real tLSval = 0.5 * tXlength - aCoordinates( 0 ) + tXcenter;

        // clean return value to return non-zero value
        return tLSval;
    }

    moris::Matrix< DDRMat > Func_Sensitivity(
            const moris::Matrix< DDRMat >& aCoordinates,
            const Vector< real >&     aGeometryParameters )
    {
        moris::Matrix< DDRMat > aReturnValue;
        return aReturnValue;
    }

    void
    HMRParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_hmr_parameter_list();

        tParameterlist( 0 )( 0 ).set( "number_of_elements_per_dimension", tNumElemsPerDim );
        tParameterlist( 0 )( 0 ).set( "domain_dimensions", tHMRDomainDimensions );
        tParameterlist( 0 )( 0 ).set( "domain_offset", "-0.0015,-0.0008" );
        tParameterlist( 0 )( 0 ).set( "domain_sidesets", "1,2,3,4" );
        tParameterlist( 0 )( 0 ).set( "lagrange_output_meshes", "0" );

        switch ( gInterpolationOrder )
        {
            case 1:
            {
                tParameterlist( 0 )( 0 ).set( "lagrange_orders", "1" );
                tParameterlist( 0 )( 0 ).set( "bspline_orders", "1" );
                break;
            }
            case 2:
            {
                tParameterlist( 0 )( 0 ).set( "lagrange_orders", "2" );
                tParameterlist( 0 )( 0 ).set( "bspline_orders", "2" );
                break;
            }
            default:
            {
                MORIS_ERROR( false, "EXA::Comsol_cut: This 2D Example can only be run with Linear or Quadratic" );
            }
        }

        tParameterlist( 0 )( 0 ).set( "lagrange_pattern", "0" );
        tParameterlist( 0 )( 0 ).set( "bspline_pattern", "0" );

        tParameterlist( 0 )( 0 ).set( "lagrange_to_bspline", "0" );

        tParameterlist( 0 )( 0 ).set( "truncate_bsplines", 1 );
        tParameterlist( 0 )( 0 ).set( "refinement_buffer", 3 );
        tParameterlist( 0 )( 0 ).set( "staircase_buffer", 3 );
        tParameterlist( 0 )( 0 ).set( "initial_refinement", "0" );
        tParameterlist( 0 )( 0 ).set( "initial_refinement_pattern", "0" );

        tParameterlist( 0 )( 0 ).set( "use_multigrid", 0 );
        tParameterlist( 0 )( 0 ).set( "severity_level", 0 );

        tParameterlist( 0 )( 0 ).set( "adaptive_refinement_level", 1 );
    }

    void
    OPTParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_opt_problem_parameter_list();

        tParameterlist( 0 )( 0 ).set( "is_optimization_problem", false );
    }

    void
    XTKParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_xtk_parameter_list();
        tParameterlist( 0 )( 0 ).set( "decompose", true );
        tParameterlist( 0 )( 0 ).set( "decomposition_type", "conformal" );
        tParameterlist( 0 )( 0 ).set( "enrich", true );
        tParameterlist( 0 )( 0 ).set( "basis_rank", "bspline" );
        tParameterlist( 0 )( 0 ).set( "enrich_mesh_indices", "0" );
        tParameterlist( 0 )( 0 ).set( "ghost_stab", true );
        tParameterlist( 0 )( 0 ).set( "multigrid", false );
        tParameterlist( 0 )( 0 ).set( "verbose", true );
        tParameterlist( 0 )( 0 ).set( "print_enriched_ig_mesh", true );
        tParameterlist( 0 )( 0 ).set( "exodus_output_XTK_ig_mesh", true );
        tParameterlist( 0 )( 0 ).set( "high_to_low_dbl_side_sets", true );
    }

    void
    GENParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
    {
        tParameterlist.resize( 3 );
        tParameterlist( 0 ).resize( 1 );

        // Main GEN parameter list
        tParameterlist( 0 )( 0 ) = prm::create_gen_parameter_list();

        // init geometry counter
        uint tGeoCounter = 0;

        // Geometry parameter lists
        tParameterlist( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        tParameterlist( 1 )( tGeoCounter ).set( "field_function_name", "Bottom_Boundary" );
        tGeoCounter++;

        tParameterlist( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        tParameterlist( 1 )( tGeoCounter ).set( "field_function_name", "Top_Boundary" );
        tGeoCounter++;

        tParameterlist( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        tParameterlist( 1 )( tGeoCounter ).set( "field_function_name", "Right_Boundary" );
        tGeoCounter++;

        tParameterlist( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        tParameterlist( 1 )( tGeoCounter ).set( "field_function_name", "Left_Boundary" );
    }

    void
    FEMParameterList( Vector< Vector< Parameter_List > >& tParameterList )
    {
        // create a cell of cell of parameter list for fem
        tParameterList.resize( 8 );

        //------------------------------------------------------------------------------
        // init property counter
        uint tPropCounter = 0;

        // create parameter list for property 1
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ) = prm::create_property_parameter_list();
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropDensity" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "0.75" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // create parameter list for property 2
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropConductivity" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "2.1e-7" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // create parameter list for property 3
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropHeatCapacity" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "2.4" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // create parameter list for property 3
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropLatentHeat" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "175" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // create parameter list for property 9
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropPCTemp" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "314.5" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropPhaseState" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "2.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // create parameter list for property 11
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropPCconst" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "3.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // create parameter list for property 4
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropImposedTemp" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Wall_Condition" );
        tPropCounter++;

        // create parameter list for property 6
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropWeightCurrent" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "100.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // create parameter list for property 7
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropWeightPrevious" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "100.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // create parameter list for property 8
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropInitialCondition" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Initial_Condition" );
        tPropCounter++;

        //------------------------------------------------------------------------------

        // init CM counter
        uint tCMCounter = 0;

        // create parameter list for constitutive model 1
        tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_name", "CMDiffusionPhaseChange" );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_type",  fem::Constitutive_Type::DIFF_LIN_ISO_PC ) ;
        tParameterList( 1 )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( 1 )( tCMCounter ).set( "properties",
                "PropConductivity,Conductivity;"
                "PropDensity,Density;"
                "PropHeatCapacity,HeatCapacity;"
                "PropLatentHeat,LatentHeat;"
                "PropPCTemp,PCTemp;"
                "PropPhaseState,PhaseStateFunction;"
                "PropPCconst,PhaseChangeConst" );
        tCMCounter++;

        //------------------------------------------------------------------------------

        // init SP counter
        uint tSPCounter = 0;

        // Ghost
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", "SPGPTemp" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type",  fem::Stabilization_Type::GHOST_DISPL ) ;
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", "0.01" );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", "PropConductivity,Material" );
        tSPCounter++;

        // GGLS parameter
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", "SPGGLSDiffusion" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type",  fem::Stabilization_Type::GGLS_DIFFUSION ) ;
        tParameterList( 2 )( tSPCounter ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties",
                "PropConductivity , Conductivity;"
                "PropDensity      , Density;"
                "PropHeatCapacity , HeatCapacity;"
                "PropLatentHeat   , LatentHeat;"
                "PropPCTemp       , PCTemp;"
                "PropPhaseState   , PhaseStateFunction;"
                "PropPCconst      , PhaseChangeConst" );
        tSPCounter++;

        // Dirichlet SP
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", "SPNitscheTemp" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type",  fem::Stabilization_Type::DIRICHLET_NITSCHE ) ;
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", "1000.0" );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", "PropConductivity,Material" );
        tSPCounter++;

        //------------------------------------------------------------------------------
        // init IWG counter
        uint tIWGCounter = 0;

        // Bulk
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGDiffusionBulk" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_BULK ) ;
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMDiffusionPhaseChange,Diffusion" );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPGGLSDiffusion,GGLSParam" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", "HMR_dummy_n_p15,HMR_dummy_c_p15" );
        tIWGCounter++;

        // Dirichlet BC
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGOutletTemp" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE ) ;
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties", "PropImposedTemp,Dirichlet" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMDiffusionPhaseChange,Diffusion" );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPNitscheTemp,DirichletNitsche" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", "iside_b0_15_b1_14" );
        tIWGCounter++;

        // Time Continuity
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGTimeContinuityTemp" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",  fem::IWG_Type::TIME_CONTINUITY_DOF ) ;
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties",
                "PropWeightCurrent,WeightCurrent;"
                "PropWeightPrevious,WeightPrevious;"
                "PropInitialCondition,InitialCondition" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", "HMR_dummy_n_p15,HMR_dummy_c_p15" );
        tParameterList( 3 )( tIWGCounter ).set( "time_continuity", true );
        tIWGCounter++;

        //------------------------------------------------------------------------------
        // init IQI counter
        uint tIQICounter = 0;

        // create parameter list for IQI 4
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkTEMP" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type",  fem::IQI_Type::DOF ) ;
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "TEMP" );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 0 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", "HMR_dummy_n_p15,HMR_dummy_c_p15" );
        tIQICounter++;

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        tParameterList( 5 ).resize( 1 );
        tParameterList( 5 )( 0 ) = prm::create_computation_parameter_list();
    }

    void
    SOLParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
    {
        tParameterlist.resize( 8 );
        for ( uint Ik = 0; Ik < 8; Ik++ )
        {
            tParameterlist( Ik ).resize( 1 );
        }

        tParameterlist( 0 )( 0 ) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL );

        tParameterlist( 1 )( 0 ) = moris::prm::create_linear_solver_parameter_list();

        tParameterlist( 2 )( 0 ) = moris::prm::create_nonlinear_algorithm_parameter_list();
        tParameterlist( 2 )( 0 ).set( "NLA_rel_res_norm_drop", 2.0e-05 );
        tParameterlist( 2 )( 0 ).set( "NLA_relaxation_parameter", 0.96 );
        tParameterlist( 2 )( 0 ).set( "NLA_max_iter", 20 );

        tParameterlist( 3 )( 0 ) = moris::prm::create_nonlinear_solver_parameter_list();
        tParameterlist( 3 )( 0 ).set( "NLA_DofTypes", "TEMP" );

        tParameterlist( 4 )( 0 ) = moris::prm::create_time_solver_algorithm_parameter_list();
        tParameterlist( 4 )( 0 ).set( "TSA_Num_Time_Steps", tStep );
        tParameterlist( 4 )( 0 ).set( "TSA_Time_Frame", tTmax );

        tParameterlist( 5 )( 0 ) = moris::prm::create_time_solver_parameter_list();
        tParameterlist( 5 )( 0 ).set( "TSA_DofTypes", "TEMP" );
        tParameterlist( 5 )( 0 ).set( "TSA_Initialize_Sol_Vec", "TEMP,0.0" );
        tParameterlist( 5 )( 0 ).set( "TSA_Output_Indices", "0" );
        tParameterlist( 5 )( 0 ).set( "TSA_Output_Criteria", "Output_Criterion" );
        tParameterlist( 5 )( 0 ).set( "TSA_time_level_per_type", "TEMP,2" );

        tParameterlist( 6 )( 0 ) = moris::prm::create_solver_warehouse_parameterlist();

        tParameterlist( 7 )( 0 ) = moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::NONE );
    }

    void
    MSIParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_msi_parameter_list();
    }

    void
    VISParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_vis_parameter_list();
        tParameterlist( 0 )( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", "Comsol_cut.exo" ) );
        tParameterlist( 0 )( 0 ).set( "Mesh_Type",  vis::VIS_Mesh_Type::STANDARD ) ;
        tParameterlist( 0 )( 0 ).set( "Set_Names", "HMR_dummy_n_p15,HMR_dummy_c_p15" );
        tParameterlist( 0 )( 0 ).set( "Field_Names", "TEMP" );
        tParameterlist( 0 )( 0 ).set( "Field_Type", "NODAL" );
        tParameterlist( 0 )( 0 ).set( "IQI_Names", "IQIBulkTEMP" );
    }

    void
    MORISGENERALParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
    {
    }

    //------------------------------------------------------------------------------
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
