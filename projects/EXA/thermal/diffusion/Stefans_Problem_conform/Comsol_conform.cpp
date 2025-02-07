/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Comsol_conform.cpp
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
#include "parameters.hpp"
#include "fn_equal_to.hpp"

#include "AztecOO.h"

#ifdef __cplusplus
extern "C" {

// global variables
extern uint gTestCaseIndex;
extern uint gInterpolationOrder;
extern bool gPrintReferenceValues;

#endif
//------------------------------------------------------------------------------
namespace moris
{

    // Geometry Parameters
    // moris::real tXlength = 0.28;
    moris::real tXlength = 0.028;
    moris::real tYlength = 0.0015;
    moris::real tXcenter = 0.5 * tXlength;
    moris::real tYcenter = 0.5 * tYlength;
    moris::real tEps     = 1.0e-4;

    // mesh
    Vector< uint > tNumElemsPerDim      = { 50, 1 };
    Vector< real > tHMRDomainDimensions = { 0.028, 0.0005 };

    // time solver parameters
    moris::sint tStep = 10;
    moris::real tTmax = 480.0;

    // ramp up of Dirichlet BC (number of time slabs to ramp up the value on the BC)
    moris::real tRampUp = 6.0;

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

    // Dummy level set function for conform mesh
    moris::real Dummy_LS(
            const moris::Matrix< DDRMat >     & aCoordinates,
            const Vector< real > & aGeometryParameters )
    {
        moris::real tLSval = aCoordinates( 0 ) + 1.0;

        // clean return value to return non-zero value
        return tLSval;
    }

    void
    HMRParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "number_of_elements_per_dimension", tNumElemsPerDim );
        aParameterLists.set( "domain_dimensions", tHMRDomainDimensions );
        aParameterLists.set( "lagrange_output_meshes", "0" );

        switch ( gInterpolationOrder )
        {
            case 1:
            {
                aParameterLists.set( "lagrange_orders", "1" );
                aParameterLists.set( "bspline_orders", "1" );
                break;
            }
            case 2:
            {
                aParameterLists.set( "lagrange_orders", "2" );
                aParameterLists.set( "bspline_orders", "2" );
                break;
            }
            default:
            {
                MORIS_ERROR( false, "EXA::Comsol_conform: This 2D Example can only be run with Linear or Quadratic" );
            }
        }

        aParameterLists.set( "lagrange_pattern", "0" );
        aParameterLists.set( "bspline_pattern", "0" );

        aParameterLists.set( "lagrange_to_bspline", "0" );

        aParameterLists.set( "refinement_buffer", 3 );
        aParameterLists.set( "staircase_buffer", 3 );
        aParameterLists.set( "initial_refinement", "0" );
        aParameterLists.set( "initial_refinement_pattern", "0" );

        aParameterLists.set( "use_multigrid", 0 );
        aParameterLists.set( "severity_level", 0 );

        aParameterLists.set( "adaptive_refinement_level", 1 );
    }

    void
    OPTParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "is_optimization_problem", false );
    }

    void
    XTKParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "decompose", true );
        aParameterLists.set( "decomposition_type", "conformal" );
        aParameterLists.set( "enrich", true );
        aParameterLists.set( "basis_rank", "bspline" );
        aParameterLists.set( "enrich_mesh_indices", "0" );
        aParameterLists.set( "ghost_stab", true );
        aParameterLists.set( "multigrid", false );
        aParameterLists.set( "verbose", true );
        aParameterLists.set( "print_enriched_ig_mesh", true );
        aParameterLists.set( "exodus_output_XTK_ig_mesh", true );
        aParameterLists.set( "high_to_low_dbl_side_sets", true );
    }

    void
    GENParameterList( Module_Parameter_Lists& aParameterLists )
    {
        // Geometry parameter lists
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists.set( "field_function_name", "Dummy_LS" );
    }

    void
    FEMParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.hack_for_legacy_fem();
        // create a cell of cell of parameter list for fem

        //------------------------------------------------------------------------------

        // create parameter list for property 1
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDensity" );
        aParameterLists.set( "function_parameters", "0.75" );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 2
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropConductivity" );
        aParameterLists.set( "function_parameters", "2.1e-7" );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 3
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropHeatCapacity" );
        aParameterLists.set( "function_parameters", "2.4" );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 3
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropLatentHeat" );
        aParameterLists.set( "function_parameters", "175" );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 9
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropPCTemp" );
        aParameterLists.set( "function_parameters", "314.5" );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropPhaseState" );
        aParameterLists.set( "function_parameters", "2.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 11
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropPCconst" );
        aParameterLists.set( "function_parameters", "3.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 4
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropImposedTemp" );
        aParameterLists.set( "value_function", "Func_Wall_Condition" );

        // create parameter list for property 6
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropWeightCurrent" );
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 7
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropWeightPrevious" );
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 8
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropInitialCondition" );
        aParameterLists.set( "value_function", "Func_Initial_Condition" );

        //------------------------------------------------------------------------------

        // create parameter list for constitutive model 1
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMDiffusionPhaseChange" );
        aParameterLists.set( "constitutive_type",  fem::Constitutive_Type::DIFF_LIN_ISO_PC ) ;
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists.set( "properties",
                "PropConductivity,Conductivity;"
                "PropDensity,Density;"
                "PropHeatCapacity,HeatCapacity;"
                "PropLatentHeat,LatentHeat;"
                "PropPCTemp,PCTemp;"
                "PropPhaseState,PhaseStateFunction;"
                "PropPCconst,PhaseChangeConst" );

        //------------------------------------------------------------------------------

        // Ghost
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPGPTemp" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::GHOST_DISPL ) ;
        aParameterLists.set( "function_parameters", "0.01" );
        aParameterLists.set( "leader_properties", "PropConductivity,Material" );

        // GGLS parameter
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPGGLSDiffusion" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::GGLS_DIFFUSION ) ;
        aParameterLists.set( "leader_dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists.set( "leader_properties",
                "PropConductivity , Conductivity;"
                "PropDensity      , Density;"
                "PropHeatCapacity , HeatCapacity;"
                "PropLatentHeat   , LatentHeat;"
                "PropPCTemp       , PCTemp;"
                "PropPhaseState   , PhaseStateFunction;"
                "PropPCconst      , PhaseChangeConst" );

        // Dirichlet SP
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPNitscheTemp" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::DIRICHLET_NITSCHE ) ;
        aParameterLists.set( "function_parameters", "1000.0" );
        aParameterLists.set( "leader_properties", "PropConductivity,Material" );

        //------------------------------------------------------------------------------
        // Bulk
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGDiffusionBulk" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_BULK ) ;
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusionPhaseChange,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPGGLSDiffusion,GGLSParam" );
        aParameterLists.set( "mesh_set_names", "HMR_dummy_n_p1" );

        // Dirichlet BC
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGOutletTemp" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE ) ;
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "leader_properties", "PropImposedTemp,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusionPhaseChange,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheTemp,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", "SideSet_4_n_p1" );

        // Time Continuity
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGTimeContinuityTemp" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::TIME_CONTINUITY_DOF ) ;
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "leader_properties",
                "PropWeightCurrent,WeightCurrent;"
                "PropWeightPrevious,WeightPrevious;"
                "PropInitialCondition,InitialCondition" );
        aParameterLists.set( "mesh_set_names", "HMR_dummy_n_p1" );
        aParameterLists.set( "time_continuity", true );

        //------------------------------------------------------------------------------
        // Temperature IQI
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkTEMP" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists.set( "dof_quantity", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "mesh_set_names", "HMR_dummy_n_p1" );

        // Max Temperature IQI
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIMaxTEMP" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::MAX_DOF ) ;
        aParameterLists.set( "dof_quantity", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "function_parameters", "313.0/30.0" );
        aParameterLists.set( "mesh_set_names", "HMR_dummy_n_p1" );

        // Latent Heat IQI
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQILatentHeatAbsorption" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::LATENT_HEAT_ABSORPTION ) ;
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "leader_properties",
                "PropDensity,Density;"
                "PropLatentHeat,LatentHeat;"
                "PropPCTemp,PCTemp;"
                "PropPhaseState,PhaseStateFunction;"
                "PropPCconst,PhaseChangeConst" );
        aParameterLists.set( "mesh_set_names", "HMR_dummy_n_p1" );

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        aParameterLists( FEM::COMPUTATION );
    }

    void
    SOLParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::AMESOS_IMPL );

        aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list();

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "NLA_rel_res_norm_drop", 1e-7 );
        aParameterLists.set( "NLA_relaxation_strategy",  sol::SolverRelaxationType::InvResNormAdaptive ) ;
        aParameterLists.set( "NLA_relaxation_parameter", 0.5 );
        aParameterLists.set( "NLA_relaxation_damping", 0.5 );
        aParameterLists.set( "NLA_max_iter", 40 );

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "NLA_DofTypes", "TEMP" );

        aParameterLists( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "TSA_Num_Time_Steps", tStep );
        aParameterLists.set( "TSA_Time_Frame", tTmax );

        aParameterLists( SOL::TIME_SOLVERS ).add_parameter_list();
        aParameterLists.set( "TSA_DofTypes", "TEMP" );
        aParameterLists.set( "TSA_Initialize_Sol_Vec", "TEMP,0.0" );
        aParameterLists.set( "TSA_Output_Indices", "0" );
        aParameterLists.set( "TSA_Output_Criteria", "Output_Criterion" );
        aParameterLists.set( "TSA_time_level_per_type", "TEMP,2" );

        aParameterLists( SOL::PRECONDITIONERS ).add_parameter_list(  sol::PreconditionerType::NONE );
    }

    void
    MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
    }

    void
    VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "File_Name", std::pair< std::string, std::string >( "./", "Comsol_conform_" + std::to_string( gTestCaseIndex ) + ".exo" ) );
        aParameterLists.set( "Mesh_Type",  vis::VIS_Mesh_Type::STANDARD ) ;
        aParameterLists.set( "Set_Names", "HMR_dummy_n_p1" );
        aParameterLists.set( "Field_Names", "TEMP,MAX_DOF,LATENT_HEAT_ABSORPTION" );
        aParameterLists.set( "Field_Type", "NODAL,GLOBAL,GLOBAL" );
        aParameterLists.set( "IQI_Names", "IQIBulkTEMP,IQIMaxTEMP,IQILatentHeatAbsorption" );
    }

    void
    MORISGENERALParameterList( Module_Parameter_Lists& aParameterLists )
    {
    }
    //------------------------------------------------------------------------------
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
