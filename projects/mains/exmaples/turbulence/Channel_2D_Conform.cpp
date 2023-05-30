/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Channel_2D_Conform.cpp
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

// Constant function for properties
void Func_Const( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
                 moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
                 moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 );
}

// Inlet velocity function
void Func_Inlet_U( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
                   moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
                   moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix.set_size( 2, 1, 0.0 );
    real tY = aFIManager->get_IP_geometry_interpolator()->valx()( 1 );
    aPropMatrix( 0 ) = 25.0 * ( 2.0 * tY - std::pow( tY, 2 ) ) / 2.0;
}

// Inlet velocity function
void Func_Inlet_V( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
                   moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
                   moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix.set_size( 1, 1, 0.0 );
    real tY = aFIManager->get_IP_geometry_interpolator()->valx()( 1 );
    aPropMatrix( 0 ) = 0.2 * ( 2.0 * tY - std::pow( tY, 2 ) ) / 2.0;
}

// Wall distance function
void Func_Wall_Distance(
    moris::Matrix< moris::DDRMat >                 & aPropMatrix,
    moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
    moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix.set_size( 1, 1, 0.0 );
    real tY = aFIManager->get_IP_geometry_interpolator()->valx()( 1 );

    // dist top
    real tDist = tY;

    // dist top
    real tDistTop = 2.0 - tY;
    if( tDist > tDistTop )
    {
        tDist = tDistTop;
    }

    if( tDist < 1E-6 )
    {
        tDist = 1E-6;
    }

    aPropMatrix( 0 ) = tDist;
}

bool Output_Criterion( moris::tsa::Time_Solver * aTimeSolver )
{
    return true;
}

moris::real Func_Plane(
    const moris::Matrix< DDRMat >     & aCoordinates,
    const moris::Cell< moris::real* > & aGeometryParameters )
{
    moris::real aReturnValue = aCoordinates( 0 ) - 100.0;
    return aReturnValue;
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

    tParameterlist( 0 )( 0 ).set( "number_of_elements_per_dimension", "20,10");
    tParameterlist( 0 )( 0 ).set( "domain_dimensions",                "16.0,2.0");
    tParameterlist( 0 )( 0 ).set( "domain_offset",                    "0.0,0.0");
    tParameterlist( 0 )( 0 ).set( "domain_sidesets",                  "1,2,3,4");
    tParameterlist( 0 )( 0 ).set( "lagrange_output_meshes",           "0");

    tParameterlist( 0 )( 0 ).set( "lagrange_orders",  std::string("1" ));
    tParameterlist( 0 )( 0 ).set( "lagrange_pattern", std::string("0" ));
    tParameterlist( 0 )( 0 ).set( "bspline_orders",   std::string("1" ));
    tParameterlist( 0 )( 0 ).set( "bspline_pattern",  std::string("0" ));

    tParameterlist( 0 )( 0 ).set( "lagrange_to_bspline", "0" );

    tParameterlist( 0 )( 0 ).set( "truncate_bsplines",  1 );
    tParameterlist( 0 )( 0 ).set( "refinement_buffer",  3 );
    tParameterlist( 0 )( 0 ).set( "staircase_buffer",   3 );
    tParameterlist( 0 )( 0 ).set( "initial_refinement", "2" );
    tParameterlist( 0 )( 0 ).set( "initial_refinement_pattern", "0" );

    tParameterlist( 0 )( 0 ).set( "use_multigrid",  0 );
    tParameterlist( 0 )( 0 ).set( "severity_level", 0 );

    tParameterlist( 0 )( 0 ).set( "adaptive_refinement_level", 0 );
}

void XTKParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
{
    tParameterlist.resize( 1 );
    tParameterlist( 0 ).resize( 1 );

    tParameterlist( 0 )( 0 ) = prm::create_xtk_parameter_list();
    tParameterlist( 0 )( 0 ).set( "decompose",           true );
    tParameterlist( 0 )( 0 ).set( "decomposition_type",  "conformal" );
    tParameterlist( 0 )( 0 ).set( "enrich",              true );
    tParameterlist( 0 )( 0 ).set( "basis_rank",          "bspline" );
    tParameterlist( 0 )( 0 ).set( "enrich_mesh_indices", "0" );
    tParameterlist( 0 )( 0 ).set( "ghost_stab",          true );
    tParameterlist( 0 )( 0 ).set( "multigrid",           false );
    tParameterlist( 0 )( 0 ).set( "print_enriched_ig_mesh", true );
    tParameterlist( 0 )( 0 ).set( "exodus_output_XTK_ig_mesh", true );
}

void GENParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
{
    tParameterlist.resize( 3 );
    tParameterlist( 0 ).resize( 1 );

    // Main GEN parameter list
    tParameterlist( 0 )( 0 ) = prm::create_gen_parameter_list();
    tParameterlist( 0 )( 0 ).set( "HMR_refinements", 0 );

     // init geometry counter
    uint tGeoCounter = 0;

    // Geometry parameter lists
    tParameterlist( 1 ).push_back( prm::create_user_defined_geometry_parameter_list() );
    tParameterlist( 1 )( tGeoCounter ).set( "field_function_name", "Func_Plane");
    tParameterlist( 1 )( tGeoCounter ).set( "sensitivity_function_name", "Func_Sensitivity");
    tParameterlist( 1 )( tGeoCounter ).set( "constant_parameters", "");
    tGeoCounter++;
}

void FEMParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterList )
{
    // create a cell of cell of parameter list for fem
    tParameterList.resize( 5 );

    //------------------------------------------------------------------------------
    // fill the property part of the parameter list

    // init property counter
    uint tPropCounter = 0;

    // create viscosity property
    tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
    tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropViscosity" );
    tParameterList( 0 )( tPropCounter ).set( "function_parameters",      "0.0018" );
    tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const" );
    tPropCounter++;

    // create density property
    tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
    tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropDensity" );
    tParameterList( 0 )( tPropCounter ).set( "function_parameters",      "1.0" );
    tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const" );
    tPropCounter++;

    // create wall distance property
    tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
    tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropWallDistance" );
    tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Wall_Distance" );
    tPropCounter++;

    // create inlet velocity property
    tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
    tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropInletU" );
    tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Inlet_U" );
    tPropCounter++;

    // create zero velocity property
    tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
    tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropZeroU" );
    tParameterList( 0 )( tPropCounter ).set( "function_parameters",      "0.0;0.0" );
    tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const" );
    tPropCounter++;

    // create init velocity property
    tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
    tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropInitialConditionU" );
    tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Inlet_U" );
    tPropCounter++;

    // create inlet viscosity property
    tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
    tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropInletV" );
    tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Inlet_V" );
    tPropCounter++;

    // create zero viscosity property
    tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
    tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropZeroV" );
    tParameterList( 0 )( tPropCounter ).set( "function_parameters",      "0.0" );
    tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const" );
    tPropCounter++;

    // create init viscosity property
    tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
    tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropInitialConditionV" );
    tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Inlet_V" );
    tPropCounter++;

    // create weight for time continuity property
    tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
    tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropWeight" );
    tParameterList( 0 )( tPropCounter ).set( "function_parameters",      "100.0" );
    tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const" );
    tPropCounter++;

    //------------------------------------------------------------------------------
    // fill the constitutive model part of the parameter list

    // init CM counter
    uint tCMCounter = 0;

    // create fluid constitutive model
    tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
    tParameterList( 1 )( tCMCounter ).set( "constitutive_name", "CMFluid" );
    tParameterList( 1 )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::FLUID_INCOMPRESSIBLE ) );
    tParameterList( 1 )( tCMCounter ).set( "dof_dependencies",  std::pair< std::string, std::string >( "VX,VY;P", "Velocity,Pressure" ) );
    tParameterList( 1 )( tCMCounter ).set( "properties",        "PropViscosity,Viscosity;PropDensity,Density" );
    tCMCounter++;

    // create turbulence constitutive model
    tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
    tParameterList( 1 )( tCMCounter ).set( "constitutive_name", "CMTurbulence" );
    tParameterList( 1 )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::FLUID_TURBULENCE ) );
    tParameterList( 1 )( tCMCounter ).set( "dof_dependencies",  std::pair< std::string, std::string >( "VX,VY;VISCOSITY", "Velocity,Viscosity" ) );
    tParameterList( 1 )( tCMCounter ).set( "properties",        "PropViscosity,Viscosity;PropDensity,Density" );
    tCMCounter++;

    //------------------------------------------------------------------------------
    // fill the stabilization parameter part of the parameter list

    // init SP counter
    uint tSPCounter = 0;

    // create SUPG stabilization parameter for Navier-Stokes
    tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
    tParameterList( 2 )( tSPCounter ).set( "stabilization_name",      "SPSUPGNS" );
    tParameterList( 2 )( tSPCounter ).set( "stabilization_type",      static_cast< uint >( fem::Stabilization_Type::INCOMPRESSIBLE_FLOW ) );
    tParameterList( 2 )( tSPCounter ).set( "function_parameters",     "36.0" );
    tParameterList( 2 )( tSPCounter ).set( "leader_properties",       "PropViscosity,Viscosity;PropDensity,Density" );
    tParameterList( 2 )( tSPCounter ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY;P", "Velocity,Pressure" ) );
    tSPCounter++;

    // create SUPG stabilization parameter for Spalart-Allmaras
    tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
    tParameterList( 2 )( tSPCounter ).set( "stabilization_name",      "SPSUPGSA" );
    tParameterList( 2 )( tSPCounter ).set( "stabilization_type",      static_cast< uint >( fem::Stabilization_Type::SUPG_SPALART_ALLMARAS_TURBULENCE ) );
    tParameterList( 2 )( tSPCounter ).set( "leader_properties",       "PropViscosity,Viscosity;PropWallDistance,WallDistance" );
    tParameterList( 2 )( tSPCounter ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY;VISCOSITY", "Velocity,Viscosity" ) );
    tSPCounter++;

    // create Nitsche stabilization parameter for velocity
    tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
    tParameterList( 2 )( tSPCounter ).set( "stabilization_name",      "SPNitscheU" );
    tParameterList( 2 )( tSPCounter ).set( "stabilization_type",      static_cast< uint >( fem::Stabilization_Type::VELOCITY_DIRICHLET_NITSCHE ) );
    tParameterList( 2 )( tSPCounter ).set( "function_parameters",     "100.0/1.0" );
    tParameterList( 2 )( tSPCounter ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
    tParameterList( 2 )( tSPCounter ).set( "leader_properties",       "PropViscosity,Viscosity;PropDensity,Density" );
    tSPCounter++;

    // create Nitsche stabilization parameter for viscosity
    tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
    tParameterList( 2 )( tSPCounter ).set( "stabilization_name",      "SPNitscheV" );
    tParameterList( 2 )( tSPCounter ).set( "stabilization_type",      static_cast< uint >( fem::Stabilization_Type::TURBULENCE_DIRICHLET_NITSCHE ) );
    tParameterList( 2 )( tSPCounter ).set( "function_parameters",     "100.0" );
    tParameterList( 2 )( tSPCounter ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VISCOSITY", "Viscosity" ) );
    tParameterList( 2 )( tSPCounter ).set( "leader_properties",       "PropViscosity,Viscosity" );
    tSPCounter++;

    //------------------------------------------------------------------------------
    // fill the IWG part of the parameter list

    // init IWG counter
    uint tIWGCounter = 0;

    // create incompressible NS velocity bulk IWG
    tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
    tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWGVelocityBulk" );
    tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_BULK ) );
    tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               "VX,VY" );
    tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies",    "VX,VY;P;VISCOSITY" );
    tParameterList( 3 )( tIWGCounter ).set( "leader_properties",          "PropDensity,Density" );
    tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid;CMTurbulence,TurbulenceFluid" );
    tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   "SPSUPGNS,IncompressibleFlow" );
    tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             "HMR_dummy_n_p0" );
    tIWGCounter++;

    // create incompressible NS pressure bulk IWG
    tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
    tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWGPressureBulk" );
    tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_BULK ) );
    tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               "P" );
    tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies",    "VX,VY;P;VISCOSITY" );
    tParameterList( 3 )( tIWGCounter ).set( "leader_properties",          "PropDensity,Density" );
    tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
    tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   "SPSUPGNS,IncompressibleFlow" );
    tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             "HMR_dummy_n_p0" );
    tIWGCounter++;

    // create SA turbulence bulk IWG
    tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
    tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWGTurbulenceBulk" );
    tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::SPALART_ALLMARAS_TURBULENCE_BULK ) );
    tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               "VISCOSITY" );
    tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies",    "VX,VY;P;VISCOSITY" );
    tParameterList( 3 )( tIWGCounter ).set( "leader_properties",          "PropViscosity,Viscosity;PropWallDistance,WallDistance" );
    tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   "SPSUPGSA,SUPG" );
    tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             "HMR_dummy_n_p0" );
    tIWGCounter++;

    // create incompressible NS velocity Dirichlet IWG for inlet
    tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
    tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWGInletVelocity" );
    tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_SYMMETRIC_NITSCHE ) );
    tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               "VX,VY" );
    tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies",    "VX,VY;P;VISCOSITY" );
    tParameterList( 3 )( tIWGCounter ).set( "leader_properties",          "PropInletU,Dirichlet" );
    tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid;CMTurbulence,TurbulenceFluid" );
    tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   "SPNitscheU,DirichletNitsche" );
    tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             "SideSet_4_n_p0" );
    tIWGCounter++;

    // create incompressible NS pressure Dirichlet IWG for inlet
    tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
    tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWGInletPressure" );
    tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_SYMMETRIC_NITSCHE ) );
    tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               "P" );
    tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies",    "VX,VY;P;VISCOSITY" );
    tParameterList( 3 )( tIWGCounter ).set( "leader_properties",          "PropInletU,Dirichlet" );
    tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
    tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             "SideSet_4_n_p0" );
    tIWGCounter++;

    // create incompressible NS velocity Dirichlet IWG for walls
    tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
    tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWGZeroVelocity" );
    tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_SYMMETRIC_NITSCHE ) );
    tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               "VX,VY" );
    tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies",    "VX,VY;P;VISCOSITY" );
    tParameterList( 3 )( tIWGCounter ).set( "leader_properties",          "PropZeroU,Dirichlet" );
    tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid;CMTurbulence,TurbulenceFluid" );
    tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   "SPNitscheU,DirichletNitsche" );
    tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             "SideSet_1_n_p0,SideSet_3_n_p0" );
    tIWGCounter++;

    // create incompressible NS pressure Dirichlet IWG for walls
    tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
    tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWGZeroPressure" );
    tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_SYMMETRIC_NITSCHE ) );
    tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               "P" );
    tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies",    "VX,VY;P;VISCOSITY" );
    tParameterList( 3 )( tIWGCounter ).set( "leader_properties",          "PropZeroU,Dirichlet" );
    tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
    tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             "SideSet_1_n_p0,SideSet_3_n_p0" );
    tIWGCounter++;

    // create viscosity Dirichlet IWG for inlet
    tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
    tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWGInletViscosity" );
    tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::SPALART_ALLMARAS_TURBULENCE_DIRICHLET_SYMMETRIC_NITSCHE ) );
    tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               "VISCOSITY" );
    tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies",    "VX,VY;P;VISCOSITY" );
    tParameterList( 3 )( tIWGCounter ).set( "leader_properties",          "PropViscosity,Viscosity;PropInletV,Dirichlet" );
    tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   "SPNitscheV,Nitsche" );
    tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             "SideSet_4_n_p0" );
    tIWGCounter++;

    // create viscosity Dirichlet IWG for walls
    tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
    tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWGZeroViscosity" );
    tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::SPALART_ALLMARAS_TURBULENCE_DIRICHLET_SYMMETRIC_NITSCHE ) );
    tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               "VISCOSITY" );
    tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies",    "VX,VY;P;VISCOSITY" );
    tParameterList( 3 )( tIWGCounter ).set( "leader_properties",          "PropViscosity,Viscosity;PropZeroV,Dirichlet" );
    tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   "SPNitscheV,Nitsche" );
    tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             "SideSet_1_n_p0,SideSet_3_n_p0" );
    tIWGCounter++;

    // create velocity time continuity IWG
    tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
    tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWGVelocityTimeContinuity" );
    tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::TIME_CONTINUITY_DOF ) );
    tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               "VX,VY" );
    tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies",    "VX,VY;P;VISCOSITY" );
    tParameterList( 3 )( tIWGCounter ).set( "leader_properties",          "PropWeight,WeightCurrent;PropWeight,WeightPrevious;PropInitialConditionU,InitialCondition" );
    tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             "HMR_dummy_n_p0" );
    tParameterList( 3 )( tIWGCounter ).set( "time_continuity",            true );
    tIWGCounter++;

    // create viscosity time continuity IWG
    tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
    tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWGViscosityTimeContinuity" );
    tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::TIME_CONTINUITY_DOF ) );
    tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               "VISCOSITY" );
    tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies",    "VX,VY;P;VISCOSITY" );
    tParameterList( 3 )( tIWGCounter ).set( "leader_properties",          "PropWeight,WeightCurrent;PropWeight,WeightPrevious;PropInitialConditionV,InitialCondition" );
    tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             "HMR_dummy_n_p0" );
    tParameterList( 3 )( tIWGCounter ).set( "time_continuity",            true );
    tIWGCounter++;

    //------------------------------------------------------------------------------
    // fill the IQI part of the parameter list

    // init IQI counter
    uint tIQICounter = 0;

    // create parameter list for IQI 1
    tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
    tParameterList( 4 )( tIQICounter ).set( "IQI_name",                   "IQIBulkVX" );
    tParameterList( 4 )( tIQICounter ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::DOF ) );
    tParameterList( 4 )( tIQICounter ).set( "IQI_output_type",            static_cast< uint >( vis::Output_Type::VX ) );
    tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies",    "VX,VY" );
    tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index",      0 );
    tParameterList( 4 )( tIQICounter ).set( "mesh_set_names",             "HMR_dummy_n_p0" );
    tIQICounter++;

    // create parameter list for IQI 2
    tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
    tParameterList( 4 )( tIQICounter ).set( "IQI_name",                   "IQIBulkVY" );
    tParameterList( 4 )( tIQICounter ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::DOF ) );
    tParameterList( 4 )( tIQICounter ).set( "IQI_output_type",            static_cast< uint >( vis::Output_Type::VY ) );
    tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies",    "VX,VY" );
    tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index",      1 );
    tParameterList( 4 )( tIQICounter ).set( "mesh_set_names",             "HMR_dummy_n_p0" );
    tIQICounter++;

    // create parameter list for IQI 3
    tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
    tParameterList( 4 )( tIQICounter ).set( "IQI_name",                   "IQIBulkP" );
    tParameterList( 4 )( tIQICounter ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::DOF ) );
    tParameterList( 4 )( tIQICounter ).set( "IQI_output_type",            static_cast< uint >( vis::Output_Type::P ) );
    tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies",    "P" );
    tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index",      0 );
    tParameterList( 4 )( tIQICounter ).set( "mesh_set_names",             "HMR_dummy_n_p0" );
    tIQICounter++;

    // create parameter list for IQI 3
    tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
    tParameterList( 4 )( tIQICounter ).set( "IQI_name",                   "IQIBulkVISCOSITY" );
    tParameterList( 4 )( tIQICounter ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::DOF ) );
    tParameterList( 4 )( tIQICounter ).set( "IQI_output_type",            static_cast< uint >( vis::Output_Type::VISCOSITY ) );
    tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies",    "VISCOSITY" );
    tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index",      0 );
    tParameterList( 4 )( tIQICounter ).set( "mesh_set_names",             "HMR_dummy_n_p0" );
    tIQICounter++;
}

void SOLParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
{
    tParameterlist.resize( 7 );
    for( uint Ik = 0; Ik < 7; Ik ++ )
    {
        tParameterlist( Ik ).resize( 1 );
    }

    tParameterlist( 0 )( 0 ) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL );

    tParameterlist( 1 )( 0 ) = moris::prm::create_linear_solver_parameter_list();

    tParameterlist( 2 )( 0 ) = moris::prm::create_nonlinear_algorithm_parameter_list();
    tParameterlist( 2 )( 0 ).set("NLA_rel_res_norm_drop", 1e-05 );
    tParameterlist( 2 )( 0 ).set("NLA_relaxation_parameter" , 0.9 );
    tParameterlist( 2 )( 0 ).set("NLA_max_iter", 20 );

    tParameterlist( 3 )( 0 ) = moris::prm::create_nonlinear_solver_parameter_list();
    tParameterlist( 3 )( 0 ).set("NLA_DofTypes"      , "VX,VY;P;VISCOSITY" );

    tParameterlist( 4 )( 0 ) = moris::prm::create_time_solver_algorithm_parameter_list();
    tParameterlist( 4 )( 0 ).set("TSA_Num_Time_Steps",     20 );
    tParameterlist( 4 )( 0 ).set("TSA_Time_Frame",         0.2 );

    tParameterlist( 5 )( 0 ) = moris::prm::create_time_solver_parameter_list();
    tParameterlist( 5 )( 0 ).set("TSA_DofTypes",            "VX,VY;P;VISCOSITY" );
    tParameterlist( 5 )( 0 ).set("TSA_Initialize_Sol_Vec",  "VX,12.5;VY,0.0;P,0.0;VISCOSITY,0.1" );
    tParameterlist( 5 )( 0 ).set("TSA_Output_Indices",      "0" );
    tParameterlist( 5 )( 0 ).set("TSA_Output_Criteria",      "Output_Criterion" );
    tParameterlist( 5 )( 0 ).set("TSA_time_level_per_type", "VX,2;VY,2;P,2;VISCOSITY,2" );

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
    tParameterlist( 0 )( 0 ).set( "File_Name"  , std::pair< std::string, std::string >( "./", "Fluid_Benchmark_2D_Channel_Turbulent_Conform.exo" ) );
    tParameterlist( 0 )( 0 ).set( "Mesh_Type"  , static_cast< uint >( vis::VIS_Mesh_Type::STANDARD ) );
    tParameterlist( 0 )( 0 ).set( "Set_Names"  , std::string( "HMR_dummy_n_p0" ) );
    tParameterlist( 0 )( 0 ).set( "Field_Names", std::string( "VX,VY,P,VISCOSITY" ) );
    tParameterlist( 0 )( 0 ).set( "Field_Type" , std::string( "NODAL,NODAL,NODAL,NODAL" ) );
    tParameterlist( 0 )( 0 ).set( "Output_Type", std::string( "VX,VY,P,VISCOSITY" ) );
}

//------------------------------------------------------------------------------
}

//------------------------------------------------------------------------------
#ifdef  __cplusplus
}
#endif

