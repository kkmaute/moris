/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Two_Channels_with_Separation_Wall_Transient.cpp
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
#endif
//------------------------------------------------------------------------------
namespace moris
{

    // Geometry Parameters
    moris::real tBottomPlane         = 0.0;                   /* Bottom plane y (m) */
    moris::real tTopPlane            = 0.41;                  /* Top plane y (m) */
    moris::real tLeftPlane           = 0.0;                   /* Top plane y (m) */
    moris::real tRightPlane          = 1.845;                 /* Top plane y (m) */
    moris::real tMid                 = 0.205;                 /* Mid y (m) */
    moris::real tWallThickness       = 0.025;                 /* Wall thickness (m) */
    moris::real tWallBottomPlane     = tMid - tWallThickness; /* Mid y (m) */
    moris::real tWallTopPlane        = tMid + tWallThickness; /* Mid y (m) */
    moris::real tRadiusTopChannel    = ( tTopPlane - tWallTopPlane ) / 2;
    moris::real tRadiusBottomChannel = ( tWallBottomPlane - tBottomPlane ) / 2;
    moris::real tYcBottomChannel     = tBottomPlane + tRadiusBottomChannel;
    moris::real tYcTopChannel        = tTopPlane - tRadiusTopChannel;

    // Property constant function
    void
    Func_Const(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    // Inlet velocity function
    void
    Func_Inlet_U(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        // unpack parameters
        real tRadiusChannel = aParameters( 0 )( 0 );
        real tYChannel      = aParameters( 1 )( 0 );
        real tDir           = aParameters( 2 )( 0 );

        // get position in space
        real tY = aFIManager->get_IP_geometry_interpolator()->valx()( 1 );

        // set size for aPropMatrix
        aPropMatrix.set_size( 2, 1, 0.0 );

        // velocity along x direction
        aPropMatrix( 0, 0 ) = -tDir * ( tY - ( tYChannel + tRadiusChannel ) ) * ( tY - ( tYChannel - tRadiusChannel ) ) / ( 2.0 * std::pow( tRadiusChannel, 2.0 ) );
    }

    // Output criterion function
    bool
    Output_Criterion( moris::tsa::Time_Solver* aTimeSolver )
    {
        return true;
    }

        // Plane geometry function
        moris::real Func_Plane(
                const moris::Matrix< DDRMat >     & aCoordinates,
                const Vector< moris::real > & aGeometryParameters )
        {
            moris::real tXNormal = ( aGeometryParameters( 0 ) );
            moris::real tXCenter = ( aGeometryParameters( 1 ) );
            moris::real tYNormal = ( aGeometryParameters( 2 ) );
            moris::real tYCenter = ( aGeometryParameters( 3 ) );

        moris::real aReturnValue =
                tXNormal * ( aCoordinates( 0 ) - tXCenter ) + tYNormal * ( aCoordinates( 1 ) - tYCenter );
        return aReturnValue;
    }

    // OPT parameter list
    void
    OPTParameterList( Vector< Vector< ParameterList > >& aParameterLists )
    {
        aParameterLists( 0 ).push_back( prm::create_opt_problem_parameter_list() );

        aParameterLists( 0 ).set( "is_optimization_problem", false );
    }

    // HMR parameter list
    void
    HMRParameterList( Vector< Vector< ParameterList > >& aParameterLists )
    {
        aParameterLists( 0 ).push_back( prm::create_hmr_parameter_list() );

        aParameterLists( 0 ).set( "number_of_elements_per_dimension", "62,31" );
        aParameterLists( 0 ).set( "domain_dimensions", "3,1" );
        aParameterLists( 0 ).set( "domain_offset", "-0.835,-0.186" );
        aParameterLists( 0 ).set( "domain_sidesets", "1,2,3,4" );
        aParameterLists( 0 ).set( "lagrange_output_meshes", "0" );

        aParameterLists( 0 ).set( "lagrange_orders", "1" );
        aParameterLists( 0 ).set( "lagrange_pattern", "0" );
        aParameterLists( 0 ).set( "bspline_orders", "1" );
        aParameterLists( 0 ).set( "bspline_pattern", "0" );

        aParameterLists( 0 ).set( "lagrange_to_bspline", "0" );

        aParameterLists( 0 ).set( "truncate_bsplines", 1 );
        aParameterLists( 0 ).set( "refinement_buffer", 3 );
        aParameterLists( 0 ).set( "staircase_buffer", 3 );
        aParameterLists( 0 ).set( "initial_refinement", "0" );
        aParameterLists( 0 ).set( "initial_refinement_pattern", "0" );

        aParameterLists( 0 ).set( "use_multigrid", 0 );
        aParameterLists( 0 ).set( "severity_level", 0 );

        aParameterLists( 0 ).set( "adaptive_refinement_level", 0 );
    }

    // XTK parameter list
    void
    XTKParameterList( Vector< Vector< ParameterList > >& aParameterLists )
    {
        aParameterLists( 0 ).push_back( prm::create_xtk_parameter_list() );
        aParameterLists( 0 ).set( "decompose", true );
        aParameterLists( 0 ).set( "decomposition_type", "conformal" );
        aParameterLists( 0 ).set( "enrich", true );
        aParameterLists( 0 ).set( "basis_rank", "bspline" );
        aParameterLists( 0 ).set( "enrich_mesh_indices", "0" );
        aParameterLists( 0 ).set( "ghost_stab", true );
        aParameterLists( 0 ).set( "multigrid", false );
        aParameterLists( 0 ).set( "high_to_low_dbl_side_sets", true );
        aParameterLists( 0 ).set( "print_enriched_ig_mesh", true );
        aParameterLists( 0 ).set( "exodus_output_XTK_ig_mesh", false );
    }

    // GEN parameter list
    void
    GENParameterList( Vector< Vector< ParameterList > >& aParameterLists )
    {
        aParameterLists( 0 ).push_back( prm::create_gen_parameter_list() );

        // Geometry parameter lists
        aParameterLists( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists( 1 ).set( "field_function_name", "Func_Plane" );
        aParameterLists( 1 ).set( "constant_parameters", 0,0,1,0 );

        aParameterLists( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists( 1 ).set( "field_function_name", "Func_Plane" );
        aParameterLists( 1 ).set( "constant_parameters", 0,0,1,0.41 );

        aParameterLists( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists( 1 ).set( "field_function_name", "Func_Plane" );
        aParameterLists( 1 ).set( "constant_parameters", 0,0,1,0.18 );

        aParameterLists( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists( 1 ).set( "field_function_name", "Func_Plane" );
        aParameterLists( 1 ).set( "constant_parameters", 0,0,1,0.23 );

        aParameterLists( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists( 1 ).set( "field_function_name", "Func_Plane" );
        aParameterLists( 1 ).set( "constant_parameters", 1,0,0,0 );

        aParameterLists( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists( 1 ).set( "field_function_name", "Func_Plane" );
        aParameterLists( 1 ).set( "constant_parameters", 1,1.845,0,0 );
    }

    // FEM parameter list
    void
    FEMParameterList( Vector< Vector< ParameterList > >& aParameterLists )
    {
        // create a cell of cell of parameter list for fem

        //------------------------------------------------------------------------------
        // fill the property part of the parameter list

        // fluid viscosity property
        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropFluidViscosity" );
        aParameterLists( 0 ).set( "function_parameters", "1.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // fluid density property
        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropFluidDensity" );
        aParameterLists( 0 ).set( "function_parameters", "1.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // fluid conductivity property
        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropFluidConductivity" );
        aParameterLists( 0 ).set( "function_parameters", "500.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // fluid heat capacity property
        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropFluidHeatCapacity" );
        aParameterLists( 0 ).set( "function_parameters", "0.001" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // solid density property
        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropSolidDensity" );
        aParameterLists( 0 ).set( "function_parameters", "2.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // solid conductivity property
        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropSolidConductivity" );
        aParameterLists( 0 ).set( "function_parameters", "1.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // solid heat capacity property
        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropSolidHeatCapacity" );
        aParameterLists( 0 ).set( "function_parameters", "0.001" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // inlet velocity top channel property
        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropInletVelocityTop" );
        aParameterLists( 0 ).set( "function_parameters", std::string( std::to_string( tRadiusTopChannel ) + "/" std::to_string( tYcTopChannel ) + "/-1.0" ) );
        aParameterLists( 0 ).set( "value_function", "Func_Inlet_U" );

        // inlet velocity bottom channel property
        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropInletVelocityBottom" );
        aParameterLists( 0 ).set( "function_parameters", std::string( std::to_string( tRadiusBottomChannel ) + "/" std::to_string( tYcBottomChannel ) + "/1.0" ) );
        aParameterLists( 0 ).set( "value_function", "Func_Inlet_U" );

        // zero velocity property
        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropZeroVelocity" );
        aParameterLists( 0 ).set( "function_parameters", "0.0;0.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // inlet temperature top channel property
        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropInletTempTop" );
        aParameterLists( 0 ).set( "function_parameters", "1.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // inlet temperature bottom channel property
        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropInletTempBottom" );
        aParameterLists( 0 ).set( "function_parameters", "0.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        //------------------------------------------------------------------------------
        // fill the constitutive model part of the parameter list

        // init CM counter
        moris::uint tCMCounter = 0;

        // fluid constitutive model
        aParameterLists( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        aParameterLists( 1 ).set( "constitutive_name", "CMFluid" );
        aParameterLists( 1 ).set( "constitutive_type",  fem::Constitutive_Type::FLUID_INCOMPRESSIBLE ) ;
        aParameterLists( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( "VX,VY;P", "Velocity,Pressure" ) );
        aParameterLists( 1 ).set( "properties", "PropFluidViscosity,Viscosity;PropFluidDensity,Density" );

        // fluid diffusion constitutive model
        aParameterLists( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        aParameterLists( 1 ).set( "constitutive_name", "CMFluidDiffusion" );
        aParameterLists( 1 ).set( "constitutive_type",  fem::Constitutive_Type::DIFF_LIN_ISO ) ;
        aParameterLists( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists( 1 ).set( "properties", "PropFluidConductivity,Conductivity;PropFluidDensity,Density;PropFluidHeatCapacity,HeatCapacity" );

        // solid diffusion constitutive model
        aParameterLists( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        aParameterLists( 1 ).set( "constitutive_name", "CMSolidDiffusion" );
        aParameterLists( 1 ).set( "constitutive_type",  fem::Constitutive_Type::DIFF_LIN_ISO ) ;
        aParameterLists( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists( 1 ).set( "properties", "PropSolidConductivity,Conductivity;PropSolidDensity,Density;PropSolidHeatCapacity,HeatCapacity" );

        //------------------------------------------------------------------------------
        // fill the stabilization parameter part of the parameter list

        // init SP counter
        moris::uint tSPCounter = 0;

        // NS SUPG stabilization parameter
        aParameterLists( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPSUPGVelocity" );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::INCOMPRESSIBLE_FLOW ) ;
        aParameterLists( 2 ).set( "function_parameters", "36.0" );
        aParameterLists( 2 ).set( "leader_properties", "PropFluidViscosity,Viscosity;PropFluidDensity,Density" );
        aParameterLists( 2 ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY;P", "Velocity,Pressure" ) );

        // advection SUPG stabilization parameter
        aParameterLists( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPSUPGTemp" );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::SUPG_ADVECTION ) ;
        aParameterLists( 2 ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
        aParameterLists( 2 ).set( "leader_properties", "PropFluidConductivity,Conductivity" );

        // nitsche for fluid velocity stabilization parameter
        aParameterLists( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPNitscheVelocity" );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::VELOCITY_DIRICHLET_NITSCHE ) ;
        aParameterLists( 2 ).set( "function_parameters", "100.0/1.0" );
        aParameterLists( 2 ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
        aParameterLists( 2 ).set( "leader_properties", "PropFluidViscosity,Viscosity;PropFluidDensity,Density" );

        // nitsche for fluid temperature stabilization parameter 4
        aParameterLists( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPNitscheFluidTemp" );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::DIRICHLET_NITSCHE ) ;
        aParameterLists( 2 ).set( "function_parameters", "100.0" );
        aParameterLists( 2 ).set( "leader_properties", "PropFluidConductivity,Material" );

        // nitsche thermal interface stabilization parameter
        aParameterLists( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPInterfaceNitsche" );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::NITSCHE_INTERFACE ) ;
        aParameterLists( 2 ).set( "function_parameters", "100.0" );
        aParameterLists( 2 ).set( "leader_properties", "PropFluidConductivity,Material" );
        aParameterLists( 2 ).set( "follower_properties", "PropSolidConductivity,Material" );

        // ghost viscous stabilization parameter
        aParameterLists( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPGPViscous" );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::VISCOUS_GHOST ) ;
        aParameterLists( 2 ).set( "function_parameters", "0.05" );
        aParameterLists( 2 ).set( "leader_properties", "PropFluidViscosity,Viscosity" );

        // ghost convective stabilization parameter
        aParameterLists( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPGPVelocity" );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::CONVECTIVE_GHOST ) ;
        aParameterLists( 2 ).set( "function_parameters", "0.05" );
        aParameterLists( 2 ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
        aParameterLists( 2 ).set( "leader_properties", "PropFluidDensity,Density" );

        // ghost fluid pressure stabilization parameter
        aParameterLists( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPGPPressure" );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::PRESSURE_GHOST ) ;
        aParameterLists( 2 ).set( "function_parameters", "0.05/1.0" );
        aParameterLists( 2 ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
        aParameterLists( 2 ).set( "leader_properties", "PropFluidViscosity,Viscosity;PropFluidDensity,Density" );

        // ghost fluid temperature stabilization parameter
        aParameterLists( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPGPFluidTemp" );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::GHOST_DISPL ) ;
        aParameterLists( 2 ).set( "function_parameters", "0.05" );
        aParameterLists( 2 ).set( "leader_properties", "PropFluidConductivity,Material" );

        // ghost solid temperature stabilization parameter
        aParameterLists( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPGPSolidTemp" );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::GHOST_DISPL ) ;
        aParameterLists( 2 ).set( "function_parameters", "0.05" );
        aParameterLists( 2 ).set( "leader_properties", "PropSolidConductivity,Material" );

        //------------------------------------------------------------------------------
        // fill the IWG part of the parameter list

        // init IWG counter
        moris::uint tIWGCounter = 0;

        // NS incompressible velocity bulk IWG
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGVelocityBulk" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_BULK ) ;
        aParameterLists( 3 ).set( "dof_residual", "VX,VY" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPSUPGVelocity,IncompressibleFlow" );
        aParameterLists( 3 ).set( "mesh_set_names", "HMR_dummy_n_p34,HMR_dummy_c_p34,HMR_dummy_n_p46,HMR_dummy_c_p46" );

        // NS incompressible pressure bulk IWG
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGPressureBulk" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_BULK ) ;
        aParameterLists( 3 ).set( "dof_residual", "P" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPSUPGVelocity,IncompressibleFlow" );
        aParameterLists( 3 ).set( "mesh_set_names", "HMR_dummy_n_p34,HMR_dummy_c_p34,HMR_dummy_n_p46,HMR_dummy_c_p46" );

        // fluid diffusion bulk IWG
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGFluidDiffusionBulk" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_BULK ) ;
        aParameterLists( 3 ).set( "dof_residual", "TEMP" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMFluidDiffusion,Diffusion" );
        aParameterLists( 3 ).set( "mesh_set_names", "HMR_dummy_n_p34,HMR_dummy_c_p34,HMR_dummy_n_p46,HMR_dummy_c_p46" );

        // solid diffusion bulk IWG
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGSolidDiffusionBulk" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_BULK ) ;
        aParameterLists( 3 ).set( "dof_residual", "TEMP" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMSolidDiffusion,Diffusion" );
        aParameterLists( 3 ).set( "mesh_set_names", "HMR_dummy_n_p42,HMR_dummy_c_p42" );

        // fluid advection bulk IWG
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGAdvectionBulk" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::ADVECTION_BULK ) ;
        aParameterLists( 3 ).set( "dof_residual", "TEMP" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMFluidDiffusion,Diffusion" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPSUPGTemp,SUPG" );
        aParameterLists( 3 ).set( "mesh_set_names", "HMR_dummy_n_p34,HMR_dummy_c_p34,HMR_dummy_n_p46,HMR_dummy_c_p46" );

        // top inlet fluid velocity IWG (velocity part)
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGInVelocityTop" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        aParameterLists( 3 ).set( "dof_residual", "VX,VY" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists( 3 ).set( "leader_properties", "PropInletVelocityTop,Dirichlet" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPNitscheVelocity,DirichletNitsche" );
        aParameterLists( 3 ).set( "mesh_set_names", "iside_b0_46_b1_47" );

        // top inlet fluid velocity IWG (pressure part)
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGInPressureTop" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        aParameterLists( 3 ).set( "dof_residual", "P" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists( 3 ).set( "leader_properties", "PropInletVelocityTop,Dirichlet" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists( 3 ).set( "mesh_set_names", "iside_b0_46_b1_47" );

        // bottom inlet fluid velocity IWG (velocity part)
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGInVelocityBottom" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        aParameterLists( 3 ).set( "dof_residual", "VX,VY" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists( 3 ).set( "leader_properties", "PropInletVelocityBottom,Dirichlet" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPNitscheVelocity,DirichletNitsche" );
        aParameterLists( 3 ).set( "mesh_set_names", "iside_b0_34_b1_32" );

        // bottom inlet fluid velocity IWG (pressure part)
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGInPressureBottom" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        aParameterLists( 3 ).set( "dof_residual", "P" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists( 3 ).set( "leader_properties", "PropInletVelocityBottom,Dirichlet" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists( 3 ).set( "mesh_set_names", "iside_b0_34_b1_32" );

        // wall fluid velocity IWG (velocity part)
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGWallVelocity" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        aParameterLists( 3 ).set( "dof_residual", "VX,VY" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists( 3 ).set( "leader_properties", "PropZeroVelocity,Dirichlet" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPNitscheVelocity,DirichletNitsche" );
        aParameterLists( 3 ).set( "mesh_set_names", "iside_b0_34_b1_2,iside_b0_34_b1_42,iside_b0_46_b1_42,iside_b0_46_b1_62" );

        // wall fluid velocity IWG (pressure part)
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGWallPressure" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        aParameterLists( 3 ).set( "dof_residual", "P" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists( 3 ).set( "leader_properties", "PropZeroVelocity,Dirichlet" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists( 3 ).set( "mesh_set_names", "iside_b0_34_b1_2,iside_b0_34_b1_42,iside_b0_46_b1_42,iside_b0_46_b1_62" );

        // top inlet fluid temperature IWG
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGInletTempTop" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        aParameterLists( 3 ).set( "dof_residual", "TEMP" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists( 3 ).set( "leader_properties", "PropInletTempTop,Dirichlet" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMFluidDiffusion,Diffusion" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPNitscheFluidTemp,DirichletNitsche" );
        aParameterLists( 3 ).set( "mesh_set_names", "iside_b0_46_b1_47" );

        // bottom inlet fluid temperature IWG
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGInletTempBottom" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        aParameterLists( 3 ).set( "dof_residual", "TEMP" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists( 3 ).set( "leader_properties", "PropInletTempBottom,Dirichlet" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMFluidDiffusion,Diffusion" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPNitscheFluidTemp,DirichletNitsche" );
        aParameterLists( 3 ).set( "mesh_set_names", "iside_b0_34_b1_32" );

        // fluid/solid thermal interface IWG
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGInterfaceFluidSolid" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE ) ;
        aParameterLists( 3 ).set( "dof_residual", "TEMP" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists( 3 ).set( "follower_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMFluidDiffusion,Diffusion" );
        aParameterLists( 3 ).set( "follower_constitutive_models", "CMSolidDiffusion,Diffusion" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPInterfaceNitsche,NitscheInterface" );
        aParameterLists( 3 ).set( "mesh_set_names", "dbl_iside_p0_34_p1_42,dbl_iside_p0_46_p1_42" );

        // ghost viscous IWG
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGGPViscous" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
        aParameterLists( 3 ).set( "dof_residual", "VX,VY" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists( 3 ).set( "follower_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPGPViscous,GhostSP" );
        aParameterLists( 3 ).set( "mesh_set_names", "ghost_p34,ghost_p46" );

        // ghost convective IWG
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGGPConvective" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
        aParameterLists( 3 ).set( "dof_residual", "VX,VY" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists( 3 ).set( "follower_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPGPVelocity,GhostSP" );
        aParameterLists( 3 ).set( "mesh_set_names", "ghost_p34,ghost_p46" );

        // ghost pressure IWG
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGGPPressure" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
        aParameterLists( 3 ).set( "dof_residual", "P" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists( 3 ).set( "follower_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPGPPressure,GhostSP" );
        aParameterLists( 3 ).set( "mesh_set_names", "ghost_p34,ghost_p46" );

        // ghost fluid temperature IWG
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGGPFluidTemp" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
        aParameterLists( 3 ).set( "dof_residual", "TEMP" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists( 3 ).set( "follower_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPGPFluidTemp,GhostSP" );
        aParameterLists( 3 ).set( "mesh_set_names", "ghost_p34,ghost_p46" );

        // ghost solid temperature IWG
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGGPSolidTemp" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
        aParameterLists( 3 ).set( "dof_residual", "TEMP" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists( 3 ).set( "follower_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPGPSolidTemp,GhostSP" );
        aParameterLists( 3 ).set( "mesh_set_names", "ghost_p42" );

        //------------------------------------------------------------------------------
        // fill the IQI part of the parameter list

        // init IQI counter
        moris::uint tIQICounter = 0;

        // VX IQI
        aParameterLists( 4 ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkVX" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists( 4 ).set( "dof_quantity", "VX,VY" );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "VX,VY" );
        aParameterLists( 4 ).set( "vectorial_field_index", 0 );
        aParameterLists( 4 ).set( "mesh_set_names", "HMR_dummy_n_p34,HMR_dummy_c_p34,HMR_dummy_n_p46,HMR_dummy_c_p46,HMR_dummy_n_p42,HMR_dummy_c_p42" );

        // VY IQI
        aParameterLists( 4 ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkVY" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists( 4 ).set( "dof_quantity", "VX,VY" );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "VX,VY" );
        aParameterLists( 4 ).set( "vectorial_field_index", 1 );
        aParameterLists( 4 ).set( "mesh_set_names", "HMR_dummy_n_p34,HMR_dummy_c_p34,HMR_dummy_n_p46,HMR_dummy_c_p46,HMR_dummy_n_p42,HMR_dummy_c_p42" );

        // P IQI
        aParameterLists( 4 ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkP" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists( 4 ).set( "dof_quantity", "P" );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "P" );
        aParameterLists( 4 ).set( "vectorial_field_index", 0 );
        aParameterLists( 4 ).set( "mesh_set_names", "HMR_dummy_n_p34,HMR_dummy_c_p34,HMR_dummy_n_p46,HMR_dummy_c_p46,HMR_dummy_n_p42,HMR_dummy_c_p42" );

        // TEMP IQI
        aParameterLists( 4 ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkTEMP" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists( 4 ).set( "dof_quantity", "TEMP" );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 4 ).set( "vectorial_field_index", 0 );
        aParameterLists( 4 ).set( "mesh_set_names", "HMR_dummy_n_p34,HMR_dummy_c_p34,HMR_dummy_n_p46,HMR_dummy_c_p46,HMR_dummy_n_p42,HMR_dummy_c_p42" );

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        aParameterLists( 5 ).push_back( prm::create_computation_parameter_list() );
    }

    // SOL parameter list
    void
    SOLParameterList( Vector< Vector< ParameterList > >& aParameterLists )
    {

        aParameterLists( 0 ).push_back( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL ) );

        aParameterLists( 1 ).push_back( moris::prm::create_linear_solver_parameter_list() );

        aParameterLists( 2 ).push_back( moris::prm::create_nonlinear_algorithm_parameter_list() );
        aParameterLists( 2 ).set( "NLA_rel_res_norm_drop", 1e-07 );

        aParameterLists( 3 ).push_back( moris::prm::create_nonlinear_solver_parameter_list() );
        aParameterLists( 3 ).set( "NLA_DofTypes", "VX,VY;P;TEMP" );

        aParameterLists( 4 ).push_back( moris::prm::create_time_solver_algorithm_parameter_list() );
        aParameterLists( 4 ).set( "TSA_Num_Time_Steps", 20 );
        aParameterLists( 4 ).set( "TSA_Time_Frame", 0.4 );

        aParameterLists( 5 ).push_back( moris::prm::create_time_solver_parameter_list() );
        aParameterLists( 5 ).set( "TSA_DofTypes", "VX,VY;P;TEMP" );
        aParameterLists( 5 ).set( "TSA_Initialize_Sol_Vec", "VX,1E-4;VY,1E-4;P,0.0;TEMP,0.0" );
        aParameterLists( 5 ).set( "TSA_Output_Indices", "0" );
        aParameterLists( 5 ).set( "TSA_Output_Criteria", "Output_Criterion" );
        aParameterLists( 5 ).set( "TSA_time_level_per_type", "VX,2;VY,2;P,2;TEMP,2" );

        aParameterLists( 6 ).push_back( moris::prm::create_solver_warehouse_parameterlist() );

        aParameterLists( 7 ).push_back( moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::NONE ) );
    }

    // MSI parametr list
    void
    MSIParameterList( Vector< Vector< ParameterList > >& aParameterLists )
    {
        aParameterLists( 0 ).push_back( prm::create_msi_parameter_list() );
    }

    // VIS parameter list
    void
    VISParameterList( Vector< Vector< ParameterList > >& aParameterLists )
    {
        aParameterLists( 0 ).push_back( prm::create_vis_parameter_list() );
        aParameterLists( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", "2_Channels_Temp.exo" ) );
        aParameterLists( 0 ).set( "Mesh_Type",  vis::VIS_Mesh_Type::STANDARD ) ;
        aParameterLists( 0 ).set( "Set_Names", "HMR_dummy_n_p34,HMR_dummy_c_p34,HMR_dummy_n_p46,HMR_dummy_c_p46,HMR_dummy_n_p42,HMR_dummy_c_p42" );
        aParameterLists( 0 ).set( "Field_Names", "VX,VY,P,TEMP" );
        aParameterLists( 0 ).set( "Field_Type", "NODAL,NODAL,NODAL,NODAL" );
        aParameterLists( 0 ).set( "IQI_Names", "IQIBulkVX,IQIBulkVY,IQIBulkP,IQIBulkTEMP" );
    }

    void
    MORISGENERALParameterList( Vector< Vector< ParameterList > >& aParameterLists )
    {
    }
    //------------------------------------------------------------------------------
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
