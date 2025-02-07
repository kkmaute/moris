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

        aParameterLists.set( "is_optimization_problem", false );
    }

    // HMR parameter list
    void
    HMRParameterList( Vector< Vector< ParameterList > >& aParameterLists )
    {
        aParameterLists( 0 ).push_back( prm::create_hmr_parameter_list() );

        aParameterLists.set( "number_of_elements_per_dimension", 62, 31 );
        aParameterLists.set( "domain_dimensions", 3.0, 1.0 );
        aParameterLists.set( "domain_offset", -0.835, -0.186 );
        aParameterLists.set( "lagrange_output_meshes", "0" );

        aParameterLists.set( "lagrange_orders", "1" );
        aParameterLists.set( "lagrange_pattern", "0" );
        aParameterLists.set( "bspline_orders", "1" );
        aParameterLists.set( "bspline_pattern", "0" );

        aParameterLists.set( "lagrange_to_bspline", "0" );

        aParameterLists.set( "refinement_buffer", 3 );
        aParameterLists.set( "staircase_buffer", 3 );
        aParameterLists.set( "initial_refinement", "0" );
        aParameterLists.set( "initial_refinement_pattern", "0" );

        aParameterLists.set( "severity_level", 0 );

        aParameterLists.set( "adaptive_refinement_level", 0 );
    }

    // XTK parameter list
    void
    XTKParameterList( Vector< Vector< ParameterList > >& aParameterLists )
    {
        aParameterLists( 0 ).push_back( prm::create_xtk_parameter_list() );
        aParameterLists.set( "decompose", true );
        aParameterLists.set( "decomposition_type", "conformal" );
        aParameterLists.set( "enrich", true );
        aParameterLists.set( "basis_rank", "bspline" );
        aParameterLists.set( "enrich_mesh_indices", "0" );
        aParameterLists.set( "ghost_stab", true );
        aParameterLists.set( "multigrid", false );
        aParameterLists.set( "high_to_low_dbl_side_sets", true );
        aParameterLists.set( "print_enriched_ig_mesh", true );
        aParameterLists.set( "exodus_output_XTK_ig_mesh", false );
    }

    // GEN parameter list
    void
    GENParameterList( Vector< Vector< ParameterList > >& aParameterLists )
    {
        aParameterLists( 0 ).push_back( prm::create_gen_parameter_list() );

        // Geometry parameter lists
        aParameterLists( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists.set( "field_function_name", "Func_Plane" );
        aParameterLists.set( "constant_parameters", 0,0,1,0 );

        aParameterLists( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists.set( "field_function_name", "Func_Plane" );
        aParameterLists.set( "constant_parameters", 0,0,1,0.41 );

        aParameterLists( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists.set( "field_function_name", "Func_Plane" );
        aParameterLists.set( "constant_parameters", 0,0,1,0.18 );

        aParameterLists( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists.set( "field_function_name", "Func_Plane" );
        aParameterLists.set( "constant_parameters", 0,0,1,0.23 );

        aParameterLists( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists.set( "field_function_name", "Func_Plane" );
        aParameterLists.set( "constant_parameters", 1,0,0,0 );

        aParameterLists( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists.set( "field_function_name", "Func_Plane" );
        aParameterLists.set( "constant_parameters", 1,1.845,0,0 );
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
        aParameterLists.set( "property_name", "PropFluidViscosity" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // fluid density property
        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropFluidDensity" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // fluid conductivity property
        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropFluidConductivity" );
        aParameterLists.set( "function_parameters", "500.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // fluid heat capacity property
        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropFluidHeatCapacity" );
        aParameterLists.set( "function_parameters", "0.001" );
        aParameterLists.set( "value_function", "Func_Const" );

        // solid density property
        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropSolidDensity" );
        aParameterLists.set( "function_parameters", "2.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // solid conductivity property
        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropSolidConductivity" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // solid heat capacity property
        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropSolidHeatCapacity" );
        aParameterLists.set( "function_parameters", "0.001" );
        aParameterLists.set( "value_function", "Func_Const" );

        // inlet velocity top channel property
        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropInletVelocityTop" );
        aParameterLists.set( "function_parameters", std::string( std::to_string( tRadiusTopChannel ) + "/" std::to_string( tYcTopChannel ) + "/-1.0" ) );
        aParameterLists.set( "value_function", "Func_Inlet_U" );

        // inlet velocity bottom channel property
        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropInletVelocityBottom" );
        aParameterLists.set( "function_parameters", std::string( std::to_string( tRadiusBottomChannel ) + "/" std::to_string( tYcBottomChannel ) + "/1.0" ) );
        aParameterLists.set( "value_function", "Func_Inlet_U" );

        // zero velocity property
        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropZeroVelocity" );
        aParameterLists.set( "function_parameters", "0.0;0.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // inlet temperature top channel property
        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropInletTempTop" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // inlet temperature bottom channel property
        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropInletTempBottom" );
        aParameterLists.set( "function_parameters", "0.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        //------------------------------------------------------------------------------
        // fill the constitutive model part of the parameter list

        // init CM counter
        moris::uint tCMCounter = 0;

        // fluid constitutive model
        aParameterLists( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        aParameterLists.set( "constitutive_name", "CMFluid" );
        aParameterLists.set( "constitutive_type",  fem::Constitutive_Type::FLUID_INCOMPRESSIBLE ) ;
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "VX,VY;P", "Velocity,Pressure" ) );
        aParameterLists.set( "properties", "PropFluidViscosity,Viscosity;PropFluidDensity,Density" );

        // fluid diffusion constitutive model
        aParameterLists( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        aParameterLists.set( "constitutive_name", "CMFluidDiffusion" );
        aParameterLists.set( "constitutive_type",  fem::Constitutive_Type::DIFF_LIN_ISO ) ;
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists.set( "properties", "PropFluidConductivity,Conductivity;PropFluidDensity,Density;PropFluidHeatCapacity,HeatCapacity" );

        // solid diffusion constitutive model
        aParameterLists( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        aParameterLists.set( "constitutive_name", "CMSolidDiffusion" );
        aParameterLists.set( "constitutive_type",  fem::Constitutive_Type::DIFF_LIN_ISO ) ;
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists.set( "properties", "PropSolidConductivity,Conductivity;PropSolidDensity,Density;PropSolidHeatCapacity,HeatCapacity" );

        //------------------------------------------------------------------------------
        // fill the stabilization parameter part of the parameter list

        // init SP counter
        moris::uint tSPCounter = 0;

        // NS SUPG stabilization parameter
        aParameterLists( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists.set( "stabilization_name", "SPSUPGVelocity" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::INCOMPRESSIBLE_FLOW ) ;
        aParameterLists.set( "function_parameters", "36.0" );
        aParameterLists.set( "leader_properties", "PropFluidViscosity,Viscosity;PropFluidDensity,Density" );
        aParameterLists.set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY;P", "Velocity,Pressure" ) );

        // advection SUPG stabilization parameter
        aParameterLists( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists.set( "stabilization_name", "SPSUPGTemp" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::SUPG_ADVECTION ) ;
        aParameterLists.set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
        aParameterLists.set( "leader_properties", "PropFluidConductivity,Conductivity" );

        // nitsche for fluid velocity stabilization parameter
        aParameterLists( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists.set( "stabilization_name", "SPNitscheVelocity" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::VELOCITY_DIRICHLET_NITSCHE ) ;
        aParameterLists.set( "function_parameters", "100.0/1.0" );
        aParameterLists.set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
        aParameterLists.set( "leader_properties", "PropFluidViscosity,Viscosity;PropFluidDensity,Density" );

        // nitsche for fluid temperature stabilization parameter 4
        aParameterLists( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists.set( "stabilization_name", "SPNitscheFluidTemp" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::DIRICHLET_NITSCHE ) ;
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "leader_properties", "PropFluidConductivity,Material" );

        // nitsche thermal interface stabilization parameter
        aParameterLists( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists.set( "stabilization_name", "SPInterfaceNitsche" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::NITSCHE_INTERFACE ) ;
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "leader_properties", "PropFluidConductivity,Material" );
        aParameterLists.set( "follower_properties", "PropSolidConductivity,Material" );

        // ghost viscous stabilization parameter
        aParameterLists( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists.set( "stabilization_name", "SPGPViscous" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::VISCOUS_GHOST ) ;
        aParameterLists.set( "function_parameters", "0.05" );
        aParameterLists.set( "leader_properties", "PropFluidViscosity,Viscosity" );

        // ghost convective stabilization parameter
        aParameterLists( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists.set( "stabilization_name", "SPGPVelocity" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::CONVECTIVE_GHOST ) ;
        aParameterLists.set( "function_parameters", "0.05" );
        aParameterLists.set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
        aParameterLists.set( "leader_properties", "PropFluidDensity,Density" );

        // ghost fluid pressure stabilization parameter
        aParameterLists( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists.set( "stabilization_name", "SPGPPressure" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::PRESSURE_GHOST ) ;
        aParameterLists.set( "function_parameters", "0.05/1.0" );
        aParameterLists.set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
        aParameterLists.set( "leader_properties", "PropFluidViscosity,Viscosity;PropFluidDensity,Density" );

        // ghost fluid temperature stabilization parameter
        aParameterLists( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists.set( "stabilization_name", "SPGPFluidTemp" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::GHOST_DISPL ) ;
        aParameterLists.set( "function_parameters", "0.05" );
        aParameterLists.set( "leader_properties", "PropFluidConductivity,Material" );

        // ghost solid temperature stabilization parameter
        aParameterLists( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists.set( "stabilization_name", "SPGPSolidTemp" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::GHOST_DISPL ) ;
        aParameterLists.set( "function_parameters", "0.05" );
        aParameterLists.set( "leader_properties", "PropSolidConductivity,Material" );

        //------------------------------------------------------------------------------
        // fill the IWG part of the parameter list

        // init IWG counter
        moris::uint tIWGCounter = 0;

        // NS incompressible velocity bulk IWG
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGVelocityBulk" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_BULK ) ;
        aParameterLists.set( "dof_residual", "VX,VY" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists.set( "stabilization_parameters", "SPSUPGVelocity,IncompressibleFlow" );
        aParameterLists.set( "mesh_set_names", "HMR_dummy_n_p34,HMR_dummy_c_p34,HMR_dummy_n_p46,HMR_dummy_c_p46" );

        // NS incompressible pressure bulk IWG
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGPressureBulk" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_BULK ) ;
        aParameterLists.set( "dof_residual", "P" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists.set( "stabilization_parameters", "SPSUPGVelocity,IncompressibleFlow" );
        aParameterLists.set( "mesh_set_names", "HMR_dummy_n_p34,HMR_dummy_c_p34,HMR_dummy_n_p46,HMR_dummy_c_p46" );

        // fluid diffusion bulk IWG
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGFluidDiffusionBulk" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_BULK ) ;
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists.set( "leader_constitutive_models", "CMFluidDiffusion,Diffusion" );
        aParameterLists.set( "mesh_set_names", "HMR_dummy_n_p34,HMR_dummy_c_p34,HMR_dummy_n_p46,HMR_dummy_c_p46" );

        // solid diffusion bulk IWG
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGSolidDiffusionBulk" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_BULK ) ;
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists.set( "leader_constitutive_models", "CMSolidDiffusion,Diffusion" );
        aParameterLists.set( "mesh_set_names", "HMR_dummy_n_p42,HMR_dummy_c_p42" );

        // fluid advection bulk IWG
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGAdvectionBulk" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::ADVECTION_BULK ) ;
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists.set( "leader_constitutive_models", "CMFluidDiffusion,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPSUPGTemp,SUPG" );
        aParameterLists.set( "mesh_set_names", "HMR_dummy_n_p34,HMR_dummy_c_p34,HMR_dummy_n_p46,HMR_dummy_c_p46" );

        // top inlet fluid velocity IWG (velocity part)
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGInVelocityTop" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        aParameterLists.set( "dof_residual", "VX,VY" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists.set( "leader_properties", "PropInletVelocityTop,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheVelocity,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", "iside_b0_46_b1_47" );

        // top inlet fluid velocity IWG (pressure part)
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGInPressureTop" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        aParameterLists.set( "dof_residual", "P" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists.set( "leader_properties", "PropInletVelocityTop,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists.set( "mesh_set_names", "iside_b0_46_b1_47" );

        // bottom inlet fluid velocity IWG (velocity part)
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGInVelocityBottom" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        aParameterLists.set( "dof_residual", "VX,VY" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists.set( "leader_properties", "PropInletVelocityBottom,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheVelocity,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", "iside_b0_34_b1_32" );

        // bottom inlet fluid velocity IWG (pressure part)
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGInPressureBottom" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        aParameterLists.set( "dof_residual", "P" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists.set( "leader_properties", "PropInletVelocityBottom,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists.set( "mesh_set_names", "iside_b0_34_b1_32" );

        // wall fluid velocity IWG (velocity part)
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGWallVelocity" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        aParameterLists.set( "dof_residual", "VX,VY" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists.set( "leader_properties", "PropZeroVelocity,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheVelocity,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", "iside_b0_34_b1_2,iside_b0_34_b1_42,iside_b0_46_b1_42,iside_b0_46_b1_62" );

        // wall fluid velocity IWG (pressure part)
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGWallPressure" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        aParameterLists.set( "dof_residual", "P" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists.set( "leader_properties", "PropZeroVelocity,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists.set( "mesh_set_names", "iside_b0_34_b1_2,iside_b0_34_b1_42,iside_b0_46_b1_42,iside_b0_46_b1_62" );

        // top inlet fluid temperature IWG
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGInletTempTop" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists.set( "leader_properties", "PropInletTempTop,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMFluidDiffusion,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheFluidTemp,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", "iside_b0_46_b1_47" );

        // bottom inlet fluid temperature IWG
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGInletTempBottom" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists.set( "leader_properties", "PropInletTempBottom,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMFluidDiffusion,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheFluidTemp,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", "iside_b0_34_b1_32" );

        // fluid/solid thermal interface IWG
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGInterfaceFluidSolid" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE ) ;
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists.set( "follower_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists.set( "leader_constitutive_models", "CMFluidDiffusion,Diffusion" );
        aParameterLists.set( "follower_constitutive_models", "CMSolidDiffusion,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPInterfaceNitsche,NitscheInterface" );
        aParameterLists.set( "mesh_set_names", "dbl_iside_p0_34_p1_42,dbl_iside_p0_46_p1_42" );

        // ghost viscous IWG
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGGPViscous" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
        aParameterLists.set( "dof_residual", "VX,VY" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists.set( "follower_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists.set( "stabilization_parameters", "SPGPViscous,GhostSP" );
        aParameterLists.set( "mesh_set_names", "ghost_p34,ghost_p46" );

        // ghost convective IWG
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGGPConvective" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
        aParameterLists.set( "dof_residual", "VX,VY" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists.set( "follower_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists.set( "stabilization_parameters", "SPGPVelocity,GhostSP" );
        aParameterLists.set( "mesh_set_names", "ghost_p34,ghost_p46" );

        // ghost pressure IWG
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGGPPressure" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
        aParameterLists.set( "dof_residual", "P" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists.set( "follower_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists.set( "stabilization_parameters", "SPGPPressure,GhostSP" );
        aParameterLists.set( "mesh_set_names", "ghost_p34,ghost_p46" );

        // ghost fluid temperature IWG
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGGPFluidTemp" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists.set( "follower_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists.set( "stabilization_parameters", "SPGPFluidTemp,GhostSP" );
        aParameterLists.set( "mesh_set_names", "ghost_p34,ghost_p46" );

        // ghost solid temperature IWG
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGGPSolidTemp" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists.set( "follower_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists.set( "stabilization_parameters", "SPGPSolidTemp,GhostSP" );
        aParameterLists.set( "mesh_set_names", "ghost_p42" );

        //------------------------------------------------------------------------------
        // fill the IQI part of the parameter list

        // init IQI counter
        moris::uint tIQICounter = 0;

        // VX IQI
        aParameterLists( 4 ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists.set( "IQI_name", "IQIBulkVX" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists.set( "dof_quantity", "VX,VY" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "mesh_set_names", "HMR_dummy_n_p34,HMR_dummy_c_p34,HMR_dummy_n_p46,HMR_dummy_c_p46,HMR_dummy_n_p42,HMR_dummy_c_p42" );

        // VY IQI
        aParameterLists( 4 ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists.set( "IQI_name", "IQIBulkVY" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists.set( "dof_quantity", "VX,VY" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY" );
        aParameterLists.set( "vectorial_field_index", 1 );
        aParameterLists.set( "mesh_set_names", "HMR_dummy_n_p34,HMR_dummy_c_p34,HMR_dummy_n_p46,HMR_dummy_c_p46,HMR_dummy_n_p42,HMR_dummy_c_p42" );

        // P IQI
        aParameterLists( 4 ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists.set( "IQI_name", "IQIBulkP" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists.set( "dof_quantity", "P" );
        aParameterLists.set( "leader_dof_dependencies", "P" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "mesh_set_names", "HMR_dummy_n_p34,HMR_dummy_c_p34,HMR_dummy_n_p46,HMR_dummy_c_p46,HMR_dummy_n_p42,HMR_dummy_c_p42" );

        // TEMP IQI
        aParameterLists( 4 ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists.set( "IQI_name", "IQIBulkTEMP" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists.set( "dof_quantity", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "mesh_set_names", "HMR_dummy_n_p34,HMR_dummy_c_p34,HMR_dummy_n_p46,HMR_dummy_c_p46,HMR_dummy_n_p42,HMR_dummy_c_p42" );

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        aParameterLists( 5 ).push_back( prm::create_computation_parameter_list() );
    }

    // SOL parameter list
    void
    SOLParameterList( Vector< Vector< ParameterList > >& aParameterLists )
    {

        aParameterLists( 0 ).push_back( add_parameter_list( sol::SolverType::AMESOS_IMPL ) );

        aParameterLists( 1 ).push_back( moris::prm::create_linear_solver_parameter_list() );

        aParameterLists( 2 ).push_back( moris::prm::create_nonlinear_algorithm_parameter_list() );
        aParameterLists.set( "NLA_rel_res_norm_drop", 1e-07 );

        aParameterLists( 3 ).push_back( moris::prm::create_nonlinear_solver_parameter_list() );
        aParameterLists.set( "NLA_DofTypes", "VX,VY;P;TEMP" );

        aParameterLists( 4 ).push_back( moris::prm::create_time_solver_algorithm_parameter_list() );
        aParameterLists.set( "TSA_Num_Time_Steps", 20 );
        aParameterLists.set( "TSA_Time_Frame", 0.4 );

        aParameterLists( 5 ).push_back( moris::prm::create_time_solver_parameter_list() );
        aParameterLists.set( "TSA_DofTypes", "VX,VY;P;TEMP" );
        aParameterLists.set( "TSA_Initialize_Sol_Vec", "VX,1E-4;VY,1E-4;P,0.0;TEMP,0.0" );
        aParameterLists.set( "TSA_Output_Indices", "0" );
        aParameterLists.set( "TSA_Output_Criteria", "Output_Criterion" );
        aParameterLists.set( "TSA_time_level_per_type", "VX,2;VY,2;P,2;TEMP,2" );

        aParameterLists( 6 ).push_back( moris::prm::create_solver_warehouse_parameterlist() );

        aParameterLists( 7 ).push_back(  sol::PreconditionerType::NONE );
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
        aParameterLists.set( "File_Name", std::pair< std::string, std::string >( "./", "2_Channels_Temp.exo" ) );
        aParameterLists.set( "Mesh_Type",  vis::VIS_Mesh_Type::STANDARD ) ;
        aParameterLists.set( "Set_Names", "HMR_dummy_n_p34,HMR_dummy_c_p34,HMR_dummy_n_p46,HMR_dummy_c_p46,HMR_dummy_n_p42,HMR_dummy_c_p42" );
        aParameterLists.set( "Field_Names", "VX,VY,P,TEMP" );
        aParameterLists.set( "Field_Type", "NODAL,NODAL,NODAL,NODAL" );
        aParameterLists.set( "IQI_Names", "IQIBulkVX,IQIBulkVY,IQIBulkP,IQIBulkTEMP" );
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
