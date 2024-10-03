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
    OPTParameterList( Vector< Vector< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).push_back( prm::create_opt_problem_parameter_list() );

        tParameterlist( 0 ).set( "is_optimization_problem", false );
    }

    // HMR parameter list
    void
    HMRParameterList( Vector< Vector< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).push_back( prm::create_hmr_parameter_list() );

        tParameterlist( 0 ).set( "number_of_elements_per_dimension", "62,31" );
        tParameterlist( 0 ).set( "domain_dimensions", "3,1" );
        tParameterlist( 0 ).set( "domain_offset", "-0.835,-0.186" );
        tParameterlist( 0 ).set( "domain_sidesets", "1,2,3,4" );
        tParameterlist( 0 ).set( "lagrange_output_meshes", "0" );

        tParameterlist( 0 ).set( "lagrange_orders", "1" );
        tParameterlist( 0 ).set( "lagrange_pattern", "0" );
        tParameterlist( 0 ).set( "bspline_orders", "1" );
        tParameterlist( 0 ).set( "bspline_pattern", "0" );

        tParameterlist( 0 ).set( "lagrange_to_bspline", "0" );

        tParameterlist( 0 ).set( "truncate_bsplines", 1 );
        tParameterlist( 0 ).set( "refinement_buffer", 3 );
        tParameterlist( 0 ).set( "staircase_buffer", 3 );
        tParameterlist( 0 ).set( "initial_refinement", "0" );
        tParameterlist( 0 ).set( "initial_refinement_pattern", "0" );

        tParameterlist( 0 ).set( "use_multigrid", 0 );
        tParameterlist( 0 ).set( "severity_level", 0 );

        tParameterlist( 0 ).set( "adaptive_refinement_level", 0 );
    }

    // XTK parameter list
    void
    XTKParameterList( Vector< Vector< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).push_back( prm::create_xtk_parameter_list() );
        tParameterlist( 0 ).set( "decompose", true );
        tParameterlist( 0 ).set( "decomposition_type", "conformal" );
        tParameterlist( 0 ).set( "enrich", true );
        tParameterlist( 0 ).set( "basis_rank", "bspline" );
        tParameterlist( 0 ).set( "enrich_mesh_indices", "0" );
        tParameterlist( 0 ).set( "ghost_stab", true );
        tParameterlist( 0 ).set( "multigrid", false );
        tParameterlist( 0 ).set( "high_to_low_dbl_side_sets", true );
        tParameterlist( 0 ).set( "print_enriched_ig_mesh", true );
        tParameterlist( 0 ).set( "exodus_output_XTK_ig_mesh", false );
    }

    // GEN parameter list
    void
    GENParameterList( Vector< Vector< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 3 );
        tParameterlist( 0 ).push_back( prm::create_gen_parameter_list() );


        // Geometry parameter lists
        tParameterlist( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        tParameterlist( 1 ).set( "field_function_name", "Func_Plane" );
        tParameterlist( 1 ).set( "constant_parameters", 0,0,1,0 );

        tParameterlist( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        tParameterlist( 1 ).set( "field_function_name", "Func_Plane" );
        tParameterlist( 1 ).set( "constant_parameters", 0,0,1,0.41 );

        tParameterlist( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        tParameterlist( 1 ).set( "field_function_name", "Func_Plane" );
        tParameterlist( 1 ).set( "constant_parameters", 0,0,1,0.18 );

        tParameterlist( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        tParameterlist( 1 ).set( "field_function_name", "Func_Plane" );
        tParameterlist( 1 ).set( "constant_parameters", 0,0,1,0.23 );

        tParameterlist( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        tParameterlist( 1 ).set( "field_function_name", "Func_Plane" );
        tParameterlist( 1 ).set( "constant_parameters", 1,0,0,0 );

        tParameterlist( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        tParameterlist( 1 ).set( "field_function_name", "Func_Plane" );
        tParameterlist( 1 ).set( "constant_parameters", 1,1.845,0,0 );
    }

    // FEM parameter list
    void
    FEMParameterList( Vector< Vector< ParameterList > >& tParameterList )
    {
        // create a cell of cell of parameter list for fem
        tParameterList.resize( 8 );

        //------------------------------------------------------------------------------
        // fill the property part of the parameter list


        // fluid viscosity property
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropFluidViscosity" );
        tParameterList( 0 ).set( "function_parameters", "1.0" );
        tParameterList( 0 ).set( "value_function", "Func_Const" );

        // fluid density property
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropFluidDensity" );
        tParameterList( 0 ).set( "function_parameters", "1.0" );
        tParameterList( 0 ).set( "value_function", "Func_Const" );

        // fluid conductivity property
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropFluidConductivity" );
        tParameterList( 0 ).set( "function_parameters", "500.0" );
        tParameterList( 0 ).set( "value_function", "Func_Const" );

        // fluid heat capacity property
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropFluidHeatCapacity" );
        tParameterList( 0 ).set( "function_parameters", "0.001" );
        tParameterList( 0 ).set( "value_function", "Func_Const" );

        // solid density property
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropSolidDensity" );
        tParameterList( 0 ).set( "function_parameters", "2.0" );
        tParameterList( 0 ).set( "value_function", "Func_Const" );

        // solid conductivity property
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropSolidConductivity" );
        tParameterList( 0 ).set( "function_parameters", "1.0" );
        tParameterList( 0 ).set( "value_function", "Func_Const" );

        // solid heat capacity property
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropSolidHeatCapacity" );
        tParameterList( 0 ).set( "function_parameters", "0.001" );
        tParameterList( 0 ).set( "value_function", "Func_Const" );

        // inlet velocity top channel property
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropInletVelocityTop" );
        tParameterList( 0 ).set( "function_parameters", std::string( std::to_string( tRadiusTopChannel ) + "/" std::to_string( tYcTopChannel ) + "/-1.0" ) );
        tParameterList( 0 ).set( "value_function", "Func_Inlet_U" );

        // inlet velocity bottom channel property
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropInletVelocityBottom" );
        tParameterList( 0 ).set( "function_parameters", std::string( std::to_string( tRadiusBottomChannel ) + "/" std::to_string( tYcBottomChannel ) + "/1.0" ) );
        tParameterList( 0 ).set( "value_function", "Func_Inlet_U" );

        // zero velocity property
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropZeroVelocity" );
        tParameterList( 0 ).set( "function_parameters", "0.0;0.0" );
        tParameterList( 0 ).set( "value_function", "Func_Const" );

        // inlet temperature top channel property
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropInletTempTop" );
        tParameterList( 0 ).set( "function_parameters", "1.0" );
        tParameterList( 0 ).set( "value_function", "Func_Const" );

        // inlet temperature bottom channel property
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropInletTempBottom" );
        tParameterList( 0 ).set( "function_parameters", "0.0" );
        tParameterList( 0 ).set( "value_function", "Func_Const" );

        //------------------------------------------------------------------------------
        // fill the constitutive model part of the parameter list

        // init CM counter
        moris::uint tCMCounter = 0;

        // fluid constitutive model
        tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 ).set( "constitutive_name", "CMFluid" );
        tParameterList( 1 ).set( "constitutive_type",  fem::Constitutive_Type::FLUID_INCOMPRESSIBLE ) ;
        tParameterList( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( "VX,VY;P", "Velocity,Pressure" ) );
        tParameterList( 1 ).set( "properties", "PropFluidViscosity,Viscosity;PropFluidDensity,Density" );

        // fluid diffusion constitutive model
        tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 ).set( "constitutive_name", "CMFluidDiffusion" );
        tParameterList( 1 ).set( "constitutive_type",  fem::Constitutive_Type::DIFF_LIN_ISO ) ;
        tParameterList( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( 1 ).set( "properties", "PropFluidConductivity,Conductivity;PropFluidDensity,Density;PropFluidHeatCapacity,HeatCapacity" );

        // solid diffusion constitutive model
        tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 ).set( "constitutive_name", "CMSolidDiffusion" );
        tParameterList( 1 ).set( "constitutive_type",  fem::Constitutive_Type::DIFF_LIN_ISO ) ;
        tParameterList( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( 1 ).set( "properties", "PropSolidConductivity,Conductivity;PropSolidDensity,Density;PropSolidHeatCapacity,HeatCapacity" );

        //------------------------------------------------------------------------------
        // fill the stabilization parameter part of the parameter list

        // init SP counter
        moris::uint tSPCounter = 0;

        // NS SUPG stabilization parameter
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 ).set( "stabilization_name", "SPSUPGVelocity" );
        tParameterList( 2 ).set( "stabilization_type",  fem::Stabilization_Type::INCOMPRESSIBLE_FLOW ) ;
        tParameterList( 2 ).set( "function_parameters", "36.0" );
        tParameterList( 2 ).set( "leader_properties", "PropFluidViscosity,Viscosity;PropFluidDensity,Density" );
        tParameterList( 2 ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY;P", "Velocity,Pressure" ) );

        // advection SUPG stabilization parameter
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 ).set( "stabilization_name", "SPSUPGTemp" );
        tParameterList( 2 ).set( "stabilization_type",  fem::Stabilization_Type::SUPG_ADVECTION ) ;
        tParameterList( 2 ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
        tParameterList( 2 ).set( "leader_properties", "PropFluidConductivity,Conductivity" );

        // nitsche for fluid velocity stabilization parameter
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 ).set( "stabilization_name", "SPNitscheVelocity" );
        tParameterList( 2 ).set( "stabilization_type",  fem::Stabilization_Type::VELOCITY_DIRICHLET_NITSCHE ) ;
        tParameterList( 2 ).set( "function_parameters", "100.0/1.0" );
        tParameterList( 2 ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
        tParameterList( 2 ).set( "leader_properties", "PropFluidViscosity,Viscosity;PropFluidDensity,Density" );

        // nitsche for fluid temperature stabilization parameter 4
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 ).set( "stabilization_name", "SPNitscheFluidTemp" );
        tParameterList( 2 ).set( "stabilization_type",  fem::Stabilization_Type::DIRICHLET_NITSCHE ) ;
        tParameterList( 2 ).set( "function_parameters", "100.0" );
        tParameterList( 2 ).set( "leader_properties", "PropFluidConductivity,Material" );

        // nitsche thermal interface stabilization parameter
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 ).set( "stabilization_name", "SPInterfaceNitsche" );
        tParameterList( 2 ).set( "stabilization_type",  fem::Stabilization_Type::NITSCHE_INTERFACE ) ;
        tParameterList( 2 ).set( "function_parameters", "100.0" );
        tParameterList( 2 ).set( "leader_properties", "PropFluidConductivity,Material" );
        tParameterList( 2 ).set( "follower_properties", "PropSolidConductivity,Material" );

        // ghost viscous stabilization parameter
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 ).set( "stabilization_name", "SPGPViscous" );
        tParameterList( 2 ).set( "stabilization_type",  fem::Stabilization_Type::VISCOUS_GHOST ) ;
        tParameterList( 2 ).set( "function_parameters", "0.05" );
        tParameterList( 2 ).set( "leader_properties", "PropFluidViscosity,Viscosity" );

        // ghost convective stabilization parameter
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 ).set( "stabilization_name", "SPGPVelocity" );
        tParameterList( 2 ).set( "stabilization_type",  fem::Stabilization_Type::CONVECTIVE_GHOST ) ;
        tParameterList( 2 ).set( "function_parameters", "0.05" );
        tParameterList( 2 ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
        tParameterList( 2 ).set( "leader_properties", "PropFluidDensity,Density" );

        // ghost fluid pressure stabilization parameter
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 ).set( "stabilization_name", "SPGPPressure" );
        tParameterList( 2 ).set( "stabilization_type",  fem::Stabilization_Type::PRESSURE_GHOST ) ;
        tParameterList( 2 ).set( "function_parameters", "0.05/1.0" );
        tParameterList( 2 ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
        tParameterList( 2 ).set( "leader_properties", "PropFluidViscosity,Viscosity;PropFluidDensity,Density" );

        // ghost fluid temperature stabilization parameter
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 ).set( "stabilization_name", "SPGPFluidTemp" );
        tParameterList( 2 ).set( "stabilization_type",  fem::Stabilization_Type::GHOST_DISPL ) ;
        tParameterList( 2 ).set( "function_parameters", "0.05" );
        tParameterList( 2 ).set( "leader_properties", "PropFluidConductivity,Material" );

        // ghost solid temperature stabilization parameter
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 ).set( "stabilization_name", "SPGPSolidTemp" );
        tParameterList( 2 ).set( "stabilization_type",  fem::Stabilization_Type::GHOST_DISPL ) ;
        tParameterList( 2 ).set( "function_parameters", "0.05" );
        tParameterList( 2 ).set( "leader_properties", "PropSolidConductivity,Material" );

        //------------------------------------------------------------------------------
        // fill the IWG part of the parameter list

        // init IWG counter
        moris::uint tIWGCounter = 0;

        // NS incompressible velocity bulk IWG
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 ).set( "IWG_name", "IWGVelocityBulk" );
        tParameterList( 3 ).set( "IWG_type",  fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_BULK ) ;
        tParameterList( 3 ).set( "dof_residual", "VX,VY" );
        tParameterList( 3 ).set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        tParameterList( 3 ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        tParameterList( 3 ).set( "stabilization_parameters", "SPSUPGVelocity,IncompressibleFlow" );
        tParameterList( 3 ).set( "mesh_set_names", "HMR_dummy_n_p34,HMR_dummy_c_p34,HMR_dummy_n_p46,HMR_dummy_c_p46" );

        // NS incompressible pressure bulk IWG
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 ).set( "IWG_name", "IWGPressureBulk" );
        tParameterList( 3 ).set( "IWG_type",  fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_BULK ) ;
        tParameterList( 3 ).set( "dof_residual", "P" );
        tParameterList( 3 ).set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        tParameterList( 3 ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        tParameterList( 3 ).set( "stabilization_parameters", "SPSUPGVelocity,IncompressibleFlow" );
        tParameterList( 3 ).set( "mesh_set_names", "HMR_dummy_n_p34,HMR_dummy_c_p34,HMR_dummy_n_p46,HMR_dummy_c_p46" );

        // fluid diffusion bulk IWG
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 ).set( "IWG_name", "IWGFluidDiffusionBulk" );
        tParameterList( 3 ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_BULK ) ;
        tParameterList( 3 ).set( "dof_residual", "TEMP" );
        tParameterList( 3 ).set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        tParameterList( 3 ).set( "leader_constitutive_models", "CMFluidDiffusion,Diffusion" );
        tParameterList( 3 ).set( "mesh_set_names", "HMR_dummy_n_p34,HMR_dummy_c_p34,HMR_dummy_n_p46,HMR_dummy_c_p46" );

        // solid diffusion bulk IWG
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 ).set( "IWG_name", "IWGSolidDiffusionBulk" );
        tParameterList( 3 ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_BULK ) ;
        tParameterList( 3 ).set( "dof_residual", "TEMP" );
        tParameterList( 3 ).set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        tParameterList( 3 ).set( "leader_constitutive_models", "CMSolidDiffusion,Diffusion" );
        tParameterList( 3 ).set( "mesh_set_names", "HMR_dummy_n_p42,HMR_dummy_c_p42" );

        // fluid advection bulk IWG
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 ).set( "IWG_name", "IWGAdvectionBulk" );
        tParameterList( 3 ).set( "IWG_type",  fem::IWG_Type::ADVECTION_BULK ) ;
        tParameterList( 3 ).set( "dof_residual", "TEMP" );
        tParameterList( 3 ).set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        tParameterList( 3 ).set( "leader_constitutive_models", "CMFluidDiffusion,Diffusion" );
        tParameterList( 3 ).set( "stabilization_parameters", "SPSUPGTemp,SUPG" );
        tParameterList( 3 ).set( "mesh_set_names", "HMR_dummy_n_p34,HMR_dummy_c_p34,HMR_dummy_n_p46,HMR_dummy_c_p46" );

        // top inlet fluid velocity IWG (velocity part)
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 ).set( "IWG_name", "IWGInVelocityTop" );
        tParameterList( 3 ).set( "IWG_type",  fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        tParameterList( 3 ).set( "dof_residual", "VX,VY" );
        tParameterList( 3 ).set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        tParameterList( 3 ).set( "leader_properties", "PropInletVelocityTop,Dirichlet" );
        tParameterList( 3 ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        tParameterList( 3 ).set( "stabilization_parameters", "SPNitscheVelocity,DirichletNitsche" );
        tParameterList( 3 ).set( "mesh_set_names", "iside_b0_46_b1_47" );

        // top inlet fluid velocity IWG (pressure part)
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 ).set( "IWG_name", "IWGInPressureTop" );
        tParameterList( 3 ).set( "IWG_type",  fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        tParameterList( 3 ).set( "dof_residual", "P" );
        tParameterList( 3 ).set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        tParameterList( 3 ).set( "leader_properties", "PropInletVelocityTop,Dirichlet" );
        tParameterList( 3 ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        tParameterList( 3 ).set( "mesh_set_names", "iside_b0_46_b1_47" );

        // bottom inlet fluid velocity IWG (velocity part)
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 ).set( "IWG_name", "IWGInVelocityBottom" );
        tParameterList( 3 ).set( "IWG_type",  fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        tParameterList( 3 ).set( "dof_residual", "VX,VY" );
        tParameterList( 3 ).set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        tParameterList( 3 ).set( "leader_properties", "PropInletVelocityBottom,Dirichlet" );
        tParameterList( 3 ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        tParameterList( 3 ).set( "stabilization_parameters", "SPNitscheVelocity,DirichletNitsche" );
        tParameterList( 3 ).set( "mesh_set_names", "iside_b0_34_b1_32" );

        // bottom inlet fluid velocity IWG (pressure part)
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 ).set( "IWG_name", "IWGInPressureBottom" );
        tParameterList( 3 ).set( "IWG_type",  fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        tParameterList( 3 ).set( "dof_residual", "P" );
        tParameterList( 3 ).set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        tParameterList( 3 ).set( "leader_properties", "PropInletVelocityBottom,Dirichlet" );
        tParameterList( 3 ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        tParameterList( 3 ).set( "mesh_set_names", "iside_b0_34_b1_32" );

        // wall fluid velocity IWG (velocity part)
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 ).set( "IWG_name", "IWGWallVelocity" );
        tParameterList( 3 ).set( "IWG_type",  fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        tParameterList( 3 ).set( "dof_residual", "VX,VY" );
        tParameterList( 3 ).set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        tParameterList( 3 ).set( "leader_properties", "PropZeroVelocity,Dirichlet" );
        tParameterList( 3 ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        tParameterList( 3 ).set( "stabilization_parameters", "SPNitscheVelocity,DirichletNitsche" );
        tParameterList( 3 ).set( "mesh_set_names", "iside_b0_34_b1_2,iside_b0_34_b1_42,iside_b0_46_b1_42,iside_b0_46_b1_62" );

        // wall fluid velocity IWG (pressure part)
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 ).set( "IWG_name", "IWGWallPressure" );
        tParameterList( 3 ).set( "IWG_type",  fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        tParameterList( 3 ).set( "dof_residual", "P" );
        tParameterList( 3 ).set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        tParameterList( 3 ).set( "leader_properties", "PropZeroVelocity,Dirichlet" );
        tParameterList( 3 ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        tParameterList( 3 ).set( "mesh_set_names", "iside_b0_34_b1_2,iside_b0_34_b1_42,iside_b0_46_b1_42,iside_b0_46_b1_62" );

        // top inlet fluid temperature IWG
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 ).set( "IWG_name", "IWGInletTempTop" );
        tParameterList( 3 ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        tParameterList( 3 ).set( "dof_residual", "TEMP" );
        tParameterList( 3 ).set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        tParameterList( 3 ).set( "leader_properties", "PropInletTempTop,Dirichlet" );
        tParameterList( 3 ).set( "leader_constitutive_models", "CMFluidDiffusion,Diffusion" );
        tParameterList( 3 ).set( "stabilization_parameters", "SPNitscheFluidTemp,DirichletNitsche" );
        tParameterList( 3 ).set( "mesh_set_names", "iside_b0_46_b1_47" );

        // bottom inlet fluid temperature IWG
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 ).set( "IWG_name", "IWGInletTempBottom" );
        tParameterList( 3 ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        tParameterList( 3 ).set( "dof_residual", "TEMP" );
        tParameterList( 3 ).set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        tParameterList( 3 ).set( "leader_properties", "PropInletTempBottom,Dirichlet" );
        tParameterList( 3 ).set( "leader_constitutive_models", "CMFluidDiffusion,Diffusion" );
        tParameterList( 3 ).set( "stabilization_parameters", "SPNitscheFluidTemp,DirichletNitsche" );
        tParameterList( 3 ).set( "mesh_set_names", "iside_b0_34_b1_32" );

        // fluid/solid thermal interface IWG
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 ).set( "IWG_name", "IWGInterfaceFluidSolid" );
        tParameterList( 3 ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE ) ;
        tParameterList( 3 ).set( "dof_residual", "TEMP" );
        tParameterList( 3 ).set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        tParameterList( 3 ).set( "follower_dof_dependencies", "VX,VY;P;TEMP" );
        tParameterList( 3 ).set( "leader_constitutive_models", "CMFluidDiffusion,Diffusion" );
        tParameterList( 3 ).set( "follower_constitutive_models", "CMSolidDiffusion,Diffusion" );
        tParameterList( 3 ).set( "stabilization_parameters", "SPInterfaceNitsche,NitscheInterface" );
        tParameterList( 3 ).set( "mesh_set_names", "dbl_iside_p0_34_p1_42,dbl_iside_p0_46_p1_42" );

        // ghost viscous IWG
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 ).set( "IWG_name", "IWGGPViscous" );
        tParameterList( 3 ).set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
        tParameterList( 3 ).set( "dof_residual", "VX,VY" );
        tParameterList( 3 ).set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        tParameterList( 3 ).set( "follower_dof_dependencies", "VX,VY;P;TEMP" );
        tParameterList( 3 ).set( "stabilization_parameters", "SPGPViscous,GhostSP" );
        tParameterList( 3 ).set( "mesh_set_names", "ghost_p34,ghost_p46" );

        // ghost convective IWG
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 ).set( "IWG_name", "IWGGPConvective" );
        tParameterList( 3 ).set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
        tParameterList( 3 ).set( "dof_residual", "VX,VY" );
        tParameterList( 3 ).set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        tParameterList( 3 ).set( "follower_dof_dependencies", "VX,VY;P;TEMP" );
        tParameterList( 3 ).set( "stabilization_parameters", "SPGPVelocity,GhostSP" );
        tParameterList( 3 ).set( "mesh_set_names", "ghost_p34,ghost_p46" );

        // ghost pressure IWG
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 ).set( "IWG_name", "IWGGPPressure" );
        tParameterList( 3 ).set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
        tParameterList( 3 ).set( "dof_residual", "P" );
        tParameterList( 3 ).set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        tParameterList( 3 ).set( "follower_dof_dependencies", "VX,VY;P;TEMP" );
        tParameterList( 3 ).set( "stabilization_parameters", "SPGPPressure,GhostSP" );
        tParameterList( 3 ).set( "mesh_set_names", "ghost_p34,ghost_p46" );

        // ghost fluid temperature IWG
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 ).set( "IWG_name", "IWGGPFluidTemp" );
        tParameterList( 3 ).set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
        tParameterList( 3 ).set( "dof_residual", "TEMP" );
        tParameterList( 3 ).set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        tParameterList( 3 ).set( "follower_dof_dependencies", "VX,VY;P;TEMP" );
        tParameterList( 3 ).set( "stabilization_parameters", "SPGPFluidTemp,GhostSP" );
        tParameterList( 3 ).set( "mesh_set_names", "ghost_p34,ghost_p46" );

        // ghost solid temperature IWG
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 ).set( "IWG_name", "IWGGPSolidTemp" );
        tParameterList( 3 ).set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
        tParameterList( 3 ).set( "dof_residual", "TEMP" );
        tParameterList( 3 ).set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        tParameterList( 3 ).set( "follower_dof_dependencies", "VX,VY;P;TEMP" );
        tParameterList( 3 ).set( "stabilization_parameters", "SPGPSolidTemp,GhostSP" );
        tParameterList( 3 ).set( "mesh_set_names", "ghost_p42" );

        //------------------------------------------------------------------------------
        // fill the IQI part of the parameter list

        // init IQI counter
        moris::uint tIQICounter = 0;

        // VX IQI
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 ).set( "IQI_name", "IQIBulkVX" );
        tParameterList( 4 ).set( "IQI_type",  fem::IQI_Type::DOF ) ;
        tParameterList( 4 ).set( "dof_quantity", "VX,VY" );
        tParameterList( 4 ).set( "leader_dof_dependencies", "VX,VY" );
        tParameterList( 4 ).set( "vectorial_field_index", 0 );
        tParameterList( 4 ).set( "mesh_set_names", "HMR_dummy_n_p34,HMR_dummy_c_p34,HMR_dummy_n_p46,HMR_dummy_c_p46,HMR_dummy_n_p42,HMR_dummy_c_p42" );

        // VY IQI
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 ).set( "IQI_name", "IQIBulkVY" );
        tParameterList( 4 ).set( "IQI_type",  fem::IQI_Type::DOF ) ;
        tParameterList( 4 ).set( "dof_quantity", "VX,VY" );
        tParameterList( 4 ).set( "leader_dof_dependencies", "VX,VY" );
        tParameterList( 4 ).set( "vectorial_field_index", 1 );
        tParameterList( 4 ).set( "mesh_set_names", "HMR_dummy_n_p34,HMR_dummy_c_p34,HMR_dummy_n_p46,HMR_dummy_c_p46,HMR_dummy_n_p42,HMR_dummy_c_p42" );

        // P IQI
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 ).set( "IQI_name", "IQIBulkP" );
        tParameterList( 4 ).set( "IQI_type",  fem::IQI_Type::DOF ) ;
        tParameterList( 4 ).set( "dof_quantity", "P" );
        tParameterList( 4 ).set( "leader_dof_dependencies", "P" );
        tParameterList( 4 ).set( "vectorial_field_index", 0 );
        tParameterList( 4 ).set( "mesh_set_names", "HMR_dummy_n_p34,HMR_dummy_c_p34,HMR_dummy_n_p46,HMR_dummy_c_p46,HMR_dummy_n_p42,HMR_dummy_c_p42" );

        // TEMP IQI
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 ).set( "IQI_name", "IQIBulkTEMP" );
        tParameterList( 4 ).set( "IQI_type",  fem::IQI_Type::DOF ) ;
        tParameterList( 4 ).set( "dof_quantity", "TEMP" );
        tParameterList( 4 ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( 4 ).set( "vectorial_field_index", 0 );
        tParameterList( 4 ).set( "mesh_set_names", "HMR_dummy_n_p34,HMR_dummy_c_p34,HMR_dummy_n_p46,HMR_dummy_c_p46,HMR_dummy_n_p42,HMR_dummy_c_p42" );

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        tParameterList( 5 ).push_back( prm::create_computation_parameter_list() );
    }

    // SOL parameter list
    void
    SOLParameterList( Vector< Vector< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 8 );


        tParameterlist( 0 ).push_back( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL ) );

        tParameterlist( 1 ).push_back( moris::prm::create_linear_solver_parameter_list() );

        tParameterlist( 2 ).push_back( moris::prm::create_nonlinear_algorithm_parameter_list() );
        tParameterlist( 2 ).set( "NLA_rel_res_norm_drop", 1e-07 );

        tParameterlist( 3 ).push_back( moris::prm::create_nonlinear_solver_parameter_list() );
        tParameterlist( 3 ).set( "NLA_DofTypes", "VX,VY;P;TEMP" );

        tParameterlist( 4 ).push_back( moris::prm::create_time_solver_algorithm_parameter_list() );
        tParameterlist( 4 ).set( "TSA_Num_Time_Steps", 20 );
        tParameterlist( 4 ).set( "TSA_Time_Frame", 0.4 );

        tParameterlist( 5 ).push_back( moris::prm::create_time_solver_parameter_list() );
        tParameterlist( 5 ).set( "TSA_DofTypes", "VX,VY;P;TEMP" );
        tParameterlist( 5 ).set( "TSA_Initialize_Sol_Vec", "VX,1E-4;VY,1E-4;P,0.0;TEMP,0.0" );
        tParameterlist( 5 ).set( "TSA_Output_Indices", "0" );
        tParameterlist( 5 ).set( "TSA_Output_Criteria", "Output_Criterion" );
        tParameterlist( 5 ).set( "TSA_time_level_per_type", "VX,2;VY,2;P,2;TEMP,2" );

        tParameterlist( 6 ).push_back( moris::prm::create_solver_warehouse_parameterlist() );

        tParameterlist( 7 ).push_back( moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::NONE ) );
    }

    // MSI parametr list
    void
    MSIParameterList( Vector< Vector< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).push_back( prm::create_msi_parameter_list() );
    }

    // VIS parameter list
    void
    VISParameterList( Vector< Vector< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).push_back( prm::create_vis_parameter_list() );
        tParameterlist( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", "2_Channels_Temp.exo" ) );
        tParameterlist( 0 ).set( "Mesh_Type",  vis::VIS_Mesh_Type::STANDARD ) ;
        tParameterlist( 0 ).set( "Set_Names", "HMR_dummy_n_p34,HMR_dummy_c_p34,HMR_dummy_n_p46,HMR_dummy_c_p46,HMR_dummy_n_p42,HMR_dummy_c_p42" );
        tParameterlist( 0 ).set( "Field_Names", "VX,VY,P,TEMP" );
        tParameterlist( 0 ).set( "Field_Type", "NODAL,NODAL,NODAL,NODAL" );
        tParameterlist( 0 ).set( "IQI_Names", "IQIBulkVX,IQIBulkVY,IQIBulkP,IQIBulkTEMP" );
    }

    void
    MORISGENERALParameterList( Vector< Vector< ParameterList > >& tParameterlist )
    {
    }
    //------------------------------------------------------------------------------
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
