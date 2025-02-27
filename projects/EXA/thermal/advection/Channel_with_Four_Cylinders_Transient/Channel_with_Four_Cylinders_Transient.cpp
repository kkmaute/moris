/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Channel_with_Four_Cylinders_Transient.cpp
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

static const bool isGhost = true;

#ifdef __cplusplus
extern "C" {
#endif
//------------------------------------------------------------------------------
namespace moris
{
    // Minimum value to be returned by level set function
    moris::real tMinLSvalue = 1.0e-6;

    // Geometry Parameters
    moris::real tPlaneBottom     = 0.0;  /* y bottom plane (m) */
    moris::real tPlaneTop        = 0.41; /* y top plane    (m) */
    moris::real tPlaneLeft       = 0.0;  /* x left plane   (m) */
    moris::real tPlaneRight      = 2.2;  /* x right plane  (m) */
    moris::real tCylinderCenterX = 0.2;
    moris::real tCylinderCenterY = 0.2;
    moris::real tCylinderRadius  = 0.05;
    moris::real tCylinderOffset  = 0.10;

    moris::sint sNumStep = 10;
    moris::real sMaxTime = 10.0;
    moris::real sNumRamp = 5.0;

    // Constant function for properties
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
        aPropMatrix.set_size( 2, 1, 0.0 );

        real tT = aFIManager->get_IP_geometry_interpolator()->valt()( 0 );
        real tY = aFIManager->get_IP_geometry_interpolator()->valx()( 1 );

        real tMagnitude = std::min( tT / ( sMaxTime / sNumStep ) / sNumRamp, 1.0 );

        aPropMatrix( 0 ) = 4.0 * tMagnitude * aParameters( 0 )( 0 ) * tY * ( 0.41 - tY ) / ( std::pow( 0.41, 2.0 ) );
    }

    void
    Func_Initial_Condition(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix.set_size( 2, 1, 0.0 );
        aPropMatrix( 0 ) = 1e-5;
    }

    bool
    Output_Criterion( moris::tsa::Time_Solver* aTimeSolver )
    {
        return true;
    }

    moris::real Func_Bottom_Plane(
            const moris::Matrix< DDRMat >& aCoordinates,
            const Vector< real >&     aGeometryParameters )
    {
        moris::real tValue = aCoordinates( 1 ) - tPlaneBottom;

        return std::abs( tValue ) < tMinLSvalue ? tMinLSvalue : tValue;
    }

    moris::real Func_Top_Plane(
            const moris::Matrix< DDRMat >& aCoordinates,
            const Vector< real >&     aGeometryParameters )
    {
        moris::real tValue = aCoordinates( 1 ) - tPlaneTop;

        return std::abs( tValue ) < tMinLSvalue ? tMinLSvalue : tValue;
    }

    moris::real Func_Left_Plane(
            const moris::Matrix< DDRMat >& aCoordinates,
            const Vector< real >&     aGeometryParameters )
    {
        moris::real tValue = aCoordinates( 0 ) - tPlaneLeft;

        return std::abs( tValue ) < tMinLSvalue ? tMinLSvalue : tValue;
    }

    moris::real Func_Right_Plane(
            const moris::Matrix< DDRMat >& aCoordinates,
            const Vector< real >&     aGeometryParameters )
    {
        moris::real tValue = aCoordinates( 0 ) - tPlaneRight;

        return std::abs( tValue ) < tMinLSvalue ? tMinLSvalue : tValue;
    }

    moris::real Func_Cylinder(
            const moris::Matrix< DDRMat >& aCoordinates,
            const Vector< real >&     aGeometryParameters )
    {
        moris::real tValue = tCylinderRadius - std::pow( std::pow( aCoordinates( 0 ) - 0.4, 2.0 ) + std::pow( aCoordinates( 1 ) - ( 0.2 - tCylinderOffset ), 2.0 ), 0.5 );

        return std::abs( tValue ) < tMinLSvalue ? tMinLSvalue : tValue;
    }

    moris::real Func_Cylinder2(
            const moris::Matrix< DDRMat >& aCoordinates,
            const Vector< real >&     aGeometryParameters )
    {
        moris::real tValue = tCylinderRadius - std::pow( std::pow( aCoordinates( 0 ) - 0.8, 2.0 ) + std::pow( aCoordinates( 1 ) - ( 0.2 + tCylinderOffset ), 2.0 ), 0.5 );

        return std::abs( tValue ) < tMinLSvalue ? tMinLSvalue : tValue;
    }

    moris::real Func_Cylinder3(
            const moris::Matrix< DDRMat >& aCoordinates,
            const Vector< real >&     aGeometryParameters )
    {
        moris::real tValue = tCylinderRadius - std::pow( std::pow( aCoordinates( 0 ) - 1.2, 2.0 ) + std::pow( aCoordinates( 1 ) - ( 0.2 - tCylinderOffset ), 2.0 ), 0.5 );

        return std::abs( tValue ) < tMinLSvalue ? tMinLSvalue : tValue;
    }

    moris::real Func_Cylinder4(
            const moris::Matrix< DDRMat >& aCoordinates,
            const Vector< real >&     aGeometryParameters )
    {
        moris::real tValue = tCylinderRadius - std::pow( std::pow( aCoordinates( 0 ) - 1.6, 2.0 ) + std::pow( aCoordinates( 1 ) - ( 0.2 + tCylinderOffset ), 2.0 ), 0.5 );

        return std::abs( tValue ) < tMinLSvalue ? tMinLSvalue : tValue;
    }

    void
    OPTParameterList( Vector< Vector< ParameterList > >& aParameterLists )
    {
        aParameterLists( 0 ).push_back( prm::create_opt_problem_parameter_list() );

        aParameterLists.set( "is_optimization_problem", false );
    }

    void
    HMRParameterList( Vector< Vector< ParameterList > >& aParameterLists )
    {
        aParameterLists( 0 ).push_back( prm::create_hmr_parameter_list() );

        aParameterLists.set( "number_of_elements_per_dimension", 22, 11 );
        aParameterLists.set( "domain_dimensions", 4.0, 2.0 );
        aParameterLists.set( "domain_offset", -1.24, -0.86 );
        aParameterLists.set( "lagrange_output_meshes", "0" );

        aParameterLists.set( "lagrange_orders", "1" );
        aParameterLists.set( "lagrange_pattern", "0" );
        aParameterLists.set( "bspline_orders", "1" );
        aParameterLists.set( "bspline_pattern", "0" );

        aParameterLists.set( "lagrange_to_bspline", "0" );

        aParameterLists.set( "refinement_buffer", 3 );
        aParameterLists.set( "staircase_buffer", 5 );
        aParameterLists.set( "initial_refinement", 4 );

        aParameterLists.set( "use_number_aura", false );

        aParameterLists.set( "adaptive_refinement_level", 1 );
    }

    void
    XTKParameterList( Vector< Vector< ParameterList > >& aParameterLists )
    {
        aParameterLists( 0 ).push_back( prm::create_xtk_parameter_list() );
        aParameterLists.set( "decompose", true );
        aParameterLists.set( "decomposition_type", "conformal" );
        aParameterLists.set( "enrich", true );
        aParameterLists.set( "basis_rank", "bspline" );
        aParameterLists.set( "enrich_mesh_indices", "0" );
        aParameterLists.set( "ghost_stab", isGhost );
        aParameterLists.set( "multigrid", false );
        aParameterLists.set( "verbose", true );
        aParameterLists.set( "print_enriched_ig_mesh", false );
        aParameterLists.set( "exodus_output_XTK_ig_mesh", false );
    }

    void
    GENParameterList( Vector< Vector< ParameterList > >& aParameterLists )
    {

        aParameterLists( 0 ).push_back( prm::create_gen_parameter_list() );

        // Geometry parameter lists
        aParameterLists( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists.set( "field_function_name", "Func_Bottom_Plane" );

        aParameterLists( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists.set( "field_function_name", "Func_Top_Plane" );

        aParameterLists( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists.set( "field_function_name", "Func_Left_Plane" );

        aParameterLists( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists.set( "field_function_name", "Func_Right_Plane" );

        aParameterLists( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists.set( "field_function_name", "Func_Cylinder" );

        aParameterLists( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists.set( "field_function_name", "Func_Cylinder2" );

        aParameterLists( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists.set( "field_function_name", "Func_Cylinder3" );

        aParameterLists( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists.set( "field_function_name", "Func_Cylinder4" );
    }

    void
    FEMParameterList( Vector< Vector< ParameterList > >& aParameterLists )
    {
        // create a cell of cell of parameter list for fem

        //------------------------------------------------------------------------------
        // fill the property part of the parameter list

        // create parameter list for property 1
        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropViscosity" );
        aParameterLists.set( "function_parameters", "0.001" );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 2
        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropDensity" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 2
        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropCapacity" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 3
        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropDirichletInletU" );
        aParameterLists.set( "function_parameters", "1.5" );
        aParameterLists.set( "value_function", "Func_Inlet_U" );

        // create parameter list for property 4
        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropDirichletZeroU" );
        aParameterLists.set( "function_parameters", "0.0;0.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 3
        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropConductivity" );
        aParameterLists.set( "function_parameters", "0.0005" );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 5
        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropInletTemp" );
        aParameterLists.set( "function_parameters", "0.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 6
        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropSideFlux" );
        aParameterLists.set( "function_parameters", "20.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 6
        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropWeightCurrent" );
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 7
        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropWeightPrevious" );
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 8
        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropInitialConditionU" );
        aParameterLists.set( "function_parameters", "1.5" );
        aParameterLists.set( "value_function", "Func_Initial_Condition" );

        // create parameter list for property 8
        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropInitialConditionTEMP" );
        aParameterLists.set( "function_parameters", "0.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        //------------------------------------------------------------------------------
        // fill the constitutive model part of the parameter list

        // create parameter list for constitutive model 1
        aParameterLists( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        aParameterLists.set( "constitutive_name", "CMFluid" );
        aParameterLists.set( "constitutive_type",  fem::Constitutive_Type::FLUID_INCOMPRESSIBLE ) ;
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "VX,VY;P", "Velocity,Pressure" ) );
        aParameterLists.set( "properties",
                "PropViscosity,Viscosity;"
                "PropDensity  ,Density" );

        // create parameter list for constitutive model 2
        aParameterLists( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        aParameterLists.set( "constitutive_name", "CMDiffusion" );
        aParameterLists.set( "constitutive_type",  fem::Constitutive_Type::DIFF_LIN_ISO ) ;
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists.set( "properties",
                "PropConductivity,Conductivity;"
                "PropDensity     ,Density;"
                "PropCapacity    ,HeatCapacity" );

        //------------------------------------------------------------------------------
        // fill the stabilization parameter part of the parameter list

        // create parameter list for stabilization parameter 1
        aParameterLists( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists.set( "stabilization_name", "SPIncFlow" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::INCOMPRESSIBLE_FLOW ) ;
        aParameterLists.set( "function_parameters", "36.0" );
        aParameterLists.set( "leader_properties", "PropViscosity,Viscosity;PropDensity,Density" );
        aParameterLists.set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY;P", "Velocity,Pressure" ) );

        // create parameter list for stabilization parameter 2
        aParameterLists( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists.set( "stabilization_name", "SPSUPGTemp" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::SUPG_ADVECTION ) ;
        aParameterLists.set( "leader_properties", "PropConductivity,Conductivity" );
        aParameterLists.set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );

        // create parameter list for stabilization parameter 2
        aParameterLists( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists.set( "stabilization_name", "SPDirichletNitscheU" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::VELOCITY_DIRICHLET_NITSCHE ) ;
        aParameterLists.set( "function_parameters", "100.0/1.0" );
        aParameterLists.set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
        aParameterLists.set( "leader_properties", "PropViscosity,Viscosity;PropDensity,Density" );

        // create parameter list for stabilization parameter 2
        aParameterLists( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists.set( "stabilization_name", "SPDirichletNitscheT" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::DIRICHLET_NITSCHE ) ;
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "leader_properties", "PropConductivity,Material" );

        // create parameter list for stabilization parameter 3
        aParameterLists( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists.set( "stabilization_name", "SPGPViscosity" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::VISCOUS_GHOST ) ;
        aParameterLists.set( "function_parameters", "0.05" );
        aParameterLists.set( "leader_properties", "PropViscosity,Viscosity" );

        // create parameter list for stabilization parameter 4
        aParameterLists( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists.set( "stabilization_name", "SPGPVelocity" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::CONVECTIVE_GHOST ) ;
        aParameterLists.set( "function_parameters", "0.05" );
        aParameterLists.set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
        aParameterLists.set( "leader_properties", "PropDensity,Density" );

        // create parameter list for stabilization parameter 5
        aParameterLists( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists.set( "stabilization_name", "SPGPPressure" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::PRESSURE_GHOST ) ;
        aParameterLists.set( "function_parameters", "0.005/1.0" );
        aParameterLists.set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
        aParameterLists.set( "leader_properties", "PropViscosity,Viscosity;PropDensity,Density" );

        // create parameter list for stabilization parameter 8
        aParameterLists( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists.set( "stabilization_name", "SPGPTemp" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::GHOST_DISPL ) ;
        aParameterLists.set( "function_parameters", "0.005" );
        aParameterLists.set( "leader_properties", "PropConductivity,Material" );

        //------------------------------------------------------------------------------
        // fill the IWG part of the parameter list

        // create parameter list for IWG 1
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGVelocityBulk" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_BULK ) ;
        aParameterLists.set( "dof_residual", "VX,VY" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists.set( "stabilization_parameters", "SPIncFlow,IncompressibleFlow" );
        aParameterLists.set( "mesh_set_names", "HMR_dummy_c_p160,HMR_dummy_n_p160" );

        // create parameter list for IWG 2
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGPressureBulk" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_BULK ) ;
        aParameterLists.set( "dof_residual", "P" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists.set( "stabilization_parameters", "SPIncFlow,IncompressibleFlow" );
        aParameterLists.set( "mesh_set_names", "HMR_dummy_c_p160,HMR_dummy_n_p160" );

        // create parameter list for IWG 3
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGDiffusionBulk" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_BULK ) ;
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusion,Diffusion" );
        aParameterLists.set( "mesh_set_names", "HMR_dummy_c_p160,HMR_dummy_n_p160" );

        // create parameter list for IWG 4
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGAdvectionBulk" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::ADVECTION_BULK ) ;
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusion,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPSUPGTemp,SUPG" );
        aParameterLists.set( "mesh_set_names", "HMR_dummy_c_p160,HMR_dummy_n_p160" );

        // create parameter list for IWG 3
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGInletVelocity" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        aParameterLists.set( "dof_residual", "VX,VY" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists.set( "leader_properties", "PropDirichletInletU,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists.set( "stabilization_parameters", "SPDirichletNitscheU,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", "iside_b0_160_b1_128" );

        // create parameter list for IWG 4
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGInletPressure" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        aParameterLists.set( "dof_residual", "P" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists.set( "leader_properties", "PropDirichletInletU,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists.set( "mesh_set_names", "iside_b0_160_b1_128" );

        // create parameter list for IWG 11
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGInletTemp" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists.set( "leader_properties", "PropInletTemp,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusion,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPDirichletNitscheT,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", "iside_b0_160_b1_128" );

        // create parameter list for IWG 11
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGCylinderFluxTemp" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_NEUMANN ) ;
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists.set( "leader_properties", "PropSideFlux,Neumann" );
        aParameterLists.set( "mesh_set_names", "iside_b0_160_b1_168,iside_b0_160_b1_164,iside_b0_160_b1_162,iside_b0_160_b1_161" );

        // create parameter list for IWG 5
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGZeroVelocity" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        aParameterLists.set( "dof_residual", "VX,VY" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists.set( "leader_properties", "PropDirichletZeroU,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists.set( "stabilization_parameters", "SPDirichletNitscheU,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", "iside_b0_160_b1_32,iside_b0_160_b1_224,iside_b0_160_b1_168,iside_b0_160_b1_164,iside_b0_160_b1_162,iside_b0_160_b1_161" );

        // create parameter list for IWG 6
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGZeroPressure" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        aParameterLists.set( "dof_residual", "P" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists.set( "leader_properties", "PropDirichletZeroU,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists.set( "mesh_set_names", "iside_b0_160_b1_32,iside_b0_160_b1_224,iside_b0_160_b1_168,iside_b0_160_b1_164,iside_b0_160_b1_162,iside_b0_160_b1_161" );

        // create parameter list for IWG 5
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGZeroVelocity" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        aParameterLists.set( "dof_residual", "VX,VY" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists.set( "leader_properties", "PropDirichletZeroU,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists.set( "stabilization_parameters", "SPDirichletNitscheU,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", "iside_b0_160_b1_32,iside_b0_160_b1_224,iside_b0_160_b1_168,iside_b0_160_b1_164,iside_b0_160_b1_162,iside_b0_160_b1_161" );

        // create parameter list for IWG 6
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGZeroPressure" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        aParameterLists.set( "dof_residual", "P" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists.set( "leader_properties", "PropDirichletZeroU,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists.set( "mesh_set_names", "iside_b0_160_b1_32,iside_b0_160_b1_224,iside_b0_160_b1_168,iside_b0_160_b1_164,iside_b0_160_b1_162,iside_b0_160_b1_161" );

        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGTimeContinuityVelocity" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::TIME_CONTINUITY_DOF ) ;
        aParameterLists.set( "dof_residual", "VX,VY" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists.set( "leader_properties",
                "PropWeightCurrent,WeightCurrent;"
                "PropWeightPrevious,WeightPrevious;"
                "PropInitialConditionU,InitialCondition" );
        aParameterLists.set( "mesh_set_names", "HMR_dummy_c_p160,HMR_dummy_n_p160" );
        aParameterLists.set( "time_continuity", true );

        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGTimeContinuityTemp" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::TIME_CONTINUITY_DOF ) ;
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
        aParameterLists.set( "leader_properties",
                "PropWeightCurrent,WeightCurrent;"
                "PropWeightPrevious,WeightPrevious;"
                "PropInitialConditionTEMP,InitialCondition" );
        aParameterLists.set( "mesh_set_names", "HMR_dummy_c_p160,HMR_dummy_n_p160" );
        aParameterLists.set( "time_continuity", true );

        if ( isGhost )
        {
            // create parameter list for IWG 7
            aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
            aParameterLists.set( "IWG_name", "IWGGPViscous" );
            aParameterLists.set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
            aParameterLists.set( "dof_residual", "VX,VY" );
            aParameterLists.set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
            aParameterLists.set( "follower_dof_dependencies", "VX,VY;P;TEMP" );
            aParameterLists.set( "stabilization_parameters", "SPGPViscosity,GhostSP" );
            aParameterLists.set( "mesh_set_names", "ghost_p160" );

            // create parameter list for IWG 8
            aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
            aParameterLists.set( "IWG_name", "IWGGPConvective" );
            aParameterLists.set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
            aParameterLists.set( "dof_residual", "VX,VY" );
            aParameterLists.set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
            aParameterLists.set( "follower_dof_dependencies", "VX,VY;P;TEMP" );
            aParameterLists.set( "stabilization_parameters", "SPGPVelocity,GhostSP" );
            aParameterLists.set( "mesh_set_names", "ghost_p160" );

            // create parameter list for IWG 9
            aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
            aParameterLists.set( "IWG_name", "IWGGPPressure" );
            aParameterLists.set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
            aParameterLists.set( "dof_residual", "P" );
            aParameterLists.set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
            aParameterLists.set( "follower_dof_dependencies", "VX,VY;P;TEMP" );
            aParameterLists.set( "stabilization_parameters", "SPGPPressure,GhostSP" );
            aParameterLists.set( "mesh_set_names", "ghost_p160" );

            // create parameter list for IWG 16
            aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
            aParameterLists.set( "IWG_name", "IWGGPTemp" );
            aParameterLists.set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
            aParameterLists.set( "dof_residual", "TEMP" );
            aParameterLists.set( "leader_dof_dependencies", "VX,VY;P;TEMP" );
            aParameterLists.set( "follower_dof_dependencies", "VX,VY;P;TEMP" );
            aParameterLists.set( "stabilization_parameters", "SPGPTemp,GhostSP" );
            aParameterLists.set( "mesh_set_names", "ghost_p160" );
            }

        //------------------------------------------------------------------------------
        // fill the IQI part of the parameter list

        // create parameter list for IQI 1
        aParameterLists( 4 ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists.set( "IQI_name", "IQIBulkVX" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists.set( "dof_quantity", "VX,VY" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "mesh_set_names", "HMR_dummy_n_p160,HMR_dummy_c_p160" );

        // create parameter list for IQI 2
        aParameterLists( 4 ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists.set( "IQI_name", "IQIBulkVY" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists.set( "dof_quantity", "VX,VY" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY" );
        aParameterLists.set( "vectorial_field_index", 1 );
        aParameterLists.set( "mesh_set_names", "HMR_dummy_n_p160,HMR_dummy_c_p160" );

        // create parameter list for IQI 3
        aParameterLists( 4 ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists.set( "IQI_name", "IQIBulkP" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists.set( "dof_quantity", "P" );
        aParameterLists.set( "leader_dof_dependencies", "P" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "mesh_set_names", "HMR_dummy_n_p160,HMR_dummy_c_p160" );

        // create parameter list for IQI 3
        aParameterLists( 4 ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists.set( "IQI_name", "IQIBulkTEMP" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists.set( "dof_quantity", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "mesh_set_names", "HMR_dummy_n_p160,HMR_dummy_c_p160" );

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        aParameterLists( 5 ).push_back( prm::create_computation_parameter_list() );
    }

    void
    SOLParameterList( Vector< Vector< ParameterList > >& aParameterLists )
    {

        aParameterLists( 0 ).push_back( add_parameter_list( sol::SolverType::BELOS_IMPL ) );
        aParameterLists.set( "ifpack_prec_type", "ILU" );

        aParameterLists( 1 ).push_back( moris::prm::create_linear_solver_parameter_list() );

                aParameterLists( 2 ).push_back( moris::prm::create_nonlinear_algorithm_parameter_list() );    // nonlinear algorithm index 0
        aParameterLists.set( "NLA_Solver_Implementation",  moris::NLA::NonlinearSolverType::NEWTON_SOLVER ) ;
        aParameterLists.set( "NLA_rel_res_norm_drop", 1e-04 );
        aParameterLists.set( "NLA_relaxation_parameter", 1.0 );
        aParameterLists.set( "NLA_max_iter", 10 );

        aParameterLists( 2 ).push_back( moris::prm::create_nonlinear_algorithm_parameter_list() );    // nonlinear algorithm index 1
        aParameterLists.set( "NLA_Solver_Implementation",  moris::NLA::NonlinearSolverType::NLBGS_SOLVER ) ;
        aParameterLists.set( "NLA_rel_res_norm_drop", 1.0 );
        aParameterLists.set( "NLA_max_iter", 1 );

        aParameterLists( 3 ).push_back( moris::prm::create_nonlinear_solver_parameter_list() );    // nonlinear solver index 0
        aParameterLists.set( "NLA_Nonlinear_solver_algorithms", "0" );             // set nonlinear algorithm with index 0
        aParameterLists.set( "NLA_Solver_Implementation",  moris::NLA::NonlinearSolverType::NEWTON_SOLVER ) ;
        aParameterLists.set( "NLA_DofTypes", "VX,VY,P" );

        aParameterLists( 3 ).push_back( moris::prm::create_nonlinear_solver_parameter_list() );    // nonlinear solver index 1
        aParameterLists.set( "NLA_Nonlinear_solver_algorithms", "0" );             // set nonlinear algorithm with index 0
        aParameterLists.set( "NLA_Solver_Implementation",  moris::NLA::NonlinearSolverType::NEWTON_SOLVER ) ;
        aParameterLists.set( "NLA_DofTypes", "TEMP" );

        aParameterLists( 3 ).push_back( moris::prm::create_nonlinear_solver_parameter_list() );    // nonlinear solver index 2
        aParameterLists.set( "NLA_Nonlinear_solver_algorithms", "1" );             // set nonlinear algorithm with index 1.
        aParameterLists.set( "NLA_Solver_Implementation",  moris::NLA::NonlinearSolverType::NLBGS_SOLVER ) ;
        aParameterLists.set( "NLA_Sub_Nonlinear_Solver", "0,1" );    // set sub nonlinear solvers with index 0 and 1
        aParameterLists.set( "NLA_DofTypes", "VX,VY,P;TEMP" );

        aParameterLists( 4 ).push_back( moris::prm::create_time_solver_algorithm_parameter_list() );
        aParameterLists.set( "TSA_Nonlinear_Solver", 2 );    // set nonlinear solver with index 2
        aParameterLists.set( "TSA_Num_Time_Steps", sNumStep );
        aParameterLists.set( "TSA_Time_Frame", sMaxTime );

        aParameterLists( 5 ).push_back( moris::prm::create_time_solver_parameter_list() );
        aParameterLists.set( "TSA_DofTypes", "VX,VY;P;TEMP" );
        aParameterLists.set( "TSA_Initialize_Sol_Vec", "VX,1.5;VY,0.0;P,0.0;TEMP,0.0" );
        aParameterLists.set( "TSA_Output_Indices", "0" );
        aParameterLists.set( "TSA_Output_Criteria", "Output_Criterion" );

        aParameterLists( 6 ).push_back( moris::prm::create_solver_warehouse_parameterlist() );

        aParameterLists( 7 ).push_back(  sol::PreconditionerType::NONE );
    }

    void
    MSIParameterList( Vector< Vector< ParameterList > >& aParameterLists )
    {
        aParameterLists( 0 ).push_back( prm::create_msi_parameter_list() );
    }

    void
    VISParameterList( Vector< Vector< ParameterList > >& aParameterLists )
    {
        aParameterLists( 0 ).push_back( prm::create_vis_parameter_list() );
        aParameterLists.set( "File_Name", std::pair< std::string, std::string >( "./", "Channel_with_Four_Cylinders_Transient.exo" ) );
        aParameterLists.set( "Mesh_Type",  vis::VIS_Mesh_Type::STANDARD ) ;
        aParameterLists.set( "Set_Names", "HMR_dummy_n_p160,HMR_dummy_c_p160" );
        aParameterLists.set( "Field_Names", "VX,VY,P,TEMP" );
        aParameterLists.set( "Field_Type", "NODAL,NODAL,NODAL,NODAL" );
        aParameterLists.set( "IQI_Names", "IQIBulkVX,IQIBulkVY,IQIBulkP,IQIBulkTEMP" );
        aParameterLists.set( "Save_Frequency", 1 );
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
