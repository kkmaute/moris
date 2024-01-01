/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Couette_Flow_Static.cpp
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
    std::string sFluid      = "HMR_dummy_n_p1,HMR_dummy_c_p1";
    std::string sFluidGhost = "ghost_p1";
    std::string sIn         = "iside_b0_1_b1_3";
    std::string sOut        = "iside_b0_1_b1_0";

    std::string sInterfaces = sIn + "," + sOut;

    std::string sAllDofTypes      = "VX,VY;P";
    std::string sVelocityDofTypes = "VX,VY";
    std::string sPressureDofTypes = "P";

    // Geometry Parameters
    Matrix< DDRMat > tCenterPoint = { { 0.0, 0.0 } }; /* Center point of the circle */
    moris::real      tRIn         = 1.0;              /* Inner circle radius (m) */
    moris::real      tROut        = 2.0;              /* Outer circle radius (m) */

    // Velocity Parameters
    moris::real tOmegaIn  = -5.0;
    moris::real tOmegaOut = 5.0;
    moris::real tKn       = 0.0;
    moris::real tb        = 1.0 / tRIn + 2.0 * tKn / std::pow( tRIn, 2.0 );
    moris::real td        = 1.0 / tROut - 2.0 * tKn / std::pow( tROut, 2.0 );
    moris::real tB        = ( tOmegaOut * tROut - tOmegaIn * tRIn * tROut / tRIn ) / ( td - tb * tROut / tRIn );
    moris::real tA        = ( tOmegaIn * tRIn - tb * tB ) / tRIn;

    // Constant function for properties
    void
    Func_Const(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    void
    ImposedVelocityFunc(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        // radius
        moris::real tx = aFIManager->get_IP_geometry_interpolator()->valx()( 0 );
        moris::real ty = aFIManager->get_IP_geometry_interpolator()->valx()( 1 );
        moris::real tR = std::pow( std::pow( tx, 2.0 ) + std::pow( ty, 2.0 ), 0.5 );

        // radial velocity
        moris::real tVTheta = tA * tR + tB / tR;

        // set size for aPropMatrix
        aPropMatrix.set_size( 2, 1, 0.0 );
        aPropMatrix( 0 ) = -tVTheta * ty / tR;
        aPropMatrix( 1 ) = tVTheta * tx / tR;
    }

    void
    AnalyticdVelocitydxFunc(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        // radius
        real tx = aFIManager->get_IP_geometry_interpolator()->valx()( 0 );
        real ty = aFIManager->get_IP_geometry_interpolator()->valx()( 1 );
        real tR = std::pow( std::pow( tx, 2.0 ) + std::pow( ty, 2.0 ), 0.5 );

        real tx2 = std::pow( tx, 2.0 );
        real ty2 = std::pow( ty, 2.0 );
        real tR2 = std::pow( tR, 2.0 );
        real tR4 = std::pow( tR, 4.0 );

        // set size for aPropMatrix
        aPropMatrix.set_size( 2, 2, 0.0 );
        aPropMatrix( 0, 0 ) = 2.0 * tB * tx * ty / tR4;
        aPropMatrix( 1, 0 ) = -tA - tB / tR2 + 2.0 * tB * ty2 / tR4;

        aPropMatrix( 0, 1 ) = tA + tB / tR2 - 2.0 * tB * tx2 / tR4;
        aPropMatrix( 1, 1 ) = -2.0 * tB * tx * ty / tR4;
    }

    bool
    Output_Criterion( moris::tsa::Time_Solver* aTimeSolver )
    {
        return true;
    }

    moris::real
    Func_Cylinder(
            const moris::Matrix< DDRMat >&     aCoordinates,
            const moris::Vector< moris::real* >& aGeometryParameters )
    {
        moris::real tR       = *( aGeometryParameters( 0 ) );
        moris::real tXCenter = *( aGeometryParameters( 1 ) );
        moris::real tYCenter = *( aGeometryParameters( 2 ) );

        moris::real aReturnValue = tR - std::pow( std::pow( aCoordinates( 0 ) - tXCenter, 2.0 ) + std::pow( aCoordinates( 1 ) - tYCenter, 2.0 ), 0.5 );

        return aReturnValue;
    }

    void
    OPTParameterList( moris::Vector< moris::Vector< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_opt_problem_parameter_list();

        tParameterlist( 0 )( 0 ).set( "is_optimization_problem", false );
    }

    void
    HMRParameterList( moris::Vector< moris::Vector< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_hmr_parameter_list();

        tParameterlist( 0 )( 0 ).set( "number_of_elements_per_dimension", "22,22" );
        tParameterlist( 0 )( 0 ).set( "domain_dimensions", "5,5" );
        tParameterlist( 0 )( 0 ).set( "domain_offset", "-2.5+0.25,-2.5+0.25" );
        tParameterlist( 0 )( 0 ).set( "domain_sidesets", "1,2,3,4" );
        tParameterlist( 0 )( 0 ).set( "lagrange_output_meshes", "0" );

        tParameterlist( 0 )( 0 ).set( "lagrange_orders", "1" );
        tParameterlist( 0 )( 0 ).set( "lagrange_pattern", "0" );
        tParameterlist( 0 )( 0 ).set( "bspline_orders", "1" );
        tParameterlist( 0 )( 0 ).set( "bspline_pattern", "0" );

        tParameterlist( 0 )( 0 ).set( "lagrange_to_bspline", "0" );

        tParameterlist( 0 )( 0 ).set( "truncate_bsplines", 1 );
        tParameterlist( 0 )( 0 ).set( "refinement_buffer", 3 );
        tParameterlist( 0 )( 0 ).set( "staircase_buffer", 3 );
        tParameterlist( 0 )( 0 ).set( "initial_refinement", "0" );
        tParameterlist( 0 )( 0 ).set( "initial_refinement_pattern", "0" );

        tParameterlist( 0 )( 0 ).set( "use_multigrid", 0 );
        tParameterlist( 0 )( 0 ).set( "severity_level", 0 );

        tParameterlist( 0 )( 0 ).set( "adaptive_refinement_level", 2 );
    }

    void
    XTKParameterList( moris::Vector< moris::Vector< ParameterList > >& tParameterlist )
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
        tParameterlist( 0 )( 0 ).set( "print_enriched_ig_mesh", true );
        tParameterlist( 0 )( 0 ).set( "exodus_output_XTK_ig_mesh", true );
    }

    void
    GENParameterList( moris::Vector< moris::Vector< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 3 );
        tParameterlist( 0 ).resize( 1 );

        // Main GEN parameter list
        tParameterlist( 0 )( 0 ) = prm::create_gen_parameter_list();

        // init geometry counter
        uint tGeoCounter = 0;

        // Geometry parameter lists
        tParameterlist( 1 ).push_back( prm::create_user_defined_geometry_parameter_list() );
        tParameterlist( 1 )( tGeoCounter ).set( "field_function_name", "Func_Cylinder" );
        tParameterlist( 1 )( tGeoCounter ).set( "constant_parameters", "1.0,0.0,0.0" );
        tGeoCounter++;

        tParameterlist( 1 ).push_back( prm::create_user_defined_geometry_parameter_list() );
        tParameterlist( 1 )( tGeoCounter ).set( "field_function_name", "Func_Cylinder" );
        tParameterlist( 1 )( tGeoCounter ).set( "constant_parameters", "2.0,0.0,0.0" );
    }

    void
    FEMParameterList( moris::Vector< moris::Vector< ParameterList > >& tParameterList )
    {
        // create a cell of cell of parameter list for fem
        tParameterList.resize( 8 );

        //------------------------------------------------------------------------------
        // fill the property part of the parameter list

        // init property counter
        uint tPropCounter = 0;

        // create parameter list for property 1
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropViscosity" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "1.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // create parameter list for property 2
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropDensity" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "1.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // create parameter list for property 3
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropDirichletInU" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "ImposedVelocityFunc" );
        tPropCounter++;

        // create parameter list for property 4
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropDirichletOutU" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "ImposedVelocityFunc" );
        tPropCounter++;

        // create parameter list for property 5
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropdUdx" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "AnalyticdVelocitydxFunc" );
        tPropCounter++;

        //------------------------------------------------------------------------------
        // fill the constitutive model part of the parameter list

        // init CM counter
        uint tCMCounter = 0;

        // create parameter list for constitutive model 1
        tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_name", "CMFluid" );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::FLUID_INCOMPRESSIBLE ) );
        tParameterList( 1 )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( "VX,VY;P", "Velocity,Pressure" ) );
        tParameterList( 1 )( tCMCounter ).set( "properties", "PropViscosity,Viscosity;PropDensity,Density" );
        tCMCounter++;

        //------------------------------------------------------------------------------
        // fill the stabilization parameter part of the parameter list

        // init SP counter
        uint tSPCounter = 0;

        // create parameter list for stabilization parameter 1
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", "SPIncFlow" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::INCOMPRESSIBLE_FLOW ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", "36.0" );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", "PropViscosity,Viscosity;PropDensity,Density" );
        tParameterList( 2 )( tSPCounter ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY;P", "Velocity,Pressure" ) );
        tSPCounter++;

        // create parameter list for stabilization parameter 2
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", "SPDirichletNitscheU" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::VELOCITY_DIRICHLET_NITSCHE ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", "100.0/1.0" );
        tParameterList( 2 )( tSPCounter ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", "PropViscosity,Viscosity;PropDensity,Density" );
        tSPCounter++;

        // create parameter list for stabilization parameter 3
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", "SPGPViscosity" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::VISCOUS_GHOST ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", "0.05" );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", "PropViscosity,Viscosity" );
        tSPCounter++;

        // create parameter list for stabilization parameter 4
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", "SPGPVelocity" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::CONVECTIVE_GHOST ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", "0.05" );
        tParameterList( 2 )( tSPCounter ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", "PropDensity,Density" );
        tSPCounter++;

        // create parameter list for stabilization parameter 5
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", "SPGPPressure" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::PRESSURE_GHOST ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", "0.005/1.0" );
        tParameterList( 2 )( tSPCounter ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", "PropViscosity,Viscosity;PropDensity,Density" );
        tSPCounter++;

        //------------------------------------------------------------------------------
        // fill the IWG part of the parameter list

        // init IWG counter
        uint tIWGCounter = 0;

        // create parameter list for IWG 1
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGVelocityBulk" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_BULK ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "VX,VY" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "VX,VY;P" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPIncFlow,IncompressibleFlow" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", sFluid );
        tIWGCounter++;

        // create parameter list for IWG 2
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGPressureBulk" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_BULK ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "P" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "VX,VY;P" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPIncFlow,IncompressibleFlow" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", sFluid );
        tIWGCounter++;

        // create parameter list for IWG 3
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGInVelocity" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_SYMMETRIC_NITSCHE ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "VX,VY" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "VX,VY;P" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties", "PropDirichletInU,Dirichlet" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPDirichletNitscheU,DirichletNitsche" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", sIn );
        tIWGCounter++;

        // create parameter list for IWG 4
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGInPressure" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_SYMMETRIC_NITSCHE ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "P" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "VX,VY;P" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties", "PropDirichletInU,Dirichlet" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", sIn );
        tIWGCounter++;

        // create parameter list for IWG 5
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGOutVelocity" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_SYMMETRIC_NITSCHE ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "VX,VY" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "VX,VY;P" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties", "PropDirichletOutU,Dirichlet" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPDirichletNitscheU,DirichletNitsche" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", sOut );
        tIWGCounter++;

        // create parameter list for IWG 6
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGOutPressure" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_SYMMETRIC_NITSCHE ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "P" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "VX,VY;P" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties", "PropDirichletOutU,Dirichlet" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", sOut );
        tIWGCounter++;

        // create parameter list for IWG 7
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGGPViscous" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::GHOST_NORMAL_FIELD ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "VX,VY" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "VX,VY;P" );
        tParameterList( 3 )( tIWGCounter ).set( "follower_dof_dependencies", "VX,VY;P" );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPGPViscosity,GhostSP" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", sFluidGhost );
        tIWGCounter++;

        // create parameter list for IWG 8
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGGPConvective" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::GHOST_NORMAL_FIELD ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "VX,VY" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "VX,VY;P" );
        tParameterList( 3 )( tIWGCounter ).set( "follower_dof_dependencies", "VX,VY;P" );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPGPVelocity,GhostSP" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", sFluidGhost );
        tIWGCounter++;

        // create parameter list for IWG 9
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGGPPressure" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::GHOST_NORMAL_FIELD ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "P" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "VX,VY;P" );
        tParameterList( 3 )( tIWGCounter ).set( "follower_dof_dependencies", "VX,VY;P" );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPGPPressure,GhostSP" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", sFluidGhost );
        tIWGCounter++;

        //------------------------------------------------------------------------------
        // fill the IQI part of the parameter list

        // init IQI counter
        uint tIQICounter = 0;

        // create parameter list for IQI 1
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkVX" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::DOF ) );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "VX,VY" );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "VX,VY" );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 0 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", sFluid );
        tIQICounter++;

        // create parameter list for IQI 2
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkVY" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::DOF ) );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "VX,VY" );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "VX,VY" );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 1 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", sFluid );
        tIQICounter++;

        // create parameter list for IQI 3
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkP" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::DOF ) );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "P" );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "P" );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 0 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", sFluid );
        tIQICounter++;

        // create parameter list for IQI 4
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkL2Error" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::L2_ERROR_ANALYTIC ) );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "VX,VY" );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "VX,VY" );
        tParameterList( 4 )( tIQICounter ).set( "leader_properties", "PropDirichletInU,L2Check" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", sFluid );
        tIQICounter++;

        // create parameter list for IQI 5
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkH1Error" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::H1_ERROR_ANALYTIC ) );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "VX,VY" );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "VX,VY" );
        tParameterList( 4 )( tIQICounter ).set( "leader_properties", "PropdUdx,H1Check" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", sFluid );
        tIQICounter++;

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        tParameterList( 5 ).resize( 1 );
        tParameterList( 5 )( 0 ) = prm::create_computation_parameter_list();
    }

    void
    SOLParameterList( moris::Vector< moris::Vector< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 8 );
        for ( uint Ik = 0; Ik < 8; Ik++ )
        {
            tParameterlist( Ik ).resize( 1 );
        }

        tParameterlist( 0 )( 0 ) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL );

        tParameterlist( 1 )( 0 ) = moris::prm::create_linear_solver_parameter_list();

        tParameterlist( 2 )( 0 ) = moris::prm::create_nonlinear_algorithm_parameter_list();
        tParameterlist( 2 )( 0 ).set( "NLA_rel_res_norm_drop", 1e-06 );
        tParameterlist( 2 )( 0 ).set( "NLA_relaxation_parameter", 1.0 );
        tParameterlist( 2 )( 0 ).set( "NLA_max_iter", 20 );

        tParameterlist( 3 )( 0 ) = moris::prm::create_nonlinear_solver_parameter_list();
        tParameterlist( 3 )( 0 ).set( "NLA_DofTypes", "VX,VY;P" );

        tParameterlist( 4 )( 0 ) = moris::prm::create_time_solver_algorithm_parameter_list();

        tParameterlist( 5 )( 0 ) = moris::prm::create_time_solver_parameter_list();
        tParameterlist( 5 )( 0 ).set( "TSA_DofTypes", "VX,VY;P" );
        tParameterlist( 5 )( 0 ).set( "TSA_Initialize_Sol_Vec", "VX,1E-4;VY,1E-4;P,0.0" );
        tParameterlist( 5 )( 0 ).set( "TSA_Output_Indices", "0" );
        tParameterlist( 5 )( 0 ).set( "TSA_Output_Criteria", "Output_Criterion" );

        tParameterlist( 6 )( 0 ) = moris::prm::create_solver_warehouse_parameterlist();

        tParameterlist( 7 )( 0 ) = moris::prm::create_preconditioner_parameter_list(sol::PreconditionerType::NONE);
    }

    void
    MSIParameterList( moris::Vector< moris::Vector< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_msi_parameter_list();
    }

    void
    VISParameterList( moris::Vector< moris::Vector< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_vis_parameter_list();
        tParameterlist( 0 )( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", "Couette_Flow_Static.exo" ) );
        tParameterlist( 0 )( 0 ).set( "Mesh_Type", static_cast< uint >( vis::VIS_Mesh_Type::STANDARD ) );
        tParameterlist( 0 )( 0 ).set( "Set_Names", sFluid + "," + sInterfaces );
        tParameterlist( 0 )( 0 ).set( "Field_Names", "VX,VY,P,L2_ERROR_ANALYTIC,H1_ERROR_ANALYTIC" );
        tParameterlist( 0 )( 0 ).set( "Field_Type", "NODAL,NODAL,NODAL,NODAL,NODAL" );
        tParameterlist( 0 )( 0 ).set( "IQI_Names", "IQIBulkVX,IQIBulkVY,IQIBulkP,IQIBulkL2Error,IQIBulkH1Error" );
    }

    void
    MORISGENERALParameterList( moris::Vector< moris::Vector< ParameterList > >& tParameterlist )
    {
    }

    //------------------------------------------------------------------------------
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif

