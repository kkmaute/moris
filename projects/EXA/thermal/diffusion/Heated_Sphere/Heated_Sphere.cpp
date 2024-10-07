/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Heated_Sphere.cpp
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
#include "cl_HMR_Element.hpp"
#include "fn_equal_to.hpp"

#include "AztecOO.h"

//---------------------------------------------------------------

// global variable for interpolation order
extern uint gInterpolationOrder;

//---------------------------------------------------------------

#ifdef __cplusplus
extern "C" {
#endif
//------------------------------------------------------------------------------
namespace moris
{
    /* ------------------------------------------------------------------------ */
    // Mesh Set Information

    std::string tInnerPhase = "HMR_dummy_n_p1,HMR_dummy_c_p1";
    std::string tOuterPhase = "HMR_dummy_n_p0,HMR_dummy_c_p0";

    std::string tInterface = "dbl_iside_p0_1_p1_0";

    std::string tSphereSurface        = "iside_b0_1_b1_0";
    std::string tBlockInnerSurface    = "iside_b0_0_b1_1";
    std::string tInterfaceSingleSides = tSphereSurface + "," + tBlockInnerSurface;

    std::string tOuterSurface = "SideSet_1_n_p0,SideSet_2_n_p0,SideSet_3_n_p0,SideSet_4_n_p0,SideSet_5_n_p0,SideSet_6_n_p0";

    std::string tInnerPhaseGhost = "ghost_p1";
    std::string tOuterPhaseGhost = "ghost_p0";

    std::string tTotalDomain   = tInnerPhase + "," + tOuterPhase;
    std::string tAllInterfaces = tOuterSurface + "," + tInterface;

    /* ------------------------------------------------------------------------ */
    // geometry parameters

    // general
    moris::real tCenterX = 0.0;
    moris::real tCenterY = 0.0;
    moris::real tCenterZ = 0.0;
    moris::real tRadius  = 0.5;

    /* ------------------------------------------------------------------------ */
    // material parameters

    // capacity
    std::string tCapacityInner = "0.0";
    std::string tCapacityOuter = "0.0";

    // density
    std::string tDensityInner = "0.0";
    std::string tDensityOuter = "0.0";

    // conductivity
    std::string tConductivityInner = "1.0";
    std::string tConductivityOuter = "0.125";

    // body flux
    std::string tHeatLoadInner = "1.0";
    std::string tHeatLoadOuter = "0.0";

    /* ------------------------------------------------------------------------ */
    // HMR parameters

    std::string tNumElemsPerDim = "8, 8, 8";
    std::string tDomainDims     = "2.0, 2.0, 2.0";
    std::string tDomainOffset   = "-1.0, -1.0, -1.0";
    std::string tDomainSidesets = "1,2,3,4,5,6";

    int tRefineBuffer = 1;

    /* ------------------------------------------------------------------------ */
    // Solver config

    moris::real tNLA_rel_res_norm_drop    = 1.0e-08;
    moris::real tNLA_relaxation_parameter = 1.0;
    int         tNLA_max_iter             = 2;

    int         tTSA_Num_Time_Steps = 1;
    moris::real tTSA_Time_Frame     = 1.0e0;

    /* ------------------------------------------------------------------------ */
    // Minimum level set value
    moris::real tMinLevs = 1.0e-8;

    /* ------------------------------------------------------------------------ */
    // Minimum level set value
    bool tUseGhost = true;

    /* ------------------------------------------------------------------------ */
    // Output Config

    std::string tOutputFileName = "Heated_Sphere.exo";

    /* ------------------------------------------------------------------------ */

    // Level set function for diamond shaped wedge
    moris::real
    Inclusion(
            const moris::Matrix< DDRMat >& aCoordinates,
            const Vector< real >&          aGeometryParameters )
    {
        // distance from sphere center
        moris::real tDx = aCoordinates( 0 ) - tCenterX;
        moris::real tDy = aCoordinates( 1 ) - tCenterY;
        moris::real tDz = aCoordinates( 2 ) - tCenterZ;

        // Compute Signed-Distance field
        moris::real tVal = tRadius - std::sqrt( tDx * tDx + tDy * tDy + tDz * tDz );

        // clean return value to return non-zero value
        return std::abs( tVal ) < tMinLevs ? tMinLevs : tVal;
    }

    /* ------------------------------------------------------------------------ */

    // Constant function for properties
    void
    Func_Const(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    void
    Func_Exact_Temperature(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        // distance from sphere center
        moris::real tDx = aFIManager->get_IP_geometry_interpolator()->valx()( 0 ) - tCenterX;
        moris::real tDy = aFIManager->get_IP_geometry_interpolator()->valx()( 1 ) - tCenterY;
        moris::real tDz = aFIManager->get_IP_geometry_interpolator()->valx()( 2 ) - tCenterZ;

        // Compute Signed-Distance field
        moris::real distance = std::sqrt( tDx * tDx + tDy * tDy + tDz * tDz );

        if ( distance <= tRadius )
        {
            aPropMatrix = { { 3.0 / 8.0 - distance * distance / 6.0 } };
        }
        else
        {
            aPropMatrix = { { ( 1.0 / 3.0 ) * ( 1.0 / distance - 1.0 ) } };
        }
    }

    /* ------------------------------------------------------------------------ */

    void
    Func_Exact_TemperatureGradient(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        // distance from sphere center
        moris::real tDx = aFIManager->get_IP_geometry_interpolator()->valx()( 0 ) - tCenterX;
        moris::real tDy = aFIManager->get_IP_geometry_interpolator()->valx()( 1 ) - tCenterY;
        moris::real tDz = aFIManager->get_IP_geometry_interpolator()->valx()( 2 ) - tCenterZ;

        // Compute Signed-Distance field
        moris::real distance = std::sqrt( tDx * tDx + tDy * tDy + tDz * tDz );

        moris::real ddistdx = tDx / distance;
        moris::real ddistdy = tDy / distance;
        moris::real ddistdz = tDz / distance;

        // set size for aPropMatrix
        aPropMatrix.set_size( 3, 1, 0.0 );

        // spatial gradients of analytic temperature distribution
        if ( distance <= tRadius )
        {
            aPropMatrix( 0, 0 ) = -2.0 / 6.0 * distance * ddistdx;
            aPropMatrix( 1, 0 ) = -2.0 / 6.0 * distance * ddistdy;
            aPropMatrix( 2, 0 ) = -2.0 / 6.0 * distance * ddistdz;
        }
        else
        {
            aPropMatrix( 0, 0 ) = -( 1.0 / 3.0 ) * ( 1.0 / distance / distance * ddistdx );
            aPropMatrix( 1, 0 ) = -( 1.0 / 3.0 ) * ( 1.0 / distance / distance * ddistdy );
            aPropMatrix( 2, 0 ) = -( 1.0 / 3.0 ) * ( 1.0 / distance / distance * ddistdz );
        }
    }

    /* ------------------------------------------------------------------------ */

    bool
    Output_Criterion( moris::tsa::Time_Solver* aTimeSolver )
    {
        return true;
    }

    /* ------------------------------------------------------------------------ */

    void
    OPTParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_opt_problem_parameter_list() );

        aParameterLists( 0 ).set( "is_optimization_problem", false );
    }

    void
    HMRParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_hmr_parameter_list() );

        aParameterLists( 0 ).set( "number_of_elements_per_dimension", tNumElemsPerDim );
        aParameterLists( 0 ).set( "domain_dimensions", tDomainDims );
        aParameterLists( 0 ).set( "domain_offset", tDomainOffset );
        aParameterLists( 0 ).set( "domain_sidesets", tDomainSidesets );
        aParameterLists( 0 ).set( "lagrange_output_meshes", "0" );

        aParameterLists( 0 ).set( "lagrange_orders", std::to_string( gInterpolationOrder ) );
        aParameterLists( 0 ).set( "lagrange_pattern", "0" );
        aParameterLists( 0 ).set( "bspline_orders", std::to_string( gInterpolationOrder ) );
        aParameterLists( 0 ).set( "bspline_pattern", "0" );

        aParameterLists( 0 ).set( "lagrange_to_bspline", "0" );

        aParameterLists( 0 ).set( "truncate_bsplines", 1 );
        aParameterLists( 0 ).set( "refinement_buffer", tRefineBuffer );
        aParameterLists( 0 ).set( "staircase_buffer", tRefineBuffer );
        aParameterLists( 0 ).set( "initial_refinement", "1" );
        aParameterLists( 0 ).set( "initial_refinement_pattern", "0" );

        aParameterLists( 0 ).set( "use_number_aura", 1 );

        aParameterLists( 0 ).set( "use_multigrid", 0 );
        aParameterLists( 0 ).set( "severity_level", 0 );
        //
        //        aParameterLists( 0 ).set( "lagrange_input_meshes", "0");
    }

    void
    XTKParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_xtk_parameter_list() );
        aParameterLists( 0 ).set( "decompose", true );
        aParameterLists( 0 ).set( "decomposition_type", "conformal" );
        aParameterLists( 0 ).set( "enrich", true );
        aParameterLists( 0 ).set( "basis_rank", "bspline" );
        aParameterLists( 0 ).set( "enrich_mesh_indices", "0" );
        aParameterLists( 0 ).set( "ghost_stab", tUseGhost );
        aParameterLists( 0 ).set( "multigrid", false );
        aParameterLists( 0 ).set( "verbose", true );
        aParameterLists( 0 ).set( "print_enriched_ig_mesh", false );
        aParameterLists( 0 ).set( "exodus_output_XTK_ig_mesh", false );
        aParameterLists( 0 ).set( "high_to_low_dbl_side_sets", true );
    }

    void
    GENParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_gen_parameter_list() );

        // Geometry parameter lists
        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists( 1 ).set( "field_function_name", "Inclusion" );
    }

    void
    FEMParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.hack_for_legacy_fem();
        // create a cell of cell of parameter list for fem

        //------------------------------------------------------------------------------

        // properties for inclusion

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropDensityInner" );
        aParameterLists( 0 ).set( "function_parameters", tDensityInner );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropCapacityInner" );
        aParameterLists( 0 ).set( "function_parameters", tCapacityInner );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropConductivityInner" );
        aParameterLists( 0 ).set( "function_parameters", tConductivityInner );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropHeatLoadInner" );
        aParameterLists( 0 ).set( "function_parameters", tHeatLoadInner );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // properties for outer material

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropDensityOuter" );
        aParameterLists( 0 ).set( "function_parameters", tDensityOuter );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropCapacityOuter" );
        aParameterLists( 0 ).set( "function_parameters", tCapacityOuter );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropConductivityOuter" );
        aParameterLists( 0 ).set( "function_parameters", tConductivityOuter );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropHeatLoadOuter" );
        aParameterLists( 0 ).set( "function_parameters", tHeatLoadOuter );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // time continuity weights
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropWeightCurrent" );
        aParameterLists( 0 ).set( "function_parameters", "100.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropWeightPrevious" );
        aParameterLists( 0 ).set( "function_parameters", "100.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // initial condition
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropInitialCondition" );
        aParameterLists( 0 ).set( "value_function", "Func_Exact_Temperature" );

        // temperature at outer surface condition
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropImposedTemperature" );
        aParameterLists( 0 ).set( "value_function", "Func_Exact_Temperature" );

        // create parameter list for property 4
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropExactTemperature" );
        aParameterLists( 0 ).set( "value_function", "Func_Exact_Temperature" );

        // create parameter list for property 5
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropExactTemperatureGradient" );
        aParameterLists( 0 ).set( "value_function", "Func_Exact_TemperatureGradient" );

        //------------------------------------------------------------------------------

        // create parameter list for constitutive model - Inclusion
        aParameterLists( 1 ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists( 1 ).set( "constitutive_name", "CMDiffusionInner" );
        aParameterLists( 1 ).set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
        aParameterLists( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists( 1 ).set( "properties",
                "PropConductivityInner , Conductivity;"
                "PropDensityInner      , Density;"
                "PropCapacityInner     , HeatCapacity" );

        // create parameter list for constitutive model - Outer Material
        aParameterLists( 1 ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists( 1 ).set( "constitutive_name", "CMDiffusionOuter" );
        aParameterLists( 1 ).set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
        aParameterLists( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists( 1 ).set( "properties",
                "PropConductivityOuter , Conductivity;"
                "PropDensityOuter      , Density;"
                "PropCapacityOuter     , HeatCapacity" );

        //------------------------------------------------------------------------------

        // create parameter list for ghost stabilization parameter for inclusion
        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPGPTempInner" );
        aParameterLists( 2 ).set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        aParameterLists( 2 ).set( "function_parameters", "0.01" );
        aParameterLists( 2 ).set( "leader_properties", "PropConductivityInner,Material" );

        // create parameter list for ghost stabilization parameter for outer material
        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPGPTempOuter" );
        aParameterLists( 2 ).set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        aParameterLists( 2 ).set( "function_parameters", "0.01" );
        aParameterLists( 2 ).set( "leader_properties", "PropConductivityOuter,Material" );

        // create parameter list for Nitsche stabilization parameter for inclusion-outer material interface
        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPInterfaceNitsche" );
        aParameterLists( 2 ).set( "stabilization_type", fem::Stabilization_Type::NITSCHE_INTERFACE );
        aParameterLists( 2 ).set( "function_parameters", "100.0" );
        aParameterLists( 2 ).set( "leader_properties", "PropConductivityInner,Material" );
        aParameterLists( 2 ).set( "follower_properties", "PropConductivityOuter,Material" );

        // create parameter list for DBC on outer surface
        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPNitscheTemp" );
        aParameterLists( 2 ).set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        aParameterLists( 2 ).set( "function_parameters", "100.0" );
        aParameterLists( 2 ).set( "leader_properties", "PropConductivityOuter,Material" );

        //------------------------------------------------------------------------------
        // create IWG for inclusion - bulk diffusion
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGDiffusionInnerBulk" );
        aParameterLists( 3 ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        aParameterLists( 3 ).set( "dof_residual", "TEMP" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMDiffusionInner,Diffusion" );
        aParameterLists( 3 ).set( "leader_properties", "PropHeatLoadInner,Load" );
        aParameterLists( 3 ).set( "mesh_set_names", tInnerPhase );

        // create IWG for outer material - bulk diffusion
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGDiffusionOuterBulk" );
        aParameterLists( 3 ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        aParameterLists( 3 ).set( "dof_residual", "TEMP" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMDiffusionOuter,Diffusion" );
        aParameterLists( 3 ).set( "leader_properties", "PropHeatLoadOuter,Load" );
        aParameterLists( 3 ).set( "mesh_set_names", tOuterPhase );

        // create parameter list for outer boundary conditions
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGOuterSurfaceTemp" );
        aParameterLists( 3 ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists( 3 ).set( "dof_residual", "TEMP" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 3 ).set( "leader_properties", "PropImposedTemperature,Dirichlet" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMDiffusionOuter,Diffusion" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPNitscheTemp,DirichletNitsche" );
        aParameterLists( 3 ).set( "mesh_set_names", tOuterSurface );

        // create parameter list for interface conditions
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGInterfaceInnerOuterTEMP" );
        aParameterLists( 3 ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE );
        aParameterLists( 3 ).set( "dof_residual", "TEMP" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 3 ).set( "follower_dof_dependencies", "TEMP" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMDiffusionInner,Diffusion" );
        aParameterLists( 3 ).set( "follower_constitutive_models", "CMDiffusionOuter,Diffusion" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPInterfaceNitsche,NitscheInterface" );
        aParameterLists( 3 ).set( "mesh_set_names", tInterface );

        if ( tUseGhost )
        {
            // create IWG for outer material - ghost
            aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists( 3 ).set( "IWG_name", "IWGGPInnerTemp" );
            aParameterLists( 3 ).set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists( 3 ).set( "dof_residual", "TEMP" );
            aParameterLists( 3 ).set( "leader_dof_dependencies", "TEMP" );
            aParameterLists( 3 ).set( "follower_dof_dependencies", "TEMP" );
            aParameterLists( 3 ).set( "stabilization_parameters", "SPGPTempInner,GhostSP" );
            aParameterLists( 3 ).set( "mesh_set_names", tInnerPhaseGhost );

            // create IWG for outer material - ghost
            aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists( 3 ).set( "IWG_name", "IWGGPInnerTemp" );
            aParameterLists( 3 ).set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists( 3 ).set( "dof_residual", "TEMP" );
            aParameterLists( 3 ).set( "leader_dof_dependencies", "TEMP" );
            aParameterLists( 3 ).set( "follower_dof_dependencies", "TEMP" );
            aParameterLists( 3 ).set( "stabilization_parameters", "SPGPTempOuter,GhostSP" );
            aParameterLists( 3 ).set( "mesh_set_names", tOuterPhaseGhost );
        }

        //        // create IWG for time continuity
        //        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        //        aParameterLists( 3 ).set( "IWG_name",                   "IWGTimeContinuityTemp") ;
        //        aParameterLists( 3 ).set( "IWG_type",                    fem::IWG_Type::TIME_CONTINUITY_DOF ) ;
        //        aParameterLists( 3 ).set( "dof_residual",               "TEMP") ;
        //        aParameterLists( 3 ).set( "leader_dof_dependencies",    "TEMP") ;
        //        aParameterLists( 3 ).set( "leader_properties",
        //                "PropWeightCurrent   ,WeightCurrent;"
        //                "PropWeightPrevious  ,WeightPrevious;"
        //                "PropInitialCondition,InitialCondition") ;
        //        aParameterLists( 3 ).set( "mesh_set_names",             tTotalDomain );
        //        aParameterLists( 3 ).set( "time_continuity",            true );
        //
        //------------------------------------------------------------------------------
        // Nodal Temperature IQI
        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkTEMP" );
        aParameterLists( 4 ).set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists( 4 ).set( "dof_quantity", "TEMP" );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 4 ).set( "vectorial_field_index", 0 );
        aParameterLists( 4 ).set( "mesh_set_names", tTotalDomain );

        // Nodal Analytic Temperature IQI
        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkTEMPAnalytic" );
        aParameterLists( 4 ).set( "IQI_type", fem::IQI_Type::PROPERTY );
        aParameterLists( 4 ).set( "leader_properties", "PropExactTemperature,Property" );
        aParameterLists( 4 ).set( "mesh_set_names", tTotalDomain );

        // L2 Error IQI
        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkL2Error" );
        aParameterLists( 4 ).set( "IQI_type", fem::IQI_Type::L2_ERROR_ANALYTIC );
        aParameterLists( 4 ).set( "dof_quantity", "TEMP" );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 4 ).set( "leader_properties", "PropExactTemperature,L2Check" );
        aParameterLists( 4 ).set( "mesh_set_names", tTotalDomain );

        // H1 Semi-Norm Error IQI
        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkH1Error" );
        aParameterLists( 4 ).set( "IQI_type", fem::IQI_Type::H1_ERROR_ANALYTIC );
        aParameterLists( 4 ).set( "dof_quantity", "TEMP" );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 4 ).set( "leader_properties", "PropExactTemperatureGradient,H1Check" );
        aParameterLists( 4 ).set( "mesh_set_names", tTotalDomain );

        // Total volume of the sphere and the block around it
        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIVolume" );
        aParameterLists( 4 ).set( "IQI_type", fem::IQI_Type::VOLUME );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 4 ).set( "mesh_set_names", tTotalDomain );

        // surface area of the sphere
        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQI_Sphere_Surface" );
        aParameterLists( 4 ).set( "IQI_type", fem::IQI_Type::VOLUME );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 4 ).set( "mesh_set_names", tSphereSurface );

        // sphere surface heat flux
        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQI_Sphere_Outer_Surface_Heat_Flux" );
        aParameterLists( 4 ).set( "IQI_type", fem::IQI_Type::TRACTION );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 4 ).set( "leader_constitutive_models", "CMDiffusionInner,TractionCM" );
        aParameterLists( 4 ).set( "mesh_set_names", tSphereSurface );
        aParameterLists( 4 ).set( "vectorial_field_index", 0 );

        // block inner surface heat flux
        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQI_Block_Inner_Surface_Heat_Flux" );
        aParameterLists( 4 ).set( "IQI_type", fem::IQI_Type::TRACTION );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 4 ).set( "leader_constitutive_models", "CMDiffusionOuter,TractionCM" );
        aParameterLists( 4 ).set( "mesh_set_names", tBlockInnerSurface );
        aParameterLists( 4 ).set( "vectorial_field_index", 0 );

        // error in heat fluxes between the two materials
        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQI_Surface_Heat_Flux_Jump" );
        aParameterLists( 4 ).set( "IQI_type", fem::IQI_Type::JUMP_TRACTION );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 4 ).set( "follower_dof_dependencies", "TEMP" );
        aParameterLists( 4 ).set( "leader_constitutive_models", "CMDiffusionInner,TractionCM" );
        aParameterLists( 4 ).set( "follower_constitutive_models", "CMDiffusionOuter,TractionCM" );
        aParameterLists( 4 ).set( "mesh_set_names", tInterface );
        aParameterLists( 4 ).set( "vectorial_field_index", 0 );

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        aParameterLists( 5 ).add_parameter_list( prm::create_computation_parameter_list() );
    }

    void
    SOLParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL ) );

        aParameterLists( 1 ).add_parameter_list( moris::prm::create_linear_solver_parameter_list() );

        aParameterLists( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
        aParameterLists( 2 ).set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        aParameterLists( 2 ).set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
        aParameterLists( 2 ).set( "NLA_max_iter", tNLA_max_iter );

        aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
        aParameterLists( 3 ).set( "NLA_DofTypes", "TEMP" );

        aParameterLists( 4 ).add_parameter_list( moris::prm::create_time_solver_algorithm_parameter_list() );
        aParameterLists( 4 ).set( "TSA_Num_Time_Steps", tTSA_Num_Time_Steps );
        aParameterLists( 4 ).set( "TSA_Time_Frame", tTSA_Time_Frame );

        aParameterLists( 5 ).add_parameter_list( moris::prm::create_time_solver_parameter_list() );
        aParameterLists( 5 ).set( "TSA_DofTypes", "TEMP" );
        aParameterLists( 5 ).set( "TSA_Initialize_Sol_Vec", "TEMP,0.0" );
        aParameterLists( 5 ).set( "TSA_Output_Indices", "0" );
        aParameterLists( 5 ).set( "TSA_Output_Criteria", "Output_Criterion" );

        aParameterLists( 6 ).add_parameter_list( moris::prm::create_solver_warehouse_parameterlist() );

        aParameterLists( 7 ).add_parameter_list( moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::NONE ) );
    }

    void
    MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_msi_parameter_list() );
        aParameterLists( 0 ).set( "order_adofs_by_host", false );
    }

    void
    VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_vis_parameter_list() );
        aParameterLists( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName ) );
        aParameterLists( 0 ).set( "Mesh_Type", vis::VIS_Mesh_Type::STANDARD );
        aParameterLists( 0 ).set( "Set_Names", tTotalDomain + "," + tAllInterfaces + "," + tInterfaceSingleSides );
        aParameterLists( 0 ).set( "Field_Names",
                "TEMP,TEMP_ANALYTIC,"
                "L2_ERROR_ANALYTIC,H1_ERROR_ANALYTIC,"
                "VOLUME,"
                "SPHERE_SURFACE_AREA,ELEMENTAL_SPHERE_SURFACE_AREA,"
                "SPHERE_SURFACE_HEAT_FLUX,TOTAL_SPHERE_SURFACE_HEAT_FLUX,"
                "BLOCK_INNER_SURFACE_HEAT_FLUX,TOTAL_BLOCK_INNER_SURFACE_HEAT_FLUX,"
                "HEAT_FLUX_JUMP" );
        aParameterLists( 0 ).set( "Field_Type",
                "NODAL,NODAL,"
                "GLOBAL,GLOBAL,"
                "GLOBAL,"
                "GLOBAL,FACETED_INT,"
                "FACETED_AVG,GLOBAL,"
                "FACETED_AVG,GLOBAL,"
                "FACETED_AVG" );
        aParameterLists( 0 ).set( "IQI_Names",
                "IQIBulkTEMP,IQIBulkTEMPAnalytic,"
                "IQIBulkL2Error,IQIBulkH1Error,"
                "IQIVolume,"
                "IQI_Sphere_Surface, IQI_Sphere_Surface,"
                "IQI_Sphere_Outer_Surface_Heat_Flux,IQI_Sphere_Outer_Surface_Heat_Flux,"
                "IQI_Block_Inner_Surface_Heat_Flux,IQI_Block_Inner_Surface_Heat_Flux,"
                "IQI_Surface_Heat_Flux_Jump" );
        aParameterLists( 0 ).set( "Save_Frequency", 1 );
    }

    void
    MORISGENERALParameterList( Module_Parameter_Lists& aParameterLists )
    {
    }

    /* ------------------------------------------------------------------------ */
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
