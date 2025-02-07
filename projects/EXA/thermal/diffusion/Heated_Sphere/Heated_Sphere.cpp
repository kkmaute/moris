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

    Vector< uint > tNumElemsPerDim = { 8, 8, 8 };
    Vector< real > tDomainDims     = { 2.0, 2.0, 2.0 };
    std::string tDomainOffset   = "-1.0, -1.0, -1.0";

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
        aParameterLists.set( "is_optimization_problem", false );
    }

    void
    HMRParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "number_of_elements_per_dimension", tNumElemsPerDim );
        aParameterLists.set( "domain_dimensions", tDomainDims );
        aParameterLists.set( "domain_offset", tDomainOffset );
        aParameterLists.set( "lagrange_output_meshes", "0" );

        aParameterLists.set( "lagrange_orders", std::to_string( gInterpolationOrder ) );
        aParameterLists.set( "lagrange_pattern", "0" );
        aParameterLists.set( "bspline_orders", std::to_string( gInterpolationOrder ) );
        aParameterLists.set( "bspline_pattern", "0" );

        aParameterLists.set( "lagrange_to_bspline", "0" );

        aParameterLists.set( "truncate_bsplines", 1 );
        aParameterLists.set( "refinement_buffer", tRefineBuffer );
        aParameterLists.set( "staircase_buffer", tRefineBuffer );
        aParameterLists.set( "initial_refinement", "1" );
        aParameterLists.set( "initial_refinement_pattern", "0" );

        aParameterLists.set( "use_number_aura", 1 );

        aParameterLists.set( "use_multigrid", 0 );
        aParameterLists.set( "severity_level", 0 );
        //
        //        aParameterLists.set( "lagrange_input_meshes", "0");
    }

    void
    XTKParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "decompose", true );
        aParameterLists.set( "decomposition_type", "conformal" );
        aParameterLists.set( "enrich", true );
        aParameterLists.set( "basis_rank", "bspline" );
        aParameterLists.set( "enrich_mesh_indices", "0" );
        aParameterLists.set( "ghost_stab", tUseGhost );
        aParameterLists.set( "multigrid", false );
        aParameterLists.set( "verbose", true );
        aParameterLists.set( "print_enriched_ig_mesh", false );
        aParameterLists.set( "exodus_output_XTK_ig_mesh", false );
        aParameterLists.set( "high_to_low_dbl_side_sets", true );
    }

    void
    GENParameterList( Module_Parameter_Lists& aParameterLists )
    {
        // Geometry parameter lists
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists.set( "field_function_name", "Inclusion" );
    }

    void
    FEMParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.hack_for_legacy_fem();
        // create a cell of cell of parameter list for fem

        //------------------------------------------------------------------------------

        // properties for inclusion

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDensityInner" );
        aParameterLists.set( "function_parameters", tDensityInner );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropCapacityInner" );
        aParameterLists.set( "function_parameters", tCapacityInner );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropConductivityInner" );
        aParameterLists.set( "function_parameters", tConductivityInner );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropHeatLoadInner" );
        aParameterLists.set( "function_parameters", tHeatLoadInner );
        aParameterLists.set( "value_function", "Func_Const" );

        // properties for outer material

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDensityOuter" );
        aParameterLists.set( "function_parameters", tDensityOuter );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropCapacityOuter" );
        aParameterLists.set( "function_parameters", tCapacityOuter );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropConductivityOuter" );
        aParameterLists.set( "function_parameters", tConductivityOuter );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropHeatLoadOuter" );
        aParameterLists.set( "function_parameters", tHeatLoadOuter );
        aParameterLists.set( "value_function", "Func_Const" );

        // time continuity weights
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropWeightCurrent" );
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropWeightPrevious" );
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // initial condition
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropInitialCondition" );
        aParameterLists.set( "value_function", "Func_Exact_Temperature" );

        // temperature at outer surface condition
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropImposedTemperature" );
        aParameterLists.set( "value_function", "Func_Exact_Temperature" );

        // create parameter list for property 4
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropExactTemperature" );
        aParameterLists.set( "value_function", "Func_Exact_Temperature" );

        // create parameter list for property 5
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropExactTemperatureGradient" );
        aParameterLists.set( "value_function", "Func_Exact_TemperatureGradient" );

        //------------------------------------------------------------------------------

        // create parameter list for constitutive model - Inclusion
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMDiffusionInner" );
        aParameterLists.set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists.set( "properties",
                "PropConductivityInner , Conductivity;"
                "PropDensityInner      , Density;"
                "PropCapacityInner     , HeatCapacity" );

        // create parameter list for constitutive model - Outer Material
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMDiffusionOuter" );
        aParameterLists.set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists.set( "properties",
                "PropConductivityOuter , Conductivity;"
                "PropDensityOuter      , Density;"
                "PropCapacityOuter     , HeatCapacity" );

        //------------------------------------------------------------------------------

        // create parameter list for ghost stabilization parameter for inclusion
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPGPTempInner" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        aParameterLists.set( "function_parameters", "0.01" );
        aParameterLists.set( "leader_properties", "PropConductivityInner,Material" );

        // create parameter list for ghost stabilization parameter for outer material
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPGPTempOuter" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        aParameterLists.set( "function_parameters", "0.01" );
        aParameterLists.set( "leader_properties", "PropConductivityOuter,Material" );

        // create parameter list for Nitsche stabilization parameter for inclusion-outer material interface
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPInterfaceNitsche" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::NITSCHE_INTERFACE );
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "leader_properties", "PropConductivityInner,Material" );
        aParameterLists.set( "follower_properties", "PropConductivityOuter,Material" );

        // create parameter list for DBC on outer surface
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPNitscheTemp" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "leader_properties", "PropConductivityOuter,Material" );

        //------------------------------------------------------------------------------
        // create IWG for inclusion - bulk diffusion
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGDiffusionInnerBulk" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusionInner,Diffusion" );
        aParameterLists.set( "leader_properties", "PropHeatLoadInner,Load" );
        aParameterLists.set( "mesh_set_names", tInnerPhase );

        // create IWG for outer material - bulk diffusion
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGDiffusionOuterBulk" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusionOuter,Diffusion" );
        aParameterLists.set( "leader_properties", "PropHeatLoadOuter,Load" );
        aParameterLists.set( "mesh_set_names", tOuterPhase );

        // create parameter list for outer boundary conditions
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGOuterSurfaceTemp" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "leader_properties", "PropImposedTemperature,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusionOuter,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheTemp,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", tOuterSurface );

        // create parameter list for interface conditions
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGInterfaceInnerOuterTEMP" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "follower_dof_dependencies", "TEMP" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusionInner,Diffusion" );
        aParameterLists.set( "follower_constitutive_models", "CMDiffusionOuter,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPInterfaceNitsche,NitscheInterface" );
        aParameterLists.set( "mesh_set_names", tInterface );

        if ( tUseGhost )
        {
            // create IWG for outer material - ghost
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "IWGGPInnerTemp" );
            aParameterLists.set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists.set( "dof_residual", "TEMP" );
            aParameterLists.set( "leader_dof_dependencies", "TEMP" );
            aParameterLists.set( "follower_dof_dependencies", "TEMP" );
            aParameterLists.set( "stabilization_parameters", "SPGPTempInner,GhostSP" );
            aParameterLists.set( "mesh_set_names", tInnerPhaseGhost );

            // create IWG for outer material - ghost
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "IWGGPInnerTemp" );
            aParameterLists.set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists.set( "dof_residual", "TEMP" );
            aParameterLists.set( "leader_dof_dependencies", "TEMP" );
            aParameterLists.set( "follower_dof_dependencies", "TEMP" );
            aParameterLists.set( "stabilization_parameters", "SPGPTempOuter,GhostSP" );
            aParameterLists.set( "mesh_set_names", tOuterPhaseGhost );
        }

        //        // create IWG for time continuity
        //        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        //        aParameterLists.set( "IWG_name",                   "IWGTimeContinuityTemp") ;
        //        aParameterLists.set( "IWG_type",                    fem::IWG_Type::TIME_CONTINUITY_DOF ) ;
        //        aParameterLists.set( "dof_residual",               "TEMP") ;
        //        aParameterLists.set( "leader_dof_dependencies",    "TEMP") ;
        //        aParameterLists.set( "leader_properties",
        //                "PropWeightCurrent   ,WeightCurrent;"
        //                "PropWeightPrevious  ,WeightPrevious;"
        //                "PropInitialCondition,InitialCondition") ;
        //        aParameterLists.set( "mesh_set_names",             tTotalDomain );
        //        aParameterLists.set( "time_continuity",            true );
        //
        //------------------------------------------------------------------------------
        // Nodal Temperature IQI
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkTEMP" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "mesh_set_names", tTotalDomain );

        // Nodal Analytic Temperature IQI
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkTEMPAnalytic" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::PROPERTY );
        aParameterLists.set( "leader_properties", "PropExactTemperature,Property" );
        aParameterLists.set( "mesh_set_names", tTotalDomain );

        // L2 Error IQI
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkL2Error" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::L2_ERROR_ANALYTIC );
        aParameterLists.set( "dof_quantity", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "leader_properties", "PropExactTemperature,L2Check" );
        aParameterLists.set( "mesh_set_names", tTotalDomain );

        // H1 Semi-Norm Error IQI
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkH1Error" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::H1_ERROR_ANALYTIC );
        aParameterLists.set( "dof_quantity", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "leader_properties", "PropExactTemperatureGradient,H1Check" );
        aParameterLists.set( "mesh_set_names", tTotalDomain );

        // Total volume of the sphere and the block around it
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIVolume" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::VOLUME );
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "mesh_set_names", tTotalDomain );

        // surface area of the sphere
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQI_Sphere_Surface" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::VOLUME );
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "mesh_set_names", tSphereSurface );

        // sphere surface heat flux
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQI_Sphere_Outer_Surface_Heat_Flux" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::TRACTION );
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusionInner,TractionCM" );
        aParameterLists.set( "mesh_set_names", tSphereSurface );
        aParameterLists.set( "vectorial_field_index", 0 );

        // block inner surface heat flux
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQI_Block_Inner_Surface_Heat_Flux" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::TRACTION );
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusionOuter,TractionCM" );
        aParameterLists.set( "mesh_set_names", tBlockInnerSurface );
        aParameterLists.set( "vectorial_field_index", 0 );

        // error in heat fluxes between the two materials
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQI_Surface_Heat_Flux_Jump" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::JUMP_TRACTION );
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "follower_dof_dependencies", "TEMP" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusionInner,TractionCM" );
        aParameterLists.set( "follower_constitutive_models", "CMDiffusionOuter,TractionCM" );
        aParameterLists.set( "mesh_set_names", tInterface );
        aParameterLists.set( "vectorial_field_index", 0 );

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
        aParameterLists.set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        aParameterLists.set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
        aParameterLists.set( "NLA_max_iter", tNLA_max_iter );

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "NLA_DofTypes", "TEMP" );

        aParameterLists( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "TSA_Num_Time_Steps", tTSA_Num_Time_Steps );
        aParameterLists.set( "TSA_Time_Frame", tTSA_Time_Frame );

        aParameterLists( SOL::TIME_SOLVERS ).add_parameter_list();
        aParameterLists.set( "TSA_DofTypes", "TEMP" );
        aParameterLists.set( "TSA_Initialize_Sol_Vec", "TEMP,0.0" );
        aParameterLists.set( "TSA_Output_Indices", "0" );
        aParameterLists.set( "TSA_Output_Criteria", "Output_Criterion" );

        aParameterLists( SOL::PRECONDITIONERS ).add_parameter_list(  sol::PreconditionerType::NONE );
    }

    void
    MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "order_adofs_by_host", false );
    }

    void
    VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName ) );
        aParameterLists.set( "Mesh_Type", vis::VIS_Mesh_Type::STANDARD );
        aParameterLists.set( "Set_Names", tTotalDomain + "," + tAllInterfaces + "," + tInterfaceSingleSides );
        aParameterLists.set( "Field_Names",
                "TEMP,TEMP_ANALYTIC,"
                "L2_ERROR_ANALYTIC,H1_ERROR_ANALYTIC,"
                "VOLUME,"
                "SPHERE_SURFACE_AREA,ELEMENTAL_SPHERE_SURFACE_AREA,"
                "SPHERE_SURFACE_HEAT_FLUX,TOTAL_SPHERE_SURFACE_HEAT_FLUX,"
                "BLOCK_INNER_SURFACE_HEAT_FLUX,TOTAL_BLOCK_INNER_SURFACE_HEAT_FLUX,"
                "HEAT_FLUX_JUMP" );
        aParameterLists.set( "Field_Type",
                "NODAL,NODAL,"
                "GLOBAL,GLOBAL,"
                "GLOBAL,"
                "GLOBAL,FACETED_INT,"
                "FACETED_AVG,GLOBAL,"
                "FACETED_AVG,GLOBAL,"
                "FACETED_AVG" );
        aParameterLists.set( "IQI_Names",
                "IQIBulkTEMP,IQIBulkTEMPAnalytic,"
                "IQIBulkL2Error,IQIBulkH1Error,"
                "IQIVolume,"
                "IQI_Sphere_Surface, IQI_Sphere_Surface,"
                "IQI_Sphere_Outer_Surface_Heat_Flux,IQI_Sphere_Outer_Surface_Heat_Flux,"
                "IQI_Block_Inner_Surface_Heat_Flux,IQI_Block_Inner_Surface_Heat_Flux,"
                "IQI_Surface_Heat_Flux_Jump" );
        aParameterLists.set( "Save_Frequency", 1 );
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
