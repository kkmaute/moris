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
#include "fn_PRM_FEM_Parameters.hpp"
#include "fn_PRM_MSI_Parameters.hpp"
#include "fn_PRM_SOL_Parameters.hpp"
#include "fn_PRM_VIS_Parameters.hpp"
#include "fn_PRM_HMR_Parameters.hpp"
#include "fn_PRM_GEN_Parameters.hpp"
#include "fn_PRM_XTK_Parameters.hpp"
#include "fn_PRM_OPT_Parameters.hpp"
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
            const moris::Matrix< DDRMat >&     aCoordinates,
            const moris::Cell< moris::real* >& aGeometryParameters )
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
    Func_Const( moris::Matrix<
                        moris::DDRMat >&                   aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    void
    Func_Exact_Temperature(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
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
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
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
    OPTParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_opt_problem_parameter_list();

        tParameterlist( 0 )( 0 ).set( "is_optimization_problem", false );
    }

    void
    HMRParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_hmr_parameter_list();

        tParameterlist( 0 )( 0 ).set( "number_of_elements_per_dimension", tNumElemsPerDim );
        tParameterlist( 0 )( 0 ).set( "domain_dimensions", tDomainDims );
        tParameterlist( 0 )( 0 ).set( "domain_offset", tDomainOffset );
        tParameterlist( 0 )( 0 ).set( "domain_sidesets", tDomainSidesets );
        tParameterlist( 0 )( 0 ).set( "lagrange_output_meshes", "0" );

        tParameterlist( 0 )( 0 ).set( "lagrange_orders", std::to_string( gInterpolationOrder ) );
        tParameterlist( 0 )( 0 ).set( "lagrange_pattern", "0" );
        tParameterlist( 0 )( 0 ).set( "bspline_orders", std::to_string( gInterpolationOrder ) );
        tParameterlist( 0 )( 0 ).set( "bspline_pattern", "0" );

        tParameterlist( 0 )( 0 ).set( "lagrange_to_bspline", "0" );

        tParameterlist( 0 )( 0 ).set( "truncate_bsplines", 1 );
        tParameterlist( 0 )( 0 ).set( "refinement_buffer", tRefineBuffer );
        tParameterlist( 0 )( 0 ).set( "staircase_buffer", tRefineBuffer );
        tParameterlist( 0 )( 0 ).set( "initial_refinement", "1" );
        tParameterlist( 0 )( 0 ).set( "initial_refinement_pattern", "0" );

        tParameterlist( 0 )( 0 ).set( "use_number_aura", 1 );

        tParameterlist( 0 )( 0 ).set( "use_multigrid", 0 );
        tParameterlist( 0 )( 0 ).set( "severity_level", 0 );
        //
        //        tParameterlist( 0 )( 0 ).set( "lagrange_input_meshes", "0");
    }

    void
    XTKParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_xtk_parameter_list();
        tParameterlist( 0 )( 0 ).set( "decompose", true );
        tParameterlist( 0 )( 0 ).set( "decomposition_type", "conformal" );
        tParameterlist( 0 )( 0 ).set( "enrich", true );
        tParameterlist( 0 )( 0 ).set( "basis_rank", "bspline" );
        tParameterlist( 0 )( 0 ).set( "enrich_mesh_indices", "0" );
        tParameterlist( 0 )( 0 ).set( "ghost_stab", tUseGhost );
        tParameterlist( 0 )( 0 ).set( "multigrid", false );
        tParameterlist( 0 )( 0 ).set( "verbose", true );
        tParameterlist( 0 )( 0 ).set( "print_enriched_ig_mesh", false );
        tParameterlist( 0 )( 0 ).set( "exodus_output_XTK_ig_mesh", false );
        tParameterlist( 0 )( 0 ).set( "high_to_low_dbl_side_sets", true );
    }

    void
    GENParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 3 );
        tParameterlist( 0 ).resize( 1 );

        // Main GEN parameter list
        tParameterlist( 0 )( 0 ) = prm::create_gen_parameter_list();

        // init geometry counter
        uint tGeoCounter = 0;

        // Geometry parameter lists
        tParameterlist( 1 ).push_back( prm::create_user_defined_geometry_parameter_list() );
        tParameterlist( 1 )( tGeoCounter ).set( "field_function_name", "Inclusion" );
    }

    void
    FEMParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterList )
    {
        // create a cell of cell of parameter list for fem
        tParameterList.resize( 8 );

        //------------------------------------------------------------------------------
        // init property counter
        uint tPropCounter = 0;

        // properties for inclusion

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ) = prm::create_property_parameter_list();
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropDensityInner" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tDensityInner );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropCapacityInner" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tCapacityInner );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropConductivityInner" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tConductivityInner );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropHeatLoadInner" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tHeatLoadInner );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // properties for outer material

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ) = prm::create_property_parameter_list();
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropDensityOuter" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tDensityOuter );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropCapacityOuter" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tCapacityOuter );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropConductivityOuter" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tConductivityOuter );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropHeatLoadOuter" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tHeatLoadOuter );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // time continuity weights
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropWeightCurrent" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "100.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropWeightPrevious" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "100.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // initial condition
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropInitialCondition" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Exact_Temperature" );
        tPropCounter++;

        // temperature at outer surface condition
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropImposedTemperature" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Exact_Temperature" );
        tPropCounter++;

        // create parameter list for property 4
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropExactTemperature" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Exact_Temperature" );
        tPropCounter++;

        // create parameter list for property 5
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropExactTemperatureGradient" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Exact_TemperatureGradient" );
        tPropCounter++;

        //------------------------------------------------------------------------------
        // init CM counter
        uint tCMCounter = 0;

        // create parameter list for constitutive model - Inclusion
        tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_name", "CMDiffusionInner" );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::DIFF_LIN_ISO ) );
        tParameterList( 1 )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( 1 )( tCMCounter ).set( "properties",
                "PropConductivityInner , Conductivity;"
                "PropDensityInner      , Density;"
                "PropCapacityInner     , HeatCapacity" );
        tCMCounter++;

        // create parameter list for constitutive model - Outer Material
        tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_name", "CMDiffusionOuter" );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::DIFF_LIN_ISO ) );
        tParameterList( 1 )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( 1 )( tCMCounter ).set( "properties",
                "PropConductivityOuter , Conductivity;"
                "PropDensityOuter      , Density;"
                "PropCapacityOuter     , HeatCapacity" );
        tCMCounter++;

        //------------------------------------------------------------------------------
        // init SP counter
        uint tSPCounter = 0;

        // create parameter list for ghost stabilization parameter for inclusion
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", "SPGPTempInner" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::GHOST_DISPL ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", "0.01" );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", "PropConductivityInner,Material" );
        tSPCounter++;

        // create parameter list for ghost stabilization parameter for outer material
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", "SPGPTempOuter" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::GHOST_DISPL ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", "0.01" );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", "PropConductivityOuter,Material" );
        tSPCounter++;

        // create parameter list for Nitsche stabilization parameter for inclusion-outer material interface
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", "SPInterfaceNitsche" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::NITSCHE_INTERFACE ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", "100.0" );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", "PropConductivityInner,Material" );
        tParameterList( 2 )( tSPCounter ).set( "follower_properties", "PropConductivityOuter,Material" );
        tSPCounter++;

        // create parameter list for DBC on outer surface
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", "SPNitscheTemp" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::DIRICHLET_NITSCHE ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", "100.0" );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", "PropConductivityOuter,Material" );
        tSPCounter++;

        //------------------------------------------------------------------------------
        // init IWG counter
        uint tIWGCounter = 0;

        // create IWG for inclusion - bulk diffusion
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGDiffusionInnerBulk" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::SPATIALDIFF_BULK ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMDiffusionInner,Diffusion" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties", "PropHeatLoadInner,Load" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tInnerPhase );
        tIWGCounter++;

        // create IWG for outer material - bulk diffusion
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGDiffusionOuterBulk" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::SPATIALDIFF_BULK ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMDiffusionOuter,Diffusion" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties", "PropHeatLoadOuter,Load" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tOuterPhase );
        tIWGCounter++;

        // create parameter list for outer boundary conditions
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGOuterSurfaceTemp" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties", "PropImposedTemperature,Dirichlet" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMDiffusionOuter,Diffusion" );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPNitscheTemp,DirichletNitsche" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tOuterSurface );
        tIWGCounter++;

        // create parameter list for interface conditions
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGInterfaceInnerOuterTEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "follower_dof_dependencies", "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMDiffusionInner,Diffusion" );
        tParameterList( 3 )( tIWGCounter ).set( "follower_constitutive_models", "CMDiffusionOuter,Diffusion" );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPInterfaceNitsche,NitscheInterface" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tInterface );
        tIWGCounter++;

        if ( tUseGhost )
        {
            // create IWG for outer material - ghost
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGGPInnerTemp" );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::GHOST_NORMAL_FIELD ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "TEMP" );
            tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "TEMP" );
            tParameterList( 3 )( tIWGCounter ).set( "follower_dof_dependencies", "TEMP" );
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPGPTempInner,GhostSP" );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tInnerPhaseGhost );
            tIWGCounter++;

            // create IWG for outer material - ghost
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGGPInnerTemp" );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::GHOST_NORMAL_FIELD ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "TEMP" );
            tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "TEMP" );
            tParameterList( 3 )( tIWGCounter ).set( "follower_dof_dependencies", "TEMP" );
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPGPTempOuter,GhostSP" );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tOuterPhaseGhost );
            tIWGCounter++;
        }

        //        // create IWG for time continuity
        //        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        //        tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWGTimeContinuityTemp") ;
        //        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::TIME_CONTINUITY_DOF ) );
        //        tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               "TEMP") ;
        //        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies",    "TEMP") ;
        //        tParameterList( 3 )( tIWGCounter ).set( "leader_properties",
        //                "PropWeightCurrent   ,WeightCurrent;"
        //                "PropWeightPrevious  ,WeightPrevious;"
        //                "PropInitialCondition,InitialCondition") ;
        //        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             tTotalDomain );
        //        tParameterList( 3 )( tIWGCounter ).set( "time_continuity",            true );
        //        tIWGCounter++;

        //------------------------------------------------------------------------------
        // init IQI counter
        uint tIQICounter = 0;

        // Nodal Temperature IQI
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkTEMP" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::DOF ) );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "TEMP" );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 0 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tTotalDomain );
        tIQICounter++;

        // Nodal Analytic Temperature IQI
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkTEMPAnalytic" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::PROPERTY ) );
        tParameterList( 4 )( tIQICounter ).set( "leader_properties", "PropExactTemperature,Property" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tTotalDomain );
        tIQICounter++;

        // L2 Error IQI
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkL2Error" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::L2_ERROR_ANALYTIC ) );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "TEMP" );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( 4 )( tIQICounter ).set( "leader_properties", "PropExactTemperature,L2Check" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tTotalDomain );
        tIQICounter++;

        // H1 Semi-Norm Error IQI
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkH1Error" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::H1_ERROR_ANALYTIC ) );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "TEMP" );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( 4 )( tIQICounter ).set( "leader_properties", "PropExactTemperatureGradient,H1Check" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tTotalDomain );
        tIQICounter++;

        // Total volume of the sphere and the block around it
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIVolume" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::VOLUME ) );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tTotalDomain );
        tIQICounter++;

        // surface area of the sphere
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQI_Sphere_Surface" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::VOLUME ) );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tSphereSurface );
        tIQICounter++;

        // sphere surface heat flux
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQI_Sphere_Outer_Surface_Heat_Flux" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", (uint)fem::IQI_Type::TRACTION );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( 4 )( tIQICounter ).set( "leader_constitutive_models", "CMDiffusionInner,TractionCM" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tSphereSurface );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 0 );
        tIQICounter++;

        // block inner surface heat flux
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQI_Block_Inner_Surface_Heat_Flux" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", (uint)fem::IQI_Type::TRACTION );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( 4 )( tIQICounter ).set( "leader_constitutive_models", "CMDiffusionOuter,TractionCM" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tBlockInnerSurface );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 0 );
        tIQICounter++;

        // error in heat fluxes between the two materials
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQI_Surface_Heat_Flux_Jump" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", (uint)fem::IQI_Type::JUMP_TRACTION );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( 4 )( tIQICounter ).set( "follower_dof_dependencies", "TEMP" );
        tParameterList( 4 )( tIQICounter ).set( "leader_constitutive_models", "CMDiffusionInner,TractionCM" );
        tParameterList( 4 )( tIQICounter ).set( "follower_constitutive_models", "CMDiffusionOuter,TractionCM" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tInterface );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 0 );
        tIQICounter++;


        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        tParameterList( 5 ).resize( 1 );
        tParameterList( 5 )( 0 ) = prm::create_computation_parameter_list();
    }

    void
    SOLParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 8 );
        for ( uint Ik = 0; Ik < 8; Ik++ )
        {
            tParameterlist( Ik ).resize( 1 );
        }

        tParameterlist( 0 )( 0 ) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL );

        tParameterlist( 1 )( 0 ) = moris::prm::create_linear_solver_parameter_list();

        tParameterlist( 2 )( 0 ) = moris::prm::create_nonlinear_algorithm_parameter_list();
        tParameterlist( 2 )( 0 ).set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        tParameterlist( 2 )( 0 ).set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
        tParameterlist( 2 )( 0 ).set( "NLA_max_iter", tNLA_max_iter );

        tParameterlist( 3 )( 0 ) = moris::prm::create_nonlinear_solver_parameter_list();
        tParameterlist( 3 )( 0 ).set( "NLA_DofTypes", "TEMP" );

        tParameterlist( 4 )( 0 ) = moris::prm::create_time_solver_algorithm_parameter_list();
        tParameterlist( 4 )( 0 ).set( "TSA_Num_Time_Steps", tTSA_Num_Time_Steps );
        tParameterlist( 4 )( 0 ).set( "TSA_Time_Frame", tTSA_Time_Frame );

        tParameterlist( 5 )( 0 ) = moris::prm::create_time_solver_parameter_list();
        tParameterlist( 5 )( 0 ).set( "TSA_DofTypes", "TEMP" );
        tParameterlist( 5 )( 0 ).set( "TSA_Initialize_Sol_Vec", "TEMP,0.0" );
        tParameterlist( 5 )( 0 ).set( "TSA_Output_Indices", "0" );
        tParameterlist( 5 )( 0 ).set( "TSA_Output_Criteria", "Output_Criterion" );

        tParameterlist( 6 )( 0 ) = moris::prm::create_solver_warehouse_parameterlist();

        tParameterlist( 7 )( 0 ) = moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::NONE );
    }

    void
    MSIParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_msi_parameter_list();
        tParameterlist( 0 )( 0 ).set( "order_adofs_by_host", false );
    }

    void
    VISParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_vis_parameter_list();
        tParameterlist( 0 )( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName ) );
        tParameterlist( 0 )( 0 ).set( "Mesh_Type", static_cast< uint >( vis::VIS_Mesh_Type::STANDARD ) );
        tParameterlist( 0 )( 0 ).set( "Set_Names", tTotalDomain + "," + tAllInterfaces + "," + tInterfaceSingleSides );
        tParameterlist( 0 )( 0 ).set( "Field_Names",
                "TEMP,TEMP_ANALYTIC,"
                "L2_ERROR_ANALYTIC,H1_ERROR_ANALYTIC,"
                "VOLUME,"
                "SPHERE_SURFACE_AREA,ELEMENTAL_SPHERE_SURFACE_AREA,"
                "SPHERE_SURFACE_HEAT_FLUX,TOTAL_SPHERE_SURFACE_HEAT_FLUX,"
                "BLOCK_INNER_SURFACE_HEAT_FLUX,TOTAL_BLOCK_INNER_SURFACE_HEAT_FLUX,"
                "HEAT_FLUX_JUMP" );
        tParameterlist( 0 )( 0 ).set( "Field_Type",
                "NODAL,NODAL,"
                "GLOBAL,GLOBAL,"
                "GLOBAL,"
                "GLOBAL,FACETED_INT,"
                "FACETED_AVG,GLOBAL,"
                "FACETED_AVG,GLOBAL,"
                "FACETED_AVG" );
        tParameterlist( 0 )( 0 ).set( "IQI_Names",
                "IQIBulkTEMP,IQIBulkTEMPAnalytic,"
                "IQIBulkL2Error,IQIBulkH1Error,"
                "IQIVolume,"
                "IQI_Sphere_Surface, IQI_Sphere_Surface,"
                "IQI_Sphere_Outer_Surface_Heat_Flux,IQI_Sphere_Outer_Surface_Heat_Flux,"
                "IQI_Block_Inner_Surface_Heat_Flux,IQI_Block_Inner_Surface_Heat_Flux,"
                "IQI_Surface_Heat_Flux_Jump" );
        tParameterlist( 0 )( 0 ).set( "Save_Frequency", 1 );
    }

    void
    MORISGENERALParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
    }

    /* ------------------------------------------------------------------------ */
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
