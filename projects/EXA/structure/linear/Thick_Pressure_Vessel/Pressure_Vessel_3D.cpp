/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Pressure_Vessel_3D.cpp
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
#include "fn_norm.hpp"

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
    // geometry parameters

    moris::real tOuterRad = 0.45;
    moris::real tInnerRad = 0.425;

    /* ------------------------------------------------------------------------ */
    // loading parameters

    std::string tInnerPressure = "-1.0";
    std::string tOuterPressure = "-2.0";

    std::string tInnerTemperature = "100.0";
    std::string tOuterTemperature = "200.0";

    /* ------------------------------------------------------------------------ */
    // material parameters

    std::string tDens = "1.0";

    std::string tEmod = "1.0";
    std::string tPois = "0.3";

    std::string tCTE     = "1.0e-5";
    std::string tRefTemp = "0.0";

    std::string tCond = "1.0";
    std::string tCap  = "0.0";

    /* ------------------------------------------------------------------------ */
    // Mesh Set Information

    std::string tVessel = "HMR_dummy_n_p1,HMR_dummy_c_p1";

    std::string tSupportSurfaceX = "SideSet_4_n_p1,SideSet_4_c_p1";
    std::string tSupportSurfaceY = "SideSet_1_n_p1,SideSet_1_c_p1";
    std::string tSupportSurfaceZ = "SideSet_5_n_p1,SideSet_5_c_p1";

    //    std::string tVessel         = "HMR_dummy_c_p2";
    //
    //    std::string tSupportSurfaceX = "SideSet_4_c_p2";
    //    std::string tSupportSurfaceY = "SideSet_1_c_p2";
    //    std::string tSupportSurfaceZ = "SideSet_5_c_p2";

    std::string tInnerPressureSurface = "iside_b0_1_b1_0";
    std::string tOuterPressureSurface = "iside_b0_1_b1_3";

    std::string tVesselGhost = "ghost_p1";

    /* ------------------------------------------------------------------------ */
    // HMR parameters

    std::string tNumElemsPerDim = "20, 20, 20";
    std::string tDomainDims     = "0.5, 0.5, 0.5";
    std::string tDomainOffset   = "0, 0, 0";
    std::string tDomainSidesets = "1,2,3,4,5,6";

    int tRefineBuffer = 1;

    /* ------------------------------------------------------------------------ */
    // Minimum level set value
    moris::real tMinLevs = 1.0e-8;

    /* ------------------------------------------------------------------------ */
    // Flag for turning on/off ghost stabilization
    bool tUseGhost = true;

    /* ------------------------------------------------------------------------ */
    // Output Config

    std::string tOutputFileName = "Pressure_Vessel_3D.exo";

    /* ------------------------------------------------------------------------ */

    void
    Func_Select_X(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        aPropMatrix.set_size( 3, 3, 0.0 );
        aPropMatrix( 0, 0 ) = 1.0;
    }

    /* ------------------------------------------------------------------------ */

    void
    Func_Select_Y(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        aPropMatrix.set_size( 3, 3, 0.0 );
        aPropMatrix( 1, 1 ) = 1.0;
    }

    /* ------------------------------------------------------------------------ */

    void
    Func_Select_Z(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        aPropMatrix.set_size( 3, 3, 0.0 );
        aPropMatrix( 2, 2 ) = 1.0;
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

    /* ------------------------------------------------------------------------ */

    bool
    Output_Criterion( moris::tsa::Time_Solver* aTimeSolver )
    {
        return true;
    }

    /* ------------------------------------------------------------------------ */

    void
    OPTParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_opt_problem_parameter_list();

        tParameterlist( 0 )( 0 ).set( "is_optimization_problem", false );
    }

    /* ------------------------------------------------------------------------ */

    void
    HMRParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
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
        tParameterlist( 0 )( 0 ).set( "initial_refinement", "0" );
        tParameterlist( 0 )( 0 ).set( "initial_refinement_pattern", "0" );

        tParameterlist( 0 )( 0 ).set( "use_number_aura", 1 );

        tParameterlist( 0 )( 0 ).set( "use_multigrid", 0 );
        tParameterlist( 0 )( 0 ).set( "severity_level", 0 );
    }

    /* ------------------------------------------------------------------------ */

    void
    XTKParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
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
        tParameterlist( 0 )( 0 ).set( "verbose", true );
        tParameterlist( 0 )( 0 ).set( "print_enriched_ig_mesh", false );
        tParameterlist( 0 )( 0 ).set( "exodus_output_XTK_ig_mesh", false );
        tParameterlist( 0 )( 0 ).set( "high_to_low_dbl_side_sets", true );
    }

    /* ------------------------------------------------------------------------ */

    void
    GENParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
    {
        tParameterlist.resize( 3 );
        tParameterlist( 0 ).resize( 1 );

        // Main GEN parameter list
        tParameterlist( 0 )( 0 ) = prm::create_gen_parameter_list();

        // init geometry counter
        uint tGeoCounter = 0;

        // Geometry parameter lists
        tParameterlist( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::SPHERE ) );
        tParameterlist( 1 )( tGeoCounter ).set( "radius", tOuterRad );
        tGeoCounter++;

        tParameterlist( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::SPHERE ) );
        tParameterlist( 1 )( tGeoCounter ).set( "radius", tInnerRad );
    }
    /* ------------------------------------------------------------------------ */

    void
    FEMParameterList( Vector< Vector< Parameter_List > >& tParameterList )
    {
        // create a cell of cell of parameter list for fem
        tParameterList.resize( 8 );

        //------------------------------------------------------------------------------
        // init property counter
        uint tPropCounter = 0;

        // properties of bars
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ) = prm::create_property_parameter_list();
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropDensity" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tDens );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropYoungs" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tEmod );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropPoisson" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tPois );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropCTE" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tCTE );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropRefTemp" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tRefTemp );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropConductivity" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tCond );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropCapacity" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tCap );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // properties of boundary conditions
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropDirichlet" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "0.0;0.0;0.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropSelectX" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "0.0;0.0;0.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Select_X" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropSelectY" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "0.0;0.0;0.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Select_Y" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropSelectZ" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "0.0;0.0;0.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Select_Z" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropInnerPressureLoad" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tInnerPressure );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropOuterPressureLoad" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tOuterPressure );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropInnerTemperature" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tInnerTemperature );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropOuterTemperature" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tOuterTemperature );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        //------------------------------------------------------------------------------
        // init CM counter
        uint tCMCounter = 0;

        // create parameter list for constitutive model 1
        tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_name", "CMStrucLinIso1" );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_type", fem::Constitutive_Type::STRUC_LIN_ISO );
        tParameterList( 1 )( tCMCounter ).set( "model_type", fem::Model_Type::FULL );
        tParameterList( 1 )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY,UZ;TEMP", "Displacement,Temperature" ) );
        tParameterList( 1 )( tCMCounter ).set( "properties",
                "PropYoungs, YoungsModulus;"
                "PropPoisson,PoissonRatio;"
                "PropCTE,    CTE;"
                "PropRefTemp,ReferenceTemperature" );
        tCMCounter++;

        // create parameter list for constitutive model - Inclusion
        tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_name", "CMDiffusion" );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
        tParameterList( 1 )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( 1 )( tCMCounter ).set( "properties",
                "PropConductivity , Conductivity;"
                "PropDensity      , Density;"
                "PropCapacity     , HeatCapacity" );
        tCMCounter++;

        //------------------------------------------------------------------------------
        // init SP counter
        uint tSPCounter = 0;

        // Nitsche stabilization parameter for structure
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", "SPNitscheStruc" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", "100.0" );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", "PropYoungs,Material" );
        tSPCounter++;

        // Nitsche stabilization parameter for thermal problem
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", "SPNitscheTemp" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", "100.0" );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", "PropConductivity,Material" );
        tSPCounter++;

        // Ghost stabilization parameter for structure
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", "SPGhostStruct" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", "0.01" );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", "PropYoungs,Material" );
        tSPCounter++;

        // Ghost stabilization parameter for thermal problem
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", "SPGhostTemp" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", "0.01" );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", "PropConductivity,Material" );
        tSPCounter++;

        //------------------------------------------------------------------------------
        // init IWG counter
        uint tIWGCounter = 0;

        // create IWG - bulk structure
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGBulkStruct" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_BULK );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "UX,UY,UZ" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY,UZ" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tVessel );
        tIWGCounter++;

        // create IWG - bulk diffusion
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGBulkTemp" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMDiffusion,Diffusion" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tVessel );
        tIWGCounter++;

        // create IWG - Dirichlet structure
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGDirichletX" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "UX,UY,UZ" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY,UZ" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties", "PropDirichlet,Dirichlet;PropSelectX,Select" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPNitscheStruc,DirichletNitsche" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tSupportSurfaceX );
        tIWGCounter++;

        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGDirichletY" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "UX,UY,UZ" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY,UZ" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties", "PropDirichlet,Dirichlet;PropSelectY,Select" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPNitscheStruc,DirichletNitsche" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tSupportSurfaceY );
        tIWGCounter++;

        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGDirichletZ" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "UX,UY,UZ" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY,UZ" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties", "PropDirichlet,Dirichlet;PropSelectZ,Select" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPNitscheStruc,DirichletNitsche" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tSupportSurfaceZ );
        tIWGCounter++;

        // create IWG - Dirichlet temp
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGDirichletTempInner" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties", "PropInnerTemperature,Dirichlet" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMDiffusion,Diffusion" );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPNitscheTemp,DirichletNitsche" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tInnerPressureSurface );
        tIWGCounter++;

        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGDirichletTempOuter" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties", "PropOuterTemperature,Dirichlet" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMDiffusion,Diffusion" );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPNitscheTemp,DirichletNitsche" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tOuterPressureSurface );
        tIWGCounter++;

        // create IWG - Neumann structure
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGNeumannInnerPressure" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_NEUMANN );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "UX,UY,UZ" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY,UZ" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties", "PropInnerPressureLoad,Pressure" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tInnerPressureSurface );
        tIWGCounter++;

        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGNeumannOuterPressure" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_NEUMANN );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "UX,UY,UZ" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY,UZ" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties", "PropOuterPressureLoad,Pressure" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tOuterPressureSurface );
        tIWGCounter++;

        if ( tUseGhost )
        {
            // create IWG - ghost structure
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGGhostStructure" );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "UX,UY,UZ" );
            tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY,UZ" );
            tParameterList( 3 )( tIWGCounter ).set( "follower_dof_dependencies", "UX,UY,UZ" );
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPGhostStruct,GhostSP" );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tVesselGhost );
            tIWGCounter++;

            // create IWG - ghost temp
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGGhostTemp" );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "TEMP" );
            tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "TEMP" );
            tParameterList( 3 )( tIWGCounter ).set( "follower_dof_dependencies", "TEMP" );
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPGhostTemp,GhostSP" );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tVesselGhost );
            tIWGCounter++;
        }

        //------------------------------------------------------------------------------
        // init IQI counter
        uint tIQICounter = 0;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkDISPX" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", fem::IQI_Type::DOF );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "UX,UY,UZ" );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "UX,UY,UZ" );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 0 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tVessel );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkDISPY" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", fem::IQI_Type::DOF );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "UX,UY,UZ" );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "UX,UY,UZ" );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 1 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tVessel );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkDISPZ" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", fem::IQI_Type::DOF );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "UX,UY,UZ" );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "UX,UY,UZ" );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 2 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tVessel );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkTEMP" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", fem::IQI_Type::DOF );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "TEMP" );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 0 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tVessel );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkStrainEnergy" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", fem::IQI_Type::STRAIN_ENERGY );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "UX,UY,UZ" );
        tParameterList( 4 )( tIQICounter ).set( "leader_constitutive_models", "CMStrucLinIso1,Elast" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tVessel );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkVolume" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", fem::IQI_Type::VOLUME );
        tParameterList( 4 )( tIQICounter ).set( "leader_properties", "PropDensity,Density" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tVessel );
        tIQICounter++;

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        tParameterList( 5 ).resize( 1 );
        tParameterList( 5 )( 0 ) = prm::create_computation_parameter_list();
    }

    /* ------------------------------------------------------------------------ */

    void
    SOLParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
    {
        tParameterlist.resize( 8 );

        tParameterlist( 0 ).push_back( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL ) );

        tParameterlist( 1 ).push_back( moris::prm::create_linear_solver_parameter_list() );

        tParameterlist( 2 ).push_back( moris::prm::create_nonlinear_algorithm_parameter_list() );
        tParameterlist( 2 )( 0 ).set( "NLA_combined_res_jac_assembly", true );

        tParameterlist( 3 ).push_back( moris::prm::create_nonlinear_solver_parameter_list() );
        tParameterlist( 3 )( 0 ).set( "NLA_DofTypes", "UX,UY,UZ;TEMP" );

        tParameterlist( 4 ).push_back( moris::prm::create_time_solver_algorithm_parameter_list() );
        tParameterlist( 4 )( 0 ).set( "TSA_Num_Time_Steps", 1 );
        tParameterlist( 4 )( 0 ).set( "TSA_Time_Frame", 1.0 );

        tParameterlist( 5 ).push_back( moris::prm::create_time_solver_parameter_list() );
        tParameterlist( 5 )( 0 ).set( "TSA_DofTypes", "UX,UY,UZ;TEMP" );
        tParameterlist( 5 )( 0 ).set( "TSA_Output_Indices", "0" );
        tParameterlist( 5 )( 0 ).set( "TSA_Output_Criteria", "Output_Criterion" );

        tParameterlist( 6 ).push_back( moris::prm::create_solver_warehouse_parameterlist() );

        tParameterlist( 7 ).push_back( moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::NONE ) );
    }

    /* ------------------------------------------------------------------------ */

    void
    MSIParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_msi_parameter_list();
    }

    /* ------------------------------------------------------------------------ */

    void
    VISParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_vis_parameter_list();
        tParameterlist( 0 )( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName ) );
        tParameterlist( 0 )( 0 ).set( "Mesh_Type", vis::VIS_Mesh_Type::STANDARD );
        tParameterlist( 0 )( 0 ).set( "Set_Names", tVessel );
        tParameterlist( 0 )( 0 ).set( "Field_Names", "UX,UY,UZ,TEMP,STRAIN_ENERGY,VOLUME" );
        tParameterlist( 0 )( 0 ).set( "Field_Type", "NODAL,NODAL,NODAL,NODAL,GLOBAL,GLOBAL" );
        tParameterlist( 0 )( 0 ).set( "IQI_Names", "IQIBulkDISPX,IQIBulkDISPY,IQIBulkDISPZ,IQIBulkTEMP,IQIBulkStrainEnergy,IQIBulkVolume" );
        tParameterlist( 0 )( 0 ).set( "Save_Frequency", 1 );
    }

    void
    MORISGENERALParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
    {
    }

    /* ------------------------------------------------------------------------ */
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
