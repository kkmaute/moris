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
#include "parameters.hpp"
#include "fn_equal_to.hpp"
#include "fn_norm.hpp"

#include "AztecOO.h"

//---------------------------------------------------------------

// global variable for interpolation order
extern uint gInterpolationOrder;

// global variable for case index
extern uint gCaseIndex;

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

    std::string tOutputFileName =
            "Pressure_Vessel_3D_Case_" + std::to_string( gCaseIndex ) + ".exo";

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
    OPTParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_opt_problem_parameter_list() );

        aParameterLists( 0 ).set( "is_optimization_problem", false );
    }

    /* ------------------------------------------------------------------------ */

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
        aParameterLists( 0 ).set( "initial_refinement", "0" );
        aParameterLists( 0 ).set( "initial_refinement_pattern", "0" );

        aParameterLists( 0 ).set( "use_number_aura", 1 );

        aParameterLists( 0 ).set( "use_multigrid", 0 );
        aParameterLists( 0 ).set( "severity_level", 0 );
    }

    /* ------------------------------------------------------------------------ */

    void
    XTKParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_xtk_parameter_list() );
        aParameterLists( 0 ).set( "decompose", true );
        aParameterLists( 0 ).set( "decomposition_type", "conformal" );
        aParameterLists( 0 ).set( "enrich", true );
        aParameterLists( 0 ).set( "basis_rank", "bspline" );
        aParameterLists( 0 ).set( "enrich_mesh_indices", "0" );
        aParameterLists( 0 ).set( "ghost_stab", true );
        aParameterLists( 0 ).set( "multigrid", false );
        aParameterLists( 0 ).set( "verbose", true );
        aParameterLists( 0 ).set( "print_enriched_ig_mesh", false );
        aParameterLists( 0 ).set( "exodus_output_XTK_ig_mesh", false );
        aParameterLists( 0 ).set( "high_to_low_dbl_side_sets", true );
    }

    /* ------------------------------------------------------------------------ */

    void
    GENParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_gen_parameter_list() );

        // Geometry parameter lists
        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::SPHERE ) );
        aParameterLists( 1 ).set( "radius", tOuterRad );

        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::SPHERE ) );
        aParameterLists( 1 ).set( "radius", tInnerRad );
    }
    /* ------------------------------------------------------------------------ */

    void
    FEMParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.hack_for_legacy_fem();
        // create a cell of cell of parameter list for fem

        //------------------------------------------------------------------------------

        // properties of bars
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropDensity" );
        aParameterLists( 0 ).set( "function_parameters", tDens );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropYoungs" );
        aParameterLists( 0 ).set( "function_parameters", tEmod );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropPoisson" );
        aParameterLists( 0 ).set( "function_parameters", tPois );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropCTE" );
        aParameterLists( 0 ).set( "function_parameters", tCTE );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropRefTemp" );
        aParameterLists( 0 ).set( "function_parameters", tRefTemp );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropConductivity" );
        aParameterLists( 0 ).set( "function_parameters", tCond );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropCapacity" );
        aParameterLists( 0 ).set( "function_parameters", tCap );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // properties of boundary conditions
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropDirichlet" );
        aParameterLists( 0 ).set( "function_parameters", "0.0;0.0;0.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropSelectX" );
        aParameterLists( 0 ).set( "function_parameters", "0.0;0.0;0.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Select_X" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropSelectY" );
        aParameterLists( 0 ).set( "function_parameters", "0.0;0.0;0.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Select_Y" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropSelectZ" );
        aParameterLists( 0 ).set( "function_parameters", "0.0;0.0;0.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Select_Z" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropInnerPressureLoad" );
        aParameterLists( 0 ).set( "function_parameters", tInnerPressure );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropOuterPressureLoad" );
        aParameterLists( 0 ).set( "function_parameters", tOuterPressure );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropInnerTemperature" );
        aParameterLists( 0 ).set( "function_parameters", tInnerTemperature );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropOuterTemperature" );
        aParameterLists( 0 ).set( "function_parameters", tOuterTemperature );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        //------------------------------------------------------------------------------

        // create parameter list for constitutive model 1
        aParameterLists( 1 ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists( 1 ).set( "constitutive_name", "CMStrucLinIso1" );
        aParameterLists( 1 ).set( "constitutive_type", fem::Constitutive_Type::STRUC_LIN_ISO );
        aParameterLists( 1 ).set( "model_type", fem::Model_Type::FULL );
        aParameterLists( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY,UZ;TEMP", "Displacement,Temperature" ) );
        aParameterLists( 1 ).set( "properties",
                "PropYoungs, YoungsModulus;"
                "PropPoisson,PoissonRatio;"
                "PropCTE,    CTE;"
                "PropRefTemp,ReferenceTemperature" );

        // create parameter list for constitutive model - Inclusion
        aParameterLists( 1 ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists( 1 ).set( "constitutive_name", "CMDiffusion" );
        aParameterLists( 1 ).set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
        aParameterLists( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists( 1 ).set( "properties",
                "PropConductivity , Conductivity;"
                "PropDensity      , Density;"
                "PropCapacity     , HeatCapacity" );

        //------------------------------------------------------------------------------

        // Nitsche stabilization parameter for structure
        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPNitscheStruc" );
        aParameterLists( 2 ).set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        aParameterLists( 2 ).set( "function_parameters", "100.0" );
        aParameterLists( 2 ).set( "leader_properties", "PropYoungs,Material" );

        // Nitsche stabilization parameter for thermal problem
        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPNitscheTemp" );
        aParameterLists( 2 ).set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        aParameterLists( 2 ).set( "function_parameters", "100.0" );
        aParameterLists( 2 ).set( "leader_properties", "PropConductivity,Material" );

        // Ghost stabilization parameter for structure
        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPGhostStruct" );
        aParameterLists( 2 ).set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        aParameterLists( 2 ).set( "function_parameters", "0.01" );
        aParameterLists( 2 ).set( "leader_properties", "PropYoungs,Material" );

        // Ghost stabilization parameter for thermal problem
        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPGhostTemp" );
        aParameterLists( 2 ).set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        aParameterLists( 2 ).set( "function_parameters", "0.01" );
        aParameterLists( 2 ).set( "leader_properties", "PropConductivity,Material" );

        //------------------------------------------------------------------------------
        // create IWG - bulk structure
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGBulkStruct" );
        aParameterLists( 3 ).set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_BULK );
        aParameterLists( 3 ).set( "dof_residual", "UX,UY,UZ" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY,UZ" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        aParameterLists( 3 ).set( "mesh_set_names", tVessel );

        // create IWG - bulk diffusion
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGBulkTemp" );
        aParameterLists( 3 ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        aParameterLists( 3 ).set( "dof_residual", "TEMP" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMDiffusion,Diffusion" );
        aParameterLists( 3 ).set( "mesh_set_names", tVessel );

        // create IWG - Dirichlet structure
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGDirichletX" );
        aParameterLists( 3 ).set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists( 3 ).set( "dof_residual", "UX,UY,UZ" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY,UZ" );
        aParameterLists( 3 ).set( "leader_properties", "PropDirichlet,Dirichlet;PropSelectX,Select" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPNitscheStruc,DirichletNitsche" );
        aParameterLists( 3 ).set( "mesh_set_names", tSupportSurfaceX );

        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGDirichletY" );
        aParameterLists( 3 ).set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists( 3 ).set( "dof_residual", "UX,UY,UZ" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY,UZ" );
        aParameterLists( 3 ).set( "leader_properties", "PropDirichlet,Dirichlet;PropSelectY,Select" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPNitscheStruc,DirichletNitsche" );
        aParameterLists( 3 ).set( "mesh_set_names", tSupportSurfaceY );

        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGDirichletZ" );
        aParameterLists( 3 ).set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists( 3 ).set( "dof_residual", "UX,UY,UZ" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY,UZ" );
        aParameterLists( 3 ).set( "leader_properties", "PropDirichlet,Dirichlet;PropSelectZ,Select" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPNitscheStruc,DirichletNitsche" );
        aParameterLists( 3 ).set( "mesh_set_names", tSupportSurfaceZ );

        // create IWG - Dirichlet temp
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGDirichletTempInner" );
        aParameterLists( 3 ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists( 3 ).set( "dof_residual", "TEMP" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 3 ).set( "leader_properties", "PropInnerTemperature,Dirichlet" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMDiffusion,Diffusion" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPNitscheTemp,DirichletNitsche" );
        aParameterLists( 3 ).set( "mesh_set_names", tInnerPressureSurface );

        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGDirichletTempOuter" );
        aParameterLists( 3 ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists( 3 ).set( "dof_residual", "TEMP" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 3 ).set( "leader_properties", "PropOuterTemperature,Dirichlet" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMDiffusion,Diffusion" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPNitscheTemp,DirichletNitsche" );
        aParameterLists( 3 ).set( "mesh_set_names", tOuterPressureSurface );

        // create IWG - Neumann structure
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGNeumannInnerPressure" );
        aParameterLists( 3 ).set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_NEUMANN );
        aParameterLists( 3 ).set( "dof_residual", "UX,UY,UZ" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY,UZ" );
        aParameterLists( 3 ).set( "leader_properties", "PropInnerPressureLoad,Pressure" );
        aParameterLists( 3 ).set( "mesh_set_names", tInnerPressureSurface );

        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGNeumannOuterPressure" );
        aParameterLists( 3 ).set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_NEUMANN );
        aParameterLists( 3 ).set( "dof_residual", "UX,UY,UZ" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY,UZ" );
        aParameterLists( 3 ).set( "leader_properties", "PropOuterPressureLoad,Pressure" );
        aParameterLists( 3 ).set( "mesh_set_names", tOuterPressureSurface );

        if ( tUseGhost )
        {
            // create IWG - ghost structure
            aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists( 3 ).set( "IWG_name", "IWGGhostStructure" );
            aParameterLists( 3 ).set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists( 3 ).set( "dof_residual", "UX,UY,UZ" );
            aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY,UZ" );
            aParameterLists( 3 ).set( "follower_dof_dependencies", "UX,UY,UZ" );
            aParameterLists( 3 ).set( "stabilization_parameters", "SPGhostStruct,GhostSP" );
            aParameterLists( 3 ).set( "mesh_set_names", tVesselGhost );

            // create IWG - ghost temp
            aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists( 3 ).set( "IWG_name", "IWGGhostTemp" );
            aParameterLists( 3 ).set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists( 3 ).set( "dof_residual", "TEMP" );
            aParameterLists( 3 ).set( "leader_dof_dependencies", "TEMP" );
            aParameterLists( 3 ).set( "follower_dof_dependencies", "TEMP" );
            aParameterLists( 3 ).set( "stabilization_parameters", "SPGhostTemp,GhostSP" );
            aParameterLists( 3 ).set( "mesh_set_names", tVesselGhost );
        }

        //------------------------------------------------------------------------------
        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkDISPX" );
        aParameterLists( 4 ).set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists( 4 ).set( "dof_quantity", "UX,UY,UZ" );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "UX,UY,UZ" );
        aParameterLists( 4 ).set( "vectorial_field_index", 0 );
        aParameterLists( 4 ).set( "mesh_set_names", tVessel );

        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkDISPY" );
        aParameterLists( 4 ).set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists( 4 ).set( "dof_quantity", "UX,UY,UZ" );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "UX,UY,UZ" );
        aParameterLists( 4 ).set( "vectorial_field_index", 1 );
        aParameterLists( 4 ).set( "mesh_set_names", tVessel );

        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkDISPZ" );
        aParameterLists( 4 ).set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists( 4 ).set( "dof_quantity", "UX,UY,UZ" );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "UX,UY,UZ" );
        aParameterLists( 4 ).set( "vectorial_field_index", 2 );
        aParameterLists( 4 ).set( "mesh_set_names", tVessel );

        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkTEMP" );
        aParameterLists( 4 ).set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists( 4 ).set( "dof_quantity", "TEMP" );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 4 ).set( "vectorial_field_index", 0 );
        aParameterLists( 4 ).set( "mesh_set_names", tVessel );

        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkStrainEnergy" );
        aParameterLists( 4 ).set( "IQI_type", fem::IQI_Type::STRAIN_ENERGY );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "UX,UY,UZ" );
        aParameterLists( 4 ).set( "leader_constitutive_models", "CMStrucLinIso1,Elast" );
        aParameterLists( 4 ).set( "mesh_set_names", tVessel );

        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkVolume" );
        aParameterLists( 4 ).set( "IQI_type", fem::IQI_Type::VOLUME );
        aParameterLists( 4 ).set( "leader_properties", "PropDensity,Density" );
        aParameterLists( 4 ).set( "mesh_set_names", tVessel );

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        aParameterLists( 5 ).add_parameter_list( prm::create_computation_parameter_list() );
    }

    /* ------------------------------------------------------------------------ */

    void
    SOLParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL ) );

        aParameterLists( 1 ).add_parameter_list( moris::prm::create_linear_solver_parameter_list() );

        aParameterLists( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
        aParameterLists( 2 ).set( "NLA_combined_res_jac_assembly", true );

        aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
        aParameterLists( 3 ).set( "NLA_DofTypes", "UX,UY,UZ;TEMP" );

        aParameterLists( 4 ).add_parameter_list( moris::prm::create_time_solver_algorithm_parameter_list() );
        aParameterLists( 4 ).set( "TSA_Num_Time_Steps", 1 );
        aParameterLists( 4 ).set( "TSA_Time_Frame", 1.0 );

        aParameterLists( 5 ).add_parameter_list( moris::prm::create_time_solver_parameter_list() );
        aParameterLists( 5 ).set( "TSA_DofTypes", "UX,UY,UZ;TEMP" );
        aParameterLists( 5 ).set( "TSA_Output_Indices", "0" );
        aParameterLists( 5 ).set( "TSA_Output_Criteria", "Output_Criterion" );

        aParameterLists( 6 ).add_parameter_list( moris::prm::create_solver_warehouse_parameterlist() );

        aParameterLists( 7 ).add_parameter_list( moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::NONE ) );
    }

    /* ------------------------------------------------------------------------ */

    void
    MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_msi_parameter_list() );
    }

    /* ------------------------------------------------------------------------ */

    void
    VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_vis_parameter_list() );
        aParameterLists( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName ) );
        aParameterLists( 0 ).set( "Mesh_Type", vis::VIS_Mesh_Type::STANDARD );
        aParameterLists( 0 ).set( "Set_Names", tVessel );
        aParameterLists( 0 ).set( "Field_Names", "UX,UY,UZ,TEMP,STRAIN_ENERGY,VOLUME" );
        aParameterLists( 0 ).set( "Field_Type", "NODAL,NODAL,NODAL,NODAL,GLOBAL,GLOBAL" );
        aParameterLists( 0 ).set( "IQI_Names", "IQIBulkDISPX,IQIBulkDISPY,IQIBulkDISPZ,IQIBulkTEMP,IQIBulkStrainEnergy,IQIBulkVolume" );
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
