/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Two_Material_Bar.cpp
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

// problem dimension: 2D or 3D
extern uint gDim;

// problem dimension: 2D or 3D
extern uint gTestCaseIndex;

//---------------------------------------------------------------

#ifdef __cplusplus
extern "C" {
#endif
//------------------------------------------------------------------------------
namespace moris
{
    /* ------------------------------------------------------------------------ */
    // Mesh Set Information

    // Phase 1: back  - Material 1
    // Phase 2: front - Material 2

    std::string tPhase1 = "HMR_dummy_n_p1,HMR_dummy_c_p1";
    std::string tPhase2 = "HMR_dummy_n_p0,HMR_dummy_c_p0";

    std::string tInterface         = "dbl_iside_p0_1_p1_0";
    std::string tInterfaceLeader   = "iside_b0_1_b1_0";
    std::string tInterfaceFollower = "iside_b0_0_b1_1";

    std::string tBackSurface = "SideSet_4_n_p1";

    std::string tFrontSurface = "SideSet_2_n_p0";

    std::string tPhase1Ghost = "ghost_p1";
    std::string tPhase2Ghost = "ghost_p0";

    std::string tTotalDomain   = tPhase1 + "," + tPhase2;
    std::string tAllInterfaces = tInterface + "," + tBackSurface + "," + tFrontSurface + "," + tInterfaceLeader + "," + tInterfaceFollower;

    /* ------------------------------------------------------------------------ */
    // geometry parameters

    // general
    moris::real sL  = 10.0;    // total length
    moris::real sL1 = 4.1;
    moris::real sL2 = sL - sL1;

    /* ------------------------------------------------------------------------ */
    // boundary condition

    // prescribed temperature
    moris::real sTpre = 0.0;

    // flux at tip fac
    moris::real sP2 = 5.0;

    /* ------------------------------------------------------------------------ */
    // material parameters

    // capacity
    std::string sCap1 = "0.0";
    std::string sCap2 = "0.0";

    // density
    std::string tDens1 = "0.0";
    std::string tDens2 = "0.0";

    // conductivity
    moris::real sK1 = 1.0;
    moris::real sK2 = 0.125;

    // body flux
    moris::real sQ1 = 1.0;
    moris::real sQ2 = 2.0;

    /* ------------------------------------------------------------------------ */
    // HMR parameters

    std::string tNumElemsPerDim = gDim == 2 ? "2,   1" : "2,   1,   1";
    std::string tDomainDims     = gDim == 2 ? std::to_string( sL ) + ", 1.0" : std::to_string( sL ) + ", 1.0, 1.0";
    std::string tDomainOffset   = gDim == 2 ? "0.0,  0.0" : "0.0,  0.0, 0.0";
    std::string tDomainSidesets = gDim == 2 ? "1,2,3,4" : "1,2,3,4,5,6";

    std::string tInterpolationOrder = "1";

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
    bool tUseGhost = false;

    /* ------------------------------------------------------------------------ */
    // Output Config

    std::string tOutputFileName = "Two_Material_Bar_" + std::to_string( gTestCaseIndex ) + ".exo";

    /* ------------------------------------------------------------------------ */
    // Level set function for diamond shaped wedge

    moris::real
    Inclusion(
            const moris::Matrix< DDRMat >& aCoordinates,
            const Vector< real >&          aGeometryParameters )
    {
        // distance from sphere center
        moris::real tX = aCoordinates( 0 );

        // Compute Signed-Distance field
        moris::real tVal = sL1 - tX;

        // clean return value to return non-zero value
        return std::abs( tVal ) < tMinLevs ? tMinLevs : tVal;
    }

    /* ------------------------------------------------------------------------ */
    // Constant function for properties

    void
    Func_Const( moris::Matrix<
                        moris::DDRMat >&              aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    /* ------------------------------------------------------------------------ */
    // Exact temperature
    // see Matlab file in test src folder

    void
    Func_Exact_Temperature(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        // distance from sphere center
        moris::real tX = aFIManager->get_IP_geometry_interpolator()->valx()( 0 );

        // Update geometry parameters
        sL2 = sL - sL1;

        if ( tX < sL1 )
        {
            aPropMatrix = { { sTpre - ( sQ1 * tX * tX ) / ( 2 * sK1 ) + ( tX * ( sP2 + sL1 * sQ1 + sL2 * sQ2 ) ) / sK1 } };
        }
        else
        {
            aPropMatrix = { { ( 2 * sL1 * sK2 * sP2 - 2 * sL1 * sK1 * sP2 + 2 * sTpre * sK1 * sK2 - sL1 * sL1 * sK1 * sQ2 + sL1 * sL1 * sK2 * sQ1 - 2 * sL1 * sL2 * sK1 * sQ2 + 2 * sL1 * sL2 * sK2 * sQ2 ) / ( 2 * sK1 * sK2 ) - ( sQ2 * tX * tX ) / ( 2 * sK2 ) + ( tX * ( sP2 + sL1 * sQ2 + sL2 * sQ2 ) ) / sK2 } };
        }
    }

    /* ------------------------------------------------------------------------ */
    // Exact temperature gradients
    // see Matlab file in test src folder

    void
    Func_Exact_TemperatureGradient(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        // distance from sphere center
        moris::real tX = aFIManager->get_IP_geometry_interpolator()->valx()( 0 );

        // Update geometry parameters
        sL2 = sL - sL1;

        // set size for aPropMatrix
        aPropMatrix.set_size( gDim, 1, 0.0 );

        // spatial gradients of analytic temperature distribution
        if ( tX <= sL1 )
        {
            aPropMatrix( 0, 0 ) = ( sP2 + sL1 * sQ1 + sL2 * sQ2 ) / sK1 - ( sQ1 * tX ) / sK1;
        }
        else
        {
            aPropMatrix( 0, 0 ) = ( sP2 + sL1 * sQ2 + sL2 * sQ2 ) / sK2 - ( sQ2 * tX ) / sK2;
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
    OPTParameterList( Vector< Submodule_Parameter_Lists >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).add_parameter_list( prm::create_opt_problem_parameter_list() );

        tParameterlist( 0 ).set( "is_optimization_problem", false );
    }

    void
    HMRParameterList( Vector< Submodule_Parameter_Lists >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).add_parameter_list( prm::create_hmr_parameter_list() );

        tParameterlist( 0 ).set( "number_of_elements_per_dimension", tNumElemsPerDim );
        tParameterlist( 0 ).set( "domain_dimensions", tDomainDims );
        tParameterlist( 0 ).set( "domain_offset", tDomainOffset );
        tParameterlist( 0 ).set( "domain_sidesets", tDomainSidesets );
        tParameterlist( 0 ).set( "lagrange_output_meshes", "0" );

        tParameterlist( 0 ).set( "lagrange_orders", std::to_string( gInterpolationOrder ) );
        tParameterlist( 0 ).set( "lagrange_pattern", "0" );
        tParameterlist( 0 ).set( "bspline_orders", std::to_string( gInterpolationOrder ) );
        tParameterlist( 0 ).set( "bspline_pattern", "0" );

        tParameterlist( 0 ).set( "lagrange_to_bspline", "0" );

        tParameterlist( 0 ).set( "truncate_bsplines", 1 );
        tParameterlist( 0 ).set( "refinement_buffer", tRefineBuffer );
        tParameterlist( 0 ).set( "staircase_buffer", tRefineBuffer );
        tParameterlist( 0 ).set( "initial_refinement", "1" );
        tParameterlist( 0 ).set( "initial_refinement_pattern", "0" );

        tParameterlist( 0 ).set( "use_number_aura", 1 );

        tParameterlist( 0 ).set( "use_multigrid", 0 );
        tParameterlist( 0 ).set( "severity_level", 0 );
    }

    void
    XTKParameterList( Vector< Submodule_Parameter_Lists >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).add_parameter_list( prm::create_xtk_parameter_list() );
        tParameterlist( 0 ).set( "decompose", true );
        tParameterlist( 0 ).set( "decomposition_type", "conformal" );
        tParameterlist( 0 ).set( "enrich", true );
        tParameterlist( 0 ).set( "basis_rank", "bspline" );
        tParameterlist( 0 ).set( "enrich_mesh_indices", "0" );
        tParameterlist( 0 ).set( "ghost_stab", tUseGhost );
        tParameterlist( 0 ).set( "multigrid", false );
        tParameterlist( 0 ).set( "verbose", true );
        tParameterlist( 0 ).set( "print_enriched_ig_mesh", false );
        tParameterlist( 0 ).set( "exodus_output_XTK_ig_mesh", true );
        tParameterlist( 0 ).set( "high_to_low_dbl_side_sets", true );
        tParameterlist( 0 ).set( "exodus_output_XTK_ip_mesh", true );
        tParameterlist( 0 ).set( "write_basis_functions", true );
        tParameterlist( 0 ).set( "write_enrichment_fields", true );
        tParameterlist( 0 ).set( "write_enrichment_fields_probe_spheres", "5.0,0.0,0.0,0.0" );
    }

    void
    GENParameterList( Vector< Submodule_Parameter_Lists >& tParameterlist )
    {
        tParameterlist.resize( 3 );
        tParameterlist( 0 ).add_parameter_list( prm::create_gen_parameter_list() );


        // Geometry parameter lists
        tParameterlist( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        tParameterlist( 1 ).set( "field_function_name", "Inclusion" );
    }

    void
    FEMParameterList( Vector< Submodule_Parameter_Lists >& tParameterList )
    {
        // create a cell of cell of parameter list for fem
        tParameterList.resize( 8 );

        //------------------------------------------------------------------------------

        // properties for inclusion

        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropDensity1" );
        tParameterList( 0 ).set( "function_parameters", tDens1 );
        tParameterList( 0 ).set( "value_function", "Func_Const" );

        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropCapacity1" );
        tParameterList( 0 ).set( "function_parameters", sCap1 );
        tParameterList( 0 ).set( "value_function", "Func_Const" );

        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropConductivity1" );
        tParameterList( 0 ).set( "function_parameters", std::to_string( sK1 ) );
        tParameterList( 0 ).set( "value_function", "Func_Const" );

        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropHeatLoad1" );
        tParameterList( 0 ).set( "function_parameters", std::to_string( sQ1 ) );
        tParameterList( 0 ).set( "value_function", "Func_Const" );

        // properties for 2 material

        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropDensity2" );
        tParameterList( 0 ).set( "function_parameters", tDens2 );
        tParameterList( 0 ).set( "value_function", "Func_Const" );

        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropCapacity2" );
        tParameterList( 0 ).set( "function_parameters", sCap2 );
        tParameterList( 0 ).set( "value_function", "Func_Const" );

        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropConductivity2" );
        tParameterList( 0 ).set( "function_parameters", std::to_string( sK2 ) );
        tParameterList( 0 ).set( "value_function", "Func_Const" );

        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropHeatLoad2" );
        tParameterList( 0 ).set( "function_parameters", std::to_string( sQ2 ) );
        tParameterList( 0 ).set( "value_function", "Func_Const" );

        // surface flux
        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropSurfaceFlux" );
        tParameterList( 0 ).set( "function_parameters", std::to_string( sP2 ) );
        tParameterList( 0 ).set( "value_function", "Func_Const" );

        // temperature at back surface
        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropImposedTemperature" );
        tParameterList( 0 ).set( "value_function", "Func_Exact_Temperature" );

        // time continuity weights
        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropWeightCurrent" );
        tParameterList( 0 ).set( "function_parameters", "100.0" );
        tParameterList( 0 ).set( "value_function", "Func_Const" );

        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropWeightPrevious" );
        tParameterList( 0 ).set( "function_parameters", "100.0" );
        tParameterList( 0 ).set( "value_function", "Func_Const" );

        // initial condition
        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropInitialCondition" );
        tParameterList( 0 ).set( "value_function", "Func_Exact_Temperature" );

        // exact temperature
        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropExactTemperature" );
        tParameterList( 0 ).set( "value_function", "Func_Exact_Temperature" );

        // exact temperature gradient
        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropExactTemperatureGradient" );
        tParameterList( 0 ).set( "value_function", "Func_Exact_TemperatureGradient" );

        //------------------------------------------------------------------------------

        // create parameter list for constitutive model - Inclusion
        tParameterList( 1 ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 ).set( "constitutive_name", "CMDiffusion1" );
        tParameterList( 1 ).set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
        tParameterList( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( 1 ).set( "properties",
                "PropConductivity1 , Conductivity;"
                "PropDensity1      , Density;"
                "PropCapacity1     , HeatCapacity" );

        // create parameter list for constitutive model - 2 Material
        tParameterList( 1 ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 ).set( "constitutive_name", "CMDiffusion2" );
        tParameterList( 1 ).set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
        tParameterList( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( 1 ).set( "properties",
                "PropConductivity2 , Conductivity;"
                "PropDensity2      , Density;"
                "PropCapacity2     , HeatCapacity" );

        //------------------------------------------------------------------------------

        // create parameter list for Nitsche stabilization parameter for inclusion-2 material interface
        tParameterList( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 ).set( "stabilization_name", "SPInterfaceNitsche" );
        tParameterList( 2 ).set( "stabilization_type", fem::Stabilization_Type::NITSCHE_INTERFACE );
        tParameterList( 2 ).set( "function_parameters", "100.0" );
        tParameterList( 2 ).set( "leader_properties", "PropConductivity1,Material" );
        tParameterList( 2 ).set( "follower_properties", "PropConductivity2,Material" );

        // create parameter list for DBC on back surface
        tParameterList( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 ).set( "stabilization_name", "SPNitscheTemp" );
        tParameterList( 2 ).set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        tParameterList( 2 ).set( "function_parameters", "100.0" );
        tParameterList( 2 ).set( "leader_properties", "PropConductivity1,Material" );

        if ( tUseGhost )
        {
            // create parameter list for ghost stabilization parameter for inclusion
            tParameterList( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
            tParameterList( 2 ).set( "stabilization_name", "SPGPTemp1" );
            tParameterList( 2 ).set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
            tParameterList( 2 ).set( "function_parameters", "0.01" );
            tParameterList( 2 ).set( "leader_properties", "PropConductivity1,Material" );

            // create parameter list for ghost stabilization parameter for 2 material
            tParameterList( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
            tParameterList( 2 ).set( "stabilization_name", "SPGPTemp2" );
            tParameterList( 2 ).set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
            tParameterList( 2 ).set( "function_parameters", "0.01" );
            tParameterList( 2 ).set( "leader_properties", "PropConductivity2,Material" );
        }

        //------------------------------------------------------------------------------
        // create IWG for inclusion - bulk diffusion
        tParameterList( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( 3 ).set( "IWG_name", "IWGDiffusion1Bulk" );
        tParameterList( 3 ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        tParameterList( 3 ).set( "dof_residual", "TEMP" );
        tParameterList( 3 ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( 3 ).set( "leader_constitutive_models", "CMDiffusion1,Diffusion" );
        tParameterList( 3 ).set( "leader_properties", "PropHeatLoad1,Load" );
        tParameterList( 3 ).set( "mesh_set_names", tPhase1 );

        // create IWG for 2 material - bulk diffusion
        tParameterList( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( 3 ).set( "IWG_name", "IWGDiffusion2Bulk" );
        tParameterList( 3 ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        tParameterList( 3 ).set( "dof_residual", "TEMP" );
        tParameterList( 3 ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( 3 ).set( "leader_constitutive_models", "CMDiffusion2,Diffusion" );
        tParameterList( 3 ).set( "leader_properties", "PropHeatLoad2,Load" );
        tParameterList( 3 ).set( "mesh_set_names", tPhase2 );

        // create IWG for Neumann boundary conditions
        tParameterList( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( 3 ).set( "IWG_name", "IWGInletFlux" );
        tParameterList( 3 ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_NEUMANN );
        tParameterList( 3 ).set( "dof_residual", "TEMP" );
        tParameterList( 3 ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( 3 ).set( "leader_properties", "PropSurfaceFlux,Neumann" );
        tParameterList( 3 ).set( "mesh_set_names", tFrontSurface );

        // create IWG for Dirichlet boundary conditions
        tParameterList( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( 3 ).set( "IWG_name", "IWG2SurfaceTemp" );
        tParameterList( 3 ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        tParameterList( 3 ).set( "dof_residual", "TEMP" );
        tParameterList( 3 ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( 3 ).set( "leader_properties", "PropImposedTemperature,Dirichlet" );
        tParameterList( 3 ).set( "leader_constitutive_models", "CMDiffusion1,Diffusion" );
        tParameterList( 3 ).set( "stabilization_parameters", "SPNitscheTemp,DirichletNitsche" );
        tParameterList( 3 ).set( "mesh_set_names", tBackSurface );

        // create parameter list for interface conditions
        tParameterList( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( 3 ).set( "IWG_name", "IWGInterface12TEMP" );
        tParameterList( 3 ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE );
        tParameterList( 3 ).set( "dof_residual", "TEMP" );
        tParameterList( 3 ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( 3 ).set( "follower_dof_dependencies", "TEMP" );
        tParameterList( 3 ).set( "leader_constitutive_models", "CMDiffusion1,Diffusion" );
        tParameterList( 3 ).set( "follower_constitutive_models", "CMDiffusion2,Diffusion" );
        tParameterList( 3 ).set( "stabilization_parameters", "SPInterfaceNitsche,NitscheInterface" );
        tParameterList( 3 ).set( "mesh_set_names", tInterface );

        if ( tUseGhost )
        {
            // create IWG for 2 material - ghost
            tParameterList( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
            tParameterList( 3 ).set( "IWG_name", "IWGGP1Temp" );
            tParameterList( 3 ).set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            tParameterList( 3 ).set( "dof_residual", "TEMP" );
            tParameterList( 3 ).set( "leader_dof_dependencies", "TEMP" );
            tParameterList( 3 ).set( "follower_dof_dependencies", "TEMP" );
            tParameterList( 3 ).set( "stabilization_parameters", "SPGPTemp1,GhostSP" );
            tParameterList( 3 ).set( "mesh_set_names", tPhase1Ghost );

            // create IWG for 2 material - ghost
            tParameterList( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
            tParameterList( 3 ).set( "IWG_name", "IWGGP1Temp" );
            tParameterList( 3 ).set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            tParameterList( 3 ).set( "dof_residual", "TEMP" );
            tParameterList( 3 ).set( "leader_dof_dependencies", "TEMP" );
            tParameterList( 3 ).set( "follower_dof_dependencies", "TEMP" );
            tParameterList( 3 ).set( "stabilization_parameters", "SPGPTemp2,GhostSP" );
            tParameterList( 3 ).set( "mesh_set_names", tPhase2Ghost );
        }

        //        // create IWG for time continuity
        //        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        //        tParameterList( 3 ).set( "IWG_name",                   "IWGTimeContinuityTemp") ;
        //        tParameterList( 3 ).set( "IWG_type",                    fem::IWG_Type::TIME_CONTINUITY_DOF ) ;
        //        tParameterList( 3 ).set( "dof_residual",               "TEMP") ;
        //        tParameterList( 3 ).set( "leader_dof_dependencies",    "TEMP") ;
        //        tParameterList( 3 ).set( "leader_properties",
        //                "PropWeightCurrent   ,WeightCurrent;"
        //                "PropWeightPrevious  ,WeightPrevious;"
        //                "PropInitialCondition,InitialCondition") ;
        //        tParameterList( 3 ).set( "mesh_set_names",             tTotalDomain );
        //        tParameterList( 3 ).set( "time_continuity",            true );
        //
        //------------------------------------------------------------------------------
        // Nodal Temperature IQI
        tParameterList( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        tParameterList( 4 ).set( "IQI_name", "IQIBulkTEMP" );
        tParameterList( 4 ).set( "IQI_type", fem::IQI_Type::DOF );
        tParameterList( 4 ).set( "dof_quantity", "TEMP" );
        tParameterList( 4 ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( 4 ).set( "vectorial_field_index", 0 );
        tParameterList( 4 ).set( "mesh_set_names", tTotalDomain );

        // create parameter list for IQI 1
        tParameterList( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        tParameterList( 4 ).set( "IQI_name", "IQIBulkTEMPAnalytic" );
        tParameterList( 4 ).set( "IQI_type", fem::IQI_Type::PROPERTY );
        tParameterList( 4 ).set( "leader_properties", "PropExactTemperature,Property" );
        tParameterList( 4 ).set( "mesh_set_names", tTotalDomain );

        // create parameter list for IQI 2
        tParameterList( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        tParameterList( 4 ).set( "IQI_name", "IQIBulkL2Error" );
        tParameterList( 4 ).set( "IQI_type", fem::IQI_Type::L2_ERROR_ANALYTIC );
        tParameterList( 4 ).set( "dof_quantity", "TEMP" );
        tParameterList( 4 ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( 4 ).set( "leader_properties", "PropExactTemperature,L2Check" );
        tParameterList( 4 ).set( "mesh_set_names", tTotalDomain );

        // create parameter list for IQI 3
        tParameterList( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        tParameterList( 4 ).set( "IQI_name", "IQIBulkH1Error" );
        tParameterList( 4 ).set( "IQI_type", fem::IQI_Type::H1_ERROR_ANALYTIC );
        tParameterList( 4 ).set( "dof_quantity", "TEMP" );
        tParameterList( 4 ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( 4 ).set( "leader_properties", "PropExactTemperatureGradient,H1Check" );
        tParameterList( 4 ).set( "mesh_set_names", tTotalDomain );

        tParameterList( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        tParameterList( 4 ).set( "IQI_name", "IQIVolume" );
        tParameterList( 4 ).set( "IQI_type", fem::IQI_Type::VOLUME );
        tParameterList( 4 ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( 4 ).set( "mesh_set_names", tTotalDomain );

        tParameterList( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        tParameterList( 4 ).set( "IQI_name", "IQI_Volume_Elem_Avg" );
        tParameterList( 4 ).set( "IQI_type", fem::IQI_Type::VOLUME );
        tParameterList( 4 ).set( "mesh_set_names", tTotalDomain );

        tParameterList( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        tParameterList( 4 ).set( "IQI_name", "IQI_Volume_Elem_Int" );
        tParameterList( 4 ).set( "IQI_type", fem::IQI_Type::VOLUME );
        tParameterList( 4 ).set( "mesh_set_names", tTotalDomain );

        tParameterList( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        tParameterList( 4 ).set( "IQI_name", "IQI_Volume_Face_Avg" );
        tParameterList( 4 ).set( "IQI_type", fem::IQI_Type::VOLUME );
        tParameterList( 4 ).set( "mesh_set_names", tInterfaceLeader );

        tParameterList( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        tParameterList( 4 ).set( "IQI_name", "IQI_Volume_Face_Int" );
        tParameterList( 4 ).set( "IQI_type", fem::IQI_Type::VOLUME );
        tParameterList( 4 ).set( "mesh_set_names", tInterfaceLeader );

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        tParameterList( 5 ).add_parameter_list( prm::create_computation_parameter_list() );
    }

    void
    SOLParameterList( Vector< Submodule_Parameter_Lists >& tParameterlist )
    {
        tParameterlist.resize( 8 );


        tParameterlist( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL ) );

        tParameterlist( 1 ).add_parameter_list( moris::prm::create_linear_solver_parameter_list() );

        tParameterlist( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
        tParameterlist( 2 ).set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        tParameterlist( 2 ).set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
        tParameterlist( 2 ).set( "NLA_max_iter", tNLA_max_iter );

        tParameterlist( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
        tParameterlist( 3 ).set( "NLA_DofTypes", "TEMP" );

        tParameterlist( 4 ).add_parameter_list( moris::prm::create_time_solver_algorithm_parameter_list() );
        tParameterlist( 4 ).set( "TSA_Num_Time_Steps", tTSA_Num_Time_Steps );
        tParameterlist( 4 ).set( "TSA_Time_Frame", tTSA_Time_Frame );

        tParameterlist( 5 ).add_parameter_list( moris::prm::create_time_solver_parameter_list() );
        tParameterlist( 5 ).set( "TSA_DofTypes", "TEMP" );
        tParameterlist( 5 ).set( "TSA_Initialize_Sol_Vec", "TEMP,0.0" );
        tParameterlist( 5 ).set( "TSA_Output_Indices", "0" );
        tParameterlist( 5 ).set( "TSA_Output_Criteria", "Output_Criterion" );

        tParameterlist( 6 ).add_parameter_list( moris::prm::create_solver_warehouse_parameterlist() );

        tParameterlist( 7 ).add_parameter_list( moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::NONE ) );
    }

    void
    MSIParameterList( Vector< Submodule_Parameter_Lists >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).add_parameter_list( prm::create_msi_parameter_list() );
        tParameterlist( 0 ).set( "order_adofs_by_host", false );
    }

    void
    VISParameterList( Vector< Submodule_Parameter_Lists >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).add_parameter_list( prm::create_vis_parameter_list() );
        tParameterlist( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName ) );
        tParameterlist( 0 ).set( "Mesh_Type", vis::VIS_Mesh_Type::STANDARD );
        tParameterlist( 0 ).set( "Set_Names", tTotalDomain + "," + tAllInterfaces );

        tParameterlist( 0 ).set( "Field_Names",
                "TEMP,TEMP_ANALYTIC,"
                "L2_ERROR_ANALYTIC,H1_ERROR_ANALYTIC,"
                "VOLUME,"
                "VOL_ELEMENTAL_AVG,VOL_ELEMENTAL_INT,"
                "VOL_FACETED_AVG,VOL_FACETED_INT" );
        tParameterlist( 0 ).set( "Field_Type",
                "NODAL,NODAL,"
                "GLOBAL,GLOBAL,"
                "GLOBAL,"
                "ELEMENTAL_AVG,ELEMENTAL_INT,"
                "FACETED_AVG,FACETED_INT" );
        tParameterlist( 0 ).set( "IQI_Names",
                "IQIBulkTEMP,IQIBulkTEMPAnalytic,"
                "IQIBulkL2Error,IQIBulkH1Error,"
                "IQIVolume,"
                "IQI_Volume_Elem_Avg,IQI_Volume_Elem_Int,"
                "IQI_Volume_Face_Avg,IQI_Volume_Face_Int" );

        tParameterlist( 0 ).set( "Save_Frequency", 1 );
    }

    void
    MORISGENERALParameterList( Vector< Submodule_Parameter_Lists >& tParameterlist )
    {
    }

    /* ------------------------------------------------------------------------ */
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
