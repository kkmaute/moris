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
#include "parameters.hpp"
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

    Vector< uint > tNumElemsPerDim = gDim == 2 ? Vector< uint >{ 2, 1 } : Vector< uint >{ 2, 1, 1 };
    std::string tDomainDims     = gDim == 2 ? std::to_string( sL ) + ", 1.0" : std::to_string( sL ) + ", 1.0, 1.0";
    std::string tDomainOffset   = gDim == 2 ? "0.0,  0.0" : "0.0,  0.0, 0.0";

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
        aParameterLists.set( "exodus_output_XTK_ig_mesh", true );
        aParameterLists.set( "high_to_low_dbl_side_sets", true );
        aParameterLists.set( "exodus_output_XTK_ip_mesh", true );
        aParameterLists.set( "write_basis_functions", true );
        aParameterLists.set( "write_enrichment_fields", true );
        aParameterLists.set( "write_enrichment_fields_probe_spheres", "5.0,0.0,0.0,0.0" );
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
        aParameterLists.set( "property_name", "PropDensity1" );
        aParameterLists.set( "function_parameters", tDens1 );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropCapacity1" );
        aParameterLists.set( "function_parameters", sCap1 );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropConductivity1" );
        aParameterLists.set( "function_parameters", std::to_string( sK1 ) );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropHeatLoad1" );
        aParameterLists.set( "function_parameters", std::to_string( sQ1 ) );
        aParameterLists.set( "value_function", "Func_Const" );

        // properties for 2 material

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDensity2" );
        aParameterLists.set( "function_parameters", tDens2 );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropCapacity2" );
        aParameterLists.set( "function_parameters", sCap2 );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropConductivity2" );
        aParameterLists.set( "function_parameters", std::to_string( sK2 ) );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropHeatLoad2" );
        aParameterLists.set( "function_parameters", std::to_string( sQ2 ) );
        aParameterLists.set( "value_function", "Func_Const" );

        // surface flux
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropSurfaceFlux" );
        aParameterLists.set( "function_parameters", std::to_string( sP2 ) );
        aParameterLists.set( "value_function", "Func_Const" );

        // temperature at back surface
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropImposedTemperature" );
        aParameterLists.set( "value_function", "Func_Exact_Temperature" );

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

        // exact temperature
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropExactTemperature" );
        aParameterLists.set( "value_function", "Func_Exact_Temperature" );

        // exact temperature gradient
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropExactTemperatureGradient" );
        aParameterLists.set( "value_function", "Func_Exact_TemperatureGradient" );

        //------------------------------------------------------------------------------

        // create parameter list for constitutive model - Inclusion
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMDiffusion1" );
        aParameterLists.set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists.set( "properties",
                "PropConductivity1 , Conductivity;"
                "PropDensity1      , Density;"
                "PropCapacity1     , HeatCapacity" );

        // create parameter list for constitutive model - 2 Material
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMDiffusion2" );
        aParameterLists.set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists.set( "properties",
                "PropConductivity2 , Conductivity;"
                "PropDensity2      , Density;"
                "PropCapacity2     , HeatCapacity" );

        //------------------------------------------------------------------------------

        // create parameter list for Nitsche stabilization parameter for inclusion-2 material interface
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPInterfaceNitsche" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::NITSCHE_INTERFACE );
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "leader_properties", "PropConductivity1,Material" );
        aParameterLists.set( "follower_properties", "PropConductivity2,Material" );

        // create parameter list for DBC on back surface
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPNitscheTemp" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "leader_properties", "PropConductivity1,Material" );

        if ( tUseGhost )
        {
            // create parameter list for ghost stabilization parameter for inclusion
            aParameterLists( FEM::STABILIZATION ).add_parameter_list();
            aParameterLists.set( "stabilization_name", "SPGPTemp1" );
            aParameterLists.set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
            aParameterLists.set( "function_parameters", "0.01" );
            aParameterLists.set( "leader_properties", "PropConductivity1,Material" );

            // create parameter list for ghost stabilization parameter for 2 material
            aParameterLists( FEM::STABILIZATION ).add_parameter_list();
            aParameterLists.set( "stabilization_name", "SPGPTemp2" );
            aParameterLists.set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
            aParameterLists.set( "function_parameters", "0.01" );
            aParameterLists.set( "leader_properties", "PropConductivity2,Material" );
        }

        //------------------------------------------------------------------------------
        // create IWG for inclusion - bulk diffusion
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGDiffusion1Bulk" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusion1,Diffusion" );
        aParameterLists.set( "leader_properties", "PropHeatLoad1,Load" );
        aParameterLists.set( "mesh_set_names", tPhase1 );

        // create IWG for 2 material - bulk diffusion
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGDiffusion2Bulk" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusion2,Diffusion" );
        aParameterLists.set( "leader_properties", "PropHeatLoad2,Load" );
        aParameterLists.set( "mesh_set_names", tPhase2 );

        // create IWG for Neumann boundary conditions
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGInletFlux" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_NEUMANN );
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "leader_properties", "PropSurfaceFlux,Neumann" );
        aParameterLists.set( "mesh_set_names", tFrontSurface );

        // create IWG for Dirichlet boundary conditions
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWG2SurfaceTemp" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "leader_properties", "PropImposedTemperature,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusion1,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheTemp,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", tBackSurface );

        // create parameter list for interface conditions
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGInterface12TEMP" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "follower_dof_dependencies", "TEMP" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusion1,Diffusion" );
        aParameterLists.set( "follower_constitutive_models", "CMDiffusion2,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPInterfaceNitsche,NitscheInterface" );
        aParameterLists.set( "mesh_set_names", tInterface );

        if ( tUseGhost )
        {
            // create IWG for 2 material - ghost
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "IWGGP1Temp" );
            aParameterLists.set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists.set( "dof_residual", "TEMP" );
            aParameterLists.set( "leader_dof_dependencies", "TEMP" );
            aParameterLists.set( "follower_dof_dependencies", "TEMP" );
            aParameterLists.set( "stabilization_parameters", "SPGPTemp1,GhostSP" );
            aParameterLists.set( "mesh_set_names", tPhase1Ghost );

            // create IWG for 2 material - ghost
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "IWGGP1Temp" );
            aParameterLists.set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists.set( "dof_residual", "TEMP" );
            aParameterLists.set( "leader_dof_dependencies", "TEMP" );
            aParameterLists.set( "follower_dof_dependencies", "TEMP" );
            aParameterLists.set( "stabilization_parameters", "SPGPTemp2,GhostSP" );
            aParameterLists.set( "mesh_set_names", tPhase2Ghost );
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

        // create parameter list for IQI 1
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkTEMPAnalytic" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::PROPERTY );
        aParameterLists.set( "leader_properties", "PropExactTemperature,Property" );
        aParameterLists.set( "mesh_set_names", tTotalDomain );

        // create parameter list for IQI 2
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkL2Error" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::L2_ERROR_ANALYTIC );
        aParameterLists.set( "dof_quantity", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "leader_properties", "PropExactTemperature,L2Check" );
        aParameterLists.set( "mesh_set_names", tTotalDomain );

        // create parameter list for IQI 3
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkH1Error" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::H1_ERROR_ANALYTIC );
        aParameterLists.set( "dof_quantity", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "leader_properties", "PropExactTemperatureGradient,H1Check" );
        aParameterLists.set( "mesh_set_names", tTotalDomain );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIVolume" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::VOLUME );
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "mesh_set_names", tTotalDomain );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQI_Volume_Elem_Avg" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::VOLUME );
        aParameterLists.set( "mesh_set_names", tTotalDomain );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQI_Volume_Elem_Int" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::VOLUME );
        aParameterLists.set( "mesh_set_names", tTotalDomain );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQI_Volume_Face_Avg" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::VOLUME );
        aParameterLists.set( "mesh_set_names", tInterfaceLeader );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQI_Volume_Face_Int" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::VOLUME );
        aParameterLists.set( "mesh_set_names", tInterfaceLeader );

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
        aParameterLists.set( "Set_Names", tTotalDomain + "," + tAllInterfaces );

        aParameterLists.set( "Field_Names",
                "TEMP,TEMP_ANALYTIC,"
                "L2_ERROR_ANALYTIC,H1_ERROR_ANALYTIC,"
                "VOLUME,"
                "VOL_ELEMENTAL_AVG,VOL_ELEMENTAL_INT,"
                "VOL_FACETED_AVG,VOL_FACETED_INT" );
        aParameterLists.set( "Field_Type",
                "NODAL,NODAL,"
                "GLOBAL,GLOBAL,"
                "GLOBAL,"
                "ELEMENTAL_AVG,ELEMENTAL_INT,"
                "FACETED_AVG,FACETED_INT" );
        aParameterLists.set( "IQI_Names",
                "IQIBulkTEMP,IQIBulkTEMPAnalytic,"
                "IQIBulkL2Error,IQIBulkH1Error,"
                "IQIVolume,"
                "IQI_Volume_Elem_Avg,IQI_Volume_Elem_Int,"
                "IQI_Volume_Face_Avg,IQI_Volume_Face_Int" );

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
