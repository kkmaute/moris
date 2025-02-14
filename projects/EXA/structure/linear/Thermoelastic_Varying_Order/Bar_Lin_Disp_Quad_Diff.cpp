/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Bar_Lin_Disp_Quad_Diff.cpp
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

//---------------------------------------------------------------

#ifdef __cplusplus
extern "C" {
#endif
//------------------------------------------------------------------------------
namespace moris
{
    // material parameters

    std::string tDens = "1.0";

    std::string tEmod = "1.0";
    std::string tPois = "0.0";

    std::string tCTE     = "0.0e-5";
    std::string tRefTemp = "0.0";

    std::string tCond = "1.0";
    std::string tCap  = "0.0";

    /* ------------------------------------------------------------------------ */
    // Mesh Set Information

    std::string tDomain1 = "HMR_dummy_n_p0,HMR_dummy_c_p0";
    std::string tDomain2 = "HMR_dummy_n_p1,HMR_dummy_c_p1";
    std::string tDomain  = "HMR_dummy_n_p0,HMR_dummy_c_p0,HMR_dummy_n_p1,HMR_dummy_c_p1";

    std::string tLeftSurf = "SideSet_4_n_p0";

    std::string tInterface = "dbl_iside_p0_0_p1_1";

    /* ------------------------------------------------------------------------ */
    // HMR parameters

    Vector< uint > tNumElemsPerDim = { 4, 2, 2 };
    Vector< real > tDomainDims     = { 2.0, 0.25, 0.25 };

    int tRefineBuffer = 0;

    /* ------------------------------------------------------------------------ */
    // Flag for turning on/off ghost stabilization
    bool tUseGhost = false;

    /* ------------------------------------------------------------------------ */
    // Output Config

    std::string tOutputFileName = "Bar_Lin_Disp_Quad_Diff.exo";

    /* ------------------------------------------------------------------------ */
    // Level set function of Sphere

    moris::real Plane(
            const moris::Matrix< DDRMat >& aCoordinates,
            const Vector< real >&     aGeometryParameters )
    {
        return aCoordinates( 0 ) - 1.2;
    }

    /* ------------------------------------------------------------------------ */

    // Constant function for properties
    void
    Func_Const(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    /* ------------------------------------------------------------------------ */

    void
    Func_Body_Load( moris::Matrix< moris::DDRMat >&        aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        moris::Matrix< DDRMat > tMat = { { 1 }, { 0 }, { 0 } };

        // moris::real tX = aFIManager->get_IG_geometry_interpolator()->valx()( 0 );
        // moris::real tValue = 3.0*std::pow(tX, 2)-tX+3;
        // aPropMatrix = aParameters( 0 )( 0 ) * tValue * tMat;

        aPropMatrix = 1.0 * tMat;
    }

    /* ------------------------------------------------------------------------ */
    void
    Func_Temp_Body_Load( moris::Matrix< moris::DDRMat >&   aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix.set_size( 1, 1, 1.0 );
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

    /* ------------------------------------------------------------------------ */

    void
    HMRParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "number_of_elements_per_dimension", tNumElemsPerDim );
        aParameterLists.set( "domain_dimensions", tDomainDims );
        aParameterLists.set( "lagrange_output_meshes", "0" );

        aParameterLists.set( "lagrange_orders", "2" );
        aParameterLists.set( "lagrange_pattern", "0" );
        aParameterLists.set( "bspline_orders", "2,1" );
        aParameterLists.set( "bspline_pattern", "0,1" );

        aParameterLists.set( "lagrange_to_bspline", "0,1" );

        aParameterLists.set( "refinement_buffer", tRefineBuffer );
        aParameterLists.set( "staircase_buffer", tRefineBuffer );
        aParameterLists.set( "initial_refinement", "0" );
    }

    /* ------------------------------------------------------------------------ */

    void
    XTKParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "decompose", true );
        aParameterLists.set( "decomposition_type", "conformal" );
        aParameterLists.set( "enrich", true );
        aParameterLists.set( "basis_rank", "bspline" );
        aParameterLists.set( "enrich_mesh_indices", "0,1" );
        aParameterLists.set( "ghost_stab", true );
        aParameterLists.set( "multigrid", false );
        aParameterLists.set( "verbose", true );
        aParameterLists.set( "print_enriched_ig_mesh", false );
        aParameterLists.set( "exodus_output_XTK_ig_mesh", true );
        aParameterLists.set( "high_to_low_dbl_side_sets", true );
    }

    /* ------------------------------------------------------------------------ */

    void
    GENParameterList( Module_Parameter_Lists& aParameterLists )
    {
        // Geometry parameter lists
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists.set( "field_function_name", "Plane" );
    }
    /* ------------------------------------------------------------------------ */

    void
    FEMParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.hack_for_legacy_fem();
        // create a cell of cell of parameter list for fem

        //------------------------------------------------------------------------------

        // properties of bars
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDensity" );
        aParameterLists.set( "function_parameters", tDens );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropYoungs" );
        aParameterLists.set( "function_parameters", tEmod );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropPoisson" );
        aParameterLists.set( "function_parameters", tPois );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropConductivity" );
        aParameterLists.set( "function_parameters", tCond );
        aParameterLists.set( "value_function", "Func_Const" );

        // properties of boundary conditions
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDirichlet" );
        aParameterLists.set( "function_parameters", "0.0;0.0;0.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropTemperature" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDispNeumann" );
        aParameterLists.set( "function_parameters", "1.0;0.0;0.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropTempNeumann" );
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropBodyLoad" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Body_Load" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropTempBodyLoad" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Temp_Body_Load" );

        //------------------------------------------------------------------------------

        // create parameter list for constitutive model 1
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMStrucLinIso1" );
        aParameterLists.set( "constitutive_type",  fem::Constitutive_Type::STRUC_LIN_ISO ) ;
        aParameterLists.set( "model_type",  fem::Model_Type::FULL ) ;
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY,UZ", "Displacement,Temperature" ) );
        aParameterLists.set( "properties",
                "PropYoungs, YoungsModulus;"
                "PropPoisson,PoissonRatio" );

        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMStrucLinIso2" );
        aParameterLists.set( "constitutive_type",  fem::Constitutive_Type::STRUC_LIN_ISO ) ;
        aParameterLists.set( "model_type",  fem::Model_Type::FULL ) ;
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY,UZ", "Displacement,Temperature" ) );
        aParameterLists.set( "properties",
                "PropYoungs, YoungsModulus;"
                "PropPoisson,PoissonRatio" );

        // create parameter list for constitutive model - Inclusion
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMDiffusion1" );
        aParameterLists.set( "constitutive_type",  fem::Constitutive_Type::DIFF_LIN_ISO ) ;
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists.set( "properties", "PropConductivity , Conductivity" );

        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMDiffusion2" );
        aParameterLists.set( "constitutive_type",  fem::Constitutive_Type::DIFF_LIN_ISO ) ;
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists.set( "properties", "PropConductivity , Conductivity" );

        //------------------------------------------------------------------------------

        // Nitsche stabilization parameter for structure
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPNitscheStruc" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::DIRICHLET_NITSCHE ) ;
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "leader_properties", "PropYoungs,Material" );

        // Nitsche stabilization parameter for thermal problem
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPNitscheTemp" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::DIRICHLET_NITSCHE ) ;
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "leader_properties", "PropConductivity,Material" );

        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPGPTemp" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::GHOST_DISPL ) ;
        aParameterLists.set( "function_parameters", "0.01" );
        aParameterLists.set( "leader_properties", "PropConductivity,Material" );

        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPInterfaceNitsche" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::NITSCHE_INTERFACE ) ;
        aParameterLists.set( "function_parameters", "10.0" );
        aParameterLists.set( "leader_properties", "PropConductivity,Material" );
        aParameterLists.set( "follower_properties", "PropConductivity,Material" );

        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPInterfaceNitscheDisp" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::NITSCHE_INTERFACE ) ;
        aParameterLists.set( "function_parameters", "10.0" );
        aParameterLists.set( "leader_properties", "PropYoungs,Material" );
        aParameterLists.set( "follower_properties", "PropYoungs,Material" );

        //------------------------------------------------------------------------------
        // create IWG - bulk structure
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGBulkStruct1" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_BULK ) ;
        aParameterLists.set( "dof_residual", "UX,UY,UZ" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY,UZ" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        aParameterLists.set( "leader_properties", "PropBodyLoad,Load" );
        aParameterLists.set( "mesh_set_names", tDomain1 );

        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGBulkStruct2" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_BULK ) ;
        aParameterLists.set( "dof_residual", "UX,UY,UZ" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY,UZ" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso2,ElastLinIso" );
        aParameterLists.set( "leader_properties", "PropBodyLoad,Load" );
        aParameterLists.set( "mesh_set_names", tDomain2 );

        // create IWG - bulk diffusion
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGBulkTemp1" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_BULK ) ;
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusion1,Diffusion" );
        aParameterLists.set( "leader_properties", "PropTempBodyLoad,Load" );
        aParameterLists.set( "mesh_set_names", tDomain1 );

        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGBulkTemp2" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_BULK ) ;
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusion2,Diffusion" );
        aParameterLists.set( "leader_properties", "PropTempBodyLoad,Load" );
        aParameterLists.set( "mesh_set_names", tDomain2 );

        // create IWG - Dirichlet structure
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGDirichlet" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE ) ;
        aParameterLists.set( "dof_residual", "UX,UY,UZ" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY,UZ;TEMP" );
        aParameterLists.set( "leader_properties", "PropDirichlet,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheStruc,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", tLeftSurf );

        // create IWG - Dirichlet temp
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGDirichletTemp" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE ) ;
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY,UZ;TEMP" );
        aParameterLists.set( "leader_properties", "PropTemperature,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusion1,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheTemp,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", tLeftSurf );

        // Interface Dirichlet BC
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGInterfaceTemp" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE ) ;
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY,UZ;TEMP" );
        aParameterLists.set( "follower_dof_dependencies", "UX,UY,UZ;TEMP" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusion1,Diffusion" );
        aParameterLists.set( "follower_constitutive_models", "CMDiffusion2,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPInterfaceNitsche ,NitscheInterface" );
        aParameterLists.set( "mesh_set_names", tInterface );

        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGInterfaceStruc" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_INTERFACE_UNSYMMETRIC_NITSCHE ) ;
        aParameterLists.set( "dof_residual", "UX,UY,UZ" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY,UZ;TEMP" );
        aParameterLists.set( "follower_dof_dependencies", "UX,UY,UZ;TEMP" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        aParameterLists.set( "follower_constitutive_models", "CMStrucLinIso2,ElastLinIso" );
        aParameterLists.set( "stabilization_parameters", "SPInterfaceNitscheDisp ,NitscheInterface" );
        aParameterLists.set( "mesh_set_names", tInterface );

        if ( tUseGhost )
        {
            // Fin Ghost
            // aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
            // aParameterLists.set( "IWG_name",                   "IWGFinGhost" );
            // aParameterLists.set( "IWG_type",                    fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
            // aParameterLists.set( "dof_residual",               "TEMP" );
            // aParameterLists.set( "leader_dof_dependencies",    "TEMP" );
            // aParameterLists.set( "follower_dof_dependencies",     "TEMP" );
            // aParameterLists.set( "stabilization_parameters",   "SPGPTemp,GhostSP" );
            // aParameterLists.set( "mesh_set_names",             tFinGhost );
            // tIWGCounter++;

            // PCM Ghost
            // aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
            // aParameterLists.set( "IWG_name",                   "IWGPcmGhost" );
            // aParameterLists.set( "IWG_type",                    fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
            // aParameterLists.set( "dof_residual",               "TEMP" );
            // aParameterLists.set( "leader_dof_dependencies",    "TEMP" );
            // aParameterLists.set( "follower_dof_dependencies",     "TEMP" );
            // aParameterLists.set( "stabilization_parameters",   "SPGPTemp,GhostSP" );
            // aParameterLists.set( "mesh_set_names",             tPcmGhost );
            // tIWGCounter++;
        }

        // create IWG - Dirichlet structure
        // aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        // aParameterLists.set( "IWG_name",                   "IWGTempNeumann");
        // aParameterLists.set( "IWG_type",                    fem::IWG_Type::SPATIALDIFF_NEUMANN ) ;
        // aParameterLists.set( "dof_residual",               "TEMP");
        // aParameterLists.set( "leader_dof_dependencies",    "UX,UY,UZ;TEMP");
        // aParameterLists.set( "leader_properties",          "PropTempNeumann,Neumann");
        // aParameterLists.set( "mesh_set_names",             tRightSurf);
        // tIWGCounter++;

        //------------------------------------------------------------------------------
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkDISPX" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists.set( "dof_quantity", "UX,UY,UZ" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY,UZ" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "mesh_set_names", tDomain );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkDISPY" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists.set( "dof_quantity", "UX,UY,UZ" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY,UZ" );
        aParameterLists.set( "vectorial_field_index", 1 );
        aParameterLists.set( "mesh_set_names", tDomain );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkDISPZ" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists.set( "dof_quantity", "UX,UY,UZ" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY,UZ" );
        aParameterLists.set( "vectorial_field_index", 2 );
        aParameterLists.set( "mesh_set_names", tDomain );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkTEMP" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists.set( "dof_quantity", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "mesh_set_names", tDomain );

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        aParameterLists( FEM::COMPUTATION );
    }

    /* ------------------------------------------------------------------------ */

    void
    SOLParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::AMESOS_IMPL );
        aParameterLists.set( "Solver_Type", "Amesos_Umfpack" );

        aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list();

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "NLA_combined_res_jac_assembly", true );

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "NLA_DofTypes", "UX,UY,UZ,TEMP" );

        aParameterLists( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "TSA_Num_Time_Steps", 1 );
        aParameterLists.set( "TSA_Time_Frame", 1.0 );

        aParameterLists( SOL::TIME_SOLVERS ).add_parameter_list();
        aParameterLists.set( "TSA_DofTypes", "UX,UY,UZ,TEMP" );
        aParameterLists.set( "TSA_Output_Indices", "0" );
        aParameterLists.set( "TSA_Output_Criteria", "Output_Criterion" );

        aParameterLists( SOL::PRECONDITIONERS ).add_parameter_list(  sol::PreconditionerType::NONE );
    }

    /* ------------------------------------------------------------------------ */

    void
    MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "UX", 1 );
        aParameterLists.set( "UY", 1 );
        aParameterLists.set( "UZ", 1 );
        aParameterLists.set( "TEMP", 0 );
    }

    /* ------------------------------------------------------------------------ */

    void
    VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName ) );
        aParameterLists.set( "Mesh_Type",  vis::VIS_Mesh_Type::STANDARD ) ;
        aParameterLists.set( "Set_Names", tDomain );
        aParameterLists.set( "Field_Names", "UX,UY,UZ,TEMP" );
        aParameterLists.set( "Field_Type", "NODAL,NODAL,NODAL,NODAL" );
        aParameterLists.set( "IQI_Names", "IQIBulkDISPX,IQIBulkDISPY,IQIBulkDISPZ,IQIBulkTEMP" );
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
