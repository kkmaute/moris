/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Single_Phase_Hollow_Cylinder_Static.cpp
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

extern uint gInterpolationOrder;

#ifdef __cplusplus
extern "C" {
#endif
//------------------------------------------------------------------------------
namespace moris
{
    // set for FEM
    std::string tSolid      = "HMR_dummy_c_p1,HMR_dummy_n_p1";
    std::string tInSide     = "iside_b0_1_b1_0";
    std::string tOutSide    = "iside_b0_1_b1_3";
    std::string tSolidGhost = "ghost_p1";

    // Geometry parameter
    moris::real tRIn     = 1.0;
    moris::real tROut    = 2.0;
    moris::real tXCenter = 0.0;
    moris::real tYCenter = 0.0;

    // Boundary conditions
    moris::real tImposedTemp      = 5.0;
    moris::real tImposedFlux      = 20.0;
    moris::real tDirichletNitsche = 100.0;
    moris::real tGhostPenalty     = 0.001;

    // Material property
    moris::real tConductivity = 1.0;

    void
    AnalyticTemperatureFunc(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        // get parameters
        real RInner  = tRIn;                                // inner radius
        real ROuter  = tROut;                               // outer radius
        real xCenter = tXCenter;                            // x coord of center
        real yCenter = tYCenter;                            // y coord of center
        real TInner  = tImposedTemp;                        // imposed temperature at inner radius
        real Q       = tImposedFlux * 2 * M_PI * ROuter;    // heat load (W)
        real kappa   = tConductivity;                       // conductivity (W/m^2)

        // get x and y coords
        real tXCoord = aFIManager->get_IP_geometry_interpolator()->valx()( 0 );
        real tYCoord = aFIManager->get_IP_geometry_interpolator()->valx()( 1 );

        // compute radius
        real R = std::sqrt( std::pow( tXCoord - xCenter, 2 ) + std::pow( tYCoord - yCenter, 2 ) );

        aPropMatrix = { { TInner + ( Q * std::log( R / RInner ) ) / ( kappa * 2 * M_PI ) } };
    }

    // Constant function for properties
    void
    Func_Const(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    void
    AnalyticdTemperaturedxFunc(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        // flux
        real tQ = tImposedFlux * 2 * M_PI * tROut;    // heat load (W)

        // radius
        real tX = aFIManager->get_IP_geometry_interpolator()->valx()( 0 );
        real tY = aFIManager->get_IP_geometry_interpolator()->valx()( 1 );
        real tR = std::pow( std::pow( tX - tXCenter, 2.0 ) + std::pow( tY - tYCenter, 2.0 ), 0.5 );

        // set size for aPropMatrix
        aPropMatrix.set_size( 2, 1, 0.0 );
        aPropMatrix( 0, 0 ) = tQ * tX / ( 2.0 * M_PI * tConductivity * std::pow( tR, 2.0 ) );
        aPropMatrix( 1, 0 ) = tQ * tY / ( 2.0 * M_PI * tConductivity * std::pow( tR, 2.0 ) );
    }

    bool
    Output_Criterion( moris::tsa::Time_Solver* aTimeSolver )
    {
        return true;
    }

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

        aParameterLists( 0 ).set( "number_of_elements_per_dimension", "21,21" );
        aParameterLists( 0 ).set( "domain_dimensions", "5,5" );
        aParameterLists( 0 ).set( "domain_offset", "-2.5+0.13,-2.5+0.18" );
        aParameterLists( 0 ).set( "domain_sidesets", "1,2,3,4" );
        aParameterLists( 0 ).set( "lagrange_output_meshes", "0" );

        aParameterLists( 0 ).set( "lagrange_orders", std::to_string( gInterpolationOrder ) );
        aParameterLists( 0 ).set( "lagrange_pattern", "0" );
        aParameterLists( 0 ).set( "bspline_orders", std::to_string( gInterpolationOrder ) );
        aParameterLists( 0 ).set( "bspline_pattern", "0" );

        aParameterLists( 0 ).set( "lagrange_to_bspline", "0" );

        aParameterLists( 0 ).set( "truncate_bsplines", 1 );
        aParameterLists( 0 ).set( "refinement_buffer", 3 );
        aParameterLists( 0 ).set( "staircase_buffer", 3 );
        aParameterLists( 0 ).set( "initial_refinement", "0" );
        aParameterLists( 0 ).set( "initial_refinement_pattern", "0" );

        aParameterLists( 0 ).set( "use_multigrid", 0 );
        aParameterLists( 0 ).set( "severity_level", 0 );

        aParameterLists( 0 ).set( "adaptive_refinement_level", 2 );
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
        aParameterLists( 0 ).set( "ghost_stab", true );
        aParameterLists( 0 ).set( "multigrid", false );
        aParameterLists( 0 ).set( "print_enriched_ig_mesh", true );
        aParameterLists( 0 ).set( "exodus_output_XTK_ig_mesh", true );
    }

    void
    GENParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_gen_parameter_list() );

        // Geometry parameter lists
        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::CIRCLE ) );
        aParameterLists( 1 ).set( "radius", 2.0 );

        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::CIRCLE ) );
        aParameterLists( 1 ).set( "radius", 1.0 );
    }

    void
    FEMParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.hack_for_legacy_fem();
        // create a cell of cell of parameter list for fem

        //------------------------------------------------------------------------------
        // fill the property part of the parameter list

        // create parameter list for property 1
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropConductivity" );
        aParameterLists( 0 ).set( "function_parameters", std::to_string( tConductivity ) );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // create parameter list for property 2
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropImposedTemp" );
        aParameterLists( 0 ).set( "function_parameters", std::to_string( tImposedTemp ) );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // create parameter list for property 3
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropImposedFlux" );
        aParameterLists( 0 ).set( "function_parameters", std::to_string( tImposedFlux ) );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // create parameter list for property 4
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropAnalyticTemp" );
        aParameterLists( 0 ).set( "value_function", "AnalyticTemperatureFunc" );

        // create parameter list for property 5
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropAnalyticdTempdx" );
        aParameterLists( 0 ).set( "value_function", "AnalyticdTemperaturedxFunc" );

        //------------------------------------------------------------------------------
        // fill the constitutive model part of the parameter list

        // create parameter list for constitutive model 1
        aParameterLists( 1 ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists( 1 ).set( "constitutive_name", "CMDiffusion" );
        aParameterLists( 1 ).set( "constitutive_type",  fem::Constitutive_Type::DIFF_LIN_ISO ) ;
        aParameterLists( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists( 1 ).set( "properties", "PropConductivity,Conductivity" );

        //------------------------------------------------------------------------------
        // fill the stabilization parameter part of the parameter list

        // create parameter list for stabilization parameter 1
        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPDirichletNitscheTemp" );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::DIRICHLET_NITSCHE ) ;
        aParameterLists( 2 ).set( "function_parameters", std::to_string( tDirichletNitsche ) );
        aParameterLists( 2 ).set( "leader_properties", "PropConductivity,Material" );

        // create parameter list for stabilization parameter 2
        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPGPTemperature" );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::GHOST_DISPL ) ;
        aParameterLists( 2 ).set( "function_parameters", std::to_string( tGhostPenalty ) );
        aParameterLists( 2 ).set( "leader_properties", "PropConductivity,Material" );

        //------------------------------------------------------------------------------
        // fill the IWG part of the parameter list

        // create parameter list for IWG 1
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGDiffusionBulk" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_BULK ) ;
        aParameterLists( 3 ).set( "dof_residual", "TEMP" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMDiffusion,Diffusion" );
        aParameterLists( 3 ).set( "mesh_set_names", tSolid );

        // create parameter list for IWG 2
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGImposedTempIn" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE ) ;
        aParameterLists( 3 ).set( "dof_residual", "TEMP" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 3 ).set( "leader_properties", "PropImposedTemp,Dirichlet" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMDiffusion,Diffusion" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPDirichletNitscheTemp,DirichletNitsche" );
        aParameterLists( 3 ).set( "mesh_set_names", tInSide );

        // create parameter list for IWG 3
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGImposedFluxOut" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_NEUMANN ) ;
        aParameterLists( 3 ).set( "dof_residual", "TEMP" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 3 ).set( "leader_properties", "PropImposedFlux,Neumann" );
        aParameterLists( 3 ).set( "mesh_set_names", tOutSide );

        // create parameter list for IWG 5
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGGhostTemp" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
        aParameterLists( 3 ).set( "dof_residual", "TEMP" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 3 ).set( "follower_dof_dependencies", "TEMP" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPGPTemperature,GhostSP" );
        aParameterLists( 3 ).set( "mesh_set_names", tSolidGhost );

        //------------------------------------------------------------------------------
        // fill the IQI part of the parameter list

        // create parameter list for IQI 1
        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkTEMP" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists( 4 ).set( "dof_quantity", "TEMP" );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 4 ).set( "vectorial_field_index", 0 );
        aParameterLists( 4 ).set( "mesh_set_names", tSolid );

        // create parameter list for IQI 1
        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkTEMPAnalytic" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::PROPERTY ) ;
        aParameterLists( 4 ).set( "leader_properties", "PropAnalyticTemp,Property" );
        aParameterLists( 4 ).set( "mesh_set_names", tSolid );

        // create parameter list for IQI 2
        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkL2Error" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::L2_ERROR_ANALYTIC ) ;
        aParameterLists( 4 ).set( "dof_quantity", "TEMP" );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 4 ).set( "leader_properties", "PropAnalyticTemp,L2Check" );
        aParameterLists( 4 ).set( "mesh_set_names", tSolid );

        // create parameter list for IQI 3
        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkH1Error" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::H1_ERROR_ANALYTIC ) ;
        aParameterLists( 4 ).set( "dof_quantity", "TEMP" );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 4 ).set( "leader_properties", "PropAnalyticdTempdx,H1Check" );
        aParameterLists( 4 ).set( "mesh_set_names", tSolid );

        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIMaxTEMP" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::MAX_DOF ) ;
        aParameterLists( 4 ).set( "dof_quantity", "TEMP" );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 4 ).set( "function_parameters", "5.0/30.0" );
        aParameterLists( 4 ).set( "mesh_set_names", tSolid );

        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIVolume" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::VOLUME ) ;
        aParameterLists( 4 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 4 ).set( "mesh_set_names", tSolid );

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
        aParameterLists( 2 ).set( "NLA_rel_res_norm_drop", 1e-07 );

        aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
        aParameterLists( 3 ).set( "NLA_DofTypes", "TEMP" );

        aParameterLists( 4 ).add_parameter_list( moris::prm::create_time_solver_algorithm_parameter_list() );

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
    }

    void
    VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_vis_parameter_list() );
        aParameterLists( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", "Single_Phase_Hollow_Cylinder_Static.exo" ) );
        aParameterLists( 0 ).set( "Mesh_Type",  vis::VIS_Mesh_Type::STANDARD ) ;
        aParameterLists( 0 ).set( "Set_Names", tSolid );
        aParameterLists( 0 ).set( "Field_Names", "TEMP,TEMP_ANALYTIC,L2_ERROR_ANALYTIC,H1_ERROR_ANALYTIC,MAX_DOF,VOLUME" );
        aParameterLists( 0 ).set( "Field_Type", "NODAL,NODAL,GLOBAL,GLOBAL,GLOBAL,GLOBAL" );
        aParameterLists( 0 ).set( "IQI_Names", "IQIBulkTEMP,IQIBulkTEMPAnalytic,IQIBulkL2Error,IQIBulkH1Error,IQIMaxTEMP,IQIVolume" );
        aParameterLists( 0 ).set( "Save_Frequency", 1 );
    }

    void
    MORISGENERALParameterList( Module_Parameter_Lists& aParameterLists )
    {
    }

    //------------------------------------------------------------------------------
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
