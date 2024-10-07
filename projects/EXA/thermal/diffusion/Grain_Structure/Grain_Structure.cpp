/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Grain_Structure.cpp
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

#ifdef  __cplusplus
extern "C"
{
#endif
//------------------------------------------------------------------------------
namespace moris
{
    bool isGhost = false;

    std::string tBulk = "HMR_dummy_n_p0,HMR_dummy_n_p1,HMR_dummy_n_p2,HMR_dummy_n_p4,HMR_dummy_n_p8,HMR_dummy_c_p0,HMR_dummy_c_p1,HMR_dummy_c_p2,HMR_dummy_c_p4,HMR_dummy_c_p8";
    std::string tRightBC = "SideSet_2_c_p0,SideSet_2_c_p4,SideSet_2_c_p8,SideSet_2_n_p0,SideSet_2_n_p4,SideSet_2_n_p8";
    std::string tLeftBC = "SideSet_4_c_p1,SideSet_4_c_p2,SideSet_4_n_p1,SideSet_4_n_p2";
    std::string tInterface = "dbl_iside_p0_0_p1_1,dbl_iside_p0_0_p1_2,dbl_iside_p0_0_p1_8,dbl_iside_p0_0_p1_4,dbl_iside_p0_1_p1_2,dbl_iside_p0_1_p1_8";

    // Constant function for properties
    void Func_Const(
            moris::Matrix< moris::DDRMat >                 & aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >  & aParameters,
            moris::fem::Field_Interpolator_Manager         * aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    bool Output_Criterion( moris::tsa::Time_Solver * aTimeSolver )
    {
        return true;
    }

    moris::real Colors_1(
                const moris::Matrix< DDRMat >     & aCoordinates,
                const Vector< real > & aGeometryParameters)
    {
        moris::real tVal = -1.0;
        if( aCoordinates(0) > 0.1117 && aCoordinates(1) <= -0.11)
        {
            tVal =1.0;
        }

        return tVal;
    }

    moris::real Colors_2(
                const moris::Matrix< DDRMat >     & aCoordinates,
                const Vector< real > & aGeometryParameters)
    {
        moris::real tVal = -1.0;
        if( std::sqrt( std::pow( aCoordinates(0) -2.0, 2) + std::pow( aCoordinates(1) -2.0, 2) ) - 1.2 <= 0.0 &&
            !(aCoordinates(0) > 0.1117 && aCoordinates(1) <= -0.11 ) )
        {
            tVal =1.0;
        }

        return tVal;
    }

       moris::real Colors_3(
                const moris::Matrix< DDRMat >     & aCoordinates,
                const Vector< real > & aGeometryParameters)
    {
        moris::real tVal = -1.0;
        if(  aCoordinates(0)-aCoordinates(1) +  1.5 <= 0.0 &&
            !( std::sqrt( std::pow( aCoordinates(0) -2.0, 2) + std::pow( aCoordinates(1) -2.0, 2) ) - 1.2 <= 0.0) &&
            !(aCoordinates(0) > 0.1117 && aCoordinates(1) <= -0.11 ))
        {
            tVal =1.0;
        }

        return tVal;
    }

       moris::real Colors_4(
                const moris::Matrix< DDRMat >     & aCoordinates,
                const Vector< real > & aGeometryParameters)
    {
        moris::real tVal = -1.0;
        if( std::sqrt( std::pow( aCoordinates(0) +2.0, 2) + std::pow( aCoordinates(1) +2.0, 2) ) - 2.6 <= 0.0 &&
            !( aCoordinates(0)-aCoordinates(1) +  1.5 <= 0.0 ) &&
            !( std::sqrt( std::pow( aCoordinates(0) -2.0, 2) + std::pow( aCoordinates(1) -2.0, 2) ) - 1.2 <= 0.0) &&
            !( aCoordinates(0) > 0.1117 && aCoordinates(1) <= -0.11 ))
        {
            tVal =1.0;
        }

        return tVal;
    }

    void OPTParameterList( Module_Parameter_Lists & aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_opt_problem_parameter_list() );

        aParameterLists( 0 ).set( "is_optimization_problem", false);
    }

    void HMRParameterList( Module_Parameter_Lists & aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_hmr_parameter_list() );

        aParameterLists( 0 ).set( "number_of_elements_per_dimension", "11,11");
        aParameterLists( 0 ).set( "domain_dimensions",                "4,4");
        aParameterLists( 0 ).set( "domain_offset",                    "-2.0,-2.0");
        aParameterLists( 0 ).set( "domain_sidesets",                  "1,2,3,4");
        aParameterLists( 0 ).set( "lagrange_output_meshes",           "0");

        aParameterLists( 0 ).set( "lagrange_orders",  "1");
        aParameterLists( 0 ).set( "lagrange_pattern", "0" );
        aParameterLists( 0 ).set( "bspline_orders",   "1");
        aParameterLists( 0 ).set( "bspline_pattern",  "0" );

        aParameterLists( 0 ).set( "lagrange_to_bspline", "0") ;

        aParameterLists( 0 ).set( "truncate_bsplines",  1 );
        aParameterLists( 0 ).set( "refinement_buffer",  1 );
        aParameterLists( 0 ).set( "staircase_buffer",   1 );
        aParameterLists( 0 ).set( "initial_refinement", "0" );
        aParameterLists( 0 ).set( "initial_refinement_pattern", "0" );

        aParameterLists( 0 ).set( "use_number_aura",    1 );

        aParameterLists( 0 ).set( "use_multigrid",  0 );
        aParameterLists( 0 ).set( "severity_level", 0 );
    }

    void XTKParameterList( Module_Parameter_Lists & aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_xtk_parameter_list() );
        aParameterLists( 0 ).set( "decompose",                     true );
        aParameterLists( 0 ).set( "decomposition_type",            "conformal") ;
        aParameterLists( 0 ).set( "enrich",                        true );
        aParameterLists( 0 ).set( "basis_rank",                    "bspline") ;
        aParameterLists( 0 ).set( "enrich_mesh_indices",           "0") ;
        aParameterLists( 0 ).set( "ghost_stab",                    isGhost );
        aParameterLists( 0 ).set( "multigrid",                     false );
        aParameterLists( 0 ).set( "verbose",                       true );
        aParameterLists( 0 ).set( "print_enriched_ig_mesh",        false );
        aParameterLists( 0 ).set( "exodus_output_XTK_ig_mesh",     true );
        aParameterLists( 0 ).set( "write_cluster_measures_to_exo", false );
    }

    void GENParameterList( Module_Parameter_Lists & aParameterLists )
    {

        aParameterLists( 0 ).add_parameter_list( prm::create_gen_parameter_list() );
        aParameterLists( 0 ).set( "output_mesh_file", "gen_output.exo" );

        // Geometry parameter lists
        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists( 1 ).set( "field_function_name", "Colors_1" );
        aParameterLists( 1 ).set( "number_of_refinements", 2 );
        aParameterLists( 1 ).set( "refinement_mesh_index", 0 );
        aParameterLists( 1 ).set( "discretization_mesh_index", -1 );
        aParameterLists( 1 ).set( "isocontour_tolerance", 1E-12 );

        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists( 1 ).set( "field_function_name", "Colors_2" );
        aParameterLists( 1 ).set( "number_of_refinements", 2 );
        aParameterLists( 1 ).set( "refinement_mesh_index", 0 );
        aParameterLists( 1 ).set( "discretization_mesh_index", -1 );
        aParameterLists( 1 ).set( "isocontour_tolerance", 1E-12 );

        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists( 1 ).set( "field_function_name", "Colors_3" );
        aParameterLists( 1 ).set( "number_of_refinements", 2 );
        aParameterLists( 1 ).set( "refinement_mesh_index", 0 );
        aParameterLists( 1 ).set( "discretization_mesh_index", -1 );
        aParameterLists( 1 ).set( "isocontour_tolerance", 1E-12 );

        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists( 1 ).set( "field_function_name", "Colors_4" );
        aParameterLists( 1 ).set( "number_of_refinements", 2 );
        aParameterLists( 1 ).set( "refinement_mesh_index", 0 );
        aParameterLists( 1 ).set( "discretization_mesh_index", -1 );
        aParameterLists( 1 ).set( "isocontour_tolerance", 1E-12 );

    }

    void FEMParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.hack_for_legacy_fem();
        // create a cell of cell of parameter list for fem

        //------------------------------------------------------------------------------
        // fill the property part of the parameter list

        // create parameter list for property 1
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name",            "PropDensity") ;
        aParameterLists( 0 ).set( "function_parameters",      "1.0") ;
        aParameterLists( 0 ).set( "value_function",           "Func_Const") ;

        // create parameter list for property 2
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name",            "PropCapacity") ;
        aParameterLists( 0 ).set( "function_parameters",      "1.0") ;
        aParameterLists( 0 ).set( "value_function",           "Func_Const") ;

        // create parameter list for property 3
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name",            "PropConductivity") ;
        aParameterLists( 0 ).set( "function_parameters",      "0.00005") ;
        aParameterLists( 0 ).set( "value_function",           "Func_Const") ;

        // create parameter list for property 5
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name",            "PropInletTemp") ;
        aParameterLists( 0 ).set( "function_parameters",      "2.0") ;
        aParameterLists( 0 ).set( "value_function",           "Func_Const") ;

        // create parameter list for property 6
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name",            "PropSideFlux") ;
        aParameterLists( 0 ).set( "function_parameters",      "2.0") ;
        aParameterLists( 0 ).set( "value_function",           "Func_Const") ;

        //------------------------------------------------------------------------------
        // fill the constitutive model part of the parameter list

        // create parameter list for constitutive model 2
        aParameterLists( 1 ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists( 1 ).set( "constitutive_name", "CMDiffusion") ;
        aParameterLists( 1 ).set( "constitutive_type",  fem::Constitutive_Type::DIFF_LIN_ISO ) ;
        aParameterLists( 1 ).set( "dof_dependencies",  std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists( 1 ).set( "properties",
                "PropConductivity,Conductivity;"
                "PropDensity     ,Density;"
                "PropCapacity    ,HeatCapacity") ;

                // create parameter list for constitutive model 2
        aParameterLists( 1 ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists( 1 ).set( "constitutive_name", "CMDiffusion_2") ;
        aParameterLists( 1 ).set( "constitutive_type",  fem::Constitutive_Type::DIFF_LIN_ISO ) ;
        aParameterLists( 1 ).set( "dof_dependencies",  std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists( 1 ).set( "properties",
                "PropConductivity,Conductivity;"
                "PropDensity     ,Density;"
                "PropCapacity    ,HeatCapacity") ;

        //------------------------------------------------------------------------------
        // fill the stabilization parameter part of the parameter list

        // create parameter list for stabilization parameter 2
        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name",      "SPDirichletNitscheT") ;
        aParameterLists( 2 ).set( "stabilization_type",       fem::Stabilization_Type::DIRICHLET_NITSCHE ) ;
        aParameterLists( 2 ).set( "function_parameters",     "100.0") ;
        aParameterLists( 2 ).set( "leader_properties",       "PropConductivity,Material") ;

        // create parameter list for stabilization parameter 8
        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name",      "SPGPTemp") ;
        aParameterLists( 2 ).set( "stabilization_type",       fem::Stabilization_Type::GHOST_DISPL ) ;
        aParameterLists( 2 ).set( "function_parameters",     "0.005") ;
        aParameterLists( 2 ).set( "leader_properties",       "PropConductivity,Material") ;

        // Temperature - Shell - PCM
        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name",  std::string("SPInterfaceNitsche") );
        aParameterLists( 2 ).set( "stabilization_type",   fem::Stabilization_Type::NITSCHE_INTERFACE ) ;
        aParameterLists( 2 ).set( "function_parameters", std::string("100.0") );
        aParameterLists( 2 ).set( "leader_properties",   std::string("PropConductivity,Material") );
        aParameterLists( 2 ).set( "follower_properties",    std::string("PropConductivity,Material") );

        //------------------------------------------------------------------------------
        // fill the IWG part of the parameter list

         // create parameter list for IWG 3
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name",                   "IWGDiffusionBulk") ;
        aParameterLists( 3 ).set( "IWG_type",                    fem::IWG_Type::SPATIALDIFF_BULK ) ;
        aParameterLists( 3 ).set( "dof_residual",               "TEMP") ;
        aParameterLists( 3 ).set( "leader_dof_dependencies",    "TEMP") ;
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMDiffusion,Diffusion") ;
        aParameterLists( 3 ).set( "mesh_set_names",             tBulk) ;

        // create parameter list for IWG 11
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name",                   "IWGInletTemp") ;
        aParameterLists( 3 ).set( "IWG_type",                    fem::IWG_Type::SPATIALDIFF_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        aParameterLists( 3 ).set( "dof_residual",               "TEMP") ;
        aParameterLists( 3 ).set( "leader_dof_dependencies",    "TEMP") ;
        aParameterLists( 3 ).set( "leader_properties",          "PropInletTemp,Dirichlet") ;
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMDiffusion,Diffusion") ;
        aParameterLists( 3 ).set( "stabilization_parameters",   "SPDirichletNitscheT,DirichletNitsche") ;
        aParameterLists( 3 ).set( "mesh_set_names",             tRightBC) ;

        // create parameter list for IWG 11
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name",                   "IWGCylinderFluxTemp") ;
        aParameterLists( 3 ).set( "IWG_type",                    fem::IWG_Type::SPATIALDIFF_NEUMANN ) ;
        aParameterLists( 3 ).set( "dof_residual",               "TEMP") ;
        aParameterLists( 3 ).set( "leader_dof_dependencies",    "TEMP") ;
        aParameterLists( 3 ).set( "leader_properties",          "PropSideFlux,Neumann") ;
        aParameterLists( 3 ).set( "mesh_set_names",             tLeftBC) ;

        // Temperature - Shell - Shell
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name",                   std::string("IWGInterfaceShellShellTEMP") );
        aParameterLists( 3 ).set( "IWG_type",                    fem::IWG_Type::SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE ) ;
        aParameterLists( 3 ).set( "dof_residual",               std::string("TEMP") );
        aParameterLists( 3 ).set( "leader_dof_dependencies",    std::string("TEMP") );
        aParameterLists( 3 ).set( "follower_dof_dependencies",     std::string("TEMP") );
        aParameterLists( 3 ).set( "leader_constitutive_models", std::string("CMDiffusion,Diffusion") );
        aParameterLists( 3 ).set( "follower_constitutive_models",  std::string("CMDiffusion_2,Diffusion") );
        aParameterLists( 3 ).set( "stabilization_parameters",
                std::string("SPInterfaceNitsche     ,NitscheInterface") );
        aParameterLists( 3 ).set( "mesh_set_names",             tInterface );

        if (isGhost)
        {
            // create parameter list for IWG 16
            aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists( 3 ).set( "IWG_name",                   "IWGGPTemp") ;
            aParameterLists( 3 ).set( "IWG_type",                    fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
            aParameterLists( 3 ).set( "dof_residual",               "TEMP") ;
            aParameterLists( 3 ).set( "leader_dof_dependencies",    "TEMP") ;
            aParameterLists( 3 ).set( "follower_dof_dependencies",     "TEMP") ;
            aParameterLists( 3 ).set( "stabilization_parameters",   "SPGPTemp,GhostSP") ;
            aParameterLists( 3 ).set( "mesh_set_names",             "ghost_p160") ;
            }

        //------------------------------------------------------------------------------
        // fill the IQI part of the parameter list

        // create parameter list for IQI 3
        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name",                   "IQIBulkTEMP") ;
        aParameterLists( 4 ).set( "IQI_type",                    fem::IQI_Type::DOF ) ;
        aParameterLists( 4 ).set( "dof_quantity",               "TEMP");
        aParameterLists( 4 ).set( "leader_dof_dependencies",    "TEMP") ;
        //aParameterLists( 4 ).set( "vectorial_field_index",      0 );
        aParameterLists( 4 ).set( "mesh_set_names",             tBulk) ;

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        aParameterLists( 5 ).add_parameter_list( prm::create_computation_parameter_list() );
    }

    void SOLParameterList( Module_Parameter_Lists & aParameterLists )
    {

        aParameterLists( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::BELOS_IMPL ) );
        aParameterLists( 0 ).set( "preconditioners", "0");

        aParameterLists( 1 ).add_parameter_list( moris::prm::create_linear_solver_parameter_list() );

        aParameterLists( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
        aParameterLists( 2 ).set("NLA_rel_res_norm_drop",    1e-04 );
        aParameterLists( 2 ).set("NLA_relaxation_parameter", 1.0  );
        aParameterLists( 2 ).set("NLA_max_iter", 1 );

        aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
        aParameterLists( 3 ).set("NLA_DofTypes"      , "TEMP") ;

        aParameterLists( 4 ).add_parameter_list( moris::prm::create_time_solver_algorithm_parameter_list() );

        aParameterLists( 5 ).add_parameter_list( moris::prm::create_time_solver_parameter_list() );
        aParameterLists( 5 ).set("TSA_DofTypes"       , "TEMP") ;
        aParameterLists( 5 ).set("TSA_Initialize_Sol_Vec" , "TEMP,0.0") ;
        aParameterLists( 5 ).set("TSA_Output_Indices" , "0") ;
        aParameterLists( 5 ).set("TSA_Output_Criteria" , "Output_Criterion") ;

        aParameterLists( 6 ).add_parameter_list( moris::prm::create_solver_warehouse_parameterlist() );

        aParameterLists( 7 ).add_parameter_list( moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::IFPACK ) );
        aParameterLists( 7 ).set( "ifpack_prec_type", "ILU" );
    }

    void MSIParameterList( Module_Parameter_Lists & aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_msi_parameter_list() );
        aParameterLists( 0 ).set("order_adofs_by_host",false);
    }

    void VISParameterList( Module_Parameter_Lists & aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_vis_parameter_list() );
        aParameterLists( 0 ).set( "File_Name"     , std::pair< std::string, std::string >( "./", "Grain_Structure.exo" ) );
        aParameterLists( 0 ).set( "Mesh_Type"     ,  vis::VIS_Mesh_Type::STANDARD ) ;
        aParameterLists( 0 ).set( "Set_Names"     ,  tBulk ) ;
        aParameterLists( 0 ).set( "Field_Names"   ,  "TEMP" ) ;
        aParameterLists( 0 ).set( "Field_Type"    ,  "NODAL" ) ;
        aParameterLists( 0 ).set( "IQI_Names"     ,  "IQIBulkTEMP" ) ;
        aParameterLists( 0 ).set( "Save_Frequency", 1 );
    }

    void MORISGENERALParameterList( Module_Parameter_Lists & aParameterLists )
    {

    }

    //------------------------------------------------------------------------------
}

//------------------------------------------------------------------------------
#ifdef  __cplusplus
}
#endif

