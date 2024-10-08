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

        aParameterLists.set( "is_optimization_problem", false);
    }

    void HMRParameterList( Module_Parameter_Lists & aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_hmr_parameter_list() );

        aParameterLists.set( "number_of_elements_per_dimension", "11,11");
        aParameterLists.set( "domain_dimensions",                "4,4");
        aParameterLists.set( "domain_offset",                    "-2.0,-2.0");
        aParameterLists.set( "domain_sidesets",                  "1,2,3,4");
        aParameterLists.set( "lagrange_output_meshes",           "0");

        aParameterLists.set( "lagrange_orders",  "1");
        aParameterLists.set( "lagrange_pattern", "0" );
        aParameterLists.set( "bspline_orders",   "1");
        aParameterLists.set( "bspline_pattern",  "0" );

        aParameterLists.set( "lagrange_to_bspline", "0") ;

        aParameterLists.set( "truncate_bsplines",  1 );
        aParameterLists.set( "refinement_buffer",  1 );
        aParameterLists.set( "staircase_buffer",   1 );
        aParameterLists.set( "initial_refinement", "0" );
        aParameterLists.set( "initial_refinement_pattern", "0" );

        aParameterLists.set( "use_number_aura",    1 );

        aParameterLists.set( "use_multigrid",  0 );
        aParameterLists.set( "severity_level", 0 );
    }

    void XTKParameterList( Module_Parameter_Lists & aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_xtk_parameter_list() );
        aParameterLists.set( "decompose",                     true );
        aParameterLists.set( "decomposition_type",            "conformal") ;
        aParameterLists.set( "enrich",                        true );
        aParameterLists.set( "basis_rank",                    "bspline") ;
        aParameterLists.set( "enrich_mesh_indices",           "0") ;
        aParameterLists.set( "ghost_stab",                    isGhost );
        aParameterLists.set( "multigrid",                     false );
        aParameterLists.set( "verbose",                       true );
        aParameterLists.set( "print_enriched_ig_mesh",        false );
        aParameterLists.set( "exodus_output_XTK_ig_mesh",     true );
        aParameterLists.set( "write_cluster_measures_to_exo", false );
    }

    void GENParameterList( Module_Parameter_Lists & aParameterLists )
    {

        aParameterLists( 0 ).add_parameter_list( prm::create_gen_parameter_list() );
        aParameterLists.set( "output_mesh_file", "gen_output.exo" );

        // Geometry parameter lists
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists.set( "field_function_name", "Colors_1" );
        aParameterLists.set( "number_of_refinements", 2 );
        aParameterLists.set( "refinement_mesh_index", 0 );
        aParameterLists.set( "discretization_mesh_index", -1 );
        aParameterLists.set( "isocontour_tolerance", 1E-12 );

        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists.set( "field_function_name", "Colors_2" );
        aParameterLists.set( "number_of_refinements", 2 );
        aParameterLists.set( "refinement_mesh_index", 0 );
        aParameterLists.set( "discretization_mesh_index", -1 );
        aParameterLists.set( "isocontour_tolerance", 1E-12 );

        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists.set( "field_function_name", "Colors_3" );
        aParameterLists.set( "number_of_refinements", 2 );
        aParameterLists.set( "refinement_mesh_index", 0 );
        aParameterLists.set( "discretization_mesh_index", -1 );
        aParameterLists.set( "isocontour_tolerance", 1E-12 );

        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists.set( "field_function_name", "Colors_4" );
        aParameterLists.set( "number_of_refinements", 2 );
        aParameterLists.set( "refinement_mesh_index", 0 );
        aParameterLists.set( "discretization_mesh_index", -1 );
        aParameterLists.set( "isocontour_tolerance", 1E-12 );

    }

    void FEMParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.hack_for_legacy_fem();
        // create a cell of cell of parameter list for fem

        //------------------------------------------------------------------------------
        // fill the property part of the parameter list

        // create parameter list for property 1
        aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name",            "PropDensity") ;
        aParameterLists.set( "function_parameters",      "1.0") ;
        aParameterLists.set( "value_function",           "Func_Const") ;

        // create parameter list for property 2
        aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name",            "PropCapacity") ;
        aParameterLists.set( "function_parameters",      "1.0") ;
        aParameterLists.set( "value_function",           "Func_Const") ;

        // create parameter list for property 3
        aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name",            "PropConductivity") ;
        aParameterLists.set( "function_parameters",      "0.00005") ;
        aParameterLists.set( "value_function",           "Func_Const") ;

        // create parameter list for property 5
        aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name",            "PropInletTemp") ;
        aParameterLists.set( "function_parameters",      "2.0") ;
        aParameterLists.set( "value_function",           "Func_Const") ;

        // create parameter list for property 6
        aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name",            "PropSideFlux") ;
        aParameterLists.set( "function_parameters",      "2.0") ;
        aParameterLists.set( "value_function",           "Func_Const") ;

        //------------------------------------------------------------------------------
        // fill the constitutive model part of the parameter list

        // create parameter list for constitutive model 2
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists.set( "constitutive_name", "CMDiffusion") ;
        aParameterLists.set( "constitutive_type",  fem::Constitutive_Type::DIFF_LIN_ISO ) ;
        aParameterLists.set( "dof_dependencies",  std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists.set( "properties",
                "PropConductivity,Conductivity;"
                "PropDensity     ,Density;"
                "PropCapacity    ,HeatCapacity") ;

                // create parameter list for constitutive model 2
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists.set( "constitutive_name", "CMDiffusion_2") ;
        aParameterLists.set( "constitutive_type",  fem::Constitutive_Type::DIFF_LIN_ISO ) ;
        aParameterLists.set( "dof_dependencies",  std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists.set( "properties",
                "PropConductivity,Conductivity;"
                "PropDensity     ,Density;"
                "PropCapacity    ,HeatCapacity") ;

        //------------------------------------------------------------------------------
        // fill the stabilization parameter part of the parameter list

        // create parameter list for stabilization parameter 2
        aParameterLists( FEM::STABILIZATION ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists.set( "stabilization_name",      "SPDirichletNitscheT") ;
        aParameterLists.set( "stabilization_type",       fem::Stabilization_Type::DIRICHLET_NITSCHE ) ;
        aParameterLists.set( "function_parameters",     "100.0") ;
        aParameterLists.set( "leader_properties",       "PropConductivity,Material") ;

        // create parameter list for stabilization parameter 8
        aParameterLists( FEM::STABILIZATION ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists.set( "stabilization_name",      "SPGPTemp") ;
        aParameterLists.set( "stabilization_type",       fem::Stabilization_Type::GHOST_DISPL ) ;
        aParameterLists.set( "function_parameters",     "0.005") ;
        aParameterLists.set( "leader_properties",       "PropConductivity,Material") ;

        // Temperature - Shell - PCM
        aParameterLists( FEM::STABILIZATION ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists.set( "stabilization_name",  std::string("SPInterfaceNitsche") );
        aParameterLists.set( "stabilization_type",   fem::Stabilization_Type::NITSCHE_INTERFACE ) ;
        aParameterLists.set( "function_parameters", std::string("100.0") );
        aParameterLists.set( "leader_properties",   std::string("PropConductivity,Material") );
        aParameterLists.set( "follower_properties",    std::string("PropConductivity,Material") );

        //------------------------------------------------------------------------------
        // fill the IWG part of the parameter list

         // create parameter list for IWG 3
        aParameterLists( FEM::IWG ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name",                   "IWGDiffusionBulk") ;
        aParameterLists.set( "IWG_type",                    fem::IWG_Type::SPATIALDIFF_BULK ) ;
        aParameterLists.set( "dof_residual",               "TEMP") ;
        aParameterLists.set( "leader_dof_dependencies",    "TEMP") ;
        aParameterLists.set( "leader_constitutive_models", "CMDiffusion,Diffusion") ;
        aParameterLists.set( "mesh_set_names",             tBulk) ;

        // create parameter list for IWG 11
        aParameterLists( FEM::IWG ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name",                   "IWGInletTemp") ;
        aParameterLists.set( "IWG_type",                    fem::IWG_Type::SPATIALDIFF_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        aParameterLists.set( "dof_residual",               "TEMP") ;
        aParameterLists.set( "leader_dof_dependencies",    "TEMP") ;
        aParameterLists.set( "leader_properties",          "PropInletTemp,Dirichlet") ;
        aParameterLists.set( "leader_constitutive_models", "CMDiffusion,Diffusion") ;
        aParameterLists.set( "stabilization_parameters",   "SPDirichletNitscheT,DirichletNitsche") ;
        aParameterLists.set( "mesh_set_names",             tRightBC) ;

        // create parameter list for IWG 11
        aParameterLists( FEM::IWG ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name",                   "IWGCylinderFluxTemp") ;
        aParameterLists.set( "IWG_type",                    fem::IWG_Type::SPATIALDIFF_NEUMANN ) ;
        aParameterLists.set( "dof_residual",               "TEMP") ;
        aParameterLists.set( "leader_dof_dependencies",    "TEMP") ;
        aParameterLists.set( "leader_properties",          "PropSideFlux,Neumann") ;
        aParameterLists.set( "mesh_set_names",             tLeftBC) ;

        // Temperature - Shell - Shell
        aParameterLists( FEM::IWG ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name",                   std::string("IWGInterfaceShellShellTEMP") );
        aParameterLists.set( "IWG_type",                    fem::IWG_Type::SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE ) ;
        aParameterLists.set( "dof_residual",               std::string("TEMP") );
        aParameterLists.set( "leader_dof_dependencies",    std::string("TEMP") );
        aParameterLists.set( "follower_dof_dependencies",     std::string("TEMP") );
        aParameterLists.set( "leader_constitutive_models", std::string("CMDiffusion,Diffusion") );
        aParameterLists.set( "follower_constitutive_models",  std::string("CMDiffusion_2,Diffusion") );
        aParameterLists.set( "stabilization_parameters",
                std::string("SPInterfaceNitsche     ,NitscheInterface") );
        aParameterLists.set( "mesh_set_names",             tInterface );

        if (isGhost)
        {
            // create parameter list for IWG 16
            aParameterLists( FEM::IWG ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists.set( "IWG_name",                   "IWGGPTemp") ;
            aParameterLists.set( "IWG_type",                    fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
            aParameterLists.set( "dof_residual",               "TEMP") ;
            aParameterLists.set( "leader_dof_dependencies",    "TEMP") ;
            aParameterLists.set( "follower_dof_dependencies",     "TEMP") ;
            aParameterLists.set( "stabilization_parameters",   "SPGPTemp,GhostSP") ;
            aParameterLists.set( "mesh_set_names",             "ghost_p160") ;
            }

        //------------------------------------------------------------------------------
        // fill the IQI part of the parameter list

        // create parameter list for IQI 3
        aParameterLists( FEM::IQI ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists.set( "IQI_name",                   "IQIBulkTEMP") ;
        aParameterLists.set( "IQI_type",                    fem::IQI_Type::DOF ) ;
        aParameterLists.set( "dof_quantity",               "TEMP");
        aParameterLists.set( "leader_dof_dependencies",    "TEMP") ;
        //aParameterLists.set( "vectorial_field_index",      0 );
        aParameterLists.set( "mesh_set_names",             tBulk) ;

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        aParameterLists( FEM::COMPUTATION ).add_parameter_list( prm::create_computation_parameter_list() );
    }

    void SOLParameterList( Module_Parameter_Lists & aParameterLists )
    {

        aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::BELOS_IMPL ) );
        aParameterLists.set( "preconditioners", "0");

        aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list( moris::prm::create_linear_solver_parameter_list() );

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
        aParameterLists.set("NLA_rel_res_norm_drop",    1e-04 );
        aParameterLists.set("NLA_relaxation_parameter", 1.0  );
        aParameterLists.set("NLA_max_iter", 1 );

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
        aParameterLists.set("NLA_DofTypes"      , "TEMP") ;

        aParameterLists( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list( moris::prm::create_time_solver_algorithm_parameter_list() );

        aParameterLists( SOL::TIME_SOLVERS ).add_parameter_list( moris::prm::create_time_solver_parameter_list() );
        aParameterLists.set("TSA_DofTypes"       , "TEMP") ;
        aParameterLists.set("TSA_Initialize_Sol_Vec" , "TEMP,0.0") ;
        aParameterLists.set("TSA_Output_Indices" , "0") ;
        aParameterLists.set("TSA_Output_Criteria" , "Output_Criterion") ;

        aParameterLists( SOL::SOLVER_WAREHOUSE ).add_parameter_list( moris::prm::create_solver_warehouse_parameterlist() );

        aParameterLists( 7 ).add_parameter_list( moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::IFPACK ) );
        aParameterLists.set( "ifpack_prec_type", "ILU" );
    }

    void MSIParameterList( Module_Parameter_Lists & aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_msi_parameter_list() );
        aParameterLists.set("order_adofs_by_host",false);
    }

    void VISParameterList( Module_Parameter_Lists & aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_vis_parameter_list() );
        aParameterLists.set( "File_Name"     , std::pair< std::string, std::string >( "./", "Grain_Structure.exo" ) );
        aParameterLists.set( "Mesh_Type"     ,  vis::VIS_Mesh_Type::STANDARD ) ;
        aParameterLists.set( "Set_Names"     ,  tBulk ) ;
        aParameterLists.set( "Field_Names"   ,  "TEMP" ) ;
        aParameterLists.set( "Field_Type"    ,  "NODAL" ) ;
        aParameterLists.set( "IQI_Names"     ,  "IQIBulkTEMP" ) ;
        aParameterLists.set( "Save_Frequency", 1 );
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

