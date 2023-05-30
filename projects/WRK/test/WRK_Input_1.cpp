/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * WRK_Input_1.cpp
 *
 */

#include <string>
#include <iostream>
#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_MSI_Equation_Object.hpp"
#include "cl_TSA_Time_Solver.hpp"
#include "cl_DLA_Solver_Interface.hpp"
#include "fn_equal_to.hpp"
#include "fn_PRM_FEM_Parameters.hpp"
#include "fn_PRM_MSI_Parameters.hpp"
#include "fn_PRM_SOL_Parameters.hpp"
#include "fn_PRM_VIS_Parameters.hpp"
#include "fn_PRM_HMR_Parameters.hpp"
#include "fn_PRM_GEN_Parameters.hpp"
#include "fn_PRM_XTK_Parameters.hpp"
#include "fn_PRM_OPT_Parameters.hpp"

#include "cl_DLA_Linear_Solver_Aztec.hpp"
#include "AztecOO.h"

#ifdef __cplusplus
extern "C" {
#endif
//------------------------------------------------------------------------------
namespace moris
{

    void
    Func1( moris::Matrix< moris::DDRMat >&                 aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    bool
    Output_Criterion( moris::tsa::Time_Solver* aTimeSolver )
    {
        return true;
    }

    moris::real
    Lvl_set_1(
            const moris::Matrix< DDRMat >&     aCoordinates,
            const moris::Cell< moris::real* >& aGeometryParameters )
    {
        return 1.01;
    }

    void
    FEMParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterList )
    {
        // create a cell of cell of parameter list for fem
        tParameterList.resize( 8 );

        //------------------------------------------------------------------------------
        // fill the property part of the parameter list
        moris::uint tNumProperties = 9;
        tParameterList( 0 ).resize( tNumProperties );

        // create parameter list for property 1
        tParameterList( 0 )( 0 ) = prm::create_property_parameter_list();
        tParameterList( 0 )( 0 ).set( "property_name", "PropConductivity1" );
        tParameterList( 0 )( 0 ).set( "function_parameters", "1.0" );
        tParameterList( 0 )( 0 ).set( "value_function", "Func1" );

        // create parameter list for property 2
        tParameterList( 0 )( 1 ) = prm::create_property_parameter_list();
        tParameterList( 0 )( 1 ).set( "property_name", "PropConductivity2" );
        tParameterList( 0 )( 1 ).set( "function_parameters", "5.0" );
        tParameterList( 0 )( 1 ).set( "value_function", "Func1" );

        // create parameter list for property 3
        tParameterList( 0 )( 2 ) = prm::create_property_parameter_list();
        tParameterList( 0 )( 2 ).set( "property_name", "PropDirichletU" );
        tParameterList( 0 )( 2 ).set( "function_parameters", "0.0;0.0" );
        tParameterList( 0 )( 2 ).set( "value_function", "Func1" );

        // create parameter list for property 4
        tParameterList( 0 )( 3 ) = prm::create_property_parameter_list();
        tParameterList( 0 )( 3 ).set( "property_name", "PropDirichletTEMP" );
        tParameterList( 0 )( 3 ).set( "function_parameters", "3.0" );
        tParameterList( 0 )( 3 ).set( "value_function", "Func1" );

        // create parameter list for property 5
        tParameterList( 0 )( 4 ) = prm::create_property_parameter_list();
        tParameterList( 0 )( 4 ).set( "property_name", "PropEMod1" );
        tParameterList( 0 )( 4 ).set( "function_parameters", "1.0" );
        tParameterList( 0 )( 4 ).set( "value_function", "Func1" );

        // create parameter list for property 6
        tParameterList( 0 )( 5 ) = prm::create_property_parameter_list();
        tParameterList( 0 )( 5 ).set( "property_name", "PropEMod2" );
        tParameterList( 0 )( 5 ).set( "function_parameters", "1.0" );
        tParameterList( 0 )( 5 ).set( "value_function", "Func1" );

        // create parameter list for property 7
        tParameterList( 0 )( 6 ) = prm::create_property_parameter_list();
        tParameterList( 0 )( 6 ).set( "property_name", "PropPoisson" );
        tParameterList( 0 )( 6 ).set( "function_parameters", "0.0" );
        tParameterList( 0 )( 6 ).set( "value_function", "Func1" );

        // create parameter list for property 8
        tParameterList( 0 )( 7 ) = prm::create_property_parameter_list();
        tParameterList( 0 )( 7 ).set( "property_name", "PropCTE" );
        tParameterList( 0 )( 7 ).set( "function_parameters", "1.0" );
        tParameterList( 0 )( 7 ).set( "value_function", "Func1" );

        // create parameter list for property 9
        tParameterList( 0 )( 8 ) = prm::create_property_parameter_list();
        tParameterList( 0 )( 8 ).set( "property_name", "PropTRef" );
        tParameterList( 0 )( 8 ).set( "function_parameters", "1.0" );
        tParameterList( 0 )( 8 ).set( "value_function", "Func1" );

        //------------------------------------------------------------------------------
        // fill the constitutive model part of the parameter list
        moris::uint tNumCMs = 4;
        tParameterList( 1 ).resize( tNumCMs );

        // create parameter list for constitutive model 1
        tParameterList( 1 )( 0 ) = prm::create_constitutive_model_parameter_list();
        tParameterList( 1 )( 0 ).set( "constitutive_name", "CMStrucLinIso1" );
        tParameterList( 1 )( 0 ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::STRUC_LIN_ISO ) );
        tParameterList( 1 )( 0 ).set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY;TEMP", "Displacement,Temperature" ) );
        tParameterList( 1 )( 0 ).set( "properties", "PropEMod1,YoungsModulus;PropPoisson,PoissonRatio;PropCTE,CTE;PropTRef,ReferenceTemperature" );
        tParameterList( 1 )( 0 ).set( "model_type", static_cast< uint >( fem::Model_Type::PLANE_STRESS ) );

        // create parameter list for constitutive model 2
        tParameterList( 1 )( 1 ) = prm::create_constitutive_model_parameter_list();
        tParameterList( 1 )( 1 ).set( "constitutive_name", "CMStrucLinIso2" );
        tParameterList( 1 )( 1 ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::STRUC_LIN_ISO ) );
        tParameterList( 1 )( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY;TEMP", "Displacement,Temperature" ) );
        tParameterList( 1 )( 1 ).set( "properties", "PropEMod2,YoungsModulus;PropPoisson,PoissonRatio;PropCTE,CTE;PropTRef,ReferenceTemperature" );
        tParameterList( 1 )( 1 ).set( "model_type", static_cast< uint >( fem::Model_Type::PLANE_STRESS ) );

        // create parameter list for constitutive model 3
        tParameterList( 1 )( 2 ) = prm::create_constitutive_model_parameter_list();
        tParameterList( 1 )( 2 ).set( "constitutive_name", "CMDiffLinIso1" );
        tParameterList( 1 )( 2 ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::DIFF_LIN_ISO ) );
        tParameterList( 1 )( 2 ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( 1 )( 2 ).set( "properties", "PropConductivity1,Conductivity" );

        // create parameter list for constitutive model 4
        tParameterList( 1 )( 3 ) = prm::create_constitutive_model_parameter_list();
        tParameterList( 1 )( 3 ).set( "constitutive_name", "CMDiffLinIso2" );
        tParameterList( 1 )( 3 ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::DIFF_LIN_ISO ) );
        tParameterList( 1 )( 3 ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( 1 )( 3 ).set( "properties", "PropConductivity2,Conductivity" );

        //------------------------------------------------------------------------------
        // fill the stabilization parameter part of the parameter list
        moris::uint tNumSPs = 2;
        tParameterList( 2 ).resize( tNumSPs );

        // create parameter list for stabilization parameter 1
        tParameterList( 2 )( 0 ) = prm::create_stabilization_parameter_parameter_list();
        tParameterList( 2 )( 0 ).set( "stabilization_name", "SPDirichletNitscheU" );
        tParameterList( 2 )( 0 ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::DIRICHLET_NITSCHE ) );
        tParameterList( 2 )( 0 ).set( "function_parameters", "100.0" );
        tParameterList( 2 )( 0 ).set( "leader_properties", "PropEMod1,Material" );

        // create parameter list for stabilization parameter 2
        tParameterList( 2 )( 1 ) = prm::create_stabilization_parameter_parameter_list();
        tParameterList( 2 )( 1 ).set( "stabilization_name", "SPDirichletNitscheTEMP" );
        tParameterList( 2 )( 1 ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::DIRICHLET_NITSCHE ) );
        tParameterList( 2 )( 1 ).set( "function_parameters", "1.0" );
        tParameterList( 2 )( 1 ).set( "leader_properties", "PropConductivity1,Material" );

        //------------------------------------------------------------------------------
        // fill the IWG part of the parameter list
        moris::uint tNumIWGs = 6;
        tParameterList( 3 ).resize( tNumIWGs );

        // create parameter list for IWG 1
        tParameterList( 3 )( 0 ) = prm::create_IWG_parameter_list();
        tParameterList( 3 )( 0 ).set( "IWG_name", "IWGBulkU_1" );
        tParameterList( 3 )( 0 ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_BULK ) );
        tParameterList( 3 )( 0 ).set( "dof_residual", "UX,UY" );
        tParameterList( 3 )( 0 ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 3 )( 0 ).set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        tParameterList( 3 )( 0 ).set( "mesh_set_names", "HMR_dummy_n_p1" );

        // create parameter list for IWG 2
        tParameterList( 3 )( 1 ) = prm::create_IWG_parameter_list();
        tParameterList( 3 )( 1 ).set( "IWG_name", "IWGBulkU_2" );
        tParameterList( 3 )( 1 ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_BULK ) );
        tParameterList( 3 )( 1 ).set( "dof_residual", "UX,UY" );
        tParameterList( 3 )( 1 ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 3 )( 1 ).set( "leader_constitutive_models", "CMStrucLinIso2,ElastLinIso" );
        tParameterList( 3 )( 1 ).set( "mesh_set_names", "" );

        // create parameter list for IWG 3
        tParameterList( 3 )( 2 ) = prm::create_IWG_parameter_list();
        tParameterList( 3 )( 2 ).set( "IWG_name", "IWGDirichletU" );
        tParameterList( 3 )( 2 ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE ) );
        tParameterList( 3 )( 2 ).set( "dof_residual", "UX,UY" );
        tParameterList( 3 )( 2 ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 3 )( 2 ).set( "leader_properties", "PropDirichletU,Dirichlet" );
        tParameterList( 3 )( 2 ).set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        tParameterList( 3 )( 2 ).set( "stabilization_parameters", "SPDirichletNitscheU,DirichletNitsche" );
        tParameterList( 3 )( 2 ).set( "mesh_set_names", "SideSet_4_n_p1" );

        // create parameter list for IWG 4
        tParameterList( 3 )( 3 ) = prm::create_IWG_parameter_list();
        tParameterList( 3 )( 3 ).set( "IWG_name", "IWGBulkTEMP_1" );
        tParameterList( 3 )( 3 ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::SPATIALDIFF_BULK ) );
        tParameterList( 3 )( 3 ).set( "dof_residual", "TEMP" );
        tParameterList( 3 )( 3 ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( 3 )( 3 ).set( "leader_constitutive_models", "CMDiffLinIso1,Diffusion" );
        tParameterList( 3 )( 3 ).set( "mesh_set_names", "HMR_dummy_n_p1" );

        // create parameter list for IWG 5
        tParameterList( 3 )( 4 ) = prm::create_IWG_parameter_list();
        tParameterList( 3 )( 4 ).set( "IWG_name", "IWGBulkTEMP_2" );
        tParameterList( 3 )( 4 ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::SPATIALDIFF_BULK ) );
        tParameterList( 3 )( 4 ).set( "dof_residual", "TEMP" );
        tParameterList( 3 )( 4 ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( 3 )( 4 ).set( "leader_constitutive_models", "CMDiffLinIso2,Diffusion" );
        tParameterList( 3 )( 4 ).set( "mesh_set_names", "" );

        // create parameter list for IWG 6
        tParameterList( 3 )( 5 ) = prm::create_IWG_parameter_list();
        tParameterList( 3 )( 5 ).set( "IWG_name", "IWGDirichletTEMP" );
        tParameterList( 3 )( 5 ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE ) );
        tParameterList( 3 )( 5 ).set( "dof_residual", "TEMP" );
        tParameterList( 3 )( 5 ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( 3 )( 5 ).set( "leader_properties", "PropDirichletTEMP,Dirichlet" );
        tParameterList( 3 )( 5 ).set( "leader_constitutive_models", "CMDiffLinIso2,Diffusion" );
        tParameterList( 3 )( 5 ).set( "stabilization_parameters", "SPDirichletNitscheTEMP,DirichletNitsche" );
        tParameterList( 3 )( 5 ).set( "mesh_set_names", "SideSet_4_n_p1" );

        //------------------------------------------------------------------------------
        // fill the IQI part of the parameter list
        moris::uint tNumIQIs = 1;
        tParameterList( 4 ).resize( tNumIQIs );

        // create parameter list for IQI 1
        tParameterList( 4 )( 0 ) = prm::create_IQI_parameter_list();
        tParameterList( 4 )( 0 ).set( "IQI_name", "IQIBulkU_1" );
        tParameterList( 4 )( 0 ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::STRAIN_ENERGY ) );
        tParameterList( 4 )( 0 ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 4 )( 0 ).set( "leader_constitutive_models", "CMStrucLinIso1,Elast" );
        tParameterList( 4 )( 0 ).set( "mesh_set_names", "HMR_dummy_n_p1" );

        // create computation parameter list
        tParameterList( 5 ).resize( 1 );
        tParameterList( 5 )( 0 ) = prm::create_computation_parameter_list();
    }

    void
    SOLParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 7 );
        for ( uint Ik = 0; Ik < 7; Ik++ )
        {
            tParameterlist( Ik ).resize( 1 );
        }

        tParameterlist( 0 )( 0 ) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AZTEC_IMPL );
        tParameterlist( 0 )( 0 ).set( "AZ_diagnostics", AZ_none );
        tParameterlist( 0 )( 0 ).set( "AZ_output", AZ_none );
        tParameterlist( 0 )( 0 ).set( "AZ_max_iter", 10000 );
        tParameterlist( 0 )( 0 ).set( "AZ_solver", AZ_gmres );
        tParameterlist( 0 )( 0 ).set( "AZ_subdomain_solve", AZ_ilu );
        tParameterlist( 0 )( 0 ).set( "AZ_graph_fill", 10 );
        tParameterlist( 0 )( 0 ).set( "ml_prec_type", "SA" );

        tParameterlist( 1 )( 0 ) = moris::prm::create_linear_solver_parameter_list();
        tParameterlist( 2 )( 0 ) = moris::prm::create_nonlinear_algorithm_parameter_list();
        tParameterlist( 3 )( 0 ) = moris::prm::create_nonlinear_solver_parameter_list();
        tParameterlist( 3 )( 0 ).set( "NLA_DofTypes", "UX,UY;TEMP" );

        tParameterlist( 4 )( 0 ) = moris::prm::create_time_solver_algorithm_parameter_list();
        tParameterlist( 5 )( 0 ) = moris::prm::create_time_solver_parameter_list();
        tParameterlist( 5 )( 0 ).set( "TSA_DofTypes", "UX,UY;TEMP" );
        tParameterlist( 5 )( 0 ).set( "TSA_Output_Indices", "0" );
        tParameterlist( 5 )( 0 ).set( "TSA_Output_Criteria", "Output_Criterion" );

        tParameterlist( 6 )( 0 ) = moris::prm::create_solver_warehouse_parameterlist();
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

        tParameterlist( 0 )( 0 ).set( "ghost_stab", true );

        tParameterlist( 0 )( 0 ).set( "multigrid", false );
    }

    void
    MSIParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_msi_parameter_list();
    }

    void
    VISParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        std::string tMorisOutput = std::getenv( "MORISOUTPUT" );

        MORIS_ERROR( tMorisOutput.size() > 0,
                "Environment variable MORISOUTPUT not set." );

        tParameterlist( 0 )( 0 ) = prm::create_vis_parameter_list();
        tParameterlist( 0 )( 0 ).set( "File_Name", std::pair< std::string, std::string >( tMorisOutput, "MDL_input_test.exo" ) );
        tParameterlist( 0 )( 0 ).set( "Set_Names", std::string( "HMR_dummy_n_p1" ) );
        tParameterlist( 0 )( 0 ).set( "Field_Names", std::string( "strain_energy_elemental,strain_energy_global,strain_energy_nodal_IP" ) );
        tParameterlist( 0 )( 0 ).set( "Field_Type", std::string( "ELEMENTAL,GLOBAL,NODAL" ) );
        tParameterlist( 0 )( 0 ).set( "IQI_Names", std::string( "IQIBulkU_1,IQIBulkU_1,IQIBulkU_1" ) );
    }

    void
    HMRParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_hmr_parameter_list();

        tParameterlist( 0 )( 0 ).set( "number_of_elements_per_dimension", std::string( "2, 1" ) );
        tParameterlist( 0 )( 0 ).set( "domain_dimensions", "2, 2" );
        tParameterlist( 0 )( 0 ).set( "domain_offset", "-1.0, -1.0" );
        tParameterlist( 0 )( 0 ).set( "domain_sidesets", "1,2,3,4" );
        tParameterlist( 0 )( 0 ).set( "lagrange_output_meshes", std::string( "0" ) );

        tParameterlist( 0 )( 0 ).set( "lagrange_orders", std::string( "1" ) );
        tParameterlist( 0 )( 0 ).set( "lagrange_pattern", std::string( "0" ) );
        tParameterlist( 0 )( 0 ).set( "bspline_orders", std::string( "1" ) );
        tParameterlist( 0 )( 0 ).set( "bspline_pattern", std::string( "0" ) );

        tParameterlist( 0 )( 0 ).set( "lagrange_to_bspline", "0" );

        tParameterlist( 0 )( 0 ).set( "truncate_bsplines", 1 );
        tParameterlist( 0 )( 0 ).set( "refinement_buffer", 3 );
        tParameterlist( 0 )( 0 ).set( "staircase_buffer", 3 );
        tParameterlist( 0 )( 0 ).set( "initial_refinement", "0" );
        tParameterlist( 0 )( 0 ).set( "initial_refinement_pattern", "0" );

        tParameterlist( 0 )( 0 ).set( "use_multigrid", 0 );
        tParameterlist( 0 )( 0 ).set( "severity_level", 2 );

        tParameterlist( 0 )( 0 ).set( "adaptive_refinement_level", 2 );
    }

    void
    GENParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 3 );
        tParameterlist( 0 ).resize( 1 );
        tParameterlist( 1 ).resize( 1 );

        // Main GEN parameter list
        tParameterlist( 0 )( 0 ) = prm::create_gen_parameter_list();

        // Geometry parameter lists
        tParameterlist( 1 )( 0 ) = prm::create_user_defined_geometry_parameter_list();
        tParameterlist( 1 )( 0 ).set( "field_function_name", "Lvl_set_1" );
        tParameterlist( 1 )( 0 ).set( "number_of_refinements", "2" );
        tParameterlist( 1 )( 0 ).set( "refinement_mesh_index", "0" );
    }

    void
    OPTParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 3 );
        tParameterlist( 0 ).resize( 1 );
        tParameterlist( 1 ).resize( 0 );
        tParameterlist( 2 ).resize( 0 );

        tParameterlist( 0 )( 0 ) = moris::prm::create_opt_problem_parameter_list();
        tParameterlist( 0 )( 0 ).set( "is_optimization_problem", false );
    }

    void
    MORISGENERALParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
    }

    //------------------------------------------------------------------------------
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
