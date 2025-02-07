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
#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_MSI_Equation_Object.hpp"
#include "cl_TSA_Time_Solver.hpp"
#include "cl_DLA_Solver_Interface.hpp"
#include "fn_equal_to.hpp"
#include "parameters.hpp"

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
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
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
            const Vector< moris::real* >& aGeometryParameters )
    {
        return 1.01;
    }

    void
    FEMParameterList( Module_Parameter_Lists& aParameterList )
    {
        aParameterList.hack_for_legacy_fem();
        // create parameter list for property 1
        aParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterList( 0 ).set( "property_name", "PropConductivity1" );
        aParameterList( 0 ).set( "function_parameters", "1.0" );
        aParameterList( 0 ).set( "value_function", "Func1" );

        // create parameter list for property 2
        aParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterList( 0 ).set( "property_name", "PropConductivity2" );
        aParameterList( 0 ).set( "function_parameters", "5.0" );
        aParameterList( 0 ).set( "value_function", "Func1" );

        // create parameter list for property 3
        aParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterList( 0 ).set( "property_name", "PropDirichletU" );
        aParameterList( 0 ).set( "function_parameters", "0.0;0.0" );
        aParameterList( 0 ).set( "value_function", "Func1" );

        // create parameter list for property 4
        aParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterList( 0 ).set( "property_name", "PropDirichletTEMP" );
        aParameterList( 0 ).set( "function_parameters", "3.0" );
        aParameterList( 0 ).set( "value_function", "Func1" );

        // create parameter list for property 5
        aParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterList( 0 ).set( "property_name", "PropEMod1" );
        aParameterList( 0 ).set( "function_parameters", "1.0" );
        aParameterList( 0 ).set( "value_function", "Func1" );

        // create parameter list for property 6
        aParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterList( 0 ).set( "property_name", "PropEMod2" );
        aParameterList( 0 ).set( "function_parameters", "1.0" );
        aParameterList( 0 ).set( "value_function", "Func1" );

        // create parameter list for property 7
        aParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterList( 0 ).set( "property_name", "PropPoisson" );
        aParameterList( 0 ).set( "function_parameters", "0.0" );
        aParameterList( 0 ).set( "value_function", "Func1" );

        // create parameter list for property 8
        aParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterList( 0 ).set( "property_name", "PropCTE" );
        aParameterList( 0 ).set( "function_parameters", "1.0" );
        aParameterList( 0 ).set( "value_function", "Func1" );

        // create parameter list for property 9
        aParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterList( 0 ).set( "property_name", "PropTRef" );
        aParameterList( 0 ).set( "function_parameters", "1.0" );
        aParameterList( 0 ).set( "value_function", "Func1" );

        //------------------------------------------------------------------------------
        // fill the constitutive model part of the parameter list

        // create parameter list for constitutive model 1
        aParameterList( 1 ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterList( 1 ).set( "constitutive_name", "CMStrucLinIso1" );
        aParameterList( 1 ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::STRUC_LIN_ISO ) );
        aParameterList( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY;TEMP", "Displacement,Temperature" ) );
        aParameterList( 1 ).set( "properties", "PropEMod1,YoungsModulus;PropPoisson,PoissonRatio;PropCTE,CTE;PropTRef,ReferenceTemperature" );
        aParameterList( 1 ).set( "model_type", static_cast< uint >( fem::Model_Type::PLANE_STRESS ) );

        // create parameter list for constitutive model 2
        aParameterList( 1 ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterList( 1 ).set( "constitutive_name", "CMStrucLinIso2" );
        aParameterList( 1 ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::STRUC_LIN_ISO ) );
        aParameterList( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY;TEMP", "Displacement,Temperature" ) );
        aParameterList( 1 ).set( "properties", "PropEMod2,YoungsModulus;PropPoisson,PoissonRatio;PropCTE,CTE;PropTRef,ReferenceTemperature" );
        aParameterList( 1 ).set( "model_type", static_cast< uint >( fem::Model_Type::PLANE_STRESS ) );

        // create parameter list for constitutive model 3
        aParameterList( 1 ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterList( 1 ).set( "constitutive_name", "CMDiffLinIso1" );
        aParameterList( 1 ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::DIFF_LIN_ISO ) );
        aParameterList( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterList( 1 ).set( "properties", "PropConductivity1,Conductivity" );

        // create parameter list for constitutive model 4
        aParameterList( 1 ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterList( 1 ).set( "constitutive_name", "CMDiffLinIso2" );
        aParameterList( 1 ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::DIFF_LIN_ISO ) );
        aParameterList( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterList( 1 ).set( "properties", "PropConductivity2,Conductivity" );

        //------------------------------------------------------------------------------
        // fill the stabilization parameter part of the parameter list

        // create parameter list for stabilization parameter 1
        aParameterList( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterList( 2 ).set( "stabilization_name", "SPDirichletNitscheU" );
        aParameterList( 2 ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::DIRICHLET_NITSCHE ) );
        aParameterList( 2 ).set( "function_parameters", "100.0" );
        aParameterList( 2 ).set( "leader_properties", "PropEMod1,Material" );

        // create parameter list for stabilization parameter 2
        aParameterList( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterList( 2 ).set( "stabilization_name", "SPDirichletNitscheTEMP" );
        aParameterList( 2 ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::DIRICHLET_NITSCHE ) );
        aParameterList( 2 ).set( "function_parameters", "1.0" );
        aParameterList( 2 ).set( "leader_properties", "PropConductivity1,Material" );

        //------------------------------------------------------------------------------
        // fill the IWG part of the parameter list

        // create parameter list for IWG 1
        aParameterList( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterList( 3 ).set( "IWG_name", "IWGBulkU_1" );
        aParameterList( 3 ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_BULK ) );
        aParameterList( 3 ).set( "dof_residual", "UX,UY" );
        aParameterList( 3 ).set( "leader_dof_dependencies", "UX,UY" );
        aParameterList( 3 ).set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        aParameterList( 3 ).set( "mesh_set_names", "HMR_dummy_n_p1" );

        // create parameter list for IWG 2
        aParameterList( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterList( 3 ).set( "IWG_name", "IWGBulkU_2" );
        aParameterList( 3 ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_BULK ) );
        aParameterList( 3 ).set( "dof_residual", "UX,UY" );
        aParameterList( 3 ).set( "leader_dof_dependencies", "UX,UY" );
        aParameterList( 3 ).set( "leader_constitutive_models", "CMStrucLinIso2,ElastLinIso" );
        aParameterList( 3 ).set( "mesh_set_names", "" );

        // create parameter list for IWG 3
        aParameterList( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterList( 3 ).set( "IWG_name", "IWGDirichletU" );
        aParameterList( 3 ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE ) );
        aParameterList( 3 ).set( "dof_residual", "UX,UY" );
        aParameterList( 3 ).set( "leader_dof_dependencies", "UX,UY" );
        aParameterList( 3 ).set( "leader_properties", "PropDirichletU,Dirichlet" );
        aParameterList( 3 ).set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        aParameterList( 3 ).set( "stabilization_parameters", "SPDirichletNitscheU,DirichletNitsche" );
        aParameterList( 3 ).set( "mesh_set_names", "SideSet_4_n_p1" );

        // create parameter list for IWG 4
        aParameterList( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterList( 3 ).set( "IWG_name", "IWGBulkTEMP_1" );
        aParameterList( 3 ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::SPATIALDIFF_BULK ) );
        aParameterList( 3 ).set( "dof_residual", "TEMP" );
        aParameterList( 3 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterList( 3 ).set( "leader_constitutive_models", "CMDiffLinIso1,Diffusion" );
        aParameterList( 3 ).set( "mesh_set_names", "HMR_dummy_n_p1" );

        // create parameter list for IWG 5
        aParameterList( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterList( 3 ).set( "IWG_name", "IWGBulkTEMP_2" );
        aParameterList( 3 ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::SPATIALDIFF_BULK ) );
        aParameterList( 3 ).set( "dof_residual", "TEMP" );
        aParameterList( 3 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterList( 3 ).set( "leader_constitutive_models", "CMDiffLinIso2,Diffusion" );
        aParameterList( 3 ).set( "mesh_set_names", "" );

        // create parameter list for IWG 6
        aParameterList( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterList( 3 ).set( "IWG_name", "IWGDirichletTEMP" );
        aParameterList( 3 ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE ) );
        aParameterList( 3 ).set( "dof_residual", "TEMP" );
        aParameterList( 3 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterList( 3 ).set( "leader_properties", "PropDirichletTEMP,Dirichlet" );
        aParameterList( 3 ).set( "leader_constitutive_models", "CMDiffLinIso2,Diffusion" );
        aParameterList( 3 ).set( "stabilization_parameters", "SPDirichletNitscheTEMP,DirichletNitsche" );
        aParameterList( 3 ).set( "mesh_set_names", "SideSet_4_n_p1" );

        //------------------------------------------------------------------------------
        // fill the IQI part of the parameter list

        // create parameter list for IQI 1
        aParameterList( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterList( 4 ).set( "IQI_name", "IQIBulkU_1" );
        aParameterList( 4 ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::STRAIN_ENERGY ) );
        aParameterList( 4 ).set( "leader_dof_dependencies", "UX,UY" );
        aParameterList( 4 ).set( "leader_constitutive_models", "CMStrucLinIso1,Elast" );
        aParameterList( 4 ).set( "mesh_set_names", "HMR_dummy_n_p1" );
    }

    void
    SOLParameterList( Module_Parameter_Lists& aParameterList )
    {
        aParameterList( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AZTEC_IMPL ) );
        aParameterList( 0 ).set( "AZ_diagnostics", AZ_none );
        aParameterList( 0 ).set( "AZ_output", AZ_none );
        aParameterList( 0 ).set( "AZ_max_iter", 10000 );
        aParameterList( 0 ).set( "AZ_solver", AZ_gmres );
        aParameterList( 0 ).set( "AZ_subdomain_solve", AZ_ilu );
        aParameterList( 0 ).set( "AZ_graph_fill", 10 );
        aParameterList( 0 ).set( "preconditioners", "0" );

        aParameterList( 1 ).add_parameter_list( moris::prm::create_linear_solver_parameter_list() );
        aParameterList( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
        aParameterList( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
        aParameterList( 3 ).set( "NLA_DofTypes", "UX,UY;TEMP" );

        aParameterList( 4 ).add_parameter_list( moris::prm::create_time_solver_algorithm_parameter_list() );
        aParameterList( 5 ).add_parameter_list( moris::prm::create_time_solver_parameter_list() );
        aParameterList( 5 ).set( "TSA_DofTypes", "UX,UY;TEMP" );
        aParameterList( 5 ).set( "TSA_Output_Indices", "0" );
        aParameterList( 5 ).set( "TSA_Output_Criteria", "Output_Criterion" );

        aParameterList( 7 ).add_parameter_list( moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::ML ) );
        aParameterList( 7 ).set( "ml_prec_type", "SA" );
    }

    void
    XTKParameterList( Module_Parameter_Lists& aParameterList )
    {
        aParameterList( 0 ).set( "decompose", true );
        aParameterList( 0 ).set( "decomposition_type", "conformal" );

        aParameterList( 0 ).set( "enrich", true );
        aParameterList( 0 ).set( "basis_rank", "bspline" );
        aParameterList( 0 ).set( "enrich_mesh_indices", "0" );

        aParameterList( 0 ).set( "ghost_stab", true );

        aParameterList( 0 ).set( "multigrid", false );
    }

    void
    MSIParameterList( Module_Parameter_Lists& aParameterList )
    {
    }

    void
    VISParameterList( Module_Parameter_Lists& aParameterList )
    {

        std::string tMorisOutput = std::getenv( "MORISOUTPUT" );

        MORIS_ERROR( tMorisOutput.size() > 0,
                "Environment variable MORISOUTPUT not set." );
        aParameterList( 0 ).set( "File_Name", std::pair< std::string, std::string >( tMorisOutput, "MDL_input_test.exo" ) );
        aParameterList( 0 ).set( "Set_Names", "HMR_dummy_n_p1" );
        aParameterList( 0 ).set( "Field_Names", "strain_energy_elemental,strain_energy_global,strain_energy_nodal_IP" );
        aParameterList( 0 ).set( "Field_Type", "ELEMENTAL_AVG,GLOBAL,NODAL" );
        aParameterList( 0 ).set( "IQI_Names", "IQIBulkU_1,IQIBulkU_1,IQIBulkU_1" );
    }

    void
    HMRParameterList( Module_Parameter_Lists& aParameterList )
    {
        aParameterList( 0 ).set( "number_of_elements_per_dimension", 2, 1 );
        aParameterList( 0 ).set( "domain_dimensions", 2.0, 2.0 );
        aParameterList( 0 ).set( "domain_offset", -1.0, -1.0 );
        aParameterList( 0 ).set( "lagrange_output_meshes", "0" );

        aParameterList( 0 ).set( "lagrange_orders", "1" );
        aParameterList( 0 ).set( "lagrange_pattern", "0" );
        aParameterList( 0 ).set( "bspline_orders", "1" );
        aParameterList( 0 ).set( "bspline_pattern", "0" );

        aParameterList( 0 ).set( "lagrange_to_bspline", "0" );

        aParameterList( 0 ).set( "refinement_buffer", 3 );
        aParameterList( 0 ).set( "staircase_buffer", 3 );
        aParameterList( 0 ).set( "initial_refinement", "0" );
        aParameterList( 0 ).set( "initial_refinement_pattern", "0" );

        aParameterList( 0 ).set( "use_multigrid", 0 );
        aParameterList( 0 ).set( "severity_level", 2 );

        aParameterList( 0 ).set( "adaptive_refinement_level", 2 );
    }

    void
    GENParameterList( Module_Parameter_Lists& aParameterList )
    {
        // Geometry parameter lists
        aParameterList( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterList( 1 ).set( "field_function_name", "Lvl_set_1" );
        aParameterList( 1 ).set( "number_of_refinements", 2u );
        aParameterList( 1 ).set( "refinement_mesh_index", 0u );
    }

    void
    OPaParameterList( Module_Parameter_Lists& aParameterList )
    {
        aParameterList( 0 ).set( "is_optimization_problem", false );
    }

    void
    MORISGENERALParameterList( Module_Parameter_Lists& aParameterList )
    {
    }

    //------------------------------------------------------------------------------
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
