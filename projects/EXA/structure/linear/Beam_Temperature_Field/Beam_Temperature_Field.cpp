/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Beam_Temperature_Field.cpp
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
#include "fn_PRM_STK_Parameters.hpp"
#include "fn_equal_to.hpp"
#include "fn_dot.hpp"
#include "fn_norm.hpp"
#include "paths.hpp"

#include "cl_Ascii.hpp"
#include "fn_Parsing_Tools.hpp"
#include "fn_assert.hpp"
#include "catch.hpp"

#include "AztecOO.h"

#ifdef __cplusplus
extern "C" {

#endif
//------------------------------------------------------------------------------
namespace moris
{
    // Input files
    std::string tPrefix         = moris::get_base_moris_dir();
    std::string tMeshFileName   = tPrefix + "/projects/EXA/structure/linear/Beam_Temperature_Field/beam_geometry_quadratic.e";
    std::string tOutputFileName = "Beam_Temperature_Field_Problem.exo";

    std::string tFieldRefPath = tPrefix + "/projects/EXA/structure/linear/Beam_Temperature_Field/beam_temperature_quadratic.e";

    /* ------------------------------------------------------------------------ */
    // material parameters

    std::string tEmod    = "1.0";
    std::string tPois    = "0.3";
    std::string tBedding = std::to_string( 1.0 * 1e-8 );
    std::string tCTE     = "0.01";
    std::string tRefTemp = "0.0";

    //------------------------------------------------------------------------------

    // FUNCTIONS

    // Constant function for properties
    void
    Func_Const(
            moris::Matrix< moris::DDRMat >                &aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > > &aParameters,
            moris::fem::Field_Interpolator_Manager        *aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    /* ------------------------------------------------------------------------ */

    void
    Func_Field(
            moris::Matrix< moris::DDRMat >                &aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > > &aParameters,
            moris::fem::Field_Interpolator_Manager        *aFIManager )
    {
        aPropMatrix = aFIManager->get_field_interpolators_for_type( mtk::Field_Type::FIELD_1 )->val();
    }

    bool
    Output_Criterion( moris::tsa::Time_Solver *aTimeSolver )
    {
        return true;
    }

    //------------------------------------------------------------------------------

    void
    OPTParameterList( Module_Parameter_Lists &aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_opt_problem_parameter_list() );

        aParameterLists.set( "is_optimization_problem", false );
        aParameterLists.set( "workflow", "STK_FEM" );
    }

    void
    STKParameterList( Module_Parameter_Lists &aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_stk_parameter_list() );
        aParameterLists.set( "input_file", tMeshFileName );
    }

    //------------------------------------------------------------------------------

    void
    FEMParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.hack_for_legacy_fem();
        // create a cell of cell of parameter list for fem

        //------------------------------------------------------------------------------

        aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropYoungs" );
        aParameterLists.set( "function_parameters", tEmod );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropPoisson" );
        aParameterLists.set( "function_parameters", tPois );
        aParameterLists.set( "value_function", "Func_Const" );

        // properties of boundary conditions
        aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropDirichlet" );
        aParameterLists.set( "function_parameters", "0.0;0.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // properties of bedding (supression for RBMs)
        aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropBedding" );
        aParameterLists.set( "function_parameters", tBedding );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropCTE" );
        aParameterLists.set( "function_parameters", tCTE );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropRefTemp" );
        aParameterLists.set( "function_parameters", tRefTemp );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropField" );
        aParameterLists.set( "value_function", "Func_Field" );
        aParameterLists.set( "field_dependencies", "FIELD_1" );

        //------------------------------------------------------------------------------

        // create parameter list for constitutive model 1
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists.set( "constitutive_name", "CMStrucLinIso1" );
        aParameterLists.set( "constitutive_type",  fem::Constitutive_Type::STRUC_LIN_ISO ) ;
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
        aParameterLists.set( "properties",
                "PropYoungs,YoungsModulus;"
                "PropPoisson,PoissonRatio;"
                "PropCTE,CTE;"
                "PropField,PropertyTemperature;"
                "PropRefTemp,ReferenceTemperature" );

        //------------------------------------------------------------------------------

        // create parameter list for stabilization parameter 1
        aParameterLists( FEM::STABILIZATION ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists.set( "stabilization_name", "SPNitsche" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::DIRICHLET_NITSCHE ) ;
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "leader_properties", "PropYoungs,Material" );

        //------------------------------------------------------------------------------
        // Bulk Terms
        aParameterLists( FEM::IWG ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGBulkU_1" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_BULK ) ;
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        aParameterLists.set( "leader_properties", "PropBedding,Bedding" );
        aParameterLists.set( "mesh_set_names", "block_main" );

        // Dirichlet Boundary
        aParameterLists( FEM::IWG ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGDirichletU" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE ) ;
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_properties", "PropDirichlet,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        aParameterLists.set( "stabilization_parameters", "SPNitsche,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", "sideset_left" );

        //------------------------------------------------------------------------------
        aParameterLists( FEM::IQI ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists.set( "IQI_name", "IQIBulkDISPX" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists.set( "dof_quantity", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "mesh_set_names", "block_main" );

        aParameterLists( FEM::IQI ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists.set( "IQI_name", "IQIBulkDISPY" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists.set( "dof_quantity", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "vectorial_field_index", 1 );
        aParameterLists.set( "mesh_set_names", "block_main" );

        // Nodal Temperature IQI
        aParameterLists( FEM::IQI ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists.set( "IQI_name", "IQITempField" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::PROPERTY ) ;
        aParameterLists.set( "leader_properties", "PropField,Property" );
        aParameterLists.set( "mesh_set_names", "block_main" );

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        aParameterLists( FEM::COMPUTATION ).add_parameter_list( prm::create_computation_parameter_list() );

        //------------------------------------------------------------------------------

        aParameterLists( FEM::FIELDS ).add_parameter_list( prm::create_fem_field_parameter_list() );
        aParameterLists.set( "field_name", "FieldNodalTEMP" );
        aParameterLists.set( "field_entity_type", "NODAL" );
        aParameterLists.set( "field_type", "FIELD_1" );
        aParameterLists.set( "field_create_from_file", tFieldRefPath );
    }

    //------------------------------------------------------------------------------

    void
    SOLParameterList( Module_Parameter_Lists &aParameterLists )
    {

        aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL ) );

        aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list( moris::prm::create_linear_solver_parameter_list() );

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
        aParameterLists.set( "NLA_rel_res_norm_drop", 1e-09 );
        aParameterLists.set( "NLA_relaxation_parameter", 1.0 );
        aParameterLists.set( "NLA_max_iter", 5 );
        aParameterLists.set( "NLA_combined_res_jac_assembly", false );

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
        aParameterLists.set( "NLA_DofTypes", "UX,UY" );

        aParameterLists( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list( moris::prm::create_time_solver_algorithm_parameter_list() );

        aParameterLists( SOL::TIME_SOLVERS ).add_parameter_list( moris::prm::create_time_solver_parameter_list() );
        aParameterLists.set( "TSA_DofTypes", "UX,UY" );
        aParameterLists.set( "TSA_Output_Indices", "0" );
        aParameterLists.set( "TSA_Output_Criteria", "Output_Criterion" );

        aParameterLists( SOL::SOLVER_WAREHOUSE ).add_parameter_list( moris::prm::create_solver_warehouse_parameterlist() );

        aParameterLists( SOL::PRECONDITIONERS ).add_parameter_list( moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::NONE ) );
    }

    //------------------------------------------------------------------------------

    void
    MSIParameterList( Module_Parameter_Lists &aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_msi_parameter_list() );
    }

    //------------------------------------------------------------------------------

    void
    VISParameterList( Module_Parameter_Lists &aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_vis_parameter_list() );
        aParameterLists.set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName ) );
        aParameterLists.set( "Mesh_Type",  vis::VIS_Mesh_Type::STANDARD ) ;
        aParameterLists.set( "Set_Names", "block_main" );
        aParameterLists.set( "Field_Names", "UX,UY,Temp" );
        aParameterLists.set( "Field_Type", "NODAL,NODAL,NODAL" );
        aParameterLists.set( "IQI_Names", "IQIBulkDISPX,IQIBulkDISPY,IQITempField" );
        aParameterLists.set( "Save_Frequency", 1 );
    }

    void
    MORISGENERALParameterList( Module_Parameter_Lists &aParameterLists )
    {
    }

    //------------------------------------------------------------------------------
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
