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
#include "fn_PRM_FEM_Parameters.hpp"
#include "fn_PRM_MSI_Parameters.hpp"
#include "fn_PRM_SOL_Parameters.hpp"
#include "fn_PRM_VIS_Parameters.hpp"
#include "fn_PRM_HMR_Parameters.hpp"
#include "fn_PRM_GEN_Parameters.hpp"
#include "fn_PRM_XTK_Parameters.hpp"
#include "fn_PRM_OPT_Parameters.hpp"
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
    OPTParameterList( Vector< Submodule_Parameter_Lists > &tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).add_parameter_list( prm::create_opt_problem_parameter_list() );

        tParameterlist( 0 ).set( "is_optimization_problem", false );
        tParameterlist( 0 ).set( "workflow", "STK_FEM" );
    }

    void
    STKParameterList( Vector< Submodule_Parameter_Lists > &tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).add_parameter_list( prm::create_stk_parameter_list() );
        tParameterlist( 0 ).set( "input_file", tMeshFileName );
    }

    //------------------------------------------------------------------------------

    void
    FEMParameterList( Vector< Submodule_Parameter_Lists > &tParameterList )
    {
        // create a cell of cell of parameter list for fem
        tParameterList.resize( 8 );

        //------------------------------------------------------------------------------

        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropYoungs" );
        tParameterList( 0 ).set( "function_parameters", tEmod );
        tParameterList( 0 ).set( "value_function", "Func_Const" );

        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropPoisson" );
        tParameterList( 0 ).set( "function_parameters", tPois );
        tParameterList( 0 ).set( "value_function", "Func_Const" );

        // properties of boundary conditions
        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropDirichlet" );
        tParameterList( 0 ).set( "function_parameters", "0.0;0.0" );
        tParameterList( 0 ).set( "value_function", "Func_Const" );

        // properties of bedding (supression for RBMs)
        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropBedding" );
        tParameterList( 0 ).set( "function_parameters", tBedding );
        tParameterList( 0 ).set( "value_function", "Func_Const" );

        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropCTE" );
        tParameterList( 0 ).set( "function_parameters", tCTE );
        tParameterList( 0 ).set( "value_function", "Func_Const" );

        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropRefTemp" );
        tParameterList( 0 ).set( "function_parameters", tRefTemp );
        tParameterList( 0 ).set( "value_function", "Func_Const" );

        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropField" );
        tParameterList( 0 ).set( "value_function", "Func_Field" );
        tParameterList( 0 ).set( "field_dependencies", "FIELD_1" );

        //------------------------------------------------------------------------------

        // create parameter list for constitutive model 1
        tParameterList( 1 ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 ).set( "constitutive_name", "CMStrucLinIso1" );
        tParameterList( 1 ).set( "constitutive_type",  fem::Constitutive_Type::STRUC_LIN_ISO ) ;
        tParameterList( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
        tParameterList( 1 ).set( "properties",
                "PropYoungs,YoungsModulus;"
                "PropPoisson,PoissonRatio;"
                "PropCTE,CTE;"
                "PropField,PropertyTemperature;"
                "PropRefTemp,ReferenceTemperature" );

        //------------------------------------------------------------------------------

        // create parameter list for stabilization parameter 1
        tParameterList( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 ).set( "stabilization_name", "SPNitsche" );
        tParameterList( 2 ).set( "stabilization_type",  fem::Stabilization_Type::DIRICHLET_NITSCHE ) ;
        tParameterList( 2 ).set( "function_parameters", "100.0" );
        tParameterList( 2 ).set( "leader_properties", "PropYoungs,Material" );

        //------------------------------------------------------------------------------
        // Bulk Terms
        tParameterList( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( 3 ).set( "IWG_name", "IWGBulkU_1" );
        tParameterList( 3 ).set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_BULK ) ;
        tParameterList( 3 ).set( "dof_residual", "UX,UY" );
        tParameterList( 3 ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 3 ).set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        tParameterList( 3 ).set( "leader_properties", "PropBedding,Bedding" );
        tParameterList( 3 ).set( "mesh_set_names", "block_main" );

        // Dirichlet Boundary
        tParameterList( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( 3 ).set( "IWG_name", "IWGDirichletU" );
        tParameterList( 3 ).set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE ) ;
        tParameterList( 3 ).set( "dof_residual", "UX,UY" );
        tParameterList( 3 ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 3 ).set( "leader_properties", "PropDirichlet,Dirichlet" );
        tParameterList( 3 ).set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        tParameterList( 3 ).set( "stabilization_parameters", "SPNitsche,DirichletNitsche" );
        tParameterList( 3 ).set( "mesh_set_names", "sideset_left" );

        //------------------------------------------------------------------------------
        tParameterList( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        tParameterList( 4 ).set( "IQI_name", "IQIBulkDISPX" );
        tParameterList( 4 ).set( "IQI_type",  fem::IQI_Type::DOF ) ;
        tParameterList( 4 ).set( "dof_quantity", "UX,UY" );
        tParameterList( 4 ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 4 ).set( "vectorial_field_index", 0 );
        tParameterList( 4 ).set( "mesh_set_names", "block_main" );

        tParameterList( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        tParameterList( 4 ).set( "IQI_name", "IQIBulkDISPY" );
        tParameterList( 4 ).set( "IQI_type",  fem::IQI_Type::DOF ) ;
        tParameterList( 4 ).set( "dof_quantity", "UX,UY" );
        tParameterList( 4 ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 4 ).set( "vectorial_field_index", 1 );
        tParameterList( 4 ).set( "mesh_set_names", "block_main" );

        // Nodal Temperature IQI
        tParameterList( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        tParameterList( 4 ).set( "IQI_name", "IQITempField" );
        tParameterList( 4 ).set( "IQI_type",  fem::IQI_Type::PROPERTY ) ;
        tParameterList( 4 ).set( "leader_properties", "PropField,Property" );
        tParameterList( 4 ).set( "mesh_set_names", "block_main" );

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        tParameterList( 5 ).add_parameter_list( prm::create_computation_parameter_list() );

        //------------------------------------------------------------------------------

        tParameterList( 6 ).add_parameter_list( prm::create_fem_field_parameter_list() );
        tParameterList( 6 ).set( "field_name", "FieldNodalTEMP" );
        tParameterList( 6 ).set( "field_entity_type", "NODAL" );
        tParameterList( 6 ).set( "field_type", "FIELD_1" );
        tParameterList( 6 ).set( "field_create_from_file", tFieldRefPath );
    }

    //------------------------------------------------------------------------------

    void
    SOLParameterList( Vector< Submodule_Parameter_Lists > &tParameterlist )
    {
        tParameterlist.resize( 8 );

        tParameterlist( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL ) );

        tParameterlist( 1 ).add_parameter_list( moris::prm::create_linear_solver_parameter_list() );

        tParameterlist( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
        tParameterlist( 2 ).set( "NLA_rel_res_norm_drop", 1e-09 );
        tParameterlist( 2 ).set( "NLA_relaxation_parameter", 1.0 );
        tParameterlist( 2 ).set( "NLA_max_iter", 5 );
        tParameterlist( 2 ).set( "NLA_combined_res_jac_assembly", false );

        tParameterlist( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
        tParameterlist( 3 ).set( "NLA_DofTypes", "UX,UY" );

        tParameterlist( 4 ).add_parameter_list( moris::prm::create_time_solver_algorithm_parameter_list() );

        tParameterlist( 5 ).add_parameter_list( moris::prm::create_time_solver_parameter_list() );
        tParameterlist( 5 ).set( "TSA_DofTypes", "UX,UY" );
        tParameterlist( 5 ).set( "TSA_Output_Indices", "0" );
        tParameterlist( 5 ).set( "TSA_Output_Criteria", "Output_Criterion" );

        tParameterlist( 6 ).add_parameter_list( moris::prm::create_solver_warehouse_parameterlist() );

        tParameterlist( 7 ).add_parameter_list( moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::NONE ) );
    }

    //------------------------------------------------------------------------------

    void
    MSIParameterList( Vector< Submodule_Parameter_Lists > &tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).add_parameter_list( prm::create_msi_parameter_list() );
    }

    //------------------------------------------------------------------------------

    void
    VISParameterList( Vector< Submodule_Parameter_Lists > &tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).add_parameter_list( prm::create_vis_parameter_list() );
        tParameterlist( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName ) );
        tParameterlist( 0 ).set( "Mesh_Type",  vis::VIS_Mesh_Type::STANDARD ) ;
        tParameterlist( 0 ).set( "Set_Names", "block_main" );
        tParameterlist( 0 ).set( "Field_Names", "UX,UY,Temp" );
        tParameterlist( 0 ).set( "Field_Type", "NODAL,NODAL,NODAL" );
        tParameterlist( 0 ).set( "IQI_Names", "IQIBulkDISPX,IQIBulkDISPY,IQITempField" );
        tParameterlist( 0 ).set( "Save_Frequency", 1 );
    }

    void
    MORISGENERALParameterList( Vector< Submodule_Parameter_Lists > &tParameterlist )
    {
    }

    //------------------------------------------------------------------------------
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
