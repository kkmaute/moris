/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Plane_Strain_Problem.cpp
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
    // Input file
    std::string tPrefix         = moris::get_base_moris_dir();
    std::string tMeshFileName   = tPrefix + "/projects/EXA/structure/linear/Plane_Strain/4elements.e";
    std::string tOutputFileName = "Plane_Strain_Problem.exo";

    /* ------------------------------------------------------------------------ */
    // material parameters

    std::string tEmod    = "1.0";
    std::string tPois    = "0.3";
    std::string tDens    = "1.0";
    std::string tBedding = std::to_string( 1.0 * 1e-6 );

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
    Func_Select_X(
            moris::Matrix< moris::DDRMat >                &aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > > &aParameters,
            moris::fem::Field_Interpolator_Manager        *aFIManager )
    {
        aPropMatrix.set_size( 2, 2, 0.0 );
        aPropMatrix( 0, 0 ) = 1.0;
    }

    /* ------------------------------------------------------------------------ */

    void
    Func_Select_Y(
            moris::Matrix< moris::DDRMat >                &aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > > &aParameters,
            moris::fem::Field_Interpolator_Manager        *aFIManager )
    {
        aPropMatrix.set_size( 2, 2, 0.0 );
        aPropMatrix( 1, 1 ) = 1.0;
    }

    bool
    Output_Criterion( moris::tsa::Time_Solver *aTimeSolver )
    {
        return true;
    }

    //------------------------------------------------------------------------------

    void
    OPTParameterList( Vector< Vector< ParameterList > > &aParameterLists )
    {
        aParameterLists( 0 ).push_back( prm::create_opt_problem_parameter_list() );

        aParameterLists( 0 ).set( "is_optimization_problem", false );
        aParameterLists( 0 ).set( "workflow", "STK_FEM" );
    }

    void
    STKParameterList( Vector< Vector< ParameterList > > &aParameterLists )
    {
        aParameterLists( 0 ).push_back( prm::create_stk_parameter_list() );
        aParameterLists( 0 ).set( "input_file", tMeshFileName );
    }

    //------------------------------------------------------------------------------

    void
    FEMParameterList( Vector< Vector< ParameterList > > &aParameterLists )
    {
        // create a cell of cell of parameter list for fem

        //------------------------------------------------------------------------------

        // properties of bars
        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropDensity" );
        aParameterLists( 0 ).set( "function_parameters", tDens );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropYoungs" );
        aParameterLists( 0 ).set( "function_parameters", tEmod );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropPoisson" );
        aParameterLists( 0 ).set( "function_parameters", tPois );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // properties of boundary conditions
        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropDirichlet" );
        aParameterLists( 0 ).set( "function_parameters", "0.0;0.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropSelectX" );
        aParameterLists( 0 ).set( "function_parameters", "0.0;0.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Select_X" );

        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropSelectY" );
        aParameterLists( 0 ).set( "function_parameters", "0.0;0.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Select_Y" );

        // create parameter list for property 10
        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropNeumann" );
        aParameterLists( 0 ).set( "function_parameters", "0.0;-1.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropThickness" );
        aParameterLists( 0 ).set( "function_parameters", "0.1" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // properties of bedding (supression for RBMs)
        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropBedding" );
        aParameterLists( 0 ).set( "function_parameters", tBedding );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        //------------------------------------------------------------------------------
        // create parameter list for constitutive model 1
        aParameterLists( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        aParameterLists( 1 ).set( "constitutive_name", "CMStrucLinIso1" );
        aParameterLists( 1 ).set( "constitutive_type",  fem::Constitutive_Type::STRUC_LIN_ISO ) ;
        aParameterLists( 1 ).set( "model_type",  fem::Model_Type::PLANE_STRAIN ) ;
        aParameterLists( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
        aParameterLists( 1 ).set( "properties",
                "PropYoungs,YoungsModulus;"
                "PropPoisson,PoissonRatio" );

        //------------------------------------------------------------------------------

        // create parameter list for stabilization parameter 1
        aParameterLists( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPNitsche" );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::DIRICHLET_NITSCHE ) ;
        aParameterLists( 2 ).set( "function_parameters", "100.0" );
        aParameterLists( 2 ).set( "leader_properties", "PropYoungs,Material" );

        //------------------------------------------------------------------------------
        // Bulk Terms
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGBulkU_1" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_BULK ) ;
        aParameterLists( 3 ).set( "dof_residual", "UX,UY" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        aParameterLists( 3 ).set( "leader_properties", "PropBedding,Bedding;PropThickness,Thickness" );
        aParameterLists( 3 ).set( "mesh_set_names", "block_main" );

        // create IWG - Dirichlet structure
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGDirichletY" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE ) ;
        aParameterLists( 3 ).set( "dof_residual", "UX,UY" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists( 3 ).set( "leader_properties", "PropDirichlet,Dirichlet;PropSelectY,Select;PropThickness,Thickness" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPNitsche,DirichletNitsche" );
        aParameterLists( 3 ).set( "mesh_set_names", "sideset_bottom" );

        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGDirichletX" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE ) ;
        aParameterLists( 3 ).set( "dof_residual", "UX,UY" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists( 3 ).set( "leader_properties", "PropDirichlet,Dirichlet;PropSelectX,Select;PropThickness,Thickness" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPNitsche,DirichletNitsche" );
        aParameterLists( 3 ).set( "mesh_set_names", "sideset_left" );

        // Neumann Boundary
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGNeumannFlux" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_NEUMANN ) ;
        aParameterLists( 3 ).set( "dof_residual", "UX,UY" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists( 3 ).set( "leader_properties", "PropNeumann,Traction;PropThickness,Thickness" );
        aParameterLists( 3 ).set( "mesh_set_names", "sideset_top" );

        //------------------------------------------------------------------------------
        aParameterLists( 4 ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkDISPX" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists( 4 ).set( "dof_quantity", "UX,UY" );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists( 4 ).set( "vectorial_field_index", 0 );
        aParameterLists( 4 ).set( "mesh_set_names", "block_main" );

        aParameterLists( 4 ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkDISPY" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists( 4 ).set( "dof_quantity", "UX,UY" );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists( 4 ).set( "vectorial_field_index", 1 );
        aParameterLists( 4 ).set( "mesh_set_names", "block_main" );

        // nodal von-mises stresses for shell
        aParameterLists( 4 ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIVMStress" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::VON_MISES_STRESS ) ;
        aParameterLists( 4 ).set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists( 4 ).set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        aParameterLists( 4 ).set( "mesh_set_names", "block_main" );

        // nodal stresses for shell
        aParameterLists( 4 ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIXStress" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::NORMAL_STRESS ) ;
        aParameterLists( 4 ).set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists( 4 ).set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        aParameterLists( 4 ).set( "vectorial_field_index", 0 );
        aParameterLists( 4 ).set( "mesh_set_names", "block_main" );

        // nodal von-mises stresses for shell
        aParameterLists( 4 ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIYStress" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::NORMAL_STRESS ) ;
        aParameterLists( 4 ).set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists( 4 ).set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        aParameterLists( 4 ).set( "vectorial_field_index", 1 );
        aParameterLists( 4 ).set( "mesh_set_names", "block_main" );

        // nodal von-mises stresses for shell
        aParameterLists( 4 ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIOOPStress" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::NORMAL_STRESS ) ;
        aParameterLists( 4 ).set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists( 4 ).set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        aParameterLists( 4 ).set( "vectorial_field_index", 2 );
        aParameterLists( 4 ).set( "mesh_set_names", "block_main" );

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        aParameterLists( 5 ).push_back( prm::create_computation_parameter_list() );
    }

    //------------------------------------------------------------------------------

    void
    SOLParameterList( Vector< Vector< ParameterList > > &aParameterLists )
    {

        aParameterLists( 0 ).push_back( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL ) );

        aParameterLists( 1 ).push_back( moris::prm::create_linear_solver_parameter_list() );

        aParameterLists( 2 ).push_back( moris::prm::create_nonlinear_algorithm_parameter_list() );
        aParameterLists( 2 ).set( "NLA_rel_res_norm_drop", 1e-09 );
        aParameterLists( 2 ).set( "NLA_relaxation_parameter", 1.0 );
        aParameterLists( 2 ).set( "NLA_max_iter", 5 );
        aParameterLists( 2 ).set( "NLA_combined_res_jac_assembly", false );

        aParameterLists( 3 ).push_back( moris::prm::create_nonlinear_solver_parameter_list() );
        aParameterLists( 3 ).set( "NLA_DofTypes", "UX,UY" );

        aParameterLists( 4 ).push_back( moris::prm::create_time_solver_algorithm_parameter_list() );

        aParameterLists( 5 ).push_back( moris::prm::create_time_solver_parameter_list() );
        aParameterLists( 5 ).set( "TSA_DofTypes", "UX,UY" );
        aParameterLists( 5 ).set( "TSA_Output_Indices", "0" );
        aParameterLists( 5 ).set( "TSA_Output_Criteria", "Output_Criterion" );

        aParameterLists( 6 ).push_back( moris::prm::create_solver_warehouse_parameterlist() );

        aParameterLists( 7 ).push_back( moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::NONE ) );
    }

    //------------------------------------------------------------------------------

    void
    MSIParameterList( Vector< Vector< ParameterList > > &aParameterLists )
    {
        aParameterLists( 0 ).push_back( prm::create_msi_parameter_list() );
    }

    //------------------------------------------------------------------------------

    void
    VISParameterList( Vector< Vector< ParameterList > > &aParameterLists )
    {
        aParameterLists( 0 ).push_back( prm::create_vis_parameter_list() );
        aParameterLists( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName ) );
        aParameterLists( 0 ).set( "Mesh_Type",  vis::VIS_Mesh_Type::STANDARD ) ;
        aParameterLists( 0 ).set( "Set_Names", "block_main" );
        aParameterLists( 0 ).set( "Field_Names", "UX,UY,STRESSVM,STRESSX,STRESSY,STRESSOOP" );
        aParameterLists( 0 ).set( "Field_Type", "NODAL,NODAL,NODAL,NODAL,NODAL,NODAL" );
        aParameterLists( 0 ).set( "IQI_Names", "IQIBulkDISPX,IQIBulkDISPY,IQIVMStress,IQIXStress,IQIYStress,IQIOOPStress" );
        aParameterLists( 0 ).set( "Save_Frequency", 1 );
    }

    void
    MORISGENERALParameterList( Vector< Vector< ParameterList > > &aParameterLists )
    {
    }

    //------------------------------------------------------------------------------
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
