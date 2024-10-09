/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Axisymmetric_Problem.cpp
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
    std::string tMeshFileName   = tPrefix + "/projects/EXA/structure/linear/Axisymmetric/cylinder_geometry.e";
    std::string tOutputFileName = "Axisymmetric_Problem.exo";

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
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    // calculates the outward normal vector for distance from rotation axis to point of interest
    // and the radius so aPropMatrix = {{2*pi*r},{r},{n1},{n2}}
    void
    Func_AxiRotation(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        MORIS_ASSERT( aParameters( 0 ).n_cols() == 2 and aParameters( 0 ).n_rows() == 2,
                "Axisymmetric rotation axis incorrectly defined. use {{x1,y1},{x2,y2}}." );

        // vector from point 1 to point 2
        moris::Matrix< moris::DDRMat > tRotVec = { //
            { aParameters( 0 )( 1, 0 ) - aParameters( 0 )( 0, 0 ) },
            { aParameters( 0 )( 1, 1 ) - aParameters( 0 )( 0, 1 ) }
        };

        // spatial location of interest
        moris::Matrix< moris::DDRMat > tX = aFIManager->get_IP_geometry_interpolator()->valx();

        // vector from point 1 to location of interest
        moris::Matrix< moris::DDRMat > tPntVec = { //
            { tX( 0 ) - aParameters( 0 )( 0, 0 ) },
            { tX( 1 ) - aParameters( 0 )( 0, 1 ) }
        };

        // minimum distance vector from line to point of interest
        moris::Matrix< moris::DDRMat > tRadVec =                 //
                tPntVec - dot( tPntVec, tRotVec ) * tRotVec /    //
                                  ( tRotVec( 0 ) * tRotVec( 0 ) + tRotVec( 1 ) * tRotVec( 1 ) );

        // init aPropMatrixSize
        aPropMatrix.set_size( 4, 1, 0.0 );

        // radius
        aPropMatrix( 1 ) = std::max( norm( tRadVec ), MORIS_REAL_EPS );

        // 2*pi*r for IWG
        aPropMatrix( 0 ) = aPropMatrix( 1 ) * 4.0 * std::acos( 0 );

        // outward radial normal vector
        aPropMatrix( { 2, 3 } ) = tRadVec / aPropMatrix( 1 );
    }

    /* ------------------------------------------------------------------------ */

    void
    Func_Select_X(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix.set_size( 2, 2, 0.0 );
        aPropMatrix( 0, 0 ) = 1.0;
    }

    /* ------------------------------------------------------------------------ */

    void
    Func_Select_Y(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix.set_size( 2, 2, 0.0 );
        aPropMatrix( 1, 1 ) = 1.0;
    }

    bool
    Output_Criterion( moris::tsa::Time_Solver* aTimeSolver )
    {
        return true;
    }

    //------------------------------------------------------------------------------

    void
    OPTParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "is_optimization_problem", false );
        aParameterLists.set( "workflow", "STK_FEM" );
    }

    void
    STKParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "input_file", tMeshFileName );
    }

    //------------------------------------------------------------------------------

    void
    FEMParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.hack_for_legacy_fem();
        // create a cell of cell of parameter list for fem

        //------------------------------------------------------------------------------

        // properties of bars
        aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropDensity" );
        aParameterLists.set( "function_parameters", tDens );
        aParameterLists.set( "value_function", "Func_Const" );

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

        aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropSelectX" );
        aParameterLists.set( "function_parameters", "0.0;0.0" );
        aParameterLists.set( "value_function", "Func_Select_X" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropSelectY" );
        aParameterLists.set( "function_parameters", "0.0;0.0" );
        aParameterLists.set( "value_function", "Func_Select_Y" );

        // create parameter list for property 10
        aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropNeumann" );
        aParameterLists.set( "function_parameters", "0.0;-1.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // create rotation axis property
        aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropRotAxis" );
        aParameterLists.set( "function_parameters", "0.0,-1.0;1.0,-1.0" );
        aParameterLists.set( "value_function", "Func_AxiRotation" );

        // properties of bedding (supression for RBMs)
        aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropBedding" );
        aParameterLists.set( "function_parameters", tBedding );
        aParameterLists.set( "value_function", "Func_Const" );

        //------------------------------------------------------------------------------
        // create parameter list for constitutive model 1
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists.set( "constitutive_name", "CMStrucLinIso1" );
        aParameterLists.set( "constitutive_type",  fem::Constitutive_Type::STRUC_LIN_ISO ) ;
        aParameterLists.set( "model_type",  fem::Model_Type::AXISYMMETRIC ) ;
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
        aParameterLists.set( "properties",
                "PropYoungs,YoungsModulus;"
                "PropPoisson,PoissonRatio;"
                "PropRotAxis,AxisymRotationAxis" );

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
        aParameterLists.set( "leader_properties", "PropBedding,Bedding;PropRotAxis,Thickness" );
        aParameterLists.set( "mesh_set_names", "block_main" );

        // create IWG - Dirichlet structure
        aParameterLists( FEM::IWG ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGDirichletX" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE ) ;
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_properties", "PropDirichlet,Dirichlet;PropSelectX,Select;PropRotAxis,Thickness" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        aParameterLists.set( "stabilization_parameters", "SPNitsche,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", "sideset_left" );

        // Neumann Boundary
        aParameterLists( FEM::IWG ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGNeumannFlux" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_NEUMANN ) ;
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_properties", "PropNeumann,Traction;PropRotAxis,Thickness" );
        aParameterLists.set( "mesh_set_names", "sideset_top" );

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

        // nodal von-mises stresses for shell
        aParameterLists( FEM::IQI ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists.set( "IQI_name", "IQIVMStress" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::VON_MISES_STRESS ) ;
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        aParameterLists.set( "mesh_set_names", "block_main" );

        // nodal stresses for shell
        aParameterLists( FEM::IQI ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists.set( "IQI_name", "IQIXStress" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::NORMAL_STRESS ) ;
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "mesh_set_names", "block_main" );

        // nodal von-mises stresses for shell
        aParameterLists( FEM::IQI ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists.set( "IQI_name", "IQIRadialStress" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::NORMAL_STRESS ) ;
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        aParameterLists.set( "vectorial_field_index", 1 );
        aParameterLists.set( "mesh_set_names", "block_main" );

        // nodal von-mises stresses for shell
        aParameterLists( FEM::IQI ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists.set( "IQI_name", "IQIOOPStress" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::NORMAL_STRESS ) ;
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        aParameterLists.set( "vectorial_field_index", 2 );
        aParameterLists.set( "mesh_set_names", "block_main" );

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        aParameterLists( FEM::COMPUTATION );
    }

    //------------------------------------------------------------------------------

    void
    SOLParameterList( Module_Parameter_Lists& aParameterLists )
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

        aParameterLists( SOL::PRECONDITIONERS ).add_parameter_list( moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::NONE ) );
    }

    //------------------------------------------------------------------------------

    void
    MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
    }

    //------------------------------------------------------------------------------

    void
    VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName ) );
        aParameterLists.set( "Mesh_Type",  vis::VIS_Mesh_Type::STANDARD ) ;
        aParameterLists.set( "Set_Names", "block_main" );
        aParameterLists.set( "Field_Names", "UX,UY,STRESSVM,STRESSX,STRESSRAD,STRESSOOP" );
        aParameterLists.set( "Field_Type", "NODAL,NODAL,NODAL,NODAL,NODAL,NODAL" );
        aParameterLists.set( "IQI_Names", "IQIBulkDISPX,IQIBulkDISPY,IQIVMStress,IQIXStress,IQIRadialStress,IQIOOPStress" );
        aParameterLists.set( "Save_Frequency", 1 );
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
