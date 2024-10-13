/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Shape_Sensitivity_Sweep.cpp
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
#include "fn_equal_to.hpp"

#include "AztecOO.h"

#ifdef __cplusplus
extern "C" {
#endif
//------------------------------------------------------------------------------
namespace moris
{
    std::string tMeshSets = "HMR_dummy_n_p1,HMR_dummy_c_p1";
    std::string tDBCSets  = "iside_b0_1_b1_3";
    std::string tNBCSets  = "SideSet_3_c_p1";    ///"iside_b0_1_b1_0";

    // FD in adjoint
    real tFEMFdEpsilon = 1.0e-7;

    // FD step size in sweep
    std::string tFDsweep = "1.0e-7";

    //--------------------------------------------------------------------------------------------------------------

    // Constant function for properties
    void
    Func_Const( moris::Matrix< moris::DDRMat >       &aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > > &aParameters,
            moris::fem::Field_Interpolator_Manager   *aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    bool
    Output_Criterion( moris::tsa::Time_Solver *aTimeSolver )
    {
        return true;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDSMat >
    get_constraint_types()
    {
        Matrix< DDSMat > tConstraintTypes( 1, 1, 1 );

        return tConstraintTypes;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_objectives( const Vector< real > &aADVs, const Vector< real > &aCriteria )
    {
        Matrix< DDRMat > tObjectives = { { aCriteria( 0 ) } };

        return tObjectives;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_constraints( const Vector< real > &aADVs, const Vector< real > &aCriteria )
    {
        Matrix< DDRMat > tConstraints = { { aCriteria( 1 ) } };

        return tConstraints;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_dobjective_dadv( const Vector< real > &aADVs, const Vector< real > &aCriteria )
    {
        Matrix< DDRMat > tDObjectiveDADV( 1, aADVs.size(), 0.0 );

        return tDObjectiveDADV;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_dobjective_dcriteria( const Vector< real > &aADVs, const Vector< real > &aCriteria )
    {
        Matrix< DDRMat > tDObjectiveDCriteria( 1, 2 );
        tDObjectiveDCriteria( 0 ) = 1;
        tDObjectiveDCriteria( 1 ) = 0;

        return tDObjectiveDCriteria;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_dconstraint_dadv( const Vector< real > &aADVs, const Vector< real > &aCriteria )
    {
        Matrix< DDRMat > tDConstraintDADV( 1, aADVs.size(), 0.0 );

        return tDConstraintDADV;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_dconstraint_dcriteria( const Vector< real > &aADVs, const Vector< real > &aCriteria )
    {
        Matrix< DDRMat > tDConstraintDCriteria( 1, 2 );
        tDConstraintDCriteria( 0 ) = 0;
        tDConstraintDCriteria( 1 ) = 1.0;

        return tDConstraintDCriteria;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    HMRParameterList( Vector< Vector< Parameter_List > > &tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_hmr_parameter_list();

        tParameterlist( 0 )( 0 ).set( "number_of_elements_per_dimension", "2, 2" );
        tParameterlist( 0 )( 0 ).set( "domain_dimensions", "2.0, 2.0" );
        tParameterlist( 0 )( 0 ).set( "domain_offset", "-1.0, -1.0" );
        tParameterlist( 0 )( 0 ).set( "domain_sidesets", "1,2,3,4" );
        tParameterlist( 0 )( 0 ).set( "lagrange_output_meshes", "0" );

        tParameterlist( 0 )( 0 ).set( "lagrange_orders", "1" );
        tParameterlist( 0 )( 0 ).set( "lagrange_pattern", "0" );
        tParameterlist( 0 )( 0 ).set( "bspline_orders", "1" );
        tParameterlist( 0 )( 0 ).set( "bspline_pattern", "0" );

        tParameterlist( 0 )( 0 ).set( "lagrange_to_bspline", "0" );

        tParameterlist( 0 )( 0 ).set( "truncate_bsplines", 1 );
        tParameterlist( 0 )( 0 ).set( "refinement_buffer", 3 );
        tParameterlist( 0 )( 0 ).set( "staircase_buffer", 3 );
        tParameterlist( 0 )( 0 ).set( "initial_refinement", "0" );
        tParameterlist( 0 )( 0 ).set( "initial_refinement_pattern", "0" );

        tParameterlist( 0 )( 0 ).set( "use_multigrid", 0 );

        tParameterlist( 0 )( 0 ).set( "adaptive_refinement_level", 1 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    XTKParameterList( Vector< Vector< Parameter_List > > &tParameterlist )
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
        tParameterlist( 0 )( 0 ).set( "verbose", false );
        tParameterlist( 0 )( 0 ).set( "print_enriched_ig_mesh", false );
        tParameterlist( 0 )( 0 ).set( "exodus_output_XTK_ig_mesh", true );
        tParameterlist( 0 )( 0 ).set( "high_to_low_dbl_side_sets", true );
    }

    //--------------------------------------------------------------------------------------------------------------
    moris::real
    Oblique_Line_Func(
            const moris::Matrix< DDRMat > &aCoordinates,
            const Vector< real >          &aGeometryParameters )
    {
        // get coordinates
        real tX = aCoordinates( 0 );
        real tY = aCoordinates( 1 );

        real tLev = 1.0 * ( tX + 0.6 ) + 1.0 * ( tY + 0.3 );

        //        if ( std::sqrt( std::pow( tX + 0.5, 2 ) + std::pow( tY + 0.5, 2 ) ) < 0.2 )
        //        {
        tLev += aGeometryParameters( 0 );
        //       }

        return tLev;
    }

    void
    Oblique_Line_Deriv(
            const moris::Matrix< moris::DDRMat >           &aCoordinates,
            const Vector< moris::Matrix< moris::DDRMat > > &aParameters,
            moris::Matrix< DDRMat >                        &aFieldSensitivity )
    {
        // get coordinates
        //        real tX = aCoordinates( 0 );
        //        real tY = aCoordinates( 1 );

        //       if ( std::sqrt( std::pow( tX + 0.5, 2 ) + std::pow( tY + 0.5, 2 ) ) < 0.2 )
        //       {
        aFieldSensitivity = { { 1.0 } };
        //        }
        //        else
        //        {
        //            aFieldSensitivity = { { 0.0 } };
        //        }
    }

    void
    GENParameterList( Vector< Vector< Parameter_List > > &tParameterlist )
    {
        tParameterlist.resize( 3 );
        tParameterlist( 0 ).resize( 1 );
        tParameterlist( 0 )( 0 ) = moris::prm::create_gen_parameter_list();
        tParameterlist( 0 )( 0 ).set( "IQI_types", "IQIBulkStrainEnergy", "IQIBulkVolume" );

        // Geometry parameter lists
        tParameterlist( 1 ).resize( 2 );

        // vertical line
        tParameterlist( 1 )( 0 ) = prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE );
        tParameterlist( 1 )( 0 ).set( "center_x", 0.8 );
        tParameterlist( 1 )( 0 ).set( "center_y", 0.3 );
        tParameterlist( 1 )( 0 ).set( "normal_x", 1.0 );
        tParameterlist( 1 )( 0 ).set( "normal_y", 0.0 );

        // oblique line
        //        tParameterlist( 1 )( 1 ) = prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE );
        //        tParameterlist( 1 )( 1 ).set( "center_x", -0.6, -0.6, -0.6 );
        //        tParameterlist( 1 )( 1 ).set( "center_y", -0.3, -0.3, -0.3 );
        //        tParameterlist( 1 )( 1 ).set( "normal_x", 1.0, 1.0, 1.0 );
        //        tParameterlist( 1 )( 1 ).set( "normal_y", 1.0, 1.0, 1.0 );

        tParameterlist( 1 )( 1 ) = prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED );
        real deps                = 0.0e-5;
        tParameterlist( 1 )( 1 ).set( "field_function_name", "Oblique_Line_Func" );
        tParameterlist( 1 )( 1 ).set( "sensitivity_function_name", "Oblique_Line_Deriv" );
        tParameterlist( 1 )( 1 ).set( "use_multilinear_interpolation", true );
        tParameterlist( 1 )( 1 ).insert( "offset_diagonal", Design_Variable( 0.0 + deps, 0.0 + deps, 0.0 + deps ) );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    FEMParameterList( Vector< Vector< Parameter_List > > &tParameterList )
    {
        // create a cell of cell of parameter list for fem
        tParameterList.resize( 8 );

        //------------------------------------------------------------------------------
        // init property counter
        uint tPropCounter = 0;

        // create parameter list for property 1
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ) = prm::create_property_parameter_list();
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropDensity" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "1.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // create parameter list for property 2
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropYoungs" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "1.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // create parameter list for property 5
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropFlux" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "10.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // create parameter list for property 4
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropDirichletU" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "0.0;0.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // create parameter list for property 10
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropTraction" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "0.0;0.1" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // create parameter list for property 7
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropPoisson" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "0.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );

        //------------------------------------------------------------------------------
        // init CM counter
        uint tCMCounter = 0;

        // create parameter list for constitutive model 1
        tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_name", "CMStrucLinIso1" );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_type", fem::Constitutive_Type::STRUC_LIN_ISO );
        tParameterList( 1 )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
        tParameterList( 1 )( tCMCounter ).set( "properties", "PropYoungs,YoungsModulus;PropPoisson,PoissonRatio" );

        //------------------------------------------------------------------------------
        // init SP counter
        uint tSPCounter = 0;

        // create parameter list for stabilization parameter 1
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", "SPNitscheTemp" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", "100.0" );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", "PropYoungs,Material" );

        //------------------------------------------------------------------------------
        // init IWG counter
        uint tIWGCounter = 0;

        // create parameter list for IWG 1
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGBulkU_1" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_BULK );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "UX,UY" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tMeshSets );
        tIWGCounter++;

        // create parameter list for IWG 2
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGDirichletU" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "UX,UY" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties", "PropDirichletU,Dirichlet" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPNitscheTemp,DirichletNitsche" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tDBCSets );
        tIWGCounter++;

        // create parameter list for IWG 3
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGTraction" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_NEUMANN );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "UX,UY" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties", "PropTraction,Traction" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tNBCSets );

        //------------------------------------------------------------------------------
        // init IQI counter
        uint tIQICounter = 0;

        // create parameter list for IQI 4
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIDispX" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", fem::IQI_Type::DOF );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "UX,UY" );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 0 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tMeshSets );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIDispY" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", fem::IQI_Type::DOF );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "UX,UY" );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 1 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tMeshSets );
        tIQICounter++;

        // create parameter list for IQI 4
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQINodeX" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", fem::IQI_Type::NODE_SENSITIVITY );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 0 );
        tParameterList( 4 )( tIQICounter ).set( "function_parameters", "0" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tMeshSets );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQINodeY" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", fem::IQI_Type::NODE_SENSITIVITY );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 1 );
        tParameterList( 4 )( tIQICounter ).set( "function_parameters", "0" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tMeshSets );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIDispXDV" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", fem::IQI_Type::DOF_SENSITIVITY );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "UX,UY" );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 0 );
        tParameterList( 4 )( tIQICounter ).set( "function_parameters", "0" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tMeshSets );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIDispYDV" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", fem::IQI_Type::DOF_SENSITIVITY );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "UX,UY" );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 1 );
        tParameterList( 4 )( tIQICounter ).set( "function_parameters", "0" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tMeshSets );
        tIQICounter++;

        // create parameter list for IQI 4
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkStrainEnergy" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", fem::IQI_Type::STRAIN_ENERGY );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 4 )( tIQICounter ).set( "leader_constitutive_models", "CMStrucLinIso1,Elast" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tMeshSets );
        tIQICounter++;

        // create parameter list for IQI 4
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkVolume" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", fem::IQI_Type::VOLUME );
        tParameterList( 4 )( tIQICounter ).set( "leader_properties", "PropDensity,Density" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tMeshSets );

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        tParameterList( 5 ).resize( 1 );
        tParameterList( 5 )( 0 ) = prm::create_computation_parameter_list();

        tParameterList( 5 )( 0 ).set( "finite_difference_scheme", fem::FDScheme_Type::POINT_3_CENTRAL );
        tParameterList( 5 )( 0 ).set( "finite_difference_perturbation_size", tFEMFdEpsilon );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    SOLParameterList( Vector< Vector< Parameter_List > > &tParameterlist )
    {
        tParameterlist.resize( 8 );
        for ( uint Ik = 0; Ik < 8; Ik++ )
        {
            tParameterlist( Ik ).resize( 1 );
        }

        tParameterlist( 0 )( 0 ) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL );
        tParameterlist( 0 )( 0 ).set( "preconditioners", "0" );

        tParameterlist( 1 )( 0 ) = moris::prm::create_linear_solver_parameter_list();

        tParameterlist( 2 )( 0 ) = moris::prm::create_nonlinear_algorithm_parameter_list();

        tParameterlist( 3 )( 0 ) = moris::prm::create_nonlinear_solver_parameter_list();
        tParameterlist( 3 )( 0 ).set( "NLA_DofTypes", "UX,UY" );

        tParameterlist( 4 )( 0 ) = moris::prm::create_time_solver_algorithm_parameter_list();
        tParameterlist( 4 )( 0 ).set( "TSA_Num_Time_Steps", 1 );
        tParameterlist( 4 )( 0 ).set( "TSA_Time_Frame", 1.0 );

        tParameterlist( 5 )( 0 ) = moris::prm::create_time_solver_parameter_list();
        tParameterlist( 5 )( 0 ).set( "TSA_DofTypes", "UX,UY" );
        tParameterlist( 5 )( 0 ).set( "TSA_Output_Indices", "0" );
        tParameterlist( 5 )( 0 ).set( "TSA_Output_Criteria", "Output_Criterion" );

        tParameterlist( 6 )( 0 ) = moris::prm::create_solver_warehouse_parameterlist();
        tParameterlist( 6 )( 0 ).set( "Sensitivity_Analysis_Type", sol::SensitivityAnalysisType::DIRECT );
        tParameterlist( 6 )( 0 ).set( "SOL_save_operator_to_matlab", "operator" );

        tParameterlist( 7 )( 0 ) = moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::IFPACK );
        tParameterlist( 7 )( 0 ).set( "ifpack_prec_type", "ILU" );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    MSIParameterList( Vector< Vector< Parameter_List > > &tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_msi_parameter_list();
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    VISParameterList( Vector< Vector< Parameter_List > > &tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_vis_parameter_list();
        tParameterlist( 0 )( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", "shape_sensitivities.exo" ) );
        tParameterlist( 0 )( 0 ).set( "Mesh_Type", vis::VIS_Mesh_Type::STANDARD );
        tParameterlist( 0 )( 0 ).set( "Set_Names", tMeshSets );
        tParameterlist( 0 )( 0 ).set( "Field_Names", "UX,UY,DPX,DPY,DUX,DUY" );
        tParameterlist( 0 )( 0 ).set( "Field_Type", "NODAL,NODAL,NODAL,NODAL,NODAL,NODAL" );
        tParameterlist( 0 )( 0 ).set( "IQI_Names", "IQIDispX,IQIDispY,IQINodeX,IQINodeY,IQIDispXDV,IQIDispYDV" );
        tParameterlist( 0 )( 0 ).set( "Analysis_Type", "FORWARD,FORWARD,SENSITIVITY,SENSITIVITY,SENSITIVITY,SENSITIVITY" );
        tParameterlist( 0 )( 0 ).set( "Save_Frequency", 1 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    OPTParameterList( Vector< Vector< Parameter_List > > &tParameterlist )
    {
        tParameterlist.resize( 3 );
        tParameterlist( 0 ).resize( 1 );
        tParameterlist( 1 ).resize( 0 );
        tParameterlist( 2 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = moris::prm::create_opt_problem_parameter_list();
        tParameterlist( 0 )( 0 ).set( "is_optimization_problem", true );
        tParameterlist( 0 )( 0 ).set( "problem", "user_defined" );
        tParameterlist( 0 )( 0 ).set( "library", "Shape_Sensitivity_Sweep.so" );

        tParameterlist( 2 )( 0 ) = moris::prm::create_sweep_parameter_list();
        tParameterlist( 2 )( 0 ).set( "hdf5_path", "shape_opt_test.hdf5" );
        tParameterlist( 2 )( 0 ).set( "evaluate_objective_gradients", true );
        tParameterlist( 2 )( 0 ).set( "evaluate_constraint_gradients", true );
        tParameterlist( 2 )( 0 ).set( "num_evaluations_per_adv", "1" );
        tParameterlist( 2 )( 0 ).set( "include_bounds", false );
        tParameterlist( 2 )( 0 ).set( "finite_difference_type", "all" );
        tParameterlist( 2 )( 0 ).set( "finite_difference_epsilons", tFDsweep );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    MORISGENERALParameterList( Vector< Vector< Parameter_List > > &tParameterlist )
    {
    }

    //------------------------------------------------------------------------------
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
