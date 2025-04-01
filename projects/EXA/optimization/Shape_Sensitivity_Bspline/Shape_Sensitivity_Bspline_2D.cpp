/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Shape_Sensitivity_Bspline.cpp
 *
 */

#include "paths.hpp"
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

// Geometry model setup
//         vertical cut         |  oblique cut
// case 0: var: 0   + analytic  |  fixed                  dArea/ds = 0.45
// case 1: fixed                |  var: 0    + analytic   dArea/ds = 1.2
// case 2: var: 0   + analytic  |  var: 1    + analytic   dArea/ds = 0.45; 1.2
// case 3: var: 1   + analytic  |  var: 0    + analytic   dArea/ds = 1.2; 0.45
// case 4: var: 0-8 + b-spline  |  var: na
// case 5: na                   |  var: 0-8  + b-spline
// case 6: var: 0-8 + b-spline  |  var: 9-15 + b-spline
extern std::shared_ptr< moris::Vector< std::string > > tGeometrySetup;

#ifdef __cplusplus
extern "C" {
#endif
//------------------------------------------------------------------------------
namespace moris
{

    std::string tMeshSets = "HMR_dummy_n_p0,HMR_dummy_c_p0";

    // Analytic function for surface mesh
    Vector< real > SM_Perturbation( const Matrix< DDRMat >& aFacetVertexCoordinates,
            const Vector< real >&                           aADVs )
    {
        // ADV only affects x coordinates
        return { aADVs( 0 ), 0.0 };
    }

    //--------------------------------------------------------------------------------------------------------------

    void SM_Perturbation_Sensitivity( const Matrix< DDRMat >& aFacetVertexCoordinates,
            const Vector< real >&                             aADVs,
            Matrix< DDRMat >&                                 aSensitivity )
    {
        // Sensitivity of the coordinate perturbation by the ADV
        aSensitivity = { { 1.0 }, { 0.0 } };
    }

    //--------------------------------------------------------------------------------------------------------------

    // Discretization scaling for surface mesh, arbitrary just to test different scalings
    Vector< real >
    Facet_Vertex_Factor( const Matrix< DDRMat >& aCoordinates )
    {
        if ( aCoordinates( 0 ) > 1.0 )
        {
            return { 0.5, 1.0 };
        }
        else
        {
            return { 1.2, 1.2 };
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    // Constant function for properties
    void
    Func_Const( moris::Matrix< moris::DDRMat >&       aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    bool
    Output_Criterion( moris::tsa::Time_Solver* aTimeSolver )
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
    compute_objectives( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tObjectives = { { aCriteria( 0 ) } };

        return tObjectives;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_constraints( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tConstraints = { { aCriteria( 1 ) } };

        return tConstraints;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_dobjective_dadv( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tDObjectiveDADV( 1, aADVs.size(), 0.0 );

        return tDObjectiveDADV;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_dobjective_dcriteria( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tDObjectiveDCriteria( 1, 2 );
        tDObjectiveDCriteria( 0 ) = 1;
        tDObjectiveDCriteria( 1 ) = 0;

        return tDObjectiveDCriteria;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_dconstraint_dadv( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tDConstraintDADV( 1, aADVs.size(), 0.0 );

        return tDConstraintDADV;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_dconstraint_dcriteria( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tDConstraintDCriteria( 1, 2 );
        tDConstraintDCriteria( 0 ) = 0;
        tDConstraintDCriteria( 1 ) = 1.0;

        return tDConstraintDCriteria;
    }

    void
    Func_Traction_U(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        if ( aFIManager->get_IG_geometry_interpolator()->valx()( 1 ) < 0.2 )
        {
            aPropMatrix = { { 0.0 },
                { aParameters( 0 )( 0 ) } };
        }
        else
        {
            aPropMatrix = { { 0.0 },
                { 0.0 } };
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    HMRParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "number_of_elements_per_dimension", "2, 2" );
        aParameterLists.set( "domain_dimensions", "2.0, 2.0" );
        aParameterLists.set( "domain_offset", "-1.0, -1.0" );
        aParameterLists.set( "domain_sidesets", "1,2,3,4" );
        aParameterLists.set( "lagrange_output_meshes", "0" );

        aParameterLists.set( "lagrange_orders", "1" );
        aParameterLists.set( "lagrange_pattern", "0" );
        aParameterLists.set( "bspline_orders", "1" );
        aParameterLists.set( "bspline_pattern", "0" );

        aParameterLists.set( "lagrange_to_bspline", "0" );

        aParameterLists.set( "truncate_bsplines", 1 );
        aParameterLists.set( "refinement_buffer", 3 );
        aParameterLists.set( "staircase_buffer", 3 );
        aParameterLists.set( "initial_refinement", "0" );
        aParameterLists.set( "initial_refinement_pattern", "0" );

        aParameterLists.set( "use_multigrid", 0 );
        // aParameterLists.set( "severity_level", 1 );

        aParameterLists.set( "adaptive_refinement_level", 1 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    XTKParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "decompose", true );
        aParameterLists.set( "decomposition_type", "conformal" );
        aParameterLists.set( "enrich", true );
        aParameterLists.set( "basis_rank", "bspline" );
        aParameterLists.set( "enrich_mesh_indices", "0" );
        aParameterLists.set( "ghost_stab", true );
        aParameterLists.set( "multigrid", false );
        aParameterLists.set( "verbose", false );
        aParameterLists.set( "print_enriched_ig_mesh", false );
        aParameterLists.set( "exodus_output_XTK_ig_mesh", true );
        aParameterLists.set( "high_to_low_dbl_side_sets", true );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    GENParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "IQI_types", "IQIBulkStrainEnergy", "IQIBulkVolume" );
        aParameterLists.set( "output_mesh_file", "gen_shape_sensitivities.exo" );

        // Geometry parameter lists

        // Create geometries based on the keys in tGeometrySetup
        for ( const auto& tGeom : *tGeometrySetup )
        {
            // Get which geometry to create
            size_t      tNamePos = tGeom.find( "-" );
            std::string tName    = tGeom.substr( 0, tNamePos );

            // Get the ADV dependency type
            size_t      tADVPos = tGeom.find( "-", tNamePos + 1 );
            std::string tADV    = tGeom.substr( tNamePos + 1, tADVPos - tNamePos - 1 );

            // Get the delaunay type
            std::string tDelaunay = tGeom.substr( tADVPos + 1, tGeom.size() - tADVPos - 1 );

            if ( tName == "VL" )    // make vertical line geometry
            {
                aParameterLists( GEN::GEOMETRIES ).add_parameter_list( gen::Field_Type::LINE );
                aParameterLists.set( "center_y", -1.0 );
                aParameterLists.set( "normal_x", 1.0 );
                aParameterLists.set( "normal_y", 0.0 );

                // Set the ADV dependency
                if ( tADV == "A" )
                {
                    aParameterLists.set( "center_x", 0.0, 0.65, 1.0 );
                }
                else if ( tADV == "D" )
                {
                    aParameterLists.set( "center_x", 0.65 );
                    aParameterLists.set( "discretization_mesh_index", 0 );
                }
            }
            else if ( tName == "OL" )    // make oblique line geometry
            {
                aParameterLists( GEN::GEOMETRIES ).add_parameter_list( gen::Field_Type::LINE );
                aParameterLists.set( "center_y", -0.52 );
                aParameterLists.set( "normal_x", .707106781 );
                aParameterLists.set( "normal_y", .707106781 );

                // Set the ADV dependency
                if ( tADV == "A" )
                {
                    aParameterLists.set( "center_x", 0.0, 0.65, 1.0 );
                }
                else if ( tADV == "D" )
                {
                    aParameterLists.set( "center_x", 0.65 );
                    aParameterLists.set( "discretization_mesh_index", 0 );
                }
            }
            else if ( tName == "SM" )    // make surface mesh geometry
            {
                aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_surface_mesh_geometry_parameter_list() );
                aParameterLists.set( "file_path", moris::get_base_moris_dir() + "projects/GEN/test/data/test_trapezoid.obj" );

                // Set the ADV dependency
                if ( tADV == "A" )
                {
                    aParameterLists.insert( "x_pos", Design_Variable( -1.0, 0.0, 1.0 ) );
                    aParameterLists.set( "field_function_name", "SM_Perturbation" );
                    aParameterLists.set( "sensitivity_function_name", "SM_Perturbation_Sensitivity" );
                }
                else if ( tADV == "D" )
                {
                    aParameterLists.set( "discretization_mesh_index", 0 );
                    aParameterLists.set( "discretization_factor_function_name", "Facet_Vertex_Factor" );
                }
            }
            else
            {
                MORIS_ERROR( false, "Geometry model not implemented in test case, available keys are VL, OL, and SM - provided key is %s", tName.c_str() );
            }

            // Set delaunay if needed
            if ( tDelaunay == "D" )
            {
                aParameterLists.set( "delaunay", true );
            }
            else
            {
                aParameterLists.set( "delaunay", false );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    FEMParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.hack_for_legacy_fem();
        // create a cell of cell of parameter list for fem

        //------------------------------------------------------------------------------

        // create parameter list for property 1
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDensity" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 2
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropYoungs" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 5
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropFlux" );
        aParameterLists.set( "function_parameters", "10.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 4
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDirichletU" );
        aParameterLists.set( "function_parameters", "0.0;0.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 10
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropTraction" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Traction_U" );

        // create parameter list for property 7
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropPoisson" );
        aParameterLists.set( "function_parameters", "0.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        //------------------------------------------------------------------------------

        // create parameter list for constitutive model 1
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMStrucLinIso1" );
        aParameterLists.set( "constitutive_type", fem::Constitutive_Type::STRUC_LIN_ISO );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
        aParameterLists.set( "properties", "PropYoungs,YoungsModulus;PropPoisson,PoissonRatio" );

        //------------------------------------------------------------------------------

        // create parameter list for stabilization parameter 1
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPNitscheTemp" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "leader_properties", "PropYoungs,Material" );

        //------------------------------------------------------------------------------
        // create parameter list for IWG 1
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGBulkU_1" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_BULK );
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        aParameterLists.set( "mesh_set_names", tMeshSets );

        // create parameter list for IWG 2
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGDirichletU" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_properties", "PropDirichletU,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheTemp,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", "SideSet_4_n_p0,SideSet_4_c_p0,SideSet_4_c_p0" );

        // create parameter list for IWG 3
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGTraction" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_NEUMANN );
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_properties", "PropTraction,Traction" );
        aParameterLists.set( "mesh_set_names", "" );

        //------------------------------------------------------------------------------
        // create parameter list for IQI 4
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIDisp" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "mesh_set_names", tMeshSets );

        // create parameter list for IQI 4
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkStrainEnergy" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::STRAIN_ENERGY );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso1,Elast" );
        aParameterLists.set( "mesh_set_names", tMeshSets );

        // create parameter list for IQI 4
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkVolume" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::VOLUME );
        aParameterLists.set( "leader_properties", "PropDensity,Density" );
        aParameterLists.set( "mesh_set_names", tMeshSets );

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        aParameterLists( FEM::COMPUTATION );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    SOLParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::AMESOS_IMPL );

        aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list();

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "NLA_combined_res_jac_assembly", false );
        aParameterLists.set( "NLA_max_iter", 1 );

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "NLA_DofTypes", "UX,UY" );

        aParameterLists( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "TSA_Num_Time_Steps", 1 );
        aParameterLists.set( "TSA_Time_Frame", 1.0 );

        aParameterLists( SOL::TIME_SOLVERS ).add_parameter_list();
        aParameterLists.set( "TSA_DofTypes", "UX,UY" );
        aParameterLists.set( "TSA_Output_Indices", "0" );
        aParameterLists.set( "TSA_Output_Criteria", "Output_Criterion" );

        aParameterLists( SOL::PRECONDITIONERS ).add_parameter_list( sol::PreconditionerType::NONE );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "File_Name", std::pair< std::string, std::string >( "./", "shape_sensitivities.exo" ) );
        aParameterLists.set( "Mesh_Type", vis::VIS_Mesh_Type::STANDARD );
        aParameterLists.set( "Set_Names", tMeshSets );
        aParameterLists.set( "Field_Names", "U" );
        aParameterLists.set( "Field_Type", "NODAL" );
        aParameterLists.set( "IQI_Names", "IQIDisp" );
        aParameterLists.set( "Save_Frequency", 1 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    OPTParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "is_optimization_problem", true );
        aParameterLists.set( "problem", "user_defined" );
        aParameterLists.set( "library", "Shape_Sensitivity_Bspline_2D.so" );

        aParameterLists( OPT::ALGORITHMS ).add_parameter_list( opt::Optimization_Algorithm_Type::SWEEP );
        aParameterLists.set( "hdf5_path", "shape_opt_test_2D.hdf5" );
        aParameterLists.set( "evaluate_objective_gradients", true );
        aParameterLists.set( "evaluate_constraint_gradients", true );
        aParameterLists.set( "num_evaluations_per_adv", "1" );
        aParameterLists.set( "include_bounds", false );
        aParameterLists.set( "finite_difference_type", "all" );
    }

    //--------------------------------------------------------------------------------------------------------------

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
