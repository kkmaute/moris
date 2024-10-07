/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Levelset_Boxbeam.cpp
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
#include "cl_HMR_Element.hpp"
#include "fn_equal_to.hpp"
#include "fn_stringify_matrix.hpp"

#include "AztecOO.h"

#ifdef __cplusplus
extern "C" {
#endif

//------------------------------------------------------------------------------

namespace moris
{
    //------------------------------------------------------------------------------
    // Main problem parameters

    real tInitialStrainEnergy = 3.26439;
    real tInitialPerimeter    = 29.0669;
    real tPerimeterPenalty    = 0.2;

    real tMaxMass = 1.0;

    real tMMAPenalty  = 5.0;
    real tMMAStepSize = 0.01;
    int  tMMAMaxIter  = 2;

    int         tInitialRef    = 2;
    std::string tLagrangeOrder = "1";
    std::string tBsplineOrder  = "1,1";

    real tElementEdgeLength = 1.0 / 15.0 / pow( 2, tInitialRef );
    real tLoadLimitY        = std::floor( 0.2 / tElementEdgeLength ) * tElementEdgeLength;

    int  tNumHoleX     = 15;
    int  tNumHoleY     = 5;
    real tHoleDiameter = 0.065;

    real tWallThickness = 0.05;

    real tBsplineLimit = 0.06;

    bool tIsOpt    = true;
    bool tUseGhost = true;

    std::string tName = "Levelset_Boxbeam";

    //------------------------------------------------------------------------------
    // Derived problem paramters

    std::string tOutputFileName = tName + ".exo";
    std::string tLibraryName    = tName + ".so";
    std::string tGENOutputFile  = "GEN_" + tName + ".exo";

    std::string tFrameSets    = "HMR_dummy_n_p2,HMR_dummy_c_p2";
    std::string tInteriorSets = "HMR_dummy_n_p1,HMR_dummy_c_p1";

    std::string tTotalDomainSets = tFrameSets + "," + tInteriorSets;

    std::string tFrameGhost    = "ghost_p2";
    std::string tInteriorGhost = "ghost_p1";

    std::string tFrameLoadSSets    = "SideSet_2_n_p2,SideSet_2_c_p2";
    std::string tFrameSupportSSets = "SideSet_4_n_p2,SideSet_4_c_p2";
    std::string tFrameFreeSSets =
            "SideSet_1_n_p2,SideSet_1_c_p2"
            "SideSet_3_n_p2,SideSet_3_c_p2";

    std::string tInterfaceVoidSSets = "iside_b0_2_b1_0,iside_b0_1_b1_0";

    std::string tFrameInteriorDSets = "dbl_iside_p0_1_p1_2";

    //------------------------------------------------------------------------------
    // Constant function for properties
    void
    Func_Const(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    //------------------------------------------------------------------------------

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
        Matrix< DDRMat > tObjectives( 1, 1 );

        real obj1 = aCriteria( 0 ) / tInitialStrainEnergy;
        real obj2 = aCriteria( 1 ) / tInitialStrainEnergy;
        real obj3 = tPerimeterPenalty * aCriteria( 3 ) / tInitialPerimeter;

        tObjectives( 0, 0 ) = obj1 + obj2 + obj3;

        std::cout << "% --------------------------------- % \n";
        std::cout << "Objective                = " << tObjectives( 0, 0 ) << " \n";
        std::cout << "Strain Energy (Frame)    = " << aCriteria( 0 ) << " ( " << obj1 / tObjectives( 0, 0 ) << " )\n";
        std::cout << "Strain Energy (Interior) = " << aCriteria( 1 ) << " ( " << obj2 / tObjectives( 0, 0 ) << " )\n";
        std::cout << "Perimeter                = " << aCriteria( 3 ) << " ( " << obj3 / tObjectives( 0, 0 ) << " )\n";
        std::cout << " \n";

        std::cout << "min ADV                  = " << aADVs.min() << " \n";
        std::cout << "max ADV                  = " << aADVs.max() << " \n"
                  << std::flush;

        return tObjectives;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_constraints( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tConstraints( 1, 1 );
        tConstraints( 0 ) = aCriteria( 2 ) / tMaxMass - 1.0;

        std::cout << "Volume     = " << aCriteria( 2 ) << " \n";
        std::cout << "Constraint = " << tConstraints( 0 ) << " \n";
        std::cout << "% --------------------------------- % \n"
                  << std::flush;

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
    compute_dobjective_dcriteria(
            const Vector< real >& aADVs,
            const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tDObjectiveDCriteria( 1, aCriteria.size(), 0.0 );

        tDObjectiveDCriteria( 0 ) = 1.0 / tInitialStrainEnergy;
        tDObjectiveDCriteria( 1 ) = 1.0 / tInitialStrainEnergy;
        tDObjectiveDCriteria( 3 ) = tPerimeterPenalty / tInitialPerimeter;

        return tDObjectiveDCriteria;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_dconstraint_dadv(
            const Vector< real >& aADVs,
            const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tDConstraintDADV( 1, aADVs.size(), 0.0 );

        return tDConstraintDADV;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_dconstraint_dcriteria(
            const Vector< real >& aADVs,
            const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tDConstraintDCriteria( 1, aCriteria.size(), 0.0 );

        tDConstraintDCriteria( 2 ) = 1.0 / tMaxMass;

        return tDConstraintDCriteria;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Func_Neumann_U(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        if ( aFIManager->get_IG_geometry_interpolator()->valx()( 1 ) < tLoadLimitY )
        {
            aPropMatrix = { { 0.0 }, { aParameters( 0 )( 0 ) } };
        }
        else
        {
            aPropMatrix = { { 0.0 }, { 0.0 } };
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    OPTParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( moris::prm::create_opt_problem_parameter_list() );
        aParameterLists( 0 ).set( "is_optimization_problem", tIsOpt );
        aParameterLists( 0 ).set( "problem", "user_defined" );
        aParameterLists( 0 ).set( "library", tLibraryName );
        aParameterLists( 0 ).set( "restart_file", "" );

        aParameterLists( 2 ).add_parameter_list( moris::prm::create_gcmma_parameter_list() );
        aParameterLists( 2 ).set( "step_size", tMMAStepSize );
        aParameterLists( 2 ).set( "penalty", tMMAPenalty );
        aParameterLists( 2 ).set( "max_its", tMMAMaxIter );    // Maximum number of iterations
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    HMRParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_hmr_parameter_list() );

        aParameterLists( 0 ).set( "number_of_elements_per_dimension", "45,15" );
        aParameterLists( 0 ).set( "domain_dimensions", "3,1" );
        aParameterLists( 0 ).set( "domain_offset", "0.0,0.0" );
        aParameterLists( 0 ).set( "domain_sidesets", "1,2,3,4" );
        aParameterLists( 0 ).set( "lagrange_output_meshes", "0" );

        aParameterLists( 0 ).set( "lagrange_orders", tLagrangeOrder );
        aParameterLists( 0 ).set( "lagrange_pattern", "0" );
        aParameterLists( 0 ).set( "bspline_orders", tBsplineOrder );
        aParameterLists( 0 ).set( "bspline_pattern", "0,0" );

        aParameterLists( 0 ).set( "lagrange_to_bspline", "0,1" );

        aParameterLists( 0 ).set( "truncate_bsplines", 1 );
        aParameterLists( 0 ).set( "refinement_buffer", 3 );
        aParameterLists( 0 ).set( "staircase_buffer", 3 );
        aParameterLists( 0 ).set( "initial_refinement", std::to_string( tInitialRef ) );

        aParameterLists( 0 ).set( "use_multigrid", 0 );
        aParameterLists( 0 ).set( "severity_level", 0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    XTKParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_xtk_parameter_list() );
        aParameterLists( 0 ).set( "decompose", true );
        aParameterLists( 0 ).set( "decomposition_type", "conformal" );
        aParameterLists( 0 ).set( "enrich", true );
        aParameterLists( 0 ).set( "basis_rank", "bspline" );
        aParameterLists( 0 ).set( "enrich_mesh_indices", "0,1" );
        aParameterLists( 0 ).set( "ghost_stab", tUseGhost );
        aParameterLists( 0 ).set( "multigrid", false );
        aParameterLists( 0 ).set( "verbose", false );
        aParameterLists( 0 ).set( "print_enriched_ig_mesh", false );
        aParameterLists( 0 ).set( "exodus_output_XTK_ig_mesh", true );
        aParameterLists( 0 ).set( "high_to_low_dbl_side_sets", true );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    GENParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists( 0 ).add_parameter_list( moris::prm::create_gen_parameter_list() );
        aParameterLists( 0 ).set( "IQI_types", "IQIBulkStrainEnergy_Frame", "IQIBulkStrainEnergy_Interior", "IQIBulkVolume_Interior", "IQIPerimeter_InterfaceVoid" );
        aParameterLists( 0 ).set( "output_mesh_file", tGENOutputFile );

        Matrix< DDUMat > tPhaseMap( 4, 1, 0 );
        tPhaseMap( 1 ) = 1;
        tPhaseMap( 2 ) = 2;
        tPhaseMap( 3 ) = 2;
        aParameterLists( 0 ).set( "phase_table", moris::ios::stringify( tPhaseMap ) );

        aParameterLists( 0 ).set( "print_phase_table", true );

        // outer frame
        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::SUPERELLIPSE ) );
        aParameterLists( 1 ).set( "center_x", 1.5 );
        aParameterLists( 1 ).set( "center_y", 0.5 );
        aParameterLists( 1 ).set( "semidiameter_x", 1.5 - tWallThickness );
        aParameterLists( 1 ).set( "semidiameter_y", 0.5 - tWallThickness );
        aParameterLists( 1 ).set( "exponent", 24.0 );

        // Geometry parameter lists
        aParameterLists( 1 ).add_parameter_list( prm::create_field_array_parameter_list( gen::Field_Type::SUPERELLIPSE ) );
        aParameterLists( 1 ).set( "semidiameter_x", tHoleDiameter );
        aParameterLists( 1 ).set( "semidiameter_y", tHoleDiameter );
        aParameterLists( 1 ).set( "exponent", 8.0 );
        aParameterLists( 1 ).set( "lower_bound_x", 0.12 );
        aParameterLists( 1 ).set( "upper_bound_x", 2.88 );
        aParameterLists( 1 ).set( "lower_bound_y", 0.12 );
        aParameterLists( 1 ).set( "upper_bound_y", 0.88 );
        aParameterLists( 1 ).set( "number_of_fields_x", tNumHoleX );
        aParameterLists( 1 ).set( "number_of_fields_y", tNumHoleY );

        if ( tIsOpt )
        {
            aParameterLists( 1 ).set( "discretization_mesh_index", 0 );
            aParameterLists( 1 ).set( "discretization_lower_bound", -tBsplineLimit );
            aParameterLists( 1 ).set( "discretization_upper_bound", tBsplineLimit );
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
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropDensity" );
        aParameterLists( 0 ).set( "function_parameters", "1.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // create parameter list for property 2
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropYoungs" );
        aParameterLists( 0 ).set( "function_parameters", "1.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // create parameter list for property 2
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropBedding" );
        aParameterLists( 0 ).set( "function_parameters", "1.0e-6" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // create parameter list for property 5
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropFlux" );
        aParameterLists( 0 ).set( "function_parameters", "10.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // create parameter list for property 4
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropDirichletU" );
        aParameterLists( 0 ).set( "function_parameters", "0.0;0.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // create parameter list for property 10
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropTraction" );
        aParameterLists( 0 ).set( "function_parameters", "1.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Neumann_U" );

        // create parameter list for property 7
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropPoisson" );
        aParameterLists( 0 ).set( "function_parameters", "0.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        //------------------------------------------------------------------------------

        // create parameter list for constitutive model 1
        aParameterLists( 1 ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists( 1 ).set( "constitutive_name", "CMStrucLinIso_Frame" );
        aParameterLists( 1 ).set( "constitutive_type",  fem::Constitutive_Type::STRUC_LIN_ISO ) ;
        aParameterLists( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
        aParameterLists( 1 ).set( "properties", "PropYoungs,YoungsModulus;PropPoisson,PoissonRatio" );

        // create parameter list for constitutive model 1
        aParameterLists( 1 ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists( 1 ).set( "constitutive_name", "CMStrucLinIso_Interior" );
        aParameterLists( 1 ).set( "constitutive_type",  fem::Constitutive_Type::STRUC_LIN_ISO ) ;
        aParameterLists( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
        aParameterLists( 1 ).set( "properties", "PropYoungs,YoungsModulus;PropPoisson,PoissonRatio" );

        //------------------------------------------------------------------------------

        // create parameter list for stabilization parameter 1
        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPNitscheDirichletBC" );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::DIRICHLET_NITSCHE ) ;
        aParameterLists( 2 ).set( "function_parameters", "100.0" );
        aParameterLists( 2 ).set( "leader_properties", "PropYoungs,Material" );

        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPNitscheFrameInteriorInterface" );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::NITSCHE_INTERFACE ) ;
        aParameterLists( 2 ).set( "function_parameters", "100.0" );
        aParameterLists( 2 ).set( "leader_properties", "PropYoungs,Material" );
        aParameterLists( 2 ).set( "follower_properties", "PropYoungs,Material" );

        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPGhost_Frame" );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::GHOST_DISPL ) ;
        aParameterLists( 2 ).set( "function_parameters", "0.01" );
        aParameterLists( 2 ).set( "leader_properties", "PropYoungs,Material" );

        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPGhost_Interior" );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::GHOST_DISPL ) ;
        aParameterLists( 2 ).set( "function_parameters", "0.01" );
        aParameterLists( 2 ).set( "leader_properties", "PropYoungs,Material" );

        //------------------------------------------------------------------------------
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGBulkU_Frame" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_BULK ) ;
        aParameterLists( 3 ).set( "dof_residual", "UX,UY" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMStrucLinIso_Frame,ElastLinIso" );
        aParameterLists( 3 ).set( "leader_properties", "PropBedding,Bedding" );
        aParameterLists( 3 ).set( "mesh_set_names", tFrameSets );

        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGBulkU_Frame" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_BULK ) ;
        aParameterLists( 3 ).set( "dof_residual", "UX,UY" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMStrucLinIso_Interior,ElastLinIso" );
        aParameterLists( 3 ).set( "leader_properties", "PropBedding,Bedding" );
        aParameterLists( 3 ).set( "mesh_set_names", tInteriorSets );

        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGDirichletU" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        aParameterLists( 3 ).set( "dof_residual", "UX,UY" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists( 3 ).set( "leader_properties", "PropDirichletU,Dirichlet" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMStrucLinIso_Frame,ElastLinIso" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPNitscheDirichletBC,DirichletNitsche" );
        aParameterLists( 3 ).set( "mesh_set_names", tFrameSupportSSets );

        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGTraction" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_NEUMANN ) ;
        aParameterLists( 3 ).set( "dof_residual", "UX,UY" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists( 3 ).set( "leader_properties", "PropTraction,Traction" );
        aParameterLists( 3 ).set( "mesh_set_names", tFrameLoadSSets );

        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGFrameInteriorInterface" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_INTERFACE_UNSYMMETRIC_NITSCHE ) ;
        aParameterLists( 3 ).set( "dof_residual", "UX,UY" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists( 3 ).set( "follower_dof_dependencies", "UX,UY" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMStrucLinIso_Frame,ElastLinIso" );
        aParameterLists( 3 ).set( "follower_constitutive_models", "CMStrucLinIso_Interior,ElastLinIso" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPNitscheFrameInteriorInterface,NitscheInterface" );
        aParameterLists( 3 ).set( "mesh_set_names", tFrameInteriorDSets );

        if ( tUseGhost )
        {
            aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists( 3 ).set( "IWG_name", "IWGGhostFrame" );
            aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
            aParameterLists( 3 ).set( "dof_residual", "UX,UY" );
            aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY" );
            aParameterLists( 3 ).set( "follower_dof_dependencies", "UX,UY" );
            aParameterLists( 3 ).set( "stabilization_parameters", "SPGhost_Frame,GhostSP" );
            aParameterLists( 3 ).set( "mesh_set_names", tFrameGhost );

            aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists( 3 ).set( "IWG_name", "IWGGhostInterior" );
            aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
            aParameterLists( 3 ).set( "dof_residual", "UX,UY" );
            aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY" );
            aParameterLists( 3 ).set( "follower_dof_dependencies", "UX,UY" );
            aParameterLists( 3 ).set( "stabilization_parameters", "SPGhost_Interior,GhostSP" );
            aParameterLists( 3 ).set( "mesh_set_names", tInteriorGhost );
        }
        //------------------------------------------------------------------------------
        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkUX" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists( 4 ).set( "dof_quantity", "UX,UY" );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists( 4 ).set( "vectorial_field_index", 0 );
        aParameterLists( 4 ).set( "mesh_set_names", tTotalDomainSets );

        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkUY" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists( 4 ).set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists( 4 ).set( "dof_quantity", "UX,UY" );
        aParameterLists( 4 ).set( "vectorial_field_index", 1 );
        aParameterLists( 4 ).set( "mesh_set_names", tTotalDomainSets );

        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkStrainEnergy_Frame" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::STRAIN_ENERGY ) ;
        aParameterLists( 4 ).set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists( 4 ).set( "leader_constitutive_models", "CMStrucLinIso_Frame,Elast" );
        aParameterLists( 4 ).set( "mesh_set_names", tFrameSets );

        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkStrainEnergy_Interior" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::STRAIN_ENERGY ) ;
        aParameterLists( 4 ).set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists( 4 ).set( "leader_constitutive_models", "CMStrucLinIso_Interior,Elast" );
        aParameterLists( 4 ).set( "mesh_set_names", tInteriorSets );

        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkVolume_Frame" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::VOLUME ) ;
        aParameterLists( 4 ).set( "leader_properties", "PropDensity,Density" );
        aParameterLists( 4 ).set( "mesh_set_names", tFrameSets );

        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkVolume_Interior" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::VOLUME ) ;
        aParameterLists( 4 ).set( "leader_properties", "PropDensity,Density" );
        aParameterLists( 4 ).set( "mesh_set_names", tInteriorSets );

        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIPerimeter_InterfaceVoid" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::VOLUME ) ;
        aParameterLists( 4 ).set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists( 4 ).set( "mesh_set_names", tInterfaceVoidSSets );

        // create computation  parameter list
        aParameterLists( 5 ).add_parameter_list( prm::create_computation_parameter_list() );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    SOLParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL ) );

        aParameterLists( 1 ).add_parameter_list( moris::prm::create_linear_solver_parameter_list() );

        aParameterLists( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
        aParameterLists( 2 ).set( "NLA_combined_res_jac_assembly", true );
        aParameterLists( 2 ).set( "NLA_rel_res_norm_drop", 1.00 );
        aParameterLists( 2 ).set( "NLA_relaxation_parameter", 1.00 );
        aParameterLists( 2 ).set( "NLA_max_iter", 1 );

        aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
        aParameterLists( 3 ).set( "NLA_DofTypes", "UX,UY" );

        aParameterLists( 4 ).add_parameter_list( moris::prm::create_time_solver_algorithm_parameter_list() );
        aParameterLists( 4 ).set( "TSA_Num_Time_Steps", 1 );
        aParameterLists( 4 ).set( "TSA_Time_Frame", 1.0 );

        aParameterLists( 5 ).add_parameter_list( moris::prm::create_time_solver_parameter_list() );
        aParameterLists( 5 ).set( "TSA_DofTypes", "UX,UY" );
        aParameterLists( 5 ).set( "TSA_Output_Indices", "0" );
        aParameterLists( 5 ).set( "TSA_Output_Criteria", "Output_Criterion" );

        aParameterLists( 6 ).add_parameter_list( moris::prm::create_solver_warehouse_parameterlist() );

        aParameterLists( 7 ).add_parameter_list( moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::NONE ) );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_msi_parameter_list() );
        aParameterLists( 0 ).set( "UX", 1 );
        aParameterLists( 0 ).set( "UY", 1 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_vis_parameter_list() );
        aParameterLists( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName ) );
        aParameterLists( 0 ).set( "Mesh_Type",  vis::VIS_Mesh_Type::STANDARD ) ;
        aParameterLists( 0 ).set( "Set_Names", tTotalDomainSets );
        aParameterLists( 0 ).set( "Field_Names", "UX,UY" );
        aParameterLists( 0 ).set( "Field_Type", "NODAL,NODAL" );
        aParameterLists( 0 ).set( "IQI_Names", "IQIBulkUX,IQIBulkUY" );
        aParameterLists( 0 ).set( "Save_Frequency", 1 );
        aParameterLists( 0 ).set( "Time_Offset", 10.0 );
    }

    void
    MORISGENERALParameterList( Module_Parameter_Lists& aParameterLists )
    {
    }

    //--------------------------------------------------------------------------------------------------------------
}    // namespace moris

//--------------------------------------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
