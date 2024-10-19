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
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
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
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
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
        aParameterLists.set( "is_optimization_problem", tIsOpt );
        aParameterLists.set( "problem", "user_defined" );
        aParameterLists.set( "library", tLibraryName );
        aParameterLists.set( "restart_file", "" );

        aParameterLists( OPT::ALGORITHMS ).add_parameter_list( opt::Optimization_Algorithm_Type::GCMMA );
        aParameterLists.set( "step_size", tMMAStepSize );
        aParameterLists.set( "penalty", tMMAPenalty );
        aParameterLists.set( "max_its", tMMAMaxIter );    // Maximum number of iterations
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    HMRParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "number_of_elements_per_dimension", "45,15" );
        aParameterLists.set( "domain_dimensions", "3,1" );
        aParameterLists.set( "domain_offset", "0.0,0.0" );
        aParameterLists.set( "domain_sidesets", "1,2,3,4" );
        aParameterLists.set( "lagrange_output_meshes", "0" );

        aParameterLists.set( "lagrange_orders", tLagrangeOrder );
        aParameterLists.set( "lagrange_pattern", "0" );
        aParameterLists.set( "bspline_orders", tBsplineOrder );
        aParameterLists.set( "bspline_pattern", "0,0" );

        aParameterLists.set( "lagrange_to_bspline", "0,1" );

        aParameterLists.set( "truncate_bsplines", 1 );
        aParameterLists.set( "refinement_buffer", 3 );
        aParameterLists.set( "staircase_buffer", 3 );
        aParameterLists.set( "initial_refinement", std::to_string( tInitialRef ) );

        aParameterLists.set( "use_multigrid", 0 );
        aParameterLists.set( "severity_level", 0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    XTKParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "decompose", true );
        aParameterLists.set( "decomposition_type", "conformal" );
        aParameterLists.set( "enrich", true );
        aParameterLists.set( "basis_rank", "bspline" );
        aParameterLists.set( "enrich_mesh_indices", "0,1" );
        aParameterLists.set( "ghost_stab", tUseGhost );
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

        aParameterLists.set( "IQI_types", "IQIBulkStrainEnergy_Frame", "IQIBulkStrainEnergy_Interior", "IQIBulkVolume_Interior", "IQIPerimeter_InterfaceVoid" );
        aParameterLists.set( "output_mesh_file", tGENOutputFile );

        Matrix< DDUMat > tPhaseMap( 4, 1, 0 );
        tPhaseMap( 1 ) = 1;
        tPhaseMap( 2 ) = 2;
        tPhaseMap( 3 ) = 2;
        aParameterLists.set( "phase_table", moris::ios::stringify( tPhaseMap ) );

        aParameterLists.set( "print_phase_table", true );

        // outer frame
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( gen::Field_Type::SUPERELLIPSE );
        aParameterLists.set( "center_x", 1.5 );
        aParameterLists.set( "center_y", 0.5 );
        aParameterLists.set( "semidiameter_x", 1.5 - tWallThickness );
        aParameterLists.set( "semidiameter_y", 0.5 - tWallThickness );
        aParameterLists.set( "exponent", 24.0 );

        // Geometry parameter lists
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_field_array_parameter_list( gen::Field_Type::SUPERELLIPSE ) );
        aParameterLists.set( "semidiameter_x", tHoleDiameter );
        aParameterLists.set( "semidiameter_y", tHoleDiameter );
        aParameterLists.set( "exponent", 8.0 );
        aParameterLists.set( "lower_bound_x", 0.12 );
        aParameterLists.set( "upper_bound_x", 2.88 );
        aParameterLists.set( "lower_bound_y", 0.12 );
        aParameterLists.set( "upper_bound_y", 0.88 );
        aParameterLists.set( "number_of_fields_x", tNumHoleX );
        aParameterLists.set( "number_of_fields_y", tNumHoleY );

        if ( tIsOpt )
        {
            aParameterLists.set( "discretization_mesh_index", 0 );
            aParameterLists.set( "discretization_lower_bound", -tBsplineLimit );
            aParameterLists.set( "discretization_upper_bound", tBsplineLimit );
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

        // create parameter list for property 2
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropBedding" );
        aParameterLists.set( "function_parameters", "1.0e-6" );
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
        aParameterLists.set( "value_function", "Func_Neumann_U" );

        // create parameter list for property 7
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropPoisson" );
        aParameterLists.set( "function_parameters", "0.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        //------------------------------------------------------------------------------

        // create parameter list for constitutive model 1
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMStrucLinIso_Frame" );
        aParameterLists.set( "constitutive_type", fem::Constitutive_Type::STRUC_LIN_ISO );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
        aParameterLists.set( "properties", "PropYoungs,YoungsModulus;PropPoisson,PoissonRatio" );

        // create parameter list for constitutive model 1
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMStrucLinIso_Interior" );
        aParameterLists.set( "constitutive_type", fem::Constitutive_Type::STRUC_LIN_ISO );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
        aParameterLists.set( "properties", "PropYoungs,YoungsModulus;PropPoisson,PoissonRatio" );

        //------------------------------------------------------------------------------

        // create parameter list for stabilization parameter 1
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPNitscheDirichletBC" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "leader_properties", "PropYoungs,Material" );

        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPNitscheFrameInteriorInterface" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::NITSCHE_INTERFACE );
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "leader_properties", "PropYoungs,Material" );
        aParameterLists.set( "follower_properties", "PropYoungs,Material" );

        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPGhost_Frame" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        aParameterLists.set( "function_parameters", "0.01" );
        aParameterLists.set( "leader_properties", "PropYoungs,Material" );

        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPGhost_Interior" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        aParameterLists.set( "function_parameters", "0.01" );
        aParameterLists.set( "leader_properties", "PropYoungs,Material" );

        //------------------------------------------------------------------------------
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGBulkU_Frame" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_BULK );
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso_Frame,ElastLinIso" );
        aParameterLists.set( "leader_properties", "PropBedding,Bedding" );
        aParameterLists.set( "mesh_set_names", tFrameSets );

        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGBulkU_Frame" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_BULK );
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso_Interior,ElastLinIso" );
        aParameterLists.set( "leader_properties", "PropBedding,Bedding" );
        aParameterLists.set( "mesh_set_names", tInteriorSets );

        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGDirichletU" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_properties", "PropDirichletU,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso_Frame,ElastLinIso" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheDirichletBC,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", tFrameSupportSSets );

        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGTraction" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_NEUMANN );
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_properties", "PropTraction,Traction" );
        aParameterLists.set( "mesh_set_names", tFrameLoadSSets );

        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGFrameInteriorInterface" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_INTERFACE_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "follower_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso_Frame,ElastLinIso" );
        aParameterLists.set( "follower_constitutive_models", "CMStrucLinIso_Interior,ElastLinIso" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheFrameInteriorInterface,NitscheInterface" );
        aParameterLists.set( "mesh_set_names", tFrameInteriorDSets );

        if ( tUseGhost )
        {
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "IWGGhostFrame" );
            aParameterLists.set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists.set( "dof_residual", "UX,UY" );
            aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
            aParameterLists.set( "follower_dof_dependencies", "UX,UY" );
            aParameterLists.set( "stabilization_parameters", "SPGhost_Frame,GhostSP" );
            aParameterLists.set( "mesh_set_names", tFrameGhost );

            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "IWGGhostInterior" );
            aParameterLists.set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists.set( "dof_residual", "UX,UY" );
            aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
            aParameterLists.set( "follower_dof_dependencies", "UX,UY" );
            aParameterLists.set( "stabilization_parameters", "SPGhost_Interior,GhostSP" );
            aParameterLists.set( "mesh_set_names", tInteriorGhost );
        }
        //------------------------------------------------------------------------------
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkUX" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "mesh_set_names", tTotalDomainSets );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkUY" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "dof_quantity", "UX,UY" );
        aParameterLists.set( "vectorial_field_index", 1 );
        aParameterLists.set( "mesh_set_names", tTotalDomainSets );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkStrainEnergy_Frame" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::STRAIN_ENERGY );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso_Frame,Elast" );
        aParameterLists.set( "mesh_set_names", tFrameSets );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkStrainEnergy_Interior" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::STRAIN_ENERGY );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso_Interior,Elast" );
        aParameterLists.set( "mesh_set_names", tInteriorSets );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkVolume_Frame" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::VOLUME );
        aParameterLists.set( "leader_properties", "PropDensity,Density" );
        aParameterLists.set( "mesh_set_names", tFrameSets );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkVolume_Interior" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::VOLUME );
        aParameterLists.set( "leader_properties", "PropDensity,Density" );
        aParameterLists.set( "mesh_set_names", tInteriorSets );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIPerimeter_InterfaceVoid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::VOLUME );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "mesh_set_names", tInterfaceVoidSSets );

        // create computation  parameter list
        aParameterLists( FEM::COMPUTATION );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    SOLParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::AMESOS_IMPL );

        aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list();

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "NLA_combined_res_jac_assembly", true );
        aParameterLists.set( "NLA_rel_res_norm_drop", 1.00 );
        aParameterLists.set( "NLA_relaxation_parameter", 1.00 );
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
        aParameterLists.set( "UX", 1 );
        aParameterLists.set( "UY", 1 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName ) );
        aParameterLists.set( "Mesh_Type", vis::VIS_Mesh_Type::STANDARD );
        aParameterLists.set( "Set_Names", tTotalDomainSets );
        aParameterLists.set( "Field_Names", "UX,UY" );
        aParameterLists.set( "Field_Type", "NODAL,NODAL" );
        aParameterLists.set( "IQI_Names", "IQIBulkUX,IQIBulkUY" );
        aParameterLists.set( "Save_Frequency", 1 );
        aParameterLists.set( "Time_Offset", 10.0 );
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
