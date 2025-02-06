/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Bedding_Sensitivity_Test.cpp
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
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosBlockGmresSolMgr.hpp"

// problem dimension: 2D or 3D
extern moris::uint gDim;

// test case index
extern moris::uint gTestCaseIndex;

#ifdef __cplusplus
extern "C" {
#endif

//------------------------------------------------------------------------------

namespace moris
{
    // forward declare the funtions for petsc and trilinos
    void create_trilinos_parameter_list( Module_Parameter_Lists& ), create_petsc_parameter_list( Module_Parameter_Lists& );
    //------------------------------------------------------------------------------
    // Main problem parameters

    std::string tName = "Bedding_Sensitivity_Test";

    bool tIs3D     = gDim == 3 ? true : false;
    bool tIsOpt    = true;
    bool tUseGhost = true;

    // hole radius
    real tOuterRadius = 0.3;
    real tInnerRadius = 0.1;

    // Nitsche interface penalty
    real tNitschePenalty = 10.0;

    // FD in adjoint
    real tFEMFdEpsilon = 1.0e-5;

    // FD step size in sweep
    std::string tFDsweep = "1.0e-5";

    // Number of constraints
    uint tNumConstraints = 3;

    // background mesh parameters
    Vector< uint > tNumElementsPerDir( 2 + tIs3D, 4 );
    std::string tDimensions        = tIs3D ? "1,1,1" : "1,1";
    std::string tOffSet            = tIs3D ? "0.0,0.0,0.0" : "0.0,0.0";

    int tRefineBuffer = 0;

    // note: pattern 0 - Levelset field  pattern 1 - displacement field
    std::string tLagrangeOrder   = "1";
    std::string tBsplineOrder    = "1";
    std::string tInitialRef      = "1";
    std::string tLagrangePattern = "0";

    int tDispOrder = 1;

    //------------------------------------------------------------------------------
    // Derived problem parameters

    std::string tProblemConfig = "_" + std::to_string( gTestCaseIndex );

    std::string tOutputFileName = tName + tProblemConfig + ".exo";
    std::string tLibraryName    = tName + ".so";
    std::string tGENOutputFile  = tName + tProblemConfig + "_GEN" + ".exo";
    std::string tHDF5FileName   = tName + tProblemConfig + "_SEN" + ".hdf5";

    std::string tMaterialSets = "HMR_dummy_n_p1,HMR_dummy_c_p1";

    std::string tTotalDomainSets = tMaterialSets;

    std::string tMaterialGhost = "ghost_p1";

    std::string tSupportSSets = "SideSet_4_n_p1,SideSet_4_c_p1";

    std::string tDofStrg = tIs3D ? "UX,UY,UZ" : "UX,UY";

    std::string tDirichletStr = tIs3D ? "0.0;0.0;0.0" : "0.0;0.0";

    std::string tBodyLoadStr = tIs3D ? "1;0.0;0.0" : "1;0.0";

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
        Matrix< DDSMat > tConstraintTypes( tNumConstraints, 1, 1 );

        return tConstraintTypes;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_objectives( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tObjectives( 1, 1 );

        real obj0 = aCriteria( 0 );
        real obj1 = aCriteria( 1 );
        real obj2 = aCriteria( 2 );

        tObjectives( 0, 0 ) = obj0 + obj1 + obj2;

        std::cout << "% --------------------------------- % \n";
        std::cout << "Objective                   = " << tObjectives( 0, 0 ) << " \n";
        std::cout << "Strain Energy               = " << aCriteria( 0 ) << " ( " << obj0 / tObjectives( 0, 0 ) << " )\n";
        std::cout << "Strain Energy with Bedding  = " << aCriteria( 1 ) << " ( " << obj1 / tObjectives( 0, 0 ) << " )\n";
        std::cout << "Volume                      = " << aCriteria( 2 ) << " ( " << obj2 / tObjectives( 0, 0 ) << " )\n";
        std::cout << " \n";

        return tObjectives;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_constraints( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tConstraints( 1, tNumConstraints );

        tConstraints( 0 ) = aCriteria( 0 );
        tConstraints( 1 ) = aCriteria( 1 );
        tConstraints( 2 ) = aCriteria( 2 );

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

        tDObjectiveDCriteria( 0 ) = 1.0;
        tDObjectiveDCriteria( 1 ) = 1.0;
        tDObjectiveDCriteria( 2 ) = 1.0;

        return tDObjectiveDCriteria;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_dconstraint_dadv(
            const Vector< real >& aADVs,
            const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tDConstraintDADV( tNumConstraints, aADVs.size(), 0.0 );

        return tDConstraintDADV;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_dconstraint_dcriteria(
            const Vector< real >& aADVs,
            const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tDConstraintDCriteria( tNumConstraints, aCriteria.size(), 0.0 );

        tDConstraintDCriteria( 0, 0 ) = 1.0;
        tDConstraintDCriteria( 1, 1 ) = 1.0;
        tDConstraintDCriteria( 2, 2 ) = 1.0;

        return tDConstraintDCriteria;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    OPTParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "is_optimization_problem", tIsOpt );
        aParameterLists.set( "problem", "user_defined" );
        aParameterLists.set( "library", tLibraryName );
        aParameterLists.set( "restart_file", "" );

        aParameterLists( OPT::ALGORITHMS ).add_parameter_list( opt::Optimization_Algorithm_Type::SWEEP );
        aParameterLists.set( "hdf5_path", tHDF5FileName );
        aParameterLists.set( "num_evaluations_per_adv", "1" );
        aParameterLists.set( "finite_difference_type", "all" );
        aParameterLists.set( "finite_difference_epsilons", tFDsweep );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    HMRParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "number_of_elements_per_dimension", tNumElementsPerDir );
        aParameterLists.set( "domain_dimensions", tDimensions );
        aParameterLists.set( "domain_offset", tOffSet );
        aParameterLists.set( "lagrange_output_meshes", "0" );

        aParameterLists.set( "lagrange_orders", tLagrangeOrder );
        aParameterLists.set( "lagrange_pattern", "0" );

        aParameterLists.set( "bspline_orders", tBsplineOrder );
        aParameterLists.set( "bspline_pattern", "0" );

        aParameterLists.set( "initial_refinement", tInitialRef );
        aParameterLists.set( "initial_refinement_pattern", "0" );
        aParameterLists.set( "use_advanced_T_matrix_scheme", 1 );

        aParameterLists.set( "lagrange_to_bspline", "0" );

        aParameterLists.set( "truncate_bsplines", 1 );
        aParameterLists.set( "refinement_buffer", tRefineBuffer );
        aParameterLists.set( "staircase_buffer", tRefineBuffer );

        aParameterLists.set( "use_number_aura", 1 );

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
        aParameterLists.set( "enrich_mesh_indices", "0" );
        aParameterLists.set( "multigrid", false );
        aParameterLists.set( "verbose", true );
        aParameterLists.set( "print_enriched_ig_mesh", false );

        aParameterLists.set( "ghost_stab", tUseGhost );

        aParameterLists.set( "exodus_output_XTK_ig_mesh", true );
        aParameterLists.set( "high_to_low_dbl_side_sets", true );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    GENParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists.set(    //
                "IQI_types",
                "IQIBulkStrainEnergy",
                "IQIBulkStrainEnergyWithBedding",
                "IQIBulkVolume" );

        aParameterLists.set( "output_mesh_file", tGENOutputFile );
        aParameterLists.set( "time_offset", 10.0 );

        Matrix< DDUMat > tPhaseMap( 4, 1, 0 );
        tPhaseMap( 0 ) = 1;
        tPhaseMap( 1 ) = 0;
        tPhaseMap( 2 ) = 1;
        tPhaseMap( 3 ) = 1;
        aParameterLists.set( "phase_table", moris::ios::stringify( tPhaseMap ) );
        aParameterLists.set( "print_phase_table", true );

        if ( tIs3D )
        {
            aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::SPHERE ) );
            aParameterLists.set( "center_z", 0.5 );
            aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::SPHERE ) );
            aParameterLists.set( "center_z", 0.5 );
        }
        else
        {
            aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::CIRCLE ) );
            aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::CIRCLE ) );
        }
        aParameterLists( GEN::GEOMETRIES )( 0 ).set( "center_x", 0.5 );
        aParameterLists( GEN::GEOMETRIES )( 0 ).set( "center_y", 0.5 );
        aParameterLists( GEN::GEOMETRIES )( 0 ).set( "radius", 0.9 * tOuterRadius, tOuterRadius, 1.1 * tOuterRadius );

        aParameterLists( GEN::GEOMETRIES )( 1 ).set( "center_x", 0.5 );
        aParameterLists( GEN::GEOMETRIES )( 1 ).set( "center_y", 0.5 );
        aParameterLists( GEN::GEOMETRIES )( 1 ).set( "radius", 0.9 * tInnerRadius, tInnerRadius, 1.1 * tInnerRadius );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    FEMParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.hack_for_legacy_fem();
        // create a cell of cell of parameter list for fem

        //------------------------------------------------------------------------------

        // create parameter list for density property
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDensity" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for Young's modulus property
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropYoungs" );
        aParameterLists.set( "function_parameters", "2.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for bedding property
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropBedding" );
        aParameterLists.set( "function_parameters", "2.0e-2" );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for DBC property
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDirichletU" );
        aParameterLists.set( "function_parameters", tDirichletStr );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for body load property
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "tPropBodyLoad" );
        aParameterLists.set( "function_parameters", tBodyLoadStr );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for Poisson ratio property
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropPoisson" );
        aParameterLists.set( "function_parameters", "0.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        //------------------------------------------------------------------------------

        // create parameter list for constitutive model
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMStrucLinIso_Material" );
        aParameterLists.set( "constitutive_type", fem::Constitutive_Type::STRUC_LIN_ISO );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( tDofStrg, "Displacement" ) );
        aParameterLists.set( "properties", "PropYoungs,YoungsModulus;PropPoisson,PoissonRatio" );

        //------------------------------------------------------------------------------

        // create parameter list for stabilization parameter 1
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPNitscheDirichletBC" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        aParameterLists.set( "function_parameters", std::to_string( tNitschePenalty ) );
        aParameterLists.set( "leader_properties", "PropYoungs,Material" );

        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", std::string( "SPGhost_Material" ) );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        aParameterLists.set( "function_parameters", std::string( "0.005" ) );
        aParameterLists.set( "leader_properties", std::string( "PropYoungs,Material" ) );

        //------------------------------------------------------------------------------
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGBulkMaterial" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_BULK );
        aParameterLists.set( "dof_residual", tDofStrg );
        aParameterLists.set( "leader_dof_dependencies", tDofStrg );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso_Material,ElastLinIso" );
        aParameterLists.set( "leader_properties", "PropBedding,Bedding;tPropBodyLoad,Load" );
        aParameterLists.set( "mesh_set_names", tMaterialSets );

        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGDirichletU" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", tDofStrg );
        aParameterLists.set( "leader_dof_dependencies", tDofStrg );
        aParameterLists.set( "leader_properties", "PropDirichletU,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso_Material,ElastLinIso" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheDirichletBC,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", tSupportSSets );

        if ( tUseGhost )
        {
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", std::string( "IWGGhostMaterial" ) );
            aParameterLists.set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists.set( "dof_residual", tDofStrg );
            aParameterLists.set( "leader_dof_dependencies", tDofStrg );
            aParameterLists.set( "follower_dof_dependencies", tDofStrg );
            aParameterLists.set( "stabilization_parameters", std::string( "SPGhost_Material,GhostSP" ) );
            aParameterLists.set( "ghost_order", (uint)tDispOrder );
            aParameterLists.set( "mesh_set_names", tMaterialGhost );
            }

        //------------------------------------------------------------------------------
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkUX" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", tDofStrg );
        aParameterLists.set( "leader_dof_dependencies", tDofStrg );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "mesh_set_names", tTotalDomainSets );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkUY" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "leader_dof_dependencies", tDofStrg );
        aParameterLists.set( "dof_quantity", tDofStrg );
        aParameterLists.set( "vectorial_field_index", 1 );
        aParameterLists.set( "mesh_set_names", tTotalDomainSets );

        if ( tIs3D )
        {
            aParameterLists( FEM::IQI ).add_parameter_list();
            aParameterLists.set( "IQI_name", "IQIBulkUZ" );
            aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
            aParameterLists.set( "leader_dof_dependencies", tDofStrg );
            aParameterLists.set( "dof_quantity", tDofStrg );
            aParameterLists.set( "vectorial_field_index", 2 );
            aParameterLists.set( "mesh_set_names", tTotalDomainSets );
            }

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkStrainEnergy" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::STRAIN_ENERGY );
        aParameterLists.set( "leader_dof_dependencies", tDofStrg );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso_Material,Elast" );
        aParameterLists.set( "mesh_set_names", tMaterialSets );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkStrainEnergyWithBedding" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::STRAIN_ENERGY );
        aParameterLists.set( "leader_dof_dependencies", tDofStrg );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso_Material,Elast" );
        aParameterLists.set( "leader_properties", "PropBedding,Bedding" );
        aParameterLists.set( "mesh_set_names", tMaterialSets );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkVolume" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::VOLUME );
        aParameterLists.set( "leader_properties", "PropDensity,Density" );
        aParameterLists.set( "mesh_set_names", tMaterialSets );

        // create computation  parameter list
        aParameterLists( FEM::COMPUTATION );
        aParameterLists.set( "print_physics_model", false );

        aParameterLists.set( "finite_difference_scheme", fem::FDScheme_Type::POINT_3_CENTRAL );
        aParameterLists.set( "finite_difference_perturbation_size", tFEMFdEpsilon );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    SOLParameterList( Module_Parameter_Lists& aParameterLists )
    {

        ( gTestCaseIndex == 21 ) ? create_petsc_parameter_list( aParameterLists ) : create_trilinos_parameter_list( aParameterLists );

        aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list();

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "NLA_combined_res_jac_assembly", true );
        aParameterLists.set( "NLA_rel_res_norm_drop", 1e-9 );
        aParameterLists.set( "NLA_relaxation_parameter", 1.00 );
        aParameterLists.set( "NLA_max_iter", 1 );

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "NLA_DofTypes", tDofStrg );

        aParameterLists( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list();

        aParameterLists( SOL::TIME_SOLVERS ).add_parameter_list();
        aParameterLists.set( "TSA_DofTypes", tDofStrg );
        aParameterLists.set( "TSA_Output_Indices", "0" );
        aParameterLists.set( "TSA_Output_Criteria", "Output_Criterion" );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "UX", 0 );
        aParameterLists.set( "UY", 0 );
        if ( tIs3D )
        {
            aParameterLists.set( "UZ", 0 );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName ) );
        aParameterLists.set( "Mesh_Type", vis::VIS_Mesh_Type::STANDARD );
        aParameterLists.set( "Set_Names", tTotalDomainSets );

        if ( tIs3D )
        {
            aParameterLists.set( "Field_Names", std::string( "UX,UY,UZ,StrainEnergy,StrainEnergyWithBedding,Volume" ) );
            aParameterLists.set( "Field_Type", std::string( "NODAL,NODAL,NODAL,GLOBAL,GLOBAL,GLOBAL" ) );
            aParameterLists.set( "IQI_Names", std::string( "IQIBulkUX,IQIBulkUY,IQIBulkUZ,IQIBulkStrainEnergy,"
                                                                    "IQIBulkStrainEnergyWithBedding,IQIBulkVolume" ) );
        }
        else
        {
            aParameterLists.set( "Field_Names", std::string( "UX,UY,StrainEnergy,StrainEnergyWithBedding,Volume" ) );
            aParameterLists.set( "Field_Type", std::string( "NODAL,NODAL,GLOBAL,GLOBAL,GLOBAL" ) );
            aParameterLists.set( "IQI_Names", std::string( "IQIBulkUX,IQIBulkUY,IQIBulkStrainEnergy,"
                                                                    "IQIBulkStrainEnergyWithBedding,IQIBulkVolume" ) );
        }

        aParameterLists.set( "Save_Frequency", 1 );
        aParameterLists.set( "Time_Offset", 10.0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    MORISGENERALParameterList( Module_Parameter_Lists& aParameterLists )
    {
    }

    //--------------------------------------------------------------------------------------------------------------
    void create_petsc_parameter_list( Module_Parameter_Lists & aParameterLists )
    {

        aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::PETSC );
        aParameterLists.set( "KSPType", "fgmres " );
        aParameterLists.set( "preconditioners", "0" );
        aParameterLists.set( "KSPTol", 1e-10 );

        aParameterLists( SOL::SOLVER_WAREHOUSE ).set( "SOL_TPL_Type", sol::MapType::Petsc );
        // aParameterLists.set( "SOL_save_operator_to_matlab", "jacc_par" );

        aParameterLists( SOL::PRECONDITIONERS ).add_parameter_list(  sol::PreconditionerType::PETSC );
        aParameterLists.set( "PCType", "mumps" );
    }

    //--------------------------------------------------------------------------------------------------------------

    void create_trilinos_parameter_list( Module_Parameter_Lists & aParameterLists )
    {
        aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::AMESOS_IMPL );

        aParameterLists( SOL::PRECONDITIONERS ).add_parameter_list(  sol::PreconditionerType::NONE );
    }

    //--------------------------------------------------------------------------------------------------------------
}    // namespace moris

//--------------------------------------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
