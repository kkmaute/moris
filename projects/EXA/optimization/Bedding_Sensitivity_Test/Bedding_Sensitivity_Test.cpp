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
#include "fn_PRM_FEM_Parameters.hpp"
#include "fn_PRM_MSI_Parameters.hpp"
#include "fn_PRM_SOL_Parameters.hpp"
#include "fn_PRM_VIS_Parameters.hpp"
#include "fn_PRM_HMR_Parameters.hpp"
#include "fn_PRM_GEN_Parameters.hpp"
#include "fn_PRM_XTK_Parameters.hpp"
#include "fn_PRM_OPT_Parameters.hpp"
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
    uint tFEMFdScheme  = static_cast< uint >( fem::FDScheme_Type::POINT_3_CENTRAL );

    // FD step size in sweep
    std::string tFDsweep = "1.0e-5";

    // Number of constraints
    uint tNumConstraints = 3;

    // background mesh parameters
    std::string tNumElementsPerDir = tIs3D ? "4,4,4" : "4,4";
    std::string tDimensions        = tIs3D ? "1,1,1" : "1,1";
    std::string tOffSet            = tIs3D ? "0.0,0.0,0.0" : "0.0,0.0";
    std::string tSideSets          = tIs3D ? "1,2,3,4,5,6" : "1,2,3,4";

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
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
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
        Matrix< DDSMat > tConstraintTypes( tNumConstraints, 1, 1 );

        return tConstraintTypes;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_objectives( Matrix< DDRMat > aADVs, Matrix< DDRMat > aCriteria )
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
    compute_constraints( Matrix< DDRMat > aADVs, Matrix< DDRMat > aCriteria )
    {
        Matrix< DDRMat > tConstraints( 1, tNumConstraints );

        tConstraints( 0 ) = aCriteria( 0 );
        tConstraints( 1 ) = aCriteria( 1 );
        tConstraints( 2 ) = aCriteria( 2 );

        return tConstraints;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_dobjective_dadv( Matrix< DDRMat > aADVs, Matrix< DDRMat > aCriteria )
    {
        Matrix< DDRMat > tDObjectiveDADV( 1, aADVs.numel(), 0.0 );

        return tDObjectiveDADV;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_dobjective_dcriteria(
            Matrix< DDRMat > aADVs,
            Matrix< DDRMat > aCriteria )
    {
        Matrix< DDRMat > tDObjectiveDCriteria( 1, aCriteria.numel(), 0.0 );

        tDObjectiveDCriteria( 0 ) = 1.0;
        tDObjectiveDCriteria( 1 ) = 1.0;
        tDObjectiveDCriteria( 2 ) = 1.0;

        return tDObjectiveDCriteria;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_dconstraint_dadv(
            Matrix< DDRMat > aADVs,
            Matrix< DDRMat > aCriteria )
    {
        Matrix< DDRMat > tDConstraintDADV( tNumConstraints, aADVs.numel(), 0.0 );

        return tDConstraintDADV;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_dconstraint_dcriteria(
            Matrix< DDRMat > aADVs,
            Matrix< DDRMat > aCriteria )
    {
        Matrix< DDRMat > tDConstraintDCriteria( tNumConstraints, aCriteria.numel(), 0.0 );

        tDConstraintDCriteria( 0, 0 ) = 1.0;
        tDConstraintDCriteria( 1, 1 ) = 1.0;
        tDConstraintDCriteria( 2, 2 ) = 1.0;

        return tDConstraintDCriteria;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    OPTParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 3 );
        tParameterlist( 0 ).resize( 1 );
        tParameterlist( 1 ).resize( 0 );
        tParameterlist( 2 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = moris::prm::create_opt_problem_parameter_list();
        tParameterlist( 0 )( 0 ).set( "is_optimization_problem", tIsOpt );
        tParameterlist( 0 )( 0 ).set( "problem", "user_defined" );
        tParameterlist( 0 )( 0 ).set( "library", tLibraryName );
        tParameterlist( 0 )( 0 ).set( "restart_file", "" );

        tParameterlist( 2 )( 0 ) = moris::prm::create_sweep_parameter_list();
        tParameterlist( 2 )( 0 ).set( "hdf5_path", tHDF5FileName );
        tParameterlist( 2 )( 0 ).set( "num_evaluations_per_adv", "1" );
        tParameterlist( 2 )( 0 ).set( "finite_difference_type", "all" );
        tParameterlist( 2 )( 0 ).set( "finite_difference_epsilons", tFDsweep );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    HMRParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_hmr_parameter_list();

        tParameterlist( 0 )( 0 ).set( "number_of_elements_per_dimension", tNumElementsPerDir );
        tParameterlist( 0 )( 0 ).set( "domain_dimensions", tDimensions );
        tParameterlist( 0 )( 0 ).set( "domain_offset", tOffSet );
        tParameterlist( 0 )( 0 ).set( "domain_sidesets", tSideSets );
        tParameterlist( 0 )( 0 ).set( "lagrange_output_meshes", "0" );

        tParameterlist( 0 )( 0 ).set( "lagrange_orders", tLagrangeOrder );
        tParameterlist( 0 )( 0 ).set( "lagrange_pattern", "0" );

        tParameterlist( 0 )( 0 ).set( "bspline_orders", tBsplineOrder );
        tParameterlist( 0 )( 0 ).set( "bspline_pattern", "0" );

        tParameterlist( 0 )( 0 ).set( "initial_refinement", tInitialRef );
        tParameterlist( 0 )( 0 ).set( "initial_refinement_pattern", "0" );
        tParameterlist( 0 )( 0 ).set( "use_advanced_T_matrix_scheme", 1 );

        tParameterlist( 0 )( 0 ).set( "lagrange_to_bspline", "0" );

        tParameterlist( 0 )( 0 ).set( "truncate_bsplines", 1 );
        tParameterlist( 0 )( 0 ).set( "refinement_buffer", tRefineBuffer );
        tParameterlist( 0 )( 0 ).set( "staircase_buffer", tRefineBuffer );

        tParameterlist( 0 )( 0 ).set( "use_number_aura", 1 );

        tParameterlist( 0 )( 0 ).set( "use_multigrid", 0 );
        tParameterlist( 0 )( 0 ).set( "severity_level", 0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    XTKParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_xtk_parameter_list();
        tParameterlist( 0 )( 0 ).set( "decompose", true );
        tParameterlist( 0 )( 0 ).set( "decomposition_type", "conformal" );
        tParameterlist( 0 )( 0 ).set( "enrich", true );
        tParameterlist( 0 )( 0 ).set( "basis_rank", "bspline" );
        tParameterlist( 0 )( 0 ).set( "enrich_mesh_indices", "0" );
        tParameterlist( 0 )( 0 ).set( "multigrid", false );
        tParameterlist( 0 )( 0 ).set( "verbose", true );
        tParameterlist( 0 )( 0 ).set( "print_enriched_ig_mesh", false );

        tParameterlist( 0 )( 0 ).set( "ghost_stab", tUseGhost );

        tParameterlist( 0 )( 0 ).set( "exodus_output_XTK_ig_mesh", true );
        tParameterlist( 0 )( 0 ).set( "high_to_low_dbl_side_sets", true );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    GENParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 3 );

        tParameterlist( 0 ).push_back( moris::prm::create_gen_parameter_list() );
        tParameterlist( 0 )( 0 ).set(    //
                "IQI_types",
                std::string( "IQIBulkStrainEnergy,"
                             "IQIBulkStrainEnergyWithBedding,"
                             "IQIBulkVolume" ) );

        tParameterlist( 0 )( 0 ).set( "output_mesh_file", tGENOutputFile );
        tParameterlist( 0 )( 0 ).set( "time_offset", 10.0 );

        Matrix< DDUMat > tPhaseMap( 4, 1, 0 );
        tPhaseMap( 0 ) = 1;
        tPhaseMap( 1 ) = 0;
        tPhaseMap( 2 ) = 1;
        tPhaseMap( 3 ) = 1;
        tParameterlist( 0 )( 0 ).set( "phase_table", moris::ios::stringify( tPhaseMap ) );

        tParameterlist( 0 )( 0 ).set( "print_phase_table", true );

        // initialize advs
        tParameterlist( 0 )( 0 ).set( "initial_advs", std::to_string( 1.0 * tOuterRadius ) + ";" + std::to_string( 1.0 * tInnerRadius ) );
        tParameterlist( 0 )( 0 ).set( "lower_bounds", std::to_string( 0.9 * tOuterRadius ) + ";" + std::to_string( 0.9 * tInnerRadius ) );
        tParameterlist( 0 )( 0 ).set( "upper_bounds", std::to_string( 1.1 * tOuterRadius ) + ";" + std::to_string( 1.1 * tInnerRadius ) );

        // init geometry counter
        uint tGeoCounter = 0;

        if ( tIs3D )
        {
            tParameterlist( 1 ).push_back( prm::create_geometry_parameter_list() );
            tParameterlist( 1 )( tGeoCounter ).set( "type", "sphere" );
            tParameterlist( 1 )( tGeoCounter ).set( "field_variable_indices", "3" );
            tParameterlist( 1 )( tGeoCounter ).set( "adv_indices", "0" );
            tParameterlist( 1 )( tGeoCounter ).set( "constant_parameters", "0.5, 0.5, 0.5" );
            tGeoCounter++;

            tParameterlist( 1 ).push_back( prm::create_geometry_parameter_list() );
            tParameterlist( 1 )( tGeoCounter ).set( "type", "sphere" );
            tParameterlist( 1 )( tGeoCounter ).set( "field_variable_indices", "3" );
            tParameterlist( 1 )( tGeoCounter ).set( "adv_indices", "1" );
            tParameterlist( 1 )( tGeoCounter ).set( "constant_parameters", "0.5, 0.5, 0.5" );
            tGeoCounter++;
        }
        else
        {
            tParameterlist( 1 ).push_back( prm::create_geometry_parameter_list() );
            tParameterlist( 1 )( tGeoCounter ).set( "type", "circle" );
            tParameterlist( 1 )( tGeoCounter ).set( "field_variable_indices", "2" );
            tParameterlist( 1 )( tGeoCounter ).set( "adv_indices", "0" );
            tParameterlist( 1 )( tGeoCounter ).set( "constant_parameters", "0.5, 0.5" );
            tGeoCounter++;

            tParameterlist( 1 ).push_back( prm::create_geometry_parameter_list() );
            tParameterlist( 1 )( tGeoCounter ).set( "type", "circle" );
            tParameterlist( 1 )( tGeoCounter ).set( "field_variable_indices", "2" );
            tParameterlist( 1 )( tGeoCounter ).set( "adv_indices", "1" );
            tParameterlist( 1 )( tGeoCounter ).set( "constant_parameters", "0.5, 0.5" );
            tGeoCounter++;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    FEMParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterList )
    {
        // create a cell of cell of parameter list for fem
        tParameterList.resize( 8 );

        //------------------------------------------------------------------------------
        // init property counter
        uint tPropCounter = 0;

        // create parameter list for density property
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ) = prm::create_property_parameter_list();
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropDensity" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "1.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // create parameter list for Young's modulus property
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropYoungs" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "2.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // create parameter list for bedding property
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropBedding" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "2.0e-3" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // create parameter list for DBC property
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropDirichletU" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tDirichletStr );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // create parameter list for body load property
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "tPropBodyLoad" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tBodyLoadStr );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // create parameter list for Poisson ratio property
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropPoisson" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "0.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        //------------------------------------------------------------------------------
        // init CM counter
        uint tCMCounter = 0;

        // create parameter list for constitutive model
        tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_name", "CMStrucLinIso_Material" );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::STRUC_LIN_ISO ) );
        tParameterList( 1 )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( tDofStrg, "Displacement" ) );
        tParameterList( 1 )( tCMCounter ).set( "properties", "PropYoungs,YoungsModulus;PropPoisson,PoissonRatio" );
        tCMCounter++;

        //------------------------------------------------------------------------------
        // init SP counter
        uint tSPCounter = 0;

        // create parameter list for stabilization parameter 1
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", "SPNitscheDirichletBC" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::DIRICHLET_NITSCHE ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", std::to_string( tNitschePenalty ) );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", "PropYoungs,Material" );
        tSPCounter++;

        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", std::string( "SPGhost_Material" ) );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::GHOST_DISPL ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", std::string( "0.005" ) );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", std::string( "PropYoungs,Material" ) );
        tSPCounter++;

        //------------------------------------------------------------------------------
        // init IWG counter
        uint tIWGCounter = 0;

        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGBulkMaterial" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_BULK ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", tDofStrg );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", tDofStrg );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMStrucLinIso_Material,ElastLinIso" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties", "PropBedding,Bedding;tPropBodyLoad,Load" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tMaterialSets );
        tIWGCounter++;

        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGDirichletU" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", tDofStrg );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", tDofStrg );
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties", "PropDirichletU,Dirichlet" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMStrucLinIso_Material,ElastLinIso" );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPNitscheDirichletBC,DirichletNitsche" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tSupportSSets );
        tIWGCounter++;

        if ( tUseGhost )
        {
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name", std::string( "IWGGhostMaterial" ) );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::GHOST_NORMAL_FIELD ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual", tDofStrg );
            tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", tDofStrg );
            tParameterList( 3 )( tIWGCounter ).set( "follower_dof_dependencies", tDofStrg );
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", std::string( "SPGhost_Material,GhostSP" ) );
            tParameterList( 3 )( tIWGCounter ).set( "ghost_order", (uint)tDispOrder );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tMaterialGhost );
            tIWGCounter++;
        }

        //------------------------------------------------------------------------------
        // init IQI counter
        uint tIQICounter = 0;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkUX" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::DOF ) );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", tDofStrg );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", tDofStrg );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 0 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tTotalDomainSets );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkUY" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::DOF ) );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", tDofStrg );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", tDofStrg );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 1 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tTotalDomainSets );
        tIQICounter++;

        if ( tIs3D )
        {
            tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
            tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkUZ" );
            tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::DOF ) );
            tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", tDofStrg );
            tParameterList( 4 )( tIQICounter ).set( "dof_quantity", tDofStrg );
            tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 2 );
            tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tTotalDomainSets );
            tIQICounter++;
        }

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkStrainEnergy" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::STRAIN_ENERGY ) );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", tDofStrg );
        tParameterList( 4 )( tIQICounter ).set( "leader_constitutive_models", "CMStrucLinIso_Material,Elast" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tMaterialSets );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkStrainEnergyWithBedding" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::STRAIN_ENERGY ) );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", tDofStrg );
        tParameterList( 4 )( tIQICounter ).set( "leader_constitutive_models", "CMStrucLinIso_Material,Elast" );
        tParameterList( 4 )( tIQICounter ).set( "leader_properties", "PropBedding,Bedding" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tMaterialSets );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkVolume" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::VOLUME ) );
        tParameterList( 4 )( tIQICounter ).set( "leader_properties", "PropDensity,Density" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tMaterialSets );
        tIQICounter++;

        // create computation  parameter list
        tParameterList( 5 ).resize( 1 );
        tParameterList( 5 )( 0 ) = prm::create_computation_parameter_list();
        tParameterList( 5 )( 0 ).set( "print_physics_model", false );

        tParameterList( 5 )( 0 ).set( "finite_difference_scheme", tFEMFdScheme );
        tParameterList( 5 )( 0 ).set( "finite_difference_perturbation_size", tFEMFdEpsilon );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    SOLParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 8 );
        for ( uint Ik = 0; Ik < 8; Ik++ )
        {
            tParameterlist( Ik ).resize( 1 );
        }

        tParameterlist( 0 )( 0 ) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL );

        tParameterlist( 1 )( 0 ) = moris::prm::create_linear_solver_parameter_list();

        tParameterlist( 2 )( 0 ) = moris::prm::create_nonlinear_algorithm_parameter_list();
        tParameterlist( 2 )( 0 ).set( "NLA_combined_res_jac_assembly", true );
        tParameterlist( 2 )( 0 ).set( "NLA_rel_res_norm_drop", 1e-9 );
        tParameterlist( 2 )( 0 ).set( "NLA_relaxation_parameter", 1.00 );
        tParameterlist( 2 )( 0 ).set( "NLA_max_iter", 1 );

        tParameterlist( 3 )( 0 ) = moris::prm::create_nonlinear_solver_parameter_list();
        tParameterlist( 3 )( 0 ).set( "NLA_DofTypes", tDofStrg );

        tParameterlist( 4 )( 0 ) = moris::prm::create_time_solver_algorithm_parameter_list();

        tParameterlist( 5 )( 0 ) = moris::prm::create_time_solver_parameter_list();
        tParameterlist( 5 )( 0 ).set( "TSA_DofTypes", tDofStrg );
        tParameterlist( 5 )( 0 ).set( "TSA_Output_Indices", "0" );
        tParameterlist( 5 )( 0 ).set( "TSA_Output_Criteria", "Output_Criterion" );

        tParameterlist( 6 )( 0 ) = moris::prm::create_solver_warehouse_parameterlist();
        
        tParameterlist( 7 )( 0 ) = moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::NONE );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    MSIParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_msi_parameter_list();
        tParameterlist( 0 )( 0 ).set( "UX", 0 );
        tParameterlist( 0 )( 0 ).set( "UY", 0 );
        if ( tIs3D )
        {
            tParameterlist( 0 )( 0 ).set( "UZ", 0 );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    VISParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_vis_parameter_list();
        tParameterlist( 0 )( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName ) );
        tParameterlist( 0 )( 0 ).set( "Mesh_Type", static_cast< uint >( vis::VIS_Mesh_Type::STANDARD ) );
        tParameterlist( 0 )( 0 ).set( "Set_Names", tTotalDomainSets );

        if ( tIs3D )
        {
            tParameterlist( 0 )( 0 ).set( "Field_Names", std::string( "UX,UY,UZ,StrainEnergy,StrainEnergyWithBedding,Volume" ) );
            tParameterlist( 0 )( 0 ).set( "Field_Type", std::string( "NODAL,NODAL,NODAL,GLOBAL,GLOBAL,GLOBAL" ) );
            tParameterlist( 0 )( 0 ).set( "IQI_Names", std::string( "IQIBulkUX,IQIBulkUY,IQIBulkUZ,IQIBulkStrainEnergy,"
                                                                    "IQIBulkStrainEnergyWithBedding,IQIBulkVolume" ) );
        }
        else
        {
            tParameterlist( 0 )( 0 ).set( "Field_Names", std::string( "UX,UY,StrainEnergy,StrainEnergyWithBedding,Volume" ) );
            tParameterlist( 0 )( 0 ).set( "Field_Type", std::string( "NODAL,NODAL,GLOBAL,GLOBAL,GLOBAL" ) );
            tParameterlist( 0 )( 0 ).set( "IQI_Names", std::string( "IQIBulkUX,IQIBulkUY,IQIBulkStrainEnergy,"
                                                                    "IQIBulkStrainEnergyWithBedding,IQIBulkVolume" ) );
        }

        tParameterlist( 0 )( 0 ).set( "Save_Frequency", 1 );
        tParameterlist( 0 )( 0 ).set( "Time_Offset", 10.0 );
    }

    void
    MORISGENERALParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
    }

    //--------------------------------------------------------------------------------------------------------------
}    // namespace moris

//--------------------------------------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
