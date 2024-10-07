/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Shape_Sensitivity_Two_Material_Sweep_Thermoelastic_Staggered.cpp
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

#include "AztecOO.h"

#ifdef __cplusplus
extern "C" {
#endif
//------------------------------------------------------------------------------
namespace moris
{
    /* ------------------------------------------------------------------------ */
    // Mesh Set Information

    // Phase 1: back  - Material 1
    // Phase 2: front - Material 2

    std::string tPhase1 = "HMR_dummy_n_p1,HMR_dummy_c_p1";

    std::string tPhase2 = "HMR_dummy_n_p0,HMR_dummy_c_p0";

    std::string tBackSurface = "SideSet_4_n_p1";

    std::string tFrontSurface = "SideSet_2_n_p1";

    std::string tPhase1Ghost = "ghost_p1";

    std::string tTotalDomain = tPhase1 + tPhase2;

    /* ------------------------------------------------------------------------ */
    // geometry parameters

    // dimensionality 2D or 3D
    uint sDim = 2;

    // general
    moris::real sL  = 6.0;    // total length
    moris::real sL1 = 4.1;
    moris::real sL2 = sL - sL1;

    /* ------------------------------------------------------------------------ */
    // material parameters

    // flux
    moris::real sP2 = 5.0;

    // capacity
    std::string sCap1 = "1.0";

    // density
    std::string tDens1 = "1.0";

    std::string tCTE1    = "1.0";
    std::string tCTE2    = "2.0";
    std::string tRefTemp = "0.0";

    std::string tEmod = "1.0";
    std::string tPois = "0.0";

    // conductivity
    moris::real sK1 = 1.0;
    moris::real sK2 = 2.0;

    // body flux
    moris::real sQ1 = 1.0;

    /* ------------------------------------------------------------------------ */
    // HMR parameters

    std::string tNumElemsPerDim = sDim == 2 ? "6,   4" : "6,   4,   4";
    std::string tDomainDims     = sDim == 2 ? std::to_string( sL ) + ", 4.0" : std::to_string( sL ) + ", 4.0, 4.0";
    std::string tDomainOffset   = sDim == 2 ? "0.0,  0.0" : "0.0,  0.0, 0.0";
    std::string tDomainSidesets = sDim == 2 ? "1,2,3,4" : "1,2,3,4,5,6";

    std::string tInterpolationOrder = "1";

    int tRefineBuffer = 1;

    int tInterfaceRefinement = 1;

    /* ------------------------------------------------------------------------ */
    // Solver config

    moris::real tNLA_rel_res_norm_drop    = 1.0e-12;
    moris::real tNLA_relaxation_parameter = 1.0;
    int         tNLA_max_iter             = 4;

    int         tTSA_Num_Time_Steps = 1;
    moris::real tTSA_Time_Frame     = 1.0e0;

    /* ------------------------------------------------------------------------ */
    // Minimum level set value
    bool tUseGhost = true;

    /* ------------------------------------------------------------------------ */
    // Output Config

    std::string tOutputFileName = "ShapeSensitivitiesTransientCircle.exo";

    /* ------------------------------------------------------------------------ */
    // Constant function for properties

    void
    Func_Const( moris::Matrix<
                        moris::DDRMat >                   &aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > > &aParameters,
            moris::fem::Field_Interpolator_Manager        *aFIManager )
    {
        aPropMatrix = aParameters( 0 );
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
        tObjectives( 0 ) = aCriteria( 0 );

        return tObjectives;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_constraints( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tConstraints( 1, 1 );
        tConstraints( 0 ) = aCriteria( 1 );

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
        tDObjectiveDCriteria( 0 ) = 1.0;
        tDObjectiveDCriteria( 1 ) = 0.0;

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
        tDConstraintDCriteria( 0 ) = 0.0;
        tDConstraintDCriteria( 1 ) = 1.0;

        return tDConstraintDCriteria;
    }

    /* ------------------------------------------------------------------------ */

    bool
    Output_Criterion( moris::tsa::Time_Solver *aTimeSolver )
    {
        return false;
    }

    /* ------------------------------------------------------------------------ */

    void
    OPTParameterList( Module_Parameter_Lists &aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( moris::prm::create_opt_problem_parameter_list() );
        aParameterLists( 0 ).set( "is_optimization_problem", true );
        aParameterLists( 0 ).set( "problem", "user_defined" );
        aParameterLists( 0 ).set( "library", "./Shape_Sensitivity_Two_Material_Sweep_Thermoelastic_Staggered.so" );

        aParameterLists( 2 ).add_parameter_list( moris::prm::create_sweep_parameter_list() );
        aParameterLists( 2 ).set( "hdf5_path", "shape_opt_test.hdf5" );
        aParameterLists( 2 ).set( "num_evaluations_per_adv", "1" );
        aParameterLists( 2 ).set( "finite_difference_type", "all" );
        aParameterLists( 2 ).set( "finite_difference_epsilons", "1e-6" );
    }

    void
    HMRParameterList( Module_Parameter_Lists &aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_hmr_parameter_list() );

        aParameterLists( 0 ).set( "number_of_elements_per_dimension", tNumElemsPerDim );
        aParameterLists( 0 ).set( "domain_dimensions", tDomainDims );
        aParameterLists( 0 ).set( "domain_offset", tDomainOffset );
        aParameterLists( 0 ).set( "domain_sidesets", tDomainSidesets );
        aParameterLists( 0 ).set( "lagrange_output_meshes", "0" );

        aParameterLists( 0 ).set( "lagrange_orders", tInterpolationOrder );
        aParameterLists( 0 ).set( "lagrange_pattern", "0" );
        aParameterLists( 0 ).set( "bspline_orders", tInterpolationOrder );
        aParameterLists( 0 ).set( "bspline_pattern", "0" );

        aParameterLists( 0 ).set( "lagrange_to_bspline", "0" );

        aParameterLists( 0 ).set( "truncate_bsplines", 1 );
        aParameterLists( 0 ).set( "refinement_buffer", tRefineBuffer );
        aParameterLists( 0 ).set( "staircase_buffer", tRefineBuffer );
        aParameterLists( 0 ).set( "initial_refinement", "1" );
        aParameterLists( 0 ).set( "initial_refinement_pattern", "0" );

        aParameterLists( 0 ).set( "use_number_aura", 1 );

        aParameterLists( 0 ).set( "use_multigrid", 0 );
        aParameterLists( 0 ).set( "severity_level", 0 );
    }

    void
    XTKParameterList( Module_Parameter_Lists &aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_xtk_parameter_list() );
        aParameterLists( 0 ).set( "decompose", true );
        aParameterLists( 0 ).set( "decomposition_type", "conformal" );
        aParameterLists( 0 ).set( "enrich", true );
        aParameterLists( 0 ).set( "basis_rank", "bspline" );
        aParameterLists( 0 ).set( "enrich_mesh_indices", "0" );
        aParameterLists( 0 ).set( "ghost_stab", tUseGhost );
        aParameterLists( 0 ).set( "multigrid", false );
        aParameterLists( 0 ).set( "verbose", true );
        aParameterLists( 0 ).set( "print_enriched_ig_mesh", false );
        aParameterLists( 0 ).set( "exodus_output_XTK_ig_mesh", true );
        aParameterLists( 0 ).set( "high_to_low_dbl_side_sets", true );
    }

    void
    GENParameterList( Module_Parameter_Lists &aParameterLists )
    {

        aParameterLists( 0 ).add_parameter_list( moris::prm::create_gen_parameter_list() );
        aParameterLists( 0 ).set( "IQI_types", "IQIBulkStrainEnergyDISP", "IQIBulkVolume" );

        // Geometry parameter lists

        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::CIRCLE ) );
        aParameterLists( 1 ).set( "center_x", 3.0, 3.0, 3.0 );
        aParameterLists( 1 ).set( "center_y", 2.21, 2.21, 2.21 );
        aParameterLists( 1 ).set( "radius", 1.4, 1.4, 1.4 );
    }

    void
    FEMParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.hack_for_legacy_fem();
        // create a cell of cell of parameter list for fem

        //------------------------------------------------------------------------------

        // properties for material 1
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropDensity1" );
        aParameterLists( 0 ).set( "function_parameters", tDens1 );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropCapacity1" );
        aParameterLists( 0 ).set( "function_parameters", sCap1 );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropConductivity1" );
        aParameterLists( 0 ).set( "function_parameters", std::to_string( sK1 ) );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropConductivity2" );
        aParameterLists( 0 ).set( "function_parameters", std::to_string( sK2 ) );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropHeatLoad1" );
        aParameterLists( 0 ).set( "function_parameters", std::to_string( sQ1 ) );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // surface flux
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropSurfaceFlux" );
        aParameterLists( 0 ).set( "function_parameters", std::to_string( sP2 ) );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // temperature at back surface
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropImposedTemperature" );
        aParameterLists( 0 ).set( "function_parameters", "0.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // properties of boundary conditions
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropDirichlet" );
        aParameterLists( 0 ).set( "function_parameters", "0.0;0.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // time continuity weights
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropWeightCurrent" );
        aParameterLists( 0 ).set( "function_parameters", "100.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropWeightPrevious" );
        aParameterLists( 0 ).set( "function_parameters", "100.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // initial condition
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropInitialCondition" );
        aParameterLists( 0 ).set( "function_parameters", "0.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropCTE1" );
        aParameterLists( 0 ).set( "function_parameters", tCTE1 );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropCTE2" );
        aParameterLists( 0 ).set( "function_parameters", tCTE2 );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropRefTemp" );
        aParameterLists( 0 ).set( "function_parameters", tRefTemp );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropYoungs" );
        aParameterLists( 0 ).set( "function_parameters", tEmod );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropPoisson" );
        aParameterLists( 0 ).set( "function_parameters", tPois );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        //------------------------------------------------------------------------------

        // create parameter list for constitutive model 1
        aParameterLists( 1 ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists( 1 ).set( "constitutive_name", "CMStrucLinIso1" );
        aParameterLists( 1 ).set( "constitutive_type",  fem::Constitutive_Type::STRUC_LIN_ISO ) ;
        aParameterLists( 1 ).set( "model_type",  fem::Model_Type::PLANE_STRESS ) ;
        aParameterLists( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY;TEMP", "Displacement,Temperature" ) );
        aParameterLists( 1 ).set( "properties",
                "PropYoungs, YoungsModulus;"
                "PropPoisson,PoissonRatio;"
                "PropCTE1,    CTE;"
                "PropRefTemp,ReferenceTemperature" );

        // create parameter list for constitutive model 1
        aParameterLists( 1 ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists( 1 ).set( "constitutive_name", "CMStrucLinIso2" );
        aParameterLists( 1 ).set( "constitutive_type",  fem::Constitutive_Type::STRUC_LIN_ISO ) ;
        aParameterLists( 1 ).set( "model_type",  fem::Model_Type::PLANE_STRESS ) ;
        aParameterLists( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY;TEMP", "Displacement,Temperature" ) );
        aParameterLists( 1 ).set( "properties",
                "PropYoungs, YoungsModulus;"
                "PropPoisson,PoissonRatio;"
                "PropCTE2,    CTE;"
                "PropRefTemp,ReferenceTemperature" );

        // create parameter list for constitutive model - Inclusion
        aParameterLists( 1 ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists( 1 ).set( "constitutive_name", "CMDiffusion1" );
        aParameterLists( 1 ).set( "constitutive_type",  fem::Constitutive_Type::DIFF_LIN_ISO ) ;
        aParameterLists( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists( 1 ).set( "properties", "PropConductivity1 , Conductivity" );

        // create parameter list for constitutive model - Inclusion
        aParameterLists( 1 ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists( 1 ).set( "constitutive_name", "CMDiffusion2" );
        aParameterLists( 1 ).set( "constitutive_type",  fem::Constitutive_Type::DIFF_LIN_ISO ) ;
        aParameterLists( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists( 1 ).set( "properties", "PropConductivity2 , Conductivity" );

        //------------------------------------------------------------------------------

        // Nitsche stabilization parameter for structure
        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPNitscheStruc" );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::DIRICHLET_NITSCHE ) ;
        aParameterLists( 2 ).set( "function_parameters", "100.0" );
        aParameterLists( 2 ).set( "leader_properties", "PropYoungs,Material" );

        // create parameter list for DBC on back surface
        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPNitscheTemp" );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::DIRICHLET_NITSCHE ) ;
        aParameterLists( 2 ).set( "function_parameters", "100.0" );
        aParameterLists( 2 ).set( "leader_properties", "PropConductivity1,Material" );

        // Dirichlet SP
        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPInterfaceNitscheDISP" );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::NITSCHE_INTERFACE ) ;
        aParameterLists( 2 ).set( "function_parameters", "100.0" );
        aParameterLists( 2 ).set( "leader_properties", "PropYoungs,Material" );
        aParameterLists( 2 ).set( "follower_properties", "PropYoungs,Material" );

        // Dirichlet SP
        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPInterfaceNitscheTEMP" );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::NITSCHE_INTERFACE ) ;
        aParameterLists( 2 ).set( "function_parameters", "100.0" );
        aParameterLists( 2 ).set( "leader_properties", "PropConductivity1,Material" );
        aParameterLists( 2 ).set( "follower_properties", "PropConductivity2,Material" );

        if ( tUseGhost )
        {
            // bulk Ghost - Shell - Temperature
            aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
            aParameterLists( 2 ).set( "stabilization_name", "SPGPTemp1" );
            aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::GHOST_DISPL ) ;
            aParameterLists( 2 ).set( "function_parameters", "0.01" );
            aParameterLists( 2 ).set( "leader_properties", "PropConductivity1,Material" );

            // bulk Ghost - PCM - Temperature
            aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
            aParameterLists( 2 ).set( "stabilization_name", "SPGPTemp2" );
            aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::GHOST_DISPL ) ;
            aParameterLists( 2 ).set( "function_parameters", "0.01" );
            aParameterLists( 2 ).set( "leader_properties", "PropConductivity2,Material" );

            // bulk Ghost - Shell - Displacements
            aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
            aParameterLists( 2 ).set( "stabilization_name", "SPGPStruct1" );
            aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::GHOST_DISPL ) ;
            aParameterLists( 2 ).set( "function_parameters", "0.01" );
            aParameterLists( 2 ).set( "leader_properties", "PropYoungs,Material" );

            // bulk Ghost - PCM - Displacements
            aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
            aParameterLists( 2 ).set( "stabilization_name", "SPGPStruct2" );
            aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::GHOST_DISPL ) ;
            aParameterLists( 2 ).set( "function_parameters", "0.01" );
            aParameterLists( 2 ).set( "leader_properties", "PropYoungs,Material" );
        }

        //------------------------------------------------------------------------------
        // create IWG - bulk structure
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGBulkStruct1" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_BULK ) ;
        aParameterLists( 3 ).set( "dof_residual", "UX,UY" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY;TEMP" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        aParameterLists( 3 ).set( "mesh_set_names", tPhase1 );

        // create IWG - bulk structure
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGBulkStruct2" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_BULK ) ;
        aParameterLists( 3 ).set( "dof_residual", "UX,UY" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY;TEMP" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMStrucLinIso2,ElastLinIso" );
        aParameterLists( 3 ).set( "mesh_set_names", tPhase2 );

        // create IWG for inclusion - bulk diffusion
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGDiffusionBulk1" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_BULK ) ;
        aParameterLists( 3 ).set( "dof_residual", "TEMP" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY;TEMP" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMDiffusion1,Diffusion" );
        // aParameterLists( 3 ).set( "leader_properties",          "PropHeatLoad1,Load");
        aParameterLists( 3 ).set( "mesh_set_names", tPhase1 );

        // create IWG for inclusion - bulk diffusion
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGDiffusionBulk2" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_BULK ) ;
        aParameterLists( 3 ).set( "dof_residual", "TEMP" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY;TEMP" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMDiffusion2,Diffusion" );
        // aParameterLists( 3 ).set( "leader_properties",          "PropHeatLoad1,Load");
        aParameterLists( 3 ).set( "mesh_set_names", tPhase2 );

        // create IWG - Dirichlet structure
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGDirichletDISP" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE ) ;
        aParameterLists( 3 ).set( "dof_residual", "UX,UY" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY;TEMP" );
        aParameterLists( 3 ).set( "leader_properties", "PropDirichlet,Dirichlet" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPNitscheStruc,DirichletNitsche" );
        aParameterLists( 3 ).set( "mesh_set_names", tBackSurface );

        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGDirichletDISP2" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE ) ;
        aParameterLists( 3 ).set( "dof_residual", "UX,UY" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY;TEMP" );
        aParameterLists( 3 ).set( "leader_properties", "PropDirichlet,Dirichlet" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPNitscheStruc,DirichletNitsche" );
        aParameterLists( 3 ).set( "mesh_set_names", tFrontSurface );

        // create IWG for Neumann boundary conditions
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGInletFlux" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_NEUMANN ) ;
        aParameterLists( 3 ).set( "dof_residual", "TEMP" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY;TEMP" );
        aParameterLists( 3 ).set( "leader_properties", "PropSurfaceFlux,Neumann" );
        aParameterLists( 3 ).set( "mesh_set_names", tFrontSurface );

        // create IWG for Dirichlet boundary conditions
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWG2SurfaceTemp" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        aParameterLists( 3 ).set( "dof_residual", "TEMP" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY;TEMP" );
        aParameterLists( 3 ).set( "leader_properties", "PropImposedTemperature,Dirichlet" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMDiffusion1,Diffusion" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPNitscheTemp,DirichletNitsche" );
        aParameterLists( 3 ).set( "mesh_set_names", tBackSurface );

        // Interface Dirichlet BC
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGInterface" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE ) ;
        aParameterLists( 3 ).set( "dof_residual", "TEMP" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY;TEMP" );
        aParameterLists( 3 ).set( "follower_dof_dependencies", "UX,UY;TEMP" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMDiffusion1,Diffusion" );
        aParameterLists( 3 ).set( "follower_constitutive_models", "CMDiffusion2,Diffusion" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPInterfaceNitscheTEMP ,NitscheInterface" );
        aParameterLists( 3 ).set( "mesh_set_names", "dbl_iside_p0_1_p1_0" );

        // Interface Dirichlet BC
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGInterface" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_INTERFACE_SYMMETRIC_NITSCHE ) ;
        aParameterLists( 3 ).set( "dof_residual", "UX,UY" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY;TEMP" );
        aParameterLists( 3 ).set( "follower_dof_dependencies", "UX,UY;TEMP" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        aParameterLists( 3 ).set( "follower_constitutive_models", "CMStrucLinIso2,ElastLinIso" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPInterfaceNitscheDISP ,NitscheInterface" );
        aParameterLists( 3 ).set( "mesh_set_names", "dbl_iside_p0_1_p1_0" );

        if ( tUseGhost )
        {
            // temperature - Shell
            aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists( 3 ).set( "IWG_name", "IWGGPTemp1" );
            aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
            aParameterLists( 3 ).set( "dof_residual", "TEMP" );
            aParameterLists( 3 ).set( "leader_dof_dependencies", "TEMP" );
            aParameterLists( 3 ).set( "follower_dof_dependencies", "TEMP" );
            aParameterLists( 3 ).set( "stabilization_parameters", "SPGPTemp1,GhostSP" );
            aParameterLists( 3 ).set( "mesh_set_names", "ghost_p1" );

            // displacements - Shell
            aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists( 3 ).set( "IWG_name", "IWGGPStruct1" );
            aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
            aParameterLists( 3 ).set( "dof_residual", "UX,UY" );
            aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY" );
            aParameterLists( 3 ).set( "follower_dof_dependencies", "UX,UY" );
            aParameterLists( 3 ).set( "stabilization_parameters", "SPGPStruct1,GhostSP" );
            aParameterLists( 3 ).set( "mesh_set_names", "ghost_p1" );

            // temperature - PCM
            aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists( 3 ).set( "IWG_name", "IWGGPTemp2" );
            aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
            aParameterLists( 3 ).set( "dof_residual", "TEMP" );
            aParameterLists( 3 ).set( "leader_dof_dependencies", "TEMP" );
            aParameterLists( 3 ).set( "follower_dof_dependencies", "TEMP" );
            aParameterLists( 3 ).set( "stabilization_parameters", "SPGPTemp2,GhostSP" );
            aParameterLists( 3 ).set( "mesh_set_names", "ghost_p0" );

            // displacements - PCM
            aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists( 3 ).set( "IWG_name", "SPGPStruct2" );
            aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
            aParameterLists( 3 ).set( "dof_residual", "UX,UY" );
            aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY" );
            aParameterLists( 3 ).set( "follower_dof_dependencies", "UX,UY" );
            aParameterLists( 3 ).set( "stabilization_parameters", "SPGPStruct2,GhostSP" );
            aParameterLists( 3 ).set( "mesh_set_names", "ghost_p0" );
        }

        //------------------------------------------------------------------------------
        // Nodal Temperature IQI
        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkTEMP" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists( 4 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 4 ).set( "vectorial_field_index", 0 );
        aParameterLists( 4 ).set( "mesh_set_names", tPhase1 );

        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkStrainEnergyDISP" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::STRAIN_ENERGY ) ;
        aParameterLists( 4 ).set( "leader_dof_dependencies", "UX,UY;TEMP" );
        aParameterLists( 4 ).set( "leader_constitutive_models", "CMStrucLinIso2,Elast" );
        aParameterLists( 4 ).set( "mesh_set_names", tPhase2 );

        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkVolume" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::VOLUME ) ;
        aParameterLists( 4 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 4 ).set( "mesh_set_names", tPhase2 );

        // Max Temperature IQI
        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIMaxTemp" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::MAX_DOF ) ;
        aParameterLists( 4 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 4 ).set( "vectorial_field_index", 0 );
        aParameterLists( 4 ).set( "function_parameters", "1.0/2.0" );
        aParameterLists( 4 ).set( "mesh_set_names", tPhase1 );

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        aParameterLists( 5 ).add_parameter_list( prm::create_computation_parameter_list() );

        aParameterLists( 5 ).set( "finite_difference_scheme", ( fem::FDScheme_Type::POINT_3_CENTRAL ) );
        aParameterLists( 5 ).set( "finite_difference_perturbation_size", 1.0e-4 );
    }

    void
    SOLParameterList( Module_Parameter_Lists &aParameterLists )
    {

        aParameterLists( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL ) );

        aParameterLists( 1 ).add_parameter_list( moris::prm::create_linear_solver_parameter_list() );

        aParameterLists( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
        aParameterLists( 2 ).set( "NLA_Solver_Implementation",  moris::NLA::NonlinearSolverType::NEWTON_SOLVER ) ;
        aParameterLists( 2 ).set( "NLA_combined_res_jac_assembly", false );
        aParameterLists( 2 ).set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        aParameterLists( 2 ).set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
        aParameterLists( 2 ).set( "NLA_max_iter", tNLA_max_iter );
        aParameterLists( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
        aParameterLists( 2 ).set( "NLA_Solver_Implementation",  moris::NLA::NonlinearSolverType::NLBGS_SOLVER ) ;

        aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
        aParameterLists( 3 ).set( "NLA_Solver_Implementation",  moris::NLA::NonlinearSolverType::NEWTON_SOLVER ) ;
        aParameterLists( 3 ).set( "NLA_DofTypes", "UX,UY" );
        aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
        aParameterLists( 3 ).set( "NLA_Solver_Implementation",  moris::NLA::NonlinearSolverType::NEWTON_SOLVER ) ;
        aParameterLists( 3 ).set( "NLA_DofTypes", "TEMP" );
        aParameterLists( 3 ).set( "NLA_Secondary_DofTypes", "UX,UY" );
        aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
        aParameterLists( 3 ).set( "NLA_Solver_Implementation",  moris::NLA::NonlinearSolverType::NLBGS_SOLVER ) ;
        aParameterLists( 3 ).set( "NLA_Sub_Nonlinear_Solver", "1,0" );
        aParameterLists( 3 ).set( "NLA_DofTypes", "UX,UY;TEMP" );
        aParameterLists( 3 ).set( "NLA_Nonlinear_solver_algorithms", "1" );

        aParameterLists( 4 ).add_parameter_list( moris::prm::create_time_solver_algorithm_parameter_list() );
        aParameterLists( 4 ).set( "TSA_Num_Time_Steps", tTSA_Num_Time_Steps );
        aParameterLists( 4 ).set( "TSA_Time_Frame", tTSA_Time_Frame );
        aParameterLists( 4 ).set( "TSA_Nonlinear_Solver", 2 );

        aParameterLists( 5 ).add_parameter_list( moris::prm::create_time_solver_parameter_list() );
        aParameterLists( 5 ).set( "TSA_DofTypes", "UX,UY;TEMP" );
        aParameterLists( 5 ).set( "TSA_Initialize_Sol_Vec", "UX,0.0;UY,0.0;TEMP,0.0" );
        aParameterLists( 5 ).set( "TSA_Output_Indices", "0" );
        aParameterLists( 5 ).set( "TSA_Output_Criteria", "Output_Criterion" );

        aParameterLists( 6 ).add_parameter_list( moris::prm::create_solver_warehouse_parameterlist() );
        aParameterLists( 6 ).set( "SOL_save_final_adjoint_vec_to_file", "Shape_Sensitivity_Two_Material_Staggered.hdf5" );

        aParameterLists( 7 ).add_parameter_list( moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::NONE ) );
    }

    void
    MSIParameterList( Module_Parameter_Lists &aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_msi_parameter_list() );
        aParameterLists( 0 ).set( "order_adofs_by_host", false );
    }

    void
    VISParameterList( Module_Parameter_Lists &aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_vis_parameter_list() );
        aParameterLists( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName ) );
        aParameterLists( 0 ).set( "Mesh_Type",  vis::VIS_Mesh_Type::STANDARD ) ;
        aParameterLists( 0 ).set( "Set_Names", tPhase1 );
        aParameterLists( 0 ).set( "Field_Names", "TEMP,MAX_DOF" );
        aParameterLists( 0 ).set( "Field_Type", "NODAL,GLOBAL" );
        aParameterLists( 0 ).set( "IQI_Names", "IQIBulkTEMP,IQIMaxTemp" );
        aParameterLists( 0 ).set( "Save_Frequency", 1 );
    }

    void
    MORISGENERALParameterList( Module_Parameter_Lists &aParameterLists )
    {
    }

    /* ------------------------------------------------------------------------ */
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
