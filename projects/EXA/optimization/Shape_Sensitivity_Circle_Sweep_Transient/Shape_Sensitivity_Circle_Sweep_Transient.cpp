/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Shape_Sensitivity_Circle_Sweep_Transient.cpp
 *
 */

#include <string>
#include <iostream>
#include "typedefs.hpp"
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

    std::string tBackSurface = "SideSet_4_n_p1";

    std::string tFrontSurface = "SideSet_2_n_p1";

    std::string tPhase1Ghost = "ghost_p1";

    std::string tTotalDomain = tPhase1;

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

    // conductivity
    moris::real sK1 = 1.0;

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

    moris::real tNLA_rel_res_norm_drop    = 1.0e-08;
    moris::real tNLA_relaxation_parameter = 1.0;
    int         tNLA_max_iter             = 2;

    int         tTSA_Num_Time_Steps = 2;
    moris::real tTSA_Time_Frame     = 1.0e0;

    /* ------------------------------------------------------------------------ */
    // Minimum level set value
    bool tUseGhost = false;

    /* ------------------------------------------------------------------------ */
    // Output Config

    std::string tOutputFileName = "ShapeSensitivitiesTransientCircle.exo";

    /* ------------------------------------------------------------------------ */
    // Constant function for properties

    void
    Func_Const( moris::Matrix<
                        moris::DDRMat >&                   aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
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
    compute_objectives( Matrix< DDRMat > aADVs, Matrix< DDRMat > aCriteria )
    {
        Matrix< DDRMat > tObjectives( 1, 1 );
        tObjectives( 0 ) = aCriteria( 0 );

        return tObjectives;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_constraints( Matrix< DDRMat > aADVs, Matrix< DDRMat > aCriteria )
    {
        Matrix< DDRMat > tConstraints( 1, 1 );
        tConstraints( 0 ) = aCriteria( 1 );

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
    compute_dobjective_dcriteria( Matrix< DDRMat > aADVs, Matrix< DDRMat > aCriteria )
    {
        Matrix< DDRMat > tDObjectiveDCriteria( 1, 2 );
        tDObjectiveDCriteria( 0 ) = 1.0;
        tDObjectiveDCriteria( 1 ) = 0.0;

        return tDObjectiveDCriteria;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_dconstraint_dadv( Matrix< DDRMat > aADVs, Matrix< DDRMat > aCriteria )
    {
        Matrix< DDRMat > tDConstraintDADV( 1, aADVs.numel(), 0.0 );

        return tDConstraintDADV;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_dconstraint_dcriteria( Matrix< DDRMat > aADVs, Matrix< DDRMat > aCriteria )
    {
        Matrix< DDRMat > tDConstraintDCriteria( 1, 2 );
        tDConstraintDCriteria( 0 ) = 0.0;
        tDConstraintDCriteria( 1 ) = 1.0;

        return tDConstraintDCriteria;
    }

    /* ------------------------------------------------------------------------ */

    bool
    Output_Criterion( moris::tsa::Time_Solver* aTimeSolver )
    {
        return false;
    }

    /* ------------------------------------------------------------------------ */

    void
    OPTParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 3 );
        tParameterlist( 0 ).resize( 1 );
        tParameterlist( 1 ).resize( 0 );
        tParameterlist( 2 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = moris::prm::create_opt_problem_parameter_list();
        tParameterlist( 0 )( 0 ).set( "is_optimization_problem", true );
        tParameterlist( 0 )( 0 ).set( "problem", "user_defined" );
        tParameterlist( 0 )( 0 ).set( "library", "./Shape_Sensitivity_Circle_Sweep_Transient.so" );

        tParameterlist( 2 )( 0 ) = moris::prm::create_sweep_parameter_list();
        tParameterlist( 2 )( 0 ).set( "hdf5_path", "shape_opt_test.hdf5" );
        tParameterlist( 2 )( 0 ).set( "num_evaluations_per_adv", "1" );
        tParameterlist( 2 )( 0 ).set( "finite_difference_type", "all" );
    }

    void
    HMRParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_hmr_parameter_list();

        tParameterlist( 0 )( 0 ).set( "number_of_elements_per_dimension", tNumElemsPerDim );
        tParameterlist( 0 )( 0 ).set( "domain_dimensions", tDomainDims );
        tParameterlist( 0 )( 0 ).set( "domain_offset", tDomainOffset );
        tParameterlist( 0 )( 0 ).set( "domain_sidesets", tDomainSidesets );
        tParameterlist( 0 )( 0 ).set( "lagrange_output_meshes", "0" );

        tParameterlist( 0 )( 0 ).set( "lagrange_orders", tInterpolationOrder );
        tParameterlist( 0 )( 0 ).set( "lagrange_pattern", "0" );
        tParameterlist( 0 )( 0 ).set( "bspline_orders", tInterpolationOrder );
        tParameterlist( 0 )( 0 ).set( "bspline_pattern", "0" );

        tParameterlist( 0 )( 0 ).set( "lagrange_to_bspline", "0" );

        tParameterlist( 0 )( 0 ).set( "truncate_bsplines", 1 );
        tParameterlist( 0 )( 0 ).set( "refinement_buffer", tRefineBuffer );
        tParameterlist( 0 )( 0 ).set( "staircase_buffer", tRefineBuffer );
        tParameterlist( 0 )( 0 ).set( "initial_refinement", "1" );
        tParameterlist( 0 )( 0 ).set( "initial_refinement_pattern", "0" );

        tParameterlist( 0 )( 0 ).set( "use_number_aura", 1 );

        tParameterlist( 0 )( 0 ).set( "use_multigrid", 0 );
        tParameterlist( 0 )( 0 ).set( "severity_level", 0 );
    }

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
        tParameterlist( 0 )( 0 ).set( "ghost_stab", tUseGhost );
        tParameterlist( 0 )( 0 ).set( "multigrid", false );
        tParameterlist( 0 )( 0 ).set( "verbose", true );
        tParameterlist( 0 )( 0 ).set( "print_enriched_ig_mesh", false );
        tParameterlist( 0 )( 0 ).set( "exodus_output_XTK_ig_mesh", true );
        tParameterlist( 0 )( 0 ).set( "high_to_low_dbl_side_sets", true );
    }

    void
    GENParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {

        tParameterlist.resize( 3 );
        tParameterlist( 0 ).resize( 1 );
        tParameterlist( 0 )( 0 ) = moris::prm::create_gen_parameter_list();
        tParameterlist( 0 )( 0 ).set( "initial_advs", "3.0; 2.21; 1.4" );
        tParameterlist( 0 )( 0 ).set( "lower_bounds", "3.0; 2.21; 1.4" );
        tParameterlist( 0 )( 0 ).set( "upper_bounds", "3.0; 2.21; 1.4" );
        tParameterlist( 0 )( 0 ).set( "IQI_types", "IQIBulkStrainEnergy,IQIBulkVolume" );

        // Geometry parameter lists
        tParameterlist( 1 ).resize( 1 );

        tParameterlist( 1 )( 0 ) = prm::create_geometry_parameter_list();
        tParameterlist( 1 )( 0 ).set( "type", "circle" );
        tParameterlist( 1 )( 0 ).set( "field_variable_indices", "0, 1, 2" );
        tParameterlist( 1 )( 0 ).set( "adv_indices", "0, 1, 2" );
    }

    void
    FEMParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterList )
    {
        // create a cell of cell of parameter list for fem
        tParameterList.resize( 8 );

        //------------------------------------------------------------------------------
        // init property counter
        uint tPropCounter = 0;

        // properties for material 1
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ) = prm::create_property_parameter_list();
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropDensity1" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tDens1 );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropCapacity1" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", sCap1 );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropConductivity1" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", std::to_string( sK1 ) );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropHeatLoad1" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", std::to_string( sQ1 ) );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // surface flux
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropSurfaceFlux" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", std::to_string( sP2 ) );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // temperature at back surface
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropImposedTemperature" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "0.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // time continuity weights
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropWeightCurrent" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "100.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropWeightPrevious" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "100.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // initial condition
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropInitialCondition" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "0.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        //------------------------------------------------------------------------------
        // init CM counter
        uint tCMCounter = 0;

        // create parameter list for constitutive model - Inclusion
        tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_name", "CMDiffusion1" );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::DIFF_LIN_ISO ) );
        tParameterList( 1 )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( 1 )( tCMCounter ).set( "properties",
                "PropConductivity1 , Conductivity;"
                "PropDensity1      , Density;"
                "PropCapacity1     , HeatCapacity" );
        tCMCounter++;

        //------------------------------------------------------------------------------
        // init SP counter
        uint tSPCounter = 0;

        // create parameter list for DBC on back surface
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", "SPNitscheTemp" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::DIRICHLET_NITSCHE ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", "10.0" );
        tParameterList( 2 )( tSPCounter ).set( "master_properties", "PropConductivity1,Material" );
        tSPCounter++;

        if ( tUseGhost )
        {
            // create parameter list for ghost stabilization parameter for inclusion
            tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
            tParameterList( 2 )( tSPCounter ).set( "stabilization_name", "SPGPTemp1" );
            tParameterList( 2 )( tSPCounter ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::GHOST_DISPL ) );
            tParameterList( 2 )( tSPCounter ).set( "function_parameters", "0.01" );
            tParameterList( 2 )( tSPCounter ).set( "master_properties", "PropConductivity1,Material" );
            tSPCounter++;
        }

        //------------------------------------------------------------------------------
        // init IWG counter
        uint tIWGCounter = 0;

        // create IWG for inclusion - bulk diffusion
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGDiffusion1Bulk" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::SPATIALDIFF_BULK ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies", "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", "CMDiffusion1,Diffusion" );
        tParameterList( 3 )( tIWGCounter ).set( "master_properties", "PropHeatLoad1,Load" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tPhase1 );
        tIWGCounter++;

        // create IWG for Neumann boundary conditions
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGInletFlux" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::SPATIALDIFF_NEUMANN ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies", "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "master_properties", "PropSurfaceFlux,Neumann" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tFrontSurface );
        tIWGCounter++;

        // create IWG for Dirichlet boundary conditions
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWG2SurfaceTemp" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies", "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "master_properties", "PropImposedTemperature,Dirichlet" );
        tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", "CMDiffusion1,Diffusion" );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPNitscheTemp,DirichletNitsche" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tBackSurface );
        tIWGCounter++;

        if ( tUseGhost )
        {
            //             create IWG for 2 material - ghost
            //             tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            //             tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWGGP1Temp") ;
            //             tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::SPATIALDIFF_GHOST ) );
            //             tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               "TEMP") ;
            //             tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    "TEMP") ;
            //             tParameterList( 3 )( tIWGCounter ).set( "slave_dof_dependencies",     "TEMP") ;
            //             tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   "SPGPTemp1,GhostDispl") ;
            //             tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             tPhase1Ghost );
            //             tIWGCounter++;
        }

        // create IWG for time continuity
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGTimeContinuityTemp" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::TIME_CONTINUITY_DOF ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies", "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "master_properties",
                "PropWeightCurrent   ,WeightCurrent;"
                "PropWeightPrevious  ,WeightPrevious;"
                "PropInitialCondition,InitialCondition" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tPhase1 );
        tParameterList( 3 )( tIWGCounter ).set( "time_continuity", true );
        tIWGCounter++;

        //------------------------------------------------------------------------------
        // init IQI counter
        uint tIQICounter = 0;

        // Nodal Temperature IQI
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkTEMP" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::DOF ) );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "TEMP" );
        tParameterList( 4 )( tIQICounter ).set( "master_dof_dependencies", "TEMP" );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 0 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tPhase1 );
        tIQICounter++;

        // create parameter list for IQI 4
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkStrainEnergy" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::STRAIN_ENERGY ) );
        tParameterList( 4 )( tIQICounter ).set( "master_dof_dependencies", "TEMP" );
        tParameterList( 4 )( tIQICounter ).set( "master_constitutive_models", "CMDiffusion1,Elast" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tPhase1 );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkVolume" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::VOLUME ) );
        tParameterList( 4 )( tIQICounter ).set( "master_dof_dependencies", "TEMP" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tPhase1 );
        tIQICounter++;

        // Max Temperature IQI
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIMaxTemp" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::MAX_DOF ) );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "TEMP" );
        tParameterList( 4 )( tIQICounter ).set( "master_dof_dependencies", "TEMP" );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 0 );
        tParameterList( 4 )( tIQICounter ).set( "function_parameters", "1.0/2.0" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tTotalDomain );
        tIQICounter++;

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        tParameterList( 5 ).resize( 1 );
        tParameterList( 5 )( 0 ) = prm::create_computation_parameter_list();
    }

    void
    SOLParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 7 );
        for ( uint Ik = 0; Ik < 7; Ik++ )
        {
            tParameterlist( Ik ).resize( 1 );
        }

        tParameterlist( 0 )( 0 ) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL );

        tParameterlist( 1 )( 0 ) = moris::prm::create_linear_solver_parameter_list();

        tParameterlist( 2 )( 0 ) = moris::prm::create_nonlinear_algorithm_parameter_list();
        tParameterlist( 2 )( 0 ).set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        tParameterlist( 2 )( 0 ).set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
        tParameterlist( 2 )( 0 ).set( "NLA_max_iter", tNLA_max_iter );
        tParameterlist( 2 )( 0 ).set( "NLA_combined_res_jac_assembly", false );

        tParameterlist( 3 )( 0 ) = moris::prm::create_nonlinear_solver_parameter_list();
        tParameterlist( 3 )( 0 ).set( "NLA_DofTypes", "TEMP" );

        tParameterlist( 4 )( 0 ) = moris::prm::create_time_solver_algorithm_parameter_list();
        tParameterlist( 4 )( 0 ).set( "TSA_Num_Time_Steps", tTSA_Num_Time_Steps );
        tParameterlist( 4 )( 0 ).set( "TSA_Time_Frame", tTSA_Time_Frame );

        tParameterlist( 5 )( 0 ) = moris::prm::create_time_solver_parameter_list();
        tParameterlist( 5 )( 0 ).set( "TSA_DofTypes", "TEMP" );
        tParameterlist( 5 )( 0 ).set( "TSA_Initialize_Sol_Vec", "TEMP,0.0" );
        tParameterlist( 5 )( 0 ).set( "TSA_Output_Indices", "0" );
        tParameterlist( 5 )( 0 ).set( "TSA_Output_Crteria", "Output_Criterion" );
        tParameterlist( 5 )( 0 ).set( "TSA_time_level_per_type", "TEMP,2" );

        tParameterlist( 6 )( 0 ) = moris::prm::create_solver_warehouse_parameterlist();
        tParameterlist( 6 )( 0 ).set( "SOL_save_final_adjoint_vec_to_file", "Shape_Sensitivity_Transient.hdf5" );
    }

    void
    MSIParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_msi_parameter_list();
        tParameterlist( 0 )( 0 ).set( "order_adofs_by_host", false );
    }

    void
    VISParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_vis_parameter_list();
        tParameterlist( 0 )( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName ) );
        tParameterlist( 0 )( 0 ).set( "Mesh_Type", static_cast< uint >( vis::VIS_Mesh_Type::STANDARD ) );
        tParameterlist( 0 )( 0 ).set( "Set_Names", tPhase1 );
        tParameterlist( 0 )( 0 ).set( "Field_Names", "TEMP,MAX_DOF" );
        tParameterlist( 0 )( 0 ).set( "Field_Type", "NODAL,GLOBAL" );
        tParameterlist( 0 )( 0 ).set( "IQI_Names", "IQIBulkTEMP,IQIMaxTemp" );
        tParameterlist( 0 )( 0 ).set( "Save_Frequency", 1 );
    }

    void
    MORISGENERALParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
    }

    /* ------------------------------------------------------------------------ */
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
