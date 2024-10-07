/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * SA_Cut_Bar_Static.cpp
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

#include "AztecOO.h"

#ifdef __cplusplus
extern "C" {

#endif
//------------------------------------------------------------------------------
namespace moris
{

    //------------------------------------------------------------------------------
    //-------------------------------- QUICK SETUP ---------------------------------
    //------------------------------------------------------------------------------

    // include/exclude Interface Nitsche Dirichlet boundary for debugging
    bool tHaveDirichlet = true;

    // include/exclude Ghost IWG for debugging
    bool tHaveGhost = true;

    // Output Config --------------------------------------------------
    // set to true for vis output, set to false for sensitivity validation
    bool        tOutputCriterion = true;
    std::string tHDF5Path        = "SA_Cut_Bar_Static.hdf5";
    std::string tLibraryName     = "SA_Cut_Bar_Static.so";
    std::string tOutputFile      = "SA_Cut_Bar_Static.exo";

    // Geometry Parameters --------------------------------------------
    moris::real tXlength        = 0.1;
    moris::real tYlength        = 0.5;
    moris::real tPcmFinRatio    = 0.9;
    moris::real tDeltaRatio     = 0.05;
    moris::real tPcmFinRatioMin = tPcmFinRatio - tDeltaRatio;
    moris::real tPcmFinRatioMax = tPcmFinRatio + tDeltaRatio;
    std::string tInterfacePos   = std::to_string( tPcmFinRatio * tYlength );

    // Solver Configuration -------------------------------------------
    moris::real tNLARelResNormDrop      = 1.0e-07;
    moris::real tNLARelaxationParameter = 1.0;
    moris::sint tNLAMaxIter             = 10;

    // material parameters --------------------------------------------

    // conductor material
    std::string tDensityFin      = "1.0";
    std::string tHeatCapFin      = "50.0";
    std::string tThermConductFin = "50.0";

    // heat storage material
    std::string tDensityPCM      = "1.0";
    std::string tHeatCapPCM      = "50.0";
    std::string tThermConductPCM = "1.0";
    std::string tLatentHeatPCM   = "1000.0";
    std::string tPCTemp          = "20.0";
    std::string tPCTempRange     = "6.0";

    // initial & boundary conditions ----------------------------------
    std::string tHeatFlux           = "1000.0";
    std::string tImposedTemperature = "10.0";

    // IQI Configuration ----------------------------------------------
    std::string tMaxTempReference = "10.0";
    std::string tMaxTempExponent  = "2.0";
    std::string tMaxTempShift     = "0.0";

    // Mesh sets ------------------------------------------------------

    // Bulk sets
    std::string tFinBulk     = "HMR_dummy_n_p3,HMR_dummy_c_p3";
    std::string tPcmBulk     = "HMR_dummy_n_p2,HMR_dummy_c_p2";
    std::string tTotalDomain = tFinBulk + "," + tPcmBulk;

    // Side sets
    std::string tFinPcmInterface       = "dbl_iside_p0_3_p1_2";
    std::string tFinNeumannInterface   = "SideSet_3_n_p3,SideSet_3_c_p3";
    std::string tPCMDirichletInterface = "SideSet_1_n_p2,SideSet_1_c_p2";

    // Ghost sets
    std::string tFinGhost = "ghost_p3";
    std::string tPcmGhost = "ghost_p2";

    // HMR parameters -------------------------------------------------
    std::string tNumElemsPerDim = "1, 40";
    std::string tDomainDims     = "0.16, 0.6";
    std::string tDomainOffset   = "-0.0342356,-0.031345";

    //------------------------------------------------------------------------------
    //-------------------------------- FUNCTIONS -----------------------------------
    //------------------------------------------------------------------------------

    /* ------------------------------------------------------------------------ */
    // PROPERTY FUNCTIONS (incl. INITIAL & BOUNDARY CONDITIONS)
    /* ------------------------------------------------------------------------ */

    // Constant function for properties
    void
    Func_Const(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    /* ------------------------------------------------------------------------ */
    // DUMMY FUNCTIONS
    /* ------------------------------------------------------------------------ */

    bool
    Output_Criterion( moris::tsa::Time_Solver* aTimeSolver )
    {
        return tOutputCriterion;
    }

    moris::Matrix< DDRMat > Func_Dummy_Sensitivity(
            const moris::Matrix< DDRMat >& aCoordinates,
            const Vector< real >&     aGeometryParameters )
    {
        moris::Matrix< DDRMat > aReturnValue = { { 0.0 } };
        return aReturnValue;
    }

    /* ------------------------------------------------------------------------ */
    // FOR SWEEP
    /* ------------------------------------------------------------------------ */

    moris::Matrix< moris::DDSMat >
    get_constraint_types()
    {
        Matrix< DDSMat > tConstraintTypes( 1, 1, 1 );
        return tConstraintTypes;
    }

    moris::Matrix< moris::DDRMat >
    compute_objectives( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        moris::Matrix< moris::DDRMat > tObjectives( 1, 1, aCriteria( 0 ) );
        return tObjectives;
    }

    moris::Matrix< moris::DDRMat >
    compute_constraints( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        moris::Matrix< moris::DDRMat > tConstraints( 1, 1, aCriteria( 1 ) );
        return tConstraints;
    }

    moris::Matrix< moris::DDRMat >
    compute_dobjective_dadv( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        moris::Matrix< moris::DDRMat > tDObjectiveDADV( 1, aADVs.size(), 0.0 );
        return tDObjectiveDADV;
    }

    moris::Matrix< moris::DDRMat >
    compute_dobjective_dcriteria( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        moris::Matrix< moris::DDRMat > tDObjectiveDCriteria( 1, 2, 0.0 );
        tDObjectiveDCriteria( 0 ) = 1.0;
        return tDObjectiveDCriteria;
    }

    moris::Matrix< moris::DDRMat >
    compute_dconstraint_dadv( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        moris::Matrix< moris::DDRMat > tDConstraintDADV( 1, aADVs.size(), 0.0 );
        return tDConstraintDADV;
    }

    moris::Matrix< moris::DDRMat >
    compute_dconstraint_dcriteria( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        moris::Matrix< moris::DDRMat > tDConstraintDCriteria( 1, 2, 0.0 );
        tDConstraintDCriteria( 1 ) = 1.0;
        return tDConstraintDCriteria;
    }

    /* ------------------------------------------------------------------------ */
    // PARAMETER LISTS
    /* ------------------------------------------------------------------------ */

    void
    HMRParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_hmr_parameter_list() );

        aParameterLists( 0 ).set( "number_of_elements_per_dimension", tNumElemsPerDim );
        aParameterLists( 0 ).set( "domain_dimensions", tDomainDims );
        aParameterLists( 0 ).set( "domain_offset", tDomainOffset );
        aParameterLists( 0 ).set( "domain_sidesets", "1,2,3,4" );
        aParameterLists( 0 ).set( "lagrange_output_meshes", "0" );

        aParameterLists( 0 ).set( "lagrange_orders", "1" );
        aParameterLists( 0 ).set( "bspline_orders", "1" );
        aParameterLists( 0 ).set( "lagrange_pattern", "0" );
        aParameterLists( 0 ).set( "bspline_pattern", "0" );

        aParameterLists( 0 ).set( "lagrange_to_bspline", "0" );

        aParameterLists( 0 ).set( "truncate_bsplines", 1 );
        aParameterLists( 0 ).set( "refinement_buffer", 0 );
        aParameterLists( 0 ).set( "staircase_buffer", 0 );
        aParameterLists( 0 ).set( "initial_refinement", "0" );
        aParameterLists( 0 ).set( "initial_refinement_pattern", "0" );

        aParameterLists( 0 ).set( "use_multigrid", 0 );
        aParameterLists( 0 ).set( "severity_level", 0 );

        aParameterLists( 0 ).set( "adaptive_refinement_level", 0 );
    }

    void
    OPTParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( moris::prm::create_opt_problem_parameter_list() );
        aParameterLists( 0 ).set( "is_optimization_problem", true );
        aParameterLists( 0 ).set( "problem", "user_defined" );
        aParameterLists( 0 ).set( "library", tLibraryName );

        aParameterLists( 2 ).add_parameter_list( moris::prm::create_sweep_parameter_list() );
        aParameterLists( 2 ).set( "print", true );
        aParameterLists( 2 ).set( "hdf5_path", tHDF5Path );
        aParameterLists( 2 ).set( "num_evaluations_per_adv", "1" );
        aParameterLists( 2 ).set( "include_bounds", false );
        aParameterLists( 2 ).set( "finite_difference_type", "all" );
    }

    void
    XTKParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_xtk_parameter_list() );
        aParameterLists( 0 ).set( "decompose", true );
        aParameterLists( 0 ).set( "decomposition_type", "conformal" );
        aParameterLists( 0 ).set( "enrich", true );
        aParameterLists( 0 ).set( "basis_rank", "bspline" );
        aParameterLists( 0 ).set( "enrich_mesh_indices", "0" );
        aParameterLists( 0 ).set( "ghost_stab", true );
        aParameterLists( 0 ).set( "multigrid", false );
        aParameterLists( 0 ).set( "verbose", true );
        aParameterLists( 0 ).set( "print_enriched_ig_mesh", true );
        aParameterLists( 0 ).set( "exodus_output_XTK_ig_mesh", true );
        aParameterLists( 0 ).set( "high_to_low_dbl_side_sets", true );
    }

    void
    GENParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_gen_parameter_list() );
        aParameterLists( 0 ).set( "IQI_types", "IQIMaxTemp", "IQIBulkVolume" );

        // Geometry parameter lists
        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
        aParameterLists( 1 ).set( "center_x", -1.0 );
        aParameterLists( 1 ).set( "center_y", 0.0 );
        aParameterLists( 1 ).set( "normal_x", 1.0 );
        aParameterLists( 1 ).set( "normal_y", 0.0 );

        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
        aParameterLists( 1 ).set( "center_x", 0.0 );
        aParameterLists( 1 ).set( "center_y", tYlength * tPcmFinRatioMin, tYlength * tPcmFinRatioMin, tYlength * tPcmFinRatioMax );
        aParameterLists( 1 ).set( "normal_x", 0.0 );
        aParameterLists( 1 ).set( "normal_y", 1.0 );
    }

    void
    FEMParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.hack_for_legacy_fem();
        // create a cell of cell of parameter list for fem

        //------------------------------------------------------------------------------

        // Density of conductor material
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropFinDensity" );
        aParameterLists( 0 ).set( "function_parameters", tDensityFin );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // Thermal conductivity of conductor material
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropFinConductivity" );
        aParameterLists( 0 ).set( "function_parameters", tThermConductFin );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // Heat capacity of conductor material
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropFinHeatCapacity" );
        aParameterLists( 0 ).set( "function_parameters", tHeatCapFin );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // Density of storage material
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropPcmDensity" );
        aParameterLists( 0 ).set( "function_parameters", tDensityPCM );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // Thermal conductivity of storage material
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropPcmConductivity" );
        aParameterLists( 0 ).set( "function_parameters", tThermConductPCM );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // Heat capacity of storage material
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropPcmHeatCapacity" );
        aParameterLists( 0 ).set( "function_parameters", tHeatCapPCM );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // Latent heat capacity of storage material
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropLatentHeat" );
        aParameterLists( 0 ).set( "function_parameters", tLatentHeatPCM );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // Melt temperature of storage material
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropPCTemp" );
        aParameterLists( 0 ).set( "function_parameters", tPCTemp );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // Phase change function
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropPhaseState" );
        aParameterLists( 0 ).set( "function_parameters", "2.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // Melting range of storage material
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropPCconst" );
        aParameterLists( 0 ).set( "function_parameters", tPCTempRange );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // Neumann BC
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropNeumannBC" );
        aParameterLists( 0 ).set( "function_parameters", tHeatFlux );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // Dirichlet BC
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropImposedTemperature" );
        aParameterLists( 0 ).set( "function_parameters", tImposedTemperature );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // Dummy latent heat for non-pc material
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropDummyLatentHeat" );
        aParameterLists( 0 ).set( "function_parameters", "0.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropDummyPCTemp" );
        aParameterLists( 0 ).set( "function_parameters", "10000.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        //------------------------------------------------------------------------------

        // constitutive model for thermal storage material
        aParameterLists( 1 ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists( 1 ).set( "constitutive_name", "CMDiffusionPcm" );
        aParameterLists( 1 ).set( "constitutive_type",  fem::Constitutive_Type::DIFF_LIN_ISO ) ;
        aParameterLists( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists( 1 ).set( "properties",
                "PropPcmConductivity,Conductivity;"
                "PropPcmDensity,Density;"
                "PropPcmHeatCapacity,HeatCapacity" );

        // constitutive model for thermal conductor material
        aParameterLists( 1 ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists( 1 ).set( "constitutive_name", "CMDiffusionFin" );
        aParameterLists( 1 ).set( "constitutive_type",  fem::Constitutive_Type::DIFF_LIN_ISO ) ;
        aParameterLists( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists( 1 ).set( "properties",
                "PropFinConductivity,Conductivity;"
                "PropFinDensity,Density;"
                "PropFinHeatCapacity,HeatCapacity" );

        //------------------------------------------------------------------------------

        // Ghost parameter for fin
        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPGPTempFin" );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::GHOST_DISPL ) ;
        aParameterLists( 2 ).set( "function_parameters", "0.01" );
        aParameterLists( 2 ).set( "leader_properties", "PropFinConductivity,Material" );

        // Ghost parameter for PCM
        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPGPTempPcm" );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::GHOST_DISPL ) ;
        aParameterLists( 2 ).set( "function_parameters", "0.01" );
        aParameterLists( 2 ).set( "leader_properties", "PropPcmConductivity,Material" );

        // GGLS parameter for fin
        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPGGLSDiffusionFin" );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::GGLS_DIFFUSION ) ;
        aParameterLists( 2 ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists( 2 ).set( "leader_properties",
                "PropFinConductivity, Conductivity;"
                "PropFinDensity     , Density;"
                "PropFinHeatCapacity, HeatCapacity;"
                "PropDummyLatentHeat, LatentHeat;"
                "PropDummyPCTemp    , PCTemp;"
                "PropPhaseState     , PhaseStateFunction;"
                "PropPCconst        , PhaseChangeConst" );

        // GGLS parameter for thermal storage material
        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPGGLSDiffusionPcm" );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::GGLS_DIFFUSION ) ;
        aParameterLists( 2 ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists( 2 ).set( "leader_properties",
                "PropPcmConductivity, Conductivity;"
                "PropPcmDensity     , Density;"
                "PropPcmHeatCapacity, HeatCapacity;"
                "PropDummyLatentHeat, LatentHeat;"
                "PropDummyPCTemp    , PCTemp;"
                "PropPhaseState     , PhaseStateFunction;"
                "PropPCconst        , PhaseChangeConst" );

        // Interface Dirichlet SP
        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPInterfaceNitsche" );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::NITSCHE_INTERFACE ) ;
        aParameterLists( 2 ).set( "function_parameters", "100.0" );
        aParameterLists( 2 ).set( "leader_properties", "PropFinConductivity,Material" );
        aParameterLists( 2 ).set( "follower_properties", "PropPcmConductivity,Material" );

        // Boundary Dirichlet SP
        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPNitscheTemp" );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::DIRICHLET_NITSCHE ) ;
        aParameterLists( 2 ).set( "function_parameters", "10.0" );
        aParameterLists( 2 ).set( "leader_properties", "PropPcmConductivity,Material" );

        //------------------------------------------------------------------------------
        // Bulk IWG for conductor material
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGDiffusionFinBulk" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_BULK ) ;
        aParameterLists( 3 ).set( "dof_residual", "TEMP" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMDiffusionFin,Diffusion" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPGGLSDiffusionFin,GGLSParam" );
        aParameterLists( 3 ).set( "mesh_set_names", tFinBulk );

        // Bulk IWG for thermal storage material
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGDiffusionPcmBulk" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_BULK ) ;
        aParameterLists( 3 ).set( "dof_residual", "TEMP" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMDiffusionPcm,Diffusion" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPGGLSDiffusionPcm,GGLSParam" );
        aParameterLists( 3 ).set( "mesh_set_names", tPcmBulk );

        // Interface Dirichlet BC
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGInterfaceFinPcm" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE ) ;
        aParameterLists( 3 ).set( "dof_residual", "TEMP" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 3 ).set( "follower_dof_dependencies", "TEMP" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMDiffusionFin,Diffusion" );
        aParameterLists( 3 ).set( "follower_constitutive_models", "CMDiffusionPcm,Diffusion" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPInterfaceNitsche ,NitscheInterface" );
        aParameterLists( 3 ).set( "mesh_set_names", tFinPcmInterface );

        // Imposed heat flux
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGInletFlux" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_NEUMANN ) ;
        aParameterLists( 3 ).set( "dof_residual", "TEMP" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 3 ).set( "leader_properties", "PropNeumannBC,Neumann" );
        aParameterLists( 3 ).set( "mesh_set_names", tFinNeumannInterface );

        // Imposed temperature
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGOutletTemp" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE ) ;
        aParameterLists( 3 ).set( "dof_residual", "TEMP" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 3 ).set( "leader_properties", "PropImposedTemperature,Dirichlet" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMDiffusionPcm,Diffusion" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPNitscheTemp,DirichletNitsche" );
        aParameterLists( 3 ).set( "mesh_set_names", tPCMDirichletInterface );

        if ( tHaveGhost )
        {
            // Fin Ghost
            aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists( 3 ).set( "IWG_name", "IWGFinGhost" );
            aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
            aParameterLists( 3 ).set( "dof_residual", "TEMP" );
            aParameterLists( 3 ).set( "leader_dof_dependencies", "TEMP" );
            aParameterLists( 3 ).set( "follower_dof_dependencies", "TEMP" );
            aParameterLists( 3 ).set( "stabilization_parameters", "SPGPTempFin,GhostSP" );
            aParameterLists( 3 ).set( "mesh_set_names", tFinGhost );

            // PCM Ghost
            aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists( 3 ).set( "IWG_name", "IWGPcmGhost" );
            aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
            aParameterLists( 3 ).set( "dof_residual", "TEMP" );
            aParameterLists( 3 ).set( "leader_dof_dependencies", "TEMP" );
            aParameterLists( 3 ).set( "follower_dof_dependencies", "TEMP" );
            aParameterLists( 3 ).set( "stabilization_parameters", "SPGPTempPcm,GhostSP" );
            aParameterLists( 3 ).set( "mesh_set_names", tPcmGhost );
            }

        //------------------------------------------------------------------------------
        // IQI - Nodal Temperature Field
        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkTEMP" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists( 4 ).set( "dof_quantity", "TEMP" );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 4 ).set( "vectorial_field_index", 0 );
        aParameterLists( 4 ).set( "mesh_set_names", tTotalDomain );

        // Volume IQI - Total Volume
        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkVolume" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::VOLUME ) ;
        aParameterLists( 4 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 4 ).set( "mesh_set_names", tPcmBulk );

        // Max Temperature IQI
        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIMaxTemp" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::MAX_DOF ) ;
        aParameterLists( 4 ).set( "dof_quantity", "TEMP" );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 4 ).set( "vectorial_field_index", 0 );
        aParameterLists( 4 ).set( "function_parameters", tMaxTempReference + "/" + tMaxTempExponent + "/" + tMaxTempShift );
        aParameterLists( 4 ).set( "mesh_set_names", tTotalDomain );

        // create computation  parameter list
        aParameterLists( 5 ).add_parameter_list( prm::create_computation_parameter_list() );
    }

    void
    SOLParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::BELOS_IMPL ) );
        aParameterLists( 0 ).set( "preconditioners", "0" );

        aParameterLists( 1 ).add_parameter_list( moris::prm::create_linear_solver_parameter_list() );

        aParameterLists( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
        aParameterLists( 2 ).set( "NLA_rel_res_norm_drop", tNLARelResNormDrop );
        aParameterLists( 2 ).set( "NLA_relaxation_parameter", tNLARelaxationParameter );
        aParameterLists( 2 ).set( "NLA_max_iter", tNLAMaxIter );
        aParameterLists( 2 ).set( "NLA_combined_res_jac_assembly", false );

        aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
        aParameterLists( 3 ).set( "NLA_DofTypes", "TEMP" );

        aParameterLists( 4 ).add_parameter_list( moris::prm::create_time_solver_algorithm_parameter_list() );

        aParameterLists( 5 ).add_parameter_list( moris::prm::create_time_solver_parameter_list() );
        aParameterLists( 5 ).set( "TSA_DofTypes", "TEMP" );
        aParameterLists( 5 ).set( "TSA_Initialize_Sol_Vec", "TEMP,0.0" );
        aParameterLists( 5 ).set( "TSA_Output_Indices", "0" );
        aParameterLists( 5 ).set( "TSA_Output_Criteria", "Output_Criterion" );
        aParameterLists( 5 ).set( "TSA_time_level_per_type", "TEMP,1" );

        aParameterLists( 6 ).add_parameter_list( moris::prm::create_solver_warehouse_parameterlist() );

        aParameterLists( 7 ).add_parameter_list( moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::IFPACK ) );
        aParameterLists( 7 ).set( "ifpack_prec_type", "ILU" );
    }

    void
    MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_msi_parameter_list() );
    }

    void
    VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_vis_parameter_list() );
        aParameterLists( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFile ) );
        aParameterLists( 0 ).set( "Mesh_Type",  vis::VIS_Mesh_Type::STANDARD ) ;

        aParameterLists( 0 ).set( "Set_Names", tTotalDomain );
        aParameterLists( 0 ).set( "Field_Names", "TEMP,MAX_DOF" );
        aParameterLists( 0 ).set( "Field_Type", "NODAL,GLOBAL" );
        aParameterLists( 0 ).set( "IQI_Names", "IQIBulkTEMP,IQIMaxTemp" );
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
