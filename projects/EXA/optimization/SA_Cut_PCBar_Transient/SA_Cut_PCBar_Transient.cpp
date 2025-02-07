/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * SA_Cut_PCBar_Transient.cpp
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

// global variables
extern uint gInterpolationOrder;
extern bool gPrintReferenceValues;

#endif
//------------------------------------------------------------------------------
namespace moris
{

    //------------------------------------------------------------------------------
    //-------------------------------- QUICK SETUP ---------------------------------
    //------------------------------------------------------------------------------

    // include/exclude Nitsche Dirichlet boundary for debugging
    bool tHaveDirichlet = true;

    // include/exclude Ghost IWG for debugging
    bool tHaveGhost = true;

    // switch between running with phase change or without
    bool tIsPhaseChange = true;

    // Output Config --------------------------------------------------
    // set to true for vis output, set to false for sensitivity validation
    bool        tOutputCriterion = true;
    std::string tHDF5Path        = "SA_Cut_PCBar_Transient.hdf5";
    std::string tLibraryName     = "SA_Cut_PCBar_Transient.so";
    std::string tOutputFile      = "SA_Cut_PCBar_Transient.exo";

    // Geometry Parameters --------------------------------------------
    moris::real tXlength        = 0.1;
    moris::real tYlength        = 0.5;
    moris::real tPcmFinRatio    = 0.9;
    moris::real tDeltaRatio     = 0.05;
    moris::real tPcmFinRatioMin = tPcmFinRatio - tDeltaRatio;
    moris::real tPcmFinRatioMax = tPcmFinRatio + tDeltaRatio;
    std::string tInterfacePos   = std::to_string( tPcmFinRatio * tYlength );

    // Solver Configuration -------------------------------------------
    moris::sint tStep = 20;
    moris::real tTmax = 0.1;

    moris::real tNLARelResNormDrop      = 5.0e-07;
    moris::real tNLARelaxationParameter = 0.5;
    moris::sint tNLAMaxIter             = 25;

    // Sensitivity Analysis Parameters
    moris::real tFEMFdEpsilon = 1.0e-8;

    // material parameters --------------------------------------------

    // conductor material
    std::string tDensityFin      = "1.0";
    std::string tHeatCapFin      = "5.0";
    std::string tThermConductFin = "50.0";

    // heat storage material
    std::string tDensityPCM      = "1.0";
    std::string tHeatCapPCM      = "50.0";
    std::string tThermConductPCM = "5.0";
    std::string tLatentHeatPCM   = "500.0";
    std::string tPCTemp          = "20.0";
    std::string tPCTempRange     = "6.0";

    // initial & boundary conditions ----------------------------------
    std::string tHeatFlux          = "1000.0";
    moris::real tNBClengthFraction = 0.1;
    moris::real tInitialTemp       = 10.0;

    // IQI Configuration ----------------------------------------------
    std::string tMaxTempReference = "1.0";
    std::string tMaxTempExponent  = "2.0";

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
    std::string tPcmGhost = "ghost_p1";

    // HMR parameters -------------------------------------------------
    // Vector< uint > tNumElemsPerDim = "2, 16";
    Vector< uint > tNumElemsPerDim = { 1, 40 };
    Vector< real > tDomainDims     = { 0.16, 0.6 };
    std::string tDomainOffset   = "-0.0342356,-0.031345";

    //------------------------------------------------------------------------------
    //-------------------------------- FUNCTIONS -----------------------------------
    //------------------------------------------------------------------------------

    /* ------------------------------------------------------------------------ */
    // PROPERTY FUNCTIONS (incl. INITIAL & BOUNDARY CONDITIONS)
    /* ------------------------------------------------------------------------ */

    void
    Func_Initial_Temperature(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        aPropMatrix = { { tInitialTemp } };
    }

    // function for concentrating the heat load at the tip
    void
    Func_Neumann_BC(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        // get coordinates
        moris::Matrix< DDRMat > tPosition  = aFIManager->get_IG_geometry_interpolator()->valx();
        moris::real             tYPosition = tPosition( 1 );

        moris::real tHeatFluxOn = 0.0;

        if ( tYPosition >= tYlength - tNBClengthFraction * tYlength )
            tHeatFluxOn = 1.0;

        // assume that heat flux is proportional to surface stagnation pressure
        aPropMatrix = tHeatFluxOn * aParameters( 0 );
    }

    // Constant function for properties
    void
    Func_Const(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
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
            const Vector< real >&          aGeometryParameters )
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
        aParameterLists.set( "number_of_elements_per_dimension", tNumElemsPerDim );
        aParameterLists.set( "domain_dimensions", tDomainDims );
        aParameterLists.set( "domain_offset", tDomainOffset );
        aParameterLists.set( "lagrange_output_meshes", "0" );

        aParameterLists.set( "lagrange_orders", "1" );
        aParameterLists.set( "bspline_orders", "1" );
        aParameterLists.set( "lagrange_pattern", "0" );
        aParameterLists.set( "bspline_pattern", "0" );

        aParameterLists.set( "lagrange_to_bspline", "0" );

        aParameterLists.set( "truncate_bsplines", 1 );
        aParameterLists.set( "refinement_buffer", 0 );
        aParameterLists.set( "staircase_buffer", 0 );
        aParameterLists.set( "initial_refinement", "0" );
        aParameterLists.set( "initial_refinement_pattern", "0" );

        aParameterLists.set( "use_multigrid", 0 );
        aParameterLists.set( "severity_level", 0 );

        aParameterLists.set( "adaptive_refinement_level", 0 );
    }

    void
    OPTParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "is_optimization_problem", true );
        aParameterLists.set( "problem", "user_defined" );
        aParameterLists.set( "library", tLibraryName );

        aParameterLists( OPT::ALGORITHMS ).add_parameter_list( opt::Optimization_Algorithm_Type::SWEEP );
        aParameterLists.set( "print", true );
        aParameterLists.set( "hdf5_path", tHDF5Path );
        aParameterLists.set( "num_evaluations_per_adv", "1" );
        aParameterLists.set( "include_bounds", false );
        aParameterLists.set( "finite_difference_type", "all" );
        aParameterLists.set( "finite_difference_epsilons", "1E-6" );
    }

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
        aParameterLists.set( "verbose", true );
        aParameterLists.set( "print_enriched_ig_mesh", true );
        aParameterLists.set( "exodus_output_XTK_ig_mesh", true );
        aParameterLists.set( "high_to_low_dbl_side_sets", true );
    }

    void
    GENParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "IQI_types", "IQIMaxTemp", "IQIBulkVolume" );

        // Geometry parameter lists
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
        aParameterLists.set( "center_x", -1.0 );
        aParameterLists.set( "center_y", 0.0 );
        aParameterLists.set( "normal_x", 1.0 );
        aParameterLists.set( "normal_y", 0.0 );

        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
        aParameterLists.set( "center_x", 0.0 );
        aParameterLists.set( "center_y", tYlength * tPcmFinRatioMin, tYlength * tPcmFinRatioMin, tYlength * tPcmFinRatioMax );
        aParameterLists.set( "normal_x", 0.0 );
        aParameterLists.set( "normal_y", 1.0 );
    }

    void
    FEMParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.hack_for_legacy_fem();
        // create a cell of cell of parameter list for fem

        //------------------------------------------------------------------------------

        // Density of conductor material
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropFinDensity" );
        aParameterLists.set( "function_parameters", tDensityFin );
        aParameterLists.set( "value_function", "Func_Const" );

        // Thermal conductivity of conductor material
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropFinConductivity" );
        aParameterLists.set( "function_parameters", tThermConductFin );
        aParameterLists.set( "value_function", "Func_Const" );

        // Heat capacity of conductor material
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropFinHeatCapacity" );
        aParameterLists.set( "function_parameters", tHeatCapFin );
        aParameterLists.set( "value_function", "Func_Const" );

        // Density of storage material
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropPcmDensity" );
        aParameterLists.set( "function_parameters", tDensityPCM );
        aParameterLists.set( "value_function", "Func_Const" );

        // Thermal conductivity of storage material
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropPcmConductivity" );
        aParameterLists.set( "function_parameters", tThermConductPCM );
        aParameterLists.set( "value_function", "Func_Const" );

        // Heat capacity of storage material
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropPcmHeatCapacity" );
        aParameterLists.set( "function_parameters", tHeatCapPCM );
        aParameterLists.set( "value_function", "Func_Const" );

        // Latent heat capacity of storage material
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropLatentHeat" );
        aParameterLists.set( "function_parameters", tLatentHeatPCM );
        aParameterLists.set( "value_function", "Func_Const" );

        // Melt temperature of storage material
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropPCTemp" );
        aParameterLists.set( "function_parameters", tPCTemp );
        aParameterLists.set( "value_function", "Func_Const" );

        // Phase change function
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropPhaseState" );
        aParameterLists.set( "function_parameters", "2.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // Melting range of storage material
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropPCconst" );
        aParameterLists.set( "function_parameters", tPCTempRange );
        aParameterLists.set( "value_function", "Func_Const" );

        // Neumann BC
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropNeumannBC" );
        aParameterLists.set( "function_parameters", tHeatFlux );
        aParameterLists.set( "value_function", "Func_Const" );

        // Initial condition (Temperature)
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropInitialCondition" );
        aParameterLists.set( "value_function", "Func_Initial_Temperature" );

        // For Time solver
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropWeightCurrent" );
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "value_function", "Func_Const" );
                aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropWeightPrevious" );
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // Dummy latent heat for non-pc material
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDummyLatentHeat" );
        aParameterLists.set( "function_parameters", "0.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDummyPCTemp" );
        aParameterLists.set( "function_parameters", "10000.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        //------------------------------------------------------------------------------

        if ( tIsPhaseChange )
        {
            // constitutive model for phase change material
            aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
            aParameterLists.set( "constitutive_name", "CMDiffusionPcm" );
            aParameterLists.set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO_PC );
            aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
            aParameterLists.set( "properties",
                    "PropPcmConductivity,Conductivity;"
                    "PropPcmDensity,Density;"
                    "PropPcmHeatCapacity,HeatCapacity;"
                    "PropLatentHeat,LatentHeat;"
                    "PropPCTemp,PCTemp;"
                    "PropPhaseState,PhaseStateFunction;"
                    "PropPCconst,PhaseChangeConst" );
            }
        else    // no phase change
        {
            // constitutive model for thermal storage material
            aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
            aParameterLists.set( "constitutive_name", "CMDiffusionPcm" );
            aParameterLists.set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
            aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
            aParameterLists.set( "properties",
                    "PropPcmConductivity,Conductivity;"
                    "PropPcmDensity,Density;"
                    "PropPcmHeatCapacity,HeatCapacity" );
            }

        // constitutive model for thermal conductor material
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMDiffusionFin" );
        aParameterLists.set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists.set( "properties",
                "PropFinConductivity,Conductivity;"
                "PropFinDensity,Density;"
                "PropFinHeatCapacity,HeatCapacity" );

        //------------------------------------------------------------------------------

        // Ghost parameter for fin
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPGPTempFin" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        aParameterLists.set( "function_parameters", "0.01" );
        aParameterLists.set( "leader_properties", "PropFinConductivity,Material" );

        // Ghost parameter for PCM
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPGPTempPcm" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        aParameterLists.set( "function_parameters", "0.01" );
        aParameterLists.set( "leader_properties", "PropPcmConductivity,Material" );

        // GGLS parameter for fin
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPGGLSDiffusionFin" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::GGLS_DIFFUSION );
        aParameterLists.set( "leader_dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists.set( "leader_properties",
                "PropFinConductivity, Conductivity;"
                "PropFinDensity     , Density;"
                "PropFinHeatCapacity, HeatCapacity;"
                "PropDummyLatentHeat, LatentHeat;"
                "PropDummyPCTemp    , PCTemp;"
                "PropPhaseState     , PhaseStateFunction;"
                "PropPCconst        , PhaseChangeConst" );

        if ( tIsPhaseChange )
        {
            // GGLS parameter for phase change material
            aParameterLists( FEM::STABILIZATION ).add_parameter_list();
            aParameterLists.set( "stabilization_name", "SPGGLSDiffusionPcm" );
            aParameterLists.set( "stabilization_type", fem::Stabilization_Type::GGLS_DIFFUSION );
            aParameterLists.set( "leader_dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
            aParameterLists.set( "leader_properties",
                    "PropPcmConductivity, Conductivity;"
                    "PropPcmDensity     , Density;"
                    "PropPcmHeatCapacity, HeatCapacity;"
                    "PropLatentHeat     , LatentHeat;"
                    "PropPCTemp         , PCTemp;"
                    "PropPhaseState     , PhaseStateFunction;"
                    "PropPCconst        , PhaseChangeConst" );
            }
        else    // no phase change
        {
            // GGLS parameter for thermal storage material
            aParameterLists( FEM::STABILIZATION ).add_parameter_list();
            aParameterLists.set( "stabilization_name", "SPGGLSDiffusionPcm" );
            aParameterLists.set( "stabilization_type", fem::Stabilization_Type::GGLS_DIFFUSION );
            aParameterLists.set( "leader_dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
            aParameterLists.set( "leader_properties",
                    "PropPcmConductivity, Conductivity;"
                    "PropPcmDensity     , Density;"
                    "PropPcmHeatCapacity, HeatCapacity;"
                    "PropDummyLatentHeat, LatentHeat;"
                    "PropDummyPCTemp    , PCTemp;"
                    "PropPhaseState     , PhaseStateFunction;"
                    "PropPCconst        , PhaseChangeConst" );
            }

        // Dirichlet SP
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPInterfaceNitsche" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::NITSCHE_INTERFACE );
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "leader_properties", "PropFinConductivity,Material" );
        aParameterLists.set( "follower_properties", "PropPcmConductivity,Material" );

        //------------------------------------------------------------------------------
        // Bulk IWG for conductor material
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGDiffusionFinBulk" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusionFin,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPGGLSDiffusionFin,GGLSParam" );
        aParameterLists.set( "mesh_set_names", tFinBulk );

        // Bulk IWG for thermal storage material
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGDiffusionPcmBulk" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusionPcm,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPGGLSDiffusionPcm,GGLSParam" );
        aParameterLists.set( "mesh_set_names", tPcmBulk );

        if ( tHaveDirichlet )
        {
            // Interface Dirichlet BC
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "IWGInterfaceFinPcm" );
            aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE );
            aParameterLists.set( "dof_residual", "TEMP" );
            aParameterLists.set( "leader_dof_dependencies", "TEMP" );
            aParameterLists.set( "follower_dof_dependencies", "TEMP" );
            aParameterLists.set( "leader_constitutive_models", "CMDiffusionFin,Diffusion" );
            aParameterLists.set( "follower_constitutive_models", "CMDiffusionPcm,Diffusion" );
            aParameterLists.set( "stabilization_parameters", "SPInterfaceNitsche ,NitscheInterface" );
            aParameterLists.set( "mesh_set_names", tFinPcmInterface );
            }

        // Imposed heat flux
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGInletFlux" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_NEUMANN );
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "leader_properties", "PropNeumannBC,Neumann" );
        aParameterLists.set( "mesh_set_names", tFinNeumannInterface );

        if ( tHaveGhost )
        {
            // Fin Ghost
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "IWGFinGhost" );
            aParameterLists.set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists.set( "dof_residual", "TEMP" );
            aParameterLists.set( "leader_dof_dependencies", "TEMP" );
            aParameterLists.set( "follower_dof_dependencies", "TEMP" );
            aParameterLists.set( "stabilization_parameters", "SPGPTempFin,GhostSP" );
            aParameterLists.set( "mesh_set_names", tFinGhost );

            // PCM Ghost
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "IWGPcmGhost" );
            aParameterLists.set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists.set( "dof_residual", "TEMP" );
            aParameterLists.set( "leader_dof_dependencies", "TEMP" );
            aParameterLists.set( "follower_dof_dependencies", "TEMP" );
            aParameterLists.set( "stabilization_parameters", "SPGPTempPcm,GhostSP" );
            aParameterLists.set( "mesh_set_names", tPcmGhost );
            }

        // Time Continuity
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGTimeContinuityTemp" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::TIME_CONTINUITY_DOF );
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "leader_properties",
                "PropWeightCurrent,WeightCurrent;"
                "PropWeightPrevious  ,WeightPrevious;"
                "PropInitialCondition,InitialCondition" );
        aParameterLists.set( "mesh_set_names", tTotalDomain );
        aParameterLists.set( "time_continuity", true );

        //------------------------------------------------------------------------------
        // IQI - Nodal Temperature Field
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkTEMP" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "mesh_set_names", tTotalDomain );

        // Volume IQI - Total Volume
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkVolume" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::VOLUME );
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "mesh_set_names", tPcmBulk );

        // Max Temperature IQI
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIMaxTemp" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::MAX_DOF );
        aParameterLists.set( "dof_quantity", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "function_parameters", tMaxTempReference + "/" + tMaxTempExponent );
        aParameterLists.set( "mesh_set_names", tTotalDomain );

        // create computation  parameter list
        aParameterLists( FEM::COMPUTATION );
        aParameterLists.set( "finite_difference_scheme", fem::FDScheme_Type::POINT_3_CENTRAL );
        aParameterLists.set( "finite_difference_perturbation_size", tFEMFdEpsilon );
    }

    void
    SOLParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::BELOS_IMPL );
        aParameterLists.set( "preconditioners", "0" );

        aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list();

        // for forward analysis
        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "NLA_rel_res_norm_drop", tNLARelResNormDrop );
        aParameterLists.set( "NLA_relaxation_parameter", tNLARelaxationParameter );
        aParameterLists.set( "NLA_relaxation_strategy", sol::SolverRelaxationType::InvResNormAdaptive );
        aParameterLists.set( "NLA_relaxation_damping", 0.5 );
        aParameterLists.set( "NLA_max_iter", tNLAMaxIter );
        aParameterLists.set( "NLA_combined_res_jac_assembly", false );

        // for adjoint analysis
        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "NLA_rel_res_norm_drop", tNLARelResNormDrop );
        aParameterLists.set( "NLA_relaxation_parameter", 1.0 );
        aParameterLists.set( "NLA_max_iter", 1 );
        aParameterLists.set( "NLA_combined_res_jac_assembly", false );

        // for forward analysis
        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "NLA_DofTypes", "TEMP" );
        aParameterLists.set( "NLA_Nonlinear_solver_algorithms", "0" );

        // for adjoint analysis
        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "NLA_DofTypes", "TEMP" );
        aParameterLists.set( "NLA_Nonlinear_solver_algorithms", "1" );

        aParameterLists( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "TSA_Num_Time_Steps", tStep );
        aParameterLists.set( "TSA_Time_Frame", tTmax );
        aParameterLists.set( "TSA_Nonlinear_Solver", 0 );                // for forward analysis
        aParameterLists.set( "TSA_Nonlinear_Sensitivity_Solver", 1 );    // for adjoint analysis

        aParameterLists( SOL::TIME_SOLVERS ).add_parameter_list();
        aParameterLists.set( "TSA_DofTypes", "TEMP" );
        aParameterLists.set( "TSA_Initialize_Sol_Vec", "TEMP,0.0" );
        aParameterLists.set( "TSA_Output_Indices", "0" );
        aParameterLists.set( "TSA_Output_Criteria", "Output_Criterion" );
        aParameterLists.set( "TSA_time_level_per_type", "TEMP,2" );

        aParameterLists( SOL::PRECONDITIONERS ).add_parameter_list(  sol::PreconditionerType::IFPACK );
        aParameterLists.set( "ifpack_prec_type", "ILU" );
    }

    void
    MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
    }

    void
    VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFile ) );
        aParameterLists.set( "Mesh_Type", vis::VIS_Mesh_Type::STANDARD );
        aParameterLists.set( "Set_Names", tTotalDomain );
        aParameterLists.set( "Field_Names", "TEMP,MAX_DOF" );
        aParameterLists.set( "Field_Type", "NODAL,GLOBAL" );
        aParameterLists.set( "IQI_Names", "IQIBulkTEMP,IQIMaxTemp" );
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
