/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Solver_Examples.cpp
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

// $MORISROOT/share/scripts/create_shared_object.sh . build_dbg Solver_Examples
// $MORISROOT/build_dbg/projects/mains/moris ./Solver_Examples.so

#ifdef __cplusplus
extern "C" {
#endif

// global variables from unit test
extern bool gHaveStaggeredFA;
extern bool gHaveStaggeredSA;
extern bool gUseMixedTimeElements;
extern bool gUseBelosWithILUT;

//------------------------------------------------------------------------------
namespace moris
{
    // bool gHaveStaggeredFA = true;
    // bool gHaveStaggeredSA = true;
    // bool gUseMixedTimeElements = false;
    // bool gUseBelosWithILUT = false;

    // file name
    std::string tName = "Solver_Examples";

    /* ------------------------------------------------------------------------ */
    // Configuration - do not change - parameters are changed later based on input above

    std::string sStructDofTypes = "UX,UY";
    std::string sAllDofTypes    = "UX,UY;TEMP";

    // time levels for structure dofs
    std::string sTLSD = "2";

    /* ------------------------------------------------------------------------ */
    // Set Names

    // Bulk sets
    std::string tPhase1      = "0";
    std::string tPhase2      = "1";
    std::string tMat1Set     = "HMR_dummy_n_p" + tPhase1 + ",HMR_dummy_c_p" + tPhase1;
    std::string tMat2Set     = "HMR_dummy_n_p" + tPhase2 + ",HMR_dummy_c_p" + tPhase2;
    std::string tTotalDomain = tMat1Set + "," + tMat2Set;

    /* ------------------------------------------------------------------------ */
    // material parameters, kg is scaled with a factor 1e-6

    // Material 1
    std::string tDensity1      = "1.0";
    std::string tHeatCapacity1 = "5.0";
    std::string tConductivity1 = "50.0";

    std::string tYoungsModulus1    = "100.0";
    std::string tPoissonRatio1     = "0.31";
    std::string tThermalExpansion1 = "1.0e-2";

    // Material 2 (PCM)
    std::string tDensity2      = "1.0";
    std::string tHeatCapacity2 = "50.0";
    std::string tConductivity2 = "5.0";
    std::string tLatentHeat2   = "500.0";
    std::string tPCTemp2       = "20.0";
    std::string tPCConst2      = "6.0";

    std::string tYoungsModulus2    = "10.0";
    std::string tPoissonRatio2     = "0.31";
    std::string tThermalExpansion2 = "1.0e-3";

    /* ------------------------------------------------------------------------ */
    // boundary conditions

    // Pressure
    std::string tLoad = "1.0";

    // Initial Temperature
    moris::real tInitialTemp   = 10.0;
    std::string tReferenceTemp = "1.0";

    // Heat flux
    std::string tHeatLoad = "1000.0";

    /* ------------------------------------------------------------------------ */
    // Solver config

    moris::real tNLA_rel_res_norm_drop    = 5.0e-06;
    moris::real tNLA_relaxation_parameter = 1.0;
    int         tNLA_max_iter             = 15;

    int         tTSA_Num_Time_Steps = 2;
    moris::real tTSA_Time_Frame     = 0.1;

    /* ------------------------------------------------------------------------ */
    // Output Config

    std::string tOutputFileName  = tName + ".exo";
    std::string tLibraryName     = tName + ".so";
    std::string tHDF5Path        = tName + ".hdf5";
    std::string tGENOutputFile   = "GEN_" + tName + ".exo";
    bool        tOutputCriterion = true;

    //------------------------------------------------------------------------------
    //-------------------------------- FUNCTIONS -----------------------------------
    //------------------------------------------------------------------------------

    /* ------------------------------------------------------------------------ */
    // GEOMETRY (LEVEL-SET) FUNCTIONS
    /* ------------------------------------------------------------------------ */

    // Interface
    //    moris::real Interface(
    //            const moris::Matrix< DDRMat >     & aCoordinates,
    //            const Vector< real > & aGeometryParameters )
    //    {
    //        // compute level set value
    //        moris::real aReturnValue = ( aCoordinates( 0 ) - *aGeometryParameters( 0 ) - 0.4 );
    //
    //        // clean return value to return non-zero value
    //        return std::abs(aReturnValue) < 1e-8 ? 1e-8 : aReturnValue;
    //    }

    /* ------------------------------------------------------------------------ */
    // PROPERTY FUNCTIONS (incl. INITIAL & BOUNDARY CONDITIONS)
    /* ------------------------------------------------------------------------ */

    // Constant function for properties
    void
    Func_Const(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    // initial temperature
    void
    Func_Initial_Condition(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        aPropMatrix = { { tInitialTemp } };
    }

    /* ------------------------------------------------------------------------ */
    // DUMMY FUNCTIONS
    /* ------------------------------------------------------------------------ */

    // Output criterion for VIS mesh
    bool
    Output_Criterion( moris::tsa::Time_Solver* aTimeSolver )
    {
        return tOutputCriterion;
    }

    // Dummy function for unused sensitivities if needed
    moris::Matrix< DDRMat > Func_Dummy_Sensitivity(
            const moris::Matrix< DDRMat >& aCoordinates,
            const Vector< real >&          aGeometryParameters )
    {
        moris::Matrix< DDRMat > aReturnValue = { { 0.0 } };
        return aReturnValue;
    }

    /* ------------------------------------------------------------------------ */
    // FOR OPTIMIZATION
    /* ------------------------------------------------------------------------ */

    moris::Matrix< moris::DDSMat >
    get_constraint_types()
    {
        Matrix< DDSMat > tConstraintTypes( 1, 1, 1 );
        return tConstraintTypes;
    }

    moris::Matrix< moris::DDRMat >
    compute_objectives(
            const Vector< real >& aADVs,
            const Vector< real >& aCriteria )
    {
        moris::Matrix< moris::DDRMat > tObjectives( 1, 1, 0.0 );

        tObjectives( 0, 0 ) = aCriteria( 0 ) + aCriteria( 1 ) + aCriteria( 2 ) + aCriteria( 3 );

        std::cout << "% --------------------------------- % \n"
                  << std::flush;
        std::cout << "Max Temp Value = " << aCriteria( 0 ) + aCriteria( 1 ) << " \n"
                  << std::flush;
        std::cout << "Strain Energy = " << aCriteria( 2 ) + aCriteria( 3 ) << " \n"
                  << std::flush;
        std::cout << "Volume 1 = " << aCriteria( 4 ) << " \n"
                  << std::flush;
        std::cout << "% --------------------------------- % \n"
                  << std::flush;
        std::cout << "Objective = " << tObjectives( 0, 0 ) << " \n"
                  << std::flush;
        std::cout << "% --------------------------------- % \n"
                  << std::flush;

        return tObjectives;
    }

    moris::Matrix< moris::DDRMat >
    compute_constraints(
            const Vector< real >& aADVs,
            const Vector< real >& aCriteria )
    {
        moris::Matrix< moris::DDRMat > tConstraints( 1, 1, 0.0 );

        tConstraints( 0, 0 ) = aCriteria( 4 );

        std::cout << "Constraint = " << tConstraints( 0, 0 ) << " \n"
                  << std::flush;
        std::cout << "% --------------------------------- % \n"
                  << std::flush;

        return tConstraints;
    }

    moris::Matrix< moris::DDRMat >
    compute_dobjective_dadv(
            const Vector< real >& aADVs,
            const Vector< real >& aCriteria )
    {
        moris::Matrix< moris::DDRMat > tDObjectiveDADV( 1, aADVs.size(), 0.0 );
        return tDObjectiveDADV;
    }

    moris::Matrix< moris::DDRMat >
    compute_dobjective_dcriteria(
            const Vector< real >& aADVs,
            const Vector< real >& aCriteria )
    {
        moris::Matrix< moris::DDRMat > tDObjectiveDCriteria;

        tDObjectiveDCriteria.set_size( 1, 5, 1.0 );
        tDObjectiveDCriteria( 4 ) = 0.0;

        return tDObjectiveDCriteria;
    }

    moris::Matrix< moris::DDRMat >
    compute_dconstraint_dadv(
            const Vector< real >& aADVs,
            const Vector< real >& aCriteria )
    {
        moris::Matrix< moris::DDRMat > tDConstraintDADV( 1, aADVs.size(), 0.0 );
        return tDConstraintDADV;
    }

    moris::Matrix< moris::DDRMat >
    compute_dconstraint_dcriteria(
            const Vector< real >& aADVs,
            const Vector< real >& aCriteria )
    {
        moris::Matrix< moris::DDRMat > tDConstraintDCriteria;
        tDConstraintDCriteria.set_size( 1, 5, 0.0 );
        tDConstraintDCriteria( 4 ) = 1.0;

        return tDConstraintDCriteria;
    }

    /* ------------------------------------------------------------------------ */
    // PARAMETER LISTS
    /* ------------------------------------------------------------------------ */

    void
    OPTParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( moris::prm::create_opt_problem_parameter_list() );
        aParameterLists( 0 ).set( "is_optimization_problem", true );
        aParameterLists( 0 ).set( "problem", "user_defined" );
        aParameterLists( 0 ).set( "library", tLibraryName );

        aParameterLists( 2 ).add_parameter_list( moris::prm::create_gcmma_parameter_list() );
        aParameterLists( 2 ).set( "max_its", 1 );
        aParameterLists( 2 ).set( "step_size", 0.2 );
    }

    /* ------------------------------------------------------------------------ */

    void
    HMRParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_hmr_parameter_list() );

        aParameterLists( 0 ).set( "number_of_elements_per_dimension", "32, 8" );
        aParameterLists( 0 ).set( "domain_dimensions", "0.8,0.2" );
        aParameterLists( 0 ).set( "domain_offset", "0.0,0.0" );
        aParameterLists( 0 ).set( "domain_sidesets", "1,2,3,4" );
        aParameterLists( 0 ).set( "lagrange_output_meshes", "0" );

        aParameterLists( 0 ).set( "lagrange_orders", "1" );
        aParameterLists( 0 ).set( "lagrange_pattern", "0" );
        aParameterLists( 0 ).set( "bspline_orders", "1" );
        aParameterLists( 0 ).set( "bspline_pattern", "0" );

        aParameterLists( 0 ).set( "truncate_bsplines", 1 );
        aParameterLists( 0 ).set( "use_number_aura", 1 );
        aParameterLists( 0 ).set( "use_multigrid", 0 );
        aParameterLists( 0 ).set( "severity_level", 0 );
    }

    /* ------------------------------------------------------------------------ */

    void
    XTKParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_xtk_parameter_list() );
        aParameterLists( 0 ).set( "decompose", true );
        aParameterLists( 0 ).set( "decomposition_type", std::string( "conformal" ) );
        aParameterLists( 0 ).set( "enrich", true );
        aParameterLists( 0 ).set( "basis_rank", std::string( "bspline" ) );
        aParameterLists( 0 ).set( "enrich_mesh_indices", std::string( "0" ) );
        aParameterLists( 0 ).set( "ghost_stab", true );
        aParameterLists( 0 ).set( "multigrid", false );
        aParameterLists( 0 ).set( "verbose", true );
        aParameterLists( 0 ).set( "print_enriched_ig_mesh", true );
        aParameterLists( 0 ).set( "exodus_output_XTK_ig_mesh", true );
        aParameterLists( 0 ).set( "high_to_low_dbl_side_sets", true );
        aParameterLists( 0 ).set( "output_path", "./" );
        aParameterLists( 0 ).set( "keep_all_opt_iters", true );
    }

    /* ------------------------------------------------------------------------ */

    void
    GENParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_gen_parameter_list() );
        aParameterLists( 0 ).set( "IQI_types", "IQIMaxTemp_1", "IQIMaxTemp_2", "IQIStrainEnergy_1", "IQIStrainEnergy_2", "IQIVolume_1" );
        aParameterLists( 0 ).set( "output_mesh_file", tGENOutputFile );

        // Interface
        //        aParameterLists( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        //        aParameterLists( 1 ).set( "field_function_name", "Interface" );
        //
        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
        aParameterLists( 1 ).set( "center_x", 0.1, 0.4, 0.7 );
        aParameterLists( 1 ).set( "center_y", 0.0 );
        aParameterLists( 1 ).set( "normal_x", 1.0 );
        aParameterLists( 1 ).set( "normal_y", 0.0 );
        aParameterLists( 1 ).set( "isocontour_tolerance", 1e-12 );    // Interface tolerance based on geometry value
    }

    /* ------------------------------------------------------------------------ */

    void
    FEMParameterList( Module_Parameter_Lists& aParameterLists )
    {
        std::string sStructDofTypes = "UX,UY";
        std::string sAllDofTypes    = "UX,UY;TEMP";

        // create a cell of cell of parameter list for fem
        uint tPropIndex = 0;
        uint tCMIndex   = 1;
        uint tSPIndex   = 2;
        uint tIWGIndex  = 3;
        uint tIQIIndex  = 4;
        uint tFEMIndex  = 5;
        // uint tFieldIndex   = 6;
        uint tPhaseIndex = 7;

        ////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////

        aParameterLists( tPhaseIndex ).add_parameter_list( prm::create_phase_parameter_list() );
        aParameterLists( tPhaseIndex ).set( "phase_name", "Phase1" );
        aParameterLists( tPhaseIndex ).set( "phase_indices", tPhase1 );

        aParameterLists( tPhaseIndex ).add_parameter_list( prm::create_phase_parameter_list() );
        aParameterLists( tPhaseIndex ).set( "phase_name", "Phase2" );
        aParameterLists( tPhaseIndex ).set( "phase_indices", tPhase2 );

        ////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////

        //------------------------------------------------------------------------------
        // MATERIAL PARAMETERS - STRUCTURE (ni-w-alloy?)
        //------------------------------------------------------------------------------

        // Density Shell
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropDensity1" );
        aParameterLists( tPropIndex ).set( "function_parameters", tDensity1 );

        // Heat Capacity Shell
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropHeatCapacity1" );
        aParameterLists( tPropIndex ).set( "function_parameters", tHeatCapacity1 );

        // Conductivity Shell
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropConductivity1" );
        aParameterLists( tPropIndex ).set( "function_parameters", tConductivity1 );

        // Youngs Modulus Shell
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropYoungsModulus1" );
        aParameterLists( tPropIndex ).set( "function_parameters", tYoungsModulus1 );

        // Poisson Ratio Shell
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropPoissonRatio1" );
        aParameterLists( tPropIndex ).set( "function_parameters", tPoissonRatio1 );

        // CTE for Shell
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropThermalExpansion1" );
        aParameterLists( tPropIndex ).set( "function_parameters", tThermalExpansion1 );

        //------------------------------------------------------------------------------
        // MATERIAL PARAMETERS - PCM (Al-Cu-Si-alloy?)
        //------------------------------------------------------------------------------

        // Density of PCM
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropDensity2" );
        aParameterLists( tPropIndex ).set( "function_parameters", tDensity2 );

        // Heat Capacity of PCM
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropHeatCapacity2" );
        aParameterLists( tPropIndex ).set( "function_parameters", tHeatCapacity2 );

        // Conductivity of PCM
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropConductivity2" );
        aParameterLists( tPropIndex ).set( "function_parameters", tConductivity2 );

        // Latent Heat of PCM
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropLatentHeat2" );
        aParameterLists( tPropIndex ).set( "function_parameters", tLatentHeat2 );

        // Phase Change Temperature of PCM
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropPCTemp2" );
        aParameterLists( tPropIndex ).set( "function_parameters", tPCTemp2 );

        // Phase Change Temperature Range of PCM
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropPCconst2" );
        aParameterLists( tPropIndex ).set( "function_parameters", tPCConst2 );

        // Cubic Phase State Function for phase change model
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropPhaseState2" );
        aParameterLists( tPropIndex ).set( "function_parameters", "2.0" );

        // Youngs Modulus for PCM
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropYoungsModulus2" );
        aParameterLists( tPropIndex ).set( "function_parameters", tYoungsModulus2 );

        // Poisson Ratio for PCM
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropPoissonRatio2" );
        aParameterLists( tPropIndex ).set( "function_parameters", tPoissonRatio2 );

        // CTE for PCM
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropThermalExpansion2" );
        aParameterLists( tPropIndex ).set( "function_parameters", tThermalExpansion2 );

        //------------------------------------------------------------------------------
        // OTHER MATERIAL PARAMETERS
        //------------------------------------------------------------------------------

        // Dummy latent heat for non-pc material
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropLatentHeat_Dummy" );
        aParameterLists( tPropIndex ).set( "function_parameters", "0.0" );

        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropPCTemp_Dummy" );
        aParameterLists( tPropIndex ).set( "function_parameters", "10000.0" );

        //------------------------------------------------------------------------------
        // BOUNDARY CONDITIONS
        //------------------------------------------------------------------------------

        // reference temperature for thermal expansion
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropReferenceTemp" );
        aParameterLists( tPropIndex ).set( "function_parameters", tReferenceTemp );

        // pressure load
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropLoad" );
        aParameterLists( tPropIndex ).set( "function_parameters", tLoad );

        // Dirichlet structure
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropDirichletStruct" );
        aParameterLists( tPropIndex ).set( "function_parameters", "0.0;0.0" );

        // heat flux from outside
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropHeatFlux" );
        aParameterLists( tPropIndex ).set( "function_parameters", tHeatLoad );
        aParameterLists( tPropIndex ).set( "value_function", "Func_Const" );

        // time continuity weights
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropWeightCurrent" );
        aParameterLists( tPropIndex ).set( "function_parameters", "100.0" );

        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropWeightPrevious" );
        aParameterLists( tPropIndex ).set( "function_parameters", "100.0" );

        // Initial Temperature
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropInitialCondition" );
        aParameterLists( tPropIndex ).set( "value_function", "Func_Initial_Condition" );

        ////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////

        //------------------------------------------------------------------------------
        // DIFFUSION
        //------------------------------------------------------------------------------

        // create parameter list for constitutive model - skin
        aParameterLists( tCMIndex ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists( tCMIndex ).set( "constitutive_name", "CMDiffusion1" );
        aParameterLists( tCMIndex ).set( "phase_name", "Phase1" );
        aParameterLists( tCMIndex ).set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
        aParameterLists( tCMIndex ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists( tCMIndex ).set( "properties", std::string( "PropConductivity1 , Conductivity;" ) + std::string( "PropDensity1      , Density;" ) + std::string( "PropHeatCapacity1 , HeatCapacity" ) );

        // diffusion with phase change - PCM - 1
        aParameterLists( tCMIndex ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists( tCMIndex ).set( "constitutive_name", "CMDiffusion2" );
        aParameterLists( tCMIndex ).set( "phase_name", "Phase2" );
        aParameterLists( tCMIndex ).set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO_PC );
        aParameterLists( tCMIndex ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists( tCMIndex ).set( "properties", std::string( "PropConductivity2, Conductivity;" ) + std::string( "PropDensity2     , Density;" ) + std::string( "PropHeatCapacity2, HeatCapacity;" ) + std::string( "PropLatentHeat2  , LatentHeat;" ) + std::string( "PropPCTemp2      , PCTemp;" ) + std::string( "PropPhaseState2  , PhaseStateFunction;" ) + std::string( "PropPCconst2     , PhaseChangeConst" ) );

        //------------------------------------------------------------------------------
        // LINEAR ELASTICITY
        //------------------------------------------------------------------------------

        // linear elasticity - skin
        aParameterLists( tCMIndex ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists( tCMIndex ).set( "constitutive_name", "CMStrucLinIso1" );
        aParameterLists( tCMIndex ).set( "phase_name", "Phase1" );
        aParameterLists( tCMIndex ).set( "model_type", fem::Model_Type::PLANE_STRESS );
        aParameterLists( tCMIndex ).set( "constitutive_type", fem::Constitutive_Type::STRUC_LIN_ISO );
        aParameterLists( tCMIndex ).set( "dof_dependencies", std::pair< std::string, std::string >( sAllDofTypes, "Displacement,Temperature" ) );
        aParameterLists( tCMIndex ).set( "properties", std::string( "PropYoungsModulus1,    YoungsModulus;" ) + std::string( "PropPoissonRatio1,     PoissonRatio;" ) + std::string( "PropThermalExpansion1, CTE;" ) + std::string( "PropReferenceTemp,     ReferenceTemperature" ) );

        // linear elasticity - fins
        aParameterLists( tCMIndex ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists( tCMIndex ).set( "constitutive_name", "CMStrucLinIso2" );
        aParameterLists( tCMIndex ).set( "phase_name", "Phase2" );
        aParameterLists( tCMIndex ).set( "model_type", fem::Model_Type::PLANE_STRESS );
        aParameterLists( tCMIndex ).set( "constitutive_type", fem::Constitutive_Type::STRUC_LIN_ISO );
        aParameterLists( tCMIndex ).set( "dof_dependencies", std::pair< std::string, std::string >( sAllDofTypes, "Displacement,Temperature" ) );
        aParameterLists( tCMIndex ).set( "properties", std::string( "PropYoungsModulus2,    YoungsModulus;" ) + std::string( "PropPoissonRatio2,     PoissonRatio;" ) + std::string( "PropThermalExpansion2, CTE;" ) + std::string( "PropReferenceTemp,          ReferenceTemperature" ) );

        ////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////

        //------------------------------------------------------------------------------
        // GGLS
        //------------------------------------------------------------------------------

        // create parameter list for GGLS stabilization parameter for Skin
        aParameterLists( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( tSPIndex ).set( "stabilization_name", std::string( "SPGGLSDiffusion1" ) );
        aParameterLists( tSPIndex ).set( "leader_phase_name", "Phase1" );
        aParameterLists( tSPIndex ).set( "stabilization_type", fem::Stabilization_Type::GGLS_DIFFUSION );
        aParameterLists( tSPIndex ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists( tSPIndex ).set( "leader_properties", std::string( "PropConductivity1    , Conductivity;" ) + std::string( "PropDensity1         , Density;" ) + std::string( "PropHeatCapacity1    , HeatCapacity;" ) + std::string( "PropLatentHeat_Dummy , LatentHeat;" ) + std::string( "PropPCTemp_Dummy     , PCTemp;" ) + std::string( "PropPhaseState2      , PhaseStateFunction;" ) + std::string( "PropPCconst2         , PhaseChangeConst" ) );

        // create parameter list for GGLS stabilization parameter for PCM
        aParameterLists( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( tSPIndex ).set( "stabilization_name", std::string( "SPGGLSDiffusion2" ) );
        aParameterLists( tSPIndex ).set( "leader_phase_name", "Phase2" );
        aParameterLists( tSPIndex ).set( "stabilization_type", fem::Stabilization_Type::GGLS_DIFFUSION );
        aParameterLists( tSPIndex ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists( tSPIndex ).set( "leader_properties", std::string( "PropConductivity2 , Conductivity;" ) + std::string( "PropDensity2      , Density;" ) + std::string( "PropHeatCapacity2 , HeatCapacity;" ) + std::string( "PropLatentHeat2   , LatentHeat;" ) + std::string( "PropPCTemp2       , PCTemp;" ) + std::string( "PropPhaseState2   , PhaseStateFunction;" ) + std::string( "PropPCconst2      , PhaseChangeConst" ) );

        //------------------------------------------------------------------------------
        // NITSCHE DIRICHLET
        //------------------------------------------------------------------------------

        // Displacements - Shell - back wall
        aParameterLists( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( tSPIndex ).set( "stabilization_name", "SPNitscheStruc" );
        aParameterLists( tSPIndex ).set( "leader_phase_name", "Phase1" );
        aParameterLists( tSPIndex ).set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        aParameterLists( tSPIndex ).set( "function_parameters", "100.0" );
        aParameterLists( tSPIndex ).set( "leader_properties", "PropYoungsModulus1,Material" );

        //------------------------------------------------------------------------------
        // NITSCHE INTERFACE
        //------------------------------------------------------------------------------

        // Temperature - Skin - PCM
        aParameterLists( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( tSPIndex ).set( "stabilization_name", std::string( "SPInterfaceNitscheTemp" ) );
        aParameterLists( tSPIndex ).set( "leader_phase_name", "Phase1" );
        aParameterLists( tSPIndex ).set( "follower_phase_name", "Phase2" );
        aParameterLists( tSPIndex ).set( "stabilization_type", fem::Stabilization_Type::NITSCHE_INTERFACE );
        aParameterLists( tSPIndex ).set( "function_parameters", std::string( "100.0" ) );
        aParameterLists( tSPIndex ).set( "leader_properties", std::string( "PropConductivity1,Material" ) );
        aParameterLists( tSPIndex ).set( "follower_properties", std::string( "PropConductivity2,Material" ) );

        // Displacements - Skin - Fins
        aParameterLists( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( tSPIndex ).set( "stabilization_name", std::string( "SPInterfaceNitscheStruct" ) );
        aParameterLists( tSPIndex ).set( "leader_phase_name", "Phase1" );
        aParameterLists( tSPIndex ).set( "follower_phase_name", "Phase2" );
        aParameterLists( tSPIndex ).set( "stabilization_type", fem::Stabilization_Type::NITSCHE_INTERFACE );
        aParameterLists( tSPIndex ).set( "function_parameters", std::string( "100.0" ) );
        aParameterLists( tSPIndex ).set( "leader_properties", std::string( "PropYoungsModulus1,Material" ) );
        aParameterLists( tSPIndex ).set( "follower_properties", std::string( "PropYoungsModulus2,Material" ) );

        //------------------------------------------------------------------------------
        // GHOST
        //------------------------------------------------------------------------------

        // bulk Ghost - 1 - Temperature
        aParameterLists( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( tSPIndex ).set( "stabilization_name", std::string( "SPGPTemp1" ) );
        aParameterLists( tSPIndex ).set( "leader_phase_name", "Phase1" );
        aParameterLists( tSPIndex ).set( "follower_phase_name", "Phase1" );
        aParameterLists( tSPIndex ).set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        aParameterLists( tSPIndex ).set( "function_parameters", std::string( "0.01" ) );
        aParameterLists( tSPIndex ).set( "leader_properties", std::string( "PropConductivity1,Material" ) );
        aParameterLists( tSPIndex ).set( "follower_properties", std::string( "PropConductivity1,Material" ) );

        // bulk Ghost - 2 - Temperature
        aParameterLists( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( tSPIndex ).set( "stabilization_name", std::string( "SPGPTemp2" ) );
        aParameterLists( tSPIndex ).set( "leader_phase_name", "Phase2" );
        aParameterLists( tSPIndex ).set( "follower_phase_name", "Phase2" );
        aParameterLists( tSPIndex ).set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        aParameterLists( tSPIndex ).set( "function_parameters", std::string( "0.01" ) );
        aParameterLists( tSPIndex ).set( "leader_properties", std::string( "PropConductivity2,Material" ) );
        aParameterLists( tSPIndex ).set( "follower_properties", std::string( "PropConductivity2,Material" ) );

        // ======= //

        // bulk Ghost - 1 - Displacements
        aParameterLists( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( tSPIndex ).set( "stabilization_name", std::string( "SPGPStruct1" ) );
        aParameterLists( tSPIndex ).set( "leader_phase_name", "Phase1" );
        aParameterLists( tSPIndex ).set( "follower_phase_name", "Phase1" );
        aParameterLists( tSPIndex ).set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        aParameterLists( tSPIndex ).set( "function_parameters", std::string( "0.01" ) );
        aParameterLists( tSPIndex ).set( "leader_properties", std::string( "PropYoungsModulus1,Material" ) );
        aParameterLists( tSPIndex ).set( "follower_properties", std::string( "PropYoungsModulus1,Material" ) );

        // bulk Ghost - 2 - Displacements
        aParameterLists( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( tSPIndex ).set( "stabilization_name", std::string( "SPGPStruct2" ) );
        aParameterLists( tSPIndex ).set( "leader_phase_name", "Phase2" );
        aParameterLists( tSPIndex ).set( "follower_phase_name", "Phase2" );
        aParameterLists( tSPIndex ).set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        aParameterLists( tSPIndex ).set( "function_parameters", std::string( "0.01" ) );
        aParameterLists( tSPIndex ).set( "leader_properties", std::string( "PropYoungsModulus2,Material" ) );
        aParameterLists( tSPIndex ).set( "follower_properties", std::string( "PropYoungsModulus2,Material" ) );

        ////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////
        //------------------------------------------------------------------------------
        // BULK IWGs
        //------------------------------------------------------------------------------

        // diffusion - Skin
        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", std::string( "IWGDiffusionBulk1" ) );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        aParameterLists( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::BULK );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "Phase1" );
        aParameterLists( tIWGIndex ).set( "dof_residual", std::string( "TEMP" ) );
        aParameterLists( tIWGIndex ).set( "leader_dof_dependencies", std::string( "TEMP" ) );
        aParameterLists( tIWGIndex ).set( "leader_constitutive_models", std::string( "CMDiffusion1,Diffusion" ) );
        aParameterLists( tIWGIndex ).set( "stabilization_parameters", std::string( "SPGGLSDiffusion1,GGLSParam" ) );

        // diffusion - PCM
        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", std::string( "IWGDiffusionBulk2" ) );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        aParameterLists( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::BULK );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "Phase2" );
        aParameterLists( tIWGIndex ).set( "dof_residual", std::string( "TEMP" ) );
        aParameterLists( tIWGIndex ).set( "leader_dof_dependencies", std::string( "TEMP" ) );
        aParameterLists( tIWGIndex ).set( "leader_constitutive_models", std::string( "CMDiffusion2,Diffusion" ) );
        aParameterLists( tIWGIndex ).set( "stabilization_parameters", std::string( "SPGGLSDiffusion2,GGLSParam" ) );

        // ======= //

        // linear elasticity - Skin
        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGStructBulk1" );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_BULK );
        aParameterLists( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::BULK );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "Phase1" );
        aParameterLists( tIWGIndex ).set( "dof_residual", sStructDofTypes );
        aParameterLists( tIWGIndex ).set( "leader_dof_dependencies", sAllDofTypes );
        aParameterLists( tIWGIndex ).set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );

        // linear elasticity - PCM
        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGStructBulk2" );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_BULK );
        aParameterLists( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::BULK );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "Phase2" );
        aParameterLists( tIWGIndex ).set( "dof_residual", sStructDofTypes );
        aParameterLists( tIWGIndex ).set( "leader_dof_dependencies", sAllDofTypes );
        aParameterLists( tIWGIndex ).set( "leader_constitutive_models", "CMStrucLinIso2,ElastLinIso" );

        //------------------------------------------------------------------------------
        // NEUMANN BCs - IWGs
        //------------------------------------------------------------------------------

        // heat flux on outside of Shell
        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", std::string( "IWGHeatFlux" ) );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_NEUMANN );
        aParameterLists( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "Phase1" );
        aParameterLists( tIWGIndex ).set( "side_ordinals", "4" );
        aParameterLists( tIWGIndex ).set( "dof_residual", "TEMP" );
        aParameterLists( tIWGIndex ).set( "leader_dof_dependencies", sAllDofTypes );
        aParameterLists( tIWGIndex ).set( "leader_properties", "PropHeatFlux,Neumann" );

        // pressure pushing on outside of Shell
        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGNeumannLoad1" );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_NEUMANN );
        aParameterLists( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "Phase1" );
        aParameterLists( tIWGIndex ).set( "side_ordinals", "3" );
        aParameterLists( tIWGIndex ).set( "dof_residual", sStructDofTypes );
        aParameterLists( tIWGIndex ).set( "leader_dof_dependencies", sAllDofTypes );
        aParameterLists( tIWGIndex ).set( "leader_properties", "PropLoad,Pressure" );

        // pressure pushing on outside of Shell
        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGNeumannLoad2" );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_NEUMANN );
        aParameterLists( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "Phase2" );
        aParameterLists( tIWGIndex ).set( "side_ordinals", "3" );
        aParameterLists( tIWGIndex ).set( "dof_residual", sStructDofTypes );
        aParameterLists( tIWGIndex ).set( "leader_dof_dependencies", sAllDofTypes );
        aParameterLists( tIWGIndex ).set( "leader_properties", "PropLoad,Pressure" );

        //------------------------------------------------------------------------------
        // DIRICHLET BCS - IWGs
        //------------------------------------------------------------------------------

        // displacements - skin - back wall
        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGDirichletStruct" );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "Phase1" );
        aParameterLists( tIWGIndex ).set( "side_ordinals", "4" );
        aParameterLists( tIWGIndex ).set( "dof_residual", sStructDofTypes );
        aParameterLists( tIWGIndex ).set( "leader_dof_dependencies", sAllDofTypes );
        aParameterLists( tIWGIndex ).set( "leader_properties", "PropDirichletStruct,Dirichlet" );
        aParameterLists( tIWGIndex ).set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        aParameterLists( tIWGIndex ).set( "stabilization_parameters", "SPNitscheStruc,DirichletNitsche" );

        //------------------------------------------------------------------------------
        // INTERFACE BCS - IWGs
        //------------------------------------------------------------------------------

        // Temperature - 1 - 2
        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", std::string( "IWGInterfaceTEMP" ) );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE );
        aParameterLists( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "Phase1" );
        aParameterLists( tIWGIndex ).set( "follower_phase_name", "Phase2" );
        aParameterLists( tIWGIndex ).set( "dof_residual", std::string( "TEMP" ) );
        aParameterLists( tIWGIndex ).set( "leader_dof_dependencies", sAllDofTypes );
        aParameterLists( tIWGIndex ).set( "follower_dof_dependencies", sAllDofTypes );
        aParameterLists( tIWGIndex ).set( "leader_constitutive_models", std::string( "CMDiffusion1,Diffusion" ) );
        aParameterLists( tIWGIndex ).set( "follower_constitutive_models", std::string( "CMDiffusion2,Diffusion" ) );
        aParameterLists( tIWGIndex ).set( "stabilization_parameters", std::string( "SPInterfaceNitscheTemp,NitscheInterface" ) );

        // Displacements - 1 - 2
        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", std::string( "IWGInterfaceStruct" ) );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_INTERFACE_SYMMETRIC_NITSCHE );
        aParameterLists( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "Phase1" );
        aParameterLists( tIWGIndex ).set( "follower_phase_name", "Phase2" );
        aParameterLists( tIWGIndex ).set( "dof_residual", sStructDofTypes );
        aParameterLists( tIWGIndex ).set( "leader_dof_dependencies", sAllDofTypes );
        aParameterLists( tIWGIndex ).set( "follower_dof_dependencies", sAllDofTypes );
        aParameterLists( tIWGIndex ).set( "leader_constitutive_models", std::string( "CMStrucLinIso1,ElastLinIso" ) );
        aParameterLists( tIWGIndex ).set( "follower_constitutive_models", std::string( "CMStrucLinIso2,ElastLinIso" ) );
        aParameterLists( tIWGIndex ).set( "stabilization_parameters", std::string( "SPInterfaceNitscheStruct,NitscheInterface" ) );

        //------------------------------------------------------------------------------
        // IWGs - GHOST
        //------------------------------------------------------------------------------

        // temperature - 1
        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", std::string( "IWGGPTemp1" ) );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
        aParameterLists( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "Phase1" );
        aParameterLists( tIWGIndex ).set( "follower_phase_name", "Phase1" );
        aParameterLists( tIWGIndex ).set( "dof_residual", std::string( "TEMP" ) );
        aParameterLists( tIWGIndex ).set( "leader_dof_dependencies", std::string( "TEMP" ) );
        aParameterLists( tIWGIndex ).set( "follower_dof_dependencies", std::string( "TEMP" ) );
        aParameterLists( tIWGIndex ).set( "stabilization_parameters", std::string( "SPGPTemp1,GhostSP" ) );

        // temperature - 2
        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", std::string( "IWGGPTemp2" ) );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
        aParameterLists( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "Phase2" );
        aParameterLists( tIWGIndex ).set( "follower_phase_name", "Phase2" );
        aParameterLists( tIWGIndex ).set( "dof_residual", std::string( "TEMP" ) );
        aParameterLists( tIWGIndex ).set( "leader_dof_dependencies", std::string( "TEMP" ) );
        aParameterLists( tIWGIndex ).set( "follower_dof_dependencies", std::string( "TEMP" ) );
        aParameterLists( tIWGIndex ).set( "stabilization_parameters", std::string( "SPGPTemp2,GhostSP" ) );

        // displacements - 1
        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", std::string( "IWGGPStruct1" ) );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
        aParameterLists( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "Phase1" );
        aParameterLists( tIWGIndex ).set( "follower_phase_name", "Phase1" );
        aParameterLists( tIWGIndex ).set( "dof_residual", sStructDofTypes );
        aParameterLists( tIWGIndex ).set( "leader_dof_dependencies", sStructDofTypes );
        aParameterLists( tIWGIndex ).set( "follower_dof_dependencies", sStructDofTypes );
        aParameterLists( tIWGIndex ).set( "stabilization_parameters", std::string( "SPGPStruct1,GhostSP" ) );

        // displacements - 2
        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", std::string( "IWGGPStruct2" ) );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
        aParameterLists( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "Phase2" );
        aParameterLists( tIWGIndex ).set( "follower_phase_name", "Phase2" );
        aParameterLists( tIWGIndex ).set( "dof_residual", sStructDofTypes );
        aParameterLists( tIWGIndex ).set( "leader_dof_dependencies", sStructDofTypes );
        aParameterLists( tIWGIndex ).set( "follower_dof_dependencies", sStructDofTypes );
        aParameterLists( tIWGIndex ).set( "stabilization_parameters", std::string( "SPGPStruct2,GhostSP" ) );

        //------------------------------------------------------------------------------
        // IWGs - TIME CONTINUITY
        //------------------------------------------------------------------------------

        // Time continuity
        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", std::string( "IWGTimeContinuityTemp2" ) );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::TIME_CONTINUITY_DOF );
        aParameterLists( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::BULK );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "Phase2" );
        aParameterLists( tIWGIndex ).set( "dof_residual", std::string( "TEMP" ) );
        aParameterLists( tIWGIndex ).set( "leader_dof_dependencies", sAllDofTypes );
        aParameterLists( tIWGIndex ).set( "leader_properties", std::string( "PropWeightCurrent,WeightCurrent;" ) + std::string( "PropWeightPrevious,WeightPrevious;" ) + std::string( "PropInitialCondition,InitialCondition" ) );
        aParameterLists( tIWGIndex ).set( "time_continuity", true );

        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", std::string( "IWGTimeContinuityTemp3" ) );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::TIME_CONTINUITY_DOF );
        aParameterLists( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::BULK );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "Phase2" );
        aParameterLists( tIWGIndex ).set( "dof_residual", std::string( "TEMP" ) );
        aParameterLists( tIWGIndex ).set( "leader_dof_dependencies", sAllDofTypes );
        aParameterLists( tIWGIndex ).set( "leader_properties", std::string( "PropWeightCurrent,WeightCurrent;" ) + std::string( "PropWeightPrevious,WeightPrevious;" ) + std::string( "PropInitialCondition,InitialCondition" ) );
        aParameterLists( tIWGIndex ).set( "time_continuity", true );

        ////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////
        // Nodal Temperature IQI
        aParameterLists( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", std::string( "IQIBulkTEMP_1" ) );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists( tIQIIndex ).set( "dof_quantity", "TEMP" );
        aParameterLists( tIQIIndex ).set( "leader_dof_dependencies", std::string( "TEMP" ) );
        aParameterLists( tIQIIndex ).set( "vectorial_field_index", 0 );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "Phase1" );
        aParameterLists( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", std::string( "IQIBulkTEMP_2" ) );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists( tIQIIndex ).set( "dof_quantity", "TEMP" );
        aParameterLists( tIQIIndex ).set( "leader_dof_dependencies", std::string( "TEMP" ) );
        aParameterLists( tIQIIndex ).set( "vectorial_field_index", 0 );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "Phase2" );

        // X-displacement
        aParameterLists( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIBulkDISPX_1" );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists( tIQIIndex ).set( "dof_quantity", sStructDofTypes );
        aParameterLists( tIQIIndex ).set( "leader_dof_dependencies", sStructDofTypes );
        aParameterLists( tIQIIndex ).set( "vectorial_field_index", 0 );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "Phase1" );
        aParameterLists( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIBulkDISPX_2" );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists( tIQIIndex ).set( "dof_quantity", sStructDofTypes );
        aParameterLists( tIQIIndex ).set( "leader_dof_dependencies", sStructDofTypes );
        aParameterLists( tIQIIndex ).set( "vectorial_field_index", 0 );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "Phase2" );

        // Y-displacement
        aParameterLists( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIBulkDISPY_1" );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists( tIQIIndex ).set( "dof_quantity", sStructDofTypes );
        aParameterLists( tIQIIndex ).set( "leader_dof_dependencies", sStructDofTypes );
        aParameterLists( tIQIIndex ).set( "vectorial_field_index", 1 );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "Phase1" );
        aParameterLists( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIBulkDISPY_2" );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists( tIQIIndex ).set( "dof_quantity", sStructDofTypes );
        aParameterLists( tIQIIndex ).set( "leader_dof_dependencies", sStructDofTypes );
        aParameterLists( tIQIIndex ).set( "vectorial_field_index", 1 );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "Phase2" );

        // ==== //

        // Max Temperature IQI
        aParameterLists( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", std::string( "IQIMaxTemp_1" ) );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::MAX_DOF );
        aParameterLists( tIQIIndex ).set( "dof_quantity", "TEMP" );
        aParameterLists( tIQIIndex ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( tIQIIndex ).set( "function_parameters", "1.0/2.0" );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "Phase1" );
        aParameterLists( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", std::string( "IQIMaxTemp_2" ) );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::MAX_DOF );
        aParameterLists( tIQIIndex ).set( "dof_quantity", "TEMP" );
        aParameterLists( tIQIIndex ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( tIQIIndex ).set( "function_parameters", "1.0/2.0" );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "Phase2" );

        // Strain Energy of Structure
        aParameterLists( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIStrainEnergy_1" );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::STRAIN_ENERGY );
        aParameterLists( tIQIIndex ).set( "leader_dof_dependencies", sAllDofTypes );
        aParameterLists( tIQIIndex ).set( "leader_constitutive_models", "CMStrucLinIso1,Elast" );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "Phase1" );
        aParameterLists( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIStrainEnergy_2" );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::STRAIN_ENERGY );
        aParameterLists( tIQIIndex ).set( "leader_dof_dependencies", sAllDofTypes );
        aParameterLists( tIQIIndex ).set( "leader_constitutive_models", "CMStrucLinIso2,Elast" );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "Phase2" );

        // Volume IQIs
        aParameterLists( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIVolume_1" );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::VOLUME );
        aParameterLists( tIQIIndex ).set( "leader_dof_dependencies", sAllDofTypes );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "Phase1" );
        aParameterLists( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIVolume_2" );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::VOLUME );
        aParameterLists( tIQIIndex ).set( "leader_dof_dependencies", sAllDofTypes );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "Phase2" );

        // create computation parameter list
        aParameterLists( tFEMIndex ).add_parameter_list( prm::create_computation_parameter_list() );
        // aParameterLists( tFEMIndex ).set( "print_physics_model", true );
    }

    void
    SOLParameterList( Module_Parameter_Lists& aParameterLists )
    {
        // time levels for structure dofs
        if ( gUseMixedTimeElements )
            sTLSD = "1";

        // ----------------------------------------------------------
        // initialize solver parameter list


        // ----------------------------------------------------------
        // linear solver algorithm

        if ( gUseBelosWithILUT )
        {
            aParameterLists( 7 ).add_parameter_list( moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::IFPACK ) );
            aParameterLists( 7 ).set( "Convergence Tolerance", 1e-12 );
            aParameterLists( 7 ).set( "ifpack_prec_type", "ILUT" );
            aParameterLists( 7 ).set( "fact: drop tolerance", 1e-10 );
            aParameterLists( 7 ).set( "fact: ilut level-of-fill", 25.0 );

            aParameterLists( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::BELOS_IMPL ) );
            aParameterLists( 0 ).set( "preconditioners", "0" );
        }
        else
        {
            aParameterLists( 7 ).add_parameter_list( moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::NONE ) );
            aParameterLists( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL ) );
        }

        // ----------------------------------------------------------
        // linear solver

        aParameterLists( 1 ).add_parameter_list( moris::prm::create_linear_solver_parameter_list() );

        // ----------------------------------------------------------

        if ( gHaveStaggeredFA )
        {
            // ----------------------------------------------------------
            // non-linear solver algorithms

            // NEWTON solver algorithm
            aParameterLists( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
            aParameterLists( 2 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
            aParameterLists( 2 ).set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
            aParameterLists( 2 ).set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
            aParameterLists( 2 ).set( "NLA_max_iter", tNLA_max_iter );
            aParameterLists( 2 ).set( "NLA_combined_res_jac_assembly", true );

            // NLBGS solver algorithm
            // NOTE: maximum iterations is set to 1 since the second (structural problem) is linear and this saves time
            aParameterLists( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
            aParameterLists( 2 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NLBGS_SOLVER );
            aParameterLists( 2 ).set( "NLA_max_iter", 1 );

            // NEWTON solver algorithm for linear problems with only one iteration (i.e. in structural part and monolythic adjoint solve)
            aParameterLists( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
            aParameterLists( 2 ).set( "NLA_rel_res_norm_drop", 1.0e-7 );
            aParameterLists( 2 ).set( "NLA_relaxation_parameter", 1.0 );
            aParameterLists( 2 ).set( "NLA_max_iter", 1 );
            aParameterLists( 2 ).set( "NLA_combined_res_jac_assembly", true );

            // ----------------------------------------------------------
            // non-linear solvers

            // NEWTON solver for (linear) structural problem and adjoint
            aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
            aParameterLists( 3 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
            aParameterLists( 3 ).set( "NLA_DofTypes", sStructDofTypes );
            aParameterLists( 3 ).set( "NLA_Nonlinear_solver_algorithms", "2" );

            // NEWTON solver for non-linear thermal problem
            aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
            aParameterLists( 3 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
            aParameterLists( 3 ).set( "NLA_DofTypes", "TEMP" );

            // NLBGS solver
            aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
            aParameterLists( 3 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NLBGS_SOLVER );
            aParameterLists( 3 ).set( "NLA_Sub_Nonlinear_Solver", "1,0" );
            aParameterLists( 3 ).set( "NLA_DofTypes", sAllDofTypes );
            aParameterLists( 3 ).set( "NLA_Nonlinear_solver_algorithms", "1" );

            // NEWTON solver for separate monolythic adjoint solve
            aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
            aParameterLists( 3 ).set( "NLA_DofTypes", sAllDofTypes );
            aParameterLists( 3 ).set( "NLA_Nonlinear_solver_algorithms", "2" );

            // ----------------------------------------------------------

            // MONOLYTHIC time solver
            aParameterLists( 4 ).add_parameter_list( moris::prm::create_time_solver_algorithm_parameter_list() );
            aParameterLists( 4 ).set( "TSA_Num_Time_Steps", tTSA_Num_Time_Steps );
            aParameterLists( 4 ).set( "TSA_Time_Frame", tTSA_Time_Frame );
            aParameterLists( 4 ).set( "TSA_Nonlinear_Solver", 2 );

            // if SA is non-staggered, use separate monolythic non-linear solver for adjoint solve
            if ( !gHaveStaggeredSA )
            {
                aParameterLists( 4 ).set( "TSA_Nonlinear_Sensitivity_Solver", 3 );
            }

            // ----------------------------------------------------------

            aParameterLists( 5 ).add_parameter_list( moris::prm::create_time_solver_parameter_list() );
            aParameterLists( 5 ).set( "TSA_DofTypes", sAllDofTypes );
            aParameterLists( 5 ).set( "TSA_Output_Indices", "0,1" );
            aParameterLists( 5 ).set( "TSA_Output_Criteria", "Output_Criterion,Output_Criterion" );
            aParameterLists( 5 ).set( "TSA_Initialize_Sol_Vec", "UX,0.0;UY,0.0;TEMP,1.0" );
            aParameterLists( 5 ).set( "TSA_time_level_per_type", "UX," + sTLSD + ";UY," + sTLSD + ";TEMP,2" );
        }    // end: staggered solver case

        // ----------------------------------------------------------

        else    // monolythic solver
        {
            aParameterLists( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
            aParameterLists( 2 ).set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
            aParameterLists( 2 ).set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
            aParameterLists( 2 ).set( "NLA_max_iter", tNLA_max_iter );
            aParameterLists( 2 ).set( "NLA_combined_res_jac_assembly", true );

            // ----------------------------------------------------------

            aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
            aParameterLists( 3 ).set( "NLA_DofTypes", "UX,UY;TEMP" );

            aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
            aParameterLists( 3 ).set( "NLA_DofTypes", "UX,UY;TEMP" );
            aParameterLists( 3 ).set( "NLA_Nonlinear_solver_algorithms", "0" );

            // ----------------------------------------------------------

            aParameterLists( 4 ).add_parameter_list( moris::prm::create_time_solver_algorithm_parameter_list() );
            aParameterLists( 4 ).set( "TSA_Num_Time_Steps", tTSA_Num_Time_Steps );
            aParameterLists( 4 ).set( "TSA_Time_Frame", tTSA_Time_Frame );

            // ----------------------------------------------------------

            aParameterLists( 5 ).add_parameter_list( moris::prm::create_time_solver_parameter_list() );
            aParameterLists( 5 ).set( "TSA_DofTypes", "UX,UY;TEMP" );
            aParameterLists( 5 ).set( "TSA_Initialize_Sol_Vec", "UX,0.0;UY,0.0;TEMP,1.0" );
            aParameterLists( 5 ).set( "TSA_Output_Indices", "0" );
            aParameterLists( 5 ).set( "TSA_Output_Criteria", "Output_Criterion" );
            aParameterLists( 5 ).set( "TSA_time_level_per_type", "UX," + sTLSD + ";UY," + sTLSD + ";TEMP,2" );
        }

        aParameterLists( 6 ).add_parameter_list( moris::prm::create_solver_warehouse_parameterlist() );
        aParameterLists( 6 ).set( "SOL_save_operator_to_matlab", "Mat.dat" );
    }

    void
    MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_msi_parameter_list() );
    }

    void
    VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        int tVisCounter = 0;

        aParameterLists( 0 ).add_parameter_list( prm::create_vis_parameter_list() );
        aParameterLists( 0 ).set( "Temp_Name", std::pair< std::string, std::string >( "./", "temp_material_1.exo" ) );
        aParameterLists( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName + "_material_1" ) );
        aParameterLists( 0 ).set( "Time_Offset", 100.0 );
        aParameterLists( 0 ).set( "Mesh_Type", vis::VIS_Mesh_Type::STANDARD );    // STANDARD_WITH_OVERLAP
        aParameterLists( 0 ).set( "Set_Names", tMat1Set );
        aParameterLists( 0 ).set( "Field_Names",
                "TEMP,UX,UY,"
                "MAX_DOF,STRAIN_ENERGY,VOLUME" );
        aParameterLists( 0 ).set( "Field_Type",
                "NODAL,NODAL,NODAL,"
                "GLOBAL,GLOBAL,GLOBAL" );
        aParameterLists( 0 ).set( "IQI_Names",
                "IQIBulkTEMP_1,IQIBulkDISPX_1,IQIBulkDISPY_1,"
                "IQIMaxTemp_1,IQIStrainEnergy_1,IQIVolume_1" );
        aParameterLists( 0 ).set( "Save_Frequency", 1 );
        aParameterLists( 0 ).set( "Output_Index", tVisCounter );
        tVisCounter++;

        aParameterLists( 0 ).add_parameter_list( prm::create_vis_parameter_list() );
        aParameterLists( 0 ).set( "Temp_Name", std::pair< std::string, std::string >( "./", "temp_material_2.exo" ) );
        aParameterLists( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName + "_material_2" ) );
        aParameterLists( 0 ).set( "Time_Offset", 100.0 );
        aParameterLists( 0 ).set( "Mesh_Type", vis::VIS_Mesh_Type::STANDARD );    // STANDARD_WITH_OVERLAP
        aParameterLists( 0 ).set( "Set_Names", tMat2Set );
        aParameterLists( 0 ).set( "Field_Names",
                "TEMP,UX,UY,"
                "MAX_DOF,STRAIN_ENERGY,VOLUME" );
        aParameterLists( 0 ).set( "Field_Type",
                "NODAL,NODAL,NODAL,"
                "GLOBAL,GLOBAL,GLOBAL" );
        aParameterLists( 0 ).set( "IQI_Names",
                "IQIBulkTEMP_2,IQIBulkDISPX_2,IQIBulkDISPY_2,"
                "IQIMaxTemp_2,IQIStrainEnergy_2,IQIVolume_2" );
        aParameterLists( 0 ).set( "Save_Frequency", 1 );
        aParameterLists( 0 ).set( "Output_Index", tVisCounter );
        tVisCounter++;
    }

    void
    MORISGENERALParameterList( Module_Parameter_Lists& aParameterLists )
    {
    }

    /* ------------------------------------------------------------------------ */
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
