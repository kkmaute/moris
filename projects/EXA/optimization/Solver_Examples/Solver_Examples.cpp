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
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    // initial temperature
    void
    Func_Initial_Condition(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
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
            const Vector< real >&     aGeometryParameters )
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
    OPTParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
    {
        tParameterlist.resize( 3 );
        tParameterlist( 0 ).resize( 1 );
        tParameterlist( 1 ).resize( 0 );
        tParameterlist( 2 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = moris::prm::create_opt_problem_parameter_list();
        tParameterlist( 0 )( 0 ).set( "is_optimization_problem", true );
        tParameterlist( 0 )( 0 ).set( "problem", "user_defined" );
        tParameterlist( 0 )( 0 ).set( "library", tLibraryName );

        tParameterlist( 2 )( 0 ) = moris::prm::create_gcmma_parameter_list();
        tParameterlist( 2 )( 0 ).set( "max_its", 1 );
        tParameterlist( 2 )( 0 ).set( "step_size", 0.2 );
    }

    /* ------------------------------------------------------------------------ */

    void
    HMRParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_hmr_parameter_list();

        tParameterlist( 0 )( 0 ).set( "number_of_elements_per_dimension", "32, 8" );
        tParameterlist( 0 )( 0 ).set( "domain_dimensions", "0.8,0.2" );
        tParameterlist( 0 )( 0 ).set( "domain_offset", "0.0,0.0" );
        tParameterlist( 0 )( 0 ).set( "domain_sidesets", "1,2,3,4" );
        tParameterlist( 0 )( 0 ).set( "lagrange_output_meshes", "0" );

        tParameterlist( 0 )( 0 ).set( "lagrange_orders", "1" );
        tParameterlist( 0 )( 0 ).set( "lagrange_pattern", "0" );
        tParameterlist( 0 )( 0 ).set( "bspline_orders", "1" );
        tParameterlist( 0 )( 0 ).set( "bspline_pattern", "0" );

        tParameterlist( 0 )( 0 ).set( "truncate_bsplines", 1 );
        tParameterlist( 0 )( 0 ).set( "use_number_aura", 1 );
        tParameterlist( 0 )( 0 ).set( "use_multigrid", 0 );
        tParameterlist( 0 )( 0 ).set( "severity_level", 0 );
    }

    /* ------------------------------------------------------------------------ */

    void
    XTKParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_xtk_parameter_list();
        tParameterlist( 0 )( 0 ).set( "decompose", true );
        tParameterlist( 0 )( 0 ).set( "decomposition_type", std::string( "conformal" ) );
        tParameterlist( 0 )( 0 ).set( "enrich", true );
        tParameterlist( 0 )( 0 ).set( "basis_rank", std::string( "bspline" ) );
        tParameterlist( 0 )( 0 ).set( "enrich_mesh_indices", std::string( "0" ) );
        tParameterlist( 0 )( 0 ).set( "ghost_stab", true );
        tParameterlist( 0 )( 0 ).set( "multigrid", false );
        tParameterlist( 0 )( 0 ).set( "verbose", true );
        tParameterlist( 0 )( 0 ).set( "print_enriched_ig_mesh", true );
        tParameterlist( 0 )( 0 ).set( "exodus_output_XTK_ig_mesh", true );
        tParameterlist( 0 )( 0 ).set( "high_to_low_dbl_side_sets", true );
        tParameterlist( 0 )( 0 ).set( "output_path", "./" );
        tParameterlist( 0 )( 0 ).set( "keep_all_opt_iters", true );
    }

    /* ------------------------------------------------------------------------ */

    void
    GENParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
    {
        tParameterlist.resize( 3 );
        tParameterlist( 0 ).resize( 1 );

        // Main GEN parameter list
        tParameterlist( 0 )( 0 ) = prm::create_gen_parameter_list();
        tParameterlist( 0 )( 0 ).set( "IQI_types", "IQIMaxTemp_1,IQIMaxTemp_2,IQIStrainEnergy_1,IQIStrainEnergy_2,IQIVolume_1" );
        tParameterlist( 0 )( 0 ).set( "output_mesh_file", tGENOutputFile );
        tParameterlist( 0 )( 0 ).set( "initial_advs", 0.4 );
        tParameterlist( 0 )( 0 ).set( "lower_bounds", 0.1 );
        tParameterlist( 0 )( 0 ).set( "upper_bounds", 0.7 );

        // init geometry counter
        uint tGeoCounter = 0;

        // Interface
        //        tParameterlist( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        //        tParameterlist( 1 )( tGeoCounter ).set( "field_function_name", "Interface" );
        //        tGeoCounter++;

        tParameterlist( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
        tParameterlist( 1 )( tGeoCounter ).set( "field_variable_indices", 0u );
        tParameterlist( 1 )( tGeoCounter ).set( "adv_indices", 0u );
        tParameterlist( 1 )( tGeoCounter ).set( "constant_parameters", 0.0, 1.0, 0.0 );
        tParameterlist( 1 )( tGeoCounter ).set( "isocontour_tolerance", 1e-12 );    // Interface tolerance based on geometry value
        tGeoCounter++;
    }

    /* ------------------------------------------------------------------------ */

    void
    FEMParameterList( Vector< Vector< Parameter_List > >& tParameterList )
    {
        std::string sStructDofTypes = "UX,UY";
        std::string sAllDofTypes    = "UX,UY;TEMP";

        // create a cell of cell of parameter list for fem
        tParameterList.resize( 9 );
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

        // phase info
        uint tPhaseCounter = 0;

        tParameterList( tPhaseIndex ).push_back( prm::create_phase_parameter_list() );
        tParameterList( tPhaseIndex )( tPhaseCounter ).set( "phase_name", "Phase1" );
        tParameterList( tPhaseIndex )( tPhaseCounter ).set( "phase_indices", tPhase1 );
        tPhaseCounter++;

        tParameterList( tPhaseIndex ).push_back( prm::create_phase_parameter_list() );
        tParameterList( tPhaseIndex )( tPhaseCounter ).set( "phase_name", "Phase2" );
        tParameterList( tPhaseIndex )( tPhaseCounter ).set( "phase_indices", tPhase2 );
        tPhaseCounter++;

        ////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////
        // init property counter
        uint tPropCounter = 0;

        //------------------------------------------------------------------------------
        // MATERIAL PARAMETERS - STRUCTURE (ni-w-alloy?)
        //------------------------------------------------------------------------------

        // Density Shell
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropDensity1" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", tDensity1 );
        tPropCounter++;

        // Heat Capacity Shell
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropHeatCapacity1" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", tHeatCapacity1 );
        tPropCounter++;

        // Conductivity Shell
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropConductivity1" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", tConductivity1 );
        tPropCounter++;

        // Youngs Modulus Shell
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropYoungsModulus1" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", tYoungsModulus1 );
        tPropCounter++;

        // Poisson Ratio Shell
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropPoissonRatio1" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", tPoissonRatio1 );
        tPropCounter++;

        // CTE for Shell
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropThermalExpansion1" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", tThermalExpansion1 );
        tPropCounter++;

        //------------------------------------------------------------------------------
        // MATERIAL PARAMETERS - PCM (Al-Cu-Si-alloy?)
        //------------------------------------------------------------------------------

        // Density of PCM
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropDensity2" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", tDensity2 );
        tPropCounter++;

        // Heat Capacity of PCM
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropHeatCapacity2" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", tHeatCapacity2 );
        tPropCounter++;

        // Conductivity of PCM
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropConductivity2" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", tConductivity2 );
        tPropCounter++;

        // Latent Heat of PCM
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropLatentHeat2" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", tLatentHeat2 );
        tPropCounter++;

        // Phase Change Temperature of PCM
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropPCTemp2" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", tPCTemp2 );
        tPropCounter++;

        // Phase Change Temperature Range of PCM
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropPCconst2" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", tPCConst2 );
        tPropCounter++;

        // Cubic Phase State Function for phase change model
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropPhaseState2" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", "2.0" );
        tPropCounter++;

        // Youngs Modulus for PCM
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropYoungsModulus2" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", tYoungsModulus2 );
        tPropCounter++;

        // Poisson Ratio for PCM
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropPoissonRatio2" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", tPoissonRatio2 );
        tPropCounter++;

        // CTE for PCM
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropThermalExpansion2" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", tThermalExpansion2 );
        tPropCounter++;

        //------------------------------------------------------------------------------
        // OTHER MATERIAL PARAMETERS
        //------------------------------------------------------------------------------

        // Dummy latent heat for non-pc material
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropLatentHeat_Dummy" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", "0.0" );
        tPropCounter++;

        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropPCTemp_Dummy" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", "10000.0" );
        tPropCounter++;

        //------------------------------------------------------------------------------
        // BOUNDARY CONDITIONS
        //------------------------------------------------------------------------------

        // reference temperature for thermal expansion
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropReferenceTemp" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", tReferenceTemp );
        tPropCounter++;

        // pressure load
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropLoad" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", tLoad );
        tPropCounter++;

        // Dirichlet structure
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropDirichletStruct" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", "0.0;0.0" );
        tPropCounter++;

        // heat flux from outside
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropHeatFlux" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", tHeatLoad );
        tParameterList( tPropIndex )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // time continuity weights
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropWeightCurrent" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", "100.0" );
        tPropCounter++;

        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropWeightPrevious" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", "100.0" );
        tPropCounter++;

        // Initial Temperature
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropInitialCondition" );
        tParameterList( tPropIndex )( tPropCounter ).set( "value_function", "Func_Initial_Condition" );
        tPropCounter++;

        ////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////
        // init CM counter
        uint tCMCounter = 0;

        //------------------------------------------------------------------------------
        // DIFFUSION
        //------------------------------------------------------------------------------

        // create parameter list for constitutive model - skin
        tParameterList( tCMIndex ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( tCMIndex )( tCMCounter ).set( "constitutive_name", "CMDiffusion1" );
        tParameterList( tCMIndex )( tCMCounter ).set( "phase_name", "Phase1" );
        tParameterList( tCMIndex )( tCMCounter ).set( "constitutive_type",  fem::Constitutive_Type::DIFF_LIN_ISO ) ;
        tParameterList( tCMIndex )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( tCMIndex )( tCMCounter ).set( "properties", std::string( "PropConductivity1 , Conductivity;" ) + std::string( "PropDensity1      , Density;" ) + std::string( "PropHeatCapacity1 , HeatCapacity" ) );
        tCMCounter++;

        // diffusion with phase change - PCM - 1
        tParameterList( tCMIndex ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( tCMIndex )( tCMCounter ).set( "constitutive_name", "CMDiffusion2" );
        tParameterList( tCMIndex )( tCMCounter ).set( "phase_name", "Phase2" );
        tParameterList( tCMIndex )( tCMCounter ).set( "constitutive_type",  fem::Constitutive_Type::DIFF_LIN_ISO_PC ) ;
        tParameterList( tCMIndex )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( tCMIndex )( tCMCounter ).set( "properties", std::string( "PropConductivity2, Conductivity;" ) + std::string( "PropDensity2     , Density;" ) + std::string( "PropHeatCapacity2, HeatCapacity;" ) + std::string( "PropLatentHeat2  , LatentHeat;" ) + std::string( "PropPCTemp2      , PCTemp;" ) + std::string( "PropPhaseState2  , PhaseStateFunction;" ) + std::string( "PropPCconst2     , PhaseChangeConst" ) );
        tCMCounter++;

        //------------------------------------------------------------------------------
        // LINEAR ELASTICITY
        //------------------------------------------------------------------------------

        // linear elasticity - skin
        tParameterList( tCMIndex ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( tCMIndex )( tCMCounter ).set( "constitutive_name", "CMStrucLinIso1" );
        tParameterList( tCMIndex )( tCMCounter ).set( "phase_name", "Phase1" );
        tParameterList( tCMIndex )( tCMCounter ).set( "model_type",  fem::Model_Type::PLANE_STRESS ) ;
        tParameterList( tCMIndex )( tCMCounter ).set( "constitutive_type",  fem::Constitutive_Type::STRUC_LIN_ISO ) ;
        tParameterList( tCMIndex )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( sAllDofTypes, "Displacement,Temperature" ) );
        tParameterList( tCMIndex )( tCMCounter ).set( "properties", std::string( "PropYoungsModulus1,    YoungsModulus;" ) + std::string( "PropPoissonRatio1,     PoissonRatio;" ) + std::string( "PropThermalExpansion1, CTE;" ) + std::string( "PropReferenceTemp,     ReferenceTemperature" ) );

        tCMCounter++;

        // linear elasticity - fins
        tParameterList( tCMIndex ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( tCMIndex )( tCMCounter ).set( "constitutive_name", "CMStrucLinIso2" );
        tParameterList( tCMIndex )( tCMCounter ).set( "phase_name", "Phase2" );
        tParameterList( tCMIndex )( tCMCounter ).set( "model_type",  fem::Model_Type::PLANE_STRESS ) ;
        tParameterList( tCMIndex )( tCMCounter ).set( "constitutive_type",  fem::Constitutive_Type::STRUC_LIN_ISO ) ;
        tParameterList( tCMIndex )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( sAllDofTypes, "Displacement,Temperature" ) );
        tParameterList( tCMIndex )( tCMCounter ).set( "properties", std::string( "PropYoungsModulus2,    YoungsModulus;" ) + std::string( "PropPoissonRatio2,     PoissonRatio;" ) + std::string( "PropThermalExpansion2, CTE;" ) + std::string( "PropReferenceTemp,          ReferenceTemperature" ) );
        tCMCounter++;

        ////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////
        // init SP counter
        uint tSPCounter = 0;

        //------------------------------------------------------------------------------
        // GGLS
        //------------------------------------------------------------------------------

        // create parameter list for GGLS stabilization parameter for Skin
        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name", std::string( "SPGGLSDiffusion1" ) );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_phase_name", "Phase1" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type",  fem::Stabilization_Type::GGLS_DIFFUSION ) ;
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_properties", std::string( "PropConductivity1    , Conductivity;" ) + std::string( "PropDensity1         , Density;" ) + std::string( "PropHeatCapacity1    , HeatCapacity;" ) + std::string( "PropLatentHeat_Dummy , LatentHeat;" ) + std::string( "PropPCTemp_Dummy     , PCTemp;" ) + std::string( "PropPhaseState2      , PhaseStateFunction;" ) + std::string( "PropPCconst2         , PhaseChangeConst" ) );
        tSPCounter++;

        // create parameter list for GGLS stabilization parameter for PCM
        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name", std::string( "SPGGLSDiffusion2" ) );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_phase_name", "Phase2" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type",  fem::Stabilization_Type::GGLS_DIFFUSION ) ;
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_properties", std::string( "PropConductivity2 , Conductivity;" ) + std::string( "PropDensity2      , Density;" ) + std::string( "PropHeatCapacity2 , HeatCapacity;" ) + std::string( "PropLatentHeat2   , LatentHeat;" ) + std::string( "PropPCTemp2       , PCTemp;" ) + std::string( "PropPhaseState2   , PhaseStateFunction;" ) + std::string( "PropPCconst2      , PhaseChangeConst" ) );
        tSPCounter++;

        //------------------------------------------------------------------------------
        // NITSCHE DIRICHLET
        //------------------------------------------------------------------------------

        // Displacements - Shell - back wall
        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name", "SPNitscheStruc" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_phase_name", "Phase1" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type",  fem::Stabilization_Type::DIRICHLET_NITSCHE ) ;
        tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters", "100.0" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_properties", "PropYoungsModulus1,Material" );
        tSPCounter++;

        //------------------------------------------------------------------------------
        // NITSCHE INTERFACE
        //------------------------------------------------------------------------------

        // Temperature - Skin - PCM
        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name", std::string( "SPInterfaceNitscheTemp" ) );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_phase_name", "Phase1" );
        tParameterList( tSPIndex )( tSPCounter ).set( "follower_phase_name", "Phase2" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type",  fem::Stabilization_Type::NITSCHE_INTERFACE ) ;
        tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters", std::string( "100.0" ) );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_properties", std::string( "PropConductivity1,Material" ) );
        tParameterList( tSPIndex )( tSPCounter ).set( "follower_properties", std::string( "PropConductivity2,Material" ) );
        tSPCounter++;

        // Displacements - Skin - Fins
        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name", std::string( "SPInterfaceNitscheStruct" ) );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_phase_name", "Phase1" );
        tParameterList( tSPIndex )( tSPCounter ).set( "follower_phase_name", "Phase2" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type",  fem::Stabilization_Type::NITSCHE_INTERFACE ) ;
        tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters", std::string( "100.0" ) );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_properties", std::string( "PropYoungsModulus1,Material" ) );
        tParameterList( tSPIndex )( tSPCounter ).set( "follower_properties", std::string( "PropYoungsModulus2,Material" ) );
        tSPCounter++;

        //------------------------------------------------------------------------------
        // GHOST
        //------------------------------------------------------------------------------

        // bulk Ghost - 1 - Temperature
        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name", std::string( "SPGPTemp1" ) );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_phase_name", "Phase1" );
        tParameterList( tSPIndex )( tSPCounter ).set( "follower_phase_name", "Phase1" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type",  fem::Stabilization_Type::GHOST_DISPL ) ;
        tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters", std::string( "0.01" ) );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_properties", std::string( "PropConductivity1,Material" ) );
        tParameterList( tSPIndex )( tSPCounter ).set( "follower_properties", std::string( "PropConductivity1,Material" ) );
        tSPCounter++;

        // bulk Ghost - 2 - Temperature
        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name", std::string( "SPGPTemp2" ) );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_phase_name", "Phase2" );
        tParameterList( tSPIndex )( tSPCounter ).set( "follower_phase_name", "Phase2" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type",  fem::Stabilization_Type::GHOST_DISPL ) ;
        tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters", std::string( "0.01" ) );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_properties", std::string( "PropConductivity2,Material" ) );
        tParameterList( tSPIndex )( tSPCounter ).set( "follower_properties", std::string( "PropConductivity2,Material" ) );
        tSPCounter++;

        // ======= //

        // bulk Ghost - 1 - Displacements
        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name", std::string( "SPGPStruct1" ) );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_phase_name", "Phase1" );
        tParameterList( tSPIndex )( tSPCounter ).set( "follower_phase_name", "Phase1" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type",  fem::Stabilization_Type::GHOST_DISPL ) ;
        tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters", std::string( "0.01" ) );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_properties", std::string( "PropYoungsModulus1,Material" ) );
        tParameterList( tSPIndex )( tSPCounter ).set( "follower_properties", std::string( "PropYoungsModulus1,Material" ) );
        tSPCounter++;

        // bulk Ghost - 2 - Displacements
        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name", std::string( "SPGPStruct2" ) );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_phase_name", "Phase2" );
        tParameterList( tSPIndex )( tSPCounter ).set( "follower_phase_name", "Phase2" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type",  fem::Stabilization_Type::GHOST_DISPL ) ;
        tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters", std::string( "0.01" ) );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_properties", std::string( "PropYoungsModulus2,Material" ) );
        tParameterList( tSPIndex )( tSPCounter ).set( "follower_properties", std::string( "PropYoungsModulus2,Material" ) );
        tSPCounter++;

        ////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////
        // init IWG counter
        uint tIWGCounter = 0;

        //------------------------------------------------------------------------------
        // BULK IWGs
        //------------------------------------------------------------------------------

        // diffusion - Skin
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", std::string( "IWGDiffusionBulk1" ) );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_BULK ) ;
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type",  fem::Element_Type::BULK ) ;
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "Phase1" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", std::string( "TEMP" ) );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_dof_dependencies", std::string( "TEMP" ) );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", std::string( "CMDiffusion1,Diffusion" ) );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", std::string( "SPGGLSDiffusion1,GGLSParam" ) );
        tIWGCounter++;

        // diffusion - PCM
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", std::string( "IWGDiffusionBulk2" ) );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_BULK ) ;
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type",  fem::Element_Type::BULK ) ;
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "Phase2" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", std::string( "TEMP" ) );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_dof_dependencies", std::string( "TEMP" ) );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", std::string( "CMDiffusion2,Diffusion" ) );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", std::string( "SPGGLSDiffusion2,GGLSParam" ) );
        tIWGCounter++;

        // ======= //

        // linear elasticity - Skin
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGStructBulk1" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_BULK ) ;
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type",  fem::Element_Type::BULK ) ;
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "Phase1" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", sStructDofTypes );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_dof_dependencies", sAllDofTypes );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        tIWGCounter++;

        // linear elasticity - PCM
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGStructBulk2" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_BULK ) ;
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type",  fem::Element_Type::BULK ) ;
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "Phase2" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", sStructDofTypes );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_dof_dependencies", sAllDofTypes );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", "CMStrucLinIso2,ElastLinIso" );
        tIWGCounter++;

        //------------------------------------------------------------------------------
        // NEUMANN BCs - IWGs
        //------------------------------------------------------------------------------

        // heat flux on outside of Shell
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", std::string( "IWGHeatFlux" ) );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_NEUMANN ) ;
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type",  fem::Element_Type::SIDESET ) ;
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "Phase1" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "side_ordinals", "4" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "TEMP" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_dof_dependencies", sAllDofTypes );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_properties", "PropHeatFlux,Neumann" );
        tIWGCounter++;

        // pressure pushing on outside of Shell
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGNeumannLoad1" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_NEUMANN ) ;
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type",  fem::Element_Type::SIDESET ) ;
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "Phase1" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "side_ordinals", "3" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", sStructDofTypes );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_dof_dependencies", sAllDofTypes );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_properties", "PropLoad,Pressure" );
        tIWGCounter++;

        // pressure pushing on outside of Shell
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGNeumannLoad2" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_NEUMANN ) ;
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type",  fem::Element_Type::SIDESET ) ;
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "Phase2" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "side_ordinals", "3" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", sStructDofTypes );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_dof_dependencies", sAllDofTypes );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_properties", "PropLoad,Pressure" );
        tIWGCounter++;

        //------------------------------------------------------------------------------
        // DIRICHLET BCS - IWGs
        //------------------------------------------------------------------------------

        // displacements - skin - back wall
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGDirichletStruct" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE ) ;
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type",  fem::Element_Type::SIDESET ) ;
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "Phase1" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "side_ordinals", "4" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", sStructDofTypes );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_dof_dependencies", sAllDofTypes );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_properties", "PropDirichletStruct,Dirichlet" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SPNitscheStruc,DirichletNitsche" );
        tIWGCounter++;

        //------------------------------------------------------------------------------
        // INTERFACE BCS - IWGs
        //------------------------------------------------------------------------------

        // Temperature - 1 - 2
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", std::string( "IWGInterfaceTEMP" ) );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE ) ;
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type",  fem::Element_Type::DOUBLE_SIDESET ) ;
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "Phase1" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "follower_phase_name", "Phase2" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", std::string( "TEMP" ) );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_dof_dependencies", sAllDofTypes );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "follower_dof_dependencies", sAllDofTypes );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", std::string( "CMDiffusion1,Diffusion" ) );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "follower_constitutive_models", std::string( "CMDiffusion2,Diffusion" ) );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", std::string( "SPInterfaceNitscheTemp,NitscheInterface" ) );
        tIWGCounter++;

        // Displacements - 1 - 2
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", std::string( "IWGInterfaceStruct" ) );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_INTERFACE_SYMMETRIC_NITSCHE ) ;
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type",  fem::Element_Type::DOUBLE_SIDESET ) ;
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "Phase1" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "follower_phase_name", "Phase2" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", sStructDofTypes );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_dof_dependencies", sAllDofTypes );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "follower_dof_dependencies", sAllDofTypes );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", std::string( "CMStrucLinIso1,ElastLinIso" ) );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "follower_constitutive_models", std::string( "CMStrucLinIso2,ElastLinIso" ) );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", std::string( "SPInterfaceNitscheStruct,NitscheInterface" ) );
        tIWGCounter++;

        //------------------------------------------------------------------------------
        // IWGs - GHOST
        //------------------------------------------------------------------------------

        // temperature - 1
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", std::string( "IWGGPTemp1" ) );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type",  fem::Element_Type::DOUBLE_SIDESET ) ;
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "Phase1" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "follower_phase_name", "Phase1" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", std::string( "TEMP" ) );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_dof_dependencies", std::string( "TEMP" ) );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "follower_dof_dependencies", std::string( "TEMP" ) );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", std::string( "SPGPTemp1,GhostSP" ) );
        tIWGCounter++;

        // temperature - 2
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", std::string( "IWGGPTemp2" ) );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type",  fem::Element_Type::DOUBLE_SIDESET ) ;
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "Phase2" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "follower_phase_name", "Phase2" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", std::string( "TEMP" ) );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_dof_dependencies", std::string( "TEMP" ) );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "follower_dof_dependencies", std::string( "TEMP" ) );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", std::string( "SPGPTemp2,GhostSP" ) );
        tIWGCounter++;

        // displacements - 1
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", std::string( "IWGGPStruct1" ) );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type",  fem::Element_Type::DOUBLE_SIDESET ) ;
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "Phase1" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "follower_phase_name", "Phase1" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", sStructDofTypes );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_dof_dependencies", sStructDofTypes );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "follower_dof_dependencies", sStructDofTypes );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", std::string( "SPGPStruct1,GhostSP" ) );
        tIWGCounter++;

        // displacements - 2
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", std::string( "IWGGPStruct2" ) );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type",  fem::Element_Type::DOUBLE_SIDESET ) ;
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "Phase2" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "follower_phase_name", "Phase2" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", sStructDofTypes );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_dof_dependencies", sStructDofTypes );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "follower_dof_dependencies", sStructDofTypes );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", std::string( "SPGPStruct2,GhostSP" ) );
        tIWGCounter++;

        //------------------------------------------------------------------------------
        // IWGs - TIME CONTINUITY
        //------------------------------------------------------------------------------

        // Time continuity
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", std::string( "IWGTimeContinuityTemp2" ) );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type",  fem::IWG_Type::TIME_CONTINUITY_DOF ) ;
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type",  fem::Element_Type::BULK ) ;
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "Phase2" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", std::string( "TEMP" ) );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_dof_dependencies", sAllDofTypes );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_properties", std::string( "PropWeightCurrent,WeightCurrent;" ) + std::string( "PropWeightPrevious,WeightPrevious;" ) + std::string( "PropInitialCondition,InitialCondition" ) );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "time_continuity", true );
        tIWGCounter++;

        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", std::string( "IWGTimeContinuityTemp3" ) );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type",  fem::IWG_Type::TIME_CONTINUITY_DOF ) ;
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type",  fem::Element_Type::BULK ) ;
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "Phase2" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", std::string( "TEMP" ) );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_dof_dependencies", sAllDofTypes );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_properties", std::string( "PropWeightCurrent,WeightCurrent;" ) + std::string( "PropWeightPrevious,WeightPrevious;" ) + std::string( "PropInitialCondition,InitialCondition" ) );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "time_continuity", true );
        tIWGCounter++;

        ////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////
        // init IQI counter
        uint tIQICounter = 0;

        // Nodal Temperature IQI
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", std::string( "IQIBulkTEMP_1" ) );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type",  fem::IQI_Type::DOF ) ;
        tParameterList( tIQIIndex )( tIQICounter ).set( "dof_quantity", "TEMP" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_dof_dependencies", std::string( "TEMP" ) );
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index", 0 );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "Phase1" );
        tIQICounter++;
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", std::string( "IQIBulkTEMP_2" ) );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type",  fem::IQI_Type::DOF ) ;
        tParameterList( tIQIIndex )( tIQICounter ).set( "dof_quantity", "TEMP" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_dof_dependencies", std::string( "TEMP" ) );
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index", 0 );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "Phase2" );
        tIQICounter++;

        // X-displacement
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIBulkDISPX_1" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type",  fem::IQI_Type::DOF ) ;
        tParameterList( tIQIIndex )( tIQICounter ).set( "dof_quantity", sStructDofTypes );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_dof_dependencies", sStructDofTypes );
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index", 0 );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "Phase1" );
        tIQICounter++;
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIBulkDISPX_2" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type",  fem::IQI_Type::DOF ) ;
        tParameterList( tIQIIndex )( tIQICounter ).set( "dof_quantity", sStructDofTypes );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_dof_dependencies", sStructDofTypes );
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index", 0 );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "Phase2" );
        tIQICounter++;

        // Y-displacement
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIBulkDISPY_1" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type",  fem::IQI_Type::DOF ) ;
        tParameterList( tIQIIndex )( tIQICounter ).set( "dof_quantity", sStructDofTypes );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_dof_dependencies", sStructDofTypes );
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index", 1 );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "Phase1" );
        tIQICounter++;
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIBulkDISPY_2" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type",  fem::IQI_Type::DOF ) ;
        tParameterList( tIQIIndex )( tIQICounter ).set( "dof_quantity", sStructDofTypes );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_dof_dependencies", sStructDofTypes );
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index", 1 );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "Phase2" );
        tIQICounter++;

        // ==== //

        // Max Temperature IQI
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", std::string( "IQIMaxTemp_1" ) );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type",  fem::IQI_Type::MAX_DOF ) ;
        tParameterList( tIQIIndex )( tIQICounter ).set( "dof_quantity", "TEMP" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "function_parameters", "1.0/2.0" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "Phase1" );
        tIQICounter++;
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", std::string( "IQIMaxTemp_2" ) );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type",  fem::IQI_Type::MAX_DOF ) ;
        tParameterList( tIQIIndex )( tIQICounter ).set( "dof_quantity", "TEMP" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "function_parameters", "1.0/2.0" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "Phase2" );
        tIQICounter++;

        // Strain Energy of Structure
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIStrainEnergy_1" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type",  fem::IQI_Type::STRAIN_ENERGY ) ;
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_dof_dependencies", sAllDofTypes );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_constitutive_models", "CMStrucLinIso1,Elast" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "Phase1" );
        tIQICounter++;
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIStrainEnergy_2" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type",  fem::IQI_Type::STRAIN_ENERGY ) ;
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_dof_dependencies", sAllDofTypes );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_constitutive_models", "CMStrucLinIso2,Elast" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "Phase2" );
        tIQICounter++;

        // Volume IQIs
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIVolume_1" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type",  fem::IQI_Type::VOLUME ) ;
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_dof_dependencies", sAllDofTypes );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "Phase1" );
        tIQICounter++;
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIVolume_2" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type",  fem::IQI_Type::VOLUME ) ;
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_dof_dependencies", sAllDofTypes );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "Phase2" );
        tIQICounter++;

        // create computation parameter list
        tParameterList( tFEMIndex ).resize( 1 );
        tParameterList( tFEMIndex )( 0 ) = prm::create_computation_parameter_list();
        // tParameterList( tFEMIndex )( 0 ).set( "print_physics_model", true );
    }

    void
    SOLParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
    {
        // time levels for structure dofs
        if ( gUseMixedTimeElements )
            sTLSD = "1";

        // ----------------------------------------------------------
        // initialize solver parameter list

        tParameterlist.resize( 8 );
        for ( uint Ik = 0; Ik < 8; Ik++ )
        {
            tParameterlist( Ik ).resize( 1 );
        }
        tParameterlist( 3 ).resize( 2 );

        // ----------------------------------------------------------
        // linear solver algorithm

        if ( gUseBelosWithILUT )
        {
            tParameterlist( 7 )( 0 ) = moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::IFPACK );
            tParameterlist( 7 )( 0 ).set( "Convergence Tolerance", 1e-12 );
            tParameterlist( 7 )( 0 ).set( "ifpack_prec_type", "ILUT" );
            tParameterlist( 7 )( 0 ).set( "fact: drop tolerance", 1e-10 );
            tParameterlist( 7 )( 0 ).set( "fact: ilut level-of-fill", 25.0 );

            tParameterlist( 0 )( 0 ) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::BELOS_IMPL );
            tParameterlist( 0 )( 0 ).set( "preconditioners", "0" );
        }
        else
        {
            tParameterlist( 7 )( 0 ) = moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::NONE );
            tParameterlist( 0 )( 0 ) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL );
        }

        // ----------------------------------------------------------
        // linear solver

        tParameterlist( 1 )( 0 ) = moris::prm::create_linear_solver_parameter_list();

        // ----------------------------------------------------------

        if ( gHaveStaggeredFA )
        {
            // ----------------------------------------------------------
            // non-linear solver algorithms

            tParameterlist( 2 ).resize( 3 );

            // NEWTON solver algorithm
            tParameterlist( 2 )( 0 ) = moris::prm::create_nonlinear_algorithm_parameter_list();
            tParameterlist( 2 )( 0 ).set( "NLA_Solver_Implementation",  moris::NLA::NonlinearSolverType::NEWTON_SOLVER ) ;
            tParameterlist( 2 )( 0 ).set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
            tParameterlist( 2 )( 0 ).set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
            tParameterlist( 2 )( 0 ).set( "NLA_max_iter", tNLA_max_iter );
            tParameterlist( 2 )( 0 ).set( "NLA_combined_res_jac_assembly", true );

            // NLBGS solver algorithm
            // NOTE: maximum iterations is set to 1 since the second (structural problem) is linear and this saves time
            tParameterlist( 2 )( 1 )                                = moris::prm::create_nonlinear_algorithm_parameter_list();
            tParameterlist( 2 )( 1 ).set( "NLA_Solver_Implementation",  moris::NLA::NonlinearSolverType::NLBGS_SOLVER ) ;
            tParameterlist( 2 )( 1 ).set( "NLA_max_iter", 1 );

            // NEWTON solver algorithm for linear problems with only one iteration (i.e. in structural part and monolythic adjoint solve)
            tParameterlist( 2 )( 2 ) = moris::prm::create_nonlinear_algorithm_parameter_list();
            tParameterlist( 2 )( 2 ).set( "NLA_rel_res_norm_drop", 1.0e-7 );
            tParameterlist( 2 )( 2 ).set( "NLA_relaxation_parameter", 1.0 );
            tParameterlist( 2 )( 2 ).set( "NLA_max_iter", 1 );
            tParameterlist( 2 )( 2 ).set( "NLA_combined_res_jac_assembly", true );

            // ----------------------------------------------------------
            // non-linear solvers

            tParameterlist( 3 ).resize( 4 );

            // NEWTON solver for (linear) structural problem and adjoint
            tParameterlist( 3 )( 0 ) = moris::prm::create_nonlinear_solver_parameter_list();
            tParameterlist( 3 )( 0 ).set( "NLA_Solver_Implementation",  moris::NLA::NonlinearSolverType::NEWTON_SOLVER ) ;
            tParameterlist( 3 )( 0 ).set( "NLA_DofTypes", sStructDofTypes );
            tParameterlist( 3 )( 0 ).set( "NLA_Nonlinear_solver_algorithms", "2" );

            // NEWTON solver for non-linear thermal problem
            tParameterlist( 3 )( 1 ) = moris::prm::create_nonlinear_solver_parameter_list();
            tParameterlist( 3 )( 1 ).set( "NLA_Solver_Implementation",  moris::NLA::NonlinearSolverType::NEWTON_SOLVER ) ;
            tParameterlist( 3 )( 1 ).set( "NLA_DofTypes", "TEMP" );

            // NLBGS solver
            tParameterlist( 3 )( 2 ) = moris::prm::create_nonlinear_solver_parameter_list();
            tParameterlist( 3 )( 2 ).set( "NLA_Solver_Implementation",  moris::NLA::NonlinearSolverType::NLBGS_SOLVER ) ;
            tParameterlist( 3 )( 2 ).set( "NLA_Sub_Nonlinear_Solver", "1,0" );
            tParameterlist( 3 )( 2 ).set( "NLA_DofTypes", sAllDofTypes );
            tParameterlist( 3 )( 2 ).set( "NLA_Nonlinear_solver_algorithms", "1" );

            // NEWTON solver for separate monolythic adjoint solve
            tParameterlist( 3 )( 3 ) = moris::prm::create_nonlinear_solver_parameter_list();
            tParameterlist( 3 )( 3 ).set( "NLA_DofTypes", sAllDofTypes );
            tParameterlist( 3 )( 3 ).set( "NLA_Nonlinear_solver_algorithms", "2" );

            // ----------------------------------------------------------

            // MONOLYTHIC time solver
            tParameterlist( 4 )( 0 ) = moris::prm::create_time_solver_algorithm_parameter_list();
            tParameterlist( 4 )( 0 ).set( "TSA_Num_Time_Steps", tTSA_Num_Time_Steps );
            tParameterlist( 4 )( 0 ).set( "TSA_Time_Frame", tTSA_Time_Frame );
            tParameterlist( 4 )( 0 ).set( "TSA_Nonlinear_solver", 2 );

            // if SA is non-staggered, use separate monolythic non-linear solver for adjoint solve
            if ( !gHaveStaggeredSA )
                tParameterlist( 4 )( 0 ).set( "TSA_nonlinear_solver_for_adjoint_solve", 3 );

            // ----------------------------------------------------------

            tParameterlist( 5 )( 0 ) = moris::prm::create_time_solver_parameter_list();
            tParameterlist( 5 )( 0 ).set( "TSA_DofTypes", sAllDofTypes );
            tParameterlist( 5 )( 0 ).set( "TSA_Output_Indices", "0,1" );
            tParameterlist( 5 )( 0 ).set( "TSA_Output_Criteria", "Output_Criterion,Output_Criterion" );
            tParameterlist( 5 )( 0 ).set( "TSA_Initialize_Sol_Vec", "UX,0.0;UY,0.0;TEMP,1.0" );
            tParameterlist( 5 )( 0 ).set( "TSA_time_level_per_type", "UX," + sTLSD + ";UY," + sTLSD + ";TEMP,2" );
        }    // end: staggered solver case

        // ----------------------------------------------------------

        else    // monolythic solver
        {
            tParameterlist( 2 )( 0 ) = moris::prm::create_nonlinear_algorithm_parameter_list();
            tParameterlist( 2 )( 0 ).set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
            tParameterlist( 2 )( 0 ).set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
            tParameterlist( 2 )( 0 ).set( "NLA_max_iter", tNLA_max_iter );
            tParameterlist( 2 )( 0 ).set( "NLA_combined_res_jac_assembly", true );

            // ----------------------------------------------------------

            tParameterlist( 3 )( 0 ) = moris::prm::create_nonlinear_solver_parameter_list();
            tParameterlist( 3 )( 0 ).set( "NLA_DofTypes", "UX,UY;TEMP" );

            tParameterlist( 3 )( 1 ) = moris::prm::create_nonlinear_solver_parameter_list();
            tParameterlist( 3 )( 1 ).set( "NLA_DofTypes", "UX,UY;TEMP" );
            tParameterlist( 3 )( 1 ).set( "NLA_Nonlinear_solver_algorithms", "0" );

            // ----------------------------------------------------------

            tParameterlist( 4 )( 0 ) = moris::prm::create_time_solver_algorithm_parameter_list();
            tParameterlist( 4 )( 0 ).set( "TSA_Num_Time_Steps", tTSA_Num_Time_Steps );
            tParameterlist( 4 )( 0 ).set( "TSA_Time_Frame", tTSA_Time_Frame );

            // ----------------------------------------------------------

            tParameterlist( 5 )( 0 ) = moris::prm::create_time_solver_parameter_list();
            tParameterlist( 5 )( 0 ).set( "TSA_DofTypes", "UX,UY;TEMP" );
            tParameterlist( 5 )( 0 ).set( "TSA_Initialize_Sol_Vec", "UX,0.0;UY,0.0;TEMP,1.0" );
            tParameterlist( 5 )( 0 ).set( "TSA_Output_Indices", "0" );
            tParameterlist( 5 )( 0 ).set( "TSA_Output_Criteria", "Output_Criterion" );
            tParameterlist( 5 )( 0 ).set( "TSA_time_level_per_type", "UX," + sTLSD + ";UY," + sTLSD + ";TEMP,2" );
        }

        tParameterlist( 6 )( 0 ) = moris::prm::create_solver_warehouse_parameterlist();
        tParameterlist( 6 )( 0 ).set( "SOL_save_operator_to_matlab", "Mat.dat" );
    }

    void
    MSIParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_msi_parameter_list();
    }

    void
    VISParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
    {
        tParameterlist.resize( 3 );
        int tVisCounter = 0;

        tParameterlist( 0 ).push_back( prm::create_vis_parameter_list() );
        tParameterlist( 0 )( tVisCounter ).set( "Temp_Name", std::pair< std::string, std::string >( "./", "temp_material_1.exo" ) );
        tParameterlist( 0 )( tVisCounter ).set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName + "_material_1" ) );
        tParameterlist( 0 )( tVisCounter ).set( "Time_Offset", 100.0 );
        tParameterlist( 0 )( tVisCounter ).set( "Mesh_Type",  vis::VIS_Mesh_Type::STANDARD ) ;    // STANDARD_WITH_OVERLAP
        tParameterlist( 0 )( tVisCounter ).set( "Set_Names", tMat1Set );
        tParameterlist( 0 )( tVisCounter ).set( "Field_Names",
                "TEMP,UX,UY,"
                "MAX_DOF,STRAIN_ENERGY,VOLUME" );
        tParameterlist( 0 )( tVisCounter ).set( "Field_Type",
                "NODAL,NODAL,NODAL,"
                "GLOBAL,GLOBAL,GLOBAL" );
        tParameterlist( 0 )( tVisCounter ).set( "IQI_Names",
                "IQIBulkTEMP_1,IQIBulkDISPX_1,IQIBulkDISPY_1,"
                "IQIMaxTemp_1,IQIStrainEnergy_1,IQIVolume_1" );
        tParameterlist( 0 )( tVisCounter ).set( "Save_Frequency", 1 );
        tParameterlist( 0 )( tVisCounter ).set( "Output_Index", tVisCounter );
        tVisCounter++;

        tParameterlist( 0 ).push_back( prm::create_vis_parameter_list() );
        tParameterlist( 0 )( tVisCounter ).set( "Temp_Name", std::pair< std::string, std::string >( "./", "temp_material_2.exo" ) );
        tParameterlist( 0 )( tVisCounter ).set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName + "_material_2" ) );
        tParameterlist( 0 )( tVisCounter ).set( "Time_Offset", 100.0 );
        tParameterlist( 0 )( tVisCounter ).set( "Mesh_Type",  vis::VIS_Mesh_Type::STANDARD ) ;    // STANDARD_WITH_OVERLAP
        tParameterlist( 0 )( tVisCounter ).set( "Set_Names", tMat2Set );
        tParameterlist( 0 )( tVisCounter ).set( "Field_Names",
                "TEMP,UX,UY,"
                "MAX_DOF,STRAIN_ENERGY,VOLUME" );
        tParameterlist( 0 )( tVisCounter ).set( "Field_Type",
                "NODAL,NODAL,NODAL,"
                "GLOBAL,GLOBAL,GLOBAL" );
        tParameterlist( 0 )( tVisCounter ).set( "IQI_Names",
                "IQIBulkTEMP_2,IQIBulkDISPX_2,IQIBulkDISPY_2,"
                "IQIMaxTemp_2,IQIStrainEnergy_2,IQIVolume_2" );
        tParameterlist( 0 )( tVisCounter ).set( "Save_Frequency", 1 );
        tParameterlist( 0 )( tVisCounter ).set( "Output_Index", tVisCounter );
        tVisCounter++;
    }

    void
    MORISGENERALParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
    {
    }

    /* ------------------------------------------------------------------------ */
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
