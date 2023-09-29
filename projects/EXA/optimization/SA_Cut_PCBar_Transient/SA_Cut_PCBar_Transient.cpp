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
#include "fn_equal_to.hpp"

#include "AztecOO.h"

#ifdef  __cplusplus
extern "C"
{

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

    //include/exclude Nitsche Dirichlet boundary for debugging
    bool tHaveDirichlet = true;

    //include/exclude Ghost IWG for debugging
    bool tHaveGhost = true;

    //switch between running with phase change or without
    bool tIsPhaseChange = true;

    // Output Config --------------------------------------------------
    // set to true for vis output, set to false for sensitivity validation
    bool tOutputCriterion    = true;
    std::string tHDF5Path    = "SA_Cut_PCBar_Transient.hdf5";
    std::string tLibraryName = "SA_Cut_PCBar_Transient.so";
    std::string tOutputFile  = "SA_Cut_PCBar_Transient.exo";

    // Geometry Parameters --------------------------------------------
    moris::real tXlength = 0.1;
    moris::real tYlength = 0.5;
    moris::real tPcmFinRatio = 0.9;
    moris::real tDeltaRatio = 0.05;
    moris::real tPcmFinRatioMin = tPcmFinRatio - tDeltaRatio;
    moris::real tPcmFinRatioMax = tPcmFinRatio + tDeltaRatio;
    std::string tInterfacePos = std::to_string( tPcmFinRatio * tYlength );

    // Solver Configuration -------------------------------------------
    moris::sint tStep = 20;
    moris::real tTmax = 0.1;

    moris::real tNLARelResNormDrop = 5.0e-07;
    moris::real tNLARelaxationParameter = 0.9;
    moris::sint tNLAMaxIter = 25;

    // Sensitivity Analysis Parameters
    moris::real tFEMFdEpsilon =  1.0e-8;
    moris::uint tFEMFdScheme  =  static_cast< uint >( fem::FDScheme_Type::POINT_3_CENTRAL);

    // material parameters --------------------------------------------

    // conductor material
    std::string tDensityFin = "1.0";
    std::string tHeatCapFin = "5.0";
    std::string tThermConductFin = "50.0";

    // heat storage material
    std::string tDensityPCM = "1.0";
    std::string tHeatCapPCM = "50.0";
    std::string tThermConductPCM = "5.0";
    std::string tLatentHeatPCM = "500.0";
    std::string tPCTemp = "20.0";
    std::string tPCTempRange = "6.0";

    // initial & boundary conditions ----------------------------------
    std::string tHeatFlux = "1000.0";
    moris::real tNBClengthFraction = 0.1;
    moris::real tInitialTemp = 10.0;

    // IQI Configuration ----------------------------------------------
    std::string tMaxTempReference = "1.0";
    std::string tMaxTempExponent = "2.0";

    // Mesh sets ------------------------------------------------------

    // Bulk sets
    std::string tFinBulk = "HMR_dummy_n_p3,HMR_dummy_c_p3";
    std::string tPcmBulk = "HMR_dummy_n_p2,HMR_dummy_c_p2";
    std::string tTotalDomain = tFinBulk + "," + tPcmBulk;

    // Side sets
    std::string tFinPcmInterface  = "dbl_iside_p0_3_p1_2";
    std::string tFinNeumannInterface  = "SideSet_3_n_p3,SideSet_3_c_p3";
    std::string tPCMDirichletInterface = "SideSet_1_n_p2,SideSet_1_c_p2";

    // Ghost sets
    std::string tFinGhost = "ghost_p3";
    std::string tPcmGhost = "ghost_p1";

    // HMR parameters -------------------------------------------------
    //std::string tNumElemsPerDim = "2, 16";
    std::string tNumElemsPerDim = "1, 40";
    std::string tDomainDims = "0.16, 0.6";
    std::string tDomainOffset = "-0.0342356,-0.031345";

    //------------------------------------------------------------------------------
    //-------------------------------- FUNCTIONS -----------------------------------
    //------------------------------------------------------------------------------

    /* ------------------------------------------------------------------------ */
    // PROPERTY FUNCTIONS (incl. INITIAL & BOUNDARY CONDITIONS)
    /* ------------------------------------------------------------------------ */

    void Func_Initial_Temperature(
            moris::Matrix< moris::DDRMat >                 & aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
            moris::fem::Field_Interpolator_Manager         * aFIManager )
    {
        aPropMatrix = {{ tInitialTemp }};
    }

    // function for concentrating the heat load at the tip
    void Func_Neumann_BC(
            moris::Matrix< moris::DDRMat >                 & aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
            moris::fem::Field_Interpolator_Manager         * aFIManager )
    {
        // get coordinates
        moris::Matrix<DDRMat> tPosition = aFIManager->get_IG_geometry_interpolator()->valx();
        moris::real tYPosition = tPosition(1);

        moris::real tHeatFluxOn = 0.0;

        if (tYPosition >=  tYlength - tNBClengthFraction * tYlength)
            tHeatFluxOn = 1.0;

        // assume that heat flux is proportional to surface stagnation pressure
        aPropMatrix = tHeatFluxOn * aParameters( 0 );
    }

    // Constant function for properties
    void Func_Const(
            moris::Matrix< moris::DDRMat >                 & aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
            moris::fem::Field_Interpolator_Manager         * aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    /* ------------------------------------------------------------------------ */
    // DUMMY FUNCTIONS
    /* ------------------------------------------------------------------------ */

    bool Output_Criterion( moris::tsa::Time_Solver * aTimeSolver )
    {
        return tOutputCriterion;
    }

    moris::Matrix< DDRMat > Func_Dummy_Sensitivity(
            const moris::Matrix< DDRMat >     & aCoordinates,
            const moris::Cell< moris::real* > & aGeometryParameters )
            {
        moris::Matrix< DDRMat > aReturnValue = {{0.0}};
        return aReturnValue;
            }

    /* ------------------------------------------------------------------------ */
    // FOR SWEEP
    /* ------------------------------------------------------------------------ */

    moris::Matrix< moris::DDSMat > get_constraint_types()
            {
        Matrix<DDSMat> tConstraintTypes( 1, 1, 1 );
        return tConstraintTypes;
            }

    moris::Matrix< moris::DDRMat > compute_objectives(Matrix<DDRMat> aADVs, moris::Matrix< moris::DDRMat > aCriteria)
            {
        moris::Matrix< moris::DDRMat > tObjectives( 1, 1, aCriteria(0) );
        return tObjectives;
            }

    moris::Matrix< moris::DDRMat > compute_constraints(Matrix<DDRMat> aADVs, moris::Matrix< moris::DDRMat > aCriteria)
            {
        moris::Matrix< moris::DDRMat > tConstraints( 1, 1, aCriteria(1) );
        return tConstraints;
            }

    moris::Matrix< moris::DDRMat > compute_dobjective_dadv(Matrix<DDRMat> aADVs, moris::Matrix< moris::DDRMat > aCriteria)
            {
        moris::Matrix< moris::DDRMat > tDObjectiveDADV( 1, aADVs.length(), 0.0 );
        return tDObjectiveDADV;
            }

    moris::Matrix< moris::DDRMat > compute_dobjective_dcriteria(moris::Matrix< moris::DDRMat > aADVs, moris::Matrix< moris::DDRMat > aCriteria)
            {
        moris::Matrix< moris::DDRMat > tDObjectiveDCriteria( 1, 2, 0.0 );
        tDObjectiveDCriteria(0) = 1.0;
        return tDObjectiveDCriteria;
            }

    moris::Matrix< moris::DDRMat > compute_dconstraint_dadv(moris::Matrix< moris::DDRMat > aADVs, moris::Matrix< moris::DDRMat > aCriteria)
            {
        moris::Matrix< moris::DDRMat > tDConstraintDADV( 1, aADVs.length(), 0.0 );
        return tDConstraintDADV;
            }

    moris::Matrix< moris::DDRMat > compute_dconstraint_dcriteria(moris::Matrix< moris::DDRMat > aADVs, moris::Matrix< moris::DDRMat > aCriteria)
            {
        moris::Matrix< moris::DDRMat > tDConstraintDCriteria( 1, 2, 0.0 );
        tDConstraintDCriteria(1) = 1.0;
        return tDConstraintDCriteria;
            }

    /* ------------------------------------------------------------------------ */
    // PARAMETER LISTS
    /* ------------------------------------------------------------------------ */

    void HMRParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_hmr_parameter_list();

        tParameterlist( 0 )( 0 ).set( "number_of_elements_per_dimension", tNumElemsPerDim );
        tParameterlist( 0 )( 0 ).set( "domain_dimensions",                tDomainDims );
        tParameterlist( 0 )( 0 ).set( "domain_offset",                    tDomainOffset );
        tParameterlist( 0 )( 0 ).set( "domain_sidesets",                  "1,2,3,4");
        tParameterlist( 0 )( 0 ).set( "lagrange_output_meshes",           "0");

        tParameterlist( 0 )( 0 ).set( "lagrange_orders",   "1" );
        tParameterlist( 0 )( 0 ).set( "bspline_orders",    "1" );
        tParameterlist( 0 )( 0 ).set( "lagrange_pattern",  "0" );
        tParameterlist( 0 )( 0 ).set( "bspline_pattern",   "0" );

        tParameterlist( 0 )( 0 ).set( "lagrange_to_bspline", "0") ;

        tParameterlist( 0 )( 0 ).set( "truncate_bsplines",  1 );
        tParameterlist( 0 )( 0 ).set( "refinement_buffer",  0 );
        tParameterlist( 0 )( 0 ).set( "staircase_buffer",   0 );
        tParameterlist( 0 )( 0 ).set( "initial_refinement", "0" );
        tParameterlist( 0 )( 0 ).set( "initial_refinement_pattern", "0" );

        tParameterlist( 0 )( 0 ).set( "use_multigrid",  0 );
        tParameterlist( 0 )( 0 ).set( "severity_level", 0 );

        tParameterlist( 0 )( 0 ).set( "adaptive_refinement_level", 0 );
    }

    void OPTParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 3 );
        tParameterlist( 0 ).resize( 1 );
        tParameterlist( 1 ).resize( 0 );
        tParameterlist( 2 ).resize( 1 );

        tParameterlist(0)(0) = moris::prm::create_opt_problem_parameter_list();
        tParameterlist(0)(0).set("is_optimization_problem", true );
        tParameterlist(0)(0).set("problem", "user_defined" );
        tParameterlist(0)(0).set("library", tLibraryName );

        tParameterlist(2)(0) = moris::prm::create_sweep_parameter_list();
        tParameterlist(2)(0).set("print", true);
        tParameterlist(2)(0).set("hdf5_path", tHDF5Path);
        tParameterlist(2)(0).set("num_evaluations_per_adv", "1");
        tParameterlist(2)(0).set("include_bounds", false);
        tParameterlist(2)(0).set("finite_difference_type", "all");
        tParameterlist(2)(0).set("finite_difference_epsilons", "1E-6");
    }

    void XTKParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_xtk_parameter_list();
        tParameterlist( 0 )( 0 ).set( "decompose",                 true );
        tParameterlist( 0 )( 0 ).set( "decomposition_type",        "conformal") ;
        tParameterlist( 0 )( 0 ).set( "enrich",                    true );
        tParameterlist( 0 )( 0 ).set( "basis_rank",                "bspline") ;
        tParameterlist( 0 )( 0 ).set( "enrich_mesh_indices",       "0") ;
        tParameterlist( 0 )( 0 ).set( "ghost_stab",                true );
        tParameterlist( 0 )( 0 ).set( "multigrid",                 false );
        tParameterlist( 0 )( 0 ).set( "verbose",                   true );
        tParameterlist( 0 )( 0 ).set( "print_enriched_ig_mesh",    true );
        tParameterlist( 0 )( 0 ).set( "exodus_output_XTK_ig_mesh", true );
        tParameterlist( 0 )( 0 ).set( "high_to_low_dbl_side_sets", true );
    }

    void GENParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 3 );
        tParameterlist( 0 ).resize( 1 );

        // Main GEN parameter list
        tParameterlist( 0 )( 0 ) = prm::create_gen_parameter_list();
        tParameterlist( 0 )( 0 ).set("initial_advs"      , std::to_string(tYlength * tPcmFinRatioMin));
        tParameterlist( 0 )( 0 ).set("lower_bounds"      , std::to_string(tYlength * tPcmFinRatioMin));
        tParameterlist( 0 )( 0 ).set("upper_bounds"      , std::to_string(tYlength * tPcmFinRatioMax));
        tParameterlist( 0 )( 0 ).set("IQI_types"         , "IQIMaxTemp,IQIBulkVolume");

        // init geometry counter
        uint tGeoCounter = 0;

        // Geometry parameter lists
        tParameterlist( 1 ).push_back( prm::create_level_set_geometry_parameter_list() );
        tParameterlist( 1 )( tGeoCounter ).set( "field_type", "plane");
        tParameterlist( 1 )( tGeoCounter ).set( "constant_parameters", "-1.0, 0.0, 1.0, 0.0");
        tGeoCounter++;

        tParameterlist( 1 ).push_back( prm::create_level_set_geometry_parameter_list() );
        tParameterlist( 1 )( tGeoCounter ).set( "field_type", "plane");
        tParameterlist( 1 )( tGeoCounter ).set( "field_variable_indices", "1");
        tParameterlist( 1 )( tGeoCounter ).set( "adv_indices", "0");
        tParameterlist( 1 )( tGeoCounter ).set( "constant_parameters", "0.0, 0.0, 1.0");
        tGeoCounter++;
    }

    void FEMParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterList )
    {
        // create a cell of cell of parameter list for fem
        tParameterList.resize( 8 );

        //------------------------------------------------------------------------------
        // init property counter
        uint tPropCounter = 0;

        // Density of conductor material
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ) = prm::create_property_parameter_list();
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropFinDensity") ;
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      tDensityFin );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const") ;
        tPropCounter++;

        // Thermal conductivity of conductor material
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropFinConductivity") ;
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      tThermConductFin );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const") ;
        tPropCounter++;

        // Heat capacity of conductor material
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropFinHeatCapacity") ;
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      tHeatCapFin );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const") ;
        tPropCounter++;

        // Density of storage material
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ) = prm::create_property_parameter_list();
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropPcmDensity") ;
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      tDensityPCM );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const") ;
        tPropCounter++;

        // Thermal conductivity of storage material
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropPcmConductivity") ;
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      tThermConductPCM );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const") ;
        tPropCounter++;

        // Heat capacity of storage material
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropPcmHeatCapacity") ;
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      tHeatCapPCM );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const") ;
        tPropCounter++;

        // Latent heat capacity of storage material
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropLatentHeat") ;
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      tLatentHeatPCM );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const") ;
        tPropCounter++;

        // Melt temperature of storage material
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropPCTemp") ;
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      tPCTemp );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const") ;
        tPropCounter++;

        // Phase change function
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropPhaseState") ;
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      "2.0") ;
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const") ;
        tPropCounter++;

        // Melting range of storage material
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropPCconst") ;
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      tPCTempRange );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const") ;
        tPropCounter++;

        // Neumann BC
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropNeumannBC") ;
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      tHeatFlux );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const") ;
        tPropCounter++;

        // Initial condition (Temperature)
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropInitialCondition") ;
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Initial_Temperature") ;
        tPropCounter++;

        // For Time solver
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropWeightCurrent") ;
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      "100.0") ;
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const") ;
        tPropCounter++;
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropWeightPrevious") ;
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      "100.0") ;
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const") ;
        tPropCounter++;

        // Dummy latent heat for non-pc material
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropDummyLatentHeat") ;
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      "0.0") ;
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const") ;
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropDummyPCTemp") ;
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      "10000.0") ;
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const") ;
        tPropCounter++;

        //------------------------------------------------------------------------------

        // init CM counter
        uint tCMCounter = 0;

        if (tIsPhaseChange)
        {
            // constitutive model for phase change material
            tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
            tParameterList( 1 )( tCMCounter ).set( "constitutive_name", "CMDiffusionPcm") ;
            tParameterList( 1 )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::DIFF_LIN_ISO_PC ) );
            tParameterList( 1 )( tCMCounter ).set( "dof_dependencies",  std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
            tParameterList( 1 )( tCMCounter ).set( "properties",
                    "PropPcmConductivity,Conductivity;"
                    "PropPcmDensity,Density;"
                    "PropPcmHeatCapacity,HeatCapacity;"
                    "PropLatentHeat,LatentHeat;"
                    "PropPCTemp,PCTemp;"
                    "PropPhaseState,PhaseStateFunction;"
                    "PropPCconst,PhaseChangeConst") ;
            tCMCounter++;
        }
        else // no phase change
        {
            // constitutive model for thermal storage material
            tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
            tParameterList( 1 )( tCMCounter ).set( "constitutive_name", "CMDiffusionPcm") ;
            tParameterList( 1 )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::DIFF_LIN_ISO ) );
            tParameterList( 1 )( tCMCounter ).set( "dof_dependencies",  std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
            tParameterList( 1 )( tCMCounter ).set( "properties",
                    "PropPcmConductivity,Conductivity;"
                    "PropPcmDensity,Density;"
                    "PropPcmHeatCapacity,HeatCapacity") ;
            tCMCounter++;
        }

        // constitutive model for thermal conductor material
        tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_name", "CMDiffusionFin") ;
        tParameterList( 1 )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::DIFF_LIN_ISO ) );
        tParameterList( 1 )( tCMCounter ).set( "dof_dependencies",  std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( 1 )( tCMCounter ).set( "properties",
                "PropFinConductivity,Conductivity;"
                "PropFinDensity,Density;"
                "PropFinHeatCapacity,HeatCapacity") ;
        tCMCounter++;

        //------------------------------------------------------------------------------

        // init SP counter
        uint tSPCounter = 0;

        // Ghost parameter for fin
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name",      "SPGPTempFin") ;
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type",      static_cast< uint >( fem::Stabilization_Type::GHOST_DISPL ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters",     "0.01") ;
        tParameterList( 2 )( tSPCounter ).set( "leader_properties",       "PropFinConductivity,Material") ;
        tSPCounter++;

        // Ghost parameter for PCM
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name",      "SPGPTempPcm") ;
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type",      static_cast< uint >( fem::Stabilization_Type::GHOST_DISPL ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters",     "0.01") ;
        tParameterList( 2 )( tSPCounter ).set( "leader_properties",       "PropPcmConductivity,Material") ;
        tSPCounter++;

        // GGLS parameter for fin
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name",      "SPGGLSDiffusionFin") ;
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type",      static_cast< uint >( fem::Stabilization_Type::GGLS_DIFFUSION ) );
        tParameterList( 2 )( tSPCounter ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties",
                "PropFinConductivity, Conductivity;"
                "PropFinDensity     , Density;"
                "PropFinHeatCapacity, HeatCapacity;"
                "PropDummyLatentHeat, LatentHeat;"
                "PropDummyPCTemp    , PCTemp;"
                "PropPhaseState     , PhaseStateFunction;"
                "PropPCconst        , PhaseChangeConst")    ;
        tSPCounter++;

        if (tIsPhaseChange)
        {
            // GGLS parameter for phase change material
            tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
            tParameterList( 2 )( tSPCounter ).set( "stabilization_name",      "SPGGLSDiffusionPcm") ;
            tParameterList( 2 )( tSPCounter ).set( "stabilization_type",      static_cast< uint >( fem::Stabilization_Type::GGLS_DIFFUSION ) );
            tParameterList( 2 )( tSPCounter ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
            tParameterList( 2 )( tSPCounter ).set( "leader_properties",
                    "PropPcmConductivity, Conductivity;"
                    "PropPcmDensity     , Density;"
                    "PropPcmHeatCapacity, HeatCapacity;"
                    "PropLatentHeat     , LatentHeat;"
                    "PropPCTemp         , PCTemp;"
                    "PropPhaseState     , PhaseStateFunction;"
                    "PropPCconst        , PhaseChangeConst")    ;
            tSPCounter++;
        }
        else // no phase change
        {
            // GGLS parameter for thermal storage material
            tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
            tParameterList( 2 )( tSPCounter ).set( "stabilization_name",      "SPGGLSDiffusionPcm") ;
            tParameterList( 2 )( tSPCounter ).set( "stabilization_type",      static_cast< uint >( fem::Stabilization_Type::GGLS_DIFFUSION ) );
            tParameterList( 2 )( tSPCounter ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
            tParameterList( 2 )( tSPCounter ).set( "leader_properties",
                    "PropPcmConductivity, Conductivity;"
                    "PropPcmDensity     , Density;"
                    "PropPcmHeatCapacity, HeatCapacity;"
                    "PropDummyLatentHeat, LatentHeat;"
                    "PropDummyPCTemp    , PCTemp;"
                    "PropPhaseState     , PhaseStateFunction;"
                    "PropPCconst        , PhaseChangeConst")    ;
            tSPCounter++;
        }

        // Dirichlet SP
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name",  "SPInterfaceNitsche") ;
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type",  static_cast< uint >( fem::Stabilization_Type::NITSCHE_INTERFACE ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", "100.0") ;
        tParameterList( 2 )( tSPCounter ).set( "leader_properties",   "PropFinConductivity,Material") ;
        tParameterList( 2 )( tSPCounter ).set( "follower_properties",    "PropPcmConductivity,Material") ;
        tSPCounter++;

        //------------------------------------------------------------------------------
        // init IWG counter
        uint tIWGCounter = 0;

        // Bulk IWG for conductor material
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWGDiffusionFinBulk") ;
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::SPATIALDIFF_BULK ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               "TEMP") ;
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies",    "TEMP") ;
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMDiffusionFin,Diffusion") ;
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   "SPGGLSDiffusionFin,GGLSParam") ;
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             tFinBulk );
        tIWGCounter++;

        // Bulk IWG for thermal storage material
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWGDiffusionPcmBulk") ;
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::SPATIALDIFF_BULK ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               "TEMP") ;
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies",    "TEMP") ;
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMDiffusionPcm,Diffusion") ;
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   "SPGGLSDiffusionPcm,GGLSParam") ;
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             tPcmBulk );
        tIWGCounter++;

        if ( tHaveDirichlet )
        {
            // Interface Dirichlet BC
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWGInterfaceFinPcm") ;
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               "TEMP") ;
            tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies",    "TEMP") ;
            tParameterList( 3 )( tIWGCounter ).set( "follower_dof_dependencies",     "TEMP") ;
            tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMDiffusionFin,Diffusion") ;
            tParameterList( 3 )( tIWGCounter ).set( "follower_constitutive_models",  "CMDiffusionPcm,Diffusion") ;
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",
                    "SPInterfaceNitsche ,NitscheInterface")   ;
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             tFinPcmInterface );
            tIWGCounter++;
        }

        // Imposed heat flux
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWGInletFlux") ;
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::SPATIALDIFF_NEUMANN ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               "TEMP") ;
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies",    "TEMP") ;
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties",          "PropNeumannBC,Neumann") ;
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             tFinNeumannInterface );
        tIWGCounter++;

        if ( tHaveGhost )
        {
            // Fin Ghost
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWGFinGhost" );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::GHOST_NORMAL_FIELD ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               "TEMP" );
            tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies",    "TEMP" );
            tParameterList( 3 )( tIWGCounter ).set( "follower_dof_dependencies",     "TEMP" );
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   "SPGPTempFin,GhostSP" );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             tFinGhost );
            tIWGCounter++;

            // PCM Ghost
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWGPcmGhost" );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::GHOST_NORMAL_FIELD ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               "TEMP" );
            tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies",    "TEMP" );
            tParameterList( 3 )( tIWGCounter ).set( "follower_dof_dependencies",     "TEMP" );
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   "SPGPTempPcm,GhostSP" );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             tPcmGhost );
            tIWGCounter++;
        }

        // Time Continuity
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWGTimeContinuityTemp") ;
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::TIME_CONTINUITY_DOF ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               "TEMP") ;
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies",    "TEMP") ;
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties",
                "PropWeightCurrent,WeightCurrent;"
                "PropWeightPrevious  ,WeightPrevious;"
                "PropInitialCondition,InitialCondition") ;
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             tTotalDomain );
        tParameterList( 3 )( tIWGCounter ).set( "time_continuity",            true );
        tIWGCounter++;

        //------------------------------------------------------------------------------
        // init IQI counter
        uint tIQICounter = 0;

        // IQI - Nodal Temperature Field
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name",                   "IQIBulkTEMP") ;
        tParameterList( 4 )( tIQICounter ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::DOF ) );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity",               "TEMP") ;
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies",    "TEMP") ;
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index",      0 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names",             tTotalDomain );
        tIQICounter++;

        // Volume IQI - Total Volume
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name",                   "IQIBulkVolume") ;
        tParameterList( 4 )( tIQICounter ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::VOLUME ) );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies",    "TEMP") ;
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names",             tPcmBulk );
        tIQICounter++;

        // Max Temperature IQI
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name",                   "IQIMaxTemp") ;
        tParameterList( 4 )( tIQICounter ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::MAX_DOF ) );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity",               "TEMP") ;
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies",    "TEMP") ;
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index",      0 );
        tParameterList( 4 )( tIQICounter ).set( "function_parameters",
                        tMaxTempReference + "/" + tMaxTempExponent ) ;
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names",             tTotalDomain );
        tIQICounter++;

        // create computation  parameter list
        tParameterList( 5 ).resize( 1 );
        tParameterList( 5 )( 0 ) = prm::create_computation_parameter_list();
        tParameterList( 5 )( 0 ).set( "finite_difference_scheme",            tFEMFdScheme  );
        tParameterList( 5 )( 0 ).set( "finite_difference_perturbation_size", tFEMFdEpsilon );
    }

    void SOLParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 7 );
        for( uint Ik = 0; Ik < 7; Ik ++)
        {
            tParameterlist( Ik ).resize( 1 );
        }

        tParameterlist( 0 )( 0 ) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::BELOS_IMPL );
        tParameterlist( 0 )( 0 ).set( "ifpack_prec_type", "ILU");

        tParameterlist( 1 )( 0 ) = moris::prm::create_linear_solver_parameter_list();

        tParameterlist( 2 ).resize(2);

        // for forward analysis
        tParameterlist( 2 )( 0 ) = moris::prm::create_nonlinear_algorithm_parameter_list();
        tParameterlist( 2 )( 0 ).set("NLA_rel_res_norm_drop", tNLARelResNormDrop );
        tParameterlist( 2 )( 0 ).set("NLA_relaxation_parameter", tNLARelaxationParameter );
        tParameterlist( 2 )( 0 ).set("NLA_relaxation_strategy",  static_cast< uint >( sol::SolverRelaxationType::InvResNormAdaptive ) );
        tParameterlist( 2 )( 0 ).set("NLA_relaxation_parameter", 0.5 );
        tParameterlist( 2 )( 0 ).set("NLA_relaxation_damping",   0.5 );
        tParameterlist( 2 )( 0 ).set("NLA_max_iter", tNLAMaxIter );
        tParameterlist( 2 )( 0 ).set("NLA_combined_res_jac_assembly", false );

        // for adjoint analysis
        tParameterlist( 2 )( 1 ) = moris::prm::create_nonlinear_algorithm_parameter_list();
        tParameterlist( 2 )( 1 ).set("NLA_rel_res_norm_drop", tNLARelResNormDrop );
        tParameterlist( 2 )( 1 ).set("NLA_relaxation_parameter", 1.0 );
        tParameterlist( 2 )( 1 ).set("NLA_max_iter", 1 );
        tParameterlist( 2 )( 1 ).set("NLA_combined_res_jac_assembly", false );

        tParameterlist( 3 ).resize(2);

        // for forward analysis
        tParameterlist( 3 )( 0 ) = moris::prm::create_nonlinear_solver_parameter_list();
        tParameterlist( 3 )( 0 ).set("NLA_DofTypes"      , "TEMP") ;
        tParameterlist( 3 )( 0 ).set("NLA_Nonlinear_solver_algorithms", "0");

        // for adjoint analysis
        tParameterlist( 3 )( 1 ) = moris::prm::create_nonlinear_solver_parameter_list();
        tParameterlist( 3 )( 1 ).set("NLA_DofTypes"      , "TEMP") ;
        tParameterlist( 3 )( 1 ).set("NLA_Nonlinear_solver_algorithms", "1");

        tParameterlist( 4 )( 0 ) = moris::prm::create_time_solver_algorithm_parameter_list();
        tParameterlist( 4 )( 0 ).set("TSA_Num_Time_Steps",     tStep );
        tParameterlist( 4 )( 0 ).set("TSA_Time_Frame",         tTmax );
        tParameterlist( 4 )( 0 ).set("TSA_Nonlinear_solver",                   0 );   // for forward analysis
        tParameterlist( 4 )( 0 ).set("TSA_nonlinear_solver_for_adjoint_solve", 1 );   // for adjoint analysis

        tParameterlist( 5 )( 0 ) = moris::prm::create_time_solver_parameter_list();
        tParameterlist( 5 )( 0 ).set("TSA_DofTypes",            "TEMP") ;
        tParameterlist( 5 )( 0 ).set("TSA_Initialize_Sol_Vec",  "TEMP,0.0") ;
        tParameterlist( 5 )( 0 ).set("TSA_Output_Indices",      "0") ;
        tParameterlist( 5 )( 0 ).set("TSA_Output_Criteria",      "Output_Criterion") ;
        tParameterlist( 5 )( 0 ).set("TSA_time_level_per_type", "TEMP,2") ;

        tParameterlist( 6 )( 0 ) = moris::prm::create_solver_warehouse_parameterlist();
    }

    void MSIParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_msi_parameter_list();
    }

    void VISParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_vis_parameter_list();
        tParameterlist( 0 )( 0 ).set( "File_Name"  , std::pair< std::string, std::string >( "./", tOutputFile ) );
        tParameterlist( 0 )( 0 ).set( "Mesh_Type"  , static_cast< uint >( vis::VIS_Mesh_Type::STANDARD ) );
        tParameterlist( 0 )( 0 ).set( "Set_Names"  , tTotalDomain );
        tParameterlist( 0 )( 0 ).set( "Field_Names",  "TEMP,MAX_DOF" ) ;
        tParameterlist( 0 )( 0 ).set( "Field_Type" ,  "NODAL,GLOBAL" ) ;
        tParameterlist( 0 )( 0 ).set( "IQI_Names",  "IQIBulkTEMP,IQIMaxTemp" ) ;
    }

    void MORISGENERALParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {

    }

    //------------------------------------------------------------------------------
}

//------------------------------------------------------------------------------
#ifdef  __cplusplus
}
#endif

