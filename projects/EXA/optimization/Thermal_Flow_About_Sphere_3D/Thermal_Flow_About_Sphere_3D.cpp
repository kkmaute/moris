/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Channel3D.cpp
 *
 */

#include <string>
#include <iostream>
#include <sstream>

#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "cl_Bitset.hpp"
#include "linalg_typedefs.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_MSI_Equation_Object.hpp"
#include "cl_TSA_Time_Solver.hpp"
#include "cl_DLA_Solver_Interface.hpp"
#include "cl_DLA_Linear_Solver_Aztec.hpp"
#include "parameters.hpp"
#include "fn_equal_to.hpp"

#include "AztecOO.h"

//---------------------------------------------------------------

#ifdef __cplusplus
extern "C" {
#endif
//------------------------------------------------------------------------------
namespace moris
{
    std::string
    moris_to_string( real tValue )
    {
        std::ostringstream streamObj;

        // Set precision
        streamObj << std::scientific;
        streamObj << std::setprecision( 15 );

        // Add value to stream
        streamObj << tValue;

        // Get string from output string stream
        return streamObj.str();
    }

    uint gInterpolationOrder = 2;

    /* ------------------------------------------------------------------------ */
    // set optimization restart iteration (if 0 no restart)
    sint tRestartId = 0;

    /* ------------------------------------------------------------------------ */

    // inlet pressure
    real tInPres = 0.001 * 101325.0;    // N/m2  - 0.6 atm

    // reference values
    real presref = tInPres;    // reference pressure
    real lenref  = 1e-3;       // reference length
    real rhoref  = 0.3289;     // reference density
    real tempref = 1.0;        // reference temperature

    // derived reference values
    real velref  = std::sqrt( presref / rhoref );
    real timeref = lenref / velref;
    real massref = rhoref * std::pow( lenref, 3.0 );

    // reynolds number based on reference values and fluid properties see below
    real reynolds = lenref * velref * 0.3289 / 43.32e-6;

    // scaling parameters: input: m -> used for computation: mm
    real tLengthScale = 1.0 / lenref;     // 1 m  = 1.0/lref;
    real tTimeScale   = 1.0 / timeref;    // 1 s  = 1.0 s
    real tMassScale   = 1.0 / massref;    // 1 kg = 1.0 kg
    real tTempScale   = 1.0 / tempref;    // 1 K  = 1.0 K

    real tForceScale    = tMassScale * tLengthScale / tTimeScale / tTimeScale;
    real tPressureScale = tForceScale / tLengthScale / tLengthScale;
    real tEnergyScale   = tForceScale * tLengthScale;
    real tPowerScale    = tEnergyScale / tTimeScale;
    real tDensityScale  = tMassScale / tLengthScale / tLengthScale / tLengthScale;

    /* ------------------------------------------------------------------------ */
    // sphere geometry parameters

    real tChannelLength = 5.00e-3 * tLengthScale;    // x-direction
    real tChannelHeight = 1.00e-3 * tLengthScale;    // y-direction
    real tChannelWidth  = 1.00e-3 * tLengthScale;    // z-direction

    real tSpherePosX = 7.50e-4 * tLengthScale;
    real tSpherePosY = 5.00e-4 * tLengthScale;
    real tSpherePosZ = 5.00e-4 * tLengthScale;

    real tSphereRadius = 2.25e-4 * tLengthScale;
    real tSphereExpon  = 2.0;

    /* ------------------------------------------------------------------------ */
    // basic geometry parameters
    real tDimX = 5.25e-3 * tLengthScale;
    real tDimY = 1.25e-3 * tLengthScale;
    real tDimZ = 1.25e-3 * tLengthScale;

    real tOffsetX = ( tChannelLength - tDimX ) / 2.0;
    real tOffsetY = ( tChannelHeight - tDimY ) / 2.0;
    real tOffsetZ = ( tChannelWidth - tDimZ ) / 2.0;

    real tApproxEleSize = 2e-4 * tLengthScale;

    /* ------------------------------------------------------------------------ */
    // background mesh
    std::string tNumElemX = moris_to_string( std::ceil( tDimX / tApproxEleSize ) );
    std::string tNumElemY = moris_to_string( std::ceil( tDimY / tApproxEleSize ) );
    std::string tNumElemZ = moris_to_string( std::ceil( tDimZ / tApproxEleSize ) );

    std::string tDomainDimX = moris_to_string( tDimX );
    std::string tDomainDimY = moris_to_string( tDimY );
    std::string tDomainDimZ = moris_to_string( tDimZ );

    std::string tDomainOffX = moris_to_string( tOffsetX );
    std::string tDomainOffY = moris_to_string( tOffsetY );
    std::string tDomainOffZ = moris_to_string( tOffsetZ );

    std::string tNumElemsPerDim = tNumElemX + "," + tNumElemY + "," + tNumElemZ;
    std::string tDomainDims     = tDomainDimX + "," + tDomainDimY + "," + tDomainDimZ;
    std::string tDomainOffset   = tDomainOffX + "," + tDomainOffY + "," + tDomainOffZ;
    std::string tDomainSidesets = "1,2,3,4,5,6";

    int tLevelsetOrder = gInterpolationOrder;
    int tDispOrder     = gInterpolationOrder;

    int tLevelsetInitialRef = 0;
    int tDispInitialRef     = 0;
    int tRefineBuffer       = 1;

    // note: pattern 0 - Levelset field  pattern 1 - displacement field
    std::string tLagrangeOrder   = std::to_string( std::max( tLevelsetOrder, tDispOrder ) );
    std::string tBsplineOrder    = std::to_string( tLevelsetOrder ) + "," + std::to_string( tDispOrder );
    std::string tInitialRef      = std::to_string( tLevelsetInitialRef ) + "," + std::to_string( tDispInitialRef );
    std::string tLagrangePattern = tLevelsetInitialRef > tDispInitialRef ? "0" : "1";

    uint tInterfaceRefinementSphere = 0;
    uint tInterfaceRefinementWalls  = 0;

    moris::real tElementEdgeLength = tApproxEleSize / ( std::pow( 2, tDispInitialRef ) );
    moris::real tGeoShift          = 0.1 * tElementEdgeLength;

    // Bspline limit
    moris::real tBSplineLimit = 5.0 * tElementEdgeLength;

    /* ------------------------------------------------------------------------ */
    // loading parameters
    real tInTemp  = 0.0;        // K     - 500 K
    real tHeatFlx = 250.0e4;    // W/m^2 - 250 W/cm^2

    std::string tInletPressure      = moris_to_string( tInPres * tPressureScale );
    std::string tInletTemperature   = moris_to_string( tInTemp * tTempScale );
    std::string tVolumetricHeatLoad = moris_to_string( tHeatFlx * tPowerScale / tLengthScale / tLengthScale / tLengthScale );

    /* ------------------------------------------------------------------------ */
    // material parameters

    std::string tFluidViscosity      = moris_to_string( 43.32e-6 * tPressureScale * tTimeScale );                //  N s/m2  air at 800 C
    std::string tFluidDensity        = moris_to_string( 0.3289 * tDensityScale );                                // kg/m3   air at 800 C and 1 atm
    std::string tFluidCapacity       = moris_to_string( 1.210e3 * tEnergyScale / tMassScale / tTempScale );      // J/kg K  cp air at 800 C and 1 atm
    std::string tFluidConductivity   = moris_to_string( 76.26e-3 * tPowerScale / tLengthScale / tTempScale );    // W/m K   air at 800 C
    std::string tFluidPressureSpring = "1e-6";

    std::string tSolidDensity      = moris_to_string( 7.0e3 * tDensityScale );                              // kg/m3
    std::string tSolidCapacity     = moris_to_string( 200.0 * tEnergyScale / tMassScale / tTempScale );     // J/kg*K
    std::string tSolidConductivity = moris_to_string( 150.0 * tPowerScale / tLengthScale / tTempScale );    // W/m K ;

    /* ------------------------------------------------------------------------ */
    // File names
    std::string tName          = "Channel3D";
    std::string tExoFile       = tName + ".exo";
    std::string tSoFile        = tName + ".so";
    std::string tHdf5File      = tName + "_SEN.hdf5";
    std::string tGENOutputFile = tName + "_GEN.exo";

    // FD in adjoint
    real tFEMFdEpsilon = 1.0e-4;

    // FD in sweep
    std::string tSweepFdEpsilon = "1.0e-4";

    // Eleement inverse estimate
    std::string tCInv = gInterpolationOrder == 2 ? "60.0" : "36.0";

    /* ------------------------------------------------------------------------ */

    // number of constraints (here: number of design criteria)
    moris::uint tNumConstraints = 9;

    // max dof IQI
    std::string tRefTemp         = moris_to_string( 250.0 * tTempScale );    // about 0.5 of inlet temperature
    std::string tMaxTempExponent = "8.0";

    /* ------------------------------------------------------------------------ */
    // solver parameters

    moris::real tNLA_rel_res_norm_drop    = 1e-8;
    moris::real tNLA_relaxation_parameter = 1.0;
    int         tNLA_max_iter             = 20;

    int         tTSA_Num_Time_Steps = 1;
    moris::real tTSA_Time_Frame     = 1.0;

    /* ------------------------------------------------------------------------ */
    // Output Config

    std::string tFluidPhase = "HMR_dummy_n_p0,HMR_dummy_c_p0";
    std::string tSolidPhase = "HMR_dummy_n_p1,HMR_dummy_c_p1";
    std::string tAllPhases  = tFluidPhase + "," + tSolidPhase;

    bool
    Output_Criterion( moris::tsa::Time_Solver* aTimeSolver )
    {
        return true;
    }

    /* ------------------------------------------------------------------------ */
    // Flag for turning on/off ghost stabilization
    bool tUseGhost = true;

    /* ------------------------------------------------------------------------ */
    // geometry parameters & LS functions

    // Minimum level set value
    moris::real tMinLevs = 1.0e-8;

    moris::real
    Func_Sphere(
            const moris::Matrix< DDRMat >& aCoordinates,
            const Vector< real >&          aGeometryParameters )
    {
        // get coordinates
        real tX = aCoordinates( 0 );
        real tY = aCoordinates( 1 );
        real tZ = aCoordinates( 2 );

        real tRadius = aGeometryParameters( 0 );

        real tReturnValue = tRadius - std::pow( std::pow( tX - tSpherePosX, tSphereExpon ) + std::pow( tY - tSpherePosY, tSphereExpon ) + std::pow( tZ - tSpherePosZ, tSphereExpon ), 1.0 / tSphereExpon );

        return std::abs( tReturnValue ) < tMinLevs ? tMinLevs : tReturnValue;
    }

    void
    Func_Sphere_Deriv(
            const moris::Matrix< moris::DDRMat >&           aCoordinates,
            const Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::Matrix< DDRMat >&                        aFieldSensitivity )
    {
        // derivative of level set function wrt sphere radius
        aFieldSensitivity = { { 1.0 } };
    }

    /* ------------------------------------------------------------------------ */
    // geometry parameters & LS functions

    moris::real
    Func_Plane(
            const moris::Matrix< DDRMat >& aCoordinates,
            const Vector< real >&          aGeometryParameters )
    {
        // get coordinates
        real tX = aCoordinates( 0 );
        real tY = aCoordinates( 1 );
        real tZ = aCoordinates( 2 );

        // get normal
        real tNx = aGeometryParameters( 0 );
        real tNy = aGeometryParameters( 1 );
        real tNz = aGeometryParameters( 2 );

        // get point on plane
        real tPx = aGeometryParameters( 3 );
        real tPy = aGeometryParameters( 4 );
        real tPz = aGeometryParameters( 5 );

        real tReturnValue = tNx * ( tPx - tX ) + tNy * ( tPy - tY ) + tNz * ( tPz - tZ );

        return std::abs( tReturnValue ) < tMinLevs ? tMinLevs : tReturnValue;
    }

    /* ------------------------------------------------------------------------ */
    // Phase assignement

    uint
    get_phase_index( const Bitset< 7 >& aGeometrySigns )
    {
        // by default in fluid
        uint tPhaseIndex = 0;

        // solid
        if ( aGeometrySigns.test( 0 ) )
        {
            return 1;
        }

        // Phase void inflow
        if ( aGeometrySigns.test( 1 ) )
        {
            return 2;
        }

        // Phase void outflow
        if ( aGeometrySigns.test( 2 ) )
        {
            return 3;
        }

        // Phase void lateral
        if ( aGeometrySigns.test( 3 )
                || aGeometrySigns.test( 4 )
                || aGeometrySigns.test( 5 )
                || aGeometrySigns.test( 6 ) )
        {
            return 4;
        }

        // Phase fluid
        return tPhaseIndex;
    }

    /* ------------------------------------------------------------------------ */

    Matrix< DDSMat >
    get_constraint_types()
    {
        Matrix< DDSMat > tConstraintTypes( tNumConstraints, 1, 1 );

        return tConstraintTypes;
    }

    /* ------------------------------------------------------------------------ */
    /*
    0 - IQIInletThermalEnergy
    1 - IQIOutletThermalEnergy
    2 - IQIInletTotalPressure
    3 - IQIOutletTotalPressure
    4 - IQIPerimeterIfc
    5 - IQIInletMassFlow
    6 - IQIOutletMassFlow
    7 - IQIMaxTemp
    8 - IQISolidVolume
     */

    Matrix< DDRMat >
    compute_objectives(
            const Vector< real >& aADVs,
            const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tObjectives = { { aCriteria( 0 ) } };

        return tObjectives;
    }

    /* ------------------------------------------------------------------------ */

    Matrix< DDRMat >
    compute_dobjective_dadv(
            const Vector< real >& aADVs,
            const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tDObjectiveDADV( 1, aADVs.size(), 0.0 );

        return tDObjectiveDADV;
    }

    /* ------------------------------------------------------------------------ */

    Matrix< DDRMat >
    compute_dobjective_dcriteria(
            const Vector< real >& aADVs,
            const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tDObjectiveDCriteria( 1, aCriteria.size(), 0.0 );

        tDObjectiveDCriteria( 0 ) = 1.0;

        return tDObjectiveDCriteria;
    }

    /* ------------------------------------------------------------------------ */

    Matrix< DDRMat >
    compute_constraints(
            const Vector< real >& aADVs,
            const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tConstraints( tNumConstraints, 1 );

        // Mass flux constraint
        tConstraints( 0 ) = aCriteria( 0 );
        tConstraints( 1 ) = aCriteria( 1 );
        tConstraints( 2 ) = aCriteria( 2 );
        tConstraints( 3 ) = aCriteria( 3 );
        tConstraints( 4 ) = aCriteria( 4 );
        tConstraints( 5 ) = aCriteria( 5 );
        tConstraints( 6 ) = aCriteria( 6 );
        tConstraints( 7 ) = aCriteria( 7 );
        tConstraints( 8 ) = aCriteria( 8 );

        return tConstraints;
    }

    /* ------------------------------------------------------------------------ */

    Matrix< DDRMat >
    compute_dconstraint_dadv(
            const Vector< real >& aADVs,
            const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tDConstraintDADV( tNumConstraints, aADVs.size(), 0.0 );

        return tDConstraintDADV;
    }

    /* ------------------------------------------------------------------------ */

    Matrix< DDRMat >
    compute_dconstraint_dcriteria(
            const Vector< real >& aADVs,
            const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tDConstraintDCriteria( tNumConstraints, aCriteria.size(), 0.0 );

        tDConstraintDCriteria( 0, 0 ) = 1.0;
        tDConstraintDCriteria( 1, 1 ) = 1.0;
        tDConstraintDCriteria( 2, 2 ) = 1.0;
        tDConstraintDCriteria( 3, 3 ) = 1.0;
        tDConstraintDCriteria( 4, 4 ) = 1.0;
        tDConstraintDCriteria( 5, 5 ) = 1.0;
        tDConstraintDCriteria( 6, 6 ) = 1.0;
        tDConstraintDCriteria( 7, 7 ) = 1.0;
        tDConstraintDCriteria( 8, 8 ) = 1.0;

        return tDConstraintDCriteria;
    }

    /* ------------------------------------------------------------------------ */

    void
    OPTParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists( 0 ).push_back( prm::create_opt_problem_parameter_list() );
        aParameterLists( 0 ).set( "is_optimization_problem", true );
        aParameterLists( 0 ).set( "problem", "user_defined" );
        aParameterLists( 0 ).set( "library", tSoFile );

        aParameterLists( 2 ).push_back( moris::prm::create_sweep_parameter_list() );
        aParameterLists( 2 ).set( "hdf5_path", tHdf5File );
        aParameterLists( 2 ).set( "num_evaluations_per_adv", "1" );
        aParameterLists( 2 ).set( "finite_difference_type", "all" );
        aParameterLists( 2 ).set( "finite_difference_epsilons", tSweepFdEpsilon );
    }

    /* ------------------------------------------------------------------------ */

    void
    HMRParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists( 0 ).push_back( prm::create_hmr_parameter_list() );

        aParameterLists( 0 ).set( "number_of_elements_per_dimension", tNumElemsPerDim );
        aParameterLists( 0 ).set( "domain_dimensions", tDomainDims );
        aParameterLists( 0 ).set( "domain_offset", tDomainOffset );
        aParameterLists( 0 ).set( "domain_sidesets", tDomainSidesets );

        aParameterLists( 0 ).set( "lagrange_output_meshes", "0" );

        aParameterLists( 0 ).set( "lagrange_orders", tLagrangeOrder );
        aParameterLists( 0 ).set( "lagrange_pattern", tLagrangePattern );

        aParameterLists( 0 ).set( "bspline_orders", tBsplineOrder );
        aParameterLists( 0 ).set( "bspline_pattern", "0,1" );

        aParameterLists( 0 ).set( "lagrange_to_bspline", "0,1" );

        aParameterLists( 0 ).set( "truncate_bsplines", 1 );
        aParameterLists( 0 ).set( "refinement_buffer", tRefineBuffer );
        aParameterLists( 0 ).set( "staircase_buffer", tRefineBuffer );

        aParameterLists( 0 ).set( "initial_refinement", tInitialRef );
        aParameterLists( 0 ).set( "initial_refinement_pattern", "0,1" );

        aParameterLists( 0 ).set( "use_number_aura", 1 );

        aParameterLists( 0 ).set( "use_multigrid", 0 );
        aParameterLists( 0 ).set( "severity_level", 0 );
    }

    /* ------------------------------------------------------------------------ */

    void
    XTKParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).push_back( prm::create_xtk_parameter_list() );
        aParameterLists( 0 ).set( "decompose", true );
        aParameterLists( 0 ).set( "decomposition_type", "conformal" );
        aParameterLists( 0 ).set( "enrich", true );
        aParameterLists( 0 ).set( "basis_rank", "bspline" );
        aParameterLists( 0 ).set( "enrich_mesh_indices", "0,1" );
        aParameterLists( 0 ).set( "ghost_stab", tUseGhost );
        aParameterLists( 0 ).set( "multigrid", false );
        aParameterLists( 0 ).set( "verbose", false );
        aParameterLists( 0 ).set( "high_to_low_dbl_side_sets", true );
        aParameterLists( 0 ).set( "print_enriched_ig_mesh", false );
        aParameterLists( 0 ).set( "exodus_output_XTK_ig_mesh", true );
    }

    /* ------------------------------------------------------------------------ */

    void
    GENParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists( 0 ).push_back( prm::create_gen_parameter_list() );
        aParameterLists( 0 ).set( "IQI_types",
                "IQIInletThermalEnergy",
                "IQIOutletThermalEnergy",
                "IQIInletTotalPressure",
                "IQIOutletTotalPressure"
                "IQIPerimeterItf",
                "IQIInletMassFlow",
                "IQIOutletMassFlow",
                "IQIMaxTemp",
                "IQISolidVolume" );

        aParameterLists( 0 ).set( "number_of_phases", 5 );
        aParameterLists( 0 ).set( "phase_function_name", "get_phase_index" );
        aParameterLists( 0 ).set( "output_mesh_file", tGENOutputFile );
        aParameterLists( 0 ).set( "time_offset", 10.0 );

        // Inclusions
        aParameterLists( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists( 1 ).set( "field_function_name", "Func_Sphere" );
        aParameterLists( 1 ).set( "sensitivity_function_name", "Func_Sphere_Deriv" );
        aParameterLists( 1 ).insert( "radius", Design_Variable( tSphereRadius * 0.9, tSphereRadius, tSphereRadius * 1.1 ) );
        aParameterLists( 1 ).set( "number_of_refinements", tInterfaceRefinementSphere );
        aParameterLists( 1 ).set( "refinement_mesh_index", 0 );
        aParameterLists( 1 ).set( "isocontour_threshold", 0.0 );
        aParameterLists( 1 ).set( "isocontour_tolerance", 1.0e-12 );
        aParameterLists( 1 ).set( "intersection_tolerance", 1.0e-12 );

        // Inlet plane
        aParameterLists( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists( 1 ).insert( "variable_1", 1.0 );
        aParameterLists( 1 ).insert( "variable_2", 0.0 );
        aParameterLists( 1 ).insert( "variable_3", 0.0 );
        aParameterLists( 1 ).insert( "variable_4", tGeoShift );
        aParameterLists( 1 ).insert( "variable_5", 0.0 );
        aParameterLists( 1 ).insert( "variable_6", 0.0 );
        aParameterLists( 1 ).set( "field_function_name", "Func_Plane" );
        aParameterLists( 1 ).set( "number_of_refinements", tInterfaceRefinementWalls );
        aParameterLists( 1 ).set( "refinement_mesh_index", 0 );
        aParameterLists( 1 ).set( "isocontour_threshold", 0.0 );
        aParameterLists( 1 ).set( "isocontour_tolerance", 1.0e-12 );
        aParameterLists( 1 ).set( "intersection_tolerance", 1.0e-12 );

        // Out plane
        aParameterLists( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists( 1 ).set( "field_function_name", "Func_Plane" );
        aParameterLists( 1 ).insert( "variable_1", -1.0 );
        aParameterLists( 1 ).insert( "variable_2", 0.0 );
        aParameterLists( 1 ).insert( "variable_3", 0.0 );
        aParameterLists( 1 ).insert( "variable_4", tChannelLength + tGeoShift );
        aParameterLists( 1 ).insert( "variable_5", 0.0 );
        aParameterLists( 1 ).insert( "variable_6", 0.0 );
        aParameterLists( 1 ).set( "number_of_refinements", tInterfaceRefinementWalls );
        aParameterLists( 1 ).set( "refinement_mesh_index", 0 );
        aParameterLists( 1 ).set( "isocontour_threshold", 0.0 );
        aParameterLists( 1 ).set( "isocontour_tolerance", 1.0e-12 );
        aParameterLists( 1 ).set( "intersection_tolerance", 1.0e-12 );

        // Lower plane
        aParameterLists( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists( 1 ).set( "field_function_name", "Func_Plane" );
        aParameterLists( 1 ).insert( "variable_1", 0.0 );
        aParameterLists( 1 ).insert( "variable_2", 1.0 );
        aParameterLists( 1 ).insert( "variable_3", 0.0 );
        aParameterLists( 1 ).insert( "variable_4", 0.0 );
        aParameterLists( 1 ).insert( "variable_5", tGeoShift );
        aParameterLists( 1 ).insert( "variable_6", 0.0 );
        aParameterLists( 1 ).set( "number_of_refinements", tInterfaceRefinementWalls );
        aParameterLists( 1 ).set( "refinement_mesh_index", 0 );
        aParameterLists( 1 ).set( "isocontour_threshold", 0.0 );
        aParameterLists( 1 ).set( "isocontour_tolerance", 1.0e-12 );
        aParameterLists( 1 ).set( "intersection_tolerance", 1.0e-12 );

        // Upper plane
        aParameterLists( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists( 1 ).set( "field_function_name", "Func_Plane" );
        aParameterLists( 1 ).insert( "variable_1", 0.0 );
        aParameterLists( 1 ).insert( "variable_2", -1.0 );
        aParameterLists( 1 ).insert( "variable_3", 0.0 );
        aParameterLists( 1 ).insert( "variable_4", 0.0 );
        aParameterLists( 1 ).insert( "variable_5", tChannelHeight + tGeoShift );
        aParameterLists( 1 ).insert( "variable_6", 0.0 );
        aParameterLists( 1 ).set( "number_of_refinements", tInterfaceRefinementWalls );
        aParameterLists( 1 ).set( "refinement_mesh_index", 0 );
        aParameterLists( 1 ).set( "isocontour_threshold", 0.0 );
        aParameterLists( 1 ).set( "isocontour_tolerance", 1.0e-12 );
        aParameterLists( 1 ).set( "intersection_tolerance", 1.0e-12 );

        // Back Plan
        aParameterLists( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists( 1 ).set( "field_function_name", "Func_Plane" );
        aParameterLists( 1 ).insert( "variable_1", 0.0 );
        aParameterLists( 1 ).insert( "variable_2", 0.0 );
        aParameterLists( 1 ).insert( "variable_3", 1.0 );
        aParameterLists( 1 ).insert( "variable_4", 0.0 );
        aParameterLists( 1 ).insert( "variable_5", 0.0 );
        aParameterLists( 1 ).insert( "variable_6", tGeoShift );
        aParameterLists( 1 ).set( "number_of_refinements", tInterfaceRefinementWalls );
        aParameterLists( 1 ).set( "refinement_mesh_index", 0 );
        aParameterLists( 1 ).set( "isocontour_threshold", 0.0 );
        aParameterLists( 1 ).set( "isocontour_tolerance", 1.0e-12 );
        aParameterLists( 1 ).set( "intersection_tolerance", 1.0e-12 );

        // Front plane
        aParameterLists( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists( 1 ).set( "field_function_name", "Func_Plane" );
        aParameterLists( 1 ).insert( "variable_1", 0.0 );
        aParameterLists( 1 ).insert( "variable_2", 0.0 );
        aParameterLists( 1 ).insert( "variable_3", -1.0 );
        aParameterLists( 1 ).insert( "variable_4", 0.0 );
        aParameterLists( 1 ).insert( "variable_5", 0.0 );
        aParameterLists( 1 ).insert( "variable_6", tChannelWidth + tGeoShift );
        aParameterLists( 1 ).set( "number_of_refinements", tInterfaceRefinementWalls );
        aParameterLists( 1 ).set( "refinement_mesh_index", 0 );
        aParameterLists( 1 ).set( "isocontour_threshold", 0.0 );
        aParameterLists( 1 ).set( "isocontour_tolerance", 1.0e-12 );
        aParameterLists( 1 ).set( "intersection_tolerance", 1.0e-12 );
    }

    /* ------------------------------------------------------------------------ */

    void
    FEMParameterList( Module_Parameter_Lists& aParameterLists )
    {
        if ( par_rank() == 0 )
        {
            std::cout << "presref        = " << presref << std::endl;
            std::cout << "lenref         = " << lenref << std::endl;
            std::cout << "rhoref         = " << rhoref << std::endl;
            std::cout << "tempref        = " << tempref << std::endl;

            std::cout << std::endl;

            std::cout << "velref         = " << velref << std::endl;
            std::cout << "timeref        = " << timeref << std::endl;
            std::cout << "massref        = " << massref << std::endl;

            std::cout << std::endl;

            std::cout << "tLengthScale   = " << tLengthScale << std::endl;
            std::cout << "tTimeScale     = " << tTimeScale << std::endl;
            std::cout << "tMassScale     = " << tMassScale << std::endl;
            std::cout << "tTempScale     = " << tTempScale << std::endl;

            std::cout << std::endl;

            std::cout << "tPressureScale = " << tPressureScale << std::endl;
            std::cout << "tEnergyScale   = " << tEnergyScale << std::endl;
            std::cout << "tPowerScale    = " << tPowerScale << std::endl;
            std::cout << "tDensityScale  = " << tDensityScale << std::endl;

            std::cout << std::endl;

            std::cout << "tInletPressure       = " << tInletPressure << std::endl;
            std::cout << "tInletTemperature    = " << tInletTemperature << std::endl;
            std::cout << "tVolumetricHeatLoad  = " << tVolumetricHeatLoad << std::endl;

            std::cout << std::endl;

            std::cout << "tFluidViscosity      = " << tFluidViscosity << std::endl;
            std::cout << "tFluidDensity        = " << tFluidDensity << std::endl;
            std::cout << "tFluidCapacity       = " << tFluidCapacity << std::endl;
            std::cout << "tFluidConductivity   = " << tFluidConductivity << std::endl;
            std::cout << "tFluidPressureSpring = " << tFluidPressureSpring << std::endl;

            std::cout << std::endl;

            std::cout << "tSolidDensity        = " << tSolidDensity << std::endl;
            std::cout << "tSolidCapacity       = " << tSolidCapacity << std::endl;
            std::cout << "tSolidConductivity   = " << tSolidConductivity << std::endl;

            std::cout << std::endl;

            std::cout << "Reynolds number     = " << 1.0 / std::stod( tFluidViscosity ) << " (" << reynolds << ")\n";
        }

        // create a cell of cell of parameter list for fem
        uint tPropIndex  = 0;
        uint tCMIndex    = 1;
        uint tSPIndex    = 2;
        uint tIWGIndex   = 3;
        uint tIQIIndex   = 4;
        uint tFEMIndex   = 5;
        uint tPhaseIndex = 7;

        //------------------------------------------------------------------------------

        aParameterLists( tPhaseIndex ).push_back( prm::create_phase_parameter_list() );
        aParameterLists( tPhaseIndex ).set( "phase_name", "PhaseFluid" );
        aParameterLists( tPhaseIndex ).set( "phase_indices", "0" );

        aParameterLists( tPhaseIndex ).push_back( prm::create_phase_parameter_list() );
        aParameterLists( tPhaseIndex ).set( "phase_name", "PhaseSolid" );
        aParameterLists( tPhaseIndex ).set( "phase_indices", "1" );

        aParameterLists( tPhaseIndex ).push_back( prm::create_phase_parameter_list() );
        aParameterLists( tPhaseIndex ).set( "phase_name", "PhaseVoidFront" );
        aParameterLists( tPhaseIndex ).set( "phase_indices", "2" );

        aParameterLists( tPhaseIndex ).push_back( prm::create_phase_parameter_list() );
        aParameterLists( tPhaseIndex ).set( "phase_name", "PhaseVoidBack" );
        aParameterLists( tPhaseIndex ).set( "phase_indices", "3" );

        aParameterLists( tPhaseIndex ).push_back( prm::create_phase_parameter_list() );
        aParameterLists( tPhaseIndex ).set( "phase_name", "PhaseVoidLateral" );
        aParameterLists( tPhaseIndex ).set( "phase_indices", "4" );

        aParameterLists( tPhaseIndex ).push_back( prm::create_phase_parameter_list() );
        aParameterLists( tPhaseIndex ).set( "phase_name", "PhaseAll" );
        aParameterLists( tPhaseIndex ).set( "phase_indices", "0,1" );

        //------------------------------------------------------------------------------
        // fill the property part of the parameter list

        // fluid properties ------------------------------------------------------------
        // create fluid viscosity property
        aParameterLists( tPropIndex ).push_back( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropFluidViscosity" );
        aParameterLists( tPropIndex ).set( "function_parameters", tFluidViscosity );

        // create fluid density property
        aParameterLists( tPropIndex ).push_back( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropFluidDensity" );
        aParameterLists( tPropIndex ).set( "function_parameters", tFluidDensity );

        // create fluid capacity property
        aParameterLists( tPropIndex ).push_back( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropFluidCapacity" );
        aParameterLists( tPropIndex ).set( "function_parameters", tFluidCapacity );

        // create fluid conductivity property
        aParameterLists( tPropIndex ).push_back( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropFluidConductivity" );
        aParameterLists( tPropIndex ).set( "function_parameters", tFluidConductivity );

        // create fluid pressure spring property
        aParameterLists( tPropIndex ).push_back( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropFluidPressureSpring" );
        aParameterLists( tPropIndex ).set( "function_parameters", tFluidPressureSpring );

        // solid properties ----------------------------------------------------------
        // create solid B density property
        aParameterLists( tPropIndex ).push_back( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropSolidDensity" );
        aParameterLists( tPropIndex ).set( "function_parameters", tSolidDensity );

        // create solid B capacity property
        aParameterLists( tPropIndex ).push_back( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropSolidCapacity" );
        aParameterLists( tPropIndex ).set( "function_parameters", tSolidCapacity );

        // create solid B conductivity property
        aParameterLists( tPropIndex ).push_back( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropSolidConductivity" );
        aParameterLists( tPropIndex ).set( "function_parameters", tSolidConductivity );

        // BC properties ---------------------------------------------------------------
        // create inlet pressure property
        aParameterLists( tPropIndex ).push_back( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropInletPressure" );
        aParameterLists( tPropIndex ).set( "function_parameters", tInletPressure );

        // create inlet temperature property
        aParameterLists( tPropIndex ).push_back( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropInletTemp" );
        aParameterLists( tPropIndex ).set( "function_parameters", tInletTemperature );

        // create wall velocity property
        aParameterLists( tPropIndex ).push_back( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropZeroU" );
        aParameterLists( tPropIndex ).set( "function_parameters", "0.0;0.0;0.0" );

        // create symmetry velocity property (x-z plane)
        aParameterLists( tPropIndex ).push_back( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropSelectUY" );
        aParameterLists( tPropIndex ).set( "function_parameters", "0.0,0.0,0.0;0.0,1.0,0.0;0.0,0.0,0.0" );

        // create symmetry velocity property (x-y plane)
        aParameterLists( tPropIndex ).push_back( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropSelectUZ" );
        aParameterLists( tPropIndex ).set( "function_parameters", "0.0,0.0,0.0;0.0,0.0,0.0;0.0,0.0,1.0" );

        // create heat load property
        aParameterLists( tPropIndex ).push_back( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropVolumetricHeatFlux" );
        aParameterLists( tPropIndex ).set( "function_parameters", tVolumetricHeatLoad );

        //------------------------------------------------------------------------------
        // fill the constitutive model part of the parameter list

        // fluid CM --------------------------------------------------------------------
        // create fluid CM
        aParameterLists( tCMIndex ).push_back( prm::create_constitutive_model_parameter_list() );
        aParameterLists( tCMIndex ).set( "constitutive_name", "CMFluid" );
        aParameterLists( tCMIndex ).set( "phase_name", "PhaseFluid" );
        aParameterLists( tCMIndex ).set( "constitutive_type", fem::Constitutive_Type::FLUID_INCOMPRESSIBLE );
        aParameterLists( tCMIndex ).set( "dof_dependencies", std::pair< std::string, std::string >( "VX,VY,VZ;P", "Velocity,Pressure" ) );
        aParameterLists( tCMIndex ).set( "properties",
                "PropFluidViscosity,Viscosity;"
                "PropFluidDensity  ,Density" );

        // create fluid diffusion CM
        aParameterLists( tCMIndex ).push_back( prm::create_constitutive_model_parameter_list() );
        aParameterLists( tCMIndex ).set( "constitutive_name", "CMFluidDiffusion" );
        aParameterLists( tCMIndex ).set( "phase_name", "PhaseFluid" );
        aParameterLists( tCMIndex ).set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
        aParameterLists( tCMIndex ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists( tCMIndex ).set( "properties",
                "PropFluidConductivity,Conductivity;"
                "PropFluidDensity     ,Density;"
                "PropFluidCapacity    ,HeatCapacity" );

        // solid CM --------------------------------------------------------------------
        // create solid B diffusion CM
        aParameterLists( tCMIndex ).push_back( prm::create_constitutive_model_parameter_list() );
        aParameterLists( tCMIndex ).set( "constitutive_name", "CMSolidDiffusion" );
        aParameterLists( tCMIndex ).set( "phase_name", "PhaseSolid" );
        aParameterLists( tCMIndex ).set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
        aParameterLists( tCMIndex ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists( tCMIndex ).set( "properties",
                "PropSolidConductivity,Conductivity;"
                "PropSolidDensity     ,Density;"
                "PropSolidCapacity    ,HeatCapacity" );

        //------------------------------------------------------------------------------
        // fill the stabilization parameter part of the parameter list

        // create SUPG/PSG fluid
        aParameterLists( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( tSPIndex ).set( "stabilization_name", "SPIncFlow" );
        aParameterLists( tSPIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tSPIndex ).set( "stabilization_type", fem::Stabilization_Type::INCOMPRESSIBLE_FLOW );
        aParameterLists( tSPIndex ).set( "function_parameters", tCInv );
        aParameterLists( tSPIndex ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY,VZ;P", "Velocity,Pressure" ) );
        aParameterLists( tSPIndex ).set( "leader_properties",
                "PropFluidViscosity,Viscosity;"
                "PropFluidDensity  ,Density" );

        // create SUPG fluid temperature advection
        aParameterLists( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( tSPIndex ).set( "stabilization_name", "SPSUPGTemp" );
        aParameterLists( tSPIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tSPIndex ).set( "stabilization_type", fem::Stabilization_Type::SUPG_ADVECTION );
        aParameterLists( tSPIndex ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY,VZ", "Velocity" ) );
        aParameterLists( tSPIndex ).set( "leader_properties", "PropFluidConductivity,Conductivity" );

        // create Nitsche for fluid velocity
        aParameterLists( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( tSPIndex ).set( "stabilization_name", "SPNitscheU" );
        aParameterLists( tSPIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tSPIndex ).set( "stabilization_type", fem::Stabilization_Type::VELOCITY_DIRICHLET_NITSCHE );
        aParameterLists( tSPIndex ).set( "function_parameters", "100.0/1.0" );
        aParameterLists( tSPIndex ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY,VZ", "Velocity" ) );
        aParameterLists( tSPIndex ).set( "leader_properties",
                "PropFluidViscosity,Viscosity;"
                "PropFluidDensity,Density" );

        // create Nitsche for fluid temperature
        aParameterLists( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( tSPIndex ).set( "stabilization_name", "SPNitscheT" );
        aParameterLists( tSPIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tSPIndex ).set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        aParameterLists( tSPIndex ).set( "function_parameters", "100.0" );
        aParameterLists( tSPIndex ).set( "leader_properties", "PropFluidConductivity,Material" );

        // create Nitsche for fluid/solid diffusion interface
        aParameterLists( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( tSPIndex ).set( "stabilization_name", "SPFSInterfaceNitsche" );
        aParameterLists( tSPIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tSPIndex ).set( "follower_phase_name", "PhaseSolid" );
        aParameterLists( tSPIndex ).set( "stabilization_type", fem::Stabilization_Type::NITSCHE_INTERFACE );
        aParameterLists( tSPIndex ).set( "function_parameters", "100.0" );
        aParameterLists( tSPIndex ).set( "leader_properties", "PropFluidConductivity,Material" );
        aParameterLists( tSPIndex ).set( "follower_properties", "PropSolidConductivity,Material" );

        // create ghost penalty viscous
        aParameterLists( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( tSPIndex ).set( "stabilization_name", "SPGPViscous" );
        aParameterLists( tSPIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tSPIndex ).set( "stabilization_type", fem::Stabilization_Type::VISCOUS_GHOST );
        aParameterLists( tSPIndex ).set( "function_parameters", "0.05" );
        aParameterLists( tSPIndex ).set( "leader_properties", "PropFluidViscosity,Viscosity" );

        // create ghost penalty convective
        aParameterLists( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( tSPIndex ).set( "stabilization_name", "SPGPVelocity" );
        aParameterLists( tSPIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tSPIndex ).set( "stabilization_type", fem::Stabilization_Type::CONVECTIVE_GHOST );
        aParameterLists( tSPIndex ).set( "function_parameters", "0.05" );
        aParameterLists( tSPIndex ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY,VZ", "Velocity" ) );
        aParameterLists( tSPIndex ).set( "leader_properties", "PropFluidDensity,Density" );

        // create ghost penalty pressure
        aParameterLists( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( tSPIndex ).set( "stabilization_name", "SPGPPressure" );
        aParameterLists( tSPIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tSPIndex ).set( "stabilization_type", fem::Stabilization_Type::PRESSURE_GHOST );
        aParameterLists( tSPIndex ).set( "function_parameters", "0.005/1.0" );
        aParameterLists( tSPIndex ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY,VZ", "Velocity" ) );
        aParameterLists( tSPIndex ).set( "leader_properties",
                "PropFluidViscosity,Viscosity;"
                "PropFluidDensity,Density" );

        // create ghost penalty fluid temperature
        aParameterLists( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( tSPIndex ).set( "stabilization_name", "SPGPFluidTemp" );
        aParameterLists( tSPIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tSPIndex ).set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        aParameterLists( tSPIndex ).set( "function_parameters", "0.05" );
        aParameterLists( tSPIndex ).set( "leader_properties", "PropFluidConductivity,Material" );

        // create ghost penalty solid temperature
        aParameterLists( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( tSPIndex ).set( "stabilization_name", "SPGPSolidTemp" );
        aParameterLists( tSPIndex ).set( "leader_phase_name", "PhaseSolid" );
        aParameterLists( tSPIndex ).set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        aParameterLists( tSPIndex ).set( "function_parameters", "0.05" );
        aParameterLists( tSPIndex ).set( "leader_properties", "PropSolidConductivity,Material" );

        //------------------------------------------------------------------------------
        // fill the IWG part of the parameter list

        // fluid bulk IWGs -------------------------------------------------------------
        // NS incompressible (velocity)
        aParameterLists( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGVelocityBulk" );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_BULK );
        aParameterLists( tIWGIndex ).set( "dof_residual", "VX,VY,VZ" );
        aParameterLists( tIWGIndex ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists( tIWGIndex ).set( "stabilization_parameters", "SPIncFlow,IncompressibleFlow" );

        // NS incompressible (pressure)
        aParameterLists( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGPressureBulk" );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_BULK );
        aParameterLists( tIWGIndex ).set( "dof_residual", "P" );
        aParameterLists( tIWGIndex ).set( "leader_properties", "PropFluidPressureSpring,PressureSpring" );
        aParameterLists( tIWGIndex ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists( tIWGIndex ).set( "stabilization_parameters", "SPIncFlow,IncompressibleFlow" );

        // diffusion
        aParameterLists( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGFluidDiffusionBulk" );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        aParameterLists( tIWGIndex ).set( "dof_residual", "TEMP" );
        aParameterLists( tIWGIndex ).set( "leader_constitutive_models", "CMFluidDiffusion,Diffusion" );

        // advection
        aParameterLists( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGFluidAdvectionBulk" );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::ADVECTION_BULK );
        aParameterLists( tIWGIndex ).set( "dof_residual", "TEMP" );
        aParameterLists( tIWGIndex ).set( "leader_constitutive_models", "CMFluidDiffusion,Diffusion" );
        aParameterLists( tIWGIndex ).set( "stabilization_parameters", "SPSUPGTemp,SUPG" );

        // solid bulk IWGs -----------------------------------------------------------
        // diffusion
        aParameterLists( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGSolidDiffusionBulk" );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseSolid" );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        aParameterLists( tIWGIndex ).set( "dof_residual", "TEMP" );
        aParameterLists( tIWGIndex ).set( "leader_properties", "PropVolumetricHeatFlux,Load" );
        aParameterLists( tIWGIndex ).set( "leader_constitutive_models", "CMSolidDiffusion,Diffusion" );

        // fluid/solid ifc IWGs ------------------------------------------------------

        // temperature
        aParameterLists( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGInterfaceFluidSolid" );
        aParameterLists( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIWGIndex ).set( "follower_phase_name", "PhaseSolid" );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_INTERFACE_UNSYMMETRIC_NITSCHE );
        aParameterLists( tIWGIndex ).set( "dof_residual", "TEMP" );
        aParameterLists( tIWGIndex ).set( "leader_constitutive_models", "CMFluidDiffusion,Diffusion" );
        aParameterLists( tIWGIndex ).set( "follower_constitutive_models", "CMSolidDiffusion,Diffusion" );
        aParameterLists( tIWGIndex ).set( "stabilization_parameters", "SPFSInterfaceNitsche,NitscheInterface" );

        // zero velocity
        aParameterLists( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGZeroVelocity" );
        aParameterLists( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIWGIndex ).set( "neighbor_phases", "PhaseSolid" );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists( tIWGIndex ).set( "dof_residual", "VX,VY,VZ" );
        aParameterLists( tIWGIndex ).set( "leader_properties", "PropZeroU,Dirichlet" );
        aParameterLists( tIWGIndex ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists( tIWGIndex ).set( "stabilization_parameters", "SPNitscheU,DirichletNitsche" );

        // zero velocity (pressure part)
        aParameterLists( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGZeroPressure" );
        aParameterLists( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIWGIndex ).set( "neighbor_phases", "PhaseSolid" );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists( tIWGIndex ).set( "dof_residual", "P" );
        aParameterLists( tIWGIndex ).set( "leader_properties", "PropZeroU,Dirichlet" );
        aParameterLists( tIWGIndex ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );

        // Laterial BC IWG (adiabatic) ----------------------------------------------------------------

        // zero velocity along lateral sides
        aParameterLists( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGZeroVelocityVoid" );
        aParameterLists( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIWGIndex ).set( "neighbor_phases", "PhaseVoidLateral" );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists( tIWGIndex ).set( "dof_residual", "VX,VY,VZ" );
        aParameterLists( tIWGIndex ).set( "leader_properties", "PropZeroU,Dirichlet" );
        aParameterLists( tIWGIndex ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists( tIWGIndex ).set( "stabilization_parameters", "SPNitscheU,DirichletNitsche" );

        // zero velocity along lateral sides (pressure part)
        aParameterLists( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGZeroPressureVoid" );
        aParameterLists( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIWGIndex ).set( "neighbor_phases", "PhaseVoidLateral" );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists( tIWGIndex ).set( "dof_residual", "P" );
        aParameterLists( tIWGIndex ).set( "leader_properties", "PropZeroU,Dirichlet" );
        aParameterLists( tIWGIndex ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );

        // Inlet BC IWG ----------------------------------------------------------------

        // inlet pressure
        aParameterLists( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGInletImposedPressure" );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIWGIndex ).set( "neighbor_phases", "PhaseVoidFront" );
        aParameterLists( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_IMPOSED_PRESSURE );
        aParameterLists( tIWGIndex ).set( "dof_residual", "VX,VY,VZ" );
        aParameterLists( tIWGIndex ).set( "leader_properties", "PropInletPressure,Pressure" );

        // inlet temperature
        aParameterLists( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGInletTemp" );
        aParameterLists( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIWGIndex ).set( "neighbor_phases", "PhaseVoidFront" );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists( tIWGIndex ).set( "dof_residual", "TEMP" );
        aParameterLists( tIWGIndex ).set( "leader_properties", "PropInletTemp,Dirichlet" );
        aParameterLists( tIWGIndex ).set( "leader_constitutive_models", "CMFluidDiffusion,Diffusion" );
        aParameterLists( tIWGIndex ).set( "stabilization_parameters", "SPNitscheT,DirichletNitsche" );

        // Ghost  ----------------------------------------------------------------

        if ( tUseGhost )
        {
            // ghost viscous
            aParameterLists( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
            aParameterLists( tIWGIndex ).set( "IWG_name", "IWGGPViscous" );
            aParameterLists( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
            aParameterLists( tIWGIndex ).set( "follower_phase_name", "PhaseFluid" );
            aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists( tIWGIndex ).set( "dof_residual", "VX,VY,VZ" );
            aParameterLists( tIWGIndex ).set( "stabilization_parameters", "SPGPViscous,GhostSP" );
            aParameterLists( tIWGIndex ).set( "ghost_order", (uint)tDispOrder );

            // ghost convective
            aParameterLists( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
            aParameterLists( tIWGIndex ).set( "IWG_name", "IWGGPConvective" );
            aParameterLists( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
            aParameterLists( tIWGIndex ).set( "follower_phase_name", "PhaseFluid" );
            aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists( tIWGIndex ).set( "dof_residual", "VX,VY,VZ" );
            aParameterLists( tIWGIndex ).set( "stabilization_parameters", "SPGPVelocity,GhostSP" );
            aParameterLists( tIWGIndex ).set( "ghost_order", (uint)tDispOrder );

            // ghost pressure
            aParameterLists( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
            aParameterLists( tIWGIndex ).set( "IWG_name", "IWGGPPressure" );
            aParameterLists( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
            aParameterLists( tIWGIndex ).set( "follower_phase_name", "PhaseFluid" );
            aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists( tIWGIndex ).set( "dof_residual", "P" );
            aParameterLists( tIWGIndex ).set( "stabilization_parameters", "SPGPPressure,GhostSP" );
            aParameterLists( tIWGIndex ).set( "ghost_order", (uint)tDispOrder );

            // ghost fluid temperature
            aParameterLists( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
            aParameterLists( tIWGIndex ).set( "IWG_name", "IWGGPFluidTemp" );
            aParameterLists( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
            aParameterLists( tIWGIndex ).set( "follower_phase_name", "PhaseFluid" );
            aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists( tIWGIndex ).set( "dof_residual", "TEMP" );
            aParameterLists( tIWGIndex ).set( "stabilization_parameters", "SPGPFluidTemp,GhostSP" );
            aParameterLists( tIWGIndex ).set( "ghost_order", (uint)tDispOrder );

            // ghost solid A temperature
            aParameterLists( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
            aParameterLists( tIWGIndex ).set( "IWG_name", "IWGGPSolidTemp" );
            aParameterLists( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseSolid" );
            aParameterLists( tIWGIndex ).set( "follower_phase_name", "PhaseSolid" );
            aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists( tIWGIndex ).set( "dof_residual", "TEMP" );
            aParameterLists( tIWGIndex ).set( "stabilization_parameters", "SPGPSolidTemp,GhostSP" );
            aParameterLists( tIWGIndex ).set( "ghost_order", (uint)tDispOrder );
        }

        //------------------------------------------------------------------------------
        // fill the IQI part of the parameter list

        // velocity VX
        aParameterLists( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIBulkVX" );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists( tIQIIndex ).set( "dof_quantity", "VX,VY,VZ" );
        aParameterLists( tIQIIndex ).set( "vectorial_field_index", 0 );

        // velocity VY
        aParameterLists( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIBulkVY" );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists( tIQIIndex ).set( "dof_quantity", "VX,VY,VZ" );
        aParameterLists( tIQIIndex ).set( "vectorial_field_index", 1 );

        // velocity VZ
        aParameterLists( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIBulkVZ" );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists( tIQIIndex ).set( "dof_quantity", "VX,VY,VZ" );
        aParameterLists( tIQIIndex ).set( "vectorial_field_index", 2 );

        // pressure
        aParameterLists( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIBulkP" );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists( tIQIIndex ).set( "dof_quantity", "P" );
        aParameterLists( tIQIIndex ).set( "vectorial_field_index", 0 );

        // temperature
        aParameterLists( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIBulkTEMP" );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseAll" );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists( tIQIIndex ).set( "dof_quantity", "TEMP" );
        aParameterLists( tIQIIndex ).set( "vectorial_field_index", 0 );

        // fluid thermal energy on inlet
        aParameterLists( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIInletThermalEnergy" );
        aParameterLists( tIQIIndex ).set( "IQI_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIQIIndex ).set( "neighbor_phases", "PhaseVoidFront" );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::THERMAL_ENERGY_CONVECTIVE_FLUX );
        aParameterLists( tIQIIndex ).set( "leader_constitutive_models", "CMFluid,Fluid" );

        // fluid thermal energy on outlet
        aParameterLists( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIOutletThermalEnergy" );
        aParameterLists( tIQIIndex ).set( "IQI_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::THERMAL_ENERGY_CONVECTIVE_FLUX );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIQIIndex ).set( "neighbor_phases", "PhaseVoidBack" );
        aParameterLists( tIQIIndex ).set( "leader_constitutive_models", "CMFluid,Fluid" );

        // fluid total pressure on inlet
        aParameterLists( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIInletTotalPressure" );
        aParameterLists( tIQIIndex ).set( "IQI_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIQIIndex ).set( "neighbor_phases", "PhaseVoidFront" );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::TOTAL_PRESSURE );
        aParameterLists( tIQIIndex ).set( "leader_constitutive_models", "CMFluid,Fluid" );

        // fluid total pressure on outlet
        aParameterLists( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIOutletTotalPressure" );
        aParameterLists( tIQIIndex ).set( "IQI_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIQIIndex ).set( "neighbor_phases", "PhaseVoidBack" );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::TOTAL_PRESSURE );
        aParameterLists( tIQIIndex ).set( "leader_constitutive_models", "CMFluid,Fluid" );

        // fluid mass flow on inlet
        aParameterLists( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIInletMassFlow" );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::MASS_FLOW );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIQIIndex ).set( "neighbor_phases", "PhaseVoidFront" );
        aParameterLists( tIQIIndex ).set( "IQI_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists( tIQIIndex ).set( "leader_constitutive_models", "CMFluid,Fluid" );

        // fluid mass flow on outlet
        aParameterLists( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIOutletMassFlow" );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIQIIndex ).set( "neighbor_phases", "PhaseVoidBack" );
        aParameterLists( tIQIIndex ).set( "IQI_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::MASS_FLOW );
        aParameterLists( tIQIIndex ).set( "leader_constitutive_models", "CMFluid,Fluid" );

        // inclusion perimeter
        aParameterLists( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIPerimeterItf" );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseSolid" );
        aParameterLists( tIQIIndex ).set( "neighbor_phases", "PhaseFluid" );
        aParameterLists( tIQIIndex ).set( "IQI_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::VOLUME );

        // max temperature in solid
        aParameterLists( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIMaxTemp" );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseSolid" );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::MAX_DOF );
        aParameterLists( tIQIIndex ).set( "function_parameters", tRefTemp + "/" + tMaxTempExponent );
        aParameterLists( tIQIIndex ).set( "dof_quantity", "TEMP" );
        aParameterLists( tIQIIndex ).set( "vectorial_field_index", 0 );

        // volume solid
        aParameterLists( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQISolidVolume" );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseSolid" );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::VOLUME );

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        aParameterLists( tFEMIndex ).push_back( prm::create_computation_parameter_list() );
        aParameterLists( tFEMIndex ).set( "print_physics_model", false );

        aParameterLists( tFEMIndex ).set( "finite_difference_scheme", fem::FDScheme_Type::POINT_3_CENTRAL );
        aParameterLists( tFEMIndex ).set( "finite_difference_perturbation_size", tFEMFdEpsilon );
    }

    void
    SOLParameterList( Module_Parameter_Lists& aParameterLists )
    {


        aParameterLists( 0 ).push_back( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL ) );
        aParameterLists( 0 ).push_back( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL ) );

#ifdef MORIS_USE_MUMPS
        aParameterLists( 0 ).set( "Solver_Type", "Amesos_Mumps" );
        aParameterLists( 0 ).set( "Solver_Type", "Amesos_Mumps" );
#else
        aParameterLists( 0 ).set( "Solver_Type", "Amesos_Superludist" );
        aParameterLists( 0 ).set( "Solver_Type", "Amesos_Superludist" );
#endif

        /*
        aParameterLists( 0 ).push_back( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::BELOS_IMPL ) );
        aParameterLists( 0 ).set( "ifpack_prec_type", "ILU");
        aParameterLists( 0 ).set( "fact: level-of-fill", 1);

        aParameterLists( 0 ).push_back( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::BELOS_IMPL ) );
        aParameterLists( 0 ).set( "ifpack_prec_type", "ILU");
        aParameterLists( 0 ).set( "fact: level-of-fill", 3);
        */

        //------------------------------------------------------------------------------

        aParameterLists( 1 ).push_back( moris::prm::create_linear_solver_parameter_list() );
        aParameterLists( 1 ).set( "DLA_Linear_solver_algorithms", "0" );

        aParameterLists( 1 ).push_back( moris::prm::create_linear_solver_parameter_list() );
        aParameterLists( 1 ).set( "DLA_Linear_solver_algorithms", "1" );

        //------------------------------------------------------------------------------

        aParameterLists( 2 ).push_back( moris::prm::create_nonlinear_algorithm_parameter_list() );    // nonlinear algorithm index 0
        aParameterLists( 2 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists( 2 ).set( "NLA_Linear_solver", 0 );
        aParameterLists( 2 ).set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        aParameterLists( 2 ).set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
        aParameterLists( 2 ).set( "NLA_max_iter", tNLA_max_iter );

        aParameterLists( 2 ).push_back( moris::prm::create_nonlinear_algorithm_parameter_list() );
        aParameterLists( 2 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists( 2 ).set( "NLA_Linear_solver", 0 );
        aParameterLists( 2 ).set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        aParameterLists( 2 ).set( "NLA_relaxation_parameter", 1.0 );
        aParameterLists( 2 ).set( "NLA_max_iter", 1 );

        aParameterLists( 2 ).push_back( moris::prm::create_nonlinear_algorithm_parameter_list() );    // nonlinear algorithm index 1
        aParameterLists( 2 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NLBGS_SOLVER );
        aParameterLists( 2 ).set( "NLA_Linear_solver", 0 );
        aParameterLists( 2 ).set( "NLA_rel_res_norm_drop", 1.0 );
        aParameterLists( 2 ).set( "NLA_max_iter", 1 );

        aParameterLists( 2 ).push_back( moris::prm::create_nonlinear_algorithm_parameter_list() );
        aParameterLists( 2 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists( 2 ).set( "NLA_Linear_solver", 1 );
        aParameterLists( 2 ).set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        aParameterLists( 2 ).set( "NLA_relaxation_parameter", 1.0 );
        aParameterLists( 2 ).set( "NLA_max_iter", 1 );

        //------------------------------------------------------------------------------

        aParameterLists( 3 ).push_back( moris::prm::create_nonlinear_solver_parameter_list() );    // 0: fluid subproblem
        aParameterLists( 3 ).set( "NLA_Nonlinear_solver_algorithms", "0" );             // set nonlinear algorithm with index 0
        aParameterLists( 3 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists( 3 ).set( "NLA_DofTypes", "VX,VY,VZ,P" );

        aParameterLists( 3 ).push_back( moris::prm::create_nonlinear_solver_parameter_list() );    // 1: thermal subproblem
        aParameterLists( 3 ).set( "NLA_Nonlinear_solver_algorithms", "1" );             // set nonlinear algorithm with index 0
        aParameterLists( 3 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists( 3 ).set( "NLA_DofTypes", "TEMP" );

        aParameterLists( 3 ).push_back( moris::prm::create_nonlinear_solver_parameter_list() );    // 2: one-way coupling via NLBGS
        aParameterLists( 3 ).set( "NLA_Nonlinear_solver_algorithms", "2" );             // set nonlinear algorithm with index 1.
        aParameterLists( 3 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NLBGS_SOLVER );
        aParameterLists( 3 ).set( "NLA_Sub_Nonlinear_Solver", "0,1" );    // set sub nonlinear solvers with index 0 and 1
        aParameterLists( 3 ).set( "NLA_DofTypes", "VX,VY,VZ,P;TEMP" );

        aParameterLists( 3 ).push_back( moris::prm::create_nonlinear_solver_parameter_list() );
        aParameterLists( 3 ).set( "NLA_DofTypes", "VX,VY,VZ,P,TEMP" );
        aParameterLists( 3 ).set( "NLA_Nonlinear_solver_algorithms", "3" );

        // ----------------------------------------------------------

        aParameterLists( 4 ).push_back( moris::prm::create_time_solver_algorithm_parameter_list() );
        aParameterLists( 4 ).set( "TSA_Nonlinear_Solver", 2 );                // using NLBGS for forward problem
        aParameterLists( 4 ).set( "TSA_Nonlinear_Sensitivity_Solver", 3 );    // using monlithic for sensitivity problem

        //------------------------------------------------------------------------------

        aParameterLists( 5 ).push_back( moris::prm::create_time_solver_parameter_list() );
        aParameterLists( 5 ).set( "TSA_DofTypes", "VX,VY,VZ;P;TEMP" );
        aParameterLists( 5 ).set( "TSA_Initialize_Sol_Vec", "VX,0.01;VY,0.0;VZ,0.0;P,0.0;TEMP,0.0" );
        aParameterLists( 5 ).set( "TSA_Output_Indices", "0" );
        aParameterLists( 5 ).set( "TSA_Output_Criteria", "Output_Criterion" );

        //------------------------------------------------------------------------------

        aParameterLists( 6 ).push_back( moris::prm::create_solver_warehouse_parameterlist() );

        aParameterLists( 7 ).push_back( moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::NONE ) );
    }

    void
    MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).push_back( prm::create_msi_parameter_list() );
        aParameterLists( 0 ).set( "VX", 1 );
        aParameterLists( 0 ).set( "VY", 1 );
        aParameterLists( 0 ).set( "VZ", 1 );
        aParameterLists( 0 ).set( "TEMP", 1 );
    }

    void
    VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).push_back( prm::create_vis_parameter_list() );
        aParameterLists( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", tExoFile ) );
        aParameterLists( 0 ).set( "Mesh_Type", vis::VIS_Mesh_Type::STANDARD );
        aParameterLists( 0 ).set( "Set_Names", tAllPhases );
        aParameterLists( 0 ).set( "Field_Names", "VX,VY,VZ,P,TEMP" );
        aParameterLists( 0 ).set( "Field_Type", "NODAL,NODAL,NODAL,NODAL,NODAL" );
        aParameterLists( 0 ).set( "IQI_Names", "IQIBulkVX,IQIBulkVY,IQIBulkVZ,IQIBulkP,IQIBulkTEMP" );
        aParameterLists( 0 ).set( "Save_Frequency", 1 );
        aParameterLists( 0 ).set( "Time_Offset", 10.0 );
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
