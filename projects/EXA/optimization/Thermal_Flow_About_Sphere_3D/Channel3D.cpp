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

    std::string tInterfaceRefinementSphere = "0";
    std::string tInterfaceRefinementWalls  = "0";

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
    uint tFEMFdScheme  = static_cast< uint >( fem::FDScheme_Type::POINT_3_CENTRAL );

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
            const moris::Matrix< DDRMat >&     aCoordinates,
            const moris::Cell< real >& aGeometryParameters )
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
            const moris::Matrix< moris::DDRMat >&                aCoordinates,
            const moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::Matrix< DDRMat >&                             aFieldSensitivity )
    {
        // derivative of level set function wrt sphere radius
        aFieldSensitivity = { { 1.0 } };
    }

    /* ------------------------------------------------------------------------ */
    // geometry parameters & LS functions

    moris::real
    Func_Plane(
            const moris::Matrix< DDRMat >&     aCoordinates,
            const moris::Cell< real >& aGeometryParameters )
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
            Matrix< DDRMat > aADVs,
            Matrix< DDRMat > aCriteria )
    {
        Matrix< DDRMat > tObjectives = { { aCriteria( 0 ) } };

        return tObjectives;
    }

    /* ------------------------------------------------------------------------ */

    Matrix< DDRMat >
    compute_dobjective_dadv(
            Matrix< DDRMat > aADVs,
            Matrix< DDRMat > aCriteria )
    {
        Matrix< DDRMat > tDObjectiveDADV( 1, aADVs.numel(), 0.0 );

        return tDObjectiveDADV;
    }

    /* ------------------------------------------------------------------------ */

    Matrix< DDRMat >
    compute_dobjective_dcriteria(
            Matrix< DDRMat > aADVs,
            Matrix< DDRMat > aCriteria )
    {
        Matrix< DDRMat > tDObjectiveDCriteria( 1, aCriteria.numel(), 0.0 );

        tDObjectiveDCriteria( 0 ) = 1.0;

        return tDObjectiveDCriteria;
    }

    /* ------------------------------------------------------------------------ */

    Matrix< DDRMat >
    compute_constraints(
            Matrix< DDRMat > aADVs,
            Matrix< DDRMat > aCriteria )
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
            Matrix< DDRMat > aADVs,
            Matrix< DDRMat > aCriteria )
    {
        Matrix< DDRMat > tDConstraintDADV( tNumConstraints, aADVs.numel(), 0.0 );

        return tDConstraintDADV;
    }

    /* ------------------------------------------------------------------------ */

    Matrix< DDRMat >
    compute_dconstraint_dcriteria(
            Matrix< DDRMat > aADVs,
            Matrix< DDRMat > aCriteria )
    {
        Matrix< DDRMat > tDConstraintDCriteria( tNumConstraints, aCriteria.numel(), 0.0 );

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
    OPTParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 3 );

        tParameterlist( 0 ).push_back( prm::create_opt_problem_parameter_list() );
        tParameterlist( 0 )( 0 ).set( "is_optimization_problem", true );
        tParameterlist( 0 )( 0 ).set( "problem", "user_defined" );
        tParameterlist( 0 )( 0 ).set( "library", tSoFile );

        tParameterlist( 1 ).resize( 0 );

        tParameterlist( 2 ).push_back( moris::prm::create_sweep_parameter_list() );
        tParameterlist( 2 )( 0 ).set( "hdf5_path", tHdf5File );
        tParameterlist( 2 )( 0 ).set( "num_evaluations_per_adv", "1" );
        tParameterlist( 2 )( 0 ).set( "finite_difference_type", "all" );
        tParameterlist( 2 )( 0 ).set( "finite_difference_epsilons", tSweepFdEpsilon );
    }

    /* ------------------------------------------------------------------------ */

    void
    HMRParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );

        tParameterlist( 0 ).push_back( prm::create_hmr_parameter_list() );

        tParameterlist( 0 )( 0 ).set( "number_of_elements_per_dimension", tNumElemsPerDim );
        tParameterlist( 0 )( 0 ).set( "domain_dimensions", tDomainDims );
        tParameterlist( 0 )( 0 ).set( "domain_offset", tDomainOffset );
        tParameterlist( 0 )( 0 ).set( "domain_sidesets", tDomainSidesets );

        tParameterlist( 0 )( 0 ).set( "lagrange_output_meshes", "0" );

        tParameterlist( 0 )( 0 ).set( "lagrange_orders", tLagrangeOrder );
        tParameterlist( 0 )( 0 ).set( "lagrange_pattern", tLagrangePattern );

        tParameterlist( 0 )( 0 ).set( "bspline_orders", tBsplineOrder );
        tParameterlist( 0 )( 0 ).set( "bspline_pattern", "0,1" );

        tParameterlist( 0 )( 0 ).set( "lagrange_to_bspline", "0,1" );

        tParameterlist( 0 )( 0 ).set( "truncate_bsplines", 1 );
        tParameterlist( 0 )( 0 ).set( "refinement_buffer", tRefineBuffer );
        tParameterlist( 0 )( 0 ).set( "staircase_buffer", tRefineBuffer );

        tParameterlist( 0 )( 0 ).set( "initial_refinement", tInitialRef );
        tParameterlist( 0 )( 0 ).set( "initial_refinement_pattern", "0,1" );

        tParameterlist( 0 )( 0 ).set( "use_number_aura", 1 );

        tParameterlist( 0 )( 0 ).set( "use_multigrid", 0 );
        tParameterlist( 0 )( 0 ).set( "severity_level", 0 );
    }

    /* ------------------------------------------------------------------------ */

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
        tParameterlist( 0 )( 0 ).set( "enrich_mesh_indices", "0,1" );
        tParameterlist( 0 )( 0 ).set( "ghost_stab", tUseGhost );
        tParameterlist( 0 )( 0 ).set( "multigrid", false );
        tParameterlist( 0 )( 0 ).set( "verbose", false );
        tParameterlist( 0 )( 0 ).set( "high_to_low_dbl_side_sets", true );
        tParameterlist( 0 )( 0 ).set( "print_enriched_ig_mesh", false );
        tParameterlist( 0 )( 0 ).set( "exodus_output_XTK_ig_mesh", true );
    }

    /* ------------------------------------------------------------------------ */

    void
    GENParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 3 );

        tParameterlist( 0 ).push_back( prm::create_gen_parameter_list() );
        tParameterlist( 0 )( 0 ).set( "IQI_types",
                "IQIInletThermalEnergy,IQIOutletThermalEnergy,"
                "IQIInletTotalPressure,IQIOutletTotalPressure,"
                "IQIPerimeterItf,"
                "IQIInletMassFlow,IQIOutletMassFlow,"
                "IQIMaxTemp,"
                "IQISolidVolume" );

        tParameterlist( 0 )( 0 ).set( "initial_advs", moris_to_string( tSphereRadius ) );
        tParameterlist( 0 )( 0 ).set( "lower_bounds", moris_to_string( tSphereRadius * 0.9 ) );
        tParameterlist( 0 )( 0 ).set( "upper_bounds", moris_to_string( tSphereRadius * 1.1 ) );

        tParameterlist( 0 )( 0 ).set( "number_of_phases", 5 );
        tParameterlist( 0 )( 0 ).set( "phase_function_name", "get_phase_index" );
        tParameterlist( 0 )( 0 ).set( "output_mesh_file", tGENOutputFile );
        tParameterlist( 0 )( 0 ).set( "time_offset", 10.0 );

        // init geometry counter
        uint tGeoCounter = 0;

        // Inclusions
        tParameterlist( 1 ).push_back( prm::create_user_defined_geometry_parameter_list() );
        tParameterlist( 1 )( tGeoCounter ).set( "field_function_name", "Func_Sphere" );
        tParameterlist( 1 )( tGeoCounter ).set( "sensitivity_function_name", "Func_Sphere_Deriv" );
        tParameterlist( 1 )( tGeoCounter ).set( "number_of_refinements", tInterfaceRefinementSphere );
        tParameterlist( 1 )( tGeoCounter ).set( "refinement_mesh_index", "0" );
        tParameterlist( 1 )( tGeoCounter ).set( "isocontour_threshold", 0.0 );
        tParameterlist( 1 )( tGeoCounter ).set( "isocontour_tolerance", 1.0e-12 );
        tParameterlist( 1 )( tGeoCounter ).set( "intersection_tolerance", 1.0e-12 );

        tParameterlist( 1 )( tGeoCounter ).set( "field_variable_indices", "0" );
        tParameterlist( 1 )( tGeoCounter ).set( "adv_indices", "0" );
        tGeoCounter++;

        // Inlet plane
        tParameterlist( 1 ).push_back( prm::create_user_defined_geometry_parameter_list() );
        tParameterlist( 1 )( tGeoCounter ).set( "field_function_name", "Func_Plane" );
        tParameterlist( 1 )( tGeoCounter ).set( "number_of_refinements", tInterfaceRefinementWalls );
        tParameterlist( 1 )( tGeoCounter ).set( "refinement_mesh_index", "0" );
        tParameterlist( 1 )( tGeoCounter ).set( "constant_parameters", "1,0,0," + moris_to_string( tGeoShift ) + ",0,0" );
        tParameterlist( 1 )( tGeoCounter ).set( "isocontour_threshold", 0.0 );
        tParameterlist( 1 )( tGeoCounter ).set( "isocontour_tolerance", 1.0e-12 );
        tParameterlist( 1 )( tGeoCounter ).set( "intersection_tolerance", 1.0e-12 );
        tGeoCounter++;

        // Out plane
        tParameterlist( 1 ).push_back( prm::create_user_defined_geometry_parameter_list() );
        tParameterlist( 1 )( tGeoCounter ).set( "field_function_name", "Func_Plane" );
        tParameterlist( 1 )( tGeoCounter ).set( "number_of_refinements", tInterfaceRefinementWalls );
        tParameterlist( 1 )( tGeoCounter ).set( "refinement_mesh_index", "0" );
        tParameterlist( 1 )( tGeoCounter ).set( "constant_parameters", "-1,0,0," + moris_to_string( tChannelLength + tGeoShift ) + ",0,0" );
        tParameterlist( 1 )( tGeoCounter ).set( "isocontour_threshold", 0.0 );
        tParameterlist( 1 )( tGeoCounter ).set( "isocontour_tolerance", 1.0e-12 );
        tParameterlist( 1 )( tGeoCounter ).set( "intersection_tolerance", 1.0e-12 );
        tGeoCounter++;

        // Lower plane
        tParameterlist( 1 ).push_back( prm::create_user_defined_geometry_parameter_list() );
        tParameterlist( 1 )( tGeoCounter ).set( "field_function_name", "Func_Plane" );
        tParameterlist( 1 )( tGeoCounter ).set( "number_of_refinements", tInterfaceRefinementWalls );
        tParameterlist( 1 )( tGeoCounter ).set( "refinement_mesh_index", "0" );
        tParameterlist( 1 )( tGeoCounter ).set( "constant_parameters", "0,1,0,0," + moris_to_string( tGeoShift ) + ",0" );
        tParameterlist( 1 )( tGeoCounter ).set( "isocontour_threshold", 0.0 );
        tParameterlist( 1 )( tGeoCounter ).set( "isocontour_tolerance", 1.0e-12 );
        tParameterlist( 1 )( tGeoCounter ).set( "intersection_tolerance", 1.0e-12 );
        tGeoCounter++;

        // Upper plane
        tParameterlist( 1 ).push_back( prm::create_user_defined_geometry_parameter_list() );
        tParameterlist( 1 )( tGeoCounter ).set( "field_function_name", "Func_Plane" );
        tParameterlist( 1 )( tGeoCounter ).set( "number_of_refinements", tInterfaceRefinementWalls );
        tParameterlist( 1 )( tGeoCounter ).set( "refinement_mesh_index", "0" );
        tParameterlist( 1 )( tGeoCounter ).set( "constant_parameters", "0,-1,0,0," + moris_to_string( tChannelHeight + tGeoShift ) + ",0" );
        tParameterlist( 1 )( tGeoCounter ).set( "isocontour_threshold", 0.0 );
        tParameterlist( 1 )( tGeoCounter ).set( "isocontour_tolerance", 1.0e-12 );
        tParameterlist( 1 )( tGeoCounter ).set( "intersection_tolerance", 1.0e-12 );
        tGeoCounter++;

        // Back Plan
        tParameterlist( 1 ).push_back( prm::create_user_defined_geometry_parameter_list() );
        tParameterlist( 1 )( tGeoCounter ).set( "field_function_name", "Func_Plane" );
        tParameterlist( 1 )( tGeoCounter ).set( "number_of_refinements", tInterfaceRefinementWalls );
        tParameterlist( 1 )( tGeoCounter ).set( "refinement_mesh_index", "0" );
        tParameterlist( 1 )( tGeoCounter ).set( "constant_parameters", "0,0,1,0,0," + moris_to_string( tGeoShift ) );
        tParameterlist( 1 )( tGeoCounter ).set( "isocontour_threshold", 0.0 );
        tParameterlist( 1 )( tGeoCounter ).set( "isocontour_tolerance", 1.0e-12 );
        tParameterlist( 1 )( tGeoCounter ).set( "intersection_tolerance", 1.0e-12 );
        tGeoCounter++;

        // Front plane
        tParameterlist( 1 ).push_back( prm::create_user_defined_geometry_parameter_list() );
        tParameterlist( 1 )( tGeoCounter ).set( "field_function_name", "Func_Plane" );
        tParameterlist( 1 )( tGeoCounter ).set( "number_of_refinements", tInterfaceRefinementWalls );
        tParameterlist( 1 )( tGeoCounter ).set( "refinement_mesh_index", "0" );
        tParameterlist( 1 )( tGeoCounter ).set( "constant_parameters", "0,0,-1,0,0," + moris_to_string( tChannelWidth + tGeoShift ) );
        tParameterlist( 1 )( tGeoCounter ).set( "isocontour_threshold", 0.0 );
        tParameterlist( 1 )( tGeoCounter ).set( "isocontour_tolerance", 1.0e-12 );
        tParameterlist( 1 )( tGeoCounter ).set( "intersection_tolerance", 1.0e-12 );
        tGeoCounter++;
    }

    /* ------------------------------------------------------------------------ */

    void
    FEMParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterList )
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
        tParameterList.resize( 9 );
        uint tPropIndex  = 0;
        uint tCMIndex    = 1;
        uint tSPIndex    = 2;
        uint tIWGIndex   = 3;
        uint tIQIIndex   = 4;
        uint tFEMIndex   = 5;
        uint tPhaseIndex = 7;

        //------------------------------------------------------------------------------
        // phase info
        uint tPhaseCounter = 0;

        tParameterList( tPhaseIndex ).push_back( prm::create_phase_parameter_list() );
        tParameterList( tPhaseIndex )( tPhaseCounter ).set( "phase_name", "PhaseFluid" );
        tParameterList( tPhaseIndex )( tPhaseCounter ).set( "phase_indices", "0" );
        tPhaseCounter++;

        tParameterList( tPhaseIndex ).push_back( prm::create_phase_parameter_list() );
        tParameterList( tPhaseIndex )( tPhaseCounter ).set( "phase_name", "PhaseSolid" );
        tParameterList( tPhaseIndex )( tPhaseCounter ).set( "phase_indices", "1" );
        tPhaseCounter++;

        tParameterList( tPhaseIndex ).push_back( prm::create_phase_parameter_list() );
        tParameterList( tPhaseIndex )( tPhaseCounter ).set( "phase_name", "PhaseVoidFront" );
        tParameterList( tPhaseIndex )( tPhaseCounter ).set( "phase_indices", "2" );
        tPhaseCounter++;

        tParameterList( tPhaseIndex ).push_back( prm::create_phase_parameter_list() );
        tParameterList( tPhaseIndex )( tPhaseCounter ).set( "phase_name", "PhaseVoidBack" );
        tParameterList( tPhaseIndex )( tPhaseCounter ).set( "phase_indices", "3" );
        tPhaseCounter++;

        tParameterList( tPhaseIndex ).push_back( prm::create_phase_parameter_list() );
        tParameterList( tPhaseIndex )( tPhaseCounter ).set( "phase_name", "PhaseVoidLateral" );
        tParameterList( tPhaseIndex )( tPhaseCounter ).set( "phase_indices", "4" );
        tPhaseCounter++;

        tParameterList( tPhaseIndex ).push_back( prm::create_phase_parameter_list() );
        tParameterList( tPhaseIndex )( tPhaseCounter ).set( "phase_name", "PhaseAll" );
        tParameterList( tPhaseIndex )( tPhaseCounter ).set( "phase_indices", "0,1" );
        tPhaseCounter++;

        //------------------------------------------------------------------------------
        // fill the property part of the parameter list

        // init property counter
        uint tPropCounter = 0;

        // fluid properties ------------------------------------------------------------
        // create fluid viscosity property
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropFluidViscosity" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", tFluidViscosity );
        tPropCounter++;

        // create fluid density property
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropFluidDensity" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", tFluidDensity );
        tPropCounter++;

        // create fluid capacity property
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropFluidCapacity" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", tFluidCapacity );
        tPropCounter++;

        // create fluid conductivity property
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropFluidConductivity" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", tFluidConductivity );

        tPropCounter++;

        // create fluid pressure spring property
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropFluidPressureSpring" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", tFluidPressureSpring );
        tPropCounter++;

        // solid properties ----------------------------------------------------------
        // create solid B density property
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropSolidDensity" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", tSolidDensity );
        tPropCounter++;

        // create solid B capacity property
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropSolidCapacity" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", tSolidCapacity );
        tPropCounter++;

        // create solid B conductivity property
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropSolidConductivity" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", tSolidConductivity );
        tPropCounter++;

        // BC properties ---------------------------------------------------------------
        // create inlet pressure property
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropInletPressure" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", tInletPressure );
        tPropCounter++;

        // create inlet temperature property
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropInletTemp" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", tInletTemperature );
        tPropCounter++;

        // create wall velocity property
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropZeroU" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", "0.0;0.0;0.0" );
        tPropCounter++;

        // create symmetry velocity property (x-z plane)
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropSelectUY" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", "0.0,0.0,0.0;0.0,1.0,0.0;0.0,0.0,0.0" );
        tPropCounter++;

        // create symmetry velocity property (x-y plane)
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropSelectUZ" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", "0.0,0.0,0.0;0.0,0.0,0.0;0.0,0.0,1.0" );
        tPropCounter++;

        // create heat load property
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropVolumetricHeatFlux" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", tVolumetricHeatLoad );
        tPropCounter++;

        //------------------------------------------------------------------------------
        // fill the constitutive model part of the parameter list

        // init CM counter
        uint tCMCounter = 0;

        // fluid CM --------------------------------------------------------------------
        // create fluid CM
        tParameterList( tCMIndex ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( tCMIndex )( tCMCounter ).set( "constitutive_name", "CMFluid" );
        tParameterList( tCMIndex )( tCMCounter ).set( "phase_name", "PhaseFluid" );
        tParameterList( tCMIndex )( tCMCounter ).set( "constitutive_type", (uint)fem::Constitutive_Type::FLUID_INCOMPRESSIBLE );
        tParameterList( tCMIndex )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( "VX,VY,VZ;P", "Velocity,Pressure" ) );
        tParameterList( tCMIndex )( tCMCounter ).set( "properties",
                "PropFluidViscosity,Viscosity;"
                "PropFluidDensity  ,Density" );
        tCMCounter++;

        // create fluid diffusion CM
        tParameterList( tCMIndex ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( tCMIndex )( tCMCounter ).set( "constitutive_name", "CMFluidDiffusion" );
        tParameterList( tCMIndex )( tCMCounter ).set( "phase_name", "PhaseFluid" );
        tParameterList( tCMIndex )( tCMCounter ).set( "constitutive_type", (uint)fem::Constitutive_Type::DIFF_LIN_ISO );
        tParameterList( tCMIndex )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( tCMIndex )( tCMCounter ).set( "properties",
                "PropFluidConductivity,Conductivity;"
                "PropFluidDensity     ,Density;"
                "PropFluidCapacity    ,HeatCapacity" );
        tCMCounter++;

        // solid CM --------------------------------------------------------------------
        // create solid B diffusion CM
        tParameterList( tCMIndex ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( tCMIndex )( tCMCounter ).set( "constitutive_name", "CMSolidDiffusion" );
        tParameterList( tCMIndex )( tCMCounter ).set( "phase_name", "PhaseSolid" );
        tParameterList( tCMIndex )( tCMCounter ).set( "constitutive_type", (uint)fem::Constitutive_Type::DIFF_LIN_ISO );
        tParameterList( tCMIndex )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( tCMIndex )( tCMCounter ).set( "properties",
                "PropSolidConductivity,Conductivity;"
                "PropSolidDensity     ,Density;"
                "PropSolidCapacity    ,HeatCapacity" );
        tCMCounter++;

        //------------------------------------------------------------------------------
        // fill the stabilization parameter part of the parameter list

        // init SP counter
        uint tSPCounter = 0;

        // create SUPG/PSG fluid
        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name", "SPIncFlow" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type", (uint)fem::Stabilization_Type::INCOMPRESSIBLE_FLOW );
        tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters", tCInv );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY,VZ;P", "Velocity,Pressure" ) );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_properties",
                "PropFluidViscosity,Viscosity;"
                "PropFluidDensity  ,Density" );
        tSPCounter++;

        // create SUPG fluid temperature advection
        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name", "SPSUPGTemp" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type", (uint)fem::Stabilization_Type::SUPG_ADVECTION );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY,VZ", "Velocity" ) );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_properties", "PropFluidConductivity,Conductivity" );
        tSPCounter++;

        // create Nitsche for fluid velocity
        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name", "SPNitscheU" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type", (uint)fem::Stabilization_Type::VELOCITY_DIRICHLET_NITSCHE );
        tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters", "100.0/1.0" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY,VZ", "Velocity" ) );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_properties",
                "PropFluidViscosity,Viscosity;"
                "PropFluidDensity,Density" );
        tSPCounter++;

        // create Nitsche for fluid temperature
        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name", "SPNitscheT" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type", (uint)fem::Stabilization_Type::DIRICHLET_NITSCHE );
        tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters", "100.0" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_properties", "PropFluidConductivity,Material" );
        tSPCounter++;

        // create Nitsche for fluid/solid diffusion interface
        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name", "SPFSInterfaceNitsche" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tSPIndex )( tSPCounter ).set( "follower_phase_name", "PhaseSolid" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type", (uint)fem::Stabilization_Type::NITSCHE_INTERFACE );
        tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters", "100.0" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_properties", "PropFluidConductivity,Material" );
        tParameterList( tSPIndex )( tSPCounter ).set( "follower_properties", "PropSolidConductivity,Material" );
        tSPCounter++;

        // create ghost penalty viscous
        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name", "SPGPViscous" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type", (uint)fem::Stabilization_Type::VISCOUS_GHOST );
        tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters", "0.05" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_properties", "PropFluidViscosity,Viscosity" );
        tSPCounter++;

        // create ghost penalty convective
        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name", "SPGPVelocity" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type", (uint)fem::Stabilization_Type::CONVECTIVE_GHOST );
        tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters", "0.05" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY,VZ", "Velocity" ) );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_properties", "PropFluidDensity,Density" );
        tSPCounter++;

        // create ghost penalty pressure
        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name", "SPGPPressure" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type", (uint)fem::Stabilization_Type::PRESSURE_GHOST );
        tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters", "0.005/1.0" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY,VZ", "Velocity" ) );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_properties",
                "PropFluidViscosity,Viscosity;"
                "PropFluidDensity,Density" );
        tSPCounter++;

        // create ghost penalty fluid temperature
        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name", "SPGPFluidTemp" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type", (uint)fem::Stabilization_Type::GHOST_DISPL );
        tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters", "0.05" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_properties", "PropFluidConductivity,Material" );
        tSPCounter++;

        // create ghost penalty solid temperature
        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name", "SPGPSolidTemp" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_phase_name", "PhaseSolid" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type", (uint)fem::Stabilization_Type::GHOST_DISPL );
        tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters", "0.05" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_properties", "PropSolidConductivity,Material" );
        tSPCounter++;

        //------------------------------------------------------------------------------
        // fill the IWG part of the parameter list

        // init IWG counter
        uint tIWGCounter = 0;

        // fluid bulk IWGs -------------------------------------------------------------
        // NS incompressible (velocity)
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGVelocityBulk" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", (uint)fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_BULK );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "VX,VY,VZ" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SPIncFlow,IncompressibleFlow" );
        tIWGCounter++;

        // NS incompressible (pressure)
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGPressureBulk" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", (uint)fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_BULK );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "P" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_properties", "PropFluidPressureSpring,PressureSpring" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SPIncFlow,IncompressibleFlow" );
        tIWGCounter++;

        // diffusion
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGFluidDiffusionBulk" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", (uint)fem::IWG_Type::SPATIALDIFF_BULK );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "TEMP" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", "CMFluidDiffusion,Diffusion" );
        tIWGCounter++;

        // advection
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGFluidAdvectionBulk" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", (uint)fem::IWG_Type::ADVECTION_BULK );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "TEMP" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", "CMFluidDiffusion,Diffusion" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SPSUPGTemp,SUPG" );
        tIWGCounter++;

        // solid bulk IWGs -----------------------------------------------------------
        // diffusion
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGSolidDiffusionBulk" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseSolid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", (uint)fem::IWG_Type::SPATIALDIFF_BULK );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "TEMP" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_properties", "PropVolumetricHeatFlux,Load" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", "CMSolidDiffusion,Diffusion" );
        tIWGCounter++;

        // fluid/solid ifc IWGs ------------------------------------------------------

        // temperature
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGInterfaceFluidSolid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", (uint)fem::Element_Type::DOUBLE_SIDESET );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "follower_phase_name", "PhaseSolid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", (uint)fem::IWG_Type::SPATIALDIFF_INTERFACE_UNSYMMETRIC_NITSCHE );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "TEMP" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", "CMFluidDiffusion,Diffusion" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "follower_constitutive_models", "CMSolidDiffusion,Diffusion" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SPFSInterfaceNitsche,NitscheInterface" );
        tIWGCounter++;

        // zero velocity
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGZeroVelocity" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", (uint)fem::Element_Type::SIDESET );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "neighbor_phases", "PhaseSolid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", (uint)fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_UNSYMMETRIC_NITSCHE );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "VX,VY,VZ" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_properties", "PropZeroU,Dirichlet" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SPNitscheU,DirichletNitsche" );
        tIWGCounter++;

        // zero velocity (pressure part)
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGZeroPressure" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", (uint)fem::Element_Type::SIDESET );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "neighbor_phases", "PhaseSolid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", (uint)fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_UNSYMMETRIC_NITSCHE );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "P" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_properties", "PropZeroU,Dirichlet" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        tIWGCounter++;

        // Laterial BC IWG (adiabatic) ----------------------------------------------------------------

        // zero velocity along lateral sides
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGZeroVelocityVoid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", (uint)fem::Element_Type::SIDESET );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "neighbor_phases", "PhaseVoidLateral" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", (uint)fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_UNSYMMETRIC_NITSCHE );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "VX,VY,VZ" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_properties", "PropZeroU,Dirichlet" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SPNitscheU,DirichletNitsche" );
        tIWGCounter++;

        // zero velocity along lateral sides (pressure part)
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGZeroPressureVoid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", (uint)fem::Element_Type::SIDESET );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "neighbor_phases", "PhaseVoidLateral" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", (uint)fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_UNSYMMETRIC_NITSCHE );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "P" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_properties", "PropZeroU,Dirichlet" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        tIWGCounter++;

        // Inlet BC IWG ----------------------------------------------------------------

        // inlet pressure
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGInletImposedPressure" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "neighbor_phases", "PhaseVoidFront" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", (uint)fem::Element_Type::SIDESET );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", (uint)fem::IWG_Type::INCOMPRESSIBLE_NS_IMPOSED_PRESSURE );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "VX,VY,VZ" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_properties", "PropInletPressure,Pressure" );
        tIWGCounter++;

        // inlet temperature
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGInletTemp" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", (uint)fem::Element_Type::SIDESET );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "neighbor_phases", "PhaseVoidFront" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", (uint)fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "TEMP" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_properties", "PropInletTemp,Dirichlet" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", "CMFluidDiffusion,Diffusion" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SPNitscheT,DirichletNitsche" );
        tIWGCounter++;

        // Ghost  ----------------------------------------------------------------

        if ( tUseGhost )
        {
            // ghost viscous
            tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGGPViscous" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", (uint)fem::Element_Type::DOUBLE_SIDESET );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFluid" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "follower_phase_name", "PhaseFluid" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", (uint)fem::IWG_Type::GHOST_NORMAL_FIELD );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "VX,VY,VZ" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SPGPViscous,GhostSP" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "ghost_order", (uint)tDispOrder );
            tIWGCounter++;

            // ghost convective
            tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGGPConvective" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", (uint)fem::Element_Type::DOUBLE_SIDESET );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFluid" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "follower_phase_name", "PhaseFluid" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", (uint)fem::IWG_Type::GHOST_NORMAL_FIELD );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "VX,VY,VZ" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SPGPVelocity,GhostSP" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "ghost_order", (uint)tDispOrder );
            tIWGCounter++;

            // ghost pressure
            tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGGPPressure" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", (uint)fem::Element_Type::DOUBLE_SIDESET );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFluid" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "follower_phase_name", "PhaseFluid" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", (uint)fem::IWG_Type::GHOST_NORMAL_FIELD );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "P" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SPGPPressure,GhostSP" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "ghost_order", (uint)tDispOrder );
            tIWGCounter++;

            // ghost fluid temperature
            tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGGPFluidTemp" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", (uint)fem::Element_Type::DOUBLE_SIDESET );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFluid" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "follower_phase_name", "PhaseFluid" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", (uint)fem::IWG_Type::GHOST_NORMAL_FIELD );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "TEMP" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SPGPFluidTemp,GhostSP" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "ghost_order", (uint)tDispOrder );
            tIWGCounter++;

            // ghost solid A temperature
            tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGGPSolidTemp" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", (uint)fem::Element_Type::DOUBLE_SIDESET );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseSolid" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "follower_phase_name", "PhaseSolid" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", (uint)fem::IWG_Type::GHOST_NORMAL_FIELD );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "TEMP" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SPGPSolidTemp,GhostSP" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "ghost_order", (uint)tDispOrder );
            tIWGCounter++;
        }

        //------------------------------------------------------------------------------
        // fill the IQI part of the parameter list

        // init IQI counter
        uint tIQICounter = 0;

        // velocity VX
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIBulkVX" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", (uint)fem::IQI_Type::DOF );
        tParameterList( tIQIIndex )( tIQICounter ).set( "dof_quantity", "VX,VY,VZ" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index", 0 );
        tIQICounter++;

        // velocity VY
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIBulkVY" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", (uint)fem::IQI_Type::DOF );
        tParameterList( tIQIIndex )( tIQICounter ).set( "dof_quantity", "VX,VY,VZ" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index", 1 );
        tIQICounter++;

        // velocity VZ
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIBulkVZ" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", (uint)fem::IQI_Type::DOF );
        tParameterList( tIQIIndex )( tIQICounter ).set( "dof_quantity", "VX,VY,VZ" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index", 2 );
        tIQICounter++;

        // pressure
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIBulkP" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", (uint)fem::IQI_Type::DOF );
        tParameterList( tIQIIndex )( tIQICounter ).set( "dof_quantity", "P" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index", 0 );
        tIQICounter++;

        // temperature
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIBulkTEMP" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseAll" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", (uint)fem::IQI_Type::DOF );
        tParameterList( tIQIIndex )( tIQICounter ).set( "dof_quantity", "TEMP" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index", 0 );
        tIQICounter++;

        // fluid thermal energy on inlet
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIInletThermalEnergy" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_bulk_type", (uint)fem::Element_Type::SIDESET );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "neighbor_phases", "PhaseVoidFront" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", (uint)fem::IQI_Type::THERMAL_ENERGY_CONVECTIVE_FLUX );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_constitutive_models", "CMFluid,Fluid" );
        tIQICounter++;

        // fluid thermal energy on outlet
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIOutletThermalEnergy" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_bulk_type", (uint)fem::Element_Type::SIDESET );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", (uint)fem::IQI_Type::THERMAL_ENERGY_CONVECTIVE_FLUX );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "neighbor_phases", "PhaseVoidBack" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_constitutive_models", "CMFluid,Fluid" );
        tIQICounter++;

        // fluid total pressure on inlet
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIInletTotalPressure" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_bulk_type", (uint)fem::Element_Type::SIDESET );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "neighbor_phases", "PhaseVoidFront" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", (uint)fem::IQI_Type::TOTAL_PRESSURE );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_constitutive_models", "CMFluid,Fluid" );
        tIQICounter++;

        // fluid total pressure on outlet
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIOutletTotalPressure" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_bulk_type", (uint)fem::Element_Type::SIDESET );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "neighbor_phases", "PhaseVoidBack" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", (uint)fem::IQI_Type::TOTAL_PRESSURE );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_constitutive_models", "CMFluid,Fluid" );
        tIQICounter++;

        // fluid mass flow on inlet
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIInletMassFlow" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", (uint)fem::IQI_Type::MASS_FLOW );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "neighbor_phases", "PhaseVoidFront" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_bulk_type", (uint)fem::Element_Type::SIDESET );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_constitutive_models", "CMFluid,Fluid" );
        tIQICounter++;

        // fluid mass flow on outlet
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIOutletMassFlow" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "neighbor_phases", "PhaseVoidBack" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_bulk_type", (uint)fem::Element_Type::SIDESET );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", (uint)fem::IQI_Type::MASS_FLOW );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_constitutive_models", "CMFluid,Fluid" );
        tIQICounter++;

        // inclusion perimeter
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIPerimeterItf" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseSolid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "neighbor_phases", "PhaseFluid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_bulk_type", (uint)fem::Element_Type::SIDESET );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", (uint)fem::IQI_Type::VOLUME );
        tIQICounter++;

        // max temperature in solid
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIMaxTemp" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseSolid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", (uint)fem::IQI_Type::MAX_DOF );
        tParameterList( tIQIIndex )( tIQICounter ).set( "function_parameters", tRefTemp + "/" + tMaxTempExponent );
        tParameterList( tIQIIndex )( tIQICounter ).set( "dof_quantity", "TEMP" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index", 0 );
        tIQICounter++;

        // volume solid
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQISolidVolume" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseSolid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", (uint)fem::IQI_Type::VOLUME );
        tIQICounter++;

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        tParameterList( tFEMIndex ).push_back( prm::create_computation_parameter_list() );
        tParameterList( tFEMIndex )( 0 ).set( "print_physics_model", false );

        tParameterList( tFEMIndex )( 0 ).set( "finite_difference_scheme", tFEMFdScheme );
        tParameterList( tFEMIndex )( 0 ).set( "finite_difference_perturbation_size", tFEMFdEpsilon );
    }

    void
    SOLParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 8 );
        for ( uint Ik = 0; Ik < 8; Ik++ )
        {
            tParameterlist( Ik ).resize( 1 );
        }

        tParameterlist( 0 ).resize( 2 );

        tParameterlist( 0 )( 0 ) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL );
        tParameterlist( 0 )( 1 ) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL );

#ifdef MORIS_USE_MUMPS
        tParameterlist( 0 )( 0 ).set( "Solver_Type", "Amesos_Mumps" );
        tParameterlist( 0 )( 1 ).set( "Solver_Type", "Amesos_Mumps" );
#else
        tParameterlist( 0 )( 0 ).set( "Solver_Type", "Amesos_Superludist" );
        tParameterlist( 0 )( 1 ).set( "Solver_Type", "Amesos_Superludist" );
#endif

        /*
        tParameterlist( 0 )( 0 ) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::BELOS_IMPL );
        tParameterlist( 0 )( 0 ).set( "ifpack_prec_type", "ILU");
        tParameterlist( 0 )( 0 ).set( "fact: level-of-fill", 1);

        tParameterlist( 0 )( 1 ) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::BELOS_IMPL );
        tParameterlist( 0 )( 1 ).set( "ifpack_prec_type", "ILU");
        tParameterlist( 0 )( 1 ).set( "fact: level-of-fill", 3);
        */

        //------------------------------------------------------------------------------

        tParameterlist( 1 ).resize( 2 );

        tParameterlist( 1 )( 0 ) = moris::prm::create_linear_solver_parameter_list();
        tParameterlist( 1 )( 0 ).set( "DLA_Linear_solver_algorithms", "0" );

        tParameterlist( 1 )( 1 ) = moris::prm::create_linear_solver_parameter_list();
        tParameterlist( 1 )( 1 ).set( "DLA_Linear_solver_algorithms", "1" );

        //------------------------------------------------------------------------------

        tParameterlist( 2 ).resize( 4 );

        tParameterlist( 2 )( 0 ) = moris::prm::create_nonlinear_algorithm_parameter_list();    // nonlinear algorithm index 0
        tParameterlist( 2 )( 0 ).set( "NLA_Solver_Implementation", (uint)moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        tParameterlist( 2 )( 0 ).set( "NLA_Linear_solver", 0 );
        tParameterlist( 2 )( 0 ).set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        tParameterlist( 2 )( 0 ).set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
        tParameterlist( 2 )( 0 ).set( "NLA_max_iter", tNLA_max_iter );

        tParameterlist( 2 )( 1 ) = moris::prm::create_nonlinear_algorithm_parameter_list();
        tParameterlist( 2 )( 1 ).set( "NLA_Solver_Implementation", (uint)moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        tParameterlist( 2 )( 1 ).set( "NLA_Linear_solver", 0 );
        tParameterlist( 2 )( 1 ).set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        tParameterlist( 2 )( 1 ).set( "NLA_relaxation_parameter", 1.0 );
        tParameterlist( 2 )( 1 ).set( "NLA_max_iter", 1 );

        tParameterlist( 2 )( 2 ) = moris::prm::create_nonlinear_algorithm_parameter_list();    // nonlinear algorithm index 1
        tParameterlist( 2 )( 2 ).set( "NLA_Solver_Implementation", (uint)moris::NLA::NonlinearSolverType::NLBGS_SOLVER );
        tParameterlist( 2 )( 2 ).set( "NLA_Linear_solver", 0 );
        tParameterlist( 2 )( 2 ).set( "NLA_rel_res_norm_drop", 1.0 );
        tParameterlist( 2 )( 2 ).set( "NLA_max_iter", 1 );

        tParameterlist( 2 )( 3 ) = moris::prm::create_nonlinear_algorithm_parameter_list();
        tParameterlist( 2 )( 3 ).set( "NLA_Solver_Implementation", (uint)moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        tParameterlist( 2 )( 3 ).set( "NLA_Linear_solver", 1 );
        tParameterlist( 2 )( 3 ).set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        tParameterlist( 2 )( 3 ).set( "NLA_relaxation_parameter", 1.0 );
        tParameterlist( 2 )( 3 ).set( "NLA_max_iter", 1 );

        //------------------------------------------------------------------------------

        tParameterlist( 3 ).resize( 4 );

        tParameterlist( 3 )( 0 ) = moris::prm::create_nonlinear_solver_parameter_list();    // 0: fluid subproblem
        tParameterlist( 3 )( 0 ).set( "NLA_Nonlinear_solver_algorithms", "0" );             // set nonlinear algorithm with index 0
        tParameterlist( 3 )( 0 ).set( "NLA_Solver_Implementation", (uint)moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        tParameterlist( 3 )( 0 ).set( "NLA_DofTypes", "VX,VY,VZ,P" );

        tParameterlist( 3 )( 1 ) = moris::prm::create_nonlinear_solver_parameter_list();    // 1: thermal subproblem
        tParameterlist( 3 )( 1 ).set( "NLA_Nonlinear_solver_algorithms", "1" );             // set nonlinear algorithm with index 0
        tParameterlist( 3 )( 1 ).set( "NLA_Solver_Implementation", (uint)moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        tParameterlist( 3 )( 1 ).set( "NLA_DofTypes", "TEMP" );

        tParameterlist( 3 )( 2 ) = moris::prm::create_nonlinear_solver_parameter_list();    // 2: one-way coupling via NLBGS
        tParameterlist( 3 )( 2 ).set( "NLA_Nonlinear_solver_algorithms", "2" );             // set nonlinear algorithm with index 1.
        tParameterlist( 3 )( 2 ).set( "NLA_Solver_Implementation", (uint)moris::NLA::NonlinearSolverType::NLBGS_SOLVER );
        tParameterlist( 3 )( 2 ).set( "NLA_Sub_Nonlinear_Solver", "0,1" );                  // set sub nonlinear solvers with index 0 and 1
        tParameterlist( 3 )( 2 ).set( "NLA_DofTypes", "VX,VY,VZ,P;TEMP" );

        tParameterlist( 3 )( 3 ) = moris::prm::create_nonlinear_solver_parameter_list();
        tParameterlist( 3 )( 3 ).set( "NLA_DofTypes", "VX,VY,VZ,P,TEMP" );
        tParameterlist( 3 )( 3 ).set( "NLA_Nonlinear_solver_algorithms", "3" );

        // ----------------------------------------------------------

        tParameterlist( 4 )( 0 ) = moris::prm::create_time_solver_algorithm_parameter_list();
        tParameterlist( 4 )( 0 ).set( "TSA_Nonlinear_solver", 2 );                      // using NLBGS for forward problem
        tParameterlist( 4 )( 0 ).set( "TSA_nonlinear_solver_for_adjoint_solve", 3 );    // using monlithic for sensitivity problem

        //------------------------------------------------------------------------------

        tParameterlist( 5 )( 0 ) = moris::prm::create_time_solver_parameter_list();
        tParameterlist( 5 )( 0 ).set( "TSA_DofTypes", "VX,VY,VZ;P;TEMP" );
        tParameterlist( 5 )( 0 ).set( "TSA_Initialize_Sol_Vec", "VX,0.01;VY,0.0;VZ,0.0;P,0.0;TEMP,0.0" );
        tParameterlist( 5 )( 0 ).set( "TSA_Output_Indices", "0" );
        tParameterlist( 5 )( 0 ).set( "TSA_Output_Criteria", "Output_Criterion" );

        //------------------------------------------------------------------------------

        tParameterlist( 6 )( 0 ) = moris::prm::create_solver_warehouse_parameterlist();

        tParameterlist( 7 )( 0 ) = moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::NONE );
    }

    void
    MSIParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_msi_parameter_list();
        tParameterlist( 0 )( 0 ).set( "VX", 1 );
        tParameterlist( 0 )( 0 ).set( "VY", 1 );
        tParameterlist( 0 )( 0 ).set( "VZ", 1 );
        tParameterlist( 0 )( 0 ).set( "TEMP", 1 );
    }

    void
    VISParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_vis_parameter_list();
        tParameterlist( 0 )( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", tExoFile ) );
        tParameterlist( 0 )( 0 ).set( "Mesh_Type", (uint)vis::VIS_Mesh_Type::STANDARD );
        tParameterlist( 0 )( 0 ).set( "Set_Names", tAllPhases );
        tParameterlist( 0 )( 0 ).set( "Field_Names", "VX,VY,VZ,P,TEMP" );
        tParameterlist( 0 )( 0 ).set( "Field_Type", "NODAL,NODAL,NODAL,NODAL,NODAL" );
        tParameterlist( 0 )( 0 ).set( "IQI_Names", "IQIBulkVX,IQIBulkVY,IQIBulkVZ,IQIBulkP,IQIBulkTEMP" );
        tParameterlist( 0 )( 0 ).set( "Save_Frequency", 1 );
        tParameterlist( 0 )( 0 ).set( "Time_Offset", 10.0 );
    }

    void
    MORISGENERALParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
    }

    //------------------------------------------------------------------------------
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
