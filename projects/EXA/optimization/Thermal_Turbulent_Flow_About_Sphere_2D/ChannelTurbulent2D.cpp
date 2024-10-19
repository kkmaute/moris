/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * ChannelTurbulent2D.cpp
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

// global variable for interpolation order
extern uint gInterpolationOrder;

// problem dimension: 2D or 3D
extern uint gTestCaseIndex;

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

    /* ------------------------------------------------------------------------ */
    // set optimization restart iteration (if 0 no restart)
    sint tRestartId = 0;

    /* ------------------------------------------------------------------------ */

    // inlet fluid and thermal properties
    real tInPres        = 0.001 * 101325.0;    // N/m2  - 0.6 atm
    real tInTemp        = 0.0;                 // K     - 500 K
    real tTurbIntensity = 3.0;

    // reference values
    real presref = tInPres;     // reference pressure
    real lenref  = 1e-3;        // reference length
    real rhoref  = 0.3289;      // reference density
    real dynvisc = 43.32e-6;    // reference kinematic visc
    real tempref = 1.0;         // reference temperature

    // derived reference values
    real velref  = std::sqrt( presref / rhoref );
    real timeref = lenref / velref;
    real massref = rhoref * std::pow( lenref, 3.0 );

    // fluid properties
    real capfld  = 1.210e3;    // J/kg K
    real prandtl = 0.7;
    real turbPr  = 0.9;

    real kinvisc = dynvisc / rhoref;
    real confld  = ( dynvisc * capfld ) / prandtl;

    // solid properties
    real rhosolid = 7.0e3;
    real consolid = 150.0;
    real capsolid = 200.0;

    // reynolds number based on reference values and fluid properties see below
    real reynolds = lenref * velref / kinvisc;

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

    real tSpherePosX = 7.50e-4 * tLengthScale;
    real tSpherePosY = 5.00e-4 * tLengthScale;

    real tSphereRadius = 2.25e-4 * tLengthScale;
    real tSphereExpon  = 2.0;

    /* ------------------------------------------------------------------------ */
    // basic geometry parameters
    real tDimX = 5.25e-3 * tLengthScale;
    real tDimY = 1.25e-3 * tLengthScale;

    real tOffsetX = ( tChannelLength - tDimX ) / 2.0;
    real tOffsetY = ( tChannelHeight - tDimY ) / 2.0;

    real tApproxEleSize = 1e-4 * tLengthScale;

    /* ------------------------------------------------------------------------ */
    // background mesh
    std::string tNumElemX = moris_to_string( std::ceil( tDimX / tApproxEleSize ) );
    std::string tNumElemY = moris_to_string( std::ceil( tDimY / tApproxEleSize ) );

    std::string tDomainDimX = moris_to_string( tDimX );
    std::string tDomainDimY = moris_to_string( tDimY );

    std::string tDomainOffX = moris_to_string( tOffsetX );
    std::string tDomainOffY = moris_to_string( tOffsetY );

    std::string tNumElemsPerDim = tNumElemX + "," + tNumElemY;
    std::string tDomainDims     = tDomainDimX + "," + tDomainDimY;
    std::string tDomainOffset   = tDomainOffX + "," + tDomainOffY;
    std::string tDomainSidesets = "1,2,3,4";

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
    real tHeatFlx = 250.0e4;    // W/m^2 - 250 W/cm^2

    std::string tInletPressure     = moris_to_string( tInPres * tPressureScale );
    std::string tInletTemperature  = moris_to_string( tInTemp * tTempScale );
    std::string tInletKinViscosity = moris_to_string( tTurbIntensity * kinvisc * tLengthScale * tLengthScale / tTimeScale );    // m2/s air at 800 C);

    std::string tVolumetricHeatLoad = moris_to_string( tHeatFlx * tPowerScale / tLengthScale / tLengthScale / tLengthScale );

    /* ------------------------------------------------------------------------ */
    // material parameters
    std::string tFluidDensity           = moris_to_string( rhoref * tDensityScale );                                // kg/m3   air at 800 C and 1 atm
    std::string tFluidDynViscosity      = moris_to_string( dynvisc * tPressureScale * tTimeScale );                 // N s/m2  air at 800 C
    std::string tFluidKinViscosity      = moris_to_string( kinvisc * tLengthScale * tLengthScale / tTimeScale );    // m2/s air at 800 C
    std::string tFluidCapacity          = moris_to_string( capfld * tEnergyScale / tMassScale / tTempScale );       // J/kg K  cp air at 800 C and 1 atm
    std::string tFluidConductivity      = moris_to_string( confld * tPowerScale / tLengthScale / tTempScale );      // W/m K   air at 800 C
    std::string tFluidTurbulencePrandtl = moris_to_string( prandtl );
    std::string tFluidPressureSpring    = "0.0";

    std::string tSolidDensity      = moris_to_string( rhosolid * tDensityScale );                              // kg/m3
    std::string tSolidCapacity     = moris_to_string( capsolid * tEnergyScale / tMassScale / tTempScale );     // J/kg*K
    std::string tSolidConductivity = moris_to_string( consolid * tPowerScale / tLengthScale / tTempScale );    // W/m K ;

    /* ------------------------------------------------------------------------ */
    // turbulence model parameters
    std::string tCMTurbFt2      = "1.0";
    std::string tCMTurbAlpha    = "10.0";
    std::string tSupgTurbPower  = "2.0";
    std::string tSupgTurbSource = "1.0";
    std::string tSupgFluidC1    = "36.0";

    /* ------------------------------------------------------------------------ */
    // material parameters for wall distance

    // conductivity
    std::string tConductivityThetaPhi = "1.0";

    // capacity
    std::string tCapacityTheta = "0.1";
    std::string tCapacityPhi   = "0.0";

    // density
    std::string tDensityTheta = "1.0";
    std::string tDensityPhi   = "0.0";

    // prescribed theta on interface
    std::string tPrescTheta = "1.0";

    // prescribed phi on interface
    std::string tPrescPhi = "0.0";

    /* ------------------------------------------------------------------------ */
    // File names
    std::string tName          = "ChannelTurbulent2D";
    std::string tProblemConfig = "_" + std::to_string( gTestCaseIndex );

    std::string tSoFile        = tName + ".so";
    std::string tExoFile       = tName + tProblemConfig + ".exo";
    std::string tHdf5File      = tName + tProblemConfig + "_SEN.hdf5";
    std::string tGENOutputFile = tName + tProblemConfig + "_GEN.exo";

    // FD in adjoint
    real tFEMFdEpsilon = 1.0e-4;

    // FD in sweep
    std::string tSweepFdEpsilon = "1.0e-4";

    // Element inverse estimate
    std::string tCInv = gInterpolationOrder == 2 ? "60.0" : "36.0";

    /* ------------------------------------------------------------------------ */

    // number of constraints (here: number of design criteria)
    moris::uint tNumConstraints = 11;

    // max dof IQI
    std::string tRefTemp         = moris_to_string( 250.0 * tTempScale );
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

    moris::real
    Func_Sphere(
            const moris::Matrix< DDRMat >& aCoordinates,
            const Vector< real >&          aGeometryParameters )
    {
        // get coordinates
        real tX = aCoordinates( 0 );
        real tY = aCoordinates( 1 );

        real tRadius = aGeometryParameters( 0 );

        real tReturnValue = tRadius - std::pow( std::pow( tX - tSpherePosX, tSphereExpon ) + std::pow( tY - tSpherePosY, tSphereExpon ), 1.0 / tSphereExpon );

        return tReturnValue;
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

        // get normal
        real tNx = aGeometryParameters( 0 );
        real tNy = aGeometryParameters( 1 );

        // get point on plane
        real tPx = aGeometryParameters( 2 );
        real tPy = aGeometryParameters( 3 );

        real tReturnValue = tNx * ( tPx - tX ) + tNy * ( tPy - tY );

        return tReturnValue;
    }

    /* ------------------------------------------------------------------------ */
    // property functions

    // inlet viscosity function
    void
    Func_Inlet_V(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        aPropMatrix.set_size( 1, 1, 0.0 );

        real        tA      = 0.05;
        moris::real tInVisc = aParameters( 0 )( 0 );
        real        tY      = M_PI * ( aFIManager->get_IP_geometry_interpolator()->valx()( 1 ) - tGeoShift ) / ( tChannelHeight );
        aPropMatrix( 0 )    = tInVisc * std::sin( tY ) / std::sqrt( tA * tA + std::sin( tY ) * std::sin( tY ) );
    }

    // wall distance function
    void
    Func_Wall_Distance( moris::Matrix< moris::DDRMat >& aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >&   aParameters,
            moris::fem::Field_Interpolator_Manager*     aFIManager )
    {
        aPropMatrix = aFIManager->get_field_interpolators_for_type( MSI::Dof_Type::PHID )->val();
    }

    // Wall distance derivative function
    void
    Func_Wall_Distance_Der( moris::Matrix< moris::DDRMat >& aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >&       aParameters,
            moris::fem::Field_Interpolator_Manager*         aFIManager )
    {
        aPropMatrix = aFIManager->get_field_interpolators_for_type( MSI::Dof_Type::PHID )->N();
    }

    /* ------------------------------------------------------------------------ */
    // Phase assignement

    uint
    get_phase_index( const Bitset< 5 >& aGeometrySigns )
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
                || aGeometrySigns.test( 4 ) )
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
    0- IQIOutletThermalEnergy
    1- IQIInletTotalPressure
    2- IQIOutletTotalPressure
    3- IQIPerimeterIfc
    4- IQIInletMassFlow
    5- IQIOutletMassFlow
    6- IQIMaxTemp
    7- IQISolidVolume
    8- IQIInletPowDisp
    9 - IQIOutletPowDisp
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

        // constraints
        tConstraints( 0 )  = aCriteria( 0 );
        tConstraints( 1 )  = aCriteria( 1 );
        tConstraints( 2 )  = aCriteria( 2 );
        tConstraints( 3 )  = aCriteria( 3 );
        tConstraints( 4 )  = aCriteria( 4 );
        tConstraints( 5 )  = aCriteria( 5 );
        tConstraints( 6 )  = aCriteria( 6 );
        tConstraints( 7 )  = aCriteria( 7 );
        tConstraints( 8 )  = aCriteria( 8 );
        tConstraints( 9 )  = aCriteria( 9 );
        tConstraints( 10 ) = aCriteria( 10 );

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

        tDConstraintDCriteria( 0, 0 )   = 1.0;
        tDConstraintDCriteria( 1, 1 )   = 1.0;
        tDConstraintDCriteria( 2, 2 )   = 1.0;
        tDConstraintDCriteria( 3, 3 )   = 1.0;
        tDConstraintDCriteria( 4, 4 )   = 1.0;
        tDConstraintDCriteria( 5, 5 )   = 1.0;
        tDConstraintDCriteria( 6, 6 )   = 1.0;
        tDConstraintDCriteria( 7, 7 )   = 1.0;
        tDConstraintDCriteria( 8, 8 )   = 1.0;
        tDConstraintDCriteria( 9, 9 )   = 1.0;
        tDConstraintDCriteria( 10, 10 ) = 1.0;

        return tDConstraintDCriteria;
    }

    /* ------------------------------------------------------------------------ */

    void
    OPTParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists.set( "is_optimization_problem", true );
        aParameterLists.set( "problem", "user_defined" );
        aParameterLists.set( "library", tSoFile );

        aParameterLists( OPT::ALGORITHMS ).add_parameter_list( opt::Optimization_Algorithm_Type::SWEEP );
        aParameterLists.set( "hdf5_path", tHdf5File );
        aParameterLists.set( "num_evaluations_per_adv", "1" );
        aParameterLists.set( "finite_difference_type", "all" );
        aParameterLists.set( "finite_difference_epsilons", tSweepFdEpsilon );
    }

    /* ------------------------------------------------------------------------ */

    void
    HMRParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists.set( "number_of_elements_per_dimension", tNumElemsPerDim );
        aParameterLists.set( "domain_dimensions", tDomainDims );
        aParameterLists.set( "domain_offset", tDomainOffset );
        aParameterLists.set( "domain_sidesets", tDomainSidesets );

        aParameterLists.set( "lagrange_output_meshes", "0" );

        aParameterLists.set( "lagrange_orders", tLagrangeOrder );
        aParameterLists.set( "lagrange_pattern", tLagrangePattern );

        aParameterLists.set( "bspline_orders", tBsplineOrder );
        aParameterLists.set( "bspline_pattern", "0,1" );

        aParameterLists.set( "lagrange_to_bspline", "0,1" );

        aParameterLists.set( "truncate_bsplines", 1 );
        aParameterLists.set( "refinement_buffer", tRefineBuffer );
        aParameterLists.set( "staircase_buffer", tRefineBuffer );

        aParameterLists.set( "initial_refinement", tInitialRef );
        aParameterLists.set( "initial_refinement_pattern", "0,1" );

        aParameterLists.set( "use_number_aura", 1 );

        aParameterLists.set( "use_multigrid", 0 );
        aParameterLists.set( "severity_level", 0 );
    }

    /* ------------------------------------------------------------------------ */

    void
    XTKParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "decompose", true );
        aParameterLists.set( "decomposition_type", "conformal" );
        aParameterLists.set( "enrich", true );
        aParameterLists.set( "basis_rank", "bspline" );
        aParameterLists.set( "enrich_mesh_indices", "0,1" );
        aParameterLists.set( "ghost_stab", tUseGhost );
        aParameterLists.set( "multigrid", false );
        aParameterLists.set( "verbose", false );
        aParameterLists.set( "high_to_low_dbl_side_sets", true );
        aParameterLists.set( "print_enriched_ig_mesh", false );
        aParameterLists.set( "exodus_output_XTK_ig_mesh", true );
    }

    /* ------------------------------------------------------------------------ */

    void
    GENParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists.set( "IQI_types",
                "IQIOutletThermalEnergy",
                "IQIInletTotalPressure",
                "IQIOutletTotalPressure",
                "IQIPerimeterItf",
                "IQIInletMassFlow",
                "IQIOutletMassFlow",
                "IQIMaxTemp",
                "IQISolidVolume",
                "IQIInletPowDisp",
                "IQIOutletPowDisp",
                "IQIVolumePowDisp" );

        aParameterLists.set( "number_of_phases", 5 );
        aParameterLists.set( "phase_function_name", "get_phase_index" );
        aParameterLists.set( "output_mesh_file", tGENOutputFile );
        aParameterLists.set( "time_offset", 10.0 );

        // Inclusions
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( gen::Field_Type::USER_DEFINED );
        aParameterLists.set( "field_function_name", "Func_Sphere" );
        aParameterLists.set( "sensitivity_function_name", "Func_Sphere_Deriv" );
        aParameterLists( 1 ).insert( "radius", Design_Variable( tSphereRadius * 0.9, tSphereRadius, tSphereRadius * 1.1 ) );
        aParameterLists.set( "number_of_refinements", tInterfaceRefinementSphere );
        aParameterLists.set( "refinement_mesh_index", 0 );
        aParameterLists.set( "isocontour_threshold", 0.0 );
        aParameterLists.set( "isocontour_tolerance", 1.0e-12 );
        aParameterLists.set( "intersection_tolerance", 1.0e-12 );

        // Inlet plane
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( gen::Field_Type::USER_DEFINED );
        aParameterLists.set( "field_function_name", "Func_Plane" );
        aParameterLists( 1 ).insert( "variable_1", 1.0 );
        aParameterLists( 1 ).insert( "variable_2", 0.0 );
        aParameterLists( 1 ).insert( "variable_3", tGeoShift );
        aParameterLists( 1 ).insert( "variable_4", 0.0 );
        aParameterLists.set( "number_of_refinements", tInterfaceRefinementWalls );
        aParameterLists.set( "refinement_mesh_index", 0 );
        aParameterLists.set( "isocontour_threshold", 0.0 );
        aParameterLists.set( "isocontour_tolerance", 1.0e-12 );
        aParameterLists.set( "intersection_tolerance", 1.0e-12 );

        // Out plane
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( gen::Field_Type::USER_DEFINED );
        aParameterLists.set( "field_function_name", "Func_Plane" );
        aParameterLists( 1 ).insert( "variable_1", -1.0 );
        aParameterLists( 1 ).insert( "variable_2", 0.0 );
        aParameterLists( 1 ).insert( "variable_3", tChannelLength + tGeoShift );
        aParameterLists( 1 ).insert( "variable_4", 0.0 );
        aParameterLists.set( "number_of_refinements", tInterfaceRefinementWalls );
        aParameterLists.set( "refinement_mesh_index", 0 );
        aParameterLists.set( "isocontour_threshold", 0.0 );
        aParameterLists.set( "isocontour_tolerance", 1.0e-12 );
        aParameterLists.set( "intersection_tolerance", 1.0e-12 );

        // Lower plane
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( gen::Field_Type::USER_DEFINED );
        aParameterLists.set( "field_function_name", "Func_Plane" );
        aParameterLists( 1 ).insert( "variable_1", 0.0 );
        aParameterLists( 1 ).insert( "variable_2", 1.0 );
        aParameterLists( 1 ).insert( "variable_3", 0.0 );
        aParameterLists( 1 ).insert( "variable_4", tGeoShift );
        aParameterLists.set( "number_of_refinements", tInterfaceRefinementWalls );
        aParameterLists.set( "refinement_mesh_index", 0 );
        aParameterLists.set( "isocontour_threshold", 0.0 );
        aParameterLists.set( "isocontour_tolerance", 1.0e-12 );
        aParameterLists.set( "intersection_tolerance", 1.0e-12 );

        // Upper plane
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( gen::Field_Type::USER_DEFINED );
        aParameterLists.set( "field_function_name", "Func_Plane" );
        aParameterLists( 1 ).insert( "variable_1", 0.0 );
        aParameterLists( 1 ).insert( "variable_2", -1.0 );
        aParameterLists( 1 ).insert( "variable_3", 0.0 );
        aParameterLists( 1 ).insert( "variable_4", tChannelHeight + tGeoShift );
        aParameterLists.set( "number_of_refinements", tInterfaceRefinementWalls );
        aParameterLists.set( "refinement_mesh_index", 0 );
        aParameterLists.set( "isocontour_threshold", 0.0 );
        aParameterLists.set( "isocontour_tolerance", 1.0e-12 );
        aParameterLists.set( "intersection_tolerance", 1.0e-12 );
    }

    /* ------------------------------------------------------------------------ */

    void
    FEMParameterList( Module_Parameter_Lists& aParameterLists )
    {
        if ( par_rank() == 0 )
        {
            std::cout << "presref        = " << presref << '\n';
            std::cout << "lenref         = " << lenref << '\n';
            std::cout << "rhoref         = " << rhoref << '\n';
            std::cout << "tempref        = " << tempref << '\n';

            std::cout << '\n';

            std::cout << "velref         = " << velref << '\n';
            std::cout << "timeref        = " << timeref << '\n';
            std::cout << "massref        = " << massref << '\n';

            std::cout << '\n';

            std::cout << "tLengthScale   = " << tLengthScale << '\n';
            std::cout << "tTimeScale     = " << tTimeScale << '\n';
            std::cout << "tMassScale     = " << tMassScale << '\n';
            std::cout << "tTempScale     = " << tTempScale << '\n';

            std::cout << '\n';

            std::cout << "tPressureScale = " << tPressureScale << '\n';
            std::cout << "tEnergyScale   = " << tEnergyScale << '\n';
            std::cout << "tPowerScale    = " << tPowerScale << '\n';
            std::cout << "tDensityScale  = " << tDensityScale << '\n';

            std::cout << '\n';

            std::cout << "tInletPressure       = " << tInletPressure << '\n';
            std::cout << "tInletTemperature    = " << tInletTemperature << '\n';
            std::cout << "tVolumetricHeatLoad  = " << tVolumetricHeatLoad << '\n';
            std::cout << "tInletKinViscosity   = " << tInletKinViscosity << '\n';

            std::cout << '\n';

            std::cout << "tFluidDensity        = " << tFluidDensity << '\n';
            std::cout << "tFluidDynViscosity   = " << tFluidDynViscosity << '\n';
            std::cout << "tFluidKinViscosity   = " << tFluidKinViscosity << '\n';
            std::cout << "tFluidCapacity       = " << tFluidCapacity << '\n';
            std::cout << "tFluidConductivity   = " << tFluidConductivity << '\n';
            std::cout << "tFluidPressureSpring = " << tFluidPressureSpring << '\n';

            std::cout << '\n';

            std::cout << "tSolidDensity        = " << tSolidDensity << '\n';
            std::cout << "tSolidCapacity       = " << tSolidCapacity << '\n';
            std::cout << "tSolidConductivity   = " << tSolidConductivity << '\n';

            std::cout << '\n';

            std::cout << "Reynolds number     = " << 1.0 / std::stod( tFluidDynViscosity ) << " (" << reynolds << ")\n";
        }

        //------------------------------------------------------------------------------

        aParameterLists( FEM::PHASES ).add_parameter_list();
        aParameterLists.set( "phase_name", "PhaseFluid" );
        aParameterLists.set( "phase_indices", "0" );

        aParameterLists( FEM::PHASES ).add_parameter_list();
        aParameterLists.set( "phase_name", "PhaseSolid" );
        aParameterLists.set( "phase_indices", "1" );

        aParameterLists( FEM::PHASES ).add_parameter_list();
        aParameterLists.set( "phase_name", "PhaseVoidFront" );
        aParameterLists.set( "phase_indices", "2" );

        aParameterLists( FEM::PHASES ).add_parameter_list();
        aParameterLists.set( "phase_name", "PhaseVoidBack" );
        aParameterLists.set( "phase_indices", "3" );

        aParameterLists( FEM::PHASES ).add_parameter_list();
        aParameterLists.set( "phase_name", "PhaseVoidLateral" );
        aParameterLists.set( "phase_indices", "4" );

        aParameterLists( FEM::PHASES ).add_parameter_list();
        aParameterLists.set( "phase_name", "PhaseAll" );
        aParameterLists.set( "phase_indices", "0,1" );

        //------------------------------------------------------------------------------
        // fill the property part of the parameter list

        // fluid properties ------------------------------------------------------------
        // create fluid viscosity property
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropFluidDynViscosity" );
        aParameterLists.set( "function_parameters", tFluidDynViscosity );

        // create fluid density property
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropFluidDensity" );
        aParameterLists.set( "function_parameters", tFluidDensity );

        // create fluid kinematic viscosity property
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropFluidKinViscosity" );
        aParameterLists.set( "function_parameters", tFluidKinViscosity );

        // create fluid capacity property
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropFluidCapacity" );
        aParameterLists.set( "function_parameters", tFluidCapacity );

        // create fluid conductivity property
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropFluidConductivity" );
        aParameterLists.set( "function_parameters", tFluidConductivity );

        // create turbulent prandtl number
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropFluidTurbPrandtl" );
        aParameterLists.set( "function_parameters", tFluidTurbulencePrandtl );

        // create fluid pressure spring property
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropFluidPressureSpring" );
        aParameterLists.set( "function_parameters", tFluidPressureSpring );

        // solid properties ----------------------------------------------------------
        // create solid B density property
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropSolidDensity" );
        aParameterLists.set( "function_parameters", tSolidDensity );

        // create solid B capacity property
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropSolidCapacity" );
        aParameterLists.set( "function_parameters", tSolidCapacity );

        // create solid B conductivity property
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropSolidConductivity" );
        aParameterLists.set( "function_parameters", tSolidConductivity );

        // BC properties ---------------------------------------------------------------
        // create inlet pressure property
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropInletPressure" );
        aParameterLists.set( "function_parameters", tInletPressure );

        // create inlet temperature property
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropInletTemp" );
        aParameterLists.set( "function_parameters", tInletTemperature );

        // create wall velocity property
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropZeroU" );
        aParameterLists.set( "function_parameters", "0.0;0.0" );

        // create symmetry velocity property (x-z plane)
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropSelectUY" );
        aParameterLists.set( "function_parameters", "0.0,0.0;0.0,1.0" );

        // create inlet kinematic viscosity property
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropInletV" );
        aParameterLists.set( "value_function", "Func_Inlet_V" );
        aParameterLists.set( "function_parameters", tInletKinViscosity );

        // create  wall kinematic viscosity property
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropZeroV" );
        aParameterLists.set( "function_parameters", "0.0" );

        // create heat load property
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropVolumetricHeatFlux" );
        aParameterLists.set( "function_parameters", tVolumetricHeatLoad );

        // Wall distance properties ----------------------------------------------------
        // create wall distance property
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropWallDistance" );
        aParameterLists.set( "dof_dependencies", "PHID" );
        aParameterLists.set( "value_function", "Func_Wall_Distance" );
        aParameterLists.set( "dof_derivative_functions", "Func_Wall_Distance_Der" );

        // create common conductivity property for theta and phi problems
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropConductivity" );
        aParameterLists.set( "function_parameters", tConductivityThetaPhi );

        // create density property for theta problem
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDensityTheta" );
        aParameterLists.set( "function_parameters", tDensityTheta );

        // create capacity property for theta problem
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropCapacityTheta" );
        aParameterLists.set( "function_parameters", tCapacityTheta );

        // create prescribed BC property for theta problem
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropPrescTheta" );
        aParameterLists.set( "function_parameters", tPrescTheta );

        // create density property for phi problem
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDensityPhi" );
        aParameterLists.set( "function_parameters", tDensityPhi );

        // create capacity property for phi problem
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropCapacityPhi" );
        aParameterLists.set( "function_parameters", tCapacityPhi );

        // create prescribed BC property for phi problem
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropPrescPhi" );
        aParameterLists.set( "function_parameters", tPrescPhi );

        // create eigenstrain BC property for phi problem
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropEigenStrainPhi" );
        aParameterLists.set( "function_parameters", "1.0" );

        // create time continuity weight current
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropWeightCurrent" );
        aParameterLists.set( "function_parameters", "10.0" );

        // create time continuity weight previous
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropWeightPrevious" );
        aParameterLists.set( "function_parameters", "10.0" );

        // create initial condition property
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropInitialCondition" );
        aParameterLists.set( "function_parameters", "0.0" );

        //------------------------------------------------------------------------------
        // fill the constitutive model part of the parameter list

        // fluid CM --------------------------------------------------------------------
        // create fluid CM
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMFluid" );
        aParameterLists.set( "phase_name", "PhaseFluid" );
        aParameterLists.set( "constitutive_type", fem::Constitutive_Type::FLUID_TURBULENCE );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "VX,VY;P;VISCOSITY", "Velocity,Pressure,Viscosity" ) );
        aParameterLists.set( "properties",
                "PropFluidDynViscosity,Viscosity;"
                "PropFluidKinViscosity,KinViscosity;"
                "PropFluidDensity  ,Density" );

        // create CM SA turbulence
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMTurbulence" );
        aParameterLists.set( "phase_name", "PhaseFluid" );
        aParameterLists.set( "constitutive_type", fem::Constitutive_Type::SPALART_ALLMARAS_TURBULENCE );
        aParameterLists.set( "function_parameters", tCMTurbFt2 + "/" + tCMTurbAlpha );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "VX,VY;VISCOSITY", "Velocity,Viscosity" ) );
        aParameterLists.set( "properties",
                "PropFluidKinViscosity,KinViscosity;"
                "PropWallDistance     ,WallDistance" );

        // create fluid diffusion CM
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMFluidDiffusion" );
        aParameterLists.set( "phase_name", "PhaseFluid" );
        aParameterLists.set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO_TURBULENCE );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists.set( "properties",
                "PropFluidConductivity,Conductivity;"
                "PropFluidDensity     ,Density;"
                "PropFluidCapacity    ,HeatCapacity;"
                "PropFluidKinViscosity,KinematicViscosity;"
                "PropFluidTurbPrandtl ,TurbulentPrandtl" );

        // solid CM --------------------------------------------------------------------
        // create solid B diffusion CM
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMSolidDiffusion" );
        aParameterLists.set( "phase_name", "PhaseSolid" );
        aParameterLists.set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists.set( "properties",
                "PropSolidConductivity,Conductivity;"
                "PropSolidDensity     ,Density;"
                "PropSolidCapacity    ,HeatCapacity" );

        // theta/phi CM --------------------------------------------------------------------
        // create fluid theta diffusion CM
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMFluidDiffusionTheta" );
        aParameterLists.set( "phase_name", "PhaseFluid" );
        aParameterLists.set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "THETA", "Temperature" ) );
        aParameterLists.set( "properties",
                "PropConductivity      , Conductivity;"
                "PropDensityTheta      , Density;"
                "PropCapacityTheta     , HeatCapacity" );

        // create solid theta diffusion CM
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMInclusionSolidDiffusionTheta" );
        aParameterLists.set( "phase_name", "PhaseSolid" );
        aParameterLists.set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "THETA", "Temperature" ) );
        aParameterLists.set( "properties",
                "PropConductivity      , Conductivity;"
                "PropDensityTheta      , Density;"
                "PropCapacityTheta     , HeatCapacity" );

        // create fluid phi diffusion CM
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMFluidDiffusionPhi" );
        aParameterLists.set( "phase_name", "PhaseFluid" );
        aParameterLists.set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "PHID;THETA", "Temperature,Theta" ) );
        aParameterLists.set( "properties",
                "PropConductivity    , Conductivity;"
                "PropDensityPhi      , Density;"
                "PropCapacityPhi     , HeatCapacity;"
                "PropEigenStrainPhi  , EigenStrain" );

        // create solid phi diffusion CM
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMInclusionSolidDiffusionPhi" );
        aParameterLists.set( "phase_name", "PhaseSolid" );
        aParameterLists.set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "PHID;THETA", "Temperature,Theta" ) );
        aParameterLists.set( "properties",
                "PropConductivity    , Conductivity;"
                "PropDensityPhi      , Density;"
                "PropCapacityPhi     , HeatCapacity;"
                "PropEigenStrainPhi  , EigenStrain" );

        //------------------------------------------------------------------------------
        // fill the stabilization parameter part of the parameter list

        // create SUPG/PSG fluid
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPIncFlow" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::INCOMPRESSIBLE_FLOW );
        aParameterLists.set( "function_parameters", tCInv + "/0.0" );
        aParameterLists.set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY;P", "Velocity,Pressure" ) );
        aParameterLists.set( "leader_properties",
                "PropFluidDynViscosity,Viscosity;"
                "PropFluidDensity  ,Density" );

        // create SUPG fluid temperature advection
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPSUPGTemp" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::SUPG_ADVECTION );
        aParameterLists.set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
        aParameterLists.set( "leader_properties", "PropFluidConductivity,Conductivity" );

        // create SUPG Spalart-Allmaras model
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPSUPGSA" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::SUPG_SPALART_ALLMARAS_TURBULENCE );
        aParameterLists.set( "leader_constitutive_models", "CMTurbulence,SpalartAllmarasTurbulence" );
        aParameterLists.set( "function_parameters", tSupgTurbPower + "/" + tSupgTurbSource );
        aParameterLists.set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY;VISCOSITY", "Velocity,Viscosity" ) );

        // create Nitsche for fluid velocity
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPNitscheU" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::VELOCITY_DIRICHLET_NITSCHE );
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
        aParameterLists.set( "leader_properties",
                "PropFluidDynViscosity,Viscosity;"
                "PropFluidDensity,Density" );

        // create Nitsche for fluid temperature
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPNitscheT" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "leader_properties", "PropFluidConductivity,Material" );

        // create Nitsche for fluid turbulent viscosity
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPNitscheV" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::TURBULENCE_DIRICHLET_NITSCHE );
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "leader_constitutive_models", "CMTurbulence,SpalartAllmarasTurbulence" );

        // create Nitsche for fluid/solid diffusion interface
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPFSInterfaceNitsche" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "follower_phase_name", "PhaseSolid" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::NITSCHE_INTERFACE );
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "leader_properties", "PropFluidConductivity,Material" );
        aParameterLists.set( "follower_properties", "PropSolidConductivity,Material" );

        // create ghost penalty viscous
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPGPViscous" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::VISCOUS_GHOST );
        aParameterLists.set( "function_parameters", "0.05" );
        aParameterLists.set( "leader_properties", "PropFluidDynViscosity,Viscosity" );

        // create ghost penalty convective
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPGPVelocity" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::CONVECTIVE_GHOST );
        aParameterLists.set( "function_parameters", "0.05" );
        aParameterLists.set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
        aParameterLists.set( "leader_properties", "PropFluidDensity,Density" );

        // create ghost penalty pressure
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPGPPressure" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::PRESSURE_GHOST );
        aParameterLists.set( "function_parameters", "0.005" );
        aParameterLists.set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
        aParameterLists.set( "leader_properties",
                "PropFluidDynViscosity,Viscosity;"
                "PropFluidDensity,Density" );

        // create ghost penalty fluid temperature
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPGPFluidTemp" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "follower_phase_name", "PhaseFluid" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        aParameterLists.set( "function_parameters", "0.05" );
        aParameterLists.set( "leader_properties", "PropFluidConductivity,Material" );
        aParameterLists.set( "follower_properties", "PropFluidConductivity,Material" );

        // create ghost fluid turbulence viscosity
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPGPViscosity" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "follower_phase_name", "PhaseFluid" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        aParameterLists.set( "function_parameters", "0.05" );
        aParameterLists.set( "leader_properties", "PropFluidKinViscosity,Material" );

        // create ghost penalty solid temperature
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPGPSolidTemp" );
        aParameterLists.set( "leader_phase_name", "PhaseSolid" );
        aParameterLists.set( "follower_phase_name", "PhaseSolid" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        aParameterLists.set( "function_parameters", "0.05" );
        aParameterLists.set( "leader_properties", "PropSolidConductivity,Material" );
        aParameterLists.set( "follower_properties", "PropSolidConductivity,Material" );

        // create ghost stabilization for theta and phi problems
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPGPFluidThetaPhi" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "follower_phase_name", "PhaseFluid" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        aParameterLists.set( "function_parameters", "0.01" );
        aParameterLists.set( "leader_properties", "PropConductivity,Material" );

        // create DBC on interface for theta problem
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPNitscheFluidThetaPhi" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "leader_properties", "PropConductivity,Material" );

        // create ghost stabilization parameter for theta and phi problems
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPGPInclusionThetaPhi" );
        aParameterLists.set( "leader_phase_name", "PhaseSolid" );
        aParameterLists.set( "follower_phase_name", "PhaseSolid" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        aParameterLists.set( "function_parameters", "0.01" );
        aParameterLists.set( "leader_properties", "PropConductivity,Material" );

        // create DBC on interface for theta problem
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPNitscheInclusionThetaPhi" );
        aParameterLists.set( "leader_phase_name", "PhaseSolid" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "leader_properties", "PropConductivity,Material" );

        //------------------------------------------------------------------------------
        // fill the IWG part of the parameter list

        // fluid bulk IWGs -------------------------------------------------------------
        // NS incompressible (velocity)
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGVelocityBulk" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_BULK );
        aParameterLists.set( "dof_residual", "VX,VY" );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists.set( "stabilization_parameters", "SPIncFlow,IncompressibleFlow" );

        // NS incompressible (pressure)
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGPressureBulk" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_BULK );
        aParameterLists.set( "dof_residual", "P" );
        aParameterLists.set( "leader_properties", "PropFluidPressureSpring,PressureSpring" );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists.set( "stabilization_parameters", "SPIncFlow,IncompressibleFlow" );

        // diffusion
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGFluidDiffusionBulk" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_constitutive_models", "CMFluidDiffusion,Diffusion" );

        // advection
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGFluidAdvectionBulk" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::ADVECTION_BULK );
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_constitutive_models", "CMFluidDiffusion,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPSUPGTemp,SUPG" );

        // turbulence
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGTurbulenceBulk" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPALART_ALLMARAS_TURBULENCE_BULK );
        aParameterLists.set( "dof_residual", "VISCOSITY" );
        aParameterLists.set( "leader_constitutive_models", "CMTurbulence,SpalartAllmarasTurbulence" );
        aParameterLists.set( "stabilization_parameters", "SPSUPGSA,SUPG" );

        // solid bulk IWGs -----------------------------------------------------------
        // diffusion
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGSolidDiffusionBulk" );
        aParameterLists.set( "leader_phase_name", "PhaseSolid" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_properties", "PropVolumetricHeatFlux,Load" );
        aParameterLists.set( "leader_constitutive_models", "CMSolidDiffusion,Diffusion" );

        // fluid/solid ifc IWGs ------------------------------------------------------

        // temperature
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGInterfaceFluidSolid" );
        aParameterLists.set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "follower_phase_name", "PhaseSolid" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_constitutive_models", "CMFluidDiffusion,Diffusion" );
        aParameterLists.set( "follower_constitutive_models", "CMSolidDiffusion,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPFSInterfaceNitsche,NitscheInterface" );

        // zero velocity
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGZeroVelocity" );
        aParameterLists.set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "neighbor_phases", "PhaseSolid" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "VX,VY" );
        aParameterLists.set( "leader_properties", "PropZeroU,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheU,DirichletNitsche" );

        // zero velocity (pressure part)
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGZeroPressure" );
        aParameterLists.set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "neighbor_phases", "PhaseSolid" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "P" );
        aParameterLists.set( "leader_properties", "PropZeroU,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );

        // zero viscosity
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGZeroViscosity" );
        aParameterLists.set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "neighbor_phases", "PhaseSolid" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPALART_ALLMARAS_TURBULENCE_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "VISCOSITY" );
        aParameterLists.set( "leader_properties", "PropZeroV,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMTurbulence,SpalartAllmarasTurbulence" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheV,Nitsche" );

        // Laterial BC IWG (adiabatic) ----------------------------------------------------------------

        // zero velocity along lateral sides
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGZeroVelocityVoid" );
        aParameterLists.set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "neighbor_phases", "PhaseVoidLateral" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "VX,VY" );
        aParameterLists.set( "leader_properties", "PropZeroU,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheU,DirichletNitsche" );

        // zero velocity along lateral sides (pressure part)
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGZeroPressureVoid" );
        aParameterLists.set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "neighbor_phases", "PhaseVoidLateral" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "P" );
        aParameterLists.set( "leader_properties", "PropZeroU,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );

        // zero viscosity along lateral sides
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGZeroViscosityVoid" );
        aParameterLists.set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "neighbor_phases", "PhaseVoidLateral" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPALART_ALLMARAS_TURBULENCE_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "VISCOSITY" );
        aParameterLists.set( "leader_properties", "PropZeroV,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMTurbulence,SpalartAllmarasTurbulence" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheV,Nitsche" );

        // Inlet BC IWG ----------------------------------------------------------------

        // inlet pressure
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGInletImposedPressure" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "neighbor_phases", "PhaseVoidFront" );
        aParameterLists.set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_IMPOSED_PRESSURE );
        aParameterLists.set( "dof_residual", "VX,VY" );
        aParameterLists.set( "leader_properties", "PropInletPressure,Pressure" );

        // inlet viscosity
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGInletViscosity" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "neighbor_phases", "PhaseVoidFront" );
        aParameterLists.set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPALART_ALLMARAS_TURBULENCE_DIRICHLET_SYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "VISCOSITY" );
        aParameterLists.set( "leader_properties", "PropInletV,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMTurbulence,SpalartAllmarasTurbulence" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheV,Nitsche" );

        // inlet temperature
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGInletTemp" );
        aParameterLists.set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "neighbor_phases", "PhaseVoidFront" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_properties", "PropInletTemp,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMFluidDiffusion,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheT,DirichletNitsche" );

        // Ghost  ----------------------------------------------------------------

        if ( tUseGhost )
        {
            // ghost viscous
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "IWGGPViscous" );
            aParameterLists.set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            aParameterLists.set( "leader_phase_name", "PhaseFluid" );
            aParameterLists.set( "follower_phase_name", "PhaseFluid" );
            aParameterLists.set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists.set( "dof_residual", "VX,VY" );
            aParameterLists.set( "stabilization_parameters", "SPGPViscous,GhostSP" );
            aParameterLists.set( "ghost_order", (uint)tDispOrder );

            // ghost convective
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "IWGGPConvective" );
            aParameterLists.set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            aParameterLists.set( "leader_phase_name", "PhaseFluid" );
            aParameterLists.set( "follower_phase_name", "PhaseFluid" );
            aParameterLists.set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists.set( "dof_residual", "VX,VY" );
            aParameterLists.set( "stabilization_parameters", "SPGPVelocity,GhostSP" );
            aParameterLists.set( "ghost_order", (uint)tDispOrder );

            // ghost pressure
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "IWGGPPressure" );
            aParameterLists.set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            aParameterLists.set( "leader_phase_name", "PhaseFluid" );
            aParameterLists.set( "follower_phase_name", "PhaseFluid" );
            aParameterLists.set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists.set( "dof_residual", "P" );
            aParameterLists.set( "stabilization_parameters", "SPGPPressure,GhostSP" );
            aParameterLists.set( "ghost_order", (uint)tDispOrder );

            // ghost fluid temperature
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "IWGGPFluidTemp" );
            aParameterLists.set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            aParameterLists.set( "leader_phase_name", "PhaseFluid" );
            aParameterLists.set( "follower_phase_name", "PhaseFluid" );
            aParameterLists.set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists.set( "dof_residual", "TEMP" );
            aParameterLists.set( "stabilization_parameters", "SPGPFluidTemp,GhostSP" );
            aParameterLists.set( "ghost_order", (uint)tDispOrder );

            // ghost fluid viscosity
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            aParameterLists.set( "IWG_name", "IWGGPFluidViscosity" );
            aParameterLists.set( "leader_phase_name", "PhaseFluid" );
            aParameterLists.set( "follower_phase_name", "PhaseFluid" );
            aParameterLists.set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists.set( "dof_residual", "VISCOSITY" );
            aParameterLists.set( "stabilization_parameters", "SPGPViscosity,GhostSP" );

            // ghost solid A temperature
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "IWGGPSolidTemp" );
            aParameterLists.set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            aParameterLists.set( "leader_phase_name", "PhaseSolid" );
            aParameterLists.set( "follower_phase_name", "PhaseSolid" );
            aParameterLists.set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists.set( "dof_residual", "TEMP" );
            aParameterLists.set( "stabilization_parameters", "SPGPSolidTemp,GhostSP" );
            aParameterLists.set( "ghost_order", (uint)tDispOrder );
        }

        // theta problem  ----------------------------------------------------------------
        // theta bulk in fluid
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGFluidDiffusionThetaBulk" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        aParameterLists.set( "dof_residual", "THETA" );
        aParameterLists.set( "leader_constitutive_models", "CMFluidDiffusionTheta,Diffusion" );

        // theta bulk in inclusion
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGInclusionSolidDiffusionThetaBulk" );
        aParameterLists.set( "leader_phase_name", "PhaseSolid" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        aParameterLists.set( "dof_residual", "THETA" );
        aParameterLists.set( "leader_constitutive_models", "CMInclusionSolidDiffusionTheta,Diffusion" );

        // create parameter list for single side interface condition
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGSurfaceOuterTheta" );
        aParameterLists.set( "leader_phase_name", "PhaseSolid" );
        aParameterLists.set( "neighbor_phases", "PhaseFluid" );
        aParameterLists.set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "THETA" );
        aParameterLists.set( "leader_properties", "PropPrescTheta,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMInclusionSolidDiffusionTheta,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheInclusionThetaPhi,DirichletNitsche" );

        // create parameter list for single side interface condition
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGSurfaceInnerTheta" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "neighbor_phases", "PhaseSolid,PhaseVoidLateral" );
        aParameterLists.set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "THETA" );
        aParameterLists.set( "leader_properties", "PropPrescTheta,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMFluidDiffusionTheta,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheFluidThetaPhi,DirichletNitsche" );

        if ( tUseGhost )
        {
            // create IWG - ghost
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "IWGGPInnerTheta" );
            aParameterLists.set( "leader_phase_name", "PhaseFluid" );
            aParameterLists.set( "follower_phase_name", "PhaseFluid" );
            aParameterLists.set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            aParameterLists.set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists.set( "dof_residual", "THETA" );
            aParameterLists.set( "stabilization_parameters", "SPGPFluidThetaPhi,GhostSP" );

            // create IWG - ghost
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "IWGGPOuterTheta" );
            aParameterLists.set( "leader_phase_name", "PhaseSolid" );
            aParameterLists.set( "follower_phase_name", "PhaseSolid" );
            aParameterLists.set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            aParameterLists.set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists.set( "dof_residual", "THETA" );
            aParameterLists.set( "stabilization_parameters", "SPGPInclusionThetaPhi,GhostSP" );
        }

        // create time side interface condition fluid
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGFluidTimeContinuityTheta" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::TIME_CONTINUITY_DOF );
        aParameterLists.set( "dof_residual", "THETA" );
        aParameterLists.set( "leader_properties",
                "PropWeightCurrent,       WeightCurrent;"
                "PropWeightPrevious,      WeightPrevious;"
                "PropInitialCondition,    InitialCondition" );

        // create time side interface condition solid
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGInclusionTimeContinuityTheta" );
        aParameterLists.set( "leader_phase_name", "PhaseSolid" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::TIME_CONTINUITY_DOF );
        aParameterLists.set( "dof_residual", "THETA" );
        aParameterLists.set( "leader_properties",
                "PropWeightCurrent,       WeightCurrent;"
                "PropWeightPrevious,      WeightPrevious;"
                "PropInitialCondition,    InitialCondition" );

        // theta problem  ----------------------------------------------------------------
        // create IWG - bulk diffusion
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGDiffusionInnerBulk" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        aParameterLists.set( "dof_residual", "PHID" );
        aParameterLists.set( "leader_constitutive_models", "CMFluidDiffusionPhi,Diffusion" );

        // create IWG - bulk diffusion
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGDiffusionOuterBulk" );
        aParameterLists.set( "leader_phase_name", "PhaseSolid" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        aParameterLists.set( "dof_residual", "PHID" );
        aParameterLists.set( "leader_constitutive_models", "CMInclusionSolidDiffusionPhi,Diffusion" );

        // create parameter list for single side interface condition
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGSurfaceInnerPhi" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "neighbor_phases", "PhaseSolid,PhaseVoidLateral" );
        aParameterLists.set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "PHID" );
        aParameterLists.set( "leader_properties", "PropPrescPhi,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMFluidDiffusionPhi,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheFluidThetaPhi,DirichletNitsche" );

        // create parameter list for single side interface condition
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGSurfaceOuterPhi" );
        aParameterLists.set( "leader_phase_name", "PhaseSolid" );
        aParameterLists.set( "neighbor_phases", "PhaseFluid" );
        aParameterLists.set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "PHID" );
        aParameterLists.set( "leader_properties", "PropPrescPhi,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMInclusionSolidDiffusionPhi,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheInclusionThetaPhi,DirichletNitsche" );

        if ( tUseGhost )
        {
            // create IWG - ghost
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "IWGGPInnerPhi" );
            aParameterLists.set( "leader_phase_name", "PhaseFluid" );
            aParameterLists.set( "follower_phase_name", "PhaseFluid" );
            aParameterLists.set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            aParameterLists.set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists.set( "dof_residual", "PHID" );
            aParameterLists.set( "stabilization_parameters", "SPGPFluidThetaPhi,GhostSP" );

            // create IWG - ghost
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "IWGGPOuterPhi" );
            aParameterLists.set( "leader_phase_name", "PhaseSolid" );
            aParameterLists.set( "follower_phase_name", "PhaseSolid" );
            aParameterLists.set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            aParameterLists.set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists.set( "dof_residual", "PHID" );
            aParameterLists.set( "stabilization_parameters", "SPGPInclusionThetaPhi,GhostSP" );
        }

        //------------------------------------------------------------------------------
        // fill the IQI part of the parameter list

        // velocity VX
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkVX" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "VX,VY" );
        aParameterLists.set( "vectorial_field_index", 0 );

        // velocity VY
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkVY" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "VX,VY" );
        aParameterLists.set( "vectorial_field_index", 1 );

        // pressure
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkP" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "P" );
        aParameterLists.set( "vectorial_field_index", 0 );

        // temperature
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkTEMP" );
        aParameterLists.set( "leader_phase_name", "PhaseAll" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "TEMP" );
        aParameterLists.set( "vectorial_field_index", 0 );

        // viscosity
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkVISCOSITY" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "VISCOSITY" );
        aParameterLists.set( "vectorial_field_index", 0 );

        // theta
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkTHETA" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "THETA" );
        aParameterLists.set( "vectorial_field_index", 0 );

        // phi
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkPHID" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "PHID" );
        aParameterLists.set( "vectorial_field_index", 0 );

        // fluid thermal energy on inlet
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIInletThermalEnergy" );
        aParameterLists.set( "IQI_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "neighbor_phases", "PhaseVoidFront" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::THERMAL_ENERGY_CONVECTIVE_FLUX );
        aParameterLists.set( "leader_constitutive_models", "CMFluidDiffusion,Diffusion" );

        // fluid thermal energy on outlet
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIOutletThermalEnergy" );
        aParameterLists.set( "IQI_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "IQI_type", fem::IQI_Type::THERMAL_ENERGY_CONVECTIVE_FLUX );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "neighbor_phases", "PhaseVoidBack" );
        aParameterLists.set( "leader_constitutive_models", "CMFluidDiffusion,Diffusion" );

        // fluid total pressure on inlet
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIInletTotalPressure" );
        aParameterLists.set( "IQI_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "neighbor_phases", "PhaseVoidFront" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::TOTAL_PRESSURE );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,Fluid" );

        // fluid total pressure on outlet
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIOutletTotalPressure" );
        aParameterLists.set( "IQI_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "neighbor_phases", "PhaseVoidBack" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::TOTAL_PRESSURE );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,Fluid" );

        // fluid mass flow on inlet
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIInletMassFlow" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::MASS_FLOW );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "neighbor_phases", "PhaseVoidFront" );
        aParameterLists.set( "IQI_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,Fluid" );

        // fluid mass flow on outlet
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIOutletMassFlow" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "neighbor_phases", "PhaseVoidBack" );
        aParameterLists.set( "IQI_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "IQI_type", fem::IQI_Type::MASS_FLOW );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,Fluid" );

        // fluid power dissipation on inlet
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIInletPowDisp" );
        aParameterLists.set( "IQI_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "neighbor_phases", "PhaseVoidFront" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::POWER_DISSIPATION );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,Fluid" );

        // fluid power dissipation on outlet
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIOutletPowDisp" );
        aParameterLists.set( "IQI_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "neighbor_phases", "PhaseVoidBack" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::POWER_DISSIPATION );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,Fluid" );

        // fluid power dissipation in volume
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIVolumePowDisp" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::POWER_DISSIPATION_BULK );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,Fluid" );

        // inclusion perimeter
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIPerimeterItf" );
        aParameterLists.set( "leader_phase_name", "PhaseSolid" );
        aParameterLists.set( "neighbor_phases", "PhaseFluid" );
        aParameterLists.set( "IQI_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "IQI_type", fem::IQI_Type::VOLUME );

        // max temperature in solid
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIMaxTemp" );
        aParameterLists.set( "leader_phase_name", "PhaseSolid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::MAX_DOF );
        aParameterLists.set( "function_parameters", tRefTemp + "/" + tMaxTempExponent );
        aParameterLists.set( "dof_quantity", "TEMP" );
        aParameterLists.set( "vectorial_field_index", 0 );

        // volume solid
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQISolidVolume" );
        aParameterLists.set( "leader_phase_name", "PhaseSolid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::VOLUME );

        // wall distance
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIWallDistance" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::PROPERTY );
        aParameterLists.set( "leader_properties", "PropWallDistance,Property" );

        // turbulent dynamic viscosity
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkTurbDynVisc" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::TURBULENT_DYNAMIC_VISCOSITY );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,Fluid_Turbulence" );

        // effective dynamic viscosity
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkEffDynVisc" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::EFFECTIVE_DYNAMIC_VISCOSITY );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,Fluid_Turbulence" );

        // effective conductivity
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkEffCond" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::EFFECTIVE_CONDUCTIVITY );
        aParameterLists.set( "leader_constitutive_models", "CMFluidDiffusion,Diffusion_Turbulence" );

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        aParameterLists( FEM::COMPUTATION );
        aParameterLists.set( "print_physics_model", false );

        aParameterLists.set( "finite_difference_scheme", fem::FDScheme_Type::POINT_3_CENTRAL );
        aParameterLists.set( "finite_difference_perturbation_size", tFEMFdEpsilon );
    }

    void
    SOLParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::AMESOS_IMPL );

#ifdef MORIS_USE_MUMPS
        aParameterLists.set( "Solver_Type", "Amesos_Mumps" );
#else
        aParameterLists.set( "Solver_Type", "Amesos_Superludist" );
#endif

        //------------------------------------------------------------------------------

        aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "DLA_Linear_solver_algorithms", "0" );

        //------------------------------------------------------------------------------

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();    // nonlinear algorithm index 0
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists.set( "NLA_Linear_solver", 0 );
        aParameterLists.set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        aParameterLists.set( "NLA_relaxation_strategy", sol::SolverRelaxationType::InvResNormAdaptive );
        aParameterLists.set( "NLA_relaxation_parameter", 0.5 );
        aParameterLists.set( "NLA_relaxation_damping", 0.5 );
        aParameterLists.set( "NLA_max_iter", 100 );

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists.set( "NLA_Linear_solver", 0 );
        aParameterLists.set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        aParameterLists.set( "NLA_relaxation_parameter", 1.0 );
        aParameterLists.set( "NLA_max_iter", 1 );

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NLBGS_SOLVER );
        aParameterLists.set( "NLA_Linear_solver", 0 );
        aParameterLists.set( "NLA_rel_res_norm_drop", 1.0 );
        aParameterLists.set( "NLA_max_iter", 1 );

        //------------------------------------------------------------------------------

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "NLA_Nonlinear_solver_algorithms", "1" );
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists.set( "NLA_DofTypes", "THETA" );

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "NLA_Nonlinear_solver_algorithms", "1" );
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists.set( "NLA_DofTypes", "PHID" );

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "NLA_Nonlinear_solver_algorithms", "0" );
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists.set( "NLA_DofTypes", "VX,VY,P,VISCOSITY" );

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "NLA_Nonlinear_solver_algorithms", "1" );
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists.set( "NLA_DofTypes", "TEMP" );

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "NLA_Nonlinear_solver_algorithms", "2" );
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NLBGS_SOLVER );
        aParameterLists.set( "NLA_Sub_Nonlinear_Solver", "0,1,2,3" );
        aParameterLists.set( "NLA_DofTypes", "THETA;PHID;VX,VY,P,VISCOSITY;TEMP" );

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "NLA_Nonlinear_solver_algorithms", "1" );
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists.set( "NLA_DofTypes", "VX,VY,P,VISCOSITY" );

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "NLA_Nonlinear_solver_algorithms", "2" );
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NLBGS_SOLVER );
        aParameterLists.set( "NLA_Sub_Nonlinear_Solver", "0,1,5,3" );
        aParameterLists.set( "NLA_DofTypes", "THETA;PHID;VX,VY,P,VISCOSITY;TEMP" );

        // ----------------------------------------------------------

        aParameterLists( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "TSA_Nonlinear_Solver", 4 );                // using NLBGS for forward problem
        aParameterLists.set( "TSA_Nonlinear_Sensitivity_Solver", 6 );    // using monlithic for sensitivity problem

        //------------------------------------------------------------------------------

        aParameterLists( SOL::TIME_SOLVERS ).add_parameter_list();
        aParameterLists.set( "TSA_DofTypes", "THETA;PHID;VX,VY,P,VISCOSITY;TEMP" );
        aParameterLists.set( "TSA_Initialize_Sol_Vec", "THETA,0.0;PHID,0.0;VX,0.01;VY,0.0;P,0.0;VISCOSITY," + tInletKinViscosity + ";TEMP,0.0" );
        aParameterLists.set( "TSA_Output_Indices", "0" );
        aParameterLists.set( "TSA_Output_Criteria", "Output_Criterion" );

        //------------------------------------------------------------------------------

        aParameterLists( SOL::PRECONDITIONERS ).add_parameter_list( sol::PreconditionerType::NONE );
    }

    void
    MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "THETA", 1 );
        aParameterLists.set( "PHID", 1 );
        aParameterLists.set( "VX", 1 );
        aParameterLists.set( "VY", 1 );
        aParameterLists.set( "P", 1 );
        aParameterLists.set( "VISCOSITY", 1 );
        aParameterLists.set( "TEMP", 1 );
    }

    void
    VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "File_Name", std::pair< std::string, std::string >( "./", tExoFile ) );
        aParameterLists.set( "Mesh_Type", vis::VIS_Mesh_Type::STANDARD );
        aParameterLists.set( "Set_Names", tAllPhases );
        aParameterLists.set( "Field_Names", "THETA,PHID,VX,VY,P,VISCOSITY,TEMP,TURBDYNVISC,EFFVISC,EFFCOND" );
        aParameterLists.set( "Field_Type", "NODAL,NODAL,NODAL,NODAL,NODAL,NODAL,NODAL,NODAL,NODAL,NODAL" );
        aParameterLists.set( "IQI_Names",
                "IQIBulkTHETA,IQIBulkPHID,IQIBulkVX,IQIBulkVY,IQIBulkP,IQIBulkVISCOSITY,IQIBulkTEMP,"
                "IQIBulkTurbDynVisc,IQIBulkEffDynVisc,IQIBulkEffCond" );
        aParameterLists.set( "Save_Frequency", 1 );
        aParameterLists.set( "Time_Offset", 10.0 );
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
