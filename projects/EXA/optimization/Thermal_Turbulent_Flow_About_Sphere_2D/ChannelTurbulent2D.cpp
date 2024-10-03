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
    OPTParameterList( Vector< Submodule_Parameter_Lists >& tParameterlist )
    {
        tParameterlist.resize( 3 );

        tParameterlist( 0 ).add_parameter_list( prm::create_opt_problem_parameter_list() );
        tParameterlist( 0 ).set( "is_optimization_problem", true );
        tParameterlist( 0 ).set( "problem", "user_defined" );
        tParameterlist( 0 ).set( "library", tSoFile );


        tParameterlist( 2 ).add_parameter_list( moris::prm::create_sweep_parameter_list() );
        tParameterlist( 2 ).set( "hdf5_path", tHdf5File );
        tParameterlist( 2 ).set( "num_evaluations_per_adv", "1" );
        tParameterlist( 2 ).set( "finite_difference_type", "all" );
        tParameterlist( 2 ).set( "finite_difference_epsilons", tSweepFdEpsilon );
    }

    /* ------------------------------------------------------------------------ */

    void
    HMRParameterList( Vector< Submodule_Parameter_Lists >& tParameterlist )
    {
        tParameterlist.resize( 1 );

        tParameterlist( 0 ).add_parameter_list( prm::create_hmr_parameter_list() );

        tParameterlist( 0 ).set( "number_of_elements_per_dimension", tNumElemsPerDim );
        tParameterlist( 0 ).set( "domain_dimensions", tDomainDims );
        tParameterlist( 0 ).set( "domain_offset", tDomainOffset );
        tParameterlist( 0 ).set( "domain_sidesets", tDomainSidesets );

        tParameterlist( 0 ).set( "lagrange_output_meshes", "0" );

        tParameterlist( 0 ).set( "lagrange_orders", tLagrangeOrder );
        tParameterlist( 0 ).set( "lagrange_pattern", tLagrangePattern );

        tParameterlist( 0 ).set( "bspline_orders", tBsplineOrder );
        tParameterlist( 0 ).set( "bspline_pattern", "0,1" );

        tParameterlist( 0 ).set( "lagrange_to_bspline", "0,1" );

        tParameterlist( 0 ).set( "truncate_bsplines", 1 );
        tParameterlist( 0 ).set( "refinement_buffer", tRefineBuffer );
        tParameterlist( 0 ).set( "staircase_buffer", tRefineBuffer );

        tParameterlist( 0 ).set( "initial_refinement", tInitialRef );
        tParameterlist( 0 ).set( "initial_refinement_pattern", "0,1" );

        tParameterlist( 0 ).set( "use_number_aura", 1 );

        tParameterlist( 0 ).set( "use_multigrid", 0 );
        tParameterlist( 0 ).set( "severity_level", 0 );
    }

    /* ------------------------------------------------------------------------ */

    void
    XTKParameterList( Vector< Submodule_Parameter_Lists >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).add_parameter_list( prm::create_xtk_parameter_list() );
        tParameterlist( 0 ).set( "decompose", true );
        tParameterlist( 0 ).set( "decomposition_type", "conformal" );
        tParameterlist( 0 ).set( "enrich", true );
        tParameterlist( 0 ).set( "basis_rank", "bspline" );
        tParameterlist( 0 ).set( "enrich_mesh_indices", "0,1" );
        tParameterlist( 0 ).set( "ghost_stab", tUseGhost );
        tParameterlist( 0 ).set( "multigrid", false );
        tParameterlist( 0 ).set( "verbose", false );
        tParameterlist( 0 ).set( "high_to_low_dbl_side_sets", true );
        tParameterlist( 0 ).set( "print_enriched_ig_mesh", false );
        tParameterlist( 0 ).set( "exodus_output_XTK_ig_mesh", true );
    }

    /* ------------------------------------------------------------------------ */

    void
    GENParameterList( Vector< Submodule_Parameter_Lists >& tParameterlist )
    {
        tParameterlist.resize( 3 );

        tParameterlist( 0 ).add_parameter_list( prm::create_gen_parameter_list() );
        tParameterlist( 0 ).set( "IQI_types",
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

        tParameterlist( 0 ).set( "number_of_phases", 5 );
        tParameterlist( 0 ).set( "phase_function_name", "get_phase_index" );
        tParameterlist( 0 ).set( "output_mesh_file", tGENOutputFile );
        tParameterlist( 0 ).set( "time_offset", 10.0 );

        // Inclusions
        tParameterlist( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        tParameterlist( 1 ).set( "field_function_name", "Func_Sphere" );
        tParameterlist( 1 ).set( "sensitivity_function_name", "Func_Sphere_Deriv" );
        tParameterlist( 1 ).insert( "radius", Design_Variable( tSphereRadius * 0.9, tSphereRadius, tSphereRadius * 1.1 ) );
        tParameterlist( 1 ).set( "number_of_refinements", tInterfaceRefinementSphere );
        tParameterlist( 1 ).set( "refinement_mesh_index", 0 );
        tParameterlist( 1 ).set( "isocontour_threshold", 0.0 );
        tParameterlist( 1 ).set( "isocontour_tolerance", 1.0e-12 );
        tParameterlist( 1 ).set( "intersection_tolerance", 1.0e-12 );

        // Inlet plane
        tParameterlist( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        tParameterlist( 1 ).set( "field_function_name", "Func_Plane" );
        tParameterlist( 1 ).insert( "variable_1", 1.0 );
        tParameterlist( 1 ).insert( "variable_2", 0.0 );
        tParameterlist( 1 ).insert( "variable_3", tGeoShift );
        tParameterlist( 1 ).insert( "variable_4", 0.0 );
        tParameterlist( 1 ).set( "number_of_refinements", tInterfaceRefinementWalls );
        tParameterlist( 1 ).set( "refinement_mesh_index", 0 );
        tParameterlist( 1 ).set( "isocontour_threshold", 0.0 );
        tParameterlist( 1 ).set( "isocontour_tolerance", 1.0e-12 );
        tParameterlist( 1 ).set( "intersection_tolerance", 1.0e-12 );

        // Out plane
        tParameterlist( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        tParameterlist( 1 ).set( "field_function_name", "Func_Plane" );
        tParameterlist( 1 ).insert( "variable_1", -1.0 );
        tParameterlist( 1 ).insert( "variable_2", 0.0 );
        tParameterlist( 1 ).insert( "variable_3", tChannelLength + tGeoShift );
        tParameterlist( 1 ).insert( "variable_4", 0.0 );
        tParameterlist( 1 ).set( "number_of_refinements", tInterfaceRefinementWalls );
        tParameterlist( 1 ).set( "refinement_mesh_index", 0 );
        tParameterlist( 1 ).set( "isocontour_threshold", 0.0 );
        tParameterlist( 1 ).set( "isocontour_tolerance", 1.0e-12 );
        tParameterlist( 1 ).set( "intersection_tolerance", 1.0e-12 );

        // Lower plane
        tParameterlist( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        tParameterlist( 1 ).set( "field_function_name", "Func_Plane" );
        tParameterlist( 1 ).insert( "variable_1", 0.0 );
        tParameterlist( 1 ).insert( "variable_2", 1.0 );
        tParameterlist( 1 ).insert( "variable_3", 0.0 );
        tParameterlist( 1 ).insert( "variable_4", tGeoShift );
        tParameterlist( 1 ).set( "number_of_refinements", tInterfaceRefinementWalls );
        tParameterlist( 1 ).set( "refinement_mesh_index", 0 );
        tParameterlist( 1 ).set( "isocontour_threshold", 0.0 );
        tParameterlist( 1 ).set( "isocontour_tolerance", 1.0e-12 );
        tParameterlist( 1 ).set( "intersection_tolerance", 1.0e-12 );

        // Upper plane
        tParameterlist( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        tParameterlist( 1 ).set( "field_function_name", "Func_Plane" );
        tParameterlist( 1 ).insert( "variable_1", 0.0 );
        tParameterlist( 1 ).insert( "variable_2", -1.0 );
        tParameterlist( 1 ).insert( "variable_3", 0.0 );
        tParameterlist( 1 ).insert( "variable_4", tChannelHeight + tGeoShift );
        tParameterlist( 1 ).set( "number_of_refinements", tInterfaceRefinementWalls );
        tParameterlist( 1 ).set( "refinement_mesh_index", 0 );
        tParameterlist( 1 ).set( "isocontour_threshold", 0.0 );
        tParameterlist( 1 ).set( "isocontour_tolerance", 1.0e-12 );
        tParameterlist( 1 ).set( "intersection_tolerance", 1.0e-12 );
    }

    /* ------------------------------------------------------------------------ */

    void
    FEMParameterList( Vector< Submodule_Parameter_Lists >& tParameterList )
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

        tParameterList( tPhaseIndex ).add_parameter_list( prm::create_phase_parameter_list() );
        tParameterList( tPhaseIndex ).set( "phase_name", "PhaseFluid" );
        tParameterList( tPhaseIndex ).set( "phase_indices", "0" );

        tParameterList( tPhaseIndex ).add_parameter_list( prm::create_phase_parameter_list() );
        tParameterList( tPhaseIndex ).set( "phase_name", "PhaseSolid" );
        tParameterList( tPhaseIndex ).set( "phase_indices", "1" );

        tParameterList( tPhaseIndex ).add_parameter_list( prm::create_phase_parameter_list() );
        tParameterList( tPhaseIndex ).set( "phase_name", "PhaseVoidFront" );
        tParameterList( tPhaseIndex ).set( "phase_indices", "2" );

        tParameterList( tPhaseIndex ).add_parameter_list( prm::create_phase_parameter_list() );
        tParameterList( tPhaseIndex ).set( "phase_name", "PhaseVoidBack" );
        tParameterList( tPhaseIndex ).set( "phase_indices", "3" );

        tParameterList( tPhaseIndex ).add_parameter_list( prm::create_phase_parameter_list() );
        tParameterList( tPhaseIndex ).set( "phase_name", "PhaseVoidLateral" );
        tParameterList( tPhaseIndex ).set( "phase_indices", "4" );

        tParameterList( tPhaseIndex ).add_parameter_list( prm::create_phase_parameter_list() );
        tParameterList( tPhaseIndex ).set( "phase_name", "PhaseAll" );
        tParameterList( tPhaseIndex ).set( "phase_indices", "0,1" );

        //------------------------------------------------------------------------------
        // fill the property part of the parameter list


        // fluid properties ------------------------------------------------------------
        // create fluid viscosity property
        tParameterList( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( tPropIndex ).set( "property_name", "PropFluidDynViscosity" );
        tParameterList( tPropIndex ).set( "function_parameters", tFluidDynViscosity );

        // create fluid density property
        tParameterList( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( tPropIndex ).set( "property_name", "PropFluidDensity" );
        tParameterList( tPropIndex ).set( "function_parameters", tFluidDensity );

        // create fluid kinematic viscosity property
        tParameterList( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( tPropIndex ).set( "property_name", "PropFluidKinViscosity" );
        tParameterList( tPropIndex ).set( "function_parameters", tFluidKinViscosity );

        // create fluid capacity property
        tParameterList( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( tPropIndex ).set( "property_name", "PropFluidCapacity" );
        tParameterList( tPropIndex ).set( "function_parameters", tFluidCapacity );

        // create fluid conductivity property
        tParameterList( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( tPropIndex ).set( "property_name", "PropFluidConductivity" );
        tParameterList( tPropIndex ).set( "function_parameters", tFluidConductivity );

        // create turbulent prandtl number
        tParameterList( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( tPropIndex ).set( "property_name", "PropFluidTurbPrandtl" );
        tParameterList( tPropIndex ).set( "function_parameters", tFluidTurbulencePrandtl );

        // create fluid pressure spring property
        tParameterList( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( tPropIndex ).set( "property_name", "PropFluidPressureSpring" );
        tParameterList( tPropIndex ).set( "function_parameters", tFluidPressureSpring );

        // solid properties ----------------------------------------------------------
        // create solid B density property
        tParameterList( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( tPropIndex ).set( "property_name", "PropSolidDensity" );
        tParameterList( tPropIndex ).set( "function_parameters", tSolidDensity );

        // create solid B capacity property
        tParameterList( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( tPropIndex ).set( "property_name", "PropSolidCapacity" );
        tParameterList( tPropIndex ).set( "function_parameters", tSolidCapacity );

        // create solid B conductivity property
        tParameterList( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( tPropIndex ).set( "property_name", "PropSolidConductivity" );
        tParameterList( tPropIndex ).set( "function_parameters", tSolidConductivity );

        // BC properties ---------------------------------------------------------------
        // create inlet pressure property
        tParameterList( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( tPropIndex ).set( "property_name", "PropInletPressure" );
        tParameterList( tPropIndex ).set( "function_parameters", tInletPressure );

        // create inlet temperature property
        tParameterList( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( tPropIndex ).set( "property_name", "PropInletTemp" );
        tParameterList( tPropIndex ).set( "function_parameters", tInletTemperature );

        // create wall velocity property
        tParameterList( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( tPropIndex ).set( "property_name", "PropZeroU" );
        tParameterList( tPropIndex ).set( "function_parameters", "0.0;0.0" );

        // create symmetry velocity property (x-z plane)
        tParameterList( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( tPropIndex ).set( "property_name", "PropSelectUY" );
        tParameterList( tPropIndex ).set( "function_parameters", "0.0,0.0;0.0,1.0" );

        // create inlet kinematic viscosity property
        tParameterList( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( tPropIndex ).set( "property_name", "PropInletV" );
        tParameterList( tPropIndex ).set( "value_function", "Func_Inlet_V" );
        tParameterList( tPropIndex ).set( "function_parameters", tInletKinViscosity );

        // create  wall kinematic viscosity property
        tParameterList( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( tPropIndex ).set( "property_name", "PropZeroV" );
        tParameterList( tPropIndex ).set( "function_parameters", "0.0" );

        // create heat load property
        tParameterList( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( tPropIndex ).set( "property_name", "PropVolumetricHeatFlux" );
        tParameterList( tPropIndex ).set( "function_parameters", tVolumetricHeatLoad );

        // Wall distance properties ----------------------------------------------------
        // create wall distance property
        tParameterList( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( tPropIndex ).set( "property_name", "PropWallDistance" );
        tParameterList( tPropIndex ).set( "dof_dependencies", "PHID" );
        tParameterList( tPropIndex ).set( "value_function", "Func_Wall_Distance" );
        tParameterList( tPropIndex ).set( "dof_derivative_functions", "Func_Wall_Distance_Der" );

        // create common conductivity property for theta and phi problems
        tParameterList( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( tPropIndex ).set( "property_name", "PropConductivity" );
        tParameterList( tPropIndex ).set( "function_parameters", tConductivityThetaPhi );

        // create density property for theta problem
        tParameterList( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( tPropIndex ).set( "property_name", "PropDensityTheta" );
        tParameterList( tPropIndex ).set( "function_parameters", tDensityTheta );

        // create capacity property for theta problem
        tParameterList( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( tPropIndex ).set( "property_name", "PropCapacityTheta" );
        tParameterList( tPropIndex ).set( "function_parameters", tCapacityTheta );

        // create prescribed BC property for theta problem
        tParameterList( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( tPropIndex ).set( "property_name", "PropPrescTheta" );
        tParameterList( tPropIndex ).set( "function_parameters", tPrescTheta );

        // create density property for phi problem
        tParameterList( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( tPropIndex ).set( "property_name", "PropDensityPhi" );
        tParameterList( tPropIndex ).set( "function_parameters", tDensityPhi );

        // create capacity property for phi problem
        tParameterList( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( tPropIndex ).set( "property_name", "PropCapacityPhi" );
        tParameterList( tPropIndex ).set( "function_parameters", tCapacityPhi );

        // create prescribed BC property for phi problem
        tParameterList( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( tPropIndex ).set( "property_name", "PropPrescPhi" );
        tParameterList( tPropIndex ).set( "function_parameters", tPrescPhi );

        // create eigenstrain BC property for phi problem
        tParameterList( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( tPropIndex ).set( "property_name", "PropEigenStrainPhi" );
        tParameterList( tPropIndex ).set( "function_parameters", "1.0" );

        // create time continuity weight current
        tParameterList( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( tPropIndex ).set( "property_name", "PropWeightCurrent" );
        tParameterList( tPropIndex ).set( "function_parameters", "10.0" );

        // create time continuity weight previous
        tParameterList( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( tPropIndex ).set( "property_name", "PropWeightPrevious" );
        tParameterList( tPropIndex ).set( "function_parameters", "10.0" );

        // create initial condition property
        tParameterList( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( tPropIndex ).set( "property_name", "PropInitialCondition" );
        tParameterList( tPropIndex ).set( "function_parameters", "0.0" );

        //------------------------------------------------------------------------------
        // fill the constitutive model part of the parameter list


        // fluid CM --------------------------------------------------------------------
        // create fluid CM
        tParameterList( tCMIndex ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        tParameterList( tCMIndex ).set( "constitutive_name", "CMFluid" );
        tParameterList( tCMIndex ).set( "phase_name", "PhaseFluid" );
        tParameterList( tCMIndex ).set( "constitutive_type", fem::Constitutive_Type::FLUID_TURBULENCE );
        tParameterList( tCMIndex ).set( "dof_dependencies", std::pair< std::string, std::string >( "VX,VY;P;VISCOSITY", "Velocity,Pressure,Viscosity" ) );
        tParameterList( tCMIndex ).set( "properties",
                "PropFluidDynViscosity,Viscosity;"
                "PropFluidKinViscosity,KinViscosity;"
                "PropFluidDensity  ,Density" );

        // create CM SA turbulence
        tParameterList( tCMIndex ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        tParameterList( tCMIndex ).set( "constitutive_name", "CMTurbulence" );
        tParameterList( tCMIndex ).set( "phase_name", "PhaseFluid" );
        tParameterList( tCMIndex ).set( "constitutive_type", fem::Constitutive_Type::SPALART_ALLMARAS_TURBULENCE );
        tParameterList( tCMIndex ).set( "function_parameters", tCMTurbFt2 + "/" + tCMTurbAlpha );
        tParameterList( tCMIndex ).set( "dof_dependencies", std::pair< std::string, std::string >( "VX,VY;VISCOSITY", "Velocity,Viscosity" ) );
        tParameterList( tCMIndex ).set( "properties",
                "PropFluidKinViscosity,KinViscosity;"
                "PropWallDistance     ,WallDistance" );

        // create fluid diffusion CM
        tParameterList( tCMIndex ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        tParameterList( tCMIndex ).set( "constitutive_name", "CMFluidDiffusion" );
        tParameterList( tCMIndex ).set( "phase_name", "PhaseFluid" );
        tParameterList( tCMIndex ).set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO_TURBULENCE );
        tParameterList( tCMIndex ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( tCMIndex ).set( "properties",
                "PropFluidConductivity,Conductivity;"
                "PropFluidDensity     ,Density;"
                "PropFluidCapacity    ,HeatCapacity;"
                "PropFluidKinViscosity,KinematicViscosity;"
                "PropFluidTurbPrandtl ,TurbulentPrandtl" );

        // solid CM --------------------------------------------------------------------
        // create solid B diffusion CM
        tParameterList( tCMIndex ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        tParameterList( tCMIndex ).set( "constitutive_name", "CMSolidDiffusion" );
        tParameterList( tCMIndex ).set( "phase_name", "PhaseSolid" );
        tParameterList( tCMIndex ).set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
        tParameterList( tCMIndex ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( tCMIndex ).set( "properties",
                "PropSolidConductivity,Conductivity;"
                "PropSolidDensity     ,Density;"
                "PropSolidCapacity    ,HeatCapacity" );

        // theta/phi CM --------------------------------------------------------------------
        // create fluid theta diffusion CM
        tParameterList( tCMIndex ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        tParameterList( tCMIndex ).set( "constitutive_name", "CMFluidDiffusionTheta" );
        tParameterList( tCMIndex ).set( "phase_name", "PhaseFluid" );
        tParameterList( tCMIndex ).set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
        tParameterList( tCMIndex ).set( "dof_dependencies", std::pair< std::string, std::string >( "THETA", "Temperature" ) );
        tParameterList( tCMIndex ).set( "properties",
                "PropConductivity      , Conductivity;"
                "PropDensityTheta      , Density;"
                "PropCapacityTheta     , HeatCapacity" );

        // create solid theta diffusion CM
        tParameterList( tCMIndex ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        tParameterList( tCMIndex ).set( "constitutive_name", "CMInclusionSolidDiffusionTheta" );
        tParameterList( tCMIndex ).set( "phase_name", "PhaseSolid" );
        tParameterList( tCMIndex ).set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
        tParameterList( tCMIndex ).set( "dof_dependencies", std::pair< std::string, std::string >( "THETA", "Temperature" ) );
        tParameterList( tCMIndex ).set( "properties",
                "PropConductivity      , Conductivity;"
                "PropDensityTheta      , Density;"
                "PropCapacityTheta     , HeatCapacity" );

        // create fluid phi diffusion CM
        tParameterList( tCMIndex ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        tParameterList( tCMIndex ).set( "constitutive_name", "CMFluidDiffusionPhi" );
        tParameterList( tCMIndex ).set( "phase_name", "PhaseFluid" );
        tParameterList( tCMIndex ).set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
        tParameterList( tCMIndex ).set( "dof_dependencies", std::pair< std::string, std::string >( "PHID;THETA", "Temperature,Theta" ) );
        tParameterList( tCMIndex ).set( "properties",
                "PropConductivity    , Conductivity;"
                "PropDensityPhi      , Density;"
                "PropCapacityPhi     , HeatCapacity;"
                "PropEigenStrainPhi  , EigenStrain" );

        // create solid phi diffusion CM
        tParameterList( tCMIndex ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        tParameterList( tCMIndex ).set( "constitutive_name", "CMInclusionSolidDiffusionPhi" );
        tParameterList( tCMIndex ).set( "phase_name", "PhaseSolid" );
        tParameterList( tCMIndex ).set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
        tParameterList( tCMIndex ).set( "dof_dependencies", std::pair< std::string, std::string >( "PHID;THETA", "Temperature,Theta" ) );
        tParameterList( tCMIndex ).set( "properties",
                "PropConductivity    , Conductivity;"
                "PropDensityPhi      , Density;"
                "PropCapacityPhi     , HeatCapacity;"
                "PropEigenStrainPhi  , EigenStrain" );

        //------------------------------------------------------------------------------
        // fill the stabilization parameter part of the parameter list


        // create SUPG/PSG fluid
        tParameterList( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex ).set( "stabilization_name", "SPIncFlow" );
        tParameterList( tSPIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tSPIndex ).set( "stabilization_type", fem::Stabilization_Type::INCOMPRESSIBLE_FLOW );
        tParameterList( tSPIndex ).set( "function_parameters", tCInv + "/0.0" );
        tParameterList( tSPIndex ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY;P", "Velocity,Pressure" ) );
        tParameterList( tSPIndex ).set( "leader_properties",
                "PropFluidDynViscosity,Viscosity;"
                "PropFluidDensity  ,Density" );

        // create SUPG fluid temperature advection
        tParameterList( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex ).set( "stabilization_name", "SPSUPGTemp" );
        tParameterList( tSPIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tSPIndex ).set( "stabilization_type", fem::Stabilization_Type::SUPG_ADVECTION );
        tParameterList( tSPIndex ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
        tParameterList( tSPIndex ).set( "leader_properties", "PropFluidConductivity,Conductivity" );

        // create SUPG Spalart-Allmaras model
        tParameterList( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex ).set( "stabilization_name", "SPSUPGSA" );
        tParameterList( tSPIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tSPIndex ).set( "stabilization_type", fem::Stabilization_Type::SUPG_SPALART_ALLMARAS_TURBULENCE );
        tParameterList( tSPIndex ).set( "leader_constitutive_models", "CMTurbulence,SpalartAllmarasTurbulence" );
        tParameterList( tSPIndex ).set( "function_parameters", tSupgTurbPower + "/" + tSupgTurbSource );
        tParameterList( tSPIndex ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY;VISCOSITY", "Velocity,Viscosity" ) );

        // create Nitsche for fluid velocity
        tParameterList( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex ).set( "stabilization_name", "SPNitscheU" );
        tParameterList( tSPIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tSPIndex ).set( "stabilization_type", fem::Stabilization_Type::VELOCITY_DIRICHLET_NITSCHE );
        tParameterList( tSPIndex ).set( "function_parameters", "100.0" );
        tParameterList( tSPIndex ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
        tParameterList( tSPIndex ).set( "leader_properties",
                "PropFluidDynViscosity,Viscosity;"
                "PropFluidDensity,Density" );

        // create Nitsche for fluid temperature
        tParameterList( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex ).set( "stabilization_name", "SPNitscheT" );
        tParameterList( tSPIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tSPIndex ).set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        tParameterList( tSPIndex ).set( "function_parameters", "100.0" );
        tParameterList( tSPIndex ).set( "leader_properties", "PropFluidConductivity,Material" );

        // create Nitsche for fluid turbulent viscosity
        tParameterList( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex ).set( "stabilization_name", "SPNitscheV" );
        tParameterList( tSPIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tSPIndex ).set( "stabilization_type", fem::Stabilization_Type::TURBULENCE_DIRICHLET_NITSCHE );
        tParameterList( tSPIndex ).set( "function_parameters", "100.0" );
        tParameterList( tSPIndex ).set( "leader_constitutive_models", "CMTurbulence,SpalartAllmarasTurbulence" );

        // create Nitsche for fluid/solid diffusion interface
        tParameterList( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex ).set( "stabilization_name", "SPFSInterfaceNitsche" );
        tParameterList( tSPIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tSPIndex ).set( "follower_phase_name", "PhaseSolid" );
        tParameterList( tSPIndex ).set( "stabilization_type", fem::Stabilization_Type::NITSCHE_INTERFACE );
        tParameterList( tSPIndex ).set( "function_parameters", "100.0" );
        tParameterList( tSPIndex ).set( "leader_properties", "PropFluidConductivity,Material" );
        tParameterList( tSPIndex ).set( "follower_properties", "PropSolidConductivity,Material" );

        // create ghost penalty viscous
        tParameterList( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex ).set( "stabilization_name", "SPGPViscous" );
        tParameterList( tSPIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tSPIndex ).set( "stabilization_type", fem::Stabilization_Type::VISCOUS_GHOST );
        tParameterList( tSPIndex ).set( "function_parameters", "0.05" );
        tParameterList( tSPIndex ).set( "leader_properties", "PropFluidDynViscosity,Viscosity" );

        // create ghost penalty convective
        tParameterList( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex ).set( "stabilization_name", "SPGPVelocity" );
        tParameterList( tSPIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tSPIndex ).set( "stabilization_type", fem::Stabilization_Type::CONVECTIVE_GHOST );
        tParameterList( tSPIndex ).set( "function_parameters", "0.05" );
        tParameterList( tSPIndex ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
        tParameterList( tSPIndex ).set( "leader_properties", "PropFluidDensity,Density" );

        // create ghost penalty pressure
        tParameterList( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex ).set( "stabilization_name", "SPGPPressure" );
        tParameterList( tSPIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tSPIndex ).set( "stabilization_type", fem::Stabilization_Type::PRESSURE_GHOST );
        tParameterList( tSPIndex ).set( "function_parameters", "0.005" );
        tParameterList( tSPIndex ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
        tParameterList( tSPIndex ).set( "leader_properties",
                "PropFluidDynViscosity,Viscosity;"
                "PropFluidDensity,Density" );

        // create ghost penalty fluid temperature
        tParameterList( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex ).set( "stabilization_name", "SPGPFluidTemp" );
        tParameterList( tSPIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tSPIndex ).set( "follower_phase_name", "PhaseFluid" );
        tParameterList( tSPIndex ).set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        tParameterList( tSPIndex ).set( "function_parameters", "0.05" );
        tParameterList( tSPIndex ).set( "leader_properties", "PropFluidConductivity,Material" );
        tParameterList( tSPIndex ).set( "follower_properties", "PropFluidConductivity,Material" );

        // create ghost fluid turbulence viscosity
        tParameterList( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex ).set( "stabilization_name", "SPGPViscosity" );
        tParameterList( tSPIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tSPIndex ).set( "follower_phase_name", "PhaseFluid" );
        tParameterList( tSPIndex ).set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        tParameterList( tSPIndex ).set( "function_parameters", "0.05" );
        tParameterList( tSPIndex ).set( "leader_properties", "PropFluidKinViscosity,Material" );

        // create ghost penalty solid temperature
        tParameterList( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex ).set( "stabilization_name", "SPGPSolidTemp" );
        tParameterList( tSPIndex ).set( "leader_phase_name", "PhaseSolid" );
        tParameterList( tSPIndex ).set( "follower_phase_name", "PhaseSolid" );
        tParameterList( tSPIndex ).set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        tParameterList( tSPIndex ).set( "function_parameters", "0.05" );
        tParameterList( tSPIndex ).set( "leader_properties", "PropSolidConductivity,Material" );
        tParameterList( tSPIndex ).set( "follower_properties", "PropSolidConductivity,Material" );

        // create ghost stabilization for theta and phi problems
        tParameterList( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex ).set( "stabilization_name", "SPGPFluidThetaPhi" );
        tParameterList( tSPIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tSPIndex ).set( "follower_phase_name", "PhaseFluid" );
        tParameterList( tSPIndex ).set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        tParameterList( tSPIndex ).set( "function_parameters", "0.01" );
        tParameterList( tSPIndex ).set( "leader_properties", "PropConductivity,Material" );

        // create DBC on interface for theta problem
        tParameterList( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex ).set( "stabilization_name", "SPNitscheFluidThetaPhi" );
        tParameterList( tSPIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tSPIndex ).set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        tParameterList( tSPIndex ).set( "function_parameters", "100.0" );
        tParameterList( tSPIndex ).set( "leader_properties", "PropConductivity,Material" );

        // create ghost stabilization parameter for theta and phi problems
        tParameterList( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex ).set( "stabilization_name", "SPGPInclusionThetaPhi" );
        tParameterList( tSPIndex ).set( "leader_phase_name", "PhaseSolid" );
        tParameterList( tSPIndex ).set( "follower_phase_name", "PhaseSolid" );
        tParameterList( tSPIndex ).set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        tParameterList( tSPIndex ).set( "function_parameters", "0.01" );
        tParameterList( tSPIndex ).set( "leader_properties", "PropConductivity,Material" );

        // create DBC on interface for theta problem
        tParameterList( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex ).set( "stabilization_name", "SPNitscheInclusionThetaPhi" );
        tParameterList( tSPIndex ).set( "leader_phase_name", "PhaseSolid" );
        tParameterList( tSPIndex ).set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        tParameterList( tSPIndex ).set( "function_parameters", "100.0" );
        tParameterList( tSPIndex ).set( "leader_properties", "PropConductivity,Material" );

        //------------------------------------------------------------------------------
        // fill the IWG part of the parameter list

        // fluid bulk IWGs -------------------------------------------------------------
        // NS incompressible (velocity)
        tParameterList( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex ).set( "IWG_name", "IWGVelocityBulk" );
        tParameterList( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex ).set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_BULK );
        tParameterList( tIWGIndex ).set( "dof_residual", "VX,VY" );
        tParameterList( tIWGIndex ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        tParameterList( tIWGIndex ).set( "stabilization_parameters", "SPIncFlow,IncompressibleFlow" );

        // NS incompressible (pressure)
        tParameterList( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex ).set( "IWG_name", "IWGPressureBulk" );
        tParameterList( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex ).set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_BULK );
        tParameterList( tIWGIndex ).set( "dof_residual", "P" );
        tParameterList( tIWGIndex ).set( "leader_properties", "PropFluidPressureSpring,PressureSpring" );
        tParameterList( tIWGIndex ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        tParameterList( tIWGIndex ).set( "stabilization_parameters", "SPIncFlow,IncompressibleFlow" );

        // diffusion
        tParameterList( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex ).set( "IWG_name", "IWGFluidDiffusionBulk" );
        tParameterList( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        tParameterList( tIWGIndex ).set( "dof_residual", "TEMP" );
        tParameterList( tIWGIndex ).set( "leader_constitutive_models", "CMFluidDiffusion,Diffusion" );

        // advection
        tParameterList( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex ).set( "IWG_name", "IWGFluidAdvectionBulk" );
        tParameterList( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex ).set( "IWG_type", fem::IWG_Type::ADVECTION_BULK );
        tParameterList( tIWGIndex ).set( "dof_residual", "TEMP" );
        tParameterList( tIWGIndex ).set( "leader_constitutive_models", "CMFluidDiffusion,Diffusion" );
        tParameterList( tIWGIndex ).set( "stabilization_parameters", "SPSUPGTemp,SUPG" );

        // turbulence
        tParameterList( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex ).set( "IWG_name", "IWGTurbulenceBulk" );
        tParameterList( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex ).set( "IWG_type", fem::IWG_Type::SPALART_ALLMARAS_TURBULENCE_BULK );
        tParameterList( tIWGIndex ).set( "dof_residual", "VISCOSITY" );
        tParameterList( tIWGIndex ).set( "leader_constitutive_models", "CMTurbulence,SpalartAllmarasTurbulence" );
        tParameterList( tIWGIndex ).set( "stabilization_parameters", "SPSUPGSA,SUPG" );

        // solid bulk IWGs -----------------------------------------------------------
        // diffusion
        tParameterList( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex ).set( "IWG_name", "IWGSolidDiffusionBulk" );
        tParameterList( tIWGIndex ).set( "leader_phase_name", "PhaseSolid" );
        tParameterList( tIWGIndex ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        tParameterList( tIWGIndex ).set( "dof_residual", "TEMP" );
        tParameterList( tIWGIndex ).set( "leader_properties", "PropVolumetricHeatFlux,Load" );
        tParameterList( tIWGIndex ).set( "leader_constitutive_models", "CMSolidDiffusion,Diffusion" );

        // fluid/solid ifc IWGs ------------------------------------------------------

        // temperature
        tParameterList( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex ).set( "IWG_name", "IWGInterfaceFluidSolid" );
        tParameterList( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
        tParameterList( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex ).set( "follower_phase_name", "PhaseSolid" );
        tParameterList( tIWGIndex ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE );
        tParameterList( tIWGIndex ).set( "dof_residual", "TEMP" );
        tParameterList( tIWGIndex ).set( "leader_constitutive_models", "CMFluidDiffusion,Diffusion" );
        tParameterList( tIWGIndex ).set( "follower_constitutive_models", "CMSolidDiffusion,Diffusion" );
        tParameterList( tIWGIndex ).set( "stabilization_parameters", "SPFSInterfaceNitsche,NitscheInterface" );

        // zero velocity
        tParameterList( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex ).set( "IWG_name", "IWGZeroVelocity" );
        tParameterList( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        tParameterList( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex ).set( "neighbor_phases", "PhaseSolid" );
        tParameterList( tIWGIndex ).set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_UNSYMMETRIC_NITSCHE );
        tParameterList( tIWGIndex ).set( "dof_residual", "VX,VY" );
        tParameterList( tIWGIndex ).set( "leader_properties", "PropZeroU,Dirichlet" );
        tParameterList( tIWGIndex ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        tParameterList( tIWGIndex ).set( "stabilization_parameters", "SPNitscheU,DirichletNitsche" );

        // zero velocity (pressure part)
        tParameterList( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex ).set( "IWG_name", "IWGZeroPressure" );
        tParameterList( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        tParameterList( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex ).set( "neighbor_phases", "PhaseSolid" );
        tParameterList( tIWGIndex ).set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_UNSYMMETRIC_NITSCHE );
        tParameterList( tIWGIndex ).set( "dof_residual", "P" );
        tParameterList( tIWGIndex ).set( "leader_properties", "PropZeroU,Dirichlet" );
        tParameterList( tIWGIndex ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );

        // zero viscosity
        tParameterList( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex ).set( "IWG_name", "IWGZeroViscosity" );
        tParameterList( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        tParameterList( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex ).set( "neighbor_phases", "PhaseSolid" );
        tParameterList( tIWGIndex ).set( "IWG_type", fem::IWG_Type::SPALART_ALLMARAS_TURBULENCE_DIRICHLET_UNSYMMETRIC_NITSCHE );
        tParameterList( tIWGIndex ).set( "dof_residual", "VISCOSITY" );
        tParameterList( tIWGIndex ).set( "leader_properties", "PropZeroV,Dirichlet" );
        tParameterList( tIWGIndex ).set( "leader_constitutive_models", "CMTurbulence,SpalartAllmarasTurbulence" );
        tParameterList( tIWGIndex ).set( "stabilization_parameters", "SPNitscheV,Nitsche" );

        // Laterial BC IWG (adiabatic) ----------------------------------------------------------------

        // zero velocity along lateral sides
        tParameterList( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex ).set( "IWG_name", "IWGZeroVelocityVoid" );
        tParameterList( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        tParameterList( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex ).set( "neighbor_phases", "PhaseVoidLateral" );
        tParameterList( tIWGIndex ).set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_UNSYMMETRIC_NITSCHE );
        tParameterList( tIWGIndex ).set( "dof_residual", "VX,VY" );
        tParameterList( tIWGIndex ).set( "leader_properties", "PropZeroU,Dirichlet" );
        tParameterList( tIWGIndex ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        tParameterList( tIWGIndex ).set( "stabilization_parameters", "SPNitscheU,DirichletNitsche" );

        // zero velocity along lateral sides (pressure part)
        tParameterList( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex ).set( "IWG_name", "IWGZeroPressureVoid" );
        tParameterList( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        tParameterList( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex ).set( "neighbor_phases", "PhaseVoidLateral" );
        tParameterList( tIWGIndex ).set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_UNSYMMETRIC_NITSCHE );
        tParameterList( tIWGIndex ).set( "dof_residual", "P" );
        tParameterList( tIWGIndex ).set( "leader_properties", "PropZeroU,Dirichlet" );
        tParameterList( tIWGIndex ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );

        // zero viscosity along lateral sides
        tParameterList( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex ).set( "IWG_name", "IWGZeroViscosityVoid" );
        tParameterList( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        tParameterList( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex ).set( "neighbor_phases", "PhaseVoidLateral" );
        tParameterList( tIWGIndex ).set( "IWG_type", fem::IWG_Type::SPALART_ALLMARAS_TURBULENCE_DIRICHLET_UNSYMMETRIC_NITSCHE );
        tParameterList( tIWGIndex ).set( "dof_residual", "VISCOSITY" );
        tParameterList( tIWGIndex ).set( "leader_properties", "PropZeroV,Dirichlet" );
        tParameterList( tIWGIndex ).set( "leader_constitutive_models", "CMTurbulence,SpalartAllmarasTurbulence" );
        tParameterList( tIWGIndex ).set( "stabilization_parameters", "SPNitscheV,Nitsche" );

        // Inlet BC IWG ----------------------------------------------------------------

        // inlet pressure
        tParameterList( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex ).set( "IWG_name", "IWGInletImposedPressure" );
        tParameterList( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex ).set( "neighbor_phases", "PhaseVoidFront" );
        tParameterList( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        tParameterList( tIWGIndex ).set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_IMPOSED_PRESSURE );
        tParameterList( tIWGIndex ).set( "dof_residual", "VX,VY" );
        tParameterList( tIWGIndex ).set( "leader_properties", "PropInletPressure,Pressure" );

        // inlet viscosity
        tParameterList( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex ).set( "IWG_name", "IWGInletViscosity" );
        tParameterList( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex ).set( "neighbor_phases", "PhaseVoidFront" );
        tParameterList( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        tParameterList( tIWGIndex ).set( "IWG_type", fem::IWG_Type::SPALART_ALLMARAS_TURBULENCE_DIRICHLET_SYMMETRIC_NITSCHE );
        tParameterList( tIWGIndex ).set( "dof_residual", "VISCOSITY" );
        tParameterList( tIWGIndex ).set( "leader_properties", "PropInletV,Dirichlet" );
        tParameterList( tIWGIndex ).set( "leader_constitutive_models", "CMTurbulence,SpalartAllmarasTurbulence" );
        tParameterList( tIWGIndex ).set( "stabilization_parameters", "SPNitscheV,Nitsche" );

        // inlet temperature
        tParameterList( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex ).set( "IWG_name", "IWGInletTemp" );
        tParameterList( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        tParameterList( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex ).set( "neighbor_phases", "PhaseVoidFront" );
        tParameterList( tIWGIndex ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        tParameterList( tIWGIndex ).set( "dof_residual", "TEMP" );
        tParameterList( tIWGIndex ).set( "leader_properties", "PropInletTemp,Dirichlet" );
        tParameterList( tIWGIndex ).set( "leader_constitutive_models", "CMFluidDiffusion,Diffusion" );
        tParameterList( tIWGIndex ).set( "stabilization_parameters", "SPNitscheT,DirichletNitsche" );

        // Ghost  ----------------------------------------------------------------

        if ( tUseGhost )
        {
            // ghost viscous
            tParameterList( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
            tParameterList( tIWGIndex ).set( "IWG_name", "IWGGPViscous" );
            tParameterList( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            tParameterList( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
            tParameterList( tIWGIndex ).set( "follower_phase_name", "PhaseFluid" );
            tParameterList( tIWGIndex ).set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            tParameterList( tIWGIndex ).set( "dof_residual", "VX,VY" );
            tParameterList( tIWGIndex ).set( "stabilization_parameters", "SPGPViscous,GhostSP" );
            tParameterList( tIWGIndex ).set( "ghost_order", (uint)tDispOrder );

            // ghost convective
            tParameterList( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
            tParameterList( tIWGIndex ).set( "IWG_name", "IWGGPConvective" );
            tParameterList( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            tParameterList( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
            tParameterList( tIWGIndex ).set( "follower_phase_name", "PhaseFluid" );
            tParameterList( tIWGIndex ).set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            tParameterList( tIWGIndex ).set( "dof_residual", "VX,VY" );
            tParameterList( tIWGIndex ).set( "stabilization_parameters", "SPGPVelocity,GhostSP" );
            tParameterList( tIWGIndex ).set( "ghost_order", (uint)tDispOrder );

            // ghost pressure
            tParameterList( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
            tParameterList( tIWGIndex ).set( "IWG_name", "IWGGPPressure" );
            tParameterList( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            tParameterList( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
            tParameterList( tIWGIndex ).set( "follower_phase_name", "PhaseFluid" );
            tParameterList( tIWGIndex ).set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            tParameterList( tIWGIndex ).set( "dof_residual", "P" );
            tParameterList( tIWGIndex ).set( "stabilization_parameters", "SPGPPressure,GhostSP" );
            tParameterList( tIWGIndex ).set( "ghost_order", (uint)tDispOrder );

            // ghost fluid temperature
            tParameterList( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
            tParameterList( tIWGIndex ).set( "IWG_name", "IWGGPFluidTemp" );
            tParameterList( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            tParameterList( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
            tParameterList( tIWGIndex ).set( "follower_phase_name", "PhaseFluid" );
            tParameterList( tIWGIndex ).set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            tParameterList( tIWGIndex ).set( "dof_residual", "TEMP" );
            tParameterList( tIWGIndex ).set( "stabilization_parameters", "SPGPFluidTemp,GhostSP" );
            tParameterList( tIWGIndex ).set( "ghost_order", (uint)tDispOrder );

            // ghost fluid viscosity
            tParameterList( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
            tParameterList( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            tParameterList( tIWGIndex ).set( "IWG_name", "IWGGPFluidViscosity" );
            tParameterList( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
            tParameterList( tIWGIndex ).set( "follower_phase_name", "PhaseFluid" );
            tParameterList( tIWGIndex ).set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            tParameterList( tIWGIndex ).set( "dof_residual", "VISCOSITY" );
            tParameterList( tIWGIndex ).set( "stabilization_parameters", "SPGPViscosity,GhostSP" );

            // ghost solid A temperature
            tParameterList( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
            tParameterList( tIWGIndex ).set( "IWG_name", "IWGGPSolidTemp" );
            tParameterList( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            tParameterList( tIWGIndex ).set( "leader_phase_name", "PhaseSolid" );
            tParameterList( tIWGIndex ).set( "follower_phase_name", "PhaseSolid" );
            tParameterList( tIWGIndex ).set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            tParameterList( tIWGIndex ).set( "dof_residual", "TEMP" );
            tParameterList( tIWGIndex ).set( "stabilization_parameters", "SPGPSolidTemp,GhostSP" );
            tParameterList( tIWGIndex ).set( "ghost_order", (uint)tDispOrder );
        }

        // theta problem  ----------------------------------------------------------------
        // theta bulk in fluid
        tParameterList( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex ).set( "IWG_name", "IWGFluidDiffusionThetaBulk" );
        tParameterList( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        tParameterList( tIWGIndex ).set( "dof_residual", "THETA" );
        tParameterList( tIWGIndex ).set( "leader_constitutive_models", "CMFluidDiffusionTheta,Diffusion" );

        // theta bulk in inclusion
        tParameterList( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex ).set( "IWG_name", "IWGInclusionSolidDiffusionThetaBulk" );
        tParameterList( tIWGIndex ).set( "leader_phase_name", "PhaseSolid" );
        tParameterList( tIWGIndex ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        tParameterList( tIWGIndex ).set( "dof_residual", "THETA" );
        tParameterList( tIWGIndex ).set( "leader_constitutive_models", "CMInclusionSolidDiffusionTheta,Diffusion" );

        // create parameter list for single side interface condition
        tParameterList( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex ).set( "IWG_name", "IWGSurfaceOuterTheta" );
        tParameterList( tIWGIndex ).set( "leader_phase_name", "PhaseSolid" );
        tParameterList( tIWGIndex ).set( "neighbor_phases", "PhaseFluid" );
        tParameterList( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        tParameterList( tIWGIndex ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        tParameterList( tIWGIndex ).set( "dof_residual", "THETA" );
        tParameterList( tIWGIndex ).set( "leader_properties", "PropPrescTheta,Dirichlet" );
        tParameterList( tIWGIndex ).set( "leader_constitutive_models", "CMInclusionSolidDiffusionTheta,Diffusion" );
        tParameterList( tIWGIndex ).set( "stabilization_parameters", "SPNitscheInclusionThetaPhi,DirichletNitsche" );

        // create parameter list for single side interface condition
        tParameterList( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex ).set( "IWG_name", "IWGSurfaceInnerTheta" );
        tParameterList( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex ).set( "neighbor_phases", "PhaseSolid,PhaseVoidLateral" );
        tParameterList( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        tParameterList( tIWGIndex ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        tParameterList( tIWGIndex ).set( "dof_residual", "THETA" );
        tParameterList( tIWGIndex ).set( "leader_properties", "PropPrescTheta,Dirichlet" );
        tParameterList( tIWGIndex ).set( "leader_constitutive_models", "CMFluidDiffusionTheta,Diffusion" );
        tParameterList( tIWGIndex ).set( "stabilization_parameters", "SPNitscheFluidThetaPhi,DirichletNitsche" );

        if ( tUseGhost )
        {
            // create IWG - ghost
            tParameterList( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
            tParameterList( tIWGIndex ).set( "IWG_name", "IWGGPInnerTheta" );
            tParameterList( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
            tParameterList( tIWGIndex ).set( "follower_phase_name", "PhaseFluid" );
            tParameterList( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            tParameterList( tIWGIndex ).set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            tParameterList( tIWGIndex ).set( "dof_residual", "THETA" );
            tParameterList( tIWGIndex ).set( "stabilization_parameters", "SPGPFluidThetaPhi,GhostSP" );

            // create IWG - ghost
            tParameterList( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
            tParameterList( tIWGIndex ).set( "IWG_name", "IWGGPOuterTheta" );
            tParameterList( tIWGIndex ).set( "leader_phase_name", "PhaseSolid" );
            tParameterList( tIWGIndex ).set( "follower_phase_name", "PhaseSolid" );
            tParameterList( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            tParameterList( tIWGIndex ).set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            tParameterList( tIWGIndex ).set( "dof_residual", "THETA" );
            tParameterList( tIWGIndex ).set( "stabilization_parameters", "SPGPInclusionThetaPhi,GhostSP" );
        }

        // create time side interface condition fluid
        tParameterList( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex ).set( "IWG_name", "IWGFluidTimeContinuityTheta" );
        tParameterList( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex ).set( "IWG_type", fem::IWG_Type::TIME_CONTINUITY_DOF );
        tParameterList( tIWGIndex ).set( "dof_residual", "THETA" );
        tParameterList( tIWGIndex ).set( "leader_properties",
                "PropWeightCurrent,       WeightCurrent;"
                "PropWeightPrevious,      WeightPrevious;"
                "PropInitialCondition,    InitialCondition" );

        // create time side interface condition solid
        tParameterList( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex ).set( "IWG_name", "IWGInclusionTimeContinuityTheta" );
        tParameterList( tIWGIndex ).set( "leader_phase_name", "PhaseSolid" );
        tParameterList( tIWGIndex ).set( "IWG_type", fem::IWG_Type::TIME_CONTINUITY_DOF );
        tParameterList( tIWGIndex ).set( "dof_residual", "THETA" );
        tParameterList( tIWGIndex ).set( "leader_properties",
                "PropWeightCurrent,       WeightCurrent;"
                "PropWeightPrevious,      WeightPrevious;"
                "PropInitialCondition,    InitialCondition" );

        // theta problem  ----------------------------------------------------------------
        // create IWG - bulk diffusion
        tParameterList( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex ).set( "IWG_name", "IWGDiffusionInnerBulk" );
        tParameterList( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        tParameterList( tIWGIndex ).set( "dof_residual", "PHID" );
        tParameterList( tIWGIndex ).set( "leader_constitutive_models", "CMFluidDiffusionPhi,Diffusion" );

        // create IWG - bulk diffusion
        tParameterList( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex ).set( "IWG_name", "IWGDiffusionOuterBulk" );
        tParameterList( tIWGIndex ).set( "leader_phase_name", "PhaseSolid" );
        tParameterList( tIWGIndex ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        tParameterList( tIWGIndex ).set( "dof_residual", "PHID" );
        tParameterList( tIWGIndex ).set( "leader_constitutive_models", "CMInclusionSolidDiffusionPhi,Diffusion" );

        // create parameter list for single side interface condition
        tParameterList( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex ).set( "IWG_name", "IWGSurfaceInnerPhi" );
        tParameterList( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex ).set( "neighbor_phases", "PhaseSolid,PhaseVoidLateral" );
        tParameterList( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        tParameterList( tIWGIndex ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        tParameterList( tIWGIndex ).set( "dof_residual", "PHID" );
        tParameterList( tIWGIndex ).set( "leader_properties", "PropPrescPhi,Dirichlet" );
        tParameterList( tIWGIndex ).set( "leader_constitutive_models", "CMFluidDiffusionPhi,Diffusion" );
        tParameterList( tIWGIndex ).set( "stabilization_parameters", "SPNitscheFluidThetaPhi,DirichletNitsche" );

        // create parameter list for single side interface condition
        tParameterList( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex ).set( "IWG_name", "IWGSurfaceOuterPhi" );
        tParameterList( tIWGIndex ).set( "leader_phase_name", "PhaseSolid" );
        tParameterList( tIWGIndex ).set( "neighbor_phases", "PhaseFluid" );
        tParameterList( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        tParameterList( tIWGIndex ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        tParameterList( tIWGIndex ).set( "dof_residual", "PHID" );
        tParameterList( tIWGIndex ).set( "leader_properties", "PropPrescPhi,Dirichlet" );
        tParameterList( tIWGIndex ).set( "leader_constitutive_models", "CMInclusionSolidDiffusionPhi,Diffusion" );
        tParameterList( tIWGIndex ).set( "stabilization_parameters", "SPNitscheInclusionThetaPhi,DirichletNitsche" );

        if ( tUseGhost )
        {
            // create IWG - ghost
            tParameterList( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
            tParameterList( tIWGIndex ).set( "IWG_name", "IWGGPInnerPhi" );
            tParameterList( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
            tParameterList( tIWGIndex ).set( "follower_phase_name", "PhaseFluid" );
            tParameterList( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            tParameterList( tIWGIndex ).set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            tParameterList( tIWGIndex ).set( "dof_residual", "PHID" );
            tParameterList( tIWGIndex ).set( "stabilization_parameters", "SPGPFluidThetaPhi,GhostSP" );

            // create IWG - ghost
            tParameterList( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
            tParameterList( tIWGIndex ).set( "IWG_name", "IWGGPOuterPhi" );
            tParameterList( tIWGIndex ).set( "leader_phase_name", "PhaseSolid" );
            tParameterList( tIWGIndex ).set( "follower_phase_name", "PhaseSolid" );
            tParameterList( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            tParameterList( tIWGIndex ).set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            tParameterList( tIWGIndex ).set( "dof_residual", "PHID" );
            tParameterList( tIWGIndex ).set( "stabilization_parameters", "SPGPInclusionThetaPhi,GhostSP" );
        }

        //------------------------------------------------------------------------------
        // fill the IQI part of the parameter list

        // velocity VX
        tParameterList( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex ).set( "IQI_name", "IQIBulkVX" );
        tParameterList( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIQIIndex ).set( "IQI_type", fem::IQI_Type::DOF );
        tParameterList( tIQIIndex ).set( "dof_quantity", "VX,VY" );
        tParameterList( tIQIIndex ).set( "vectorial_field_index", 0 );

        // velocity VY
        tParameterList( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex ).set( "IQI_name", "IQIBulkVY" );
        tParameterList( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIQIIndex ).set( "IQI_type", fem::IQI_Type::DOF );
        tParameterList( tIQIIndex ).set( "dof_quantity", "VX,VY" );
        tParameterList( tIQIIndex ).set( "vectorial_field_index", 1 );

        // pressure
        tParameterList( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex ).set( "IQI_name", "IQIBulkP" );
        tParameterList( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIQIIndex ).set( "IQI_type", fem::IQI_Type::DOF );
        tParameterList( tIQIIndex ).set( "dof_quantity", "P" );
        tParameterList( tIQIIndex ).set( "vectorial_field_index", 0 );

        // temperature
        tParameterList( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex ).set( "IQI_name", "IQIBulkTEMP" );
        tParameterList( tIQIIndex ).set( "leader_phase_name", "PhaseAll" );
        tParameterList( tIQIIndex ).set( "IQI_type", fem::IQI_Type::DOF );
        tParameterList( tIQIIndex ).set( "dof_quantity", "TEMP" );
        tParameterList( tIQIIndex ).set( "vectorial_field_index", 0 );

        // viscosity
        tParameterList( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex ).set( "IQI_name", "IQIBulkVISCOSITY" );
        tParameterList( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIQIIndex ).set( "IQI_type", fem::IQI_Type::DOF );
        tParameterList( tIQIIndex ).set( "dof_quantity", "VISCOSITY" );
        tParameterList( tIQIIndex ).set( "vectorial_field_index", 0 );

        // theta
        tParameterList( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex ).set( "IQI_name", "IQIBulkTHETA" );
        tParameterList( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIQIIndex ).set( "IQI_type", fem::IQI_Type::DOF );
        tParameterList( tIQIIndex ).set( "dof_quantity", "THETA" );
        tParameterList( tIQIIndex ).set( "vectorial_field_index", 0 );

        // phi
        tParameterList( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex ).set( "IQI_name", "IQIBulkPHID" );
        tParameterList( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIQIIndex ).set( "IQI_type", fem::IQI_Type::DOF );
        tParameterList( tIQIIndex ).set( "dof_quantity", "PHID" );
        tParameterList( tIQIIndex ).set( "vectorial_field_index", 0 );

        // fluid thermal energy on inlet
        tParameterList( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex ).set( "IQI_name", "IQIInletThermalEnergy" );
        tParameterList( tIQIIndex ).set( "IQI_bulk_type", fem::Element_Type::SIDESET );
        tParameterList( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIQIIndex ).set( "neighbor_phases", "PhaseVoidFront" );
        tParameterList( tIQIIndex ).set( "IQI_type", fem::IQI_Type::THERMAL_ENERGY_CONVECTIVE_FLUX );
        tParameterList( tIQIIndex ).set( "leader_constitutive_models", "CMFluidDiffusion,Diffusion" );

        // fluid thermal energy on outlet
        tParameterList( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex ).set( "IQI_name", "IQIOutletThermalEnergy" );
        tParameterList( tIQIIndex ).set( "IQI_bulk_type", fem::Element_Type::SIDESET );
        tParameterList( tIQIIndex ).set( "IQI_type", fem::IQI_Type::THERMAL_ENERGY_CONVECTIVE_FLUX );
        tParameterList( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIQIIndex ).set( "neighbor_phases", "PhaseVoidBack" );
        tParameterList( tIQIIndex ).set( "leader_constitutive_models", "CMFluidDiffusion,Diffusion" );

        // fluid total pressure on inlet
        tParameterList( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex ).set( "IQI_name", "IQIInletTotalPressure" );
        tParameterList( tIQIIndex ).set( "IQI_bulk_type", fem::Element_Type::SIDESET );
        tParameterList( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIQIIndex ).set( "neighbor_phases", "PhaseVoidFront" );
        tParameterList( tIQIIndex ).set( "IQI_type", fem::IQI_Type::TOTAL_PRESSURE );
        tParameterList( tIQIIndex ).set( "leader_constitutive_models", "CMFluid,Fluid" );

        // fluid total pressure on outlet
        tParameterList( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex ).set( "IQI_name", "IQIOutletTotalPressure" );
        tParameterList( tIQIIndex ).set( "IQI_bulk_type", fem::Element_Type::SIDESET );
        tParameterList( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIQIIndex ).set( "neighbor_phases", "PhaseVoidBack" );
        tParameterList( tIQIIndex ).set( "IQI_type", fem::IQI_Type::TOTAL_PRESSURE );
        tParameterList( tIQIIndex ).set( "leader_constitutive_models", "CMFluid,Fluid" );

        // fluid mass flow on inlet
        tParameterList( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex ).set( "IQI_name", "IQIInletMassFlow" );
        tParameterList( tIQIIndex ).set( "IQI_type", fem::IQI_Type::MASS_FLOW );
        tParameterList( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIQIIndex ).set( "neighbor_phases", "PhaseVoidFront" );
        tParameterList( tIQIIndex ).set( "IQI_bulk_type", fem::Element_Type::SIDESET );
        tParameterList( tIQIIndex ).set( "leader_constitutive_models", "CMFluid,Fluid" );

        // fluid mass flow on outlet
        tParameterList( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex ).set( "IQI_name", "IQIOutletMassFlow" );
        tParameterList( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIQIIndex ).set( "neighbor_phases", "PhaseVoidBack" );
        tParameterList( tIQIIndex ).set( "IQI_bulk_type", fem::Element_Type::SIDESET );
        tParameterList( tIQIIndex ).set( "IQI_type", fem::IQI_Type::MASS_FLOW );
        tParameterList( tIQIIndex ).set( "leader_constitutive_models", "CMFluid,Fluid" );

        // fluid power dissipation on inlet
        tParameterList( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex ).set( "IQI_name", "IQIInletPowDisp" );
        tParameterList( tIQIIndex ).set( "IQI_bulk_type", fem::Element_Type::SIDESET );
        tParameterList( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIQIIndex ).set( "neighbor_phases", "PhaseVoidFront" );
        tParameterList( tIQIIndex ).set( "IQI_type", fem::IQI_Type::POWER_DISSIPATION );
        tParameterList( tIQIIndex ).set( "leader_constitutive_models", "CMFluid,Fluid" );

        // fluid power dissipation on outlet
        tParameterList( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex ).set( "IQI_name", "IQIOutletPowDisp" );
        tParameterList( tIQIIndex ).set( "IQI_bulk_type", fem::Element_Type::SIDESET );
        tParameterList( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIQIIndex ).set( "neighbor_phases", "PhaseVoidBack" );
        tParameterList( tIQIIndex ).set( "IQI_type", fem::IQI_Type::POWER_DISSIPATION );
        tParameterList( tIQIIndex ).set( "leader_constitutive_models", "CMFluid,Fluid" );

        // fluid power dissipation in volume
        tParameterList( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex ).set( "IQI_name", "IQIVolumePowDisp" );
        tParameterList( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIQIIndex ).set( "IQI_type", fem::IQI_Type::POWER_DISSIPATION_BULK );
        tParameterList( tIQIIndex ).set( "leader_constitutive_models", "CMFluid,Fluid" );

        // inclusion perimeter
        tParameterList( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex ).set( "IQI_name", "IQIPerimeterItf" );
        tParameterList( tIQIIndex ).set( "leader_phase_name", "PhaseSolid" );
        tParameterList( tIQIIndex ).set( "neighbor_phases", "PhaseFluid" );
        tParameterList( tIQIIndex ).set( "IQI_bulk_type", fem::Element_Type::SIDESET );
        tParameterList( tIQIIndex ).set( "IQI_type", fem::IQI_Type::VOLUME );

        // max temperature in solid
        tParameterList( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex ).set( "IQI_name", "IQIMaxTemp" );
        tParameterList( tIQIIndex ).set( "leader_phase_name", "PhaseSolid" );
        tParameterList( tIQIIndex ).set( "IQI_type", fem::IQI_Type::MAX_DOF );
        tParameterList( tIQIIndex ).set( "function_parameters", tRefTemp + "/" + tMaxTempExponent );
        tParameterList( tIQIIndex ).set( "dof_quantity", "TEMP" );
        tParameterList( tIQIIndex ).set( "vectorial_field_index", 0 );

        // volume solid
        tParameterList( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex ).set( "IQI_name", "IQISolidVolume" );
        tParameterList( tIQIIndex ).set( "leader_phase_name", "PhaseSolid" );
        tParameterList( tIQIIndex ).set( "IQI_type", fem::IQI_Type::VOLUME );

        // wall distance
        tParameterList( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex ).set( "IQI_name", "IQIWallDistance" );
        tParameterList( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIQIIndex ).set( "IQI_type", fem::IQI_Type::PROPERTY );
        tParameterList( tIQIIndex ).set( "leader_properties", "PropWallDistance,Property" );

        // turbulent dynamic viscosity
        tParameterList( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex ).set( "IQI_name", "IQIBulkTurbDynVisc" );
        tParameterList( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIQIIndex ).set( "IQI_type", fem::IQI_Type::TURBULENT_DYNAMIC_VISCOSITY );
        tParameterList( tIQIIndex ).set( "leader_constitutive_models", "CMFluid,Fluid_Turbulence" );

        // effective dynamic viscosity
        tParameterList( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex ).set( "IQI_name", "IQIBulkEffDynVisc" );
        tParameterList( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIQIIndex ).set( "IQI_type", fem::IQI_Type::EFFECTIVE_DYNAMIC_VISCOSITY );
        tParameterList( tIQIIndex ).set( "leader_constitutive_models", "CMFluid,Fluid_Turbulence" );

        // effective conductivity
        tParameterList( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex ).set( "IQI_name", "IQIBulkEffCond" );
        tParameterList( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIQIIndex ).set( "IQI_type", fem::IQI_Type::EFFECTIVE_CONDUCTIVITY );
        tParameterList( tIQIIndex ).set( "leader_constitutive_models", "CMFluidDiffusion,Diffusion_Turbulence" );

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        tParameterList( tFEMIndex ).add_parameter_list( prm::create_computation_parameter_list() );
        tParameterList( tFEMIndex ).set( "print_physics_model", false );

        tParameterList( tFEMIndex ).set( "finite_difference_scheme", fem::FDScheme_Type::POINT_3_CENTRAL );
        tParameterList( tFEMIndex ).set( "finite_difference_perturbation_size", tFEMFdEpsilon );
    }

    void
    SOLParameterList( Vector< Submodule_Parameter_Lists >& tParameterlist )
    {
        tParameterlist.resize( 8 );

        tParameterlist( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL ) );

#ifdef MORIS_USE_MUMPS
        tParameterlist( 0 ).set( "Solver_Type", "Amesos_Mumps" );
#else
        tParameterlist( 0 ).set( "Solver_Type", "Amesos_Superludist" );
#endif

        //------------------------------------------------------------------------------


        tParameterlist( 1 ).add_parameter_list( moris::prm::create_linear_solver_parameter_list() );
        tParameterlist( 1 ).set( "DLA_Linear_solver_algorithms", "0" );

        //------------------------------------------------------------------------------


        tParameterlist( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );    // nonlinear algorithm index 0
        tParameterlist( 2 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        tParameterlist( 2 ).set( "NLA_Linear_solver", 0 );
        tParameterlist( 2 ).set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        tParameterlist( 2 ).set( "NLA_relaxation_strategy", sol::SolverRelaxationType::InvResNormAdaptive );
        tParameterlist( 2 ).set( "NLA_relaxation_parameter", 0.5 );
        tParameterlist( 2 ).set( "NLA_relaxation_damping", 0.5 );
        tParameterlist( 2 ).set( "NLA_max_iter", 100 );

        tParameterlist( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
        tParameterlist( 2 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        tParameterlist( 2 ).set( "NLA_Linear_solver", 0 );
        tParameterlist( 2 ).set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        tParameterlist( 2 ).set( "NLA_relaxation_parameter", 1.0 );
        tParameterlist( 2 ).set( "NLA_max_iter", 1 );

        tParameterlist( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
        tParameterlist( 2 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NLBGS_SOLVER );
        tParameterlist( 2 ).set( "NLA_Linear_solver", 0 );
        tParameterlist( 2 ).set( "NLA_rel_res_norm_drop", 1.0 );
        tParameterlist( 2 ).set( "NLA_max_iter", 1 );

        //------------------------------------------------------------------------------


        tParameterlist( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
        tParameterlist( 3 ).set( "NLA_Nonlinear_solver_algorithms", "1" );
        tParameterlist( 3 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        tParameterlist( 3 ).set( "NLA_DofTypes", "THETA" );

        tParameterlist( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
        tParameterlist( 3 ).set( "NLA_Nonlinear_solver_algorithms", "1" );
        tParameterlist( 3 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        tParameterlist( 3 ).set( "NLA_DofTypes", "PHID" );

        tParameterlist( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
        tParameterlist( 3 ).set( "NLA_Nonlinear_solver_algorithms", "0" );
        tParameterlist( 3 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        tParameterlist( 3 ).set( "NLA_DofTypes", "VX,VY,P,VISCOSITY" );

        tParameterlist( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
        tParameterlist( 3 ).set( "NLA_Nonlinear_solver_algorithms", "1" );
        tParameterlist( 3 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        tParameterlist( 3 ).set( "NLA_DofTypes", "TEMP" );

        tParameterlist( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
        tParameterlist( 3 ).set( "NLA_Nonlinear_solver_algorithms", "2" );
        tParameterlist( 3 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NLBGS_SOLVER );
        tParameterlist( 3 ).set( "NLA_Sub_Nonlinear_Solver", "0,1,2,3" );
        tParameterlist( 3 ).set( "NLA_DofTypes", "THETA;PHID;VX,VY,P,VISCOSITY;TEMP" );

        tParameterlist( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
        tParameterlist( 3 ).set( "NLA_Nonlinear_solver_algorithms", "1" );
        tParameterlist( 3 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        tParameterlist( 3 ).set( "NLA_DofTypes", "VX,VY,P,VISCOSITY" );

        tParameterlist( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
        tParameterlist( 3 ).set( "NLA_Nonlinear_solver_algorithms", "2" );
        tParameterlist( 3 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NLBGS_SOLVER );
        tParameterlist( 3 ).set( "NLA_Sub_Nonlinear_Solver", "0,1,5,3" );
        tParameterlist( 3 ).set( "NLA_DofTypes", "THETA;PHID;VX,VY,P,VISCOSITY;TEMP" );

        // ----------------------------------------------------------


        tParameterlist( 4 ).add_parameter_list( moris::prm::create_time_solver_algorithm_parameter_list() );
        tParameterlist( 4 ).set( "TSA_Nonlinear_Solver", 4 );                // using NLBGS for forward problem
        tParameterlist( 4 ).set( "TSA_Nonlinear_Sensitivity_Solver", 6 );    // using monlithic for sensitivity problem

        //------------------------------------------------------------------------------


        tParameterlist( 5 ).add_parameter_list( moris::prm::create_time_solver_parameter_list() );
        tParameterlist( 5 ).set( "TSA_DofTypes", "THETA;PHID;VX,VY,P,VISCOSITY;TEMP" );
        tParameterlist( 5 ).set( "TSA_Initialize_Sol_Vec", "THETA,0.0;PHID,0.0;VX,0.01;VY,0.0;P,0.0;VISCOSITY," + tInletKinViscosity + ";TEMP,0.0" );
        tParameterlist( 5 ).set( "TSA_Output_Indices", "0" );
        tParameterlist( 5 ).set( "TSA_Output_Criteria", "Output_Criterion" );

        //------------------------------------------------------------------------------


        tParameterlist( 6 ).add_parameter_list( moris::prm::create_solver_warehouse_parameterlist() );

        tParameterlist( 7 ).add_parameter_list( moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::NONE ) );
    }

    void
    MSIParameterList( Vector< Submodule_Parameter_Lists >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).add_parameter_list( prm::create_msi_parameter_list() );
        tParameterlist( 0 ).set( "THETA", 1 );
        tParameterlist( 0 ).set( "PHID", 1 );
        tParameterlist( 0 ).set( "VX", 1 );
        tParameterlist( 0 ).set( "VY", 1 );
        tParameterlist( 0 ).set( "P", 1 );
        tParameterlist( 0 ).set( "VISCOSITY", 1 );
        tParameterlist( 0 ).set( "TEMP", 1 );
    }

    void
    VISParameterList( Vector< Submodule_Parameter_Lists >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).add_parameter_list( prm::create_vis_parameter_list() );
        tParameterlist( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", tExoFile ) );
        tParameterlist( 0 ).set( "Mesh_Type", vis::VIS_Mesh_Type::STANDARD );
        tParameterlist( 0 ).set( "Set_Names", tAllPhases );
        tParameterlist( 0 ).set( "Field_Names", "THETA,PHID,VX,VY,P,VISCOSITY,TEMP,TURBDYNVISC,EFFVISC,EFFCOND" );
        tParameterlist( 0 ).set( "Field_Type", "NODAL,NODAL,NODAL,NODAL,NODAL,NODAL,NODAL,NODAL,NODAL,NODAL" );
        tParameterlist( 0 ).set( "IQI_Names",
                "IQIBulkTHETA,IQIBulkPHID,IQIBulkVX,IQIBulkVY,IQIBulkP,IQIBulkVISCOSITY,IQIBulkTEMP,"
                "IQIBulkTurbDynVisc,IQIBulkEffDynVisc,IQIBulkEffCond" );
        tParameterlist( 0 ).set( "Save_Frequency", 1 );
        tParameterlist( 0 ).set( "Time_Offset", 10.0 );
    }

    void
    MORISGENERALParameterList( Vector< Submodule_Parameter_Lists >& tParameterlist )
    {
    }

    //------------------------------------------------------------------------------
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
