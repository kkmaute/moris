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
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include <cmath>

// global variable for interpolation order
extern uint gInterpolationOrder;

extern uint gTestIndex;

#ifdef __cplusplus
extern "C" {
#endif
//------------------------------------------------------------------------------
namespace moris
{
    // funciton declarations forward
    void create_trilinos_solver_parameter_list( Module_Parameter_Lists& );
    void create_petsc_solver_parameter_list( Module_Parameter_Lists& );

    /* ------------------------------------------------------------------------ */
    // function to convert real value into string

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
    // problem type: forward only or optimization
    bool tIsOptimization = false;

    /* ------------------------------------------------------------------------ */
    // Interpolation order

    uint tInterpolationOrder = 2;

    uint tSolverType = 3;

    /* ------------------------------------------------------------------------ */
    // geometry parameters

    // dimension of domain
    real tDimX = 1.00;    // x-direction
    real tDimY = 1.00;    // y-direction

    // position and size of sphere
    real tPlanePos = 0.52;

    // approximate element size
    real tApproxEleSize = 1.0 / 4.0;

    /* ------------------------------------------------------------------------ */
    // loading parameters

    // precsribed temperature
    real tInTemp = 0.0;

    // prescribed heat flux - volumetric
    real tHeatFlx = 0.0;

    // prescribed heat flux - surface
    std::string tSurfaceHeatLoad = "0.0";

    /* ------------------------------------------------------------------------ */
    // material parameters

    std::string tDensity      = moris_to_string( 1.0 );
    std::string tCapacity     = moris_to_string( 1.0 );
    std::string tConductivity = moris_to_string( 1.916229656078740e-01 );

    /* ------------------------------------------------------------------------ */
    // parameters for weak enforcement of boundary conditions and ghost stabilization

    // Nitsche penalty
    std::string tNitschePenThermal = "10.0";

    // Flag to turn on/off ghost
    bool tUseGhost = false;

    /* ------------------------------------------------------------------------ */
    // parameters for transient simulations

    // number of time steps and total simulation time
    int         tTimeSteps = 0;
    moris::real tMaxTime   = 1.0;

    // time penalty parameter
    real tTimePenalty = 1.0;

    // time step size
    moris::real tDeltaTime = tMaxTime / (real)tTimeSteps;

    // number of steps to linearly increase heat laod
    moris::real tRampTime = 10;    // tMaxTime/10.0;

    // flag to turn on/off transient simulation
    bool tUseTimeContinuity = tTimeSteps > 1 ? true : false;

    /* ------------------------------------------------------------------------ */
    // Nonlinear solver parameters

    // maximum number of Newton steps
    int tNLA_max_iter = 1;

    // required drop of resdiual
    moris::real tNLA_rel_res_norm_drop = 1e-12;

    // relaxation parameter
    moris::real tNLA_relaxation_parameter = 1.0;

    /* ------------------------------------------------------------------------ */
    // File names

    std::string tName          = "Laplace";
    std::string tExoFile       = tName + ".exo";
    std::string tSoFile        = tName + ".so";
    std::string tHdf5File      = tName + "_SEN.hdf5";
    std::string tGENOutputFile = tName + "_GEN.exo";

    /* ------------------------------------------------------------------------ */
    //     // background mesh_sol

    // number of elements
    std::string tNumElemX = moris_to_string( std::ceil( tDimX / tApproxEleSize ) );
    std::string tNumElemY = moris_to_string( std::ceil( tDimY / tApproxEleSize ) );

    // dimension of background mesh
    std::string tDomainDimX = moris_to_string( tDimX );
    std::string tDomainDimY = moris_to_string( tDimY );

    // offeset of background mesh
    std::string tDomainOffX = moris_to_string( 0.0001 * ( -3.0 ) );
    std::string tDomainOffY = moris_to_string( 0.0 );

    // setting up information for HMR parameter list
    std::string tNumElemsPerDim = tNumElemX + "," + tNumElemY;
    std::string tDomainDims     = tDomainDimX + "," + tDomainDimY;
    std::string tDomainOffset   = tDomainOffX + "," + tDomainOffY;
    std::string tDomainSidesets = "1,2,3,4";

    int tLevelsetOrder = tInterpolationOrder;
    int tDispOrder     = tInterpolationOrder;

    int tLevelsetInitialRef = 0;
    int tDispInitialRef     = 0;
    int tRefineBuffer       = 1;

    std::string tLagrangeOrder   = std::to_string( std::max( tLevelsetOrder, tDispOrder ) );
    std::string tBsplineOrder    = std::to_string( tLevelsetOrder ) + "," + std::to_string( tDispOrder );
    std::string tInitialRef      = std::to_string( tLevelsetInitialRef ) + "," + std::to_string( tDispInitialRef );
    std::string tLagrangePattern = "0";

    uint tInterfaceRefinementSphere = 0;
    uint tInterfaceRefinementWalls  = 0;

    moris::real tElementEdgeLength = tApproxEleSize / ( std::pow( 2, tDispInitialRef ) );
    moris::real tGeoShift          = 0.1 * tElementEdgeLength;

    /* ------------------------------------------------------------------------ */
    // Optimization parameters

    // Bspline limit
    moris::real tBSplineLimit = 5.0 * tElementEdgeLength;

    // FD in adjoint
    real tFEMFdEpsilon = 1.0e-4;

    // FD in sweep
    std::string tSweepFdEpsilon = "1.0e-4";

    // number of constraints (here: number of design criteria)
    moris::uint tNumConstraints = 7;

    /* ------------------------------------------------------------------------ */
    // Output Config

    std::string tSolidPhase         = "HMR_dummy_n_p0,HMR_dummy_c_p0";
    std::string tVoidPhase          = "HMR_dummy_n_p1,HMR_dummy_c_p1";
    std::string tAllPhases          = tSolidPhase + "," + tVoidPhase;
    std::string tOuterSurfaceSolid  = "SideSet_1_n_p0,SideSet_2_n_p0,SideSet_3_n_p0,SideSet_4_n_p0,SideSet_1_c_p0,SideSet_2_c_p0,SideSet_3_c_p0,SideSet_4_c_p0";
    std::string tOuterSurfaceVoid   = "SideSet_1_n_p1,SideSet_2_n_p1,SideSet_3_n_p1,SideSet_4_n_p1,SideSet_1_c_p1,SideSet_2_c_p1,SideSet_3_c_p1,SideSet_4_c_p1";
    std::string tAllPhaseInterfaces = tAllPhases + "," + tOuterSurfaceSolid + "," + tOuterSurfaceVoid;

    /* ------------------------------------------------------------------------ */
    // geometry parameters & LS functions

    // sphere level set function
    moris::real
    Func_Sphere(
            const moris::Matrix< DDRMat >& aCoordinates,
            const Vector< real >&          aGeometryParameters )
    {
        // get coordinates
        real tX = aCoordinates( 0 );

        // compute level set value
        real tLS = tX - tPlanePos;

        // return the level set value
        return tLS;
    }

    /* ------------------------------------------------------------------------ */
    // geometry parameters & LS functions

    // plane function
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

    void
    Func_Anal_Sol(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        // intermeidate variables measuring the length of the solid domain
        real tA = tPlanePos;
        real tB = 1.0;

        aPropMatrix.set_size( 1, 1, 0.0 );

        // get coordinates
        real tX = aFIManager->get_IP_geometry_interpolator()->valx()( 0 );
        real tY = aFIManager->get_IP_geometry_interpolator()->valx()( 1 );

        // initailize the temprture
        real tTemp = 0;
        for ( int n = 1; n < 200; n++ )
        {
            // cast to real
            n = (real)n;

            // compute the coefficient
            real tCoeff = 2 / tB / std::cosh( n * M_PI * tA / tB ) * ( tB * tB * ( 2 * tB + ( -2 * tB + ( -1 + tB ) * n * n * M_PI * M_PI ) * std::cos( n * M_PI ) ) ) / ( std::pow( n * M_PI, 3.0 ) );

            tTemp += tCoeff * std::sin( n * M_PI * tY / tB ) * std::cosh( n * M_PI * ( tX ) / tB );
        }

        aPropMatrix( 0, 0 ) = tTemp;
    }

    /* ------------------------------------------------------------------------ */
    void
    Func_HeatLoad(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        // get the constant heat flux and set it as output
        aPropMatrix.set_size( 1, 1 );

        aPropMatrix( 0, 0 ) = tHeatFlx;
    }

    void
    Func_TempDistro(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        aPropMatrix.set_size( 1, 1 );

        const Matrix< DDRMat > tCoord = aFIManager->get_IP_geometry_interpolator()->valx();

        // the temp distrubution is y*(1-y) applied on the right side of the domain
        aPropMatrix( 0, 0 ) = tCoord( 1 ) * ( 1 - tCoord( 1 ) );
    }

    /* ------------------------------------------------------------------------ */
    // stored thermal energy

    void
    Func_StoredThermalEnergy(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        aPropMatrix.set_size( 1, 1 );

        auto tFITemp = aFIManager->get_field_interpolators_for_type( MSI::Dof_Type::TEMP );

        aPropMatrix( 0, 0 ) = std::stod( tDensity ) * std::stod( tCapacity ) * tFITemp->gradt( 1 )( 0 );
    }

    /* ------------------------------------------------------------------------ */
    // Phase assignement

    uint
    get_phase_index( const Bitset< 1 >& aGeometrySigns )
    {
        // by default solid
        uint tPhaseIndex = 0;

        // Phase void inside sphere
        if ( aGeometrySigns.test( 0 ) )
        {
            return 1;
        }

        // Phase
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

        return tDConstraintDCriteria;
    }

    /* ------------------------------------------------------------------------ */
    // function to dtermine when to oupt results

    bool
    Output_Criterion( moris::tsa::Time_Solver* aTimeSolver )
    {
        return true;
    }

    /* ------------------------------------------------------------------------ */

    void
    OPTParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists( 0 ).add_parameter_list( prm::create_opt_problem_parameter_list() );
        aParameterLists( 0 ).set( "is_optimization_problem", tIsOptimization );
        aParameterLists( 0 ).set( "problem", "user_defined" );
        aParameterLists( 0 ).set( "library", tSoFile );

        aParameterLists( 2 ).add_parameter_list( moris::prm::create_sweep_parameter_list() );
        aParameterLists( 2 ).set( "hdf5_path", tHdf5File );
        aParameterLists( 2 ).set( "num_evaluations_per_adv", "1" );
        aParameterLists( 2 ).set( "finite_difference_type", "all" );
        aParameterLists( 2 ).set( "finite_difference_epsilons", tSweepFdEpsilon );
    }

    /* ------------------------------------------------------------------------ */

    void
    HMRParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists( 0 ).add_parameter_list( prm::create_hmr_parameter_list() );

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

        aParameterLists( 0 ).set( "basis_function_vtk_file", "basisinhmr.vtk" );
        aParameterLists( 0 ).set( "write_lagrange_output_mesh_to_exodus", "lagrangehmr.exo" );
    }

    /* ------------------------------------------------------------------------ */

    void
    XTKParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_xtk_parameter_list() );
        aParameterLists( 0 ).set( "decompose", true );
        aParameterLists( 0 ).set( "decomposition_type", "conformal" );
        aParameterLists( 0 ).set( "enrich", true );
        aParameterLists( 0 ).set( "use_SPG_based_enrichment", true );
        aParameterLists( 0 ).set( "basis_rank", "bspline" );
        aParameterLists( 0 ).set( "enrich_mesh_indices", "0" );
        aParameterLists( 0 ).set( "ghost_stab", tUseGhost );
        aParameterLists( 0 ).set( "multigrid", false );
        aParameterLists( 0 ).set( "verbose", false );
        aParameterLists( 0 ).set( "high_to_low_dbl_side_sets", true );
        aParameterLists( 0 ).set( "print_enriched_ig_mesh", true );
        aParameterLists( 0 ).set( "exodus_output_XTK_ig_mesh", true );
        aParameterLists( 0 ).set( "exodus_output_XTK_ip_mesh", true );
        aParameterLists( 0 ).set( "activate_basis_agglomeration", true );
        aParameterLists( 0 ).set( "write_enrichment_fields", true );
        aParameterLists( 0 ).set( "write_enrichment_fields_probe_spheres", "1.0,0.5,0.5,1.0" );
    }

    /* ------------------------------------------------------------------------ */

    void
    GENParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists( 0 ).add_parameter_list( prm::create_gen_parameter_list() );
        aParameterLists( 0 ).set( "IQI_types",
                "IQIInputThermalEnergy",
                "IQIDiffusiveLower",
                "IQIDiffusiveUpper",
                "IQIDiffusiveFront",
                "IQIDiffusiveBack",
                "IQIStoredThermalEnergy",
                "IQIInputThermalEnergySurface" );

        aParameterLists( 0 ).set( "number_of_phases", 2 );
        aParameterLists( 0 ).set( "phase_function_name", "get_phase_index" );
        aParameterLists( 0 ).set( "output_mesh_file", tGENOutputFile );
        aParameterLists( 0 ).set( "time_offset", 10.0 );

        // Inclusions
        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists( 1 ).set( "field_function_name", "Func_Sphere" );
        aParameterLists( 1 ).set( "number_of_refinements", tInterfaceRefinementSphere );
        aParameterLists( 1 ).set( "refinement_mesh_index", 0 );
        aParameterLists( 1 ).set( "isocontour_threshold", 0.0 );
        aParameterLists( 1 ).set( "isocontour_tolerance", 1.0e-12 );
        aParameterLists( 1 ).set( "intersection_tolerance", 1.0e-12 );
        aParameterLists( 1 ).set( "use_multilinear_interpolation", true );
        aParameterLists( 1 ).set( "discretization_mesh_index", 0 );
        aParameterLists( 1 ).set( "discretization_lower_bound", -2.0 );
        aParameterLists( 1 ).set( "discretization_upper_bound", 2.0 );

    }

    /* ------------------------------------------------------------------------ */

    void
    FEMParameterList( Module_Parameter_Lists& aParameterLists )
    {

        // create a cell of cell of parameter list for fem
        uint tPropIndex  = 0;
        uint tCMIndex    = 1;
        uint tSPIndex    = 2;
        uint tIWGIndex   = 3;
        uint tIQIIndex   = 4;
        uint tFEMIndex   = 5;
        uint tPhaseIndex = 7;

        //------------------------------------------------------------------------------

        aParameterLists( tPhaseIndex ).add_parameter_list( prm::create_phase_parameter_list() );
        aParameterLists( tPhaseIndex ).set( "phase_name", "PhaseSolid" );
        aParameterLists( tPhaseIndex ).set( "phase_indices", "0" );

        aParameterLists( tPhaseIndex ).add_parameter_list( prm::create_phase_parameter_list() );
        aParameterLists( tPhaseIndex ).set( "phase_name", "PhaseVoid" );
        aParameterLists( tPhaseIndex ).set( "phase_indices", "1" );

        //------------------------------------------------------------------------------
        // fill the property part of the parameter list

        // soild properties ------------------------------------------------------------

        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropDensity" );
        aParameterLists( tPropIndex ).set( "function_parameters", tDensity );

        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropCapacity" );
        aParameterLists( tPropIndex ).set( "function_parameters", tCapacity );

        // create  conductivity property
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropConductivity" );
        aParameterLists( tPropIndex ).set( "function_parameters", tConductivity );

        // BC properties ---------------------------------------------------------------

        // create inlet temperature property
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropInletTemp" );
        aParameterLists( tPropIndex ).set( "function_parameters", moris_to_string( tInTemp ) );

        // create inlet temperature property
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropInletTemp2" );
        aParameterLists( tPropIndex ).set( "function_parameters", moris_to_string( 1.0 ) );
        aParameterLists( tPropIndex ).set( "value_function", "Func_TempDistro" );

        // create inlet temperature property
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropAnalSol" );
        aParameterLists( tPropIndex ).set( "function_parameters", moris_to_string( 1.0 ) );
        aParameterLists( tPropIndex ).set( "value_function", "Func_Anal_Sol" );

        // create heat load property
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropVolumetricHeatFlux" );
        aParameterLists( tPropIndex ).set( "value_function", "Func_HeatLoad" );

        // create heat load property
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropSurfaceHeatFlux" );
        aParameterLists( tPropIndex ).set( "function_parameters", tSurfaceHeatLoad );

        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropInitialTemp" );
        aParameterLists( tPropIndex ).set( "function_parameters", moris_to_string( tInTemp ) );

        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropTimeContinuity" );
        aParameterLists( tPropIndex ).set( "function_parameters", moris_to_string( tTimePenalty * std::stod( tDensity ) * std::stod( tCapacity ) ) );

        // Stored Thermal Energy
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropStoredThermalEnergy" );
        aParameterLists( tPropIndex ).set( "value_function", "Func_StoredThermalEnergy" );

        //------------------------------------------------------------------------------
        // fill the constitutive model part of the parameter list

        //  CM --------------------------------------------------------------------
        // create  diffusion CM
        aParameterLists( tCMIndex ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists( tCMIndex ).set( "constitutive_name", "CMDiffusion" );
        aParameterLists( tCMIndex ).set( "phase_name", "PhaseSolid" );
        aParameterLists( tCMIndex ).set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
        aParameterLists( tCMIndex ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists( tCMIndex ).set( "properties",
                "PropConductivity,Conductivity;"
                "PropDensity     ,Density;"
                "PropCapacity    ,HeatCapacity" );

        aParameterLists( tCMIndex ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists( tCMIndex ).set( "constitutive_name", "CMDiffusionL2" );
        aParameterLists( tCMIndex ).set( "phase_name", "PhaseSolid" );
        aParameterLists( tCMIndex ).set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
        aParameterLists( tCMIndex ).set( "dof_dependencies", std::pair< std::string, std::string >( "L2", "Temperature" ) );
        aParameterLists( tCMIndex ).set( "properties",
                "PropConductivity,Conductivity;"
                "PropDensity     ,Density;"
                "PropCapacity    ,HeatCapacity" );

        //------------------------------------------------------------------------------
        // fill the stabilization parameter part of the parameter list

        // create Nitsche for  temperature
        aParameterLists( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( tSPIndex ).set( "stabilization_name", "SPNitsche" );
        aParameterLists( tSPIndex ).set( "leader_phase_name", "PhaseSolid" );
        aParameterLists( tSPIndex ).set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        aParameterLists( tSPIndex ).set( "function_parameters", tNitschePenThermal );
        aParameterLists( tSPIndex ).set( "leader_properties", "PropConductivity,Material" );

        // create ghost penalty  temperature
        aParameterLists( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( tSPIndex ).set( "stabilization_name", "SPGPTemp" );
        aParameterLists( tSPIndex ).set( "leader_phase_name", "PhaseSolid" );
        aParameterLists( tSPIndex ).set( "follower_phase_name", "PhaseSolid" );
        aParameterLists( tSPIndex ).set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        aParameterLists( tSPIndex ).set( "function_parameters", "0.05" );
        aParameterLists( tSPIndex ).set( "leader_properties", "PropConductivity,Material" );

        //------------------------------------------------------------------------------
        // fill the IWG part of the parameter list

        //  bulk IWGs -------------------------------------------------------------
        // diffusion
        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGDiffusionBulk" );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseSolid" );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        aParameterLists( tIWGIndex ).set( "dof_residual", "TEMP" );
        aParameterLists( tIWGIndex ).set( "leader_properties", "PropVolumetricHeatFlux,Load" );
        aParameterLists( tIWGIndex ).set( "leader_constitutive_models", "CMDiffusion,Diffusion" );

        //--------------------------------------------------------------------------------
        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGDiffusionBulk" );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseSolid" );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        aParameterLists( tIWGIndex ).set( "dof_residual", "L2" );
        aParameterLists( tIWGIndex ).set( "leader_properties", "PropVolumetricHeatFlux,Load" );
        aParameterLists( tIWGIndex ).set( "leader_constitutive_models", "CMDiffusionL2,Diffusion" );

        // Inlet BC IWG ----------------------------------------------------------------

        //
        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGInletTempTopBottom" );
        aParameterLists( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseSolid" );
        aParameterLists( tIWGIndex ).set( "side_ordinals", "1,3" );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists( tIWGIndex ).set( "dof_residual", "TEMP" );
        aParameterLists( tIWGIndex ).set( "leader_properties", "PropInletTemp,Dirichlet" );
        aParameterLists( tIWGIndex ).set( "leader_constitutive_models", "CMDiffusion,Diffusion" );
        aParameterLists( tIWGIndex ).set( "stabilization_parameters", "SPNitsche,DirichletNitsche" );

        //--------------------------------------------------------------------------------

        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGInletTempRight" );
        aParameterLists( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseSolid" );
        aParameterLists( tIWGIndex ).set( "neighbor_phases", "PhaseVoid" );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists( tIWGIndex ).set( "dof_residual", "TEMP" );
        aParameterLists( tIWGIndex ).set( "leader_properties", "PropInletTemp2,Dirichlet" );
        aParameterLists( tIWGIndex ).set( "leader_constitutive_models", "CMDiffusion,Diffusion" );
        aParameterLists( tIWGIndex ).set( "stabilization_parameters", "SPNitsche,DirichletNitsche" );

        // Inlet BC IWG ----------------------------------------------------------------

        //
        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGInletTempTopBottom" );
        aParameterLists( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseSolid" );
        aParameterLists( tIWGIndex ).set( "side_ordinals", "1,3" );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists( tIWGIndex ).set( "dof_residual", "L2" );
        aParameterLists( tIWGIndex ).set( "leader_properties", "PropInletTemp,Dirichlet" );
        aParameterLists( tIWGIndex ).set( "leader_constitutive_models", "CMDiffusionL2,Diffusion" );
        aParameterLists( tIWGIndex ).set( "stabilization_parameters", "SPNitsche,DirichletNitsche" );

        //--------------------------------------------------------------------------------

        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGInletTempRight" );
        aParameterLists( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseSolid" );
        aParameterLists( tIWGIndex ).set( "neighbor_phases", "PhaseVoid" );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists( tIWGIndex ).set( "dof_residual", "L2" );
        aParameterLists( tIWGIndex ).set( "leader_properties", "PropInletTemp2,Dirichlet" );
        aParameterLists( tIWGIndex ).set( "leader_constitutive_models", "CMDiffusionL2,Diffusion" );
        aParameterLists( tIWGIndex ).set( "stabilization_parameters", "SPNitsche,DirichletNitsche" );
        // Outlet BC IWG ----------------------------------------------------------------

        // Time continuity ----------------------------------------------------------------

        if ( tUseTimeContinuity )
        {
            // Time continuity temperature
            aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists( tIWGIndex ).set( "IWG_name", "IWGTimeContinuity" );
            aParameterLists( tIWGIndex ).set( "IWG_type", ( fem::IWG_Type::TIME_CONTINUITY_DOF ) );
            aParameterLists( tIWGIndex ).set( "dof_residual", "TEMP" );
            aParameterLists( tIWGIndex ).set( "leader_properties",
                    "PropTimeContinuity,WeightCurrent;"
                    "PropTimeContinuity,WeightPrevious;"
                    "PropInitialTemp,InitialCondition" );
            aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseSolid" );
            aParameterLists( tIWGIndex ).set( "time_continuity", true );
            }

        // Ghost  ----------------------------------------------------------------

        if ( tUseGhost )
        {
            // ghost  temperature
            aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists( tIWGIndex ).set( "IWG_name", "IWGGPTemp" );
            aParameterLists( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseSolid" );
            aParameterLists( tIWGIndex ).set( "follower_phase_name", "PhaseSolid" );
            aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists( tIWGIndex ).set( "dof_residual", "TEMP" );
            aParameterLists( tIWGIndex ).set( "stabilization_parameters", "SPGPTemp,GhostSP" );
            aParameterLists( tIWGIndex ).set( "ghost_order", (uint)tDispOrder );

            aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists( tIWGIndex ).set( "IWG_name", "IWGGPTemp" );
            aParameterLists( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseSolid" );
            aParameterLists( tIWGIndex ).set( "follower_phase_name", "PhaseSolid" );
            aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists( tIWGIndex ).set( "dof_residual", "L2" );
            aParameterLists( tIWGIndex ).set( "stabilization_parameters", "SPGPTemp,GhostSP" );
            aParameterLists( tIWGIndex ).set( "ghost_order", (uint)tDispOrder );
        }

        //------------------------------------------------------------------------------
        // fill the IQI part of the parameter list

        // input thermal energy
        aParameterLists( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIAnalTemp" );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseSolid" );
        aParameterLists( tIQIIndex ).set( "IQI_type", ( fem::IQI_Type::PROPERTY ) );
        aParameterLists( tIQIIndex ).set( "leader_properties", "PropAnalSol,Property" );
        aParameterLists( tIQIIndex ).set( "vectorial_field_index", 0 );

        // temperature
        aParameterLists( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIBulkTEMP" );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseSolid" );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists( tIQIIndex ).set( "dof_quantity", "TEMP" );
        aParameterLists( tIQIIndex ).set( "vectorial_field_index", 0 );

        aParameterLists( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIBulkL2" );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseSolid" );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists( tIQIIndex ).set( "dof_quantity", "L2" );
        aParameterLists( tIQIIndex ).set( "vectorial_field_index", 0 );

        // input thermal energy
        aParameterLists( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIInputThermalEnergy" );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseSolid" );
        aParameterLists( tIQIIndex ).set( "IQI_type", ( fem::IQI_Type::PROPERTY ) );
        aParameterLists( tIQIIndex ).set( "leader_properties", "PropVolumetricHeatFlux,Property" );

        aParameterLists( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIBulkL2Error" );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::L2_ERROR_ANALYTIC );
        aParameterLists( tIQIIndex ).set( "dof_quantity", "TEMP" );
        aParameterLists( tIQIIndex ).set( "leader_properties", "PropAnalSol,L2Check" );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseSolid" );

        // input thermal energy
        aParameterLists( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIInputThermalEnergySurface" );
        aParameterLists( tIQIIndex ).set( "IQI_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseSolid" );
        aParameterLists( tIQIIndex ).set( "side_ordinals", "2" );
        aParameterLists( tIQIIndex ).set( "IQI_type", ( fem::IQI_Type::PROPERTY ) );
        aParameterLists( tIQIIndex ).set( "leader_properties", "PropSurfaceHeatFlux,Property" );

        //  diffusive flux
        aParameterLists( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIDiffusiveLower" );
        aParameterLists( tIQIIndex ).set( "IQI_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseSolid" );
        aParameterLists( tIQIIndex ).set( "side_ordinals", "1" );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::THERMAL_ENERGY_DIFFUSIVE_FLUX );
        aParameterLists( tIQIIndex ).set( "leader_constitutive_models", "CMDiffusion,Diffusion" );

        //  diffusive flux
        aParameterLists( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIDiffusiveUpper" );
        aParameterLists( tIQIIndex ).set( "IQI_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseSolid" );
        aParameterLists( tIQIIndex ).set( "side_ordinals", "3" );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::THERMAL_ENERGY_DIFFUSIVE_FLUX );
        aParameterLists( tIQIIndex ).set( "leader_constitutive_models", "CMDiffusion,Diffusion" );

        //  diffusive flux
        aParameterLists( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIDiffusiveFront" );
        aParameterLists( tIQIIndex ).set( "IQI_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseSolid" );
        aParameterLists( tIQIIndex ).set( "side_ordinals", "4" );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::THERMAL_ENERGY_DIFFUSIVE_FLUX );
        aParameterLists( tIQIIndex ).set( "leader_constitutive_models", "CMDiffusion,Diffusion" );

        //  diffusive flux
        aParameterLists( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIDiffusiveBack" );
        aParameterLists( tIQIIndex ).set( "IQI_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseSolid" );
        aParameterLists( tIQIIndex ).set( "side_ordinals", "2" );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::THERMAL_ENERGY_DIFFUSIVE_FLUX );
        aParameterLists( tIQIIndex ).set( "leader_constitutive_models", "CMDiffusion,Diffusion" );

        // stored thermal energy
        aParameterLists( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIStoredThermalEnergy" );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseSolid" );
        aParameterLists( tIQIIndex ).set( "IQI_type", ( fem::IQI_Type::PROPERTY ) );
        aParameterLists( tIQIIndex ).set( "leader_properties", "PropStoredThermalEnergy,Property" );

        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkEigen" );
        aParameterLists( 4 ).set( "IQI_type", ( fem::IQI_Type::EIGEN_VECTOR ) );
        aParameterLists( 4 ).set( "function_parameters", "0" );
        aParameterLists( 4 ).set( "dof_quantity", "TEMP" );
        aParameterLists( 4 ).set( "leader_phase_name", "PhaseSolid" );
        aParameterLists( 4 ).set( "vectorial_field_index", 0 );

        if ( gTestIndex == 1 )
        {
            aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
            aParameterLists( 4 ).set( "IQI_name", "IQIBulkEigenVal" );
            aParameterLists( 4 ).set( "IQI_type", ( fem::IQI_Type::EIGEN_VALUE ) );
            aParameterLists( 4 ).set( "function_parameters", "0" );
            aParameterLists( 4 ).set( "dof_quantity", "TEMP" );
            aParameterLists( 4 ).set( "leader_phase_name", "PhaseSolid" );
            aParameterLists( 4 ).set( "vectorial_field_index", 0 );
            }

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        aParameterLists( tFEMIndex ).add_parameter_list( prm::create_computation_parameter_list() );
        aParameterLists( tFEMIndex ).set( "print_physics_model", false );

        aParameterLists( tFEMIndex ).set( "finite_difference_scheme", fem::FDScheme_Type::POINT_3_CENTRAL );
        aParameterLists( tFEMIndex ).set( "finite_difference_perturbation_size", tFEMFdEpsilon );
    }

    void
    SOLParameterList( Module_Parameter_Lists& aParameterLists )
    {
        std::cout << "gTestIndex_" + std::to_string( par_rank() ) + ": " << gTestIndex << '\n';
        gTestIndex == 0 ? create_trilinos_solver_parameter_list( aParameterLists ) : create_petsc_solver_parameter_list( aParameterLists );
    }

    void
    MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_msi_parameter_list() );
        aParameterLists( 0 ).set( "TEMP", 0 );
        aParameterLists( 0 ).set( "number_eigen_vectors", gTestIndex == 0 ? 5 : 1 );
        aParameterLists( 0 ).set( "order_adofs_by_host", true );
    }

    void
    VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_vis_parameter_list() );
        aParameterLists( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", tExoFile ) );
        aParameterLists( 0 ).set( "Mesh_Type", vis::VIS_Mesh_Type::STANDARD );
        aParameterLists( 0 ).set( "Set_Names", tAllPhaseInterfaces );
        aParameterLists( 0 ).set( "Field_Names", "TEMP,HEAT,TEMP_A,L2Nodal,L2Glob,DiffLower,EigenVec,EigenVal,TEMPL" );
        aParameterLists( 0 ).set( "Field_Type", "NODAL,NODAL,NODAL,NODAL,GLOBAL,GLOBAL,NODAL,NODAL,NODAL" );
        aParameterLists( 0 ).set( "IQI_Names", "IQIBulkTEMP,IQIInputThermalEnergy,IQIAnalTemp,IQIBulkL2Error,IQIBulkL2Error,IQIDiffusiveLower,IQIBulkEigen,IQIBulkEigenVal,IQIBulkL2" );
        aParameterLists( 0 ).set( "Save_Frequency", 1 );
        aParameterLists( 0 ).set( "Time_Offset", 10.0 );
    }

    void
    MORISGENERALParameterList( Module_Parameter_Lists& aParameterLists )
    {
    }

    //------------------------------------------------------------------------------

    void
    create_petsc_solver_parameter_list( Module_Parameter_Lists& aParameterLists )
    {

        // 5 linear solvers, 1 linear solver to find the displacement,
        // 2 slpec solver and the associated linear solver object
        aParameterLists( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::PETSC ) );
        aParameterLists( 0 ).set( "KSPType", "fgmres" );
        aParameterLists( 0 ).set( "preconditioners", "0" );    // 10 shift_invert

        aParameterLists( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::PETSC ) );
        aParameterLists( 0 ).set( "KSPType", "preonly" );
        aParameterLists( 0 ).set( "preconditioners", "0" );    // 10 shift_invert

        aParameterLists( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::PETSC ) );
        aParameterLists( 0 ).set( "KSPType", "preonly" );
        aParameterLists( 0 ).set( "preconditioners", "1" );    // 10 shift_invert

        // find max eigen value
        aParameterLists( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::SLEPC_SOLVER ) );
        aParameterLists( 0 ).set( "Eigen_Algorithm", "power" );
        aParameterLists( 0 ).set( "Which", std::string( "LM" ) );
        aParameterLists( 0 ).set( "Num_Eig_Vals", 1 );
        aParameterLists( 0 ).set( "sub_linear_solver", "1" );    // 10 shift_invert
        aParameterLists( 0 ).set( "is_symmetric", false );       // 10 shift_invert
        aParameterLists( 0 ).set( "Update_Flag", true );         // 10 shift_invert

        // find min eigen value
        aParameterLists( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::SLEPC_SOLVER ) );
        aParameterLists( 0 ).set( "Eigen_Algorithm", "power" );
        aParameterLists( 0 ).set( "Which", std::string( "LM" ) );
        aParameterLists( 0 ).set( "Num_Eig_Vals", 1 );
        aParameterLists( 0 ).set( "STType", "shift_invert" );
        aParameterLists( 0 ).set( "sub_linear_solver", "2" );    // 10 shift_invert
        aParameterLists( 0 ).set( "is_symmetric", false );       // 10 shift_invert
        aParameterLists( 0 ).set( "Update_Flag", false );        // 10 shift_invert

        // precondioerr
        aParameterLists( 7 ).add_parameter_list( moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::PETSC ) );
        aParameterLists( 7 ).set( "PCType", "none" );

        // Ifpack precondioner for the eigen solve
        aParameterLists( 7 ).add_parameter_list( moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::PETSC ) );
        aParameterLists( 7 ).set( "PCType", "mumps" );

        // **
        aParameterLists( 1 ).add_parameter_list( moris::prm::create_linear_solver_parameter_list() );    // linear problem index 0
        aParameterLists( 1 ).set( "DLA_Linear_solver_algorithms", "0" );

        aParameterLists( 1 ).add_parameter_list( moris::prm::create_linear_solver_parameter_list() );    // eigen problem index 1
        aParameterLists( 1 ).set( "DLA_Linear_solver_algorithms", "3" );
        aParameterLists( 1 ).set( "RHS_Matrix_Type", "IdentityMat" );

        aParameterLists( 1 ).add_parameter_list( moris::prm::create_linear_solver_parameter_list() );    // eigen problem index 1
        aParameterLists( 1 ).set( "DLA_Linear_solver_algorithms", "4" );
        aParameterLists( 1 ).set( "RHS_Matrix_Type", "IdentityMat" );

        //------------------------------------------------------------------------------

        aParameterLists( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );    // nonlinear algorithm index 0
        aParameterLists( 2 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists( 2 ).set( "NLA_Linear_solver", 0 );
        aParameterLists( 2 ).set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        aParameterLists( 2 ).set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
        aParameterLists( 2 ).set( "NLA_max_iter", tNLA_max_iter );

        aParameterLists( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
        aParameterLists( 2 ).set( "NLA_Linear_solver", 1 );
        aParameterLists( 2 ).set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        aParameterLists( 2 ).set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
        aParameterLists( 2 ).set( "NLA_max_iter", tNLA_max_iter );
        aParameterLists( 2 ).set( "NLA_is_eigen_problem", true );

        aParameterLists( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
        aParameterLists( 2 ).set( "NLA_Linear_solver", 2 );
        aParameterLists( 2 ).set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        aParameterLists( 2 ).set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
        aParameterLists( 2 ).set( "NLA_max_iter", tNLA_max_iter );
        aParameterLists( 2 ).set( "NLA_is_eigen_problem", true );

        aParameterLists( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );    // nonlinear algorithm index 1
        aParameterLists( 2 ).set( "NLA_Solver_Implementation", static_cast< uint >( moris::NLA::NonlinearSolverType::NLBGS_SOLVER ) );
        aParameterLists( 2 ).set( "NLA_rel_res_norm_drop", 1.0e-9 );
        aParameterLists( 2 ).set( "NLA_max_iter", 1 );

        //------------------------------------------------------------------------------

        aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );    // 1: thermal subproblem
        aParameterLists( 3 ).set( "NLA_Nonlinear_solver_algorithms", "0,1" );           // set nonlinear algorithm with index 0
        aParameterLists( 3 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists( 3 ).set( "NLA_DofTypes", "TEMP" );

        aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );    // 1: thermal subproblem
        aParameterLists( 3 ).set( "NLA_Nonlinear_solver_algorithms", "0" );             // set nonlinear algorithm with index 0
        aParameterLists( 3 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists( 3 ).set( "NLA_DofTypes", "L2" );

        aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );    // nonlinear solver index 2
        aParameterLists( 3 ).set( "NLA_Solver_Implementation", static_cast< uint >( moris::NLA::NonlinearSolverType::NLBGS_SOLVER ) );
        aParameterLists( 3 ).set( "NLA_Nonlinear_solver_algorithms", "3" );    // set nonlinear algorithm with index 1.
        aParameterLists( 3 ).set( "NLA_Sub_Nonlinear_Solver", "0,1" );         // set sub nonlinear solvers with index 0 and 1
        aParameterLists( 3 ).set( "NLA_DofTypes", "TEMP;L2" );

        // ----------------------------------------------------------

        aParameterLists( 4 ).add_parameter_list( moris::prm::create_time_solver_algorithm_parameter_list() );
        aParameterLists( 4 ).set( "TSA_Nonlinear_Solver", 2 );                // using NLBGS for forward problem
        aParameterLists( 4 ).set( "TSA_Nonlinear_Sensitivity_Solver", 2 );    // using monlithic for sensitivity problem

        if ( tUseTimeContinuity )
        {
            aParameterLists( 4 ).set( "TSA_Num_Time_Steps", tTimeSteps );
            aParameterLists( 4 ).set( "TSA_Time_Frame", tMaxTime );
        }

        //------------------------------------------------------------------------------

        aParameterLists( 5 ).add_parameter_list( moris::prm::create_time_solver_parameter_list() );
        aParameterLists( 5 ).set( "TSA_DofTypes", "TEMP;L2" );
        // aParameterLists( 5 ).set( "TSA_Initialize_Sol_Vec", "TEMP,0.0" );
        aParameterLists( 5 ).set( "TSA_Output_Indices", "0" );
        aParameterLists( 5 ).set( "TSA_Output_Criteria", "Output_Criterion" );

        if ( tUseTimeContinuity )
        {
            aParameterLists( 5 ).set( "TSA_time_level_per_type", "TEMP,2" );
        }
        else
        {
            aParameterLists( 5 ).set( "TSA_time_level_per_type", "TEMP,1" );
        }

        //------------------------------------------------------------------------------

        aParameterLists( 6 ).add_parameter_list( moris::prm::create_solver_warehouse_parameterlist() );
        aParameterLists( 6 ).set( "SOL_save_operator_to_matlab", "Jacobian_petsc" );
        aParameterLists( 6 ).set( "SOL_TPL_Type", static_cast< uint >( sol::MapType::Petsc ) );
    }

    //------------------------------------------------------------------------------
    void
    create_trilinos_solver_parameter_list( Module_Parameter_Lists& aParameterLists )
    {


        sint tVerbosity = Belos::Errors + Belos::Warnings + Belos::IterationDetails + Belos::TimingDetails + Belos::StatusTestDetails + Belos::FinalSummary;

        switch ( tSolverType )
        {
            case 0:    // Amesos - Pardiso
                aParameterLists( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL ) );

                //                 aParameterLists( 0 ).set( "Solver_Type", "Amesos_Mumps");
                break;

            case 1:    // Amesos - Mumps
                aParameterLists( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL ) );

                aParameterLists( 0 ).set( "Solver_Type", "Amesos_Mumps" );
                break;

            case 2:    // Belos
                aParameterLists( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::BELOS_IMPL ) );

                // Solver type: GMRES, Flexible GMRES, Block CG , PseudoBlockCG, Stochastic CG, Recycling GMRES, Recycling CG, MINRES, LSQR, TFQMR
                //              Pseudoblock TFQMR, Seed GMRES, Seed CG
                aParameterLists( 0 ).set( "Solver Type", "Flexible GMRES" );

                // Diagnostics: Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::StatusTestDetails
                aParameterLists( 0 ).set( "Verbosity", tVerbosity );

                // Maximum number of blocks in Krylov factorization
                aParameterLists( 0 ).set( "Num Blocks", 250 );

                // Block size to be used by iterative solver
                aParameterLists( 0 ).set( "Block Size", 1 );

                // Allowable Belos solver iterations
                aParameterLists( 0 ).set( "Maximum Iterations", 500 );

                // Allowable Belos solver iterations
                aParameterLists( 0 ).set( "Maximum Restarts", 2 );

                // Convergence criteria
                aParameterLists( 0 ).set( "Convergence Tolerance", 1e-9 );

                // Left or right preconditioner
                aParameterLists( 0 ).set( "Left-right Preconditioner", "right" );

                // Preconditioner

                // ifpack - ILU
                // aParameterLists( 0 ).set( "ifpack_prec_type",  "ILU");
                // aParameterLists( 0 ).set( "fact: level-of-fill",  5 );

                // ifpack - ILUT
                aParameterLists( 0 ).set( "ifpack_prec_type", "ILUT" );
                aParameterLists( 0 ).set( "fact: ilut level-of-fill", 5.0 );
                aParameterLists( 0 ).set( "fact: drop tolerance", 1e-12 );

                // ifpack with direct solve
                // aParameterLists( 0 ).set( "ifpack_prec_type",  "amesos");
                // aParameterLists( 0 ).set( "amesos: solver type", "Amesos_Pardiso");

                // AMG with defaults for non-symmetric system
                // aParameterLists( 0 ).set( "ml_prec_type",  "NSSA");
                // aParameterLists( 0 ).set( "PDE equations", 3);

                break;

            case 3:    // AZTEC
                aParameterLists( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AZTEC_IMPL ) );

                // options are: AZ_gmres, AZ_gmres_condnum, AZ_cg, AZ_cg_condnum, AZ_cgs, AZ_tfqmr, AZ_bicgstab
                aParameterLists( 0 ).set( "AZ_solver", AZ_gmres );

                // Allowable Aztec solver iterations
                aParameterLists( 0 ).set( "AZ_max_iter", 500 );

                // Allowable Aztec iterative residual
                aParameterLists( 0 ).set( "rel_residual", 1e-16 );

                // set Az_conv -convergence criteria
                // options are AZ_r0, AZ_rhs, AZ_Anorm, AZ_noscaled, AZ_sol
                aParameterLists( 0 ).set( "AZ_conv", AZ_r0 );

                // set Az_diagnostic parameters
                // Set whether or not diagnostics for every linear iteration are printed or not. options are AZ_all, AZ_none
                aParameterLists( 0 ).set( "AZ_diagnostics", AZ_all );

                // set AZ_output options
                // options are AZ_all, AZ_none, AZ_warnings, AZ_last, AZ_summary
                aParameterLists( 0 ).set( "AZ_output", AZ_all );

                // Determines the submatrices factored with the domain decomposition algorithms
                // Option to specify with how many rows from other processors each processor\u2019s local submatrix is augmented.
                aParameterLists( 0 ).set( "AZ_overlap", 1 );

                // Determines how overlapping subdomain results are combined when different processors have computed different values for the same unknown.
                // Options are AZ_standard, AZ_symmetric
                aParameterLists( 0 ).set( "AZ_type_overlap", AZ_standard );

                // Determines whether RCM reordering will be done in conjunction with domain decomposition incomplete factorizations.
                // Option to enable (=1) or disable (=0) the Reverse Cuthill\u2013McKee (RCM) algorithm to reorder system equations for smaller bandwidth
                aParameterLists( 0 ).set( "AZ_reorder", 1 );

                // Use preconditioner from a previous Iterate() call
                // Option are AZ_calc, AZ_recalc, AZ_reuse
                aParameterLists( 0 ).set( "AZ_pre_calc", AZ_calc );

                // Determines  whether  matrix  factorization  information will be kept after this solve
                // for example for preconditioner_recalculation
                aParameterLists( 0 ).set( "AZ_keep_info", 0 );

                //--------------------------GMRES specific solver parameters--------------------------------------------------------------------------
                // Set AZ_kspace
                // Krylov subspace size for restarted GMRES
                // Setting mKrylovSpace larger improves the robustness, decreases iteration count, but increases memory consumption.
                // For very difficult problems, set it equal to the maximum number of iterations.
                aParameterLists( 0 ).set( "AZ_kspace", 250 );

                // Set AZ_orthog
                // AZ_classic or AZ_modified
                aParameterLists( 0 ).set( "AZ_orthog", AZ_classic );

                // Set AZ_rthresh
                // Parameter used to modify the relative magnitude of the diagonal entries of the matrix that is used to compute
                // any of the incomplete factorization preconditioners
                aParameterLists( 0 ).set( "AZ_rthresh", 0.0 );

                // Set AZ_athresh
                // Parameter used to modify the absolute magnitude of the diagonal entries of the matrix that is used to compute
                // any of the incomplete factorization preconditioners
                aParameterLists( 0 ).set( "AZ_athresh", 0.0 );

                //--------------------------Preconsitioner specific parameters--------------------------------------------------------------------------
                // Determine which preconditioner is used
                // Options are AZ_none, AZ_Jacobi, AZ_sym_GS, AZ_Neumann, AZ_ls, AZ_dom_decomp,
                aParameterLists( 0 ).set( "AZ_precond", AZ_dom_decomp );

                // Set preconditioner subdomain solve - direct solve or incomplete
                // Options are AZ_lu, AZ_ilut, AZ_ilu, AZ_rilu, AZ_bilu, AZ_icc
                aParameterLists( 0 ).set( "AZ_subdomain_solve", AZ_ilut );

                // Set preconditioner polynomial order - polynomial preconditioning, Gauss-Seidel, Jacobi
                aParameterLists( 0 ).set( "AZ_poly_ord", 3 );

                // Set drop tolerance - for LU, ILUT
                aParameterLists( 0 ).set( "AZ_drop", 1.0e-12 );

                // Set level of graph fill in - for ilu(k), icc(k), bilu(k)
                aParameterLists( 0 ).set( "AZ_graph_fill", 3 );

                // Set ilut fill
                aParameterLists( 0 ).set( "AZ_ilut_fill", 4.0 );

                // Set Damping or relaxation parameter used for RILU
                aParameterLists( 0 ).set( "AZ_omega", 1.0 );

                // use the correct precondioiner
                aParameterLists( 0 ).set( "preconditioners", "0" );

                // Preconditioner using ifpack
                // aParameterLists( 0 ).set( "ifpack_prec_type"    , "ILU");
                // aParameterLists( 0 ).set( "fact: level-of-fill" ,  3       );
                // aParameterLists( 0 ).set( "fact: drop tolerance",  1.0e-2 );
                // aParameterLists( 0 ).set( "prec_reuse"          ,  false );
                break;
        }

        // eigenvalue problem to find the largest eigenvalues
        aParameterLists( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::EIGEN_SOLVER ) );
        aParameterLists( 0 ).set( "Eigen_Algorithm", "EIGALG_BLOCK_KRYLOV_SCHUR" );
        aParameterLists( 0 ).set( "Verbosity", false );
        aParameterLists( 0 ).set( "Which", "SM" );
        aParameterLists( 0 ).set( "Block_Size", 5 );          // Block Size should be same as Number of Eigen values
        aParameterLists( 0 ).set( "NumFreeDofs", 1000 );      // For 2D problem of rectangular elements number of free dofs = 2*node_x*node_y
        aParameterLists( 0 ).set( "Num_Eig_Vals", 5 );        // Number of Eigen values should be same as Block Size
        aParameterLists( 0 ).set( "Num_Blocks", 3 );          // Number of Blocks should satisfy : Num_Blocks*Block_Size < InitVec Length
        aParameterLists( 0 ).set( "MaxSubSpaceDims", 15 );    // Max Subspace Dimension = 3*Block_Size*Num_Eig_Vals
        aParameterLists( 0 ).set( "Initial_Guess", 0 );
        aParameterLists( 0 ).set( "MaxRestarts", 100 );
        aParameterLists( 0 ).set( "Convergence_Tolerance", 1e-02 );
        aParameterLists( 0 ).set( "Relative_Convergence_Tolerance", true );
        aParameterLists( 0 ).set( "preconditioners", "1" );
        aParameterLists( 0 ).set( "preconditioners_linear_operator", "0" );

        // eigenvalue problem to find the smallest eigenvalues
        aParameterLists( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::EIGEN_SOLVER ) );
        aParameterLists( 0 ).set( "Eigen_Algorithm", "EIGALG_BLOCK_KRYLOV_SCHUR_AMESOS" );
        aParameterLists( 0 ).set( "Verbosity", false );
        aParameterLists( 0 ).set( "Which", "LM" );
        aParameterLists( 0 ).set( "Block_Size", 5 );          // Block Size should be same as Number of Eigen values
        aParameterLists( 0 ).set( "NumFreeDofs", 1000 );      // For 2D problem of rectangular elements number of free dofs = 2*node_x*node_y
        aParameterLists( 0 ).set( "Num_Eig_Vals", 5 );        // Number of Eigen values should be same as Block Size
        aParameterLists( 0 ).set( "Num_Blocks", 3 );          // Number of Blocks should satisfy : Num_Blocks*Block_Size < InitVec Length
        aParameterLists( 0 ).set( "MaxSubSpaceDims", 15 );    // Max Subspace Dimension = 3*Block_Size*Num_Eig_Vals
        aParameterLists( 0 ).set( "Initial_Guess", 0 );
        aParameterLists( 0 ).set( "MaxRestarts", 100 );
        aParameterLists( 0 ).set( "Convergence_Tolerance", 1e-02 );
        aParameterLists( 0 ).set( "Relative_Convergence_Tolerance", true );
        aParameterLists( 0 ).set( "preconditioners", "1" );

        // ML precondioner for the linear solver
        aParameterLists( 7 ).add_parameter_list( moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::ML ) );
        aParameterLists( 7 ).set( "ml_prec_type", "NSSA" );
        aParameterLists( 7 ).set( "PDE equations", 1 );

        // Ifpack precondioner for the eigen solve
        aParameterLists( 7 ).add_parameter_list( moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::IFPACK ) );
        aParameterLists( 7 ).set( "ifpack_prec_type", "ILU" );
        aParameterLists( 7 ).set( "fact: level-of-fill", 2 );
        aParameterLists( 7 ).set( "fact: drop tolerance", 1.0e-2 );
        aParameterLists( 7 ).set( "prec_reuse", false );

        //------------------------------------------------------------------------------

        aParameterLists( 1 ).add_parameter_list( moris::prm::create_linear_solver_parameter_list() );    // linear problem index 0
        aParameterLists( 1 ).set( "DLA_Linear_solver_algorithms", "0" );
        aParameterLists( 1 ).set( "DLA_operator_condition_number_with_moris", "dense" );
        aParameterLists( 1 ).set( "DLA_prec_operator_condition_number_with_moris", "dense" );

        aParameterLists( 1 ).add_parameter_list( moris::prm::create_linear_solver_parameter_list() );    // eigen problem index 1
        aParameterLists( 1 ).set( "DLA_Linear_solver_algorithms", "1" );
        aParameterLists( 1 ).set( "RHS_Matrix_Type", "IdentityMat" );

        aParameterLists( 1 ).add_parameter_list( moris::prm::create_linear_solver_parameter_list() );    // eigen problem index 1
        aParameterLists( 1 ).set( "DLA_Linear_solver_algorithms", "2" );
        aParameterLists( 1 ).set( "RHS_Matrix_Type", "IdentityMat" );

        //------------------------------------------------------------------------------

        aParameterLists( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );    // nonlinear algorithm index 0
        aParameterLists( 2 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists( 2 ).set( "NLA_Linear_solver", 0 );
        aParameterLists( 2 ).set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        aParameterLists( 2 ).set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
        aParameterLists( 2 ).set( "NLA_max_iter", tNLA_max_iter );

        aParameterLists( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
        aParameterLists( 2 ).set( "NLA_Linear_solver", 1 );
        aParameterLists( 2 ).set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        aParameterLists( 2 ).set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
        aParameterLists( 2 ).set( "NLA_max_iter", tNLA_max_iter );
        aParameterLists( 2 ).set( "NLA_is_eigen_problem", true );

        aParameterLists( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
        aParameterLists( 2 ).set( "NLA_Linear_solver", 2 );
        aParameterLists( 2 ).set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        aParameterLists( 2 ).set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
        aParameterLists( 2 ).set( "NLA_max_iter", tNLA_max_iter );
        aParameterLists( 2 ).set( "NLA_is_eigen_problem", true );

        //------------------------------------------------------------------------------

        aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );    // 1: thermal subproblem
        aParameterLists( 3 ).set( "NLA_Nonlinear_solver_algorithms", "0,1,2" );         // set nonlinear algorithm with index 0
        aParameterLists( 3 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists( 3 ).set( "NLA_DofTypes", "TEMP" );

        // ----------------------------------------------------------

        aParameterLists( 4 ).add_parameter_list( moris::prm::create_time_solver_algorithm_parameter_list() );
        aParameterLists( 4 ).set( "TSA_Nonlinear_Solver", 0 );                // using NLBGS for forward problem
        aParameterLists( 4 ).set( "TSA_Nonlinear_Sensitivity_Solver", 0 );    // using monlithic for sensitivity problem

        if ( tUseTimeContinuity )
        {
            aParameterLists( 4 ).set( "TSA_Num_Time_Steps", tTimeSteps );
            aParameterLists( 4 ).set( "TSA_Time_Frame", tMaxTime );
        }

        //------------------------------------------------------------------------------

        aParameterLists( 5 ).add_parameter_list( moris::prm::create_time_solver_parameter_list() );
        aParameterLists( 5 ).set( "TSA_DofTypes", "TEMP" );
        aParameterLists( 5 ).set( "TSA_Initialize_Sol_Vec", "TEMP,0.0" );
        aParameterLists( 5 ).set( "TSA_Output_Indices", "0" );
        aParameterLists( 5 ).set( "TSA_Output_Criteria", "Output_Criterion" );

        if ( tUseTimeContinuity )
        {
            aParameterLists( 5 ).set( "TSA_time_level_per_type", "TEMP,2" );
        }
        else
        {
            aParameterLists( 5 ).set( "TSA_time_level_per_type", "TEMP,1" );
        }

        //------------------------------------------------------------------------------

        aParameterLists( 6 ).add_parameter_list( moris::prm::create_solver_warehouse_parameterlist() );
        aParameterLists( 6 ).set( "SOL_save_operator_to_matlab", "heat" );
    }

    //------------------------------------------------------------------------------
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
