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
    uint tNumElemX = std::ceil( tDimX / tApproxEleSize );
    uint tNumElemY = std::ceil( tDimY / tApproxEleSize );

    int tLevelsetOrder = tInterpolationOrder;
    int tDispOrder     = tInterpolationOrder;

    int tLevelsetInitialRef = 0;
    int tDispInitialRef     = 0;
    int tRefineBuffer       = 1;

    std::string tLagrangeOrder   = std::to_string( std::max( tLevelsetOrder, tDispOrder ) );
    std::string tBsplineOrder    = std::to_string( tLevelsetOrder ) + "," + std::to_string( tDispOrder );
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

        aParameterLists.set( "is_optimization_problem", tIsOptimization );
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

        aParameterLists.set( "number_of_elements_per_dimension", tNumElemX, tNumElemY );
        aParameterLists.set( "domain_dimensions", tDimX, tDimY );
        aParameterLists.set( "domain_offset", -0.0003, 0.0 );

        aParameterLists.set( "lagrange_output_meshes", "0" );

        aParameterLists.set( "lagrange_orders", tLagrangeOrder );
        aParameterLists.set( "lagrange_pattern", tLagrangePattern );

        aParameterLists.set( "bspline_orders", tBsplineOrder );
        aParameterLists.set( "bspline_pattern", "0,1" );

        aParameterLists.set( "lagrange_to_bspline", "0,1" );

        aParameterLists.set( "refinement_buffer", tRefineBuffer );
        aParameterLists.set( "staircase_buffer", tRefineBuffer );

        aParameterLists.set( "initial_refinement", tLevelsetInitialRef, tDispInitialRef );

        aParameterLists.set( "basis_function_vtk_file", "basisinhmr.vtk" );
        aParameterLists.set( "lagrange_mesh_output_file_name", "lagrangehmr.exo" );
    }

    /* ------------------------------------------------------------------------ */

    void
    XTKParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "decompose", true );
        aParameterLists.set( "decomposition_type", "conformal" );
        aParameterLists.set( "use_SPG_based_enrichment", true );
        aParameterLists.set( "enrich_mesh_indices", "0" );
        aParameterLists.set( "ghost_stab", tUseGhost );
        aParameterLists.set( "multigrid", false );
        aParameterLists.set( "verbose", false );
        aParameterLists.set( "high_to_low_dbl_side_sets", true );
        aParameterLists.set( "print_enriched_ig_mesh", true );
        aParameterLists.set( "exodus_output_XTK_ig_mesh", true );
        aParameterLists.set( "exodus_output_XTK_ip_mesh", true );
        aParameterLists.set( "activate_basis_agglomeration", true );
        aParameterLists.set( "write_enrichment_fields", true );
        aParameterLists.set( "write_enrichment_fields_probe_spheres", "1.0,0.5,0.5,1.0" );
    }

    /* ------------------------------------------------------------------------ */

    void
    GENParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists.set( "IQI_types",
                "IQIInputThermalEnergy",
                "IQIDiffusiveLower",
                "IQIDiffusiveUpper",
                "IQIDiffusiveFront",
                "IQIDiffusiveBack",
                "IQIStoredThermalEnergy",
                "IQIInputThermalEnergySurface" );

        aParameterLists.set( "number_of_phases", 2 );
        aParameterLists.set( "phase_function_name", "get_phase_index" );
        aParameterLists.set( "output_mesh_file", tGENOutputFile );
        aParameterLists.set( "time_offset", 10.0 );

        // Inclusions
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists.set( "field_function_name", "Func_Sphere" );
        aParameterLists.set( "number_of_refinements", tInterfaceRefinementSphere );
        aParameterLists.set( "refinement_mesh_index", 0 );
        aParameterLists.set( "isocontour_threshold", 0.0 );
        aParameterLists.set( "isocontour_tolerance", 1.0e-12 );
        aParameterLists.set( "intersection_tolerance", 1.0e-12 );
        aParameterLists.set( "use_multilinear_interpolation", true );
        aParameterLists.set( "discretization_mesh_index", 0 );
        aParameterLists.set( "discretization_lower_bound", -2.0 );
        aParameterLists.set( "discretization_upper_bound", 2.0 );

    }

    /* ------------------------------------------------------------------------ */

    void
    FEMParameterList( Module_Parameter_Lists& aParameterLists )
    {

         //------------------------------------------------------------------------------

        aParameterLists( FEM::PHASES ).add_parameter_list();
        aParameterLists.set( "phase_name", "PhaseSolid" );
        aParameterLists.set( "phase_indices", "0" );

        aParameterLists( FEM::PHASES ).add_parameter_list();
        aParameterLists.set( "phase_name", "PhaseVoid" );
        aParameterLists.set( "phase_indices", "1" );

        //------------------------------------------------------------------------------
        // fill the property part of the parameter list

        // soild properties ------------------------------------------------------------

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDensity" );
        aParameterLists.set( "function_parameters", tDensity );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropCapacity" );
        aParameterLists.set( "function_parameters", tCapacity );

        // create  conductivity property
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropConductivity" );
        aParameterLists.set( "function_parameters", tConductivity );

        // BC properties ---------------------------------------------------------------

        // create inlet temperature property
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropInletTemp" );
        aParameterLists.set( "function_parameters", moris_to_string( tInTemp ) );

        // create inlet temperature property
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropInletTemp2" );
        aParameterLists.set( "function_parameters", moris_to_string( 1.0 ) );
        aParameterLists.set( "value_function", "Func_TempDistro" );

        // create inlet temperature property
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropAnalSol" );
        aParameterLists.set( "function_parameters", moris_to_string( 1.0 ) );
        aParameterLists.set( "value_function", "Func_Anal_Sol" );

        // create heat load property
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropVolumetricHeatFlux" );
        aParameterLists.set( "value_function", "Func_HeatLoad" );

        // create heat load property
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropSurfaceHeatFlux" );
        aParameterLists.set( "function_parameters", tSurfaceHeatLoad );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropInitialTemp" );
        aParameterLists.set( "function_parameters", moris_to_string( tInTemp ) );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropTimeContinuity" );
        aParameterLists.set( "function_parameters", moris_to_string( tTimePenalty * std::stod( tDensity ) * std::stod( tCapacity ) ) );

        // Stored Thermal Energy
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropStoredThermalEnergy" );
        aParameterLists.set( "value_function", "Func_StoredThermalEnergy" );

        //------------------------------------------------------------------------------
        // fill the constitutive model part of the parameter list

        //  CM --------------------------------------------------------------------
        // create  diffusion CM
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMDiffusion" );
        aParameterLists.set( "phase_name", "PhaseSolid" );
        aParameterLists.set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists.set( "properties",
                "PropConductivity,Conductivity;"
                "PropDensity     ,Density;"
                "PropCapacity    ,HeatCapacity" );

        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMDiffusionL2" );
        aParameterLists.set( "phase_name", "PhaseSolid" );
        aParameterLists.set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "L2", "Temperature" ) );
        aParameterLists.set( "properties",
                "PropConductivity,Conductivity;"
                "PropDensity     ,Density;"
                "PropCapacity    ,HeatCapacity" );

        //------------------------------------------------------------------------------
        // fill the stabilization parameter part of the parameter list

        // create Nitsche for  temperature
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPNitsche" );
        aParameterLists.set( "leader_phase_name", "PhaseSolid" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        aParameterLists.set( "function_parameters", tNitschePenThermal );
        aParameterLists.set( "leader_properties", "PropConductivity,Material" );

        // create ghost penalty  temperature
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPGPTemp" );
        aParameterLists.set( "leader_phase_name", "PhaseSolid" );
        aParameterLists.set( "follower_phase_name", "PhaseSolid" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        aParameterLists.set( "function_parameters", "0.05" );
        aParameterLists.set( "leader_properties", "PropConductivity,Material" );

        //------------------------------------------------------------------------------
        // fill the IWG part of the parameter list

        //  bulk IWGs -------------------------------------------------------------
        // diffusion
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGDiffusionBulk" );
        aParameterLists.set( "leader_phase_name", "PhaseSolid" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_properties", "PropVolumetricHeatFlux,Load" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusion,Diffusion" );

        //--------------------------------------------------------------------------------
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGDiffusionBulk" );
        aParameterLists.set( "leader_phase_name", "PhaseSolid" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        aParameterLists.set( "dof_residual", "L2" );
        aParameterLists.set( "leader_properties", "PropVolumetricHeatFlux,Load" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusionL2,Diffusion" );

        // Inlet BC IWG ----------------------------------------------------------------

        //
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGInletTempTopBottom" );
        aParameterLists.set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "leader_phase_name", "PhaseSolid" );
        aParameterLists.set( "side_ordinals", "1,3" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_properties", "PropInletTemp,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusion,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPNitsche,DirichletNitsche" );

        //--------------------------------------------------------------------------------

        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGInletTempRight" );
        aParameterLists.set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "leader_phase_name", "PhaseSolid" );
        aParameterLists.set( "neighbor_phases", "PhaseVoid" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_properties", "PropInletTemp2,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusion,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPNitsche,DirichletNitsche" );

        // Inlet BC IWG ----------------------------------------------------------------

        //
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGInletTempTopBottom" );
        aParameterLists.set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "leader_phase_name", "PhaseSolid" );
        aParameterLists.set( "side_ordinals", "1,3" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "L2" );
        aParameterLists.set( "leader_properties", "PropInletTemp,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusionL2,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPNitsche,DirichletNitsche" );

        //--------------------------------------------------------------------------------

        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGInletTempRight" );
        aParameterLists.set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "leader_phase_name", "PhaseSolid" );
        aParameterLists.set( "neighbor_phases", "PhaseVoid" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "L2" );
        aParameterLists.set( "leader_properties", "PropInletTemp2,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusionL2,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPNitsche,DirichletNitsche" );
        // Outlet BC IWG ----------------------------------------------------------------

        // Time continuity ----------------------------------------------------------------

        if ( tUseTimeContinuity )
        {
            // Time continuity temperature
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "IWGTimeContinuity" );
            aParameterLists.set( "IWG_type", ( fem::IWG_Type::TIME_CONTINUITY_DOF ) );
            aParameterLists.set( "dof_residual", "TEMP" );
            aParameterLists.set( "leader_properties",
                    "PropTimeContinuity,WeightCurrent;"
                    "PropTimeContinuity,WeightPrevious;"
                    "PropInitialTemp,InitialCondition" );
            aParameterLists.set( "leader_phase_name", "PhaseSolid" );
            aParameterLists.set( "time_continuity", true );
            }

        // Ghost  ----------------------------------------------------------------

        if ( tUseGhost )
        {
            // ghost  temperature
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "IWGGPTemp" );
            aParameterLists.set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            aParameterLists.set( "leader_phase_name", "PhaseSolid" );
            aParameterLists.set( "follower_phase_name", "PhaseSolid" );
            aParameterLists.set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists.set( "dof_residual", "TEMP" );
            aParameterLists.set( "stabilization_parameters", "SPGPTemp,GhostSP" );
            aParameterLists.set( "ghost_order", (uint)tDispOrder );

            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "IWGGPTemp" );
            aParameterLists.set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            aParameterLists.set( "leader_phase_name", "PhaseSolid" );
            aParameterLists.set( "follower_phase_name", "PhaseSolid" );
            aParameterLists.set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists.set( "dof_residual", "L2" );
            aParameterLists.set( "stabilization_parameters", "SPGPTemp,GhostSP" );
            aParameterLists.set( "ghost_order", (uint)tDispOrder );
        }

        //------------------------------------------------------------------------------
        // fill the IQI part of the parameter list

        // input thermal energy
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIAnalTemp" );
        aParameterLists.set( "leader_phase_name", "PhaseSolid" );
        aParameterLists.set( "IQI_type", ( fem::IQI_Type::PROPERTY ) );
        aParameterLists.set( "leader_properties", "PropAnalSol,Property" );
        aParameterLists.set( "vectorial_field_index", 0 );

        // temperature
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkTEMP" );
        aParameterLists.set( "leader_phase_name", "PhaseSolid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "TEMP" );
        aParameterLists.set( "vectorial_field_index", 0 );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkL2" );
        aParameterLists.set( "leader_phase_name", "PhaseSolid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "L2" );
        aParameterLists.set( "vectorial_field_index", 0 );

        // input thermal energy
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIInputThermalEnergy" );
        aParameterLists.set( "leader_phase_name", "PhaseSolid" );
        aParameterLists.set( "IQI_type", ( fem::IQI_Type::PROPERTY ) );
        aParameterLists.set( "leader_properties", "PropVolumetricHeatFlux,Property" );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkL2Error" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::L2_ERROR_ANALYTIC );
        aParameterLists.set( "dof_quantity", "TEMP" );
        aParameterLists.set( "leader_properties", "PropAnalSol,L2Check" );
        aParameterLists.set( "leader_phase_name", "PhaseSolid" );

        // input thermal energy
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIInputThermalEnergySurface" );
        aParameterLists.set( "IQI_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "leader_phase_name", "PhaseSolid" );
        aParameterLists.set( "side_ordinals", "2" );
        aParameterLists.set( "IQI_type", ( fem::IQI_Type::PROPERTY ) );
        aParameterLists.set( "leader_properties", "PropSurfaceHeatFlux,Property" );

        //  diffusive flux
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIDiffusiveLower" );
        aParameterLists.set( "IQI_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "leader_phase_name", "PhaseSolid" );
        aParameterLists.set( "side_ordinals", "1" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::THERMAL_ENERGY_DIFFUSIVE_FLUX );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusion,Diffusion" );

        //  diffusive flux
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIDiffusiveUpper" );
        aParameterLists.set( "IQI_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "leader_phase_name", "PhaseSolid" );
        aParameterLists.set( "side_ordinals", "3" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::THERMAL_ENERGY_DIFFUSIVE_FLUX );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusion,Diffusion" );

        //  diffusive flux
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIDiffusiveFront" );
        aParameterLists.set( "IQI_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "leader_phase_name", "PhaseSolid" );
        aParameterLists.set( "side_ordinals", "4" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::THERMAL_ENERGY_DIFFUSIVE_FLUX );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusion,Diffusion" );

        //  diffusive flux
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIDiffusiveBack" );
        aParameterLists.set( "IQI_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "leader_phase_name", "PhaseSolid" );
        aParameterLists.set( "side_ordinals", "2" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::THERMAL_ENERGY_DIFFUSIVE_FLUX );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusion,Diffusion" );

        // stored thermal energy
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIStoredThermalEnergy" );
        aParameterLists.set( "leader_phase_name", "PhaseSolid" );
        aParameterLists.set( "IQI_type", ( fem::IQI_Type::PROPERTY ) );
        aParameterLists.set( "leader_properties", "PropStoredThermalEnergy,Property" );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkEigen" );
        aParameterLists.set( "IQI_type", ( fem::IQI_Type::EIGEN_VECTOR ) );
        aParameterLists.set( "function_parameters", "0" );
        aParameterLists.set( "dof_quantity", "TEMP" );
        aParameterLists.set( "leader_phase_name", "PhaseSolid" );
        aParameterLists.set( "vectorial_field_index", 0 );

        if ( gTestIndex == 1 )
        {
            aParameterLists( FEM::IQI ).add_parameter_list();
            aParameterLists.set( "IQI_name", "IQIBulkEigenVal" );
            aParameterLists.set( "IQI_type", ( fem::IQI_Type::EIGEN_VALUE ) );
            aParameterLists.set( "function_parameters", "0" );
            aParameterLists.set( "dof_quantity", "TEMP" );
            aParameterLists.set( "leader_phase_name", "PhaseSolid" );
            aParameterLists.set( "vectorial_field_index", 0 );
            }

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
        std::cout << "gTestIndex_" + std::to_string( par_rank() ) + ": " << gTestIndex << '\n';
        gTestIndex == 0 ? create_trilinos_solver_parameter_list( aParameterLists ) : create_petsc_solver_parameter_list( aParameterLists );
    }

    void
    MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "TEMP", 0 );
        aParameterLists.set( "number_eigen_vectors", gTestIndex == 0 ? 5 : 1 );
        aParameterLists.set( "order_adofs_by_host", true );
    }

    void
    VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "File_Name", std::pair< std::string, std::string >( "./", tExoFile ) );
        aParameterLists.set( "Mesh_Type", vis::VIS_Mesh_Type::STANDARD );
        aParameterLists.set( "Set_Names", tAllPhaseInterfaces );
        aParameterLists.set( "Field_Names", "TEMP,HEAT,TEMP_A,L2Nodal,L2Glob,DiffLower,EigenVec,EigenVal,TEMPL" );
        aParameterLists.set( "Field_Type", "NODAL,NODAL,NODAL,NODAL,GLOBAL,GLOBAL,NODAL,NODAL,NODAL" );
        aParameterLists.set( "IQI_Names", "IQIBulkTEMP,IQIInputThermalEnergy,IQIAnalTemp,IQIBulkL2Error,IQIBulkL2Error,IQIDiffusiveLower,IQIBulkEigen,IQIBulkEigenVal,IQIBulkL2" );
        aParameterLists.set( "Save_Frequency", 1 );
        aParameterLists.set( "Time_Offset", 10.0 );
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
        aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::PETSC );
        aParameterLists.set( "KSPType", "fgmres" );
        aParameterLists.set( "preconditioners", "0" );    // 10 shift_invert

        aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::PETSC );
        aParameterLists.set( "KSPType", "preonly" );
        aParameterLists.set( "preconditioners", "0" );    // 10 shift_invert

        aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::PETSC );
        aParameterLists.set( "KSPType", "preonly" );
        aParameterLists.set( "preconditioners", "1" );    // 10 shift_invert

        // find max eigen value
        aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::SLEPC_SOLVER );
        aParameterLists.set( "Eigen_Algorithm", "power" );
        aParameterLists.set( "Which", std::string( "LM" ) );
        aParameterLists.set( "Num_Eig_Vals", 1 );
        aParameterLists.set( "sub_linear_solver", "1" );    // 10 shift_invert
        aParameterLists.set( "is_symmetric", false );       // 10 shift_invert
        aParameterLists.set( "Update_Flag", true );         // 10 shift_invert

        // find min eigen value
        aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::SLEPC_SOLVER );
        aParameterLists.set( "Eigen_Algorithm", "power" );
        aParameterLists.set( "Which", std::string( "LM" ) );
        aParameterLists.set( "Num_Eig_Vals", 1 );
        aParameterLists.set( "STType", "shift_invert" );
        aParameterLists.set( "sub_linear_solver", "2" );    // 10 shift_invert
        aParameterLists.set( "is_symmetric", false );       // 10 shift_invert
        aParameterLists.set( "Update_Flag", false );        // 10 shift_invert

        // precondioerr
        aParameterLists( SOL::PRECONDITIONERS ).add_parameter_list(  sol::PreconditionerType::PETSC );
        aParameterLists.set( "PCType", "none" );

        // Ifpack precondioner for the eigen solve
        aParameterLists( SOL::PRECONDITIONERS ).add_parameter_list(  sol::PreconditionerType::PETSC );
        aParameterLists.set( "PCType", "mumps" );

        // **
        aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list();    // linear problem index 0
        aParameterLists.set( "DLA_Linear_solver_algorithms", "0" );

        aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list();    // eigen problem index 1
        aParameterLists.set( "DLA_Linear_solver_algorithms", "3" );
        aParameterLists.set( "RHS_Matrix_Type", "IdentityMat" );

        aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list();    // eigen problem index 1
        aParameterLists.set( "DLA_Linear_solver_algorithms", "4" );
        aParameterLists.set( "RHS_Matrix_Type", "IdentityMat" );

        //------------------------------------------------------------------------------

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();    // nonlinear algorithm index 0
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists.set( "NLA_Linear_solver", 0 );
        aParameterLists.set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        aParameterLists.set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
        aParameterLists.set( "NLA_max_iter", tNLA_max_iter );

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "NLA_Linear_solver", 1 );
        aParameterLists.set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        aParameterLists.set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
        aParameterLists.set( "NLA_max_iter", tNLA_max_iter );
        aParameterLists.set( "NLA_is_eigen_problem", true );

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "NLA_Linear_solver", 2 );
        aParameterLists.set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        aParameterLists.set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
        aParameterLists.set( "NLA_max_iter", tNLA_max_iter );
        aParameterLists.set( "NLA_is_eigen_problem", true );

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();    // nonlinear algorithm index 1
        aParameterLists.set( "NLA_Solver_Implementation", static_cast< uint >( moris::NLA::NonlinearSolverType::NLBGS_SOLVER ) );
        aParameterLists.set( "NLA_rel_res_norm_drop", 1.0e-9 );
        aParameterLists.set( "NLA_max_iter", 1 );

        //------------------------------------------------------------------------------

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();    // 1: thermal subproblem
        aParameterLists.set( "NLA_Nonlinear_solver_algorithms", "0,1" );           // set nonlinear algorithm with index 0
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists.set( "NLA_DofTypes", "TEMP" );

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();    // 1: thermal subproblem
        aParameterLists.set( "NLA_Nonlinear_solver_algorithms", "0" );             // set nonlinear algorithm with index 0
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists.set( "NLA_DofTypes", "L2" );

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();    // nonlinear solver index 2
        aParameterLists.set( "NLA_Solver_Implementation", static_cast< uint >( moris::NLA::NonlinearSolverType::NLBGS_SOLVER ) );
        aParameterLists.set( "NLA_Nonlinear_solver_algorithms", "3" );    // set nonlinear algorithm with index 1.
        aParameterLists.set( "NLA_Sub_Nonlinear_Solver", "0,1" );         // set sub nonlinear solvers with index 0 and 1
        aParameterLists.set( "NLA_DofTypes", "TEMP;L2" );

        // ----------------------------------------------------------

        aParameterLists( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "TSA_Nonlinear_Solver", 2 );                // using NLBGS for forward problem
        aParameterLists.set( "TSA_Nonlinear_Sensitivity_Solver", 2 );    // using monlithic for sensitivity problem

        if ( tUseTimeContinuity )
        {
            aParameterLists.set( "TSA_Num_Time_Steps", tTimeSteps );
            aParameterLists.set( "TSA_Time_Frame", tMaxTime );
        }

        //------------------------------------------------------------------------------

        aParameterLists( SOL::TIME_SOLVERS ).add_parameter_list();
        aParameterLists.set( "TSA_DofTypes", "TEMP;L2" );
        // aParameterLists.set( "TSA_Initialize_Sol_Vec", "TEMP,0.0" );
        aParameterLists.set( "TSA_Output_Indices", "0" );
        aParameterLists.set( "TSA_Output_Criteria", "Output_Criterion" );

        if ( tUseTimeContinuity )
        {
            aParameterLists.set( "TSA_time_level_per_type", "TEMP,2" );
        }
        else
        {
            aParameterLists.set( "TSA_time_level_per_type", "TEMP,1" );
        }

        //------------------------------------------------------------------------------

        aParameterLists( SOL::SOLVER_WAREHOUSE ).set( "SOL_save_operator_to_matlab", "Jacobian_petsc" );
        aParameterLists.set( "SOL_TPL_Type", static_cast< uint >( sol::MapType::Petsc ) );
    }

    //------------------------------------------------------------------------------
    void
    create_trilinos_solver_parameter_list( Module_Parameter_Lists& aParameterLists )
    {


        sint tVerbosity = Belos::Errors + Belos::Warnings + Belos::IterationDetails + Belos::TimingDetails + Belos::StatusTestDetails + Belos::FinalSummary;

        switch ( tSolverType )
        {
            case 0:    // Amesos - Pardiso
                aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::AMESOS_IMPL );

                //                 aParameterLists.set( "Solver_Type", "Amesos_Mumps");
                break;

            case 1:    // Amesos - Mumps
                aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::AMESOS_IMPL );

                aParameterLists.set( "Solver_Type", "Amesos_Mumps" );
                break;

            case 2:    // Belos
                aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::BELOS_IMPL );

                // Solver type: GMRES, Flexible GMRES, Block CG , PseudoBlockCG, Stochastic CG, Recycling GMRES, Recycling CG, MINRES, LSQR, TFQMR
                //              Pseudoblock TFQMR, Seed GMRES, Seed CG
                aParameterLists.set( "Solver Type", "Flexible GMRES" );

                // Diagnostics: Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::StatusTestDetails
                aParameterLists.set( "Verbosity", tVerbosity );

                // Maximum number of blocks in Krylov factorization
                aParameterLists.set( "Num Blocks", 250 );

                // Block size to be used by iterative solver
                aParameterLists.set( "Block Size", 1 );

                // Allowable Belos solver iterations
                aParameterLists.set( "Maximum Iterations", 500 );

                // Allowable Belos solver iterations
                aParameterLists.set( "Maximum Restarts", 2 );

                // Convergence criteria
                aParameterLists.set( "Convergence Tolerance", 1e-9 );

                // Left or right preconditioner
                aParameterLists.set( "Left-right Preconditioner", "right" );

                // Preconditioner

                // ifpack - ILU
                // aParameterLists.set( "ifpack_prec_type",  "ILU");
                // aParameterLists.set( "fact: level-of-fill",  5 );

                // ifpack - ILUT
                aParameterLists.set( "ifpack_prec_type", "ILUT" );
                aParameterLists.set( "fact: ilut level-of-fill", 5.0 );
                aParameterLists.set( "fact: drop tolerance", 1e-12 );

                // ifpack with direct solve
                // aParameterLists.set( "ifpack_prec_type",  "amesos");
                // aParameterLists.set( "amesos: solver type", "Amesos_Pardiso");

                // AMG with defaults for non-symmetric system
                // aParameterLists.set( "ml_prec_type",  "NSSA");
                // aParameterLists.set( "PDE equations", 3);

                break;

            case 3:    // AZTEC
                aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::AZTEC_IMPL );

                // options are: AZ_gmres, AZ_gmres_condnum, AZ_cg, AZ_cg_condnum, AZ_cgs, AZ_tfqmr, AZ_bicgstab
                aParameterLists.set( "AZ_solver", AZ_gmres );

                // Allowable Aztec solver iterations
                aParameterLists.set( "AZ_max_iter", 500 );

                // Allowable Aztec iterative residual
                aParameterLists.set( "rel_residual", 1e-16 );

                // set Az_conv -convergence criteria
                // options are AZ_r0, AZ_rhs, AZ_Anorm, AZ_noscaled, AZ_sol
                aParameterLists.set( "AZ_conv", AZ_r0 );

                // set Az_diagnostic parameters
                // Set whether or not diagnostics for every linear iteration are printed or not. options are AZ_all, AZ_none
                aParameterLists.set( "AZ_diagnostics", AZ_all );

                // set AZ_output options
                // options are AZ_all, AZ_none, AZ_warnings, AZ_last, AZ_summary
                aParameterLists.set( "AZ_output", AZ_all );

                // Determines the submatrices factored with the domain decomposition algorithms
                // Option to specify with how many rows from other processors each processor\u2019s local submatrix is augmented.
                aParameterLists.set( "AZ_overlap", 1 );

                // Determines how overlapping subdomain results are combined when different processors have computed different values for the same unknown.
                // Options are AZ_standard, AZ_symmetric
                aParameterLists.set( "AZ_type_overlap", AZ_standard );

                // Determines whether RCM reordering will be done in conjunction with domain decomposition incomplete factorizations.
                // Option to enable (=1) or disable (=0) the Reverse Cuthill\u2013McKee (RCM) algorithm to reorder system equations for smaller bandwidth
                aParameterLists.set( "AZ_reorder", 1 );

                // Use preconditioner from a previous Iterate() call
                // Option are AZ_calc, AZ_recalc, AZ_reuse
                aParameterLists.set( "AZ_pre_calc", AZ_calc );

                // Determines  whether  matrix  factorization  information will be kept after this solve
                // for example for preconditioner_recalculation
                aParameterLists.set( "AZ_keep_info", 0 );

                //--------------------------GMRES specific solver parameters--------------------------------------------------------------------------
                // Set AZ_kspace
                // Krylov subspace size for restarted GMRES
                // Setting mKrylovSpace larger improves the robustness, decreases iteration count, but increases memory consumption.
                // For very difficult problems, set it equal to the maximum number of iterations.
                aParameterLists.set( "AZ_kspace", 250 );

                // Set AZ_orthog
                // AZ_classic or AZ_modified
                aParameterLists.set( "AZ_orthog", AZ_classic );

                // Set AZ_rthresh
                // Parameter used to modify the relative magnitude of the diagonal entries of the matrix that is used to compute
                // any of the incomplete factorization preconditioners
                aParameterLists.set( "AZ_rthresh", 0.0 );

                // Set AZ_athresh
                // Parameter used to modify the absolute magnitude of the diagonal entries of the matrix that is used to compute
                // any of the incomplete factorization preconditioners
                aParameterLists.set( "AZ_athresh", 0.0 );

                //--------------------------Preconsitioner specific parameters--------------------------------------------------------------------------
                // Determine which preconditioner is used
                // Options are AZ_none, AZ_Jacobi, AZ_sym_GS, AZ_Neumann, AZ_ls, AZ_dom_decomp,
                aParameterLists.set( "AZ_precond", AZ_dom_decomp );

                // Set preconditioner subdomain solve - direct solve or incomplete
                // Options are AZ_lu, AZ_ilut, AZ_ilu, AZ_rilu, AZ_bilu, AZ_icc
                aParameterLists.set( "AZ_subdomain_solve", AZ_ilut );

                // Set preconditioner polynomial order - polynomial preconditioning, Gauss-Seidel, Jacobi
                aParameterLists.set( "AZ_poly_ord", 3 );

                // Set drop tolerance - for LU, ILUT
                aParameterLists.set( "AZ_drop", 1.0e-12 );

                // Set level of graph fill in - for ilu(k), icc(k), bilu(k)
                aParameterLists.set( "AZ_graph_fill", 3 );

                // Set ilut fill
                aParameterLists.set( "AZ_ilut_fill", 4.0 );

                // Set Damping or relaxation parameter used for RILU
                aParameterLists.set( "AZ_omega", 1.0 );

                // use the correct precondioiner
                aParameterLists.set( "preconditioners", "0" );

                // Preconditioner using ifpack
                // aParameterLists.set( "ifpack_prec_type"    , "ILU");
                // aParameterLists.set( "fact: level-of-fill" ,  3       );
                // aParameterLists.set( "fact: drop tolerance",  1.0e-2 );
                // aParameterLists.set( "prec_reuse"          ,  false );
                break;
        }

        // eigenvalue problem to find the largest eigenvalues
        aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::EIGEN_SOLVER );
        aParameterLists.set( "Eigen_Algorithm", "EIGALG_BLOCK_KRYLOV_SCHUR" );
        aParameterLists.set( "Verbosity", false );
        aParameterLists.set( "Which", "SM" );
        aParameterLists.set( "Block_Size", 5 );          // Block Size should be same as Number of Eigen values
        aParameterLists.set( "NumFreeDofs", 1000 );      // For 2D problem of rectangular elements number of free dofs = 2*node_x*node_y
        aParameterLists.set( "Num_Eig_Vals", 5 );        // Number of Eigen values should be same as Block Size
        aParameterLists.set( "Num_Blocks", 3 );          // Number of Blocks should satisfy : Num_Blocks*Block_Size < InitVec Length
        aParameterLists.set( "MaxSubSpaceDims", 15 );    // Max Subspace Dimension = 3*Block_Size*Num_Eig_Vals
        aParameterLists.set( "Initial_Guess", 0 );
        aParameterLists.set( "MaxRestarts", 100 );
        aParameterLists.set( "Convergence_Tolerance", 1e-02 );
        aParameterLists.set( "Relative_Convergence_Tolerance", true );
        aParameterLists.set( "preconditioners", "1" );
        aParameterLists.set( "preconditioners_linear_operator", "0" );

        // eigenvalue problem to find the smallest eigenvalues
        aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::EIGEN_SOLVER );
        aParameterLists.set( "Eigen_Algorithm", "EIGALG_BLOCK_KRYLOV_SCHUR_AMESOS" );
        aParameterLists.set( "Verbosity", false );
        aParameterLists.set( "Which", "LM" );
        aParameterLists.set( "Block_Size", 5 );          // Block Size should be same as Number of Eigen values
        aParameterLists.set( "NumFreeDofs", 1000 );      // For 2D problem of rectangular elements number of free dofs = 2*node_x*node_y
        aParameterLists.set( "Num_Eig_Vals", 5 );        // Number of Eigen values should be same as Block Size
        aParameterLists.set( "Num_Blocks", 3 );          // Number of Blocks should satisfy : Num_Blocks*Block_Size < InitVec Length
        aParameterLists.set( "MaxSubSpaceDims", 15 );    // Max Subspace Dimension = 3*Block_Size*Num_Eig_Vals
        aParameterLists.set( "Initial_Guess", 0 );
        aParameterLists.set( "MaxRestarts", 100 );
        aParameterLists.set( "Convergence_Tolerance", 1e-02 );
        aParameterLists.set( "Relative_Convergence_Tolerance", true );
        aParameterLists.set( "preconditioners", "1" );

        // ML precondioner for the linear solver
        aParameterLists( SOL::PRECONDITIONERS ).add_parameter_list(  sol::PreconditionerType::ML );
        aParameterLists.set( "ml_prec_type", "NSSA" );
        aParameterLists.set( "PDE equations", 1 );

        // Ifpack precondioner for the eigen solve
        aParameterLists( SOL::PRECONDITIONERS ).add_parameter_list(  sol::PreconditionerType::IFPACK );
        aParameterLists.set( "ifpack_prec_type", "ILU" );
        aParameterLists.set( "fact: level-of-fill", 2 );
        aParameterLists.set( "fact: drop tolerance", 1.0e-2 );
        aParameterLists.set( "prec_reuse", false );

        //------------------------------------------------------------------------------

        aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list();    // linear problem index 0
        aParameterLists.set( "DLA_Linear_solver_algorithms", "0" );
        aParameterLists.set( "DLA_operator_condition_number_with_moris", "dense" );
        aParameterLists.set( "DLA_prec_operator_condition_number_with_moris", "dense" );

        aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list();    // eigen problem index 1
        aParameterLists.set( "DLA_Linear_solver_algorithms", "1" );
        aParameterLists.set( "RHS_Matrix_Type", "IdentityMat" );

        aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list();    // eigen problem index 1
        aParameterLists.set( "DLA_Linear_solver_algorithms", "2" );
        aParameterLists.set( "RHS_Matrix_Type", "IdentityMat" );

        //------------------------------------------------------------------------------

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();    // nonlinear algorithm index 0
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists.set( "NLA_Linear_solver", 0 );
        aParameterLists.set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        aParameterLists.set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
        aParameterLists.set( "NLA_max_iter", tNLA_max_iter );

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "NLA_Linear_solver", 1 );
        aParameterLists.set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        aParameterLists.set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
        aParameterLists.set( "NLA_max_iter", tNLA_max_iter );
        aParameterLists.set( "NLA_is_eigen_problem", true );

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "NLA_Linear_solver", 2 );
        aParameterLists.set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        aParameterLists.set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
        aParameterLists.set( "NLA_max_iter", tNLA_max_iter );
        aParameterLists.set( "NLA_is_eigen_problem", true );

        //------------------------------------------------------------------------------

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();    // 1: thermal subproblem
        aParameterLists.set( "NLA_Nonlinear_solver_algorithms", "0,1,2" );         // set nonlinear algorithm with index 0
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists.set( "NLA_DofTypes", "TEMP" );

        // ----------------------------------------------------------

        aParameterLists( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "TSA_Nonlinear_Solver", 0 );                // using NLBGS for forward problem
        aParameterLists.set( "TSA_Nonlinear_Sensitivity_Solver", 0 );    // using monlithic for sensitivity problem

        if ( tUseTimeContinuity )
        {
            aParameterLists.set( "TSA_Num_Time_Steps", tTimeSteps );
            aParameterLists.set( "TSA_Time_Frame", tMaxTime );
        }

        //------------------------------------------------------------------------------

        aParameterLists( SOL::TIME_SOLVERS ).add_parameter_list();
        aParameterLists.set( "TSA_DofTypes", "TEMP" );
        aParameterLists.set( "TSA_Initialize_Sol_Vec", "TEMP,0.0" );
        aParameterLists.set( "TSA_Output_Indices", "0" );
        aParameterLists.set( "TSA_Output_Criteria", "Output_Criterion" );

        if ( tUseTimeContinuity )
        {
            aParameterLists.set( "TSA_time_level_per_type", "TEMP,2" );
        }
        else
        {
            aParameterLists.set( "TSA_time_level_per_type", "TEMP,1" );
        }

        //------------------------------------------------------------------------------

        aParameterLists( SOL::SOLVER_WAREHOUSE ).set( "SOL_save_operator_to_matlab", "heat" );
    }

    //------------------------------------------------------------------------------
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
