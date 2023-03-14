#include <string>
#include <iostream>
#include <sstream>

#include "typedefs.hpp"
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
#include "fn_PRM_STK_Parameters.hpp"
#include "fn_equal_to.hpp"
#include "paths.hpp"

#include "AztecOO.h"

//---------------------------------------------------------------

#ifdef __cplusplus
extern "C" {
#endif
//------------------------------------------------------------------------------
namespace moris
{
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

    uint tInterpolationOrder = 1;

    /* ------------------------------------------------------------------------ */
    // geometry parameters

    // dimension of domain
    real tDimX = 1.00;    // x-direction
    real tDimY = 1.00;    // y-direction

    // position and size of sphere
    real tSpherePosX = 0;
    real tSpherePosY = 0;

    // radius of sphere
    real tSphereRadius = 0.31;
    real tSphereExpon  = 2.0;

    // interpolation of geometry with element
    bool tUseAnalyticGeometry = true;

    // approximate element size
    real tApproxEleSize = 1.0 / 4.0;

    /* ------------------------------------------------------------------------ */
    // loading parameters

    // precsribed temperature
    real tInTemp = 0.0;

    // prescribed heat flux - volumetric
    real tHeatFlx = 100.0;

    // prescribed heat flux - surface
    std::string tSurfaceHeatLoad = "100.0";

    /* ------------------------------------------------------------------------ */
    // material parameters

    std::string tDensity      = moris_to_string( 1.0 );
    std::string tCapacity     = moris_to_string( 1.0 );
    std::string tConductivity = moris_to_string( 1.916229656078740e-01 );

    /* ------------------------------------------------------------------------ */
    // parameters for weak enforcement of boundary conditions and ghost stabilization

    // Nitsche penalty
    std::string tNitschePenThermal = "100.0";

    // Flag to turn on/off ghost
    bool tUseGhost = true;

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
    int tNLA_max_iter = 20;

    // required drop of resdiual
    moris::real tNLA_rel_res_norm_drop = 1e-12;

    // relaxation parameter
    moris::real tNLA_relaxation_parameter = 1.0;

    /* ------------------------------------------------------------------------ */
    // File names

    std::string tName          = "STK_XTK_HeatConduction";
    std::string tExoFile       = tName + ".exo";
    std::string tSoFile        = tName + ".so";
    std::string tHdf5File      = tName + "_SEN.hdf5";
    std::string tGENOutputFile = tName + "_GEN.exo";

    // interpolation orders for level set and temperature fields
    int tLevelsetOrder = 2;
    int tDispOrder     = 1;

    // refinement for level set and temperature fields
    int tLevelsetInitialRef = 0;
    int tDispInitialRef     = 2;

    // refinement buffer
    int tRefineBuffer = 0;

    // refinement around sphere
    std::string tInterfaceRefinementSphere = "0";

    // automatic setup of Lagrange and Bspline discretizations
    std::string tLagrangeOrder   = std::to_string( std::max( tLevelsetOrder, tDispOrder ) );
    std::string tBsplineOrder    = std::to_string( tLevelsetOrder ) + "," + std::to_string( tDispOrder );
    std::string tInitialRef      = std::to_string( tLevelsetInitialRef ) + "," + std::to_string( tDispInitialRef );
    std::string tLagrangePattern = tLevelsetInitialRef > tDispInitialRef ? "0" : "1";

    // length of element edge
    moris::real tElementEdgeLength = tApproxEleSize / ( std::pow( 2, tDispInitialRef ) );

    // shift factor
    moris::real tGeoShift = 0.1 * tElementEdgeLength;

    /* ------------------------------------------------------------------------ */
    // Optimization parameters

    // Bspline limit
    moris::real tBSplineLimit = 5.0 * tElementEdgeLength;

    // FD in adjoint
    real tFEMFdEpsilon = 1.0e-4;
    uint tFEMFdScheme  = static_cast< uint >( fem::FDScheme_Type::POINT_3_CENTRAL );

    // FD in sweep
    std::string tSweepFdEpsilon = "1.0e-4";

    // number of constraints (here: number of design criteria)
    moris::uint tNumConstraints = 7;

    /* ------------------------------------------------------------------------ */
    // Output Config
    std::string tSolidPhase = "hmr_dummy_c_p0,hmr_dummy_n_p0";
    std::string tVoidPhase  = "hmr_dummy_c_p1,hmr_dummy_n_p1";
    std::string mesh_sets   = tSolidPhase + "," + tVoidPhase;
    std::string inlet_set   = "sideset_1_n_p0";
    std::string outlet_set  = "iside_b0_0_b1_1";
    std::string wall_sets   = "sideset_1_n_p0,sideset_2_n_p0";
    std::string tInterface  = "dbl_iside_p0_1_p1_0";

    /* ------------------------------------------------------------------------ */
    // geometry parameters & LS functions

    // sphere level set function
    moris::real
    Func_Sphere(
            const moris::Matrix< DDRMat >&     aCoordinates,
            const moris::Cell< moris::real* >& aGeometryParameters )
    {
        // get coordinates
        real tX = aCoordinates( 0 );
        real tY = aCoordinates( 1 );

        real tRadius = tUseAnalyticGeometry ? *aGeometryParameters( 0 ) : tSphereRadius;

        real tReturnValue = tRadius - std::pow( std::pow( tX - tSpherePosX, tSphereExpon ) + std::pow( tY - tSpherePosY, tSphereExpon ), 1.0 / tSphereExpon );

        return tReturnValue;
    }

    void
    Func_Sphere_Deriv(
            const moris::Matrix< moris::DDRMat >&                aCoordinates,
            const moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::Matrix< DDRMat >&                             aFieldSensitivity )
    {
        if ( tUseAnalyticGeometry )
        {
            // derivative of level set function wrt sphere radius
            aFieldSensitivity = { { 1.0 } };
        }
    }

    /* ------------------------------------------------------------------------ */
    // geometry parameters & LS functions

    // plane function
    moris::real
    Func_Plane(
            const moris::Matrix< DDRMat >&     aCoordinates,
            const moris::Cell< moris::real* >& aGeometryParameters )
    {
        // get coordinates
        real tX = aCoordinates( 0 );
        real tY = aCoordinates( 1 );

        // get normal
        real tNx = *aGeometryParameters( 0 );
        real tNy = *aGeometryParameters( 1 );

        // get point on plane
        real tPx = *aGeometryParameters( 2 );
        real tPy = *aGeometryParameters( 3 );

        real tReturnValue = tNx * ( tPx - tX ) + tNy * ( tPy - tY );

        return tReturnValue;
    }

    /* ------------------------------------------------------------------------ */
    void
    Func_HeatLoad(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix.set_size( 1, 1 );

        const Matrix< DDRMat > tCoord = aFIManager->get_IP_geometry_interpolator()->valx();

        aPropMatrix( 0, 0 ) = tHeatFlx;    // std::max(0.0,std::min(1.0,tTime/tRampTime)) * std::exp(-20*dist*dist) * tHeatFlx;
    }

    /* ------------------------------------------------------------------------ */
    // stored thermal energy

    void
    Func_StoredThermalEnergy(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
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
    OPTParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 3 );

        tParameterlist( 0 ).push_back( prm::create_opt_problem_parameter_list() );
        tParameterlist( 0 )( 0 ).set( "is_optimization_problem", tIsOptimization );
        tParameterlist( 0 )( 0 ).set( "workflow", "STK_XTK" );
        tParameterlist( 0 )( 0 ).set( "problem", "user_defined" );
        tParameterlist( 0 )( 0 ).set( "library", tSoFile );

        tParameterlist( 1 ).resize( 0 );

        tParameterlist( 2 ).push_back( moris::prm::create_sweep_parameter_list() );
        tParameterlist( 2 )( 0 ).set( "hdf5_path", tHdf5File );
        tParameterlist( 2 )( 0 ).set( "num_evaluations_per_adv", "1" );
        tParameterlist( 2 )( 0 ).set( "finite_difference_type", "all" );
        tParameterlist( 2 )( 0 ).set( "finite_difference_epsilons", tSweepFdEpsilon );
    }

    void
    STKParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        std::string tPrefix       = moris::get_base_moris_dir();
        std::string tMeshFileName = tPrefix + "/projects/EXA/structure/linear/STK_XTK_HeatConduction/STK_XTK_HeatConduction.g";

        tParameterlist( 0 )( 0 ) = prm::create_stk_parameter_list();
        tParameterlist( 0 )( 0 ).set( "input_file", tMeshFileName );

        gLogger.set_severity_level( 0 );
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
        tParameterlist( 0 )( 0 ).set( "basis_rank", "node" );
        tParameterlist( 0 )( 0 ).set( "enrich_mesh_indices", "0,1" );
        tParameterlist( 0 )( 0 ).set( "ghost_stab", tUseGhost );
        tParameterlist( 0 )( 0 ).set( "multigrid", false );
        tParameterlist( 0 )( 0 ).set( "verbose", true );
        tParameterlist( 0 )( 0 ).set( "verbose_level", moris::uint( 1 ) );
        tParameterlist( 0 )( 0 ).set( "high_to_low_dbl_side_sets", true );
        tParameterlist( 0 )( 0 ).set( "print_enriched_ig_mesh", true );
        tParameterlist( 0 )( 0 ).set( "exodus_output_XTK_ig_mesh", true );
    }

    /* ------------------------------------------------------------------------ */

    void
    GENParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 3 );

        tParameterlist( 0 ).push_back( prm::create_gen_parameter_list() );
        tParameterlist( 0 )( 0 ).set( "isocontour_threshold", 0.0 );
        tParameterlist( 0 )( 0 ).set( "isocontour_tolerance", 1.0e-12 );
        tParameterlist( 0 )( 0 ).set( "intersection_tolerance", 1.0e-12 );
        tParameterlist( 0 )( 0 ).set( "IQI_types",
                "IQIInputThermalEnergy,"
                "IQIDiffusiveLower,"
                "IQIDiffusiveUpper,"
                "IQIDiffusiveFront,"
                "IQIDiffusiveBack,"
                "IQIStoredThermalEnergy,"
                "IQIInputThermalEnergySurface" );

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
        tParameterlist( 1 )( tGeoCounter ).set( "multilinear_intersections", false );

        if ( tUseAnalyticGeometry )
        {
            tParameterlist( 0 )( 0 ).set( "initial_advs", moris_to_string( tSphereRadius ) );
            tParameterlist( 0 )( 0 ).set( "lower_bounds", moris_to_string( tSphereRadius * 0.9 ) );
            tParameterlist( 0 )( 0 ).set( "upper_bounds", moris_to_string( tSphereRadius * 1.1 ) );

            tParameterlist( 1 )( tGeoCounter ).set( "field_variable_indices", "0" );
            tParameterlist( 1 )( tGeoCounter ).set( "adv_indices", "0" );
        }
        else
        {
            tParameterlist( 1 )( tGeoCounter ).set( "discretization_mesh_index", 0 );
            tParameterlist( 1 )( tGeoCounter ).set( "discretization_lower_bound", -2.0 );
            tParameterlist( 1 )( tGeoCounter ).set( "discretization_upper_bound", 2.0 );
        }
        tGeoCounter++;
    }

    /* ------------------------------------------------------------------------ */

    void
    FEMParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterList )
    {

        // create a cell of cell of parameter list for fem
        tParameterList.resize( 8 );
        uint tPropIndex = 0;
        uint tCMIndex   = 1;
        uint tSPIndex   = 2;
        uint tIWGIndex  = 3;
        uint tIQIIndex  = 4;
        uint tFEMIndex  = 5;


        //------------------------------------------------------------------------------
        // fill the property part of the parameter list

        // init property counter
        uint tPropCounter = 0;

        // soild properties ------------------------------------------------------------

        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropDensity" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", tDensity );
        tPropCounter++;

        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropCapacity" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", tCapacity );
        tPropCounter++;

        // create  conductivity property
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropConductivity" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", tConductivity );
        tPropCounter++;

        // create  conductivity property
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropConductivityVoid" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", "10.0" );
        tPropCounter++;

        // BC properties ---------------------------------------------------------------

        // create inlet temperature property
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropInletTemp" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", moris_to_string( tInTemp ) );
        tPropCounter++;

        // create heat load property
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropVolumetricHeatFlux" );
        tParameterList( tPropIndex )( tPropCounter ).set( "value_function", "Func_HeatLoad" );
        tPropCounter++;

        // create heat load property
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropSurfaceHeatFlux" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", tSurfaceHeatLoad );
        tPropCounter++;

        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropInitialTemp" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", moris_to_string( tInTemp ) );
        tPropCounter++;

        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropTimeContinuity" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", moris_to_string( tTimePenalty * std::stod( tDensity ) * std::stod( tCapacity ) ) );
        tPropCounter++;

        // Stored Thermal Energy
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropStoredThermalEnergy" );
        tParameterList( tPropIndex )( tPropCounter ).set( "value_function", "Func_StoredThermalEnergy" );
        tPropCounter++;

        //------------------------------------------------------------------------------
        // fill the constitutive model part of the parameter list

        // init CM counter
        uint tCMCounter = 0;

        //  CM --------------------------------------------------------------------
        // create  diffusion CM
        tParameterList( tCMIndex ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( tCMIndex )( tCMCounter ).set( "constitutive_name", "CMDiffusion" );
        tParameterList( tCMIndex )( tCMCounter ).set( "constitutive_type", (uint)fem::Constitutive_Type::DIFF_LIN_ISO );
        tParameterList( tCMIndex )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( tCMIndex )( tCMCounter ).set( "properties",
                "PropConductivity,Conductivity;"
                "PropDensity     ,Density;"
                "PropCapacity    ,HeatCapacity" );
        tCMCounter++;

        tParameterList( tCMIndex ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( tCMIndex )( tCMCounter ).set( "constitutive_name", "CMDiffusionVoid" );
        tParameterList( tCMIndex )( tCMCounter ).set( "constitutive_type", (uint)fem::Constitutive_Type::DIFF_LIN_ISO );
        tParameterList( tCMIndex )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( tCMIndex )( tCMCounter ).set( "properties",
                "PropConductivityVoid,Conductivity;"
                "PropDensity     ,Density;"
                "PropCapacity    ,HeatCapacity" );
        tCMCounter++;

        //------------------------------------------------------------------------------
        // fill the stabilization parameter part of the parameter list

        uint tSPCounter = 0;

        // create Nitsche for  temperature
        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name", "SPNitscheT" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type", (uint)fem::Stabilization_Type::DIRICHLET_NITSCHE );
        tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters", tNitschePenThermal );
        tParameterList( tSPIndex )( tSPCounter ).set( "master_properties", "PropConductivity,Material" );
        tSPCounter++;

        // create parameter list for Nitsche stabilization parameter for inclusion-2 material interface
        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name", "SPInterfaceNitsche" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type", (uint)fem::Stabilization_Type::NITSCHE_INTERFACE );
        tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters", tNitschePenThermal );
        tParameterList( tSPIndex )( tSPCounter ).set( "master_properties", "PropConductivity,Material" );
        tParameterList( tSPIndex )( tSPCounter ).set( "slave_properties", "PropConductivityVoid,Material" );
        tSPCounter++;

        // create ghost penalty  temperature
        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name", "SPGPTemp" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type", (uint)fem::Stabilization_Type::GHOST_DISPL );
        tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters", "0.05" );
        tParameterList( tSPIndex )( tSPCounter ).set( "master_properties", "PropConductivity,Material" );
        tSPCounter++;

        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name", "SPGPTempVoid" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type", (uint)fem::Stabilization_Type::GHOST_DISPL );
        tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters", "0.05" );
        tParameterList( tSPIndex )( tSPCounter ).set( "master_properties", "PropConductivityVoid,Material" );
        tSPCounter++;

        //------------------------------------------------------------------------------
        // fill the IWG part of the parameter list

        // init IWG counter
        uint tIWGCounter = 0;

        //  bulk IWGs -------------------------------------------------------------
        // diffusion
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGDiffusionBulk" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", (uint)fem::IWG_Type::SPATIALDIFF_BULK );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "TEMP" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "master_dof_dependencies", "TEMP" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "master_properties", "PropVolumetricHeatFlux,Load" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "master_constitutive_models", "CMDiffusion,Diffusion" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "mesh_set_names", mesh_sets );
        tIWGCounter++;

        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGDiffusionVoid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", (uint)fem::IWG_Type::SPATIALDIFF_BULK );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "TEMP" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "master_dof_dependencies", "TEMP" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "master_properties", "PropVolumetricHeatFlux,Load" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "master_constitutive_models", "CMDiffusionVoid,Diffusion" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "mesh_set_names", tVoidPhase );
        tIWGCounter++;

        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGInterfaceSolidVoid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", (uint)fem::Element_Type::DOUBLE_SIDESET );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", (uint)fem::IWG_Type::SPATIALDIFF_INTERFACE_UNSYMMETRIC_NITSCHE );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "TEMP" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "master_dof_dependencies", "TEMP" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "master_constitutive_models", "CMDiffusion,Diffusion" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "slave_constitutive_models", "CMDiffusionVoid,Diffusion" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SPInterfaceNitsche,NitscheInterface" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "mesh_set_names", tInterface );
        tIWGCounter++;


        // Inlet BC IWG ----------------------------------------------------------------

        // inlet temperature
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGInletTemp" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", (uint)fem::Element_Type::SIDESET );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", (uint)fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "TEMP" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "master_dof_dependencies", "TEMP" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "master_properties", "PropInletTemp,Dirichlet" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "master_constitutive_models", "CMDiffusion,Diffusion" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SPNitscheT,DirichletNitsche" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "mesh_set_names", inlet_set );
        tIWGCounter++;

        // Outlet BC IWG ----------------------------------------------------------------

        // heat flux around inclusion
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGIntefaceHeat" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", (uint)fem::Element_Type::SIDESET );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", (uint)fem::IWG_Type::SPATIALDIFF_NEUMANN );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "TEMP" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "master_dof_dependencies", "TEMP" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "master_properties", "PropSurfaceHeatFlux,Neumann" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "mesh_set_names", outlet_set );
        tIWGCounter++;

        // Time continuity ----------------------------------------------------------------

        if ( tUseTimeContinuity )
        {
            // Time continuity temperature
            tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGTimeContinuity" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", (uint)( fem::IWG_Type::TIME_CONTINUITY_DOF ) );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "TEMP" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "master_dof_dependencies", "TEMP" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "master_properties",
                    "PropTimeContinuity,WeightCurrent;"
                    "PropTimeContinuity,WeightPrevious;"
                    "PropInitialTemp,InitialCondition" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "time_continuity", true );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "mesh_set_names", mesh_sets );
            tIWGCounter++;
        }

        // Ghost  ----------------------------------------------------------------

        if ( tUseGhost )
        {
            // ghost  temperature
            tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGGPTemp" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", (uint)fem::IWG_Type::GHOST_NORMAL_FIELD );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "TEMP" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "master_dof_dependencies", "TEMP" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "slave_dof_dependencies", "TEMP" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SPGPTemp,GhostSP" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "ghost_order", (uint)tDispOrder );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "mesh_set_names", "ghost_p0" );
            tIWGCounter++;

            // ghost  temperature
            tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGGPTempVoid" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", (uint)fem::Element_Type::DOUBLE_SIDESET );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", (uint)fem::IWG_Type::GHOST_NORMAL_FIELD );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "TEMP" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "master_dof_dependencies", "TEMP" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "slave_dof_dependencies", "TEMP" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SPGPTempVoid,GhostSP" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "ghost_order", (uint)tDispOrder );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "mesh_set_names", "ghost_p1" );
            tIWGCounter++;
        }

        //------------------------------------------------------------------------------
        // fill the IQI part of the parameter list

        // init IQI counter
        uint tIQICounter = 0;

        // temperature
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIBulkTEMP" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", (uint)fem::IQI_Type::DOF );
        tParameterList( tIQIIndex )( tIQICounter ).set( "dof_quantity", "TEMP" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "master_dof_dependencies", "TEMP" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index", 0 );
        tParameterList( tIQIIndex )( tIQICounter ).set( "mesh_set_names", mesh_sets );
        tIQICounter++;

        // input thermal energy
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIInputThermalEnergy" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", (uint)( fem::IQI_Type::PROPERTY ) );
        tParameterList( tIQIIndex )( tIQICounter ).set( "master_properties", "PropVolumetricHeatFlux,Property" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "mesh_set_names", tSolidPhase );
        tIQICounter++;

        // input thermal energy
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIInputThermalEnergySurface" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", (uint)( fem::IQI_Type::PROPERTY ) );
        tParameterList( tIQIIndex )( tIQICounter ).set( "master_properties", "PropSurfaceHeatFlux,Property" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "mesh_set_names", wall_sets );
        tIQICounter++;

        //  diffusive flux
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIDiffusiveLower" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_bulk_type", (uint)fem::Element_Type::SIDESET );
        tParameterList( tIQIIndex )( tIQICounter ).set( "side_ordinals", "1" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", (uint)fem::IQI_Type::THERMAL_ENERGY_DIFFUSIVE_FLUX );
        tParameterList( tIQIIndex )( tIQICounter ).set( "master_constitutive_models", "CMDiffusion,Diffusion" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "mesh_set_names", wall_sets );
        tIQICounter++;

        //  diffusive flux
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIDiffusiveUpper" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_bulk_type", (uint)fem::Element_Type::SIDESET );
        tParameterList( tIQIIndex )( tIQICounter ).set( "side_ordinals", "3" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", (uint)fem::IQI_Type::THERMAL_ENERGY_DIFFUSIVE_FLUX );
        tParameterList( tIQIIndex )( tIQICounter ).set( "master_constitutive_models", "CMDiffusion,Diffusion" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "mesh_set_names", wall_sets );
        tIQICounter++;

        //  diffusive flux
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIDiffusiveFront" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_bulk_type", (uint)fem::Element_Type::SIDESET );
        tParameterList( tIQIIndex )( tIQICounter ).set( "side_ordinals", "4" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", (uint)fem::IQI_Type::THERMAL_ENERGY_DIFFUSIVE_FLUX );
        tParameterList( tIQIIndex )( tIQICounter ).set( "master_constitutive_models", "CMDiffusion,Diffusion" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "mesh_set_names", wall_sets );
        tIQICounter++;

        //  diffusive flux
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIDiffusiveBack" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_bulk_type", (uint)fem::Element_Type::SIDESET );
        tParameterList( tIQIIndex )( tIQICounter ).set( "side_ordinals", "2" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", (uint)fem::IQI_Type::THERMAL_ENERGY_DIFFUSIVE_FLUX );
        tParameterList( tIQIIndex )( tIQICounter ).set( "master_constitutive_models", "CMDiffusion,Diffusion" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "mesh_set_names", wall_sets );
        tIQICounter++;

        // stored thermal energy
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIStoredThermalEnergy" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", (uint)( fem::IQI_Type::PROPERTY ) );
        tParameterList( tIQIIndex )( tIQICounter ).set( "master_properties", "PropStoredThermalEnergy,Property" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "mesh_set_names", wall_sets );
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
        tParameterlist.resize( 7 );
        for ( uint Ik = 0; Ik < 7; Ik++ )
        {
            tParameterlist( Ik ).resize( 1 );
        }

        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL );
#ifdef MORIS_USE_MUMPS
        tParameterlist( 0 )( 0 ).set( "Solver_Type", "Amesos_Mumps" );
#else
        tParameterlist( 0 )( 0 ).set( "Solver_Type", "Amesos_Superludist" );
#endif

        //------------------------------------------------------------------------------

        tParameterlist( 1 ).resize( 1 );

        tParameterlist( 1 )( 0 ) = moris::prm::create_linear_solver_parameter_list();
        tParameterlist( 1 )( 0 ).set( "DLA_Linear_solver_algorithms", "0" );

        //------------------------------------------------------------------------------

        tParameterlist( 2 ).resize( 1 );

        tParameterlist( 2 )( 0 ) = moris::prm::create_nonlinear_algorithm_parameter_list();    // nonlinear algorithm index 0
        tParameterlist( 2 )( 0 ).set( "NLA_Solver_Implementation", (uint)moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        tParameterlist( 2 )( 0 ).set( "NLA_Linear_solver", 0 );
        tParameterlist( 2 )( 0 ).set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        tParameterlist( 2 )( 0 ).set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
        tParameterlist( 2 )( 0 ).set( "NLA_max_iter", tNLA_max_iter );

        //------------------------------------------------------------------------------

        tParameterlist( 3 ).resize( 1 );
        tParameterlist( 3 )( 0 ) = moris::prm::create_nonlinear_solver_parameter_list();    // 1: thermal subproblem
        tParameterlist( 3 )( 0 ).set( "NLA_Nonlinear_solver_algorithms", "0" );             // set nonlinear algorithm with index 0
        tParameterlist( 3 )( 0 ).set( "NLA_Solver_Implementation", (uint)moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        tParameterlist( 3 )( 0 ).set( "NLA_DofTypes", "TEMP" );

        // ----------------------------------------------------------

        tParameterlist( 4 )( 0 ) = moris::prm::create_time_solver_algorithm_parameter_list();
        tParameterlist( 4 )( 0 ).set( "TSA_Nonlinear_solver", 0 );                      // using NLBGS for forward problem
        tParameterlist( 4 )( 0 ).set( "TSA_nonlinear_solver_for_adjoint_solve", 0 );    // using monlithic for sensitivity problem

        if ( tUseTimeContinuity )
        {
            tParameterlist( 4 )( 0 ).set( "TSA_Num_Time_Steps", tTimeSteps );
            tParameterlist( 4 )( 0 ).set( "TSA_Time_Frame", tMaxTime );
        }

        //------------------------------------------------------------------------------

        tParameterlist( 5 )( 0 ) = moris::prm::create_time_solver_parameter_list();
        tParameterlist( 5 )( 0 ).set( "TSA_DofTypes", "TEMP" );
        tParameterlist( 5 )( 0 ).set( "TSA_Initialize_Sol_Vec", "TEMP,0.0" );
        tParameterlist( 5 )( 0 ).set( "TSA_Output_Indices", "0" );
        tParameterlist( 5 )( 0 ).set( "TSA_Output_Criteria", "Output_Criterion" );

        if ( tUseTimeContinuity )
        {
            tParameterlist( 5 )( 0 ).set( "TSA_time_level_per_type", "TEMP,2" );
        }
        else
        {
            tParameterlist( 5 )( 0 ).set( "TSA_time_level_per_type", "TEMP,1" );
        }

        //------------------------------------------------------------------------------

        tParameterlist( 6 )( 0 ) = moris::prm::create_solver_warehouse_parameterlist();
    }

    void
    MSIParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_msi_parameter_list();
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
        tParameterlist( 0 )( 0 ).set( "Set_Names", mesh_sets );
        tParameterlist( 0 )( 0 ).set( "Field_Names", "TEMP,HEAT" );
        tParameterlist( 0 )( 0 ).set( "Field_Type", "NODAL,NODAL" );
        tParameterlist( 0 )( 0 ).set( "IQI_Names", "IQIBulkTEMP,IQIInputThermalEnergy" );
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
