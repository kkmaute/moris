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
#include "fn_equal_to.hpp"

#include "AztecOO.h"
#include "AztecOO.h"
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include <cmath>

// global variable for interpolation order
extern uint gInterpolationOrder;


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

    uint tInterpolationOrder = 2;
    
    uint tSolverType = 3;
    
    /* ------------------------------------------------------------------------ */
    // geometry parameters

    // dimension of domain
    real tDimX = 1.00;// x-direction
    real tDimY = 1.00;// y-direction

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
    moris::real tRampTime = 10;// tMaxTime/10.0;

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
    std::string tDomainOffX = moris_to_string( 0.0001*(-3.0) );
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

    std::string tInterfaceRefinementSphere = "0";
    std::string tInterfaceRefinementWalls  = "0";

    moris::real tElementEdgeLength = tApproxEleSize / ( std::pow( 2, tDispInitialRef ) );
    moris::real tGeoShift          = 0.1 * tElementEdgeLength;

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

    std::string tSolidPhase = "HMR_dummy_n_p0,HMR_dummy_c_p0";
    std::string tVoidPhase  = "HMR_dummy_n_p1,HMR_dummy_c_p1";
    std::string tAllPhases  = tSolidPhase + "," + tVoidPhase;
    std::string tOuterSurfaceSolid = "SideSet_1_n_p0,SideSet_2_n_p0,SideSet_3_n_p0,SideSet_4_n_p0,SideSet_1_c_p0,SideSet_2_c_p0,SideSet_3_c_p0,SideSet_4_c_p0";
    std::string tOuterSurfaceVoid = "SideSet_1_n_p1,SideSet_2_n_p1,SideSet_3_n_p1,SideSet_4_n_p1,SideSet_1_c_p1,SideSet_2_c_p1,SideSet_3_c_p1,SideSet_4_c_p1";
    std::string tAllPhaseInterfaces = tAllPhases + "," + tOuterSurfaceSolid + "," + tOuterSurfaceVoid;

    /* ------------------------------------------------------------------------ */
    // geometry parameters & LS functions

    // sphere level set function
    moris::real
    Func_Sphere(
            const moris::Matrix< DDRMat >&     aCoordinates,
            const moris::Cell< real >& aGeometryParameters )
    {
        // get coordinates
        real tX = aCoordinates( 0 );

        // compute level set value
        real tLS = tX - tPlanePos ;
        
        // return the level set value
        return tLS;
    }

    /* ------------------------------------------------------------------------ */
    // geometry parameters & LS functions

    // plane function
    moris::real
    Func_Plane(
            const moris::Matrix< DDRMat >&     aCoordinates,
            const moris::Cell< real >& aGeometryParameters )
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
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        // intermeidate variables measuring the length of the solid domain
        real tA = tPlanePos; 
        real tB = 1.0; 
        
        aPropMatrix.set_size( 1, 1 , 0.0);

        // get coordinates
        real tX = aFIManager->get_IP_geometry_interpolator()->valx()(0);
        real tY = aFIManager->get_IP_geometry_interpolator()->valx()(1);

        //initailize the temprture
        real tTemp = 0; 
        for (int n= 1 ; n < 200 ; n++ )
        {
            // cast to real
            n = (real) n; 
            
            // compute the coefficient
            real tCoeff = 2/tB / std::cosh(n*M_PI*tA/tB) * (tB*tB *(2*tB + (-2*tB + (-1 + tB)* n*n *M_PI*M_PI) *std::cos(n *M_PI) ))/(std::pow(n *M_PI,3.0)) ;
            
            tTemp += tCoeff* std::sin(n*M_PI*tY/tB) * std::cosh(n*M_PI*(tX)/tB) ; 

        }
        
        aPropMatrix( 0, 0 ) =  tTemp ; 
    }

    /* ------------------------------------------------------------------------ */
    void
    Func_HeatLoad(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        // get the constant heat flux and set it as output
        aPropMatrix.set_size( 1, 1 );

        aPropMatrix( 0, 0 ) = tHeatFlx;
    }
    
    void
    Func_TempDistro(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix.set_size( 1, 1 );

        const Matrix< DDRMat > tCoord = aFIManager->get_IP_geometry_interpolator()->valx();

        // the temp distrubution is y*(1-y) applied on the right side of the domain
        aPropMatrix( 0, 0 ) = tCoord(1)*(1-tCoord(1));
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
        
        tParameterlist( 0 )( 0 ).set( "basis_function_vtk_file", "basisinhmr.vtk" );
        tParameterlist( 0 )( 0 ).set( "write_lagrange_output_mesh_to_exodus", "lagrangehmr.exo" );
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
        tParameterlist( 0 )( 0 ).set( "use_SPG_based_enrichment", true );
        tParameterlist( 0 )( 0 ).set( "basis_rank", "bspline" );
        tParameterlist( 0 )( 0 ).set( "enrich_mesh_indices", "0" );
        tParameterlist( 0 )( 0 ).set( "ghost_stab", tUseGhost );
        tParameterlist( 0 )( 0 ).set( "multigrid", false );
        tParameterlist( 0 )( 0 ).set( "verbose", false );
        tParameterlist( 0 )( 0 ).set( "high_to_low_dbl_side_sets", true );
        tParameterlist( 0 )( 0 ).set( "print_enriched_ig_mesh", true );
        tParameterlist( 0 )( 0 ).set( "exodus_output_XTK_ig_mesh", true );
        tParameterlist( 0 )( 0 ).set( "exodus_output_XTK_ip_mesh", true );
        tParameterlist( 0 )( 0 ).set( "activate_basis_agglomeration", true );  
        tParameterlist( 0 )( 0 ).set( "write_enrichment_fields", true );
        tParameterlist( 0 )( 0 ).set( "write_enrichment_fields_probe_spheres", "1.0,0.5,0.5,1.0" );
        

    }

    /* ------------------------------------------------------------------------ */

    void
    GENParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 3 );

        tParameterlist( 0 ).push_back( prm::create_gen_parameter_list() );
        tParameterlist( 0 )( 0 ).set( "IQI_types",
                "IQIInputThermalEnergy,"
                "IQIDiffusiveLower,"
                "IQIDiffusiveUpper,"
                "IQIDiffusiveFront,"
                "IQIDiffusiveBack,"
                "IQIStoredThermalEnergy,"
                "IQIInputThermalEnergySurface" );

        tParameterlist( 0 )( 0 ).set( "number_of_phases", 2 );
        tParameterlist( 0 )( 0 ).set( "phase_function_name", "get_phase_index" );
        tParameterlist( 0 )( 0 ).set( "output_mesh_file", tGENOutputFile );
        tParameterlist( 0 )( 0 ).set( "time_offset", 10.0 );
        
        // init geometry counter
        uint tGeoCounter = 0;

        // Inclusions
        tParameterlist( 1 ).push_back( prm::create_user_defined_geometry_parameter_list() );
        tParameterlist( 1 )( tGeoCounter ).set( "field_function_name", "Func_Sphere" );
        tParameterlist( 1 )( tGeoCounter ).set( "number_of_refinements", tInterfaceRefinementSphere );
        tParameterlist( 1 )( tGeoCounter ).set( "refinement_mesh_index", "0" );
        tParameterlist( 1 )( tGeoCounter ).set( "isocontour_threshold", 0.0 );
        tParameterlist( 1 )( tGeoCounter ).set( "isocontour_tolerance", 1.0e-12 );
        tParameterlist( 1 )( tGeoCounter ).set( "intersection_tolerance", 1.0e-12 );
        tParameterlist( 1 )( tGeoCounter ).set( "multilinear_intersections", true );
        tParameterlist( 1 )( tGeoCounter ).set( "discretization_mesh_index", 0 );
        tParameterlist( 1 )( tGeoCounter ).set( "discretization_lower_bound", -2.0 );
        tParameterlist( 1 )( tGeoCounter ).set( "discretization_upper_bound", 2.0 );

        tGeoCounter++;
    }

    /* ------------------------------------------------------------------------ */

    void
    FEMParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterList )
    {

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
        tParameterList( tPhaseIndex )( tPhaseCounter ).set( "phase_name", "PhaseSolid" );
        tParameterList( tPhaseIndex )( tPhaseCounter ).set( "phase_indices", "0" );
        tPhaseCounter++;

        tParameterList( tPhaseIndex ).push_back( prm::create_phase_parameter_list() );
        tParameterList( tPhaseIndex )( tPhaseCounter ).set( "phase_name", "PhaseVoid" );
        tParameterList( tPhaseIndex )( tPhaseCounter ).set( "phase_indices", "1" );
        tPhaseCounter++;

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

        // BC properties ---------------------------------------------------------------

        // create inlet temperature property
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropInletTemp" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", moris_to_string( tInTemp ) );
        tPropCounter++;
        
               // create inlet temperature property
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropInletTemp2" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", moris_to_string( 1.0 ) );
        tParameterList( tPropIndex )( tPropCounter ).set( "value_function", "Func_TempDistro" );
        tPropCounter++;
        
                 // create inlet temperature property
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropAnalSol" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", moris_to_string( 1.0 ) );
        tParameterList( tPropIndex )( tPropCounter ).set( "value_function", "Func_Anal_Sol" );
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
        tParameterList( tCMIndex )( tCMCounter ).set( "phase_name", "PhaseSolid" );
        tParameterList( tCMIndex )( tCMCounter ).set( "constitutive_type", (uint)fem::Constitutive_Type::DIFF_LIN_ISO );
        tParameterList( tCMIndex )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( tCMIndex )( tCMCounter ).set( "properties",
                "PropConductivity,Conductivity;"
                "PropDensity     ,Density;"
                "PropCapacity    ,HeatCapacity" );
        tCMCounter++;

        //------------------------------------------------------------------------------
        // fill the stabilization parameter part of the parameter list

        uint tSPCounter = 0;

        // create Nitsche for  temperature
        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name", "SPNitsche" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_phase_name", "PhaseSolid" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type", (uint)fem::Stabilization_Type::DIRICHLET_NITSCHE );
        tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters", tNitschePenThermal );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_properties", "PropConductivity,Material" );
        tSPCounter++;

        // create ghost penalty  temperature
        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name", "SPGPTemp" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_phase_name", "PhaseSolid" );
        tParameterList( tSPIndex )( tSPCounter ).set( "follower_phase_name", "PhaseSolid" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type", (uint)fem::Stabilization_Type::GHOST_DISPL );
        tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters", "0.05" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_properties", "PropConductivity,Material" );
        tSPCounter++;

        //------------------------------------------------------------------------------
        // fill the IWG part of the parameter list

        // init IWG counter
        uint tIWGCounter = 0;

        //  bulk IWGs -------------------------------------------------------------
        // diffusion
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGDiffusionBulk" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseSolid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", (uint)fem::IWG_Type::SPATIALDIFF_BULK );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "TEMP" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_properties", "PropVolumetricHeatFlux,Load" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", "CMDiffusion,Diffusion" );
        tIWGCounter++;

        // Inlet BC IWG ----------------------------------------------------------------

        //   
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );         
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGInletTempTopBottom" );         
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", (uint)fem::Element_Type::SIDESET );         
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseSolid" );         
        tParameterList( tIWGIndex )( tIWGCounter ).set( "side_ordinals", "1,3" );         
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", (uint)fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );         
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "TEMP" );         
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_properties", "PropInletTemp,Dirichlet" );         
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", "CMDiffusion,Diffusion" );         
        tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SPNitsche,DirichletNitsche" );         
        tIWGCounter++;    
        
        //--------------------------------------------------------------------------------
        
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );        
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGInletTempRight" );         
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", (uint)fem::Element_Type::SIDESET );         
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseSolid" );        
        tParameterList( tIWGIndex )( tIWGCounter ).set( "neighbor_phases", "PhaseVoid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", (uint)fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );       
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "TEMP" );        
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_properties", "PropInletTemp2,Dirichlet" );       
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", "CMDiffusion,Diffusion" );         
        tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SPNitsche,DirichletNitsche" );         
        tIWGCounter++;          
        
        // Outlet BC IWG ----------------------------------------------------------------         


        // Time continuity ----------------------------------------------------------------

        if ( tUseTimeContinuity )
        {
            // Time continuity temperature
            tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGTimeContinuity" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", (uint)( fem::IWG_Type::TIME_CONTINUITY_DOF ) );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "TEMP" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_properties",
                    "PropTimeContinuity,WeightCurrent;"
                    "PropTimeContinuity,WeightPrevious;"
                    "PropInitialTemp,InitialCondition" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseSolid" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "time_continuity", true );
            tIWGCounter++;
        }

        // Ghost  ----------------------------------------------------------------

        if ( tUseGhost )
        {
            // ghost  temperature
            tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGGPTemp" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", (uint)fem::Element_Type::DOUBLE_SIDESET );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseSolid" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "follower_phase_name", "PhaseSolid" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", (uint)fem::IWG_Type::GHOST_NORMAL_FIELD );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "TEMP" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SPGPTemp,GhostSP" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "ghost_order", (uint)tDispOrder );
            tIWGCounter++;
        }

        //------------------------------------------------------------------------------
        // fill the IQI part of the parameter list

        // init IQI counter
        uint tIQICounter = 0;
        
                       // input thermal energy
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIAnalTemp" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseSolid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", (uint)( fem::IQI_Type::PROPERTY ) );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_properties", "PropAnalSol,Property" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index", 0 );
        tIQICounter++;
        
        // temperature
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIBulkTEMP" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseSolid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", (uint)fem::IQI_Type::DOF );
        tParameterList( tIQIIndex )( tIQICounter ).set( "dof_quantity", "TEMP" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index", 0 );
        tIQICounter++;

        // input thermal energy
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIInputThermalEnergy" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseSolid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", (uint)( fem::IQI_Type::PROPERTY ) );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_properties", "PropVolumetricHeatFlux,Property" );
        tIQICounter++;
        
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIBulkL2Error" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::L2_ERROR_ANALYTIC ) );
        tParameterList( tIQIIndex )( tIQICounter ).set( "dof_quantity", "TEMP" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_properties", "PropAnalSol,L2Check" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseSolid" );
        tIQICounter++;
        
 

        // input thermal energy
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIInputThermalEnergySurface" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_bulk_type", (uint)fem::Element_Type::SIDESET );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseSolid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "side_ordinals", "2" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", (uint)( fem::IQI_Type::PROPERTY ) );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_properties", "PropSurfaceHeatFlux,Property" );
        tIQICounter++;

        //  diffusive flux
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIDiffusiveLower" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_bulk_type", (uint)fem::Element_Type::SIDESET );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseSolid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "side_ordinals", "1" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", (uint)fem::IQI_Type::THERMAL_ENERGY_DIFFUSIVE_FLUX );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_constitutive_models", "CMDiffusion,Diffusion" );
        tIQICounter++;

        //  diffusive flux
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIDiffusiveUpper" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_bulk_type", (uint)fem::Element_Type::SIDESET );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseSolid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "side_ordinals", "3" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", (uint)fem::IQI_Type::THERMAL_ENERGY_DIFFUSIVE_FLUX );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_constitutive_models", "CMDiffusion,Diffusion" );
        tIQICounter++;

        //  diffusive flux
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIDiffusiveFront" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_bulk_type", (uint)fem::Element_Type::SIDESET );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseSolid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "side_ordinals", "4" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", (uint)fem::IQI_Type::THERMAL_ENERGY_DIFFUSIVE_FLUX );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_constitutive_models", "CMDiffusion,Diffusion" );
        tIQICounter++;

        //  diffusive flux
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIDiffusiveBack" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_bulk_type", (uint)fem::Element_Type::SIDESET );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseSolid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "side_ordinals", "2" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", (uint)fem::IQI_Type::THERMAL_ENERGY_DIFFUSIVE_FLUX );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_constitutive_models", "CMDiffusion,Diffusion" );
        tIQICounter++;

        // stored thermal energy
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIStoredThermalEnergy" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseSolid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", (uint)( fem::IQI_Type::PROPERTY ) );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_properties", "PropStoredThermalEnergy,Property" );
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

        sint tVerbosity = Belos::Errors + Belos::Warnings + Belos::IterationDetails + Belos::TimingDetails + Belos::StatusTestDetails + Belos::FinalSummary;

         switch ( tSolverType )
         {   
            case 0: // Amesos - Pardiso
                tParameterlist( 0 ) ( 0 ) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL );

//                 tParameterlist( 0 )( 0 ).set( "Solver_Type", "Amesos_Mumps");
                break;

            case 1: // Amesos - Mumps
                tParameterlist( 0 )( 0 ) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL );

                tParameterlist( 0 )( 0 ).set( "Solver_Type", "Amesos_Mumps");
                break;

            case 2: // Belos
                tParameterlist( 0 )( 0 ) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::BELOS_IMPL );

                // Solver type: GMRES, Flexible GMRES, Block CG , PseudoBlockCG, Stochastic CG, Recycling GMRES, Recycling CG, MINRES, LSQR, TFQMR
                //              Pseudoblock TFQMR, Seed GMRES, Seed CG
                tParameterlist( 0 )( 0 ).set( "Solver Type" ,  "Flexible GMRES" );

                // Diagnostics: Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::StatusTestDetails
                tParameterlist( 0 )( 0 ).set( "Verbosity" , tVerbosity );

                // Maximum number of blocks in Krylov factorization
                tParameterlist( 0 )( 0 ).set( "Num Blocks", 250   );

                // Block size to be used by iterative solver
                tParameterlist( 0 )( 0 ).set( "Block Size", 1   );

                // Allowable Belos solver iterations
                tParameterlist( 0 )( 0 ).set( "Maximum Iterations" , 500 );

                // Allowable Belos solver iterations
                tParameterlist( 0 )( 0 ).set( "Maximum Restarts" , 2 );

                // Convergence criteria
                tParameterlist( 0 )( 0 ).set( "Convergence Tolerance" ,  1e-9 );

                // Left or right preconditioner
                tParameterlist( 0 )( 0 ).set( "Left-right Preconditioner" , "right" );

                // Preconditioner

                // ifpack - ILU
                //tParameterlist( 0 )( 0 ).set( "ifpack_prec_type",  "ILU");
                //tParameterlist( 0 )( 0 ).set( "fact: level-of-fill",  5 );

                // ifpack - ILUT
                tParameterlist( 0 )( 0 ).set( "ifpack_prec_type",  "ILUT");
                tParameterlist( 0 )( 0 ).set( "fact: ilut level-of-fill", 5.0 );
                tParameterlist( 0 )( 0 ).set( "fact: drop tolerance", 1e-12 );

                // ifpack with direct solve
                //tParameterlist( 0 )( 0 ).set( "ifpack_prec_type",  "amesos");
                //tParameterlist( 0 )( 0 ).set( "amesos: solver type", "Amesos_Pardiso");

                // AMG with defaults for non-symmetric system
                //tParameterlist( 0 )( 0 ).set( "ml_prec_type",  "NSSA");
                //tParameterlist( 0 )( 0 ).set( "PDE equations", 3);

                break;

            case 3: // AZTEC
                tParameterlist( 0 )( 0 ) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AZTEC_IMPL );

                //options are: AZ_gmres, AZ_gmres_condnum, AZ_cg, AZ_cg_condnum, AZ_cgs, AZ_tfqmr, AZ_bicgstab
                tParameterlist( 0 )( 0 ).set( "AZ_solver" ,  AZ_gmres );

                // Allowable Aztec solver iterations
                tParameterlist( 0 )( 0 ).set( "AZ_max_iter", 500   );

                // Allowable Aztec iterative residual
                tParameterlist( 0 )( 0 ).set( "rel_residual" , 1e-16 );

                // set Az_conv -convergence criteria
                // options are AZ_r0, AZ_rhs, AZ_Anorm, AZ_noscaled, AZ_sol
                tParameterlist( 0 )( 0 ).set( "AZ_conv" ,  AZ_r0 );

                // set Az_diagnostic parameters
                // Set whether or not diagnostics for every linear iteration are printed or not. options are AZ_all, AZ_none
                tParameterlist( 0 )( 0 ).set( "AZ_diagnostics" ,  AZ_all );

                // set AZ_output options
                // options are AZ_all, AZ_none, AZ_warnings, AZ_last, AZ_summary
                tParameterlist( 0 )( 0 ).set( "AZ_output" ,  AZ_all );

                // Determines the submatrices factored with the domain decomposition algorithms
                // Option to specify with how many rows from other processors each processor\u2019s local submatrix is augmented.
                tParameterlist( 0 )( 0 ).set( "AZ_overlap" , 1 );

                // Determines how overlapping subdomain results are combined when different processors have computed different values for the same unknown.
                // Options are AZ_standard, AZ_symmetric
                tParameterlist( 0 )( 0 ).set( "AZ_type_overlap" , AZ_standard );

                // Determines whether RCM reordering will be done in conjunction with domain decomposition incomplete factorizations.
                // Option to enable (=1) or disable (=0) the Reverse Cuthill\u2013McKee (RCM) algorithm to reorder system equations for smaller bandwidth
                tParameterlist( 0 )( 0 ).set( "AZ_reorder" , 1 );

                // Use preconditioner from a previous Iterate() call
                // Option are AZ_calc, AZ_recalc, AZ_reuse
                tParameterlist( 0 )( 0 ).set( "AZ_pre_calc" , AZ_calc );

                // Determines  whether  matrix  factorization  information will be kept after this solve
                // for example for preconditioner_recalculation
                tParameterlist( 0 )( 0 ).set( "AZ_keep_info" , 0 );

                //--------------------------GMRES specific solver parameters--------------------------------------------------------------------------
                // Set AZ_kspace
                // Krylov subspace size for restarted GMRES
                // Setting mKrylovSpace larger improves the robustness, decreases iteration count, but increases memory consumption.
                // For very difficult problems, set it equal to the maximum number of iterations.
                tParameterlist( 0 )( 0 ).set( "AZ_kspace" ,250 );

                // Set AZ_orthog
                //AZ_classic or AZ_modified
                tParameterlist( 0 )( 0 ).set( "AZ_orthog" , AZ_classic );

                // Set AZ_rthresh
                // Parameter used to modify the relative magnitude of the diagonal entries of the matrix that is used to compute
                // any of the incomplete factorization preconditioners
                tParameterlist( 0 )( 0 ).set( "AZ_rthresh" ,  0.0 );

                // Set AZ_athresh
                // Parameter used to modify the absolute magnitude of the diagonal entries of the matrix that is used to compute
                // any of the incomplete factorization preconditioners
                tParameterlist( 0 )( 0 ).set( "AZ_athresh" ,  0.0 );

                //--------------------------Preconsitioner specific parameters--------------------------------------------------------------------------
                // Determine which preconditioner is used
                // Options are AZ_none, AZ_Jacobi, AZ_sym_GS, AZ_Neumann, AZ_ls, AZ_dom_decomp,
                tParameterlist( 0 )( 0 ).set( "AZ_precond" ,  AZ_dom_decomp );

                // Set preconditioner subdomain solve - direct solve or incomplete
                // Options are AZ_lu, AZ_ilut, AZ_ilu, AZ_rilu, AZ_bilu, AZ_icc
                tParameterlist( 0 )( 0 ).set( "AZ_subdomain_solve" ,  AZ_ilut );

                // Set preconditioner polynomial order - polynomial preconditioning, Gauss-Seidel, Jacobi
                tParameterlist( 0 )( 0 ).set( "AZ_poly_ord" ,  3 );

                // Set drop tolerance - for LU, ILUT
                tParameterlist( 0 )( 0 ).set(  "AZ_drop" ,  1.0e-12 );

                // Set level of graph fill in - for ilu(k), icc(k), bilu(k)
                tParameterlist( 0 )( 0 ).set( "AZ_graph_fill" ,  3 );

                // Set ilut fill
                tParameterlist( 0 )( 0 ).set( "AZ_ilut_fill" ,  4.0 );

                // Set Damping or relaxation parameter used for RILU
                tParameterlist( 0 )( 0 ).set( "AZ_omega" ,  1.0 );

                // Preconditioner using ifpack
                //tParameterlist( 0 )( 0 ).set( "ifpack_prec_type"    , "ILU");
                //tParameterlist( 0 )( 0 ).set( "fact: level-of-fill" ,  3       );
                //tParameterlist( 0 )( 0 ).set( "fact: drop tolerance",  1.0e-2 );
                //tParameterlist( 0 )( 0 ).set( "prec_reuse"          ,  false );

                // AMG with defaults for non-symmetric system
                tParameterlist( 0 )( 0 ).set( "ml_prec_type",  "NSSA");
                tParameterlist( 0 )( 0 ).set( "PDE equations", 1);

                break;
        }

        //------------------------------------------------------------------------------

        tParameterlist( 1 ).resize( 1 );

        tParameterlist( 1 )( 0 ) = moris::prm::create_linear_solver_parameter_list();
        tParameterlist( 1 )( 0 ).set( "DLA_Linear_solver_algorithms", "0" );

        //------------------------------------------------------------------------------

        tParameterlist( 2 ).resize( 1 );

        tParameterlist( 2 )( 0 ) = moris::prm::create_nonlinear_algorithm_parameter_list();// nonlinear algorithm index 0
        tParameterlist( 2 )( 0 ).set( "NLA_Solver_Implementation", (uint)moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        tParameterlist( 2 )( 0 ).set( "NLA_Linear_solver", 0 );
        tParameterlist( 2 )( 0 ).set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        tParameterlist( 2 )( 0 ).set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
        tParameterlist( 2 )( 0 ).set( "NLA_max_iter", tNLA_max_iter );

        //------------------------------------------------------------------------------

        tParameterlist( 3 ).resize( 1 );
        tParameterlist( 3 )( 0 ) = moris::prm::create_nonlinear_solver_parameter_list();// 1: thermal subproblem
        tParameterlist( 3 )( 0 ).set( "NLA_Nonlinear_solver_algorithms", "0" );// set nonlinear algorithm with index 0
        tParameterlist( 3 )( 0 ).set( "NLA_Solver_Implementation", (uint)moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        tParameterlist( 3 )( 0 ).set( "NLA_DofTypes", "TEMP" );

        // ----------------------------------------------------------

        tParameterlist( 4 )( 0 ) = moris::prm::create_time_solver_algorithm_parameter_list();
        tParameterlist( 4 )( 0 ).set( "TSA_Nonlinear_solver", 0 );// using NLBGS for forward problem
        tParameterlist( 4 )( 0 ).set( "TSA_nonlinear_solver_for_adjoint_solve", 0 );// using monlithic for sensitivity problem

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
        tParameterlist( 6 )( 0 ).set( "SOL_save_operator_to_matlab", "heat" );
//         tParameterlist( 6 )( 0 ).set( "SOL_TPL_Type",  (uint)( sol::MapType::Petsc ) );
    }

    void
    MSIParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_msi_parameter_list();
        tParameterlist( 0 )( 0 ).set( "TEMP", 0 );
    }

    void
    VISParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_vis_parameter_list();
        tParameterlist( 0 )( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", tExoFile ) );
        tParameterlist( 0 )( 0 ).set( "Mesh_Type", (uint)vis::VIS_Mesh_Type::STANDARD );
        tParameterlist( 0 )( 0 ).set( "Set_Names", tAllPhaseInterfaces );
        tParameterlist( 0 )( 0 ).set( "Field_Names", "TEMP,HEAT,TEMP_A,L2Nodal,L2Glob,DiffLower" );
        tParameterlist( 0 )( 0 ).set( "Field_Type", "NODAL,NODAL,NODAL,NODAL,GLOBAL,GLOBAL" );
        tParameterlist( 0 )( 0 ).set( "IQI_Names", "IQIBulkTEMP,IQIInputThermalEnergy,IQIAnalTemp,IQIBulkL2Error,IQIBulkL2Error,IQIDiffusiveLower" );
        tParameterlist( 0 )( 0 ).set( "Save_Frequency", 1 );
        tParameterlist( 0 )( 0 ).set( "Time_Offset", 10.0 );
    }

    void
    MORISGENERALParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
    }

    //------------------------------------------------------------------------------
}// namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
