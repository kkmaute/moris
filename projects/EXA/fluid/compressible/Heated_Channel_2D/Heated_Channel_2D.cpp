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
#endif
//------------------------------------------------------------------------------
namespace moris
{
    // problem size
    moris::real tChannelLength = 100.0; // in meters
    moris::real tChannelHeight = 1.0;   // in meters

    // mesh
    moris::uint tIpOrder = 2;     // polynomial order for interpolation
    moris::uint tNumXElems = 100; // number of elements in x-direction
    moris::uint tNumYElems = 1;   // number of elements in y-direction

    // Material Parameters
    moris::real tViscosity    = 1.716e-5; /* Dynamic Viscosity mu () */
    moris::real tHeatCapacity = 0.718e+3; /* Heat Capacity Cv () */
    moris::real tGasConstant  = 287.058;  /* Specific Gas Constant R () */
    moris::real tConductivity = 24.35e-3; /* Thermal Conductivity kappa() */

    //------------------------------------------------------------------------------

    // Heat Loading
    moris::real tPeakHeatLoad = 1.0e5; /* maximum volumetric heat load in center of channel in W/m^3 */
    moris::real tQfac = 200.0;         /* exponential decay factor for heat load distribution */

    // Initial Conditions
    moris::real tInitialTemperature  = 273.0;  /* Initial temperature 273 K */
    moris::real tInitialDensity      = 1.0;    /* Initial density 1.0 kg/m^3 */
    moris::real tInitialPressure     = tGasConstant * tInitialDensity * tInitialTemperature;  /* Initial pressure */

    // Nitsche Penalty for Dirichlet BCs
    moris::real tNitscheGammaP    = 0.0;    /* Penalty for pressure Dirichlet BCs */
    moris::real tNitscheGammaVX   = 100.0;  /* Penalty for x-vel Dirichlet BCs */
    moris::real tNitscheGammaVY   = 100.0;  /* Penalty for y-vel Dirichlet BCs */
    moris::real tNitscheGammaT    = 100.0;  /* Penalty for temperature Dirichlet BCs */

    // Time Continuity Weights
    moris::real tTimePenalty = 1000.0;

    // related parameters
    moris::real tGamma = 1.0 + tGasConstant / tHeatCapacity;
    moris::real tCs = std::sqrt( tGamma * tGasConstant * tInitialTemperature );

    //------------------------------------------------------------------------------

    // transient configuration
    int tNumTimeSteps = 10; // number of elements in time dimension
    moris::real tTimeStepSize = 20.0 * tChannelLength / 100.0 / tCs / 2.0;
    moris::real tTimeFrame = (real) tNumTimeSteps * tTimeStepSize;  // resulting duration
    moris::real tTCWeight = tTimePenalty / tTimeStepSize;

    // Newton configuration
    moris::real tNewtonRelaxation = 1.0;
    moris::real tNewtonTolerance = 1.0e-9;
    int tMaxNewtonSteps = 10;

    // stabilization
    bool tHaveGLS = false;

    // BC configuration
    bool tUseUpwindForPressureBC = true;
    bool tHaveFixedEnds = false; // close off ends of channel, impose zero velocity
    bool tHaveTopBottomBCs = false;

    // switch for time continuity  - for debugging
    bool tHaveTimeContinuity = true;

    // switch for bulk contribution  - for debugging
    bool tHaveBulk = true;

    // write jacobian to file - for debugging
    bool tWriteJacToMatlab = false;

    // order DoFs in global system by host - for debugging
    bool tOrderAdofsByHost = false;

    // use only Lagrange shape functions (and not B-Splines)
    bool tUseLagrange = true;

    //------------------------------------------------------------------------------

    // convert Nitsche penalty to string
    std::string sNitscheGammas = 
            ios::stringify( tNitscheGammaP ) + ";" +
            ios::stringify( tNitscheGammaVX ) + ";" +
            ios::stringify( tNitscheGammaVY ) + ";" +
            ios::stringify( tNitscheGammaT  );

    // property string, do not change
    std::string tPropertyString = "PropViscosity,DynamicViscosity;PropConductivity,ThermalConductivity";

    // bulk phase
    std::string sFluid = "HMR_dummy_c_p0,HMR_dummy_n_p0";

    // Constant function for properties
    void Func_Const(
            moris::Matrix< moris::DDRMat >                 & aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
            moris::fem::Field_Interpolator_Manager         * aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    //------------------------------------------------------------------------------

    // heat load distribution
	void Func_Heat_Load_Distribution(
            moris::Matrix< moris::DDRMat >                 & aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
            moris::fem::Field_Interpolator_Manager         * aFIManager )
    {
        // get element size in x-direction
        real tHx = tChannelLength / tNumXElems;

        // get x-coordinate
        real tX = aFIManager->get_IP_geometry_interpolator()->valx()( 0 );

        // get element number in x-direction
        real tElemNumber = std::floor( tX / tHx ); 

        // get element left and right node positions
        real tLeftNodePos = ( tElemNumber - 1.0 ) * tHx;
        real tRightNodePos = tElemNumber * tHx;

        // compute heat loads at left and right nodes
        real tQl = std::exp( -1.0 * tQfac * std::pow( 2.0 *  tLeftNodePos / tChannelLength - 1.0, 2.0 ) );
        real tQr = std::exp( -1.0 * tQfac * std::pow( 2.0 * tRightNodePos / tChannelLength - 1.0, 2.0 ) );

        // compute linear interpolation
        real tAlpha = ( ( tX - tLeftNodePos ) / ( tRightNodePos - tLeftNodePos ) );
        real tQ = tQl + tAlpha * ( tQr - tQl );
        //real tQ = std::exp( -1.0 * tQfac * std::pow( 2.0 * tX / tChannelLength - 1.0, 2.0 ) );

        // return value
        aPropMatrix = tQ * aParameters( 0 );
    }

    //------------------------------------------------------------------------------

    // initial pressure
	void Func_Initial_Pressure(
            moris::Matrix< moris::DDRMat >                 & aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
            moris::fem::Field_Interpolator_Manager         * aFIManager )
    {
        aPropMatrix = { { tInitialPressure } };
    }

    // initial velocity
	void Func_Initial_Velocity(
            moris::Matrix< moris::DDRMat >                 & aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
            moris::fem::Field_Interpolator_Manager         * aFIManager )
    {
        aPropMatrix = { { 0.0 } , { 0.0 } };
    }

    // initial pressure
	void Func_Initial_Temperature(
            moris::Matrix< moris::DDRMat >                 & aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
            moris::fem::Field_Interpolator_Manager         * aFIManager )
    {
        aPropMatrix = { { tInitialTemperature } };
    }

    //------------------------------------------------------------------------------

    // local mach Number
	void Func_Mach_Number(
            moris::Matrix< moris::DDRMat >                 & aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
            moris::fem::Field_Interpolator_Manager         * aFIManager )
    {
        // get field values
        real tTemp = aFIManager->get_field_interpolators_for_type( MSI::Dof_Type::TEMP )->val()( 0 );
        real tVx   = aFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX   )->val()( 0 );
        real tVy   = aFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX   )->val()( 1 );
        real tAbsVel = std::sqrt( tVx * tVx + tVy * tVy );

        // get ratio of specific heats
        real tKappa = ( tHeatCapacity + tGasConstant ) / tHeatCapacity;

        // compute and return Mach number
        real tMa = tAbsVel / std::sqrt( tKappa * tGasConstant * tTemp );
        aPropMatrix = { { tMa } };
    }

    // local reynolds Number
	void Func_Reynolds_Number(
            moris::Matrix< moris::DDRMat >                 & aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
            moris::fem::Field_Interpolator_Manager         * aFIManager )
    {
        // get field values
        real tTemp = aFIManager->get_field_interpolators_for_type( MSI::Dof_Type::TEMP )->val()( 0 );
        real tPres = aFIManager->get_field_interpolators_for_type( MSI::Dof_Type::P    )->val()( 0 );
        real tVx   = aFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX   )->val()( 0 );
        real tVy   = aFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX   )->val()( 1 );
        real tAbsVel = std::sqrt( tVx * tVx + tVy * tVy );

        // compute local density
        real tRho = tPres / ( tGasConstant * tTemp );

        // compute and return Reynolds number
        real tRe = tRho * tAbsVel * tChannelHeight / tViscosity;
        aPropMatrix = { { tRe } };
    }

    //------------------------------------------------------------------------------

    // Output criterion function
    bool Output_Criterion( moris::tsa::Time_Solver * aTimeSolver )
    {
        return true;
    }

    moris::real Func_Dummy_Plane(
            const moris::Matrix< DDRMat >     & aCoordinates,
            const moris::Cell< moris::real* > & aGeometryParameters )
    {
        moris::real aReturnValue = aCoordinates( 1 ) - 10000; //tPlaneBottom - 0.01;
        return aReturnValue;
    }

    //------------------------------------------------------------------------------

    void OPTParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );
        tParameterlist( 0 )( 0 ) = prm::create_opt_problem_parameter_list();
        tParameterlist( 0 )( 0 ).set( "is_optimization_problem", false);
    }

    //------------------------------------------------------------------------------

    void HMRParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
// print values
std::cout << "Time step size: " << tTimeStepSize << " \n" << std::flush;
std::cout << "Time continuity weight: " << tTCWeight << " \n" << std::flush;

        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_hmr_parameter_list();

        tParameterlist( 0 )( 0 ).set( "number_of_elements_per_dimension", ios::stringify( tNumXElems ) + "," + ios::stringify( tNumYElems ) );
        tParameterlist( 0 )( 0 ).set( "domain_dimensions",                ios::stringify( tChannelLength ) + "," + ios::stringify( tChannelHeight ) );
        tParameterlist( 0 )( 0 ).set( "domain_offset",                    "0.0," + ios::stringify( tChannelHeight / -2.0 ) );
        tParameterlist( 0 )( 0 ).set( "domain_sidesets",                  "1,2,3,4");
        tParameterlist( 0 )( 0 ).set( "lagrange_output_meshes",           "0");

        tParameterlist( 0 )( 0 ).set( "lagrange_orders",  ios::stringify( tIpOrder ) );
        tParameterlist( 0 )( 0 ).set( "lagrange_pattern", "0" );
        tParameterlist( 0 )( 0 ).set( "bspline_orders",   ios::stringify( tIpOrder ) );
        tParameterlist( 0 )( 0 ).set( "bspline_pattern",  "0" );
        if ( tUseLagrange )
        {
            tParameterlist( 0 )( 0 ).set( "lagrange_to_bspline", "-1") ;
        }

        tParameterlist( 0 )( 0 ).set( "truncate_bsplines",  1 );
        tParameterlist( 0 )( 0 ).set( "use_number_aura", 1);
        tParameterlist( 0 )( 0 ).set( "use_multigrid",  0 );
        tParameterlist( 0 )( 0 ).set( "severity_level", 0 );
    }

    //------------------------------------------------------------------------------

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
        tParameterlist( 0 )( 0 ).set( "ghost_stab",                false );
        tParameterlist( 0 )( 0 ).set( "multigrid",                 false );
        tParameterlist( 0 )( 0 ).set( "print_enriched_ig_mesh",    true );
        tParameterlist( 0 )( 0 ).set( "exodus_output_XTK_ig_mesh", true );
    }

    //------------------------------------------------------------------------------

    void GENParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 3 );
        tParameterlist( 0 ).resize( 1 );

        // Main GEN parameter list
        tParameterlist( 0 )( 0 ) = prm::create_gen_parameter_list();

        // init geometry counter
        uint tGeoCounter = 0;

        // Dummy plane
        tParameterlist( 1 ).push_back( prm::create_user_defined_geometry_parameter_list() );
        tParameterlist( 1 )( tGeoCounter ).set( "field_function_name", "Func_Dummy_Plane");
    }

    //------------------------------------------------------------------------------

    void FEMParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterList )
    {
        // create a cell of cell of parameter list for fem
        tParameterList.resize( 8 );

        // create a cell of cell of parameter list for fem
        tParameterList.resize( 9 );
        uint tPropIndex  = 0;
        uint tCMIndex    = 1;
        uint tSPIndex    = 2;
        uint tIWGIndex   = 3;
        uint tIQIIndex   = 4;
        uint tFEMIndex   = 5;
        uint tPhaseIndex = 7;
        uint tMMIndex    = 8;

        //------------------------------------------------------------------------------
        // phase info
        uint tPhaseCounter = 0;

        tParameterList( tPhaseIndex ).push_back( prm::create_phase_parameter_list() );
        tParameterList( tPhaseIndex )( tPhaseCounter ).set( "phase_name",       "PhaseFluid" );
        tParameterList( tPhaseIndex )( tPhaseCounter ).set( "phase_indices",    "0"  );
        tPhaseCounter++;

        //------------------------------------------------------------------------------
        // fill the property part of the parameter list

        // init property counter
        uint tPropCounter = 0;

        // Dynamic Viscosity mu
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name",            "PropViscosity") ;
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters",      ios::stringify( tViscosity ) );
        tParameterList( tPropIndex )( tPropCounter ).set( "value_function",           "Func_Const") ;
        tPropCounter++;

        // Heat Capacity Cv
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name",            "PropHeatCapacity") ;
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters",      ios::stringify( tHeatCapacity ) );
        tParameterList( tPropIndex )( tPropCounter ).set( "value_function",           "Func_Const") ;
        tPropCounter++;

        // Specific Gas Constant R
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name",            "PropGasConstant") ;
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters",      ios::stringify( tGasConstant ) );
        tParameterList( tPropIndex )( tPropCounter ).set( "value_function",           "Func_Const") ;
        tPropCounter++;

        // Thermal Conductivity kappa
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name",            "PropConductivity") ;
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters",      ios::stringify( tConductivity ) );
        tParameterList( tPropIndex )( tPropCounter ).set( "value_function",           "Func_Const") ;
        tPropCounter++;

        // Volumetric Heat load
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name",            "PropHeatLoad") ;
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters",      ios::stringify( tPeakHeatLoad ) );
        tParameterList( tPropIndex )( tPropCounter ).set( "value_function",           "Func_Heat_Load_Distribution") ;
        tPropCounter++;

        // velocity for no-slip
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name",            "PropZeroU") ;
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters",      "0.0;0.0") ;
        tParameterList( tPropIndex )( tPropCounter ).set( "value_function",           "Func_Const") ;
        tPropCounter++;

        // Outlet pressure BCs
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name",            "PropPressureBC" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters",      std::to_string( tInitialPressure ) );
        tParameterList( tPropIndex )( tPropCounter ).set( "value_function",           "Func_Const" );
        tPropCounter++;

        // Initial Pressure
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name",            "PropInitialPressure") ;
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters",      "1.0") ;
        tParameterList( tPropIndex )( tPropCounter ).set( "value_function",           "Func_Initial_Pressure") ;
        tPropCounter++;

        // Initial Velocity
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name",            "PropInitialVelocity") ;
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters",      "1.0") ;
        tParameterList( tPropIndex )( tPropCounter ).set( "value_function",           "Func_Initial_Velocity") ;
        tPropCounter++;

        // Initial Temperature
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name",            "PropInitialTemperature") ;
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters",      "1.0") ;
        tParameterList( tPropIndex )( tPropCounter ).set( "value_function",           "Func_Initial_Temperature") ;
        tPropCounter++;

        // Mach Number
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name",            "PropMachNumber") ;
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters",      "1.0") ;
        tParameterList( tPropIndex )( tPropCounter ).set( "value_function",           "Func_Mach_Number") ;
        tPropCounter++;

        // Reynolds Number
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name",            "PropReynoldsNumber") ;
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters",      "1.0") ;
        tParameterList( tPropIndex )( tPropCounter ).set( "value_function",           "Func_Reynolds_Number") ;
        tPropCounter++;

        // create upwind weight factor
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name",            "PropUpwind") ;
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters",      "1.0") ;
        tParameterList( tPropIndex )( tPropCounter ).set( "value_function",           "Func_Const") ;
        tPropCounter++;

        // select matrix for y-direction
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name",            "PropSelectY");
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters",      "0.0,0.0;0.0,1.0");
        tPropCounter++;

        // time continuity weights        
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name",            "PropWeightCurrent" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters",      std::to_string( tTCWeight ) );
        tParameterList( tPropIndex )( tPropCounter ).set( "value_function",           "Func_Const" );
        tPropCounter++;

        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name",            "PropWeightPrevious" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters",      std::to_string( tTCWeight ) );
        tParameterList( tPropIndex )( tPropCounter ).set( "value_function",           "Func_Const" );
        tPropCounter++;

        // dummy property
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name",            "PropDummy") ;
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters",      "1.0") ;
        tParameterList( tPropIndex )( tPropCounter ).set( "value_function",           "Func_Const") ;
        tPropCounter++;

        // zero property
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name",            "PropZero") ;
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters",      "0.0") ;
        tParameterList( tPropIndex )( tPropCounter ).set( "value_function",           "Func_Const") ;
        tPropCounter++;
                
        // zero property
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name",            "PropVectorZero") ;
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters",      "0.0;0.0") ;
        tParameterList( tPropIndex )( tPropCounter ).set( "value_function",           "Func_Const") ;
        tPropCounter++;

        //------------------------------------------------------------------------------
        // fill the material model part of the parameter list

        // init CM counter
        uint tMMCounter = 0;

        // create fluid constitutive model
        tParameterList( tMMIndex ).push_back( prm::create_material_model_parameter_list() );
        tParameterList( tMMIndex )( tMMCounter ).set( "material_name", "MMFluid" );
        tParameterList( tMMIndex )( tMMCounter ).set( "phase_name",    "PhaseFluid") ;
        tParameterList( tMMIndex )( tMMCounter ).set( "material_type", (uint) fem::Material_Type::PERFECT_GAS );
        tParameterList( tMMIndex )( tMMCounter ).set( "dof_dependencies",  std::pair< std::string, std::string >( 
                                                                           "P;TEMP", "Pressure,Temperature" ) );
        tParameterList( tMMIndex )( tMMCounter ).set( "properties", "PropHeatCapacity,IsochoricHeatCapacity;"
                                                                    "PropGasConstant,SpecificGasConstant"    );
        tMMCounter++;

        //------------------------------------------------------------------------------
        // fill the constitutive model part of the parameter list

        // init CM counter
        uint tCMCounter = 0;

        // create fluid constitutive model
        tParameterList( tCMIndex ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( tCMIndex )( tCMCounter ).set( "constitutive_name", "CMFluid" );
        tParameterList( tCMIndex )( tCMCounter ).set( "phase_name",        "PhaseFluid") ;
        tParameterList( tCMIndex )( tCMCounter ).set( "constitutive_type", (uint) fem::Constitutive_Type::FLUID_COMPRESSIBLE_NEWTONIAN );
        tParameterList( tCMIndex )( tCMCounter ).set( "dof_dependencies",  std::pair< std::string, std::string >( 
                                                                           "P;VX,VY;TEMP", "Pressure,Velocity,Temperature" ) );
        tParameterList( tCMIndex )( tCMCounter ).set( "properties", "PropViscosity,DynamicViscosity;"
                                                                    "PropConductivity,ThermalConductivity"    );
        tParameterList( tCMIndex )( tCMCounter ).set( "material_model", "MMFluid,ThermodynamicMaterialModel" );
        tCMCounter++;

        //------------------------------------------------------------------------------
        // fill the stabilization parameter part of the parameter list

        // init SP counter
        uint tSPCounter = 0;

        // create NITSCHE SP
        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name",      "NitscheSP" );
        tParameterList( tSPIndex )( tSPCounter ).set( "master_phase_name",       "PhaseFluid" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type",      (uint) fem::Stabilization_Type::COMPRESSIBLE_DIRICHLET_NITSCHE );
        tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters",     sNitscheGammas );
        tParameterList( tSPIndex )( tSPCounter ).set( "master_properties",  "PropViscosity,DynamicViscosity;"
                                                                            "PropConductivity,ThermalConductivity" );
        tSPCounter++;

        // create DUMMY SP for GLS (simply has value 1.0 everywhere)
        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name",      "DummySP" );
        tParameterList( tSPIndex )( tSPCounter ).set( "master_phase_name",       "PhaseFluid" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type",      (uint) fem::Stabilization_Type::DIRICHLET_NITSCHE );
        tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters",     "1.0" );
        tParameterList( tSPIndex )( tSPCounter ).set( "master_properties",       "PropDummy,Material" );
        tSPCounter++;


        //------------------------------------------------------------------------------
        // fill the IWG part of the parameter list

        // init IWG counter
        uint tIWGCounter = 0;

        if ( tHaveBulk )
        {
            // bulk IWG
            tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name",                   "IWGBulk" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "master_phase_name",          "PhaseFluid" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type",                   (uint) fem::IWG_Type::COMPRESSIBLE_NS_BULK );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual",               "P;VX,VY;TEMP" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "master_properties",          "PropViscosity,DynamicViscosity;"
                                                                                        "PropConductivity,ThermalConductivity;"
                                                                                        "PropHeatLoad,BodyHeatLoad" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "master_material_model",      "MMFluid,FluidMM" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "master_constitutive_models", "CMFluid,FluidCM" );
            if ( tHaveGLS )
            {
                tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters",   "DummySP,GLS" );
            }
            tIWGCounter++;
        }

        // Boundary IWG for top and bottom
        // Boundary IWG inlet
        // Boundary IWG outlet
        // tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        // tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name",                   "IWGBoundaryOutlet" );
        // tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type",              (uint) fem::Element_Type::SIDESET );
        // tParameterList( tIWGIndex )( tIWGCounter ).set( "master_phase_name",          "PhaseFluid" );
        // tParameterList( tIWGIndex )( tIWGCounter ).set( "side_ordinals",              "2" );
        // tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type",                   (uint) fem::IWG_Type::COMPRESSIBLE_NS_BOUNDARY );
        // tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual",               "P;VX,VY;TEMP" );
        // tParameterList( tIWGIndex )( tIWGCounter ).set( "master_properties",          "PropInitialPressure,Pressure;"
        //                                                                               "PropViscosity,DynamicViscosity;" 
        //                                                                               "PropConductivity,ThermalConductivity" );
        // tParameterList( tIWGIndex )( tIWGCounter ).set( "master_material_model",      "MMFluid,FluidMM" );
        // tParameterList( tIWGIndex )( tIWGCounter ).set( "master_constitutive_models", "CMFluid,FluidCM" );
        // tIWGCounter++;

        // Nitsche IWG for top and bottom
        if ( tHaveTopBottomBCs )
        {
            tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name",                   "IWGNitscheSides" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type",              (uint) fem::Element_Type::SIDESET );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "master_phase_name",          "PhaseFluid" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "side_ordinals",              "1,3" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type",                   (uint) fem::IWG_Type::COMPRESSIBLE_NS_DIRICHLET_SYMMETRIC_NITSCHE );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual",               "P;VX,VY;TEMP" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "master_properties",          "PropZeroU,PrescribedVelocity;"
                                                                                          "PropSelectY,SelectVelocity;"
                                                                                          "PropViscosity,DynamicViscosity;"
                                                                                          "PropConductivity,ThermalConductivity" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "master_material_model",      "MMFluid,FluidMM" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "master_constitutive_models", "CMFluid,FluidCM" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters",   "NitscheSP,NitschePenaltyParameter" );
            tIWGCounter++;
        }

        // Nitsche IWGs for Outlets
        if ( tHaveFixedEnds )
        {
            tPropertyString = tPropertyString + ";PropZeroU,PrescribedVelocity";
        }
        else
        {
            tPropertyString = tPropertyString + ";PropInitialPressure,PrescribedDof1";

            if( tUseUpwindForPressureBC )
            {
                tPropertyString = tPropertyString + ";PropUpwind,PressureUpwind";
            }
        }

        // Nitsche IWGs for Outlets
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name",                   "IWGNitscheOutlets" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type",              (uint) fem::Element_Type::SIDESET );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "master_phase_name",          "PhaseFluid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "side_ordinals",              "2,4" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type",                   (uint) fem::IWG_Type::COMPRESSIBLE_NS_DIRICHLET_SYMMETRIC_NITSCHE );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual",               "P;VX,VY;TEMP" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "master_properties",          tPropertyString );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "master_material_model",      "MMFluid,FluidMM" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "master_constitutive_models", "CMFluid,FluidCM" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters",   "NitscheSP,NitschePenaltyParameter" );
        tIWGCounter++;

        if ( tHaveTimeContinuity )
        {
            // Time continuity for Pressure
            tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name",                   "IWGTimeContinuityPressure" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type",              ( uint ) fem::Element_Type::BULK );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "master_phase_name",          "PhaseFluid" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type",                   ( uint ) fem::IWG_Type::TIME_CONTINUITY_DOF );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual",               "P" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "master_dof_dependencies",    "P;VX,VY;TEMP" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "master_properties",          "PropWeightCurrent,WeightCurrent;"
                                                                                          "PropWeightPrevious,WeightPrevious;"
                                                                                          "PropInitialPressure,InitialCondition" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "time_continuity",            true );
            tIWGCounter++;

            // Time continuity for Velocity
            tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name",                   "IWGTimeContinuityVelocity" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type",              ( uint ) fem::Element_Type::BULK );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "master_phase_name",          "PhaseFluid" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type",                   ( uint ) fem::IWG_Type::TIME_CONTINUITY_DOF );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual",               "VX,VY" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "master_dof_dependencies",    "P;VX,VY;TEMP" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "master_properties",          "PropWeightCurrent,WeightCurrent;"
                                                                                          "PropWeightPrevious,WeightPrevious;"
                                                                                          "PropInitialVelocity,InitialCondition" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "time_continuity",            true );
            tIWGCounter++;

            // Time continuity for Temperature
            tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name",                   "IWGTimeContinuityTemp" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type",              ( uint ) fem::Element_Type::BULK );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "master_phase_name",          "PhaseFluid" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type",                   ( uint ) fem::IWG_Type::TIME_CONTINUITY_DOF );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual",               "TEMP" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "master_dof_dependencies",    "P;VX,VY;TEMP" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "master_properties",          "PropWeightCurrent,WeightCurrent;"
                                                                                          "PropWeightPrevious,WeightPrevious;"
                                                                                          "PropInitialTemperature,InitialCondition" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "time_continuity",            true );
            tIWGCounter++;
        }

        //------------------------------------------------------------------------------
        // fill the IQI part of the parameter list

        // init IQI counter
        uint tIQICounter = 0;

        // pressure
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name",               "IQIBulkP" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "master_phase_name",      "PhaseFluid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type",               (uint) fem::IQI_Type::DOF );
        tParameterList( tIQIIndex )( tIQICounter ).set( "dof_quantity",           "P" );
        tIQICounter++;

        // velocity VX
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name",               "IQIBulkVX" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "master_phase_name",      "PhaseFluid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type",               (uint) fem::IQI_Type::DOF );
        tParameterList( tIQIIndex )( tIQICounter ).set( "dof_quantity",           "VX,VY" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index",  0 );
        tIQICounter++;

        // velocity VY
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name",               "IQIBulkVY" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "master_phase_name",      "PhaseFluid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type",               (uint) fem::IQI_Type::DOF );
        tParameterList( tIQIIndex )( tIQICounter ).set( "dof_quantity",           "VX,VY" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index",  1 );
        tIQICounter++;

        // temperature
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name",               "IQIBulkTEMP" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "master_phase_name",      "PhaseFluid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type",               (uint) fem::IQI_Type::DOF );
        tParameterList( tIQIIndex )( tIQICounter ).set( "dof_quantity",           "TEMP" );
        tIQICounter++;

        // local Mach number
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name",               "IQIMachNumber" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "master_phase_name",      "PhaseFluid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type",               (uint) fem::IQI_Type::PROPERTY );
        tParameterList( tIQIIndex )( tIQICounter ).set( "master_properties",      "PropMachNumber,Property");
        tIQICounter++;

        // local Reynolds number
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name",               "IQIReynoldsNumber" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "master_phase_name",      "PhaseFluid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type",               (uint) fem::IQI_Type::PROPERTY );
        tParameterList( tIQIIndex )( tIQICounter ).set( "master_properties",      "PropReynoldsNumber,Property");
        tIQICounter++;

        // heat load distribution
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name",               "IQIHeatLoad" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "master_phase_name",      "PhaseFluid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type",               (uint) fem::IQI_Type::PROPERTY );
        tParameterList( tIQIIndex )( tIQICounter ).set( "master_properties",      "PropHeatLoad,Property");
        tIQICounter++;

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        tParameterList( tFEMIndex ).resize( 1 );
        tParameterList( tFEMIndex )( 0 ) = prm::create_computation_parameter_list();

        //tParameterList( tFEMIndex )( 0 ).set( "finite_difference_scheme",            tFEMFdScheme  );
        //tParameterList( tFEMIndex )( 0 ).set( "finite_difference_perturbation_size", tFEMFdEpsilon );
    }

    //------------------------------------------------------------------------------

    void SOLParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 7 );
        for( uint Ik = 0; Ik < 7; Ik ++ )
        {
            tParameterlist( Ik ).resize( 1 );
        }

        tParameterlist( 0 )( 0 ) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL );

        tParameterlist( 1 )( 0 ) = moris::prm::create_linear_solver_parameter_list();

        tParameterlist( 2 )( 0 ) = moris::prm::create_nonlinear_algorithm_parameter_list();
        tParameterlist( 2 )( 0 ).set("NLA_rel_res_norm_drop",    tNewtonTolerance );
        tParameterlist( 2 )( 0 ).set("NLA_relaxation_parameter", tNewtonRelaxation );
        tParameterlist( 2 )( 0 ).set("NLA_max_iter",             tMaxNewtonSteps );

        tParameterlist( 3 )( 0 ) = moris::prm::create_nonlinear_solver_parameter_list();
        tParameterlist( 3 )( 0 ).set("NLA_DofTypes", "P;VX,VY;TEMP" );

        tParameterlist( 4 )( 0 ) = moris::prm::create_time_solver_algorithm_parameter_list();
        //tParameterlist( 4 )( 0 ).set("TSA_Nonlinear_solver", 2);       
        tParameterlist( 4 )( 0 ).set("TSA_Num_Time_Steps", tNumTimeSteps );
        tParameterlist( 4 )( 0 ).set("TSA_Time_Frame",     tTimeFrame );

        tParameterlist( 5 )( 0 ) = moris::prm::create_time_solver_parameter_list();
        tParameterlist( 5 )( 0 ).set("TSA_DofTypes"           , "P;VX,VY;TEMP" );
        tParameterlist( 5 )( 0 ).set("TSA_Initialize_Sol_Vec" , "P," + ios::stringify( tInitialPressure ) + 
                                                                ";VX,0.0;VY,0.0;TEMP," + ios::stringify( tInitialTemperature ) );
        tParameterlist( 5 )( 0 ).set("TSA_Output_Indices"     , "0" );
        tParameterlist( 5 )( 0 ).set("TSA_Output_Crteria"     , "Output_Criterion" );
        tParameterlist( 5 )( 0 ).set("TSA_time_level_per_type", "P,2;VX,2;VY,2;TEMP,2" );

        tParameterlist( 6 )( 0 ) = moris::prm::create_solver_warehouse_parameterlist();
        if ( tWriteJacToMatlab )
        {
            tParameterlist( 6 )( 0 ).set( "SOL_save_operator_to_matlab", "Matlab_2e.dat" );
        }
    }

    //------------------------------------------------------------------------------

    void MSIParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_msi_parameter_list();
        tParameterlist( 0 )( 0 ).set( "order_adofs_by_host", tOrderAdofsByHost );

    }

    //------------------------------------------------------------------------------

    void VISParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_vis_parameter_list();
        tParameterlist( 0 )( 0 ).set( "File_Name"     , std::pair< std::string, std::string >( "./", "Heated_Channel_2D.exo" ) );
        tParameterlist( 0 )( 0 ).set( "Mesh_Type"     , static_cast< uint >( vis::VIS_Mesh_Type::STANDARD ) );
        tParameterlist( 0 )( 0 ).set( "Set_Names"     , sFluid );
        tParameterlist( 0 )( 0 ).set( "Field_Names"   , "P,VX,VY,TEMP,Ma,Re,Q" );
        tParameterlist( 0 )( 0 ).set( "Field_Type"    , "NODAL,NODAL,NODAL,NODAL,NODAL,NODAL,NODAL" );
        tParameterlist( 0 )( 0 ).set( "IQI_Names"     , "IQIBulkP,IQIBulkVX,IQIBulkVY,IQIBulkTEMP,IQIMachNumber,IQIReynoldsNumber,IQIHeatLoad" ) ;
        tParameterlist( 0 )( 0 ).set( "Save_Frequency", 1 );
        tParameterlist( 0 )( 0 ).set( "Time_Offset"   , 10.0 );
    }

    //------------------------------------------------------------------------------

    void MORISGENERALParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {

    }

    //------------------------------------------------------------------------------
}

//------------------------------------------------------------------------------
#ifdef  __cplusplus
}
#endif
