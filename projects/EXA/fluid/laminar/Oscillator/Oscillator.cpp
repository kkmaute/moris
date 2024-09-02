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
#include "fn_PRM_STK_Parameters.hpp"
#include "paths.hpp"

#include "fn_equal_to.hpp"

#include "AztecOO.h"

extern moris::Logger gLogger;

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

    real tLengthScale    = 1.0;
    real tVelocityScale  = 1.0;
    real tViscosityScale = 1.0;
    real tDensityScale   = 1.0;

    real tThickness = 3e-3;       // the thickness of the 2D domain is 3mm
    real tViscosity = 1.12e-3;    // the viscosity of water is 0.00112Pa*s
    real tDensity   = 999;        // the density of water is 999kg/m^3

    real tReynolds = tDensityScale * tLengthScale * tVelocityScale / tViscosityScale;     // the Reynolds number is Re = rho*L*U/mu
    real tDarcy    = 10 * std::pow( tLengthScale, 2.0 ) / std::pow( tThickness, 2.0 );    // the height-based Darcy number is Da_h = 10*(L/h)^2

    real tTimeEnd   = 0.010;                            // the final time: original: 26.0
    real tTimeStep  = 0.005;                            // time step in seconds (COMSOL uses between 0.001 and 0.005 for 1m/s)
    real tTimeRamp  = 0.25;                             // ramp velocity up over 1s
    real tTimeScale = tLengthScale / tVelocityScale;    // the time scale is L/U

    std::string tMu    = moris_to_string( ( tViscosity / tViscosityScale ) / tReynolds );                 // the effective viscosity is mu for dimensional and 1/Re for non-dimensional
    std::string tRho   = moris_to_string( tDensity / tDensityScale );                                     // the effective density is rho for dimensional and 1 for non-dimensional
    std::string tAlpha = moris_to_string( ( tViscosity / tViscosityScale ) / ( tReynolds * tDarcy ) );    // the effective impermeability is 10*mu/h^2 for dimensional and 1/(Re*Da_h) for non-dimensional

    // prescribed inletVelocity
    real tVelIn = 5e-1 / tVelocityScale;    // always 1 since velocity scale is always inlet velocity

    /* ------------------------------------------------------------------------ */
    // parameters for weak enforcement of boundary conditions and ghost stabilization

    // Nitsche penalty
    std::string tNitsche = "50.0";

    /* ------------------------------------------------------------------------ */
    // Nonlinear solver parameters

    // maximum number of Newton steps
    int tNLA_max_iter = 20;

    // required drop of residual
    moris::real tNLA_rel_res_norm_drop = 1e-3;

    // relaxation parameter
    moris::real tNLA_relaxation_parameter = 1.0;

    /* ------------------------------------------------------------------------ */
    // parameters for transient simulations

    // number of time steps and total simulation time
    moris::real tMaxTime   = tTimeEnd / tTimeScale;
    int         tTimeSteps = tTimeEnd / tTimeStep;

    // time penalty parameter
    real tTimePenalty = 1.0;

    // time step size
    moris::real tDeltaTime = tMaxTime / (real)tTimeSteps;

    // number of steps to linearly increase inlet velocity
    moris::real tRampTime = tTimeRamp / tTimeScale;

    // flag to turn on/off transient simulation
    // bool tUseTimeContinuity = tTimeSteps > 1 ? true : false;
    bool tUseTimeContinuity = true;

    /* ------------------------------------------------------------------------ */
    // File names

    std::string tName    = "Oscillator";
    std::string tExoFile = tName + ".exo";
    std::string tSoFile  = tName + ".so";

    /* ------------------------------------------------------------------------ */
    // Output Config

    std::string mesh_sets =
            "block_1,block_2,block_3,block_4,block_5,block_6,block_7,block_8,block_9,block_10,"
            "block_11,block_12,block_13,block_14,block_15,block_16,block_17,block_18,block_19,block_20";

    std::string inlet_sets  = "surface_1,surface_2";
    std::string outlet_sets = "surface_3,surface_4";
    std::string wall_sets   = "surface_5,surface_6";

    /* ------------------------------------------------------------------------ */
    // Constant function for properties
    void
    Func_Const(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    /* ------------------------------------------------------------------------ */
    // Ramping of inlet velocity
    void
    Func_InletVel(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        real tHeightChannel = 3e-3 / tLengthScale;

        // get position in space
        real tY = aFIManager->get_IP_geometry_interpolator()->valx()( 1 );

        // velocity along x direction
        real tVelFuncY = 1.0 - 4.0 / tHeightChannel / tHeightChannel * tY * tY;

        // get position in time
        moris::real tTime = aFIManager->get_IP_geometry_interpolator()->valt()( 0 );

        // time ramping function
        real tBeta     = 4;
        real tEta      = 0.5;
        real tTimeRel  = tTime / tRampTime;
        real tVelFuncT = std::max( 0.0, std::min( 1.0, ( std::tanh( tBeta * tEta ) + std::tanh( tBeta * ( tTimeRel - tEta ) ) ) / ( std::tanh( tBeta * tEta ) + std::tanh( tBeta * ( 1.0 - tEta ) ) ) ) );    // smooth ramping

        aPropMatrix = { { tVelFuncT * tVelFuncY * tVelIn }, { 0.0 } };    // PARABOLIC RAMPING
    }

    /* ------------------------------------------------------------------------ */
    // function to determine when to output results

    bool
    Output_Criterion( moris::tsa::Time_Solver* aTimeSolver )
    {
        return true;
    }

    /* ------------------------------------------------------------------------ */

    void
    OPTParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_opt_problem_parameter_list();

        tParameterlist( 0 )( 0 ).set( "is_optimization_problem", false );
        tParameterlist( 0 )( 0 ).set( "workflow", "STK_FEM" );
    }

    /* ------------------------------------------------------------------------ */

    void
    STKParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        std::string tPrefix       = moris::get_base_moris_dir();
        std::string tMeshFileName = tPrefix + "/projects/EXA/fluid/laminar/Oscillator/Oscillator.g";

        tParameterlist( 0 )( 0 ) = prm::create_stk_parameter_list();
        tParameterlist( 0 )( 0 ).set( "input_file", tMeshFileName );

        gLogger.set_severity_level( 0 );
    }

    /* ------------------------------------------------------------------------ */

    void
    FEMParameterList( Vector< Vector< Parameter_List > >& tParameterList )
    {
        if ( par_rank() == 0 )
        {
            std::cout << "===============================================================" << '\n';
            std::cout << "Reynolds number = " << moris_to_string( tReynolds ) << '\n';
            std::cout << "Time scale = " << moris_to_string( tTimeScale ) << '\n';
            std::cout << "Effective viscosity = " << tMu << '\n';
            std::cout << "Effective density = " << tRho << '\n';
            std::cout << "===============================================================" << '\n';
        }

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

        // solid properties ------------------------------------------------------------

        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropDensity" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", tRho );
        tPropCounter++;

        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropViscosity" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", tMu );
        tPropCounter++;

        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropTimeContinuity" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", moris_to_string( tTimePenalty * std::stod( tRho ) ) );
        tPropCounter++;

        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropInitialVel" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", "0.0;0.0" );
        tPropCounter++;

        // BC properties ---------------------------------------------------------------

        // create inlet velocity property
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropInletVel" );
        tParameterList( tPropIndex )( tPropCounter ).set( "value_function", "Func_InletVel" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", "1.0" );
        tPropCounter++;

        // create no-slip velocity property
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropNoSlip" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", "0.0;0.0" );
        tPropCounter++;

        // create outlet pressure property
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropBackFlowPrevention" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", "1.0" );
        tPropCounter++;

        //------------------------------------------------------------------------------
        // fill the constitutive model part of the parameter list

        // init CM counter
        uint tCMCounter = 0;

        // create fluid constitutive model
        tParameterList( tCMIndex ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( tCMIndex )( tCMCounter ).set( "constitutive_name", "CMFluid" );
        tParameterList( tCMIndex )( tCMCounter ).set( "constitutive_type", fem::Constitutive_Type::FLUID_INCOMPRESSIBLE );
        tParameterList( tCMIndex )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( "VX,VY;P", "Velocity,Pressure" ) );
        tParameterList( tCMIndex )( tCMCounter ).set( "properties",
                "PropViscosity,Viscosity;"
                "PropDensity,Density" );
        tCMCounter++;

        //------------------------------------------------------------------------------
        // fill the stabilization parameter part of the parameter list

        uint tSPCounter = 0;

        // create SUPG stabilization parameter for Navier-Stokes
        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name", "SPSUPGNS" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type", fem::Stabilization_Type::INCOMPRESSIBLE_FLOW );
        tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters", "36.0 / 1.0" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_properties", "PropViscosity,Viscosity;PropDensity,Density" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY;P", "Velocity,Pressure" ) );
        tSPCounter++;

        // create Nitsche stabilization parameter for velocity
        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name", "SPNitscheU" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type", fem::Stabilization_Type::VELOCITY_DIRICHLET_NITSCHE );
        tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters", tNitsche + "/1.0" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_properties", "PropViscosity,Viscosity;PropDensity,Density" );
        tSPCounter++;

        //------------------------------------------------------------------------------
        // fill the IWG part of the parameter list

        // init IWG counter
        uint tIWGCounter = 0;

        //  bulk IWGs -------------------------------------------------------------
        /* For pure fluid phase */
        // incompressible NS velocity bulk IWG
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGVelocityBulk" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_BULK );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "VX,VY" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_dof_dependencies", "VX,VY;P" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SPSUPGNS,IncompressibleFlow" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "mesh_set_names", mesh_sets );
        tIWGCounter++;

        // incompressible NS pressure bulk IWG
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGPressureBulk" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_BULK );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "P" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_dof_dependencies", "VX,VY;P" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SPSUPGNS,IncompressibleFlow" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "mesh_set_names", mesh_sets );
        tIWGCounter++;

        // Walls IWGs -----------------------------------------------i

        // velocity Dirichlet IWG for interface
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGInterfaceFluidSolidVelocity" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_UNSYMMETRIC_NITSCHE );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "VX,VY" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_dof_dependencies", "VX,VY;P" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_properties", "PropNoSlip,Dirichlet" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SPNitscheU,DirichletNitsche" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "mesh_set_names", wall_sets );
        tIWGCounter++;

        // pressure contribution from velocity Dirichlet IWG for interface
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGInterfaceFluidSolidPressure" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_UNSYMMETRIC_NITSCHE );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "P" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_dof_dependencies", "VX,VY;P" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_properties", "PropNoSlip,Dirichlet" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "mesh_set_names", wall_sets );
        tIWGCounter++;

        // Inlet BC IWG ----------------------------------------------------------------

        // inlet velocity
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGInletVelocity" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_UNSYMMETRIC_NITSCHE );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "VX,VY" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_dof_dependencies", "VX,VY;P" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_properties", "PropInletVel,Dirichlet" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SPNitscheU,DirichletNitsche" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "mesh_set_names", inlet_sets );
        tIWGCounter++;

        // pressure contribution from inlet velocity
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGInletPressure" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_UNSYMMETRIC_NITSCHE );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "P" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_dof_dependencies", "VX,VY;P" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_properties", "PropInletVel,Dirichlet" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "mesh_set_names", inlet_sets );
        tIWGCounter++;

        // Outlet BC IWG ----------------------------------------------------------------

        // create inlet backflow prevention BC
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGBackflowPrevention" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_IMPOSED_PRESSURE );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "VX,VY" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_dof_dependencies", "VX,VY;P" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_properties", "PropDensity,Density;PropBackFlowPrevention,BackFlowPrevention" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "mesh_set_names", outlet_sets );
        tIWGCounter++;

        // Time continuity ----------------------------------------------------------------
        if ( tUseTimeContinuity )
        {
            // Time continuity
            tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGTimeContinuity" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", ( fem::IWG_Type::TIME_CONTINUITY_DOF ) );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "VX,VY" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_dof_dependencies", "VX,VY;P" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_properties",
                    "PropTimeContinuity,WeightCurrent;"
                    "PropTimeContinuity,WeightPrevious;"
                    "PropInitialVel,InitialCondition" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "time_continuity", true );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "mesh_set_names", mesh_sets );
            tIWGCounter++;
        }

        //------------------------------------------------------------------------------
        // fill the IQI part of the parameter list

        // init IQI counter
        uint tIQICounter = 0;

        // x velocity
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIBulkVX" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", fem::IQI_Type::DOF );
        tParameterList( tIQIIndex )( tIQICounter ).set( "dof_quantity", "VX,VY" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_dof_dependencies", "VX,VY;P" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index", 0 );
        tParameterList( tIQIIndex )( tIQICounter ).set( "mesh_set_names", mesh_sets );
        tIQICounter++;

        // y velocity
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIBulkVY" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", fem::IQI_Type::DOF );
        tParameterList( tIQIIndex )( tIQICounter ).set( "dof_quantity", "VX,VY" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_dof_dependencies", "VX,VY;P" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index", 1 );
        tParameterList( tIQIIndex )( tIQICounter ).set( "mesh_set_names", mesh_sets );
        tIQICounter++;

        // pressure
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIBulkP" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", fem::IQI_Type::DOF );
        tParameterList( tIQIIndex )( tIQICounter ).set( "dof_quantity", "P" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_dof_dependencies", "VX,VY;P" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index", 0 );
        tParameterList( tIQIIndex )( tIQICounter ).set( "mesh_set_names", mesh_sets );
        tIQICounter++;

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        tParameterList( tFEMIndex ).push_back( prm::create_computation_parameter_list() );
    }

    void
    SOLParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
    {
        tParameterlist.resize( 8 );
        for ( uint Ik = 0; Ik < 8; Ik++ )
        {
            tParameterlist( Ik ).resize( 1 );
        }

        tParameterlist( 0 )( 0 ) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL );

#ifdef MORIS_USE_MUMPS
        tParameterlist( 0 )( 0 ).set( "Solver_Type", "Amesos_Mumps" );
#else
        tParameterlist( 0 )( 0 ).set( "Solver_Type", "Amesos_Superludist" );
#endif
        //------------------------------------------------------------------------------

        tParameterlist( 1 )( 0 ) = moris::prm::create_linear_solver_parameter_list();
        tParameterlist( 1 )( 0 ).set( "DLA_Linear_solver_algorithms", "0" );

        //------------------------------------------------------------------------------

        tParameterlist( 2 )( 0 ) = moris::prm::create_nonlinear_algorithm_parameter_list();    // nonlinear algorithm index 0
        tParameterlist( 2 )( 0 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        tParameterlist( 2 )( 0 ).set( "NLA_Linear_solver", 0 );
        tParameterlist( 2 )( 0 ).set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        tParameterlist( 2 )( 0 ).set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
        tParameterlist( 2 )( 0 ).set( "NLA_max_iter", tNLA_max_iter );

        //------------------------------------------------------------------------------

        tParameterlist( 3 )( 0 ) = moris::prm::create_nonlinear_solver_parameter_list();    // 1: thermal subproblem
        tParameterlist( 3 )( 0 ).set( "NLA_Nonlinear_solver_algorithms", "0" );             // set nonlinear algorithm with index 0
        tParameterlist( 3 )( 0 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        tParameterlist( 3 )( 0 ).set( "NLA_DofTypes", "VX,VY;P" );

        // ----------------------------------------------------------

        tParameterlist( 4 )( 0 ) = moris::prm::create_time_solver_algorithm_parameter_list();
        tParameterlist( 4 )( 0 ).set( "TSA_Nonlinear_Solver", 0 );
        tParameterlist( 4 )( 0 ).set( "TSA_Nonlinear_Sensitivity_Solver", 0 );
        if ( tUseTimeContinuity )
        {
            tParameterlist( 4 )( 0 ).set( "TSA_Num_Time_Steps", tTimeSteps );
            tParameterlist( 4 )( 0 ).set( "TSA_Time_Frame", tMaxTime );
        }

        //------------------------------------------------------------------------------

        tParameterlist( 5 )( 0 ) = moris::prm::create_time_solver_parameter_list();
        tParameterlist( 5 )( 0 ).set( "TSA_DofTypes", "VX,VY;P" );
        //        tParameterlist( 5 )( 0 ).set( "TSA_Initialize_Sol_Vec", "VX,0.0;VY,0.0;P,0.0" );
        tParameterlist( 5 )( 0 ).set( "TSA_Initialize_Sol_Vec", "VX," + moris_to_string( tVelIn ) + ";VY,0.0;P,0.0" );
        tParameterlist( 5 )( 0 ).set( "TSA_Output_Indices", "0" );
        tParameterlist( 5 )( 0 ).set( "TSA_Output_Criteria", "Output_Criterion" );

        if ( tUseTimeContinuity )
        {
            tParameterlist( 5 )( 0 ).set( "TSA_time_level_per_type", "VX,2;VY,2;P,2" );
        }
        else
        {
            tParameterlist( 5 )( 0 ).set( "TSA_time_level_per_type", "VX,1;VY,1;P,1" );
        }

        //------------------------------------------------------------------------------

        tParameterlist( 6 )( 0 ) = moris::prm::create_solver_warehouse_parameterlist();

        tParameterlist( 7 )( 0 ) = moris::prm::create_preconditioner_parameter_list(sol::PreconditionerType::NONE);
    }

    void
    MSIParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_msi_parameter_list();
    }

    void
    VISParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_vis_parameter_list();
        tParameterlist( 0 )( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", tExoFile ) );
        tParameterlist( 0 )( 0 ).set( "Mesh_Type", vis::VIS_Mesh_Type::STANDARD );
        tParameterlist( 0 )( 0 ).set( "Set_Names", mesh_sets );
        tParameterlist( 0 )( 0 ).set( "Field_Names", "VX,VY,P" );
        tParameterlist( 0 )( 0 ).set( "Field_Type", "NODAL,NODAL,NODAL" );
        tParameterlist( 0 )( 0 ).set( "IQI_Names", "IQIBulkVX,IQIBulkVY,IQIBulkP" );
        tParameterlist( 0 )( 0 ).set( "Save_Frequency", 1 );
        tParameterlist( 0 )( 0 ).set( "Time_Offset", 10.0 );
    }

    void
    MORISGENERALParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
    {
    }

    //------------------------------------------------------------------------------
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
