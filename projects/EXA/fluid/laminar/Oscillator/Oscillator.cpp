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
    OPTParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_opt_problem_parameter_list() );

        aParameterLists.set( "is_optimization_problem", false );
        aParameterLists.set( "workflow", "STK_FEM" );
    }

    /* ------------------------------------------------------------------------ */

    void
    STKParameterList( Module_Parameter_Lists& aParameterLists )
    {

        std::string tPrefix       = moris::get_base_moris_dir();
        std::string tMeshFileName = tPrefix + "/projects/EXA/fluid/laminar/Oscillator/Oscillator.g";

        aParameterLists( 0 ).add_parameter_list( prm::create_stk_parameter_list() );;
        aParameterLists.set( "input_file", tMeshFileName );

        gLogger.set_severity_level( 0 );
    }

    /* ------------------------------------------------------------------------ */

    void
    FEMParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.hack_for_legacy_fem();
        if ( par_rank() == 0 )
        {
            std::cout << "===============================================================" << '\n';
            std::cout << "Reynolds number = " << moris_to_string( tReynolds ) << '\n';
            std::cout << "Time scale = " << moris_to_string( tTimeScale ) << '\n';
            std::cout << "Effective viscosity = " << tMu << '\n';
            std::cout << "Effective density = " << tRho << '\n';
            std::cout << "===============================================================" << '\n';
        }


        //------------------------------------------------------------------------------
        // fill the property part of the parameter list

        // solid properties ------------------------------------------------------------

        aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropDensity" );
        aParameterLists.set( "function_parameters", tRho );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropViscosity" );
        aParameterLists.set( "function_parameters", tMu );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropTimeContinuity" );
        aParameterLists.set( "function_parameters", moris_to_string( tTimePenalty * std::stod( tRho ) ) );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropInitialVel" );
        aParameterLists.set( "function_parameters", "0.0;0.0" );

        // BC properties ---------------------------------------------------------------

        // create inlet velocity property
        aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropInletVel" );
        aParameterLists.set( "value_function", "Func_InletVel" );
        aParameterLists.set( "function_parameters", "1.0" );

        // create no-slip velocity property
        aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropNoSlip" );
        aParameterLists.set( "function_parameters", "0.0;0.0" );

        // create outlet pressure property
        aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropBackFlowPrevention" );
        aParameterLists.set( "function_parameters", "1.0" );

        //------------------------------------------------------------------------------
        // fill the constitutive model part of the parameter list

        // create fluid constitutive model
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists.set( "constitutive_name", "CMFluid" );
        aParameterLists.set( "constitutive_type", fem::Constitutive_Type::FLUID_INCOMPRESSIBLE );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "VX,VY;P", "Velocity,Pressure" ) );
        aParameterLists.set( "properties",
                "PropViscosity,Viscosity;"
                "PropDensity,Density" );

        //------------------------------------------------------------------------------
        // fill the stabilization parameter part of the parameter list

        // create SUPG stabilization parameter for Navier-Stokes
        aParameterLists( FEM::STABILIZATION ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists.set( "stabilization_name", "SPSUPGNS" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::INCOMPRESSIBLE_FLOW );
        aParameterLists.set( "function_parameters", "36.0 / 1.0" );
        aParameterLists.set( "leader_properties", "PropViscosity,Viscosity;PropDensity,Density" );
        aParameterLists.set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY;P", "Velocity,Pressure" ) );

        // create Nitsche stabilization parameter for velocity
        aParameterLists( FEM::STABILIZATION ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists.set( "stabilization_name", "SPNitscheU" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::VELOCITY_DIRICHLET_NITSCHE );
        aParameterLists.set( "function_parameters", tNitsche + "/1.0" );
        aParameterLists.set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
        aParameterLists.set( "leader_properties", "PropViscosity,Viscosity;PropDensity,Density" );

        //------------------------------------------------------------------------------
        // fill the IWG part of the parameter list

        //  bulk IWGs -------------------------------------------------------------
        /* For pure fluid phase */
        // incompressible NS velocity bulk IWG
        aParameterLists( FEM::IWG ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGVelocityBulk" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_BULK );
        aParameterLists.set( "dof_residual", "VX,VY" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY;P" );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists.set( "stabilization_parameters", "SPSUPGNS,IncompressibleFlow" );
        aParameterLists.set( "mesh_set_names", mesh_sets );

        // incompressible NS pressure bulk IWG
        aParameterLists( FEM::IWG ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGPressureBulk" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_BULK );
        aParameterLists.set( "dof_residual", "P" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY;P" );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists.set( "stabilization_parameters", "SPSUPGNS,IncompressibleFlow" );
        aParameterLists.set( "mesh_set_names", mesh_sets );

        // Walls IWGs -----------------------------------------------i

        // velocity Dirichlet IWG for interface
        aParameterLists( FEM::IWG ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGInterfaceFluidSolidVelocity" );
        aParameterLists.set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "VX,VY" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY;P" );
        aParameterLists.set( "leader_properties", "PropNoSlip,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheU,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", wall_sets );

        // pressure contribution from velocity Dirichlet IWG for interface
        aParameterLists( FEM::IWG ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGInterfaceFluidSolidPressure" );
        aParameterLists.set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "P" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY;P" );
        aParameterLists.set( "leader_properties", "PropNoSlip,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists.set( "mesh_set_names", wall_sets );

        // Inlet BC IWG ----------------------------------------------------------------

        // inlet velocity
        aParameterLists( FEM::IWG ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGInletVelocity" );
        aParameterLists.set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "VX,VY" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY;P" );
        aParameterLists.set( "leader_properties", "PropInletVel,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheU,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", inlet_sets );

        // pressure contribution from inlet velocity
        aParameterLists( FEM::IWG ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGInletPressure" );
        aParameterLists.set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "P" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY;P" );
        aParameterLists.set( "leader_properties", "PropInletVel,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists.set( "mesh_set_names", inlet_sets );

        // Outlet BC IWG ----------------------------------------------------------------

        // create inlet backflow prevention BC
        aParameterLists( FEM::IWG ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGBackflowPrevention" );
        aParameterLists.set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_IMPOSED_PRESSURE );
        aParameterLists.set( "dof_residual", "VX,VY" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY;P" );
        aParameterLists.set( "leader_properties", "PropDensity,Density;PropBackFlowPrevention,BackFlowPrevention" );
        aParameterLists.set( "mesh_set_names", outlet_sets );

        // Time continuity ----------------------------------------------------------------
        if ( tUseTimeContinuity )
        {
            // Time continuity
            aParameterLists( FEM::IWG ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists.set( "IWG_name", "IWGTimeContinuity" );
            aParameterLists.set( "IWG_type", ( fem::IWG_Type::TIME_CONTINUITY_DOF ) );
            aParameterLists.set( "dof_residual", "VX,VY" );
            aParameterLists.set( "leader_dof_dependencies", "VX,VY;P" );
            aParameterLists.set( "leader_properties",
                    "PropTimeContinuity,WeightCurrent;"
                    "PropTimeContinuity,WeightPrevious;"
                    "PropInitialVel,InitialCondition" );
            aParameterLists.set( "time_continuity", true );
            aParameterLists.set( "mesh_set_names", mesh_sets );
            }

        //------------------------------------------------------------------------------
        // fill the IQI part of the parameter list

        // x velocity
        aParameterLists( FEM::IQI ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists.set( "IQI_name", "IQIBulkVX" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "VX,VY" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY;P" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "mesh_set_names", mesh_sets );

        // y velocity
        aParameterLists( FEM::IQI ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists.set( "IQI_name", "IQIBulkVY" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "VX,VY" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY;P" );
        aParameterLists.set( "vectorial_field_index", 1 );
        aParameterLists.set( "mesh_set_names", mesh_sets );

        // pressure
        aParameterLists( FEM::IQI ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists.set( "IQI_name", "IQIBulkP" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "P" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY;P" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "mesh_set_names", mesh_sets );

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        aParameterLists( FEM::COMPUTATION ).add_parameter_list( prm::create_computation_parameter_list() );
    }

    void
    SOLParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL ) );

#ifdef MORIS_USE_MUMPS
        aParameterLists.set( "Solver_Type", "Amesos_Mumps" );
#else
        aParameterLists.set( "Solver_Type", "Amesos_Superludist" );
#endif
        //------------------------------------------------------------------------------

        aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list( moris::prm::create_linear_solver_parameter_list() );
        aParameterLists.set( "DLA_Linear_solver_algorithms", "0" );

        //------------------------------------------------------------------------------

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );    // nonlinear algorithm index 0
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists.set( "NLA_Linear_solver", 0 );
        aParameterLists.set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        aParameterLists.set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
        aParameterLists.set( "NLA_max_iter", tNLA_max_iter );

        //------------------------------------------------------------------------------

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );    // 1: thermal subproblem
        aParameterLists.set( "NLA_Nonlinear_solver_algorithms", "0" );             // set nonlinear algorithm with index 0
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists.set( "NLA_DofTypes", "VX,VY;P" );

        // ----------------------------------------------------------

        aParameterLists( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list( moris::prm::create_time_solver_algorithm_parameter_list() );
        aParameterLists.set( "TSA_Nonlinear_Solver", 0 );
        aParameterLists.set( "TSA_Nonlinear_Sensitivity_Solver", 0 );
        if ( tUseTimeContinuity )
        {
            aParameterLists.set( "TSA_Num_Time_Steps", tTimeSteps );
            aParameterLists.set( "TSA_Time_Frame", tMaxTime );
        }

        //------------------------------------------------------------------------------

        aParameterLists( SOL::TIME_SOLVERS ).add_parameter_list( moris::prm::create_time_solver_parameter_list() );
        aParameterLists.set( "TSA_DofTypes", "VX,VY;P" );
        //        aParameterLists.set( "TSA_Initialize_Sol_Vec", "VX,0.0;VY,0.0;P,0.0" );
        aParameterLists.set( "TSA_Initialize_Sol_Vec", "VX," + moris_to_string( tVelIn ) + ";VY,0.0;P,0.0" );
        aParameterLists.set( "TSA_Output_Indices", "0" );
        aParameterLists.set( "TSA_Output_Criteria", "Output_Criterion" );

        if ( tUseTimeContinuity )
        {
            aParameterLists.set( "TSA_time_level_per_type", "VX,2;VY,2;P,2" );
        }
        else
        {
            aParameterLists.set( "TSA_time_level_per_type", "VX,1;VY,1;P,1" );
        }

        //------------------------------------------------------------------------------

        aParameterLists( SOL::SOLVER_WAREHOUSE ).add_parameter_list( moris::prm::create_solver_warehouse_parameterlist() );

        aParameterLists( SOL::PRECONDITIONERS ).add_parameter_list( moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::NONE ) );
    }

    void
    MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_msi_parameter_list() );
    }

    void
    VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_vis_parameter_list() );
        aParameterLists.set( "File_Name", std::pair< std::string, std::string >( "./", tExoFile ) );
        aParameterLists.set( "Mesh_Type", vis::VIS_Mesh_Type::STANDARD );
        aParameterLists.set( "Set_Names", mesh_sets );
        aParameterLists.set( "Field_Names", "VX,VY,P" );
        aParameterLists.set( "Field_Type", "NODAL,NODAL,NODAL" );
        aParameterLists.set( "IQI_Names", "IQIBulkVX,IQIBulkVY,IQIBulkP" );
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
