/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Heated_Bubble_2D.cpp
 *
 */

#include <string>
#include <iostream>
#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_MSI_Equation_Object.hpp"
#include "cl_TSA_Time_Solver.hpp"
#include "cl_DLA_Solver_Interface.hpp"
#include "cl_DLA_Linear_Solver_Aztec.hpp"
#include "parameters.hpp"
#include "fn_equal_to.hpp"

#include "AztecOO.h"

#ifdef __cplusplus
extern "C" {
#endif
//------------------------------------------------------------------------------
namespace moris
{
    // problem size
    moris::real tChannelLength = 50.0;    // in meters
    moris::real tChannelHeight = 50.0;    // in meters

    // mesh
    moris::uint tIpOrder   = 2;     // polynomial order for interpolation
    moris::uint tNumXElems = 50;    // number of elements in x-direction
    moris::uint tNumYElems = 50;    // number of elements in y-direction

    // Material Parameters
    moris::real tViscosity    = 1.716e-5; /* Dynamic Viscosity mu () */
    moris::real tHeatCapacity = 0.718e+3; /* Heat Capacity Cv () */
    moris::real tGasConstant  = 287.058;  /* Specific Gas Constant R () */
    moris::real tConductivity = 24.35e-3; /* Thermal Conductivity kappa() */

    //------------------------------------------------------------------------------

    // Heat Loading
    moris::real tPeakHeatLoad = 1.0e5; /* maximum volumetric heat load in center of channel in W/m^3 */
    moris::real tQfac         = 200.0; /* exponential decay factor for heat load distribution */

    // Initial Conditions
    moris::real tInitialTemperature = 273.0;                                                /* Initial temperature 273 K */
    moris::real tInitialDensity     = 1.0;                                                  /* Initial density 1.0 kg/m^3 */
    moris::real tInitialPressure    = tGasConstant * tInitialDensity * tInitialTemperature; /* Initial pressure */

    // Nitsche Penalty for Dirichlet BCs
    moris::real tNitscheGammaP  = 0.0;        /* Penalty for pressure Dirichlet BCs */
    moris::real tNitscheGammaVX = 10000000.0; /* Penalty for x-vel Dirichlet BCs */
    moris::real tNitscheGammaVY = 10000000.0; /* Penalty for y-vel Dirichlet BCs */
    moris::real tNitscheGammaT  = 0.0;        /* Penalty for temperature Dirichlet BCs */

    // Time Continuity Weights
    moris::real tTimePenalty = 10.0;

    // related parameters
    moris::real tGamma = 1.0 + tGasConstant / tHeatCapacity;
    moris::real tCs    = std::sqrt( tGamma * tGasConstant * tInitialTemperature );

    //------------------------------------------------------------------------------

    // transient configuration
    int         tNumTimeSteps = 30;    // number of elements in time dimension
    moris::real tTimeStepSize = 20.0 * tChannelLength / 100.0 / tCs / 2.0;
    moris::real tTimeFrame    = (real)tNumTimeSteps * tTimeStepSize;    // resulting duration
    moris::real tTCWeight     = tTimePenalty / tTimeStepSize;

    // Newton configuration
    moris::real tNewtonRelaxation = 1.0;
    moris::real tNewtonTolerance  = 1.0e-7;
    int         tMaxNewtonSteps   = 10;

    // stabilization
    bool tHaveGLS = true;

    // BC configuration
    bool tUseUpwindForPressureBC = true;     // fix pressure at ends of channel
    bool tHaveFixedEnds          = false;    // close off ends of channel, impose zero velocity
    bool tHaveTopBottomBCs       = false;    // impose zero velocity in normal direction on channel sides

    // switch for time continuity - for debugging
    bool tHaveTimeContinuity = true;

    // switch for bulk contribution - for debugging
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
            ios::stringify( tNitscheGammaP ) + ";" + ios::stringify( tNitscheGammaVX ) + ";" + ios::stringify( tNitscheGammaVY ) + ";" + ios::stringify( tNitscheGammaT );

    // property string, do not change
    std::string tPropertyString = "PropViscosity,DynamicViscosity;PropConductivity,ThermalConductivity";

    // bulk phase
    std::string sFluid = "HMR_dummy_c_p0,HMR_dummy_n_p0";

    // Constant function for properties
    void
    Func_Const(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    //------------------------------------------------------------------------------

    // heat load distribution
    void
    Func_Heat_Load_Distribution(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        // get coordinates
        real tX = aFIManager->get_IP_geometry_interpolator()->valx()( 0 );
        real tY = aFIManager->get_IP_geometry_interpolator()->valx()( 1 );

        // get radius
        real tR = std::sqrt( tX * tX + tY * tY );

        // compute Heat Load
        real tQ = std::exp( -1.0 * tQfac * std::pow( 2.0 * tR / tChannelLength, 2.0 ) );

        // return value
        aPropMatrix = tQ * aParameters( 0 );
    }

    //------------------------------------------------------------------------------

    // initial pressure
    void
    Func_Initial_Pressure(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        aPropMatrix = { { tInitialPressure } };
    }

    // initial velocity
    void
    Func_Initial_Velocity(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        aPropMatrix = { { 0.0 }, { 0.0 } };
    }

    // initial pressure
    void
    Func_Initial_Temperature(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        aPropMatrix = { { tInitialTemperature } };
    }

    //------------------------------------------------------------------------------

    // local mach Number
    void
    Func_Mach_Number(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        // get field values
        real tTemp   = aFIManager->get_field_interpolators_for_type( MSI::Dof_Type::TEMP )->val()( 0 );
        real tVx     = aFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX )->val()( 0 );
        real tVy     = aFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX )->val()( 1 );
        real tAbsVel = std::sqrt( tVx * tVx + tVy * tVy );

        // get ratio of specific heats
        real tKappa = ( tHeatCapacity + tGasConstant ) / tHeatCapacity;

        // compute and return Mach number
        real tMa    = tAbsVel / std::sqrt( tKappa * tGasConstant * tTemp );
        aPropMatrix = { { tMa } };
    }

    // local reynolds Number
    void
    Func_Reynolds_Number(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        // get field values
        real tTemp   = aFIManager->get_field_interpolators_for_type( MSI::Dof_Type::TEMP )->val()( 0 );
        real tPres   = aFIManager->get_field_interpolators_for_type( MSI::Dof_Type::P )->val()( 0 );
        real tVx     = aFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX )->val()( 0 );
        real tVy     = aFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX )->val()( 1 );
        real tAbsVel = std::sqrt( tVx * tVx + tVy * tVy );

        // compute local density
        real tRho = tPres / ( tGasConstant * tTemp );

        // compute and return Reynolds number
        real tRe    = tRho * tAbsVel * tChannelHeight / tViscosity;
        aPropMatrix = { { tRe } };
    }

    //------------------------------------------------------------------------------

    // Output criterion function
    bool
    Output_Criterion( moris::tsa::Time_Solver* aTimeSolver )
    {
        return true;
    }

    moris::real
    Func_Dummy_Plane(
            const moris::Matrix< DDRMat >& aCoordinates,
            const Vector< real >&          aGeometryParameters )
    {
        moris::real aReturnValue = aCoordinates( 1 ) - 10000;    // tPlaneBottom - 0.01;
        return aReturnValue;
    }

    //------------------------------------------------------------------------------

    void
    OPTParameterList( Vector< Vector< ParameterList > >& aParameterLists )
    {
        aParameterLists( 0 ).push_back( prm::create_opt_problem_parameter_list() );
        aParameterLists.set( "is_optimization_problem", false );
    }

    //------------------------------------------------------------------------------

    void
    HMRParameterList( Vector< Vector< ParameterList > >& aParameterLists )
    {
        // print values
        std::cout << "Time step size: " << tTimeStepSize << " \n"
                  << std::flush;
        std::cout << "Time continuity weight: " << tTCWeight << " \n"
                  << std::flush;

        aParameterLists( 0 ).push_back( prm::create_hmr_parameter_list() );

        std::string tOffSet = ios::stringify( tChannelLength / -2.0 ) + "," + ios::stringify( tChannelHeight / -2.0 );

        aParameterLists.set( "number_of_elements_per_dimension", ios::stringify( tNumXElems ) + "," + ios::stringify( tNumYElems ) );
        aParameterLists.set( "domain_dimensions", ios::stringify( tChannelLength ) + "," + ios::stringify( tChannelHeight ) );
        aParameterLists.set( "domain_offset", tOffSet );
        aParameterLists.set( "lagrange_output_meshes", "0" );

        aParameterLists.set( "lagrange_orders", ios::stringify( tIpOrder ) );
        aParameterLists.set( "lagrange_pattern", "0" );
        aParameterLists.set( "bspline_orders", ios::stringify( tIpOrder ) );
        aParameterLists.set( "bspline_pattern", "0" );
        if ( tUseLagrange )
        {
            aParameterLists.set( "lagrange_to_bspline", "-1" );
        }

        aParameterLists.set( "truncate_bsplines", 1 );
        aParameterLists.set( "use_number_aura", 1 );
        aParameterLists.set( "use_multigrid", 0 );
        aParameterLists.set( "severity_level", 0 );
    }

    //------------------------------------------------------------------------------

    void
    XTKParameterList( Vector< Vector< ParameterList > >& aParameterLists )
    {
        aParameterLists( 0 ).push_back( prm::create_xtk_parameter_list() );
        aParameterLists.set( "decompose", true );
        aParameterLists.set( "decomposition_type", "conformal" );
        aParameterLists.set( "enrich", true );
        aParameterLists.set( "basis_rank", "bspline" );
        aParameterLists.set( "enrich_mesh_indices", "0" );
        aParameterLists.set( "ghost_stab", false );
        aParameterLists.set( "multigrid", false );
        aParameterLists.set( "print_enriched_ig_mesh", true );
        aParameterLists.set( "exodus_output_XTK_ig_mesh", true );
    }

    //------------------------------------------------------------------------------

    void
    GENParameterList( Vector< Vector< ParameterList > >& aParameterLists )
    {
        aParameterLists( 0 ).push_back( prm::create_gen_parameter_list() );

        // Dummy plane
        aParameterLists( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists.set( "field_function_name", "Func_Dummy_Plane" );
    }

    //------------------------------------------------------------------------------

    void
    FEMParameterList( Vector< Vector< ParameterList > >& aParameterLists )
    {
        // create a cell of cell of parameter list for fem


        //------------------------------------------------------------------------------

        aParameterLists( tPhaseIndex ).push_back( prm::create_phase_parameter_list() );
        aParameterLists.set( "phase_name", "PhaseFluid" );
        aParameterLists.set( "phase_indices", "0" );

        //------------------------------------------------------------------------------
        // fill the property part of the parameter list

        // Dynamic Viscosity mu
        aParameterLists( tPropIndex ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropViscosity" );
        aParameterLists.set( "function_parameters", ios::stringify( tViscosity ) );
        aParameterLists.set( "value_function", "Func_Const" );

        // Heat Capacity Cv
        aParameterLists( tPropIndex ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropHeatCapacity" );
        aParameterLists.set( "function_parameters", ios::stringify( tHeatCapacity ) );
        aParameterLists.set( "value_function", "Func_Const" );

        // Specific Gas Constant R
        aParameterLists( tPropIndex ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropGasConstant" );
        aParameterLists.set( "function_parameters", ios::stringify( tGasConstant ) );
        aParameterLists.set( "value_function", "Func_Const" );

        // Thermal Conductivity kappa
        aParameterLists( tPropIndex ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropConductivity" );
        aParameterLists.set( "function_parameters", ios::stringify( tConductivity ) );
        aParameterLists.set( "value_function", "Func_Const" );

        // Volumetric Heat load
        aParameterLists( tPropIndex ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropHeatLoad" );
        aParameterLists.set( "function_parameters", ios::stringify( tPeakHeatLoad ) );
        aParameterLists.set( "value_function", "Func_Heat_Load_Distribution" );    // Func_Heat_Load_Distribution

        // velocity for no-slip
        aParameterLists( tPropIndex ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropZeroU" );
        aParameterLists.set( "function_parameters", "0.0;0.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // Outlet pressure BCs
        aParameterLists( tPropIndex ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropPressureBC" );
        aParameterLists.set( "function_parameters", std::to_string( tInitialPressure ) );
        aParameterLists.set( "value_function", "Func_Const" );

        // Initial Pressure
        aParameterLists( tPropIndex ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropInitialPressure" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Initial_Pressure" );

        // Initial Velocity
        aParameterLists( tPropIndex ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropInitialVelocity" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Initial_Velocity" );

        // Initial Temperature
        aParameterLists( tPropIndex ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropInitialTemperature" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Initial_Temperature" );

        // Mach Number
        aParameterLists( tPropIndex ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropMachNumber" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Mach_Number" );

        // Reynolds Number
        aParameterLists( tPropIndex ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropReynoldsNumber" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Reynolds_Number" );

        // create upwind weight factor
        aParameterLists( tPropIndex ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropUpwind" );
        aParameterLists.set( "function_parameters", "3.774074960608205e-03" );    // 3.774074960608205e-03
        aParameterLists.set( "value_function", "Func_Const" );

        // select matrix for y-direction
        aParameterLists( tPropIndex ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropSelectY" );
        aParameterLists.set( "function_parameters", "0.0,0.0;0.0,1.0" );

        // time continuity weights
        aParameterLists( tPropIndex ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropWeightCurrent" );
        aParameterLists.set( "function_parameters", std::to_string( tTCWeight ) );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( tPropIndex ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropWeightPrevious" );
        aParameterLists.set( "function_parameters", std::to_string( tTCWeight ) );
        aParameterLists.set( "value_function", "Func_Const" );

        // dummy property
        aParameterLists( tPropIndex ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropDummy" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // zero property
        aParameterLists( tPropIndex ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropZero" );
        aParameterLists.set( "function_parameters", "0.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // zero property
        aParameterLists( tPropIndex ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropVectorZero" );
        aParameterLists.set( "function_parameters", "0.0;0.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        //------------------------------------------------------------------------------
        // fill the material model part of the parameter list

        // init CM counter
        uint tMMCounter = 0;

        // create fluid constitutive model
        aParameterLists( tMMIndex ).push_back( prm::create_material_model_parameter_list() );
        aParameterLists.set( "material_name", "MMFluid" );
        aParameterLists.set( "phase_name", "PhaseFluid" );
        aParameterLists.set( "material_type", fem::Material_Type::PERFECT_GAS );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "P;TEMP", "Pressure,Temperature" ) );
        aParameterLists.set( "properties",
                "PropHeatCapacity,IsochoricHeatCapacity;"
                "PropGasConstant,SpecificGasConstant" );
        tMMCounter++;

        //------------------------------------------------------------------------------
        // fill the constitutive model part of the parameter list

        // create fluid constitutive model
        aParameterLists( tCMIndex ).push_back( prm::create_constitutive_model_parameter_list() );
        aParameterLists.set( "constitutive_name", "CMFluid" );
        aParameterLists.set( "phase_name", "PhaseFluid" );
        aParameterLists.set( "constitutive_type", fem::Constitutive_Type::FLUID_COMPRESSIBLE_NEWTONIAN );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "P;VX,VY;TEMP", "Pressure,Velocity,Temperature" ) );
        aParameterLists.set( "properties",
                "PropViscosity,DynamicViscosity;"
                "PropConductivity,ThermalConductivity" );
        aParameterLists.set( "material_model", "MMFluid,ThermodynamicMaterialModel" );

        //------------------------------------------------------------------------------
        // fill the stabilization parameter part of the parameter list

        // create NITSCHE SP
        aParameterLists( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists.set( "stabilization_name", "NitscheSP" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::COMPRESSIBLE_DIRICHLET_NITSCHE );
        aParameterLists.set( "function_parameters", sNitscheGammas );
        aParameterLists.set( "leader_properties",
                "PropViscosity,DynamicViscosity;"
                "PropConductivity,ThermalConductivity" );

        // create DUMMY SP for GLS (simply has value 1.0 everywhere)
        aParameterLists( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists.set( "stabilization_name", "DummySP" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "leader_properties", "PropDummy,Material" );

        //------------------------------------------------------------------------------
        // fill the IWG part of the parameter list

        if ( tHaveBulk )
        {
            // bulk IWG
            aParameterLists( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
            aParameterLists.set( "IWG_name", "IWGBulk" );
            aParameterLists.set( "leader_phase_name", "PhaseFluid" );
            aParameterLists.set( "IWG_type", fem::IWG_Type::COMPRESSIBLE_NS_BULK );
            aParameterLists.set( "dof_residual", "P;VX,VY;TEMP" );
            aParameterLists.set( "leader_properties",
                    "PropViscosity,DynamicViscosity;"
                    "PropConductivity,ThermalConductivity;"
                    "PropHeatLoad,BodyHeatLoad" );
            aParameterLists.set( "leader_material_model", "MMFluid,FluidMM" );
            aParameterLists.set( "leader_constitutive_models", "CMFluid,FluidCM" );
            if ( tHaveGLS )
            {
                aParameterLists.set( "stabilization_parameters", "DummySP,GLS" );
            }
            }

        // Boundary IWG for top and bottom
        // Boundary IWG inlet
        // Boundary IWG outlet
        // aParameterLists( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        // aParameterLists.set( "IWG_name",                   "IWGBoundaryOutlet" );
        // aParameterLists.set( "IWG_bulk_type",               fem::Element_Type::SIDESET );
        // aParameterLists.set( "leader_phase_name",          "PhaseFluid" );
        // aParameterLists.set( "side_ordinals",              "2" );
        // aParameterLists.set( "IWG_type",                    fem::IWG_Type::COMPRESSIBLE_NS_BOUNDARY );
        // aParameterLists.set( "dof_residual",               "P;VX,VY;TEMP" );
        // aParameterLists.set( "leader_properties",          "PropInitialPressure,Pressure;"
        //                                                                               "PropViscosity,DynamicViscosity;"
        //                                                                               "PropConductivity,ThermalConductivity" );
        // aParameterLists.set( "leader_material_model",      "MMFluid,FluidMM" );
        // aParameterLists.set( "leader_constitutive_models", "CMFluid,FluidCM" );
        // tIWGCounter++;

        // Nitsche IWG for top and bottom
        if ( tHaveTopBottomBCs )
        {
            aParameterLists( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
            aParameterLists.set( "IWG_name", "IWGNitscheSides" );
            aParameterLists.set( "IWG_bulk_type", fem::Element_Type::SIDESET );
            aParameterLists.set( "leader_phase_name", "PhaseFluid" );
            aParameterLists.set( "side_ordinals", "1,3" );
            aParameterLists.set( "IWG_type", fem::IWG_Type::COMPRESSIBLE_NS_DIRICHLET_UNSYMMETRIC_NITSCHE );
            aParameterLists.set( "dof_residual", "P;VX,VY;TEMP" );
            aParameterLists.set( "leader_properties",
                    "PropZeroU,PrescribedVelocity;"
                    "PropSelectY,SelectVelocity;"
                    "PropViscosity,DynamicViscosity;"
                    "PropConductivity,ThermalConductivity" );
            aParameterLists.set( "leader_material_model", "MMFluid,FluidMM" );
            aParameterLists.set( "leader_constitutive_models", "CMFluid,FluidCM" );
            aParameterLists.set( "stabilization_parameters", "NitscheSP,NitschePenaltyParameter" );
            }

        // Nitsche IWGs for Outlets
        if ( tHaveFixedEnds )
        {
            tPropertyString = tPropertyString + ";PropZeroU,PrescribedVelocity";
        }
        else
        {
            tPropertyString = tPropertyString + ";PropInitialPressure,PrescribedDof1";

            if ( tUseUpwindForPressureBC )
            {
                tPropertyString = tPropertyString + ";PropUpwind,PressureUpwind";
            }
        }

        // Nitsche IWGs for Outlets
        aParameterLists( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGNitscheOutlets" );
        aParameterLists.set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "side_ordinals", "1,2,3,4" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::COMPRESSIBLE_NS_DIRICHLET_SYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "P;VX,VY;TEMP" );
        aParameterLists.set( "leader_properties", tPropertyString );
        aParameterLists.set( "leader_material_model", "MMFluid,FluidMM" );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,FluidCM" );
        aParameterLists.set( "stabilization_parameters", "NitscheSP,NitschePenaltyParameter" );

        if ( tHaveTimeContinuity )
        {
            // Time continuity for Pressure
            aParameterLists( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
            aParameterLists.set( "IWG_name", "IWGTimeContinuityPressure" );
            aParameterLists.set( "IWG_bulk_type", fem::Element_Type::BULK );
            aParameterLists.set( "leader_phase_name", "PhaseFluid" );
            aParameterLists.set( "IWG_type", fem::IWG_Type::TIME_CONTINUITY_DOF );
            aParameterLists.set( "dof_residual", "P" );
            aParameterLists.set( "leader_dof_dependencies", "P;VX,VY;TEMP" );
            aParameterLists.set( "leader_properties",
                    "PropWeightCurrent,WeightCurrent;"
                    "PropWeightPrevious,WeightPrevious;"
                    "PropInitialPressure,InitialCondition" );
            aParameterLists.set( "time_continuity", true );

            // Time continuity for Velocity
            aParameterLists( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
            aParameterLists.set( "IWG_name", "IWGTimeContinuityVelocity" );
            aParameterLists.set( "IWG_bulk_type", fem::Element_Type::BULK );
            aParameterLists.set( "leader_phase_name", "PhaseFluid" );
            aParameterLists.set( "IWG_type", fem::IWG_Type::TIME_CONTINUITY_DOF );
            aParameterLists.set( "dof_residual", "VX,VY" );
            aParameterLists.set( "leader_dof_dependencies", "P;VX,VY;TEMP" );
            aParameterLists.set( "leader_properties",
                    "PropWeightCurrent,WeightCurrent;"
                    "PropWeightPrevious,WeightPrevious;"
                    "PropInitialVelocity,InitialCondition" );
            aParameterLists.set( "time_continuity", true );

            // Time continuity for Temperature
            aParameterLists( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
            aParameterLists.set( "IWG_name", "IWGTimeContinuityTemp" );
            aParameterLists.set( "IWG_bulk_type", fem::Element_Type::BULK );
            aParameterLists.set( "leader_phase_name", "PhaseFluid" );
            aParameterLists.set( "IWG_type", fem::IWG_Type::TIME_CONTINUITY_DOF );
            aParameterLists.set( "dof_residual", "TEMP" );
            aParameterLists.set( "leader_dof_dependencies", "P;VX,VY;TEMP" );
            aParameterLists.set( "leader_properties",
                    "PropWeightCurrent,WeightCurrent;"
                    "PropWeightPrevious,WeightPrevious;"
                    "PropInitialTemperature,InitialCondition" );
            aParameterLists.set( "time_continuity", true );
            }

        //------------------------------------------------------------------------------
        // fill the IQI part of the parameter list

        // pressure
        aParameterLists( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists.set( "IQI_name", "IQIBulkP" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "P" );

        // velocity VX
        aParameterLists( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists.set( "IQI_name", "IQIBulkVX" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "VX,VY" );
        aParameterLists.set( "vectorial_field_index", 0 );

        // velocity VY
        aParameterLists( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists.set( "IQI_name", "IQIBulkVY" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "VX,VY" );
        aParameterLists.set( "vectorial_field_index", 1 );

        // temperature
        aParameterLists( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists.set( "IQI_name", "IQIBulkTEMP" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "TEMP" );

        // local Mach number
        aParameterLists( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists.set( "IQI_name", "IQIMachNumber" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::PROPERTY );
        aParameterLists.set( "leader_properties", "PropMachNumber,Property" );

        // local Reynolds number
        aParameterLists( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists.set( "IQI_name", "IQIReynoldsNumber" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::PROPERTY );
        aParameterLists.set( "leader_properties", "PropReynoldsNumber,Property" );

        // heat load distribution
        aParameterLists( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists.set( "IQI_name", "IQIHeatLoad" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::PROPERTY );
        aParameterLists.set( "leader_properties", "PropHeatLoad,Property" );

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        aParameterLists( tFEMIndex ).push_back( prm::create_computation_parameter_list() );

        // aParameterLists.set( "finite_difference_scheme",            tFEMFdScheme  );
        // aParameterLists.set( "finite_difference_perturbation_size", tFEMFdEpsilon );
    }

    //------------------------------------------------------------------------------

    void
    SOLParameterList( Vector< Vector< ParameterList > >& aParameterLists )
    {

        aParameterLists( 0 ).push_back( add_parameter_list( sol::SolverType::AMESOS_IMPL ) );

        aParameterLists( 1 ).push_back( moris::prm::create_linear_solver_parameter_list() );

        aParameterLists( 2 ).push_back( moris::prm::create_nonlinear_algorithm_parameter_list() );
        aParameterLists.set( "NLA_rel_res_norm_drop", tNewtonTolerance );
        aParameterLists.set( "NLA_relaxation_parameter", tNewtonRelaxation );
        aParameterLists.set( "NLA_max_iter", tMaxNewtonSteps );

        aParameterLists( 3 ).push_back( moris::prm::create_nonlinear_solver_parameter_list() );
        aParameterLists.set( "NLA_DofTypes", "P;VX,VY;TEMP" );

        aParameterLists( 4 ).push_back( moris::prm::create_time_solver_algorithm_parameter_list() );
        // aParameterLists.set("TSA_Nonlinear_Solver", 2);
        aParameterLists.set( "TSA_Num_Time_Steps", tNumTimeSteps );
        aParameterLists.set( "TSA_Time_Frame", tTimeFrame );

        aParameterLists( 5 ).push_back( moris::prm::create_time_solver_parameter_list() );
        aParameterLists.set( "TSA_DofTypes", "P;VX,VY;TEMP" );
        aParameterLists.set( "TSA_Initialize_Sol_Vec", "P," + ios::stringify( tInitialPressure ) + ";VX,0.0;VY,0.0;TEMP," + ios::stringify( tInitialTemperature ) );
        // aParameterLists.set("TSA_Initialize_Sol_Vec" , "P,1.0;VX,1.0;VY,0.0;TEMP,1.0" );
        aParameterLists.set( "TSA_Output_Indices", "0" );
        aParameterLists.set( "TSA_Output_Criteria", "Output_Criterion" );
        aParameterLists.set( "TSA_time_level_per_type", "P,2;VX,2;VY,2;TEMP,2" );

        aParameterLists( 6 ).push_back( moris::prm::create_solver_warehouse_parameterlist() );
        if ( tWriteJacToMatlab )
        {
            aParameterLists.set( "SOL_save_operator_to_matlab", "Jacobian.dat" );
        }

        aParameterLists( 7 ).push_back(  sol::PreconditionerType::NONE );
    }

    //------------------------------------------------------------------------------

    void
    MSIParameterList( Vector< Vector< ParameterList > >& aParameterLists )
    {
        aParameterLists( 0 ).push_back( prm::create_msi_parameter_list() );
        aParameterLists.set( "order_adofs_by_host", tOrderAdofsByHost );
    }

    //------------------------------------------------------------------------------

    void
    VISParameterList( Vector< Vector< ParameterList > >& aParameterLists )
    {
        aParameterLists( 0 ).push_back( prm::create_vis_parameter_list() );
        aParameterLists.set( "File_Name", std::pair< std::string, std::string >( "./", "Heated_Bubble_2D.exo" ) );
        aParameterLists.set( "Mesh_Type", vis::VIS_Mesh_Type::STANDARD );
        aParameterLists.set( "Set_Names", sFluid );
        aParameterLists.set( "Field_Names", "P,VX,VY,TEMP,Ma,Re,Q" );
        aParameterLists.set( "Field_Type", "NODAL,NODAL,NODAL,NODAL,NODAL,NODAL,NODAL" );
        aParameterLists.set( "IQI_Names", "IQIBulkP,IQIBulkVX,IQIBulkVY,IQIBulkTEMP,IQIMachNumber,IQIReynoldsNumber,IQIHeatLoad" );
        aParameterLists.set( "Save_Frequency", 1 );
        aParameterLists.set( "Time_Offset", 10.0 );
    }

    //------------------------------------------------------------------------------

    void
    MORISGENERALParameterList( Vector< Vector< ParameterList > >& aParameterLists )
    {
    }

    //------------------------------------------------------------------------------
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
