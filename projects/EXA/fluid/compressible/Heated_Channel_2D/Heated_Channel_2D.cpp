/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Heated_Channel_2D.cpp
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
    moris::real tChannelLength = 100.0;    // in meters
    moris::real tChannelHeight = 1.0;      // in meters

    // mesh
    moris::uint tIpOrder   = 2;      // polynomial order for interpolation
    moris::uint tNumXElems = 100;    // number of elements in x-direction
    moris::uint tNumYElems = 1;      // number of elements in y-direction

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
    moris::real tNitscheGammaP  = 0.0;   /* Penalty for pressure Dirichlet BCs */
    moris::real tNitscheGammaVX = 100.0; /* Penalty for x-vel Dirichlet BCs */
    moris::real tNitscheGammaVY = 100.0; /* Penalty for y-vel Dirichlet BCs */
    moris::real tNitscheGammaT  = 0.0;   /* Penalty for temperature Dirichlet BCs */

    // Time Continuity Weights
    moris::real tTimePenalty = 0.1;

    // related parameters
    moris::real tGamma = 1.0 + tGasConstant / tHeatCapacity;
    moris::real tCs    = std::sqrt( tGamma * tGasConstant * tInitialTemperature );

    //------------------------------------------------------------------------------

    // transient configuration
    int         tNumTimeSteps = 10;    // number of elements in time dimension
    moris::real tTimeStepSize = 20.0 * tChannelLength / 100.0 / tCs / 2.0;
    moris::real tTimeFrame    = (real)tNumTimeSteps * tTimeStepSize;    // resulting duration
    moris::real tTCWeight     = tTimePenalty / tTimeStepSize;

    // Newton configuration
    moris::real tNewtonRelaxation = 1.0;
    moris::real tNewtonTolerance  = 1.0e-9;
    int         tMaxNewtonSteps   = 10;

    // stabilization
    bool tHaveGLS = false;

    // BC configuration
    bool tUseUpwindForPressureBC = true;     // fix pressure at ends of channel
    bool tHaveFixedEnds          = false;    // close off ends of channel, impose zero velocity
    bool tHaveTopBottomBCs       = true;     // impose zero velocity in normal direction on channel sides

    // switch for time continuity - for debugging
    bool tHaveTimeContinuity = true;

    // switch for bulk contribution - for debugging
    bool tHaveBulk = true;

    // write jacobian to file - for debugging
    bool tWriteJacToMatlab = false;

    // write solution vectors to file in TSA - for debugging
    bool tWriteSolVecToHDF5 = false;

    // write LHS to file in DLA - for debugging
    bool tWriteLhsToHDF5 = false;

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
        // get x-coordinate
        real tX = aFIManager->get_IP_geometry_interpolator()->valx()( 0 );

        // compute heat load function
        real tQ = std::exp( -1.0 * tQfac * std::pow( 2.0 * tX / tChannelLength - 1.0, 2.0 ) );

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
    OPTParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "is_optimization_problem", false );
    }

    //------------------------------------------------------------------------------

    void
    HMRParameterList( Module_Parameter_Lists& aParameterLists )
    {
        // print values
        std::cout << "Time step size: " << tTimeStepSize << " \n"
                  << std::flush;
        std::cout << "Time continuity weight: " << tTCWeight << " \n"
                  << std::flush;

        aParameterLists.set( "number_of_elements_per_dimension", ios::stringify( tNumXElems ) + "," + ios::stringify( tNumYElems ) );
        aParameterLists.set( "domain_dimensions", ios::stringify( tChannelLength ) + "," + ios::stringify( tChannelHeight ) );
        aParameterLists.set( "domain_offset", "0.0,0.0" );
        aParameterLists.set( "domain_sidesets", "1,2,3,4" );
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
    XTKParameterList( Module_Parameter_Lists& aParameterLists )
    {
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
    GENParameterList( Module_Parameter_Lists& aParameterLists )
    {
        // Dummy plane
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( gen::Field_Type::USER_DEFINED );
        aParameterLists.set( "field_function_name", "Func_Dummy_Plane" );
    }

    //------------------------------------------------------------------------------

    void
    FEMParameterList( Module_Parameter_Lists& aParameterLists )
    {
        //------------------------------------------------------------------------------

        aParameterLists( FEM::PHASES ).add_parameter_list();
        aParameterLists.set( "phase_name", "PhaseFluid" );
        aParameterLists.set( "phase_indices", "0" );

        //------------------------------------------------------------------------------
        // fill the property part of the parameter list

        // Dynamic Viscosity mu
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropViscosity" );
        aParameterLists.set( "function_parameters", ios::stringify( tViscosity ) );
        aParameterLists.set( "value_function", "Func_Const" );

        // Heat Capacity Cv
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropHeatCapacity" );
        aParameterLists.set( "function_parameters", ios::stringify( tHeatCapacity ) );
        aParameterLists.set( "value_function", "Func_Const" );

        // Specific Gas Constant R
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropGasConstant" );
        aParameterLists.set( "function_parameters", ios::stringify( tGasConstant ) );
        aParameterLists.set( "value_function", "Func_Const" );

        // Thermal Conductivity kappa
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropConductivity" );
        aParameterLists.set( "function_parameters", ios::stringify( tConductivity ) );
        aParameterLists.set( "value_function", "Func_Const" );

        // Volumetric Heat load
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropHeatLoad" );
        aParameterLists.set( "function_parameters", ios::stringify( tPeakHeatLoad ) );
        aParameterLists.set( "value_function", "Func_Heat_Load_Distribution" );    // Func_Heat_Load_Distribution

        // velocity for no-slip
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropZeroU" );
        aParameterLists.set( "function_parameters", "0.0;0.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // Outlet pressure BCs
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropPressureBC" );
        aParameterLists.set( "function_parameters", std::to_string( tInitialPressure ) );
        aParameterLists.set( "value_function", "Func_Const" );

        // Initial Pressure
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropInitialPressure" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Initial_Pressure" );

        // Initial Velocity
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropInitialVelocity" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Initial_Velocity" );

        // Initial Temperature
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropInitialTemperature" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Initial_Temperature" );

        // Mach Number
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropMachNumber" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Mach_Number" );

        // Reynolds Number
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropReynoldsNumber" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Reynolds_Number" );

        // create upwind weight factor
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropUpwind" );
        aParameterLists.set( "function_parameters", "1.0" );    // 3.774074960608205e-03
        aParameterLists.set( "value_function", "Func_Const" );

        // select matrix for y-direction
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropSelectY" );
        aParameterLists.set( "function_parameters", "0.0,0.0;0.0,1.0" );

        // time continuity weights
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropWeightCurrent" );
        aParameterLists.set( "function_parameters", std::to_string( tTCWeight ) );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropWeightPrevious" );
        aParameterLists.set( "function_parameters", std::to_string( tTCWeight ) );
        aParameterLists.set( "value_function", "Func_Const" );

        // dummy property
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDummy" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // zero property
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropZero" );
        aParameterLists.set( "function_parameters", "0.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // zero property
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropVectorZero" );
        aParameterLists.set( "function_parameters", "0.0;0.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        //------------------------------------------------------------------------------
        // fill the material model part of the parameter list

        // init CM counter
        uint tMMCounter = 0;

        // create fluid constitutive model
        aParameterLists( FEM::MATERIAL_MODELS ).add_parameter_list();
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
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
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
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "NitscheSP" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::COMPRESSIBLE_DIRICHLET_NITSCHE );
        aParameterLists.set( "function_parameters", sNitscheGammas );
        aParameterLists.set( "leader_properties",
                "PropViscosity,DynamicViscosity;"
                "PropConductivity,ThermalConductivity" );

        // create DUMMY SP for GLS (simply has value 1.0 everywhere)
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
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
            aParameterLists( FEM::IWG ).add_parameter_list();
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
            aParameterLists( FEM::IWG ).add_parameter_list();
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
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGNitscheOutlets" );
        aParameterLists.set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "side_ordinals", "2,4" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::COMPRESSIBLE_NS_DIRICHLET_SYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "P;VX,VY;TEMP" );
        aParameterLists.set( "leader_properties", tPropertyString );
        aParameterLists.set( "leader_material_model", "MMFluid,FluidMM" );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,FluidCM" );
        aParameterLists.set( "stabilization_parameters", "NitscheSP,NitschePenaltyParameter" );

        if ( tHaveTimeContinuity )
        {
            // Time continuity for Pressure
            aParameterLists( FEM::IWG ).add_parameter_list();
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
            aParameterLists( FEM::IWG ).add_parameter_list();
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
            aParameterLists( FEM::IWG ).add_parameter_list();
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
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkP" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "P" );

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

        // temperature
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkTEMP" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "TEMP" );

        // local Mach number
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIMachNumber" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::PROPERTY );
        aParameterLists.set( "leader_properties", "PropMachNumber,Property" );

        // local Reynolds number
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIReynoldsNumber" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::PROPERTY );
        aParameterLists.set( "leader_properties", "PropReynoldsNumber,Property" );

        // heat load distribution
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIHeatLoad" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::PROPERTY );
        aParameterLists.set( "leader_properties", "PropHeatLoad,Property" );

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        aParameterLists( FEM::COMPUTATION );

        // aParameterLists.set( "finite_difference_scheme",            tFEMFdScheme  );
        // aParameterLists.set( "finite_difference_perturbation_size", tFEMFdEpsilon );
    }

    //------------------------------------------------------------------------------

    void
    SOLParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::AMESOS_IMPL );

        aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list();
        if ( tWriteLhsToHDF5 )
        {
            aParameterLists.set( "DLA_LHS_output_filename", "LHS" );
        }

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "NLA_rel_res_norm_drop", tNewtonTolerance );
        aParameterLists.set( "NLA_relaxation_parameter", tNewtonRelaxation );
        aParameterLists.set( "NLA_max_iter", tMaxNewtonSteps );

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "NLA_DofTypes", "P;VX,VY;TEMP" );

        aParameterLists( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list();
        // aParameterLists.set("TSA_Nonlinear_Solver", 2);
        aParameterLists.set( "TSA_Num_Time_Steps", tNumTimeSteps );
        aParameterLists.set( "TSA_Time_Frame", tTimeFrame );

        aParameterLists( SOL::TIME_SOLVERS ).add_parameter_list();
        aParameterLists.set( "TSA_DofTypes", "P;VX,VY;TEMP" );
        aParameterLists.set( "TSA_Initialize_Sol_Vec", "P," + ios::stringify( tInitialPressure ) + ";VX,0.0;VY,0.0;TEMP," + ios::stringify( tInitialTemperature ) );
        // aParameterLists.set("TSA_Initialize_Sol_Vec" , "P,1.0;VX,1.0;VY,0.0;TEMP,1.0" );
        // aParameterLists.set("TSA_Initialize_Sol_Vec" , "InitSolVec.hdf5" );
        aParameterLists.set( "TSA_Output_Indices", "0" );
        aParameterLists.set( "TSA_Output_Criteria", "Output_Criterion" );
        aParameterLists.set( "TSA_time_level_per_type", "P,2;VX,2;VY,2;TEMP,2" );

        aParameterLists( SOL::SOLVER_WAREHOUSE );
        if ( tWriteJacToMatlab )
        {
            aParameterLists.set( "SOL_save_operator_to_matlab", "Heated_Channel" );
        }
        if ( tWriteSolVecToHDF5 )
        {
            aParameterLists.set( "TSA_Save_Sol_Vecs_to_file", "SolVec" );
        }

        aParameterLists( SOL::PRECONDITIONERS ).add_parameter_list( sol::PreconditionerType::NONE );
    }

    //------------------------------------------------------------------------------

    void
    MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "order_adofs_by_host", tOrderAdofsByHost );
    }

    //------------------------------------------------------------------------------

    void
    VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "File_Name", std::pair< std::string, std::string >( "./", "Heated_Channel_2D.exo" ) );
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
    MORISGENERALParameterList( Module_Parameter_Lists& aParameterLists )
    {
    }

    //------------------------------------------------------------------------------
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
