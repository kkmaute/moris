/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Channel_2D_Compressible.cpp
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

#ifdef  __cplusplus
extern "C"
{
#endif
//------------------------------------------------------------------------------
namespace moris
{
    // problem size
    moris::real tChannelLength = 1.0; // in meters
    moris::real tChannelHeight = 0.1;  // in meters
    bool tAllowSlip = false;

    // mesh
    moris::uint tIpOrder = 2;    // polynomial order for interpolation
    moris::uint tNumXElems = 1; // number of elements in x-direction
    moris::uint tNumYElems = 1; // number of elements in y-direction

    moris::real tElenX = tChannelLength/tNumXElems;
    moris::real tElenY = tChannelHeight/tNumYElems;

    // pseudo-transient - static problem run in transient
    bool tIsPseudoTransient = false;
    moris::real tPseudoTimeFrame = 1.0;

    // stabilization
    bool tHaveGLS = true;

    // shape functions used
    bool tUseLagrange = true;

    // write jacobian to file - for debugging
    bool tWriteJacAndResToMatlab = true;

    // write solution vectors to file in TSA - for debugging
    bool tWriteSolVecToHDF5 = true;

    // write LHS to file in DLA - for debugging
    bool tWriteLhsToHDF5 = true;

    // turn contributions on/off - for debugging
    bool tHaveBulk = false;
    bool tHaveTopBottomVelNitsche = true;
    bool tHavePressureInletBC = false;
    bool tHavePressureOutletBC = false;
    bool tHaveTemperatureInletBC = false;

    //------------------------------------------------------------------------------

    // Material Parameters
    moris::real tViscosity    = 1.716e-3; /* Dynamic Viscosity mu () */
    moris::real tHeatCapacity = 0.718e+3; /* Heat Capacity Cv () */
    moris::real tGasConstant  = 287.06;   /* Specific Gas Constant R () */
    moris::real tConductivity = 24.35e-3; /* Thermal Conductivity kappa() */

    // Boundary Conditions
    moris::real tInletVelocity     = 1.0;  /* Inlet velocity 10 m/s  () */
    moris::real tInletTemperature  = 300.0;  /* Inlet temperature 300 K  () */
    moris::real tOutletPressure    = 1.0e5;  /* Outlet pressure 1 bar() */

    moris::real tPrsFac = 1.0e-4;
    moris::real tInletPressure     = (1.0 + tPrsFac) * tOutletPressure;  /* Inlet pressure 1+eps bar() */

    moris::real tInitialDensity  = tOutletPressure/tGasConstant/tInletTemperature;
    moris::real tInitialXVelocity = std::sqrt(tPrsFac*tOutletPressure/tInitialDensity);
    moris::real tInitialYVelocity = 0.0;

    // Nitsche Penalty
    moris::real tNitscheGamma = 1.0e2 * std::sqrt( tElenX*tElenX + tElenY*tElenY );  /* Penalty for Dirichlet BC */

    // scale Nitsche penalty factors with material
    bool tAutoScaleNitschePenalty = true;

    //------------------------------------------------------------------------------

    // convert Nitsche penalty to string
    std::string sNitscheGammas =
            ios::stringify( 0.0 ) + ";" +
            ios::stringify( tNitscheGamma ) + ";" +
            ios::stringify( tNitscheGamma ) + ";" +
            ios::stringify( tNitscheGamma );

    // bulk phase
    std::string sFluid = "HMR_dummy_c_p0,HMR_dummy_n_p0";

    // Constant function for properties
    void Func_Const(
            moris::Matrix< moris::DDRMat >                 & aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >  & aParameters,
            moris::fem::Field_Interpolator_Manager         * aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    //------------------------------------------------------------------------------

    // initial pressure
	void Func_Initial_Pressure(
            moris::Matrix< moris::DDRMat >                 & aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >  & aParameters,
            moris::fem::Field_Interpolator_Manager         * aFIManager )
    {
        aPropMatrix = { { tOutletPressure } };
    }

    // initial velocity
	void Func_Initial_Velocity(
            moris::Matrix< moris::DDRMat >                 & aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >  & aParameters,
            moris::fem::Field_Interpolator_Manager         * aFIManager )
    {
        aPropMatrix = { { tInitialXVelocity } , { tInitialYVelocity } };
    }

    // initial pressure
	void Func_Initial_Temperature(
            moris::Matrix< moris::DDRMat >                 & aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >  & aParameters,
            moris::fem::Field_Interpolator_Manager         * aFIManager )
    {
        aPropMatrix = { { tInletTemperature } };
    }

    //------------------------------------------------------------------------------

    // local mach Number
	void Func_Mach_Number(
            moris::Matrix< moris::DDRMat >                 & aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >  & aParameters,
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
            Vector< moris::Matrix< moris::DDRMat > >  & aParameters,
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
            const Vector< real > & aGeometryParameters )
    {
        moris::real aReturnValue = aCoordinates( 1 ) - 10000; //tPlaneBottom - 0.01;
        return aReturnValue;
    }

    //------------------------------------------------------------------------------

    void OPTParameterList( Vector< Vector< ParameterList > > & aParameterLists )
    {
        aParameterLists( 0 ).push_back( prm::create_opt_problem_parameter_list() );
        aParameterLists( 0 ).set( "is_optimization_problem", false);
    }

    //------------------------------------------------------------------------------

    void HMRParameterList( Vector< Vector< ParameterList > > & aParameterLists )
    {
        aParameterLists( 0 ).push_back( prm::create_hmr_parameter_list() );

        aParameterLists( 0 ).set( "number_of_elements_per_dimension", ios::stringify( tNumXElems ) + "," + ios::stringify( tNumYElems ) );
        aParameterLists( 0 ).set( "domain_dimensions",                ios::stringify( tChannelLength ) + "," + ios::stringify( tChannelHeight ) );
        aParameterLists( 0 ).set( "domain_offset",                    "0.0," + ios::stringify( tChannelHeight / -2.0 ) );
        aParameterLists( 0 ).set( "domain_sidesets",                  "1,2,3,4");
        aParameterLists( 0 ).set( "lagrange_output_meshes",           "0");

        aParameterLists( 0 ).set( "lagrange_orders",  ios::stringify( tIpOrder ) );
        aParameterLists( 0 ).set( "lagrange_pattern", "0" );
        aParameterLists( 0 ).set( "bspline_orders",   ios::stringify( tIpOrder ) );
        aParameterLists( 0 ).set( "bspline_pattern",  "0" );
        if ( tUseLagrange )
        {
            aParameterLists( 0 ).set( "lagrange_to_bspline", "-1") ;
        }

        aParameterLists( 0 ).set( "truncate_bsplines", 1 );
        aParameterLists( 0 ).set( "use_number_aura",   1 );
        aParameterLists( 0 ).set( "use_multigrid",     0 );
        aParameterLists( 0 ).set( "severity_level",    0 );
    }

    //------------------------------------------------------------------------------

    void XTKParameterList( Vector< Vector< ParameterList > > & aParameterLists )
    {
        aParameterLists( 0 ).push_back( prm::create_xtk_parameter_list() );
        aParameterLists( 0 ).set( "decompose",                 true );
        aParameterLists( 0 ).set( "decomposition_type",        "conformal") ;
        aParameterLists( 0 ).set( "enrich",                    true );
        aParameterLists( 0 ).set( "basis_rank",                "bspline") ;
        aParameterLists( 0 ).set( "enrich_mesh_indices",       "0") ;
        aParameterLists( 0 ).set( "ghost_stab",                false );
        aParameterLists( 0 ).set( "multigrid",                 false );
        aParameterLists( 0 ).set( "print_enriched_ig_mesh",    true );
        aParameterLists( 0 ).set( "exodus_output_XTK_ig_mesh", true );
    }

    //------------------------------------------------------------------------------

    void GENParameterList( Vector< Vector< ParameterList > > & aParameterLists )
    {
        aParameterLists( 0 ).push_back( prm::create_gen_parameter_list() );

        // Dummy plane
        aParameterLists( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists( 1 ).set( "field_function_name", "Func_Dummy_Plane");
    }

    //------------------------------------------------------------------------------

    void FEMParameterList( Vector< Vector< ParameterList > > & aParameterLists )
    {
        // create a cell of cell of parameter list for fem

        // create a cell of cell of parameter list for fem
        uint tPropIndex  = 0;
        uint tCMIndex    = 1;
        uint tSPIndex    = 2;
        uint tIWGIndex   = 3;
        uint tIQIIndex   = 4;
        uint tFEMIndex   = 5;
        uint tPhaseIndex = 7;
        uint tMMIndex    = 8;

        //------------------------------------------------------------------------------

        aParameterLists( tPhaseIndex ).push_back( prm::create_phase_parameter_list() );
        aParameterLists( tPhaseIndex ).set( "phase_name",       "PhaseFluid" );
        aParameterLists( tPhaseIndex ).set( "phase_indices",    "0"  );

        //------------------------------------------------------------------------------
        // fill the property part of the parameter list

        // Dynamic Viscosity mu
        aParameterLists( tPropIndex ).push_back( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name",            "PropViscosity") ;
        aParameterLists( tPropIndex ).set( "function_parameters",      ios::stringify( tViscosity ) );
        aParameterLists( tPropIndex ).set( "value_function",           "Func_Const") ;

        // Heat Capacity Cv
        aParameterLists( tPropIndex ).push_back( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name",            "PropHeatCapacity") ;
        aParameterLists( tPropIndex ).set( "function_parameters",      ios::stringify( tHeatCapacity ) );
        aParameterLists( tPropIndex ).set( "value_function",           "Func_Const") ;

        // Specific Gas Constant R
        aParameterLists( tPropIndex ).push_back( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name",            "PropGasConstant") ;
        aParameterLists( tPropIndex ).set( "function_parameters",      ios::stringify( tGasConstant ) );
        aParameterLists( tPropIndex ).set( "value_function",           "Func_Const") ;

        // Thermal Conductivity kappa
        aParameterLists( tPropIndex ).push_back( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name",            "PropConductivity") ;
        aParameterLists( tPropIndex ).set( "function_parameters",      ios::stringify( tConductivity ) );
        aParameterLists( tPropIndex ).set( "value_function",           "Func_Const") ;

        // velocity for no-slip
        aParameterLists( tPropIndex ).push_back( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name",            "PropZeroU") ;
        aParameterLists( tPropIndex ).set( "function_parameters",      "0.0;0.0") ;
        aParameterLists( tPropIndex ).set( "value_function",           "Func_Const") ;

        // Outlet pressure
        aParameterLists( tPropIndex ).push_back( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name",            "PropOutletPressure") ;
        aParameterLists( tPropIndex ).set( "function_parameters",      ios::stringify( tOutletPressure ) );
        aParameterLists( tPropIndex ).set( "value_function",           "Func_Const") ;

        // Inlet pressure
        aParameterLists( tPropIndex ).push_back( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name",            "PropInletPressure") ;

        aParameterLists( tPropIndex ).set( "function_parameters",      ios::stringify( tInletPressure ) );
        aParameterLists( tPropIndex ).set( "value_function",           "Func_Const" );

        // inlet temperature
        aParameterLists( tPropIndex ).push_back( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name",            "PropInletTemperature" );
        aParameterLists( tPropIndex ).set( "function_parameters",      ios::stringify( tInletTemperature ) );
        aParameterLists( tPropIndex ).set( "value_function",           "Func_Const" );

        // Initial Pressure
        aParameterLists( tPropIndex ).push_back( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name",            "PropInitialPressure") ;
        aParameterLists( tPropIndex ).set( "function_parameters",      "1.0") ;
        aParameterLists( tPropIndex ).set( "value_function",           "Func_Initial_Pressure") ;

        // Initial Velocity
        aParameterLists( tPropIndex ).push_back( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name",            "PropInitialVelocity") ;
        aParameterLists( tPropIndex ).set( "function_parameters",      "1.0") ;
        aParameterLists( tPropIndex ).set( "value_function",           "Func_Initial_Velocity") ;

        // Initial Temperature
        aParameterLists( tPropIndex ).push_back( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name",            "PropInitialTemperature") ;
        aParameterLists( tPropIndex ).set( "function_parameters",      "1.0") ;
        aParameterLists( tPropIndex ).set( "value_function",           "Func_Initial_Temperature") ;

        // Mach Number
        aParameterLists( tPropIndex ).push_back( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name",            "PropMachNumber") ;
        aParameterLists( tPropIndex ).set( "function_parameters",      "1.0") ;
        aParameterLists( tPropIndex ).set( "value_function",           "Func_Mach_Number") ;

        // Reynolds Number
        aParameterLists( tPropIndex ).push_back( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name",            "PropReynoldsNumber") ;
        aParameterLists( tPropIndex ).set( "function_parameters",      "1.0") ;
        aParameterLists( tPropIndex ).set( "value_function",           "Func_Reynolds_Number") ;

        // create upwind weight factor
        aParameterLists( tPropIndex ).push_back( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name",            "PropUpwind") ;
        aParameterLists( tPropIndex ).set( "function_parameters",      "1.0") ;
        aParameterLists( tPropIndex ).set( "value_function",           "Func_Const") ;

        // select matrix for y-direction
        aParameterLists( tPropIndex ).push_back( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name",            "PropSelectY");
        aParameterLists( tPropIndex ).set( "function_parameters",      "0.0,0.0;0.0,1.0");

        // dummy property
        aParameterLists( tPropIndex ).push_back( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name",            "PropDummy") ;
        aParameterLists( tPropIndex ).set( "function_parameters",      "1.0") ;
        aParameterLists( tPropIndex ).set( "value_function",           "Func_Const") ;

        // zero property
        aParameterLists( tPropIndex ).push_back( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name",            "PropZero") ;
        aParameterLists( tPropIndex ).set( "function_parameters",      "0.0") ;
        aParameterLists( tPropIndex ).set( "value_function",           "Func_Const") ;

        // zero property
        aParameterLists( tPropIndex ).push_back( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name",            "PropVectorZero") ;
        aParameterLists( tPropIndex ).set( "function_parameters",      "0.0;0.0") ;
        aParameterLists( tPropIndex ).set( "value_function",           "Func_Const") ;

        //------------------------------------------------------------------------------
        // fill the material model part of the parameter list

        // init CM counter
        uint tMMCounter = 0;

        // create fluid constitutive model
        aParameterLists( tMMIndex ).push_back( prm::create_material_model_parameter_list() );
        aParameterLists( tMMIndex ).set( "material_name", "MMFluid" );
        aParameterLists( tMMIndex ).set( "phase_name",    "PhaseFluid") ;
        aParameterLists( tMMIndex ).set( "material_type",  fem::Material_Type::PERFECT_GAS );
        aParameterLists( tMMIndex ).set( "dof_dependencies",  std::pair< std::string, std::string >(
                                                                           "P;TEMP", "Pressure,Temperature" ) );
        aParameterLists( tMMIndex ).set( "properties", "PropHeatCapacity,IsochoricHeatCapacity;"
                                                                    "PropGasConstant,SpecificGasConstant"    );
        tMMCounter++;

        //------------------------------------------------------------------------------
        // fill the constitutive model part of the parameter list

        // create fluid constitutive model
        aParameterLists( tCMIndex ).push_back( prm::create_constitutive_model_parameter_list() );
        aParameterLists( tCMIndex ).set( "constitutive_name", "CMFluid" );
        aParameterLists( tCMIndex ).set( "phase_name",        "PhaseFluid") ;
        aParameterLists( tCMIndex ).set( "constitutive_type",  fem::Constitutive_Type::FLUID_COMPRESSIBLE_NEWTONIAN );
        aParameterLists( tCMIndex ).set( "dof_dependencies",  std::pair< std::string, std::string >(
                                                                           "P;VX,VY;TEMP", "Pressure,Velocity,Temperature" ) );
        aParameterLists( tCMIndex ).set( "properties", "PropViscosity,DynamicViscosity;"
                                                                    "PropConductivity,ThermalConductivity"    );
        aParameterLists( tCMIndex ).set( "material_model", "MMFluid,ThermodynamicMaterialModel" );

        //------------------------------------------------------------------------------
        // fill the stabilization parameter part of the parameter list

        // create NITSCHE SP
        aParameterLists( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( tSPIndex ).set( "stabilization_name",      "NitscheSP" );
        aParameterLists( tSPIndex ).set( "leader_phase_name",       "PhaseFluid" );
        aParameterLists( tSPIndex ).set( "stabilization_type",       fem::Stabilization_Type::COMPRESSIBLE_DIRICHLET_NITSCHE );
        aParameterLists( tSPIndex ).set( "function_parameters",     sNitscheGammas );
        if ( tAutoScaleNitschePenalty )
        {
            aParameterLists( tSPIndex ).set( "leader_properties",  "PropViscosity,DynamicViscosity;"
                                                                                "PropConductivity,ThermalConductivity" );
        }
        else
        {
            aParameterLists( tSPIndex ).set( "leader_properties",  "PropDummy,DynamicViscosity;"
                                                                                "PropDummy,ThermalConductivity" );
        }

        // create DUMMY SP for GLS
        aParameterLists( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( tSPIndex ).set( "stabilization_name",      "DummySP" );
        aParameterLists( tSPIndex ).set( "leader_phase_name",       "PhaseFluid" );
        aParameterLists( tSPIndex ).set( "stabilization_type",       fem::Stabilization_Type::DIRICHLET_NITSCHE );
        aParameterLists( tSPIndex ).set( "function_parameters",     "1.0" );
        aParameterLists( tSPIndex ).set( "leader_properties",       "PropDummy,Material" );

        //------------------------------------------------------------------------------
        // fill the IWG part of the parameter list

        // bulk IWG
        if( tHaveBulk )
        {
            aParameterLists( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
            aParameterLists( tIWGIndex ).set( "IWG_name",                   "IWGBulk" );
            aParameterLists( tIWGIndex ).set( "leader_phase_name",          "PhaseFluid" );
            aParameterLists( tIWGIndex ).set( "IWG_type",                    fem::IWG_Type::COMPRESSIBLE_NS_BULK );
            aParameterLists( tIWGIndex ).set( "dof_residual",               "P;VX,VY;TEMP" );
            aParameterLists( tIWGIndex ).set( "leader_properties",          "PropViscosity,DynamicViscosity;"
                                                                                          "PropConductivity,ThermalConductivity" );
            aParameterLists( tIWGIndex ).set( "leader_material_model",      "MMFluid,FluidMM" );
            aParameterLists( tIWGIndex ).set( "leader_constitutive_models", "CMFluid,FluidCM" );
            if ( tHaveGLS )
            {
                aParameterLists( tIWGIndex ).set( "stabilization_parameters",   "DummySP,GLS" );
            }
            }

        // Nitsche IWG for top and bottom
        if( tHaveTopBottomVelNitsche )
        {
            aParameterLists( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
            aParameterLists( tIWGIndex ).set( "IWG_name",                   "IWGNitscheSides" );
            aParameterLists( tIWGIndex ).set( "IWG_bulk_type",               fem::Element_Type::SIDESET );
            aParameterLists( tIWGIndex ).set( "leader_phase_name",          "PhaseFluid" );
            aParameterLists( tIWGIndex ).set( "side_ordinals",              "1,3" );
            aParameterLists( tIWGIndex ).set( "IWG_type",                    fem::IWG_Type::COMPRESSIBLE_NS_DIRICHLET_SYMMETRIC_NITSCHE );
            aParameterLists( tIWGIndex ).set( "dof_residual",               "P;VX,VY;TEMP" );
            if ( tAllowSlip )
            {
                aParameterLists( tIWGIndex ).set( "leader_properties",          "PropZeroU,PrescribedVelocity;"
                                                                                            "PropSelectY,SelectVelocity;"
                                                                                            "PropViscosity,DynamicViscosity;"
                                                                                            "PropConductivity,ThermalConductivity" );
            }
            else
            {
                aParameterLists( tIWGIndex ).set( "leader_properties",          "PropZeroU,PrescribedVelocity;"
                                                                                            "PropViscosity,DynamicViscosity;"
                                                                                            "PropConductivity,ThermalConductivity" );
            }
            aParameterLists( tIWGIndex ).set( "leader_material_model",      "MMFluid,FluidMM" );
            aParameterLists( tIWGIndex ).set( "leader_constitutive_models", "CMFluid,FluidCM" );
            aParameterLists( tIWGIndex ).set( "stabilization_parameters",   "NitscheSP,NitschePenaltyParameter" );
            }

        // Upwind pressure inlet IWG
        if( tHavePressureInletBC )
        {
            aParameterLists( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
            aParameterLists( tIWGIndex ).set( "IWG_name",                   "IWGNitscheInlet" );
            aParameterLists( tIWGIndex ).set( "IWG_bulk_type",               fem::Element_Type::SIDESET );
            aParameterLists( tIWGIndex ).set( "leader_phase_name",          "PhaseFluid" );
            aParameterLists( tIWGIndex ).set( "side_ordinals",              "4" );
            aParameterLists( tIWGIndex ).set( "IWG_type",                    fem::IWG_Type::COMPRESSIBLE_NS_DIRICHLET_UNSYMMETRIC_NITSCHE );
            aParameterLists( tIWGIndex ).set( "dof_residual",               "P;VX,VY;TEMP" );
            aParameterLists( tIWGIndex ).set( "leader_properties",          "PropInletPressure,PrescribedDof1;"
                                                                                          "PropViscosity,DynamicViscosity;"
                                                                                          "PropConductivity,ThermalConductivity" );
            aParameterLists( tIWGIndex ).set( "leader_material_model",      "MMFluid,FluidMM" );
            aParameterLists( tIWGIndex ).set( "leader_constitutive_models", "CMFluid,FluidCM" );
            aParameterLists( tIWGIndex ).set( "stabilization_parameters",   "NitscheSP,NitschePenaltyParameter" );
            }

        // Temperature Nitsche inlet IWG
        if( tHaveTemperatureInletBC )
        {
            aParameterLists( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
            aParameterLists( tIWGIndex ).set( "IWG_name",                   "IWGNitscheInlet" );
            aParameterLists( tIWGIndex ).set( "IWG_bulk_type",               fem::Element_Type::SIDESET );
            aParameterLists( tIWGIndex ).set( "leader_phase_name",          "PhaseFluid" );
            aParameterLists( tIWGIndex ).set( "side_ordinals",              "4" );
            aParameterLists( tIWGIndex ).set( "IWG_type",                    fem::IWG_Type::COMPRESSIBLE_NS_DIRICHLET_UNSYMMETRIC_NITSCHE );
            aParameterLists( tIWGIndex ).set( "dof_residual",               "P;VX,VY;TEMP" );
            aParameterLists( tIWGIndex ).set( "leader_properties",          "PropInletTemperature,PrescribedDof3;"
                                                                                          "PropViscosity,DynamicViscosity;"
                                                                                          "PropConductivity,ThermalConductivity" );
            aParameterLists( tIWGIndex ).set( "leader_material_model",      "MMFluid,FluidMM" );
            aParameterLists( tIWGIndex ).set( "leader_constitutive_models", "CMFluid,FluidCM" );
            aParameterLists( tIWGIndex ).set( "stabilization_parameters",   "NitscheSP,NitschePenaltyParameter" );
            }

        // Nitsche IWG outlet
        if( tHavePressureOutletBC )
        {
            aParameterLists( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
            aParameterLists( tIWGIndex ).set( "IWG_name",                   "IWGNitscheOutlet" );
            aParameterLists( tIWGIndex ).set( "IWG_bulk_type",               fem::Element_Type::SIDESET );
            aParameterLists( tIWGIndex ).set( "leader_phase_name",          "PhaseFluid" );
            aParameterLists( tIWGIndex ).set( "side_ordinals",              "2" );
            aParameterLists( tIWGIndex ).set( "IWG_type",                    fem::IWG_Type::COMPRESSIBLE_NS_DIRICHLET_SYMMETRIC_NITSCHE );
            aParameterLists( tIWGIndex ).set( "dof_residual",               "P;VX,VY;TEMP" );
            aParameterLists( tIWGIndex ).set( "leader_properties",          "PropOutletPressure,PrescribedDof1;"
                                                                                          "PropUpwind,PressureUpwind;"
                                                                                          "PropViscosity,DynamicViscosity;"
                                                                                          "PropConductivity,ThermalConductivity" );
            aParameterLists( tIWGIndex ).set( "leader_material_model",      "MMFluid,FluidMM" );
            aParameterLists( tIWGIndex ).set( "leader_constitutive_models", "CMFluid,FluidCM" );
            aParameterLists( tIWGIndex ).set( "stabilization_parameters",   "NitscheSP,NitschePenaltyParameter" );
            }

        //------------------------------------------------------------------------------
        // fill the IQI part of the parameter list

        // pressure
        aParameterLists( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name",               "IQIBulkP" );
        aParameterLists( tIQIIndex ).set( "leader_phase_name",      "PhaseFluid" );
        aParameterLists( tIQIIndex ).set( "IQI_type",                fem::IQI_Type::DOF );
        aParameterLists( tIQIIndex ).set( "dof_quantity",           "P" );

        // velocity VX
        aParameterLists( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name",               "IQIBulkVX" );
        aParameterLists( tIQIIndex ).set( "leader_phase_name",      "PhaseFluid" );
        aParameterLists( tIQIIndex ).set( "IQI_type",                fem::IQI_Type::DOF );
        aParameterLists( tIQIIndex ).set( "dof_quantity",           "VX,VY" );
        aParameterLists( tIQIIndex ).set( "vectorial_field_index",  0 );

        // velocity VY
        aParameterLists( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name",               "IQIBulkVY" );
        aParameterLists( tIQIIndex ).set( "leader_phase_name",      "PhaseFluid" );
        aParameterLists( tIQIIndex ).set( "IQI_type",                fem::IQI_Type::DOF );
        aParameterLists( tIQIIndex ).set( "dof_quantity",           "VX,VY" );
        aParameterLists( tIQIIndex ).set( "vectorial_field_index",  1 );

        // temperature
        aParameterLists( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name",               "IQIBulkTEMP" );
        aParameterLists( tIQIIndex ).set( "leader_phase_name",      "PhaseFluid" );
        aParameterLists( tIQIIndex ).set( "IQI_type",                fem::IQI_Type::DOF );
        aParameterLists( tIQIIndex ).set( "dof_quantity",           "TEMP" );

        // local Mach number
        aParameterLists( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name",               "IQIMachNumber" );
        aParameterLists( tIQIIndex ).set( "leader_phase_name",      "PhaseFluid" );
        aParameterLists( tIQIIndex ).set( "IQI_type",                fem::IQI_Type::PROPERTY );
        aParameterLists( tIQIIndex ).set( "leader_properties",      "PropMachNumber,Property");

        // local Reynolds number
        aParameterLists( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name",               "IQIReynoldsNumber" );
        aParameterLists( tIQIIndex ).set( "leader_phase_name",      "PhaseFluid" );
        aParameterLists( tIQIIndex ).set( "IQI_type",                fem::IQI_Type::PROPERTY );
        aParameterLists( tIQIIndex ).set( "leader_properties",      "PropReynoldsNumber,Property");

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        aParameterLists( tFEMIndex ).push_back( prm::create_computation_parameter_list() );

        //aParameterLists( tFEMIndex ).set( "finite_difference_scheme",            tFEMFdScheme  );
        //aParameterLists( tFEMIndex ).set( "finite_difference_perturbation_size", tFEMFdEpsilon );
    }

    //------------------------------------------------------------------------------

    void SOLParameterList( Vector< Vector< ParameterList > > & aParameterLists )
    {

        aParameterLists( 0 ).push_back( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL ) );

        aParameterLists( 1 ).push_back( moris::prm::create_linear_solver_parameter_list() );
        if ( tWriteLhsToHDF5 )
        {
            aParameterLists( 1 ).set("DLA_LHS_output_filename", "LHS" );
        }

        aParameterLists( 2 ).push_back( moris::prm::create_nonlinear_algorithm_parameter_list() );
        aParameterLists( 2 ).set("NLA_rel_res_norm_drop",    2.0e-06 );
        aParameterLists( 2 ).set("NLA_relaxation_parameter", 1.0 );
        aParameterLists( 2 ).set("NLA_max_iter",             10 );

        aParameterLists( 3 ).push_back( moris::prm::create_nonlinear_solver_parameter_list() );
        aParameterLists( 3 ).set("NLA_DofTypes", "P;VX,VY;TEMP") ;

        aParameterLists( 4 ).push_back( moris::prm::create_time_solver_algorithm_parameter_list() );

        aParameterLists( 5 ).push_back( moris::prm::create_time_solver_parameter_list() );
        aParameterLists( 5 ).set("TSA_DofTypes"           , "P;VX,VY;TEMP" );
        aParameterLists( 5 ).set("TSA_Output_Indices"     , "0" );
        aParameterLists( 5 ).set("TSA_Output_Criteria"     , "Output_Criterion" );

        // for pseudo-transient case, have a time step
        if ( tIsPseudoTransient )
        {
            aParameterLists( 4 ).set("TSA_Num_Time_Steps", 1 );
            aParameterLists( 4 ).set("TSA_Time_Frame",     tPseudoTimeFrame );
        }

        aParameterLists( 5 ).set("TSA_Initialize_Sol_Vec",  "P," + ios::stringify( tOutletPressure ) +
                                                                ";VX," + ios::stringify( tInitialXVelocity ) +
                                                                ";VY," + ios::stringify( tInitialYVelocity ) +
                                                                ";TEMP," + ios::stringify( tInletTemperature ) );

        aParameterLists( 6 ).push_back( moris::prm::create_solver_warehouse_parameterlist() );
        if ( tWriteJacAndResToMatlab )
        {
            aParameterLists( 6 ).set( "SOL_save_operator_to_matlab", "Channel_Compressible" );
        }
        if ( tWriteSolVecToHDF5 )
        {
            aParameterLists( 6 ).set( "TSA_Save_Sol_Vecs_to_file", "SolVec" );
        }

        aParameterLists( 7 ).push_back( moris::prm::create_preconditioner_parameter_list(sol::PreconditionerType::NONE) );
    }

    //------------------------------------------------------------------------------

    void MSIParameterList( Vector< Vector< ParameterList > > & aParameterLists )
    {
        aParameterLists( 0 ).push_back( prm::create_msi_parameter_list() );
    }

    //------------------------------------------------------------------------------

    void VISParameterList( Vector< Vector< ParameterList > > & aParameterLists )
    {
        aParameterLists( 0 ).push_back( prm::create_vis_parameter_list() );
        aParameterLists( 0 ).set( "File_Name"     , std::pair< std::string, std::string >( "./", "Channel_2D_Compressible.exo" ) );
        aParameterLists( 0 ).set( "Mesh_Type"     ,  vis::VIS_Mesh_Type::STANDARD ) ;
        aParameterLists( 0 ).set( "Set_Names"     , sFluid );
        aParameterLists( 0 ).set( "Field_Names"   , "P,VX,VY,TEMP,InletVelX,InletVelY,Ma,Re" );
        aParameterLists( 0 ).set( "Field_Type"    , "NODAL,NODAL,NODAL,NODAL,NODAL,NODAL,NODAL,NODAL" );
        aParameterLists( 0 ).set( "IQI_Names"     , "IQIBulkP,IQIBulkVX,IQIBulkVY,IQIBulkTEMP,IQIinletVelX,IQIinletVelY,IQIMachNumber,IQIReynoldsNumber" ) ;
        aParameterLists( 0 ).set( "Save_Frequency", 1 );
        aParameterLists( 0 ).set( "Time_Offset"   , 10.0 );
    }

    //------------------------------------------------------------------------------

    void MORISGENERALParameterList( Vector< Vector< ParameterList > > & aParameterLists )
    {

    }

    //------------------------------------------------------------------------------
}

//------------------------------------------------------------------------------
#ifdef  __cplusplus
}
#endif

