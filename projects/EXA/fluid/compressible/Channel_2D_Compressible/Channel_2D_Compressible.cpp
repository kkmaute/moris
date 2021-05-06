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
    // bulk phase
    std::string sFluid = "HMR_dummy_c_p0,HMR_dummy_n_p0";

    //Material Parameters
    moris::real tViscosity    = 1.716e-5; /* Dynamic Viscosity mu () */
    moris::real tHeatCapacity = 0.718e+3; /* Heat Capacity Cv () */
    moris::real tGasConstant  = 287.06;   /* Specific Gas Constant R () */
    moris::real tConductivity = 24.35e-3; /* Thermal Conductivity kappa() */

    // Boundary Conditions
    moris::real tInletVelocity     =  10.0;  /* Inlet velocity 10 m/s  () */
    moris::real tInletTemperature  = 300.0;  /* Inlet temperature 300 K  () */
    moris::real tOutletPressure    = 1.0e5;  /* Outlet pressure 1 bar() */

    // Nitsche Penalty
    moris::real tNitscheGamma   = 100.0;  /* Penalty for Dirichlet BC */

    // Constant function for properties
    void Func_Const(
            moris::Matrix< moris::DDRMat >                 & aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
            moris::fem::Field_Interpolator_Manager         * aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

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

    void OPTParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );
        tParameterlist( 0 )( 0 ) = prm::create_opt_problem_parameter_list();
        tParameterlist( 0 )( 0 ).set( "is_optimization_problem", false);
    }

    void HMRParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_hmr_parameter_list();

        tParameterlist( 0 )( 0 ).set( "number_of_elements_per_dimension", "100,20");
        tParameterlist( 0 )( 0 ).set( "domain_dimensions",                "10.0,2.0");
        tParameterlist( 0 )( 0 ).set( "domain_offset",                    "0.0,-1.5");
        tParameterlist( 0 )( 0 ).set( "domain_sidesets",                  "1,2,3,4");
        tParameterlist( 0 )( 0 ).set( "lagrange_output_meshes",           "0");

        tParameterlist( 0 )( 0 ).set( "lagrange_orders",  "2" );
        tParameterlist( 0 )( 0 ).set( "lagrange_pattern", "0" );
        tParameterlist( 0 )( 0 ).set( "bspline_orders",   "2" );
        tParameterlist( 0 )( 0 ).set( "bspline_pattern",  "0" );

        tParameterlist( 0 )( 0 ).set( "truncate_bsplines",  1 );
        tParameterlist( 0 )( 0 ).set( "use_number_aura", 1);
        tParameterlist( 0 )( 0 ).set( "use_multigrid",  0 );
        tParameterlist( 0 )( 0 ).set( "severity_level", 0 );
    }

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
        tParameterlist( 0 )( 0 ).set( "print_enriched_ig_mesh",    false );
        tParameterlist( 0 )( 0 ).set( "exodus_output_XTK_ig_mesh", false );
    }

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
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters",      ios::stringify(tViscosity) );
        tParameterList( tPropIndex )( tPropCounter ).set( "value_function",           "Func_Const") ;
        tPropCounter++;

        // Heat Capacity Cv
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name",            "PropHeatCapacity") ;
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters",      ios::stringify(tHeatCapacity) );
        tParameterList( tPropIndex )( tPropCounter ).set( "value_function",           "Func_Const") ;
        tPropCounter++;

        // Specific Gas Constant R
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name",            "PropGasConstant") ;
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters",      ios::stringify(tGasConstant) );
        tParameterList( tPropIndex )( tPropCounter ).set( "value_function",           "Func_Const") ;
        tPropCounter++;

        // Thermal Conductivity kappa
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name",            "PropConductivity") ;
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters",      ios::stringify(tConductivity) );
        tParameterList( tPropIndex )( tPropCounter ).set( "value_function",           "Func_Const") ;
        tPropCounter++;

        // velocity for no-slip
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name",            "PropZeroU") ;
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters",      "0.0;0.0") ;
        tParameterList( tPropIndex )( tPropCounter ).set( "value_function",           "Func_Const") ;
        tPropCounter++;

        // Inlet Velocity
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name",            "PropInletU") ;
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters",      ios::stringify(tInletVelocity) + ",0.0" );
        tParameterList( tPropIndex )( tPropCounter ).set( "value_function",           "Func_Const") ;
        tPropCounter++;

        // Outlet pressure
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name",            "PropOutletPressure") ;
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters",      std::to_string(tOutletPressure) );
        tParameterList( tPropIndex )( tPropCounter ).set( "value_function",           "Func_Const") ;
        tPropCounter++;

        // inlet temperature
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name",            "PropInletTemperature") ;
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters",      std::to_string(tInletTemperature) );
        tParameterList( tPropIndex )( tPropCounter ).set( "value_function",           "Func_Const") ;
        tPropCounter++;

        // create upwind weight factor
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name",            "PropUpwind") ;
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters",      "1.0") ;
        tParameterList( tPropIndex )( tPropCounter ).set( "value_function",           "Func_Const") ;
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
        // tParameterList( tSPIndex )( tSPCounter ).set( "master_dof_dependencies", std::pair< std::string, std::string >(
        //                                                                          "P;VX,VY;T", "Pressure,Velocity,Temperature" ) );
        tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters",     ios::stringify( tNitscheGamma ) );
        tParameterList( tSPIndex )( tSPCounter ).set( "master_properties",  "PropViscosity,DynamicViscosity;"
                                                                            "PropConductivity,ThermalConductivity" );
        tSPCounter++;

        // create DUMMY SP for GLS
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

        // bulk IWG
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name",                   "IWGBulk" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "master_phase_name",          "PhaseFluid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type",                   (uint) fem::IWG_Type::COMPRESSIBLE_NS_BULK );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual",               "P;VX,VY;T" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "master_material_model",      "MMFluid,FluidMM" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "master_constitutive_models", "CMFluid,FluidCM" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters",   "DummySP,GLS" );
        tIWGCounter++;

        // Boundary IWG for top and bottom
        // Boundary IWG inlet
        // Boundary IWG outlet
        // tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        // tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name",                   "IWGBoundaryOutlet" );
        // tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type",              (uint) fem::Element_Type::SIDESET );
        // tParameterList( tIWGIndex )( tIWGCounter ).set( "master_phase_name",          "PhaseFluid" );
        // tParameterList( tIWGIndex )( tIWGCounter ).set( "side_ordinals",              "2" );
        // tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type",                   (uint) fem::IWG_Type::COMPRESSIBLE_NS_BOUNDARY );
        // tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual",               "P;VX,VY;T" );
        // tParameterList( tIWGIndex )( tIWGCounter ).set( "master_properties",          "PropOutletPressure,Pressure;"
        //                                                                               "PropViscosity,DynamicViscosity;" 
        //                                                                               "PropConductivity,ThermalConductivity" );
        // tParameterList( tIWGIndex )( tIWGCounter ).set( "master_material_model",     "MMFluid,FluidMM" );
        // tParameterList( tIWGIndex )( tIWGCounter ).set( "master_constitutive_models", "CMFluid,FluidCM" );
        // tIWGCounter++;

        // Nitsche IWG for top and bottom
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name",                   "IWGNitscheSides" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type",              (uint) fem::Element_Type::SIDESET );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "master_phase_name",          "PhaseFluid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "side_ordinals",              "1,3" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type",                   (uint) fem::IWG_Type::COMPRESSIBLE_NS_DIRICHLET_SYMMETRIC_NITSCHE );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual",               "P;VX,VY;T" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "master_properties",          "PropZeroU,PrescribedVelocity;"
                                                                                      "PropUpwind,Upwind;"
                                                                                      "PropViscosity,DynamicViscosity;"
                                                                                      "PropConductivity,ThermalConductivity" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "master_material_model",     "MMFluid,FluidMM" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "master_constitutive_models", "CMFluid,FluidCM" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters",   "NitscheSP,NitschePenaltyParameter" );
        tIWGCounter++;

        // Nitsche IWG inlet
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name",                   "IWGNitscheInlet" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type",              (uint) fem::Element_Type::SIDESET );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "master_phase_name",          "PhaseFluid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "side_ordinals",              "4" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type",                   (uint) fem::IWG_Type::COMPRESSIBLE_NS_DIRICHLET_SYMMETRIC_NITSCHE );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual",               "P;VX,VY;T" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "master_properties",          "PropInletU,PrescribedVelocity;"
                                                                                      "PropInletTemperature,PrescribedDof3;"
                                                                                      "PropUpwind,Upwind;"
                                                                                      "PropViscosity,DynamicViscosity;"
                                                                                      "PropConductivity,ThermalConductivity" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "master_material_model",     "MMFluid,FluidMM" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "master_constitutive_models", "CMFluid,FluidCM" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters",   "NitscheSP,NitschePenaltyParameter" );
        tIWGCounter++;

        // Nitsche IWG outlet
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name",                   "IWGNitscheOutlet" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type",              (uint) fem::Element_Type::SIDESET );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "master_phase_name",          "PhaseFluid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "side_ordinals",              "2" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type",                   (uint) fem::IWG_Type::COMPRESSIBLE_NS_DIRICHLET_SYMMETRIC_NITSCHE );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual",               "P;VX,VY;T" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "master_properties",          "PropOutletPressure,PrescribedDof1;"
                                                                                      "PropUpwind,Upwind;"
                                                                                      "PropViscosity,DynamicViscosity;"
                                                                                      "PropConductivity,ThermalConductivity" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "master_material_model",     "MMFluid,FluidMM" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "master_constitutive_models", "CMFluid,FluidCM" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters",   "NitscheSP,NitschePenaltyParameter" );
        tIWGCounter++;


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

        // pressure
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name",               "IQIBulkTEMP" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "master_phase_name",      "PhaseFluid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type",               (uint) fem::IQI_Type::DOF );
        tParameterList( tIQIIndex )( tIQICounter ).set( "dof_quantity",           "TEMP" );
        tIQICounter++;

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        tParameterList( tFEMIndex ).resize( 1 );
        tParameterList( tFEMIndex )( 0 ) = prm::create_computation_parameter_list();

        //tParameterList( tFEMIndex )( 0 ).set( "finite_difference_scheme",            tFEMFdScheme  );
        //tParameterList( tFEMIndex )( 0 ).set( "finite_difference_perturbation_size", tFEMFdEpsilon );
    }

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
        tParameterlist( 2 )( 0 ).set("NLA_rel_res_norm_drop",    1e-06 );
        tParameterlist( 2 )( 0 ).set("NLA_relaxation_parameter", 0.9 );
        tParameterlist( 2 )( 0 ).set("NLA_max_iter",             50 );

        tParameterlist( 3 )( 0 ) = moris::prm::create_nonlinear_solver_parameter_list();
        tParameterlist( 3 )( 0 ).set("NLA_DofTypes", "P;VX,VY;TEMP") ;

        tParameterlist( 4 )( 0 ) = moris::prm::create_time_solver_algorithm_parameter_list();

        tParameterlist( 5 )( 0 ) = moris::prm::create_time_solver_parameter_list();
        tParameterlist( 5 )( 0 ).set("TSA_DofTypes",            "P;VX,VY;TEMP") ;
        tParameterlist( 5 )( 0 ).set("TSA_Initialize_Sol_Vec",  "P,1e5;VX,10.0;VY,0.0;TEMP,300.0") ;
        tParameterlist( 5 )( 0 ).set("TSA_Output_Indices",      "0") ;
        tParameterlist( 5 )( 0 ).set("TSA_Output_Crteria",      "Output_Criterion") ;

        tParameterlist( 6 )( 0 ) = moris::prm::create_solver_warehouse_parameterlist();
    }

    void MSIParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_msi_parameter_list();
        tParameterlist( 0 )( 0 ).set("P"    , 1);
        tParameterlist( 0 )( 0 ).set("VX"   , 1);
        tParameterlist( 0 )( 0 ).set("VY"   , 1);
        tParameterlist( 0 )( 0 ).set("TEMP" , 1);
    }

    void VISParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_vis_parameter_list();
        tParameterlist( 0 )( 0 ).set( "File_Name"     , std::pair< std::string, std::string >( "./", "Channel_2D_Compressible.exo" ) );
        tParameterlist( 0 )( 0 ).set( "Mesh_Type"     , static_cast< uint >( vis::VIS_Mesh_Type::STANDARD ) );
        tParameterlist( 0 )( 0 ).set( "Set_Names"     , sFluid );
        tParameterlist( 0 )( 0 ).set( "Field_Names"   , "P,VX,VY,TEMP" );
        tParameterlist( 0 )( 0 ).set( "Field_Type"    , "NODAL,NODAL,NODAL,NODAL" );
        tParameterlist( 0 )( 0 ).set( "IQI_Names"     , "IQIBulkP,IQIBulkVX,IQIBulkVY,IQIBulkTEMP" ) ;
        tParameterlist( 0 )( 0 ).set( "Save_Frequency", 1 );
        tParameterlist( 0 )( 0 ).set( "Time_Offset"   , 10.0 );
    }

    //------------------------------------------------------------------------------
}

//------------------------------------------------------------------------------
#ifdef  __cplusplus
}
#endif
