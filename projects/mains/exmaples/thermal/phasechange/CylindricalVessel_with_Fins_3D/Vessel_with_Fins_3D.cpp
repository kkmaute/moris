#include <string>
#include <iostream>
#include <math.h>

#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_MSI_Equation_Object.hpp"
#include "cl_TSA_Time_Solver.hpp"
#include "cl_DLA_Solver_Interface.hpp"
#include "cl_DLA_Linear_Solver_Aztec.hpp"
#include "cl_PRM_FEM_Parameters.hpp"
#include "cl_PRM_MSI_Parameters.hpp"
#include "cl_PRM_SOL_Parameters.hpp"
#include "cl_PRM_VIS_Parameters.hpp"
#include "cl_PRM_HMR_Parameters.hpp"
#include "fn_PRM_GEN_Parameters.hpp"
#include "cl_PRM_XTK_Parameters.hpp"
#include "cl_PRM_OPT_Parameters.hpp"
#include "fn_equal_to.hpp"

#include "AztecOO.h"


#ifdef  __cplusplus
extern "C"
{
#endif
//------------------------------------------------------------------------------
namespace moris
{

    // Geometry Parameters
    moris::real tCenterX         = 0.5;   /* y bottom plane (m) */
    moris::real tCenterY         = 0.5;   /* y bottom plane (m) */
    moris::real tOuterRad        = 0.4;   /* y top plane    (m) */
    moris::real tInnerRad        = 0.385; /* y top plane    (m) */
    moris::real tFinRatlen       = 0.50;
    moris::real tFinWidth        = 0.01;
    moris::real tFinExp          = 24.0;

    // Constant function for properties
    void Func_Const(
            moris::Matrix< moris::DDRMat >                 & aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
            moris::fem::Field_Interpolator_Manager         * aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    void Func_Initial_Condition(
            moris::Matrix< moris::DDRMat >                 & aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
            moris::fem::Field_Interpolator_Manager         * aFIManager )
    {
        aPropMatrix = {{ 489.0 }};
    }

    bool Output_Criterion( moris::tsa::Time_Solver * aTimeSolver )
    {
        return true;
    }

    moris::real Outer_Ring(
            const moris::Matrix< DDRMat >     & aCoordinates,
            const moris::Cell< moris::real* > & aGeometryParameters )
    {
        moris::real aReturnValue = tOuterRad - std::sqrt( std::pow(aCoordinates(0) - tCenterX,2.0) + std::pow(aCoordinates(1) - tCenterY,2.0));
        return std::abs(aReturnValue)<1e-8 ? 1e-8 : aReturnValue;
    }

    moris::real Inner_Ring(
            const moris::Matrix< DDRMat >     & aCoordinates,
            const moris::Cell< moris::real* > & aGeometryParameters )
    {
        moris::real aReturnValue = tInnerRad - std::sqrt( std::pow(aCoordinates(0) - tCenterX,2.0) + std::pow(aCoordinates(1) - tCenterY,2.0));
        return std::abs(aReturnValue)<1e-8 ? 1e-8 : aReturnValue;
    }

    moris::real Fins(
            const moris::Matrix< DDRMat >     & aCoordinates,
            const moris::Cell< moris::real* > & aGeometryParameters )
    {
        moris::real tFinLen = tFinRatlen * tOuterRad;
        moris::real tFinRad = tOuterRad - tFinLen/2.0;

        moris::real tInc    = 2.0*M_PI/6.0;
        moris::real tAlpha  = 0.0;

        moris::real aReturnValue = -1e12;

        while ( tAlpha < 2*M_PI )
        {
            moris::Matrix< DDRMat >  tFinCenter = { {tFinRad * std::cos(tAlpha) + tCenterX },{tFinRad *std::sin(tAlpha) + tCenterY},{0}};

            moris::Matrix< DDRMat > tTmatTrans = {{std::cos(tAlpha),std::sin(tAlpha),0},{-std::sin(tAlpha),std::cos(tAlpha),0},{0,0,1}};

            moris::Matrix< DDRMat > tLocalCoordinates = tTmatTrans * trans(aCoordinates);
            moris::Matrix< DDRMat > tLocalFinCenter   = tTmatTrans * tFinCenter;

            moris::real tPhi =  1.0 - std::pow( std::pow( (tLocalCoordinates(0) - tLocalFinCenter(0))/tFinLen,   tFinExp)
            +std::pow( (tLocalCoordinates(1) - tLocalFinCenter(1))/tFinWidth, tFinExp), 1.0/tFinExp);

            aReturnValue = std::max(aReturnValue,tPhi);

            tAlpha += tInc;
        }

        return std::abs(aReturnValue)<1e-8 ? 1e-8 : aReturnValue;

    }

    moris::Matrix< DDRMat > Func_Sensitivity(
            const moris::Matrix< DDRMat >     & aCoordinates,
            const moris::Cell< moris::real* > & aGeometryParameters )
    {
        moris::Matrix< DDRMat > aReturnValue;
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

        tParameterlist( 0 )( 0 ).set( "number_of_elements_per_dimension", std::string("32,32,4"));
        tParameterlist( 0 )( 0 ).set( "domain_dimensions",                std::string("1.0,1.0,0.1"));
        tParameterlist( 0 )( 0 ).set( "domain_offset",                    std::string("0,0,0") );
        tParameterlist( 0 )( 0 ).set( "domain_sidesets",                  std::string("1,2,3,4,5,6"));
        tParameterlist( 0 )( 0 ).set( "lagrange_output_meshes",           std::string("0"));

        tParameterlist( 0 )( 0 ).set( "lagrange_orders",  std::string( "1" ));
        tParameterlist( 0 )( 0 ).set( "lagrange_pattern", std::string( "0" ));
        tParameterlist( 0 )( 0 ).set( "bspline_orders",   std::string( "1" ));
        tParameterlist( 0 )( 0 ).set( "bspline_pattern",  std::string( "0" ));

        tParameterlist( 0 )( 0 ).set( "lagrange_to_bspline", std::string("0") );

        tParameterlist( 0 )( 0 ).set( "truncate_bsplines",  1 );
        tParameterlist( 0 )( 0 ).set( "refinement_buffer",  2 );
        tParameterlist( 0 )( 0 ).set( "staircase_buffer",   2 );
        tParameterlist( 0 )( 0 ).set( "initial_refinement", 0 );
		
        tParameterlist( 0 )( 0 ).set( "use_number_aura", 1);
        
		tParameterlist( 0 )( 0 ).set( "use_multigrid",  0 );
        tParameterlist( 0 )( 0 ).set( "severity_level", 0 );

        tParameterlist( 0 )( 0 ).set( "adaptive_refinement_level", 4 );
    }

    void XTKParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_xtk_parameter_list();
        tParameterlist( 0 )( 0 ).set( "decompose",                 true );
        tParameterlist( 0 )( 0 ).set( "decomposition_type",        std::string("conformal") );
        tParameterlist( 0 )( 0 ).set( "enrich",                    true );
        tParameterlist( 0 )( 0 ).set( "basis_rank",                std::string("bspline") );
        tParameterlist( 0 )( 0 ).set( "enrich_mesh_indices",       std::string("0") );
        tParameterlist( 0 )( 0 ).set( "ghost_stab",                true );
        tParameterlist( 0 )( 0 ).set( "multigrid",                 false );
        tParameterlist( 0 )( 0 ).set( "verbose",                   true );
        tParameterlist( 0 )( 0 ).set( "print_enriched_ig_mesh",    true );
        tParameterlist( 0 )( 0 ).set( "exodus_output_XTK_ig_mesh", false );
        tParameterlist( 0 )( 0 ).set( "high_to_low_dbl_side_sets", true );
    }

    void GENParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 3 );
        tParameterlist( 0 ).resize( 1 );

        // Main GEN parameter list
        tParameterlist( 0 )( 0 ) = prm::create_gen_parameter_list();
        tParameterlist( 0 )( 0 ).set( "HMR_refinements", 2 );

        // init geometry counter
        uint tGeoCounter = 0;

        // Geometry parameter lists
        tParameterlist( 1 ).push_back( prm::create_user_defined_geometry_parameter_list() );
        tParameterlist( 1 )( tGeoCounter ).set( "field_function_name", "Outer_Ring");
        tParameterlist( 1 )( tGeoCounter ).set( "sensitivity_function_name", "Func_Sensitivity");
        tParameterlist( 1 )( tGeoCounter ).set( "constant_parameters", "");
        tGeoCounter++;

        tParameterlist( 1 ).push_back( prm::create_user_defined_geometry_parameter_list() );
        tParameterlist( 1 )( tGeoCounter ).set( "field_function_name", "Inner_Ring");
        tParameterlist( 1 )( tGeoCounter ).set( "sensitivity_function_name", "Func_Sensitivity");
        tParameterlist( 1 )( tGeoCounter ).set( "constant_parameters", "");
        tGeoCounter++;

        tParameterlist( 1 ).push_back( prm::create_user_defined_geometry_parameter_list() );
        tParameterlist( 1 )( tGeoCounter ).set( "field_function_name", "Fins");
        tParameterlist( 1 )( tGeoCounter ).set( "sensitivity_function_name", "Func_Sensitivity");
        tParameterlist( 1 )( tGeoCounter ).set( "constant_parameters", "");
        tGeoCounter++;
    }

    void FEMParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterList )
    {
        // create a cell of cell of parameter list for fem
        tParameterList.resize( 5 );

        //------------------------------------------------------------------------------
        // init property counter
        uint tPropCounter = 0;

        // create parameter list for outer ring

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ) = prm::create_property_parameter_list();
        tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropDensity_Outer") );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      std::string("1.0") );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Const") );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropHeatCapacity_Outer") );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      std::string("1.0") );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Const") );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropConductivity_Outer") );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      std::string("0.01") );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Const") );
        tPropCounter++;

        // create parameter list for inner circle

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ) = prm::create_property_parameter_list();
        tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropDensity_Inner") );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      std::string("1.0") );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Const") );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropHeatCapacity_Inner") );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      std::string("1.0") );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Const") );
        tPropCounter++;


        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropConductivity_Inner") );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      std::string("0.01") );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Const") );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropLatentHeat_Inner") );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      std::string("250.0") );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Const") );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropPCTemp_Inner") );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      std::string("500.") );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Const") );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropPhaseState_Inner") );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      std::string("2.0") );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Const") );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropPCconst_Inner") );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      std::string("20.0") );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Const") );
        tPropCounter++;

        // flux on outer ring

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropImposedFlux") );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      std::string("100.0") );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Const") );
        tPropCounter++;

        // algorithmic parameters

        // create parameter list for property 6
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropWeightCurrent") );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      std::string("100.0") );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Const") );
        tPropCounter++;

        // create parameter list for property 7
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropWeightPrevious") );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      std::string("100.0") );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Const") );
        tPropCounter++;

        // create parameter list for property 8
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropInitialCondition") );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Initial_Condition") );
        tPropCounter++;

        //------------------------------------------------------------------------------
        // init CM counter
        uint tCMCounter = 0;

        // create parameter list for constitutive model - outer ring
        tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_name", std::string("CMDiffusion_Outer") );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::DIFF_LIN_ISO ) );
        tParameterList( 1 )( tCMCounter ).set( "dof_dependencies",  std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( 1 )( tCMCounter ).set( "properties",        
	        std::string("PropConductivity_Outer , Conductivity;") +
	        std::string("PropDensity_Outer      , Density;")      +
		std::string("PropHeatCapacity_Outer     , Heat_Capacity") );
        tCMCounter++;
 
        // create parameter list for constitutive model - outer ring
        tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_name", std::string("CMDiffusion_Outer_Dummy") );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::DIFF_LIN_ISO ) );
        tParameterList( 1 )( tCMCounter ).set( "dof_dependencies",  std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( 1 )( tCMCounter ).set( "properties",        
	        std::string("PropConductivity_Outer , Conductivity;") +
	        std::string("PropDensity_Outer      , Density;")      +
		std::string("PropHeatCapacity_Outer     , Heat_Capacity") );
        tCMCounter++;

        // create parameter list for constitutive model - inner circle
        tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_name", std::string("CMDiffusion_Inner") );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::DIFF_LIN_ISO_PC ) );
        tParameterList( 1 )( tCMCounter ).set( "dof_dependencies",  std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( 1 )( tCMCounter ).set( "properties",
                std::string("PropConductivity_Inner, Conductivity;")         +
                std::string("PropDensity_Inner     , Density;")              +
                std::string("PropHeatCapacity_Inner, Heat_Capacity;")        +
                std::string("PropLatentHeat_Inner  , Latent_Heat;")          +
                std::string("PropPCTemp_Inner      , PC_Temp;")              +
                std::string("PropPhaseState_Inner  , Phase_State_Function;") +
                std::string("PropPCconst_Inner     , Phase_Change_Const")    );
        tCMCounter++;

        //------------------------------------------------------------------------------
        // init SP counter
        uint tSPCounter = 0;

        // create parameter list for stabilization parameter 2
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name",      std::string("SPGPTemp_Outer") );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type",      static_cast< uint >( fem::Stabilization_Type::GHOST_DISPL ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters",     std::string("0.01") );
        tParameterList( 2 )( tSPCounter ).set( "master_properties",       std::string("PropConductivity_Outer,Material") );
        tSPCounter++;

        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name",      std::string("SPGPTemp_Inner") );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type",      static_cast< uint >( fem::Stabilization_Type::GHOST_DISPL ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters",     std::string("0.01") );
        tParameterList( 2 )( tSPCounter ).set( "master_properties",       std::string("PropConductivity_Inner,Material") );
        tSPCounter++;

        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name",      std::string("SPGGLSDiffusion_Inner") );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type",      static_cast< uint >( fem::Stabilization_Type::GGLS_DIFFUSION_PC ) );
        tParameterList( 2 )( tSPCounter ).set( "master_dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( 2 )( tSPCounter ).set( "master_properties",
                std::string("PropConductivity_Inner , Conductivity;")         +
                std::string("PropDensity_Inner      , Density;")              +
                std::string("PropHeatCapacity_Inner , Heat_Capacity;")        +
                std::string("PropLatentHeat_Inner   , Latent_Heat;")          +
                std::string("PropPCTemp_Inner       , PC_Temp;")              +
                std::string("PropPhaseState_Inner   , Phase_State_Function;") +
                std::string("PropPCconst_Inner      , Phase_Change_Const")    );
        tSPCounter++;

        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name",  std::string("SPInterfaceNitscheOuter") );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type",  static_cast< uint >( fem::Stabilization_Type::NITSCHE_INTERFACE ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", std::string("100.0") );
        tParameterList( 2 )( tSPCounter ).set( "master_properties",   std::string("PropConductivity_Outer,Material") );
        tParameterList( 2 )( tSPCounter ).set( "slave_properties",    std::string("PropConductivity_Outer,Material") );
        tSPCounter++;

        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name",  std::string("SPInterfaceNitscheOuterInner") );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type",  static_cast< uint >( fem::Stabilization_Type::NITSCHE_INTERFACE ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", std::string("100.0") );
        tParameterList( 2 )( tSPCounter ).set( "master_properties",   std::string("PropConductivity_Outer,Material") );
        tParameterList( 2 )( tSPCounter ).set( "slave_properties",    std::string("PropConductivity_Inner,Material") );
        tSPCounter++;

        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name",  std::string("SPInterfaceMasterWeightOuter") );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type",  static_cast< uint >( fem::Stabilization_Type::MASTER_WEIGHT_INTERFACE ) );
        tParameterList( 2 )( tSPCounter ).set( "master_properties",   std::string("PropConductivity_Outer,Material") );
        tParameterList( 2 )( tSPCounter ).set( "slave_properties",    std::string("PropConductivity_Outer,Material") );
        tSPCounter++;

        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name",  std::string("SPInterfaceMasterWeightOuterInner") );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type",  static_cast< uint >( fem::Stabilization_Type::MASTER_WEIGHT_INTERFACE ) );
        tParameterList( 2 )( tSPCounter ).set( "master_properties",   std::string("PropConductivity_Outer,Material") );
        tParameterList( 2 )( tSPCounter ).set( "slave_properties",    std::string("PropConductivity_Inner,Material") );
        tSPCounter++;

        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name",  std::string("SPInterfaceSlaveWeightOuter") );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type",  static_cast< uint >( fem::Stabilization_Type::SLAVE_WEIGHT_INTERFACE ) );
        tParameterList( 2 )( tSPCounter ).set( "master_properties",   std::string("PropConductivity_Outer,Material") );
        tParameterList( 2 )( tSPCounter ).set( "slave_properties",    std::string("PropConductivity_Outer,Material") );
        tSPCounter++;

        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name",  std::string("SPInterfaceSlaveWeightOuterInner") );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type",  static_cast< uint >( fem::Stabilization_Type::SLAVE_WEIGHT_INTERFACE ) );
        tParameterList( 2 )( tSPCounter ).set( "master_properties",   std::string("PropConductivity_Outer,Material") );
        tParameterList( 2 )( tSPCounter ).set( "slave_properties",    std::string("PropConductivity_Inner,Material") );
        tSPCounter++;

        //------------------------------------------------------------------------------
        // init IWG counter
        uint tIWGCounter = 0;

        // create parameter  for outer ring
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   std::string("IWGDiffusionBulk") );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::SPATIALDIFF_BULK ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               std::string("TEMP") );
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    std::string("TEMP") );
        tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", std::string("CMDiffusion_Outer,Diffusion") );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             std::string("HMR_dummy_n_p4,HMR_dummy_n_p5,HMR_dummy_n_p7,HMR_dummy_c_p4,HMR_dummy_c_p5,HMR_dummy_c_p7") );
        tIWGCounter++;

        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   std::string("IWGDiffusionBulk") );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::SPATIALDIFF_BULK ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               std::string("TEMP") );
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    std::string("TEMP") );
        tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", std::string("CMDiffusion_Inner,Diffusion") );
	    tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   std::string("SPGGLSDiffusion_Inner,GGLS_Param") );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             std::string("HMR_dummy_n_p6,HMR_dummy_c_p6") );
        tIWGCounter++;
		
		tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   std::string("IWGInletFlux") );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::SPATIALDIFF_NEUMANN ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               std::string("TEMP") );
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    std::string("TEMP") );
        tParameterList( 3 )( tIWGCounter ).set( "master_properties",          std::string("PropImposedFlux,Neumann") );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             std::string("iside_b0_4_b1_0,iside_b0_5_b1_1") );
        tIWGCounter++;

        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   std::string("IWGInterfaceTEMPOuter") );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::SPATIALDIFF_INTERFACE ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               std::string("TEMP") );
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    std::string("TEMP") );
        tParameterList( 3 )( tIWGCounter ).set( "slave_dof_dependencies",     std::string("TEMP") );
        tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", std::string("CMDiffusion_Outer,Diffusion") );
        tParameterList( 3 )( tIWGCounter ).set( "slave_constitutive_models",  std::string("CMDiffusion_Outer_Dummy,Diffusion") );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",
                std::string("SPInterfaceNitscheOuter     , NitscheInterface;")      +
                std::string("SPInterfaceMasterWeightOuter, MasterWeightInterface;") +
                std::string("SPInterfaceSlaveWeightOuter , SlaveWeightInterface") );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             std::string("dbl_iside_p0_4_p1_5,dbl_iside_p0_5_p1_7") );
        tIWGCounter++;

        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   std::string("IWGInterfaceTEMPOuterInner") );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::SPATIALDIFF_INTERFACE ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               std::string("TEMP") );
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    std::string("TEMP") );
        tParameterList( 3 )( tIWGCounter ).set( "slave_dof_dependencies",     std::string("TEMP") );
        tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", std::string("CMDiffusion_Outer,Diffusion") );
        tParameterList( 3 )( tIWGCounter ).set( "slave_constitutive_models",  std::string("CMDiffusion_Inner,Diffusion") );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",
                std::string("SPInterfaceNitscheOuterInner     , NitscheInterface;")      +
                std::string("SPInterfaceMasterWeightOuterInner, MasterWeightInterface;") +
                std::string("SPInterfaceSlaveWeightOuterInner , SlaveWeightInterface") );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             std::string("dbl_iside_p0_7_p1_6,dbl_iside_p0_4_p1_6") );
        tIWGCounter++;

        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   std::string("IWGGPTempOuter") );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::SPATIALDIFF_GHOST ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               std::string("TEMP") );
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    std::string("TEMP") );
        tParameterList( 3 )( tIWGCounter ).set( "slave_dof_dependencies",     std::string("TEMP") );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   std::string("SPGPTemp_Outer,GhostDispl") );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             std::string("ghost_p4,ghost_p5,ghost_p7") );
        tIWGCounter++;

        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   std::string("IWGGPTempInner") );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::SPATIALDIFF_GHOST ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               std::string("TEMP") );
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    std::string("TEMP") );
        tParameterList( 3 )( tIWGCounter ).set( "slave_dof_dependencies",     std::string("TEMP") );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   std::string("SPGPTemp_Inner,GhostDispl") );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             std::string("ghost_p6") );
        tIWGCounter++;

        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   std::string("IWGTimeContinuityTemp") );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::TIME_CONTINUITY_DOF ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               std::string("TEMP") );
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    std::string("TEMP") );
        tParameterList( 3 )( tIWGCounter ).set( "master_properties",          
	        std::string("PropWeightCurrent   , WeightCurrent;")   +
	        std::string("PropWeightPrevious  , WeightPrevious;")  + 
		    std::string("PropInitialCondition, InitialCondition") );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",
                std::string("HMR_dummy_n_p4,HMR_dummy_c_p4,") +
                std::string("HMR_dummy_n_p5,HMR_dummy_c_p5,") +
                std::string("HMR_dummy_n_p6,HMR_dummy_c_p6,") +
                std::string("HMR_dummy_n_p7,HMR_dummy_c_p7")  );
        tParameterList( 3 )( tIWGCounter ).set( "time_continuity",            true );
        tIWGCounter++;

        //------------------------------------------------------------------------------
        // init IQI counter
        uint tIQICounter = 0;

        // create parameter list for IQI 4
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name",                   std::string("IQIBulkTEMP") );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::DOF ) );
        tParameterList( 4 )( tIQICounter ).set( "IQI_output_type",            static_cast< uint >( vis::Output_Type::TEMP ) );
        tParameterList( 4 )( tIQICounter ).set( "master_dof_dependencies",    std::string("TEMP") );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index",      0 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names",
                std::string("HMR_dummy_n_p4,HMR_dummy_c_p4,") +
                std::string("HMR_dummy_n_p5,HMR_dummy_c_p5,") +
                std::string("HMR_dummy_n_p6,HMR_dummy_c_p6,") +
                std::string("HMR_dummy_n_p7,HMR_dummy_c_p7")  );
        tIQICounter++;
    }

    void SOLParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 7 );
        for( uint Ik = 0; Ik < 7; Ik ++)
        {
            tParameterlist( Ik ).resize( 1 );
        }

        tParameterlist( 0 )( 0 ) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL );

        tParameterlist( 1 )( 0 ) = moris::prm::create_linear_solver_parameter_list();

        tParameterlist( 2 )( 0 ) = moris::prm::create_nonlinear_algorithm_parameter_list();
        tParameterlist( 2 )( 0 ).set("NLA_rel_res_norm_drop",    1e-04 );
        tParameterlist( 2 )( 0 ).set("NLA_relaxation_parameter", 1.00  );
        tParameterlist( 2 )( 0 ).set("NLA_max_iter", 20 );

        tParameterlist( 3 )( 0 ) = moris::prm::create_nonlinear_solver_parameter_list();
        tParameterlist( 3 )( 0 ).set("NLA_DofTypes"      , std::string("TEMP") );

        tParameterlist( 4 )( 0 ) = moris::prm::create_time_solver_algorithm_parameter_list();
        tParameterlist( 4 )( 0 ).set("TSA_Num_Time_Steps",     50 );
        tParameterlist( 4 )( 0 ).set("TSA_Time_Frame",         0.5 );

        tParameterlist( 5 )( 0 ) = moris::prm::create_time_solver_parameter_list();
        tParameterlist( 5 )( 0 ).set("TSA_DofTypes",           std::string("TEMP") );
        tParameterlist( 5 )( 0 ).set("TSA_Initialize_Sol_Vec", std::string("TEMP,0.0") );
        tParameterlist( 5 )( 0 ).set("TSA_Output_Indices",     std::string("0") );
        tParameterlist( 5 )( 0 ).set("TSA_Output_Crteria",     std::string("Output_Criterion") );
        tParameterlist( 5 )( 0 ).set("TSA_time_level_per_type",std::string("TEMP,2") );

        tParameterlist( 6 )( 0 ) = moris::prm::create_solver_warehouse_parameterlist();
    }

    void MSIParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_msi_parameter_list();
    }

    void VISParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_vis_parameter_list();
        tParameterlist( 0 )( 0 ).set( "File_Name"  , std::pair< std::string, std::string >( "./", "Vessel_with_Fins_3D.exo" ) );
        tParameterlist( 0 )( 0 ).set( "Mesh_Type"  , static_cast< uint >( vis::VIS_Mesh_Type::STANDARD ) );
        tParameterlist( 0 )( 0 ).set( "Set_Names"  ,
                std::string("HMR_dummy_n_p4,HMR_dummy_c_p4,") +
                std::string("HMR_dummy_n_p5,HMR_dummy_c_p5,") +
                std::string("HMR_dummy_n_p6,HMR_dummy_c_p6,") +
                std::string("HMR_dummy_n_p7,HMR_dummy_c_p7")  );
        tParameterlist( 0 )( 0 ).set( "Field_Names", std::string( "TEMP" ) );
        tParameterlist( 0 )( 0 ).set( "Field_Type" , std::string( "NODAL" ) );
        tParameterlist( 0 )( 0 ).set( "Output_Type", std::string( "TEMP" ) );
    }

    //------------------------------------------------------------------------------
}

//------------------------------------------------------------------------------
#ifdef  __cplusplus
}
#endif
