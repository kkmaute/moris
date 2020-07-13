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
    moris::real tXlength = 0.28;
    moris::real tYlength = 0.0015;
    moris::real tXcenter = 0.5 * tXlength;
    moris::real tYcenter = 0.5 * tYlength;
    moris::real tEps = 1.0e-4;


    // time solver parameters
    moris::sint tStep = 10;
    moris::real tTmax = 480.0;

    // ramp up of Dirichlet BC (number of time slabs to ramp up the value on the BC)
    moris::real tRampUp = 6.0;


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
        aPropMatrix = {{ 313.0 }};
    }

    void Func_Wall_Condition(
            moris::Matrix< moris::DDRMat >                 & aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
            moris::fem::Field_Interpolator_Manager 	* aFIManager )
    {
        real tT = aFIManager->get_IP_geometry_interpolator()->valt()( 0 );

        real value = std::min(350.0, 313.0+(350-313)*tT/(tRampUp*tTmax/tStep));

        aPropMatrix      = {{ 0.0 }};
        aPropMatrix(0,0) = value;
    }

    bool Output_Criterion( moris::tsa::Time_Solver * aTimeSolver )
    {
        return true;
    }

    // Dummy level set function for conform mesh
    moris::real Dummy_LS(
            const moris::Matrix< DDRMat >     & aCoordinates,
            const moris::Cell< moris::real* > & aGeometryParameters )
    {
        moris::real tLSval = aCoordinates(0) + 1.0  ;

        // clean return value to return non-zero value
        return std::abs(tLSval) < 1.0e-8  ?  1.0e-8 : tLSval;
    }

    moris::Matrix< DDRMat > Func_Sensitivity(
            const moris::Matrix< DDRMat >     & aCoordinates,
            const moris::Cell< moris::real* > & aGeometryParameters )
        {
        moris::Matrix< DDRMat > aReturnValue;
        return aReturnValue;
        }

    void HMRParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_hmr_parameter_list();

        tParameterlist( 0 )( 0 ).set( "number_of_elements_per_dimension", std::string("500,1"));
        tParameterlist( 0 )( 0 ).set( "domain_dimensions",                std::string("0.28,0.0005"));
        tParameterlist( 0 )( 0 ).set( "domain_offset",                    std::string("0.0,0.0") );
        tParameterlist( 0 )( 0 ).set( "domain_sidesets",                  std::string("1,2,3,4"));
        tParameterlist( 0 )( 0 ).set( "lagrange_output_meshes",           std::string("0"));

        tParameterlist( 0 )( 0 ).set( "lagrange_orders",  std::string( "1" ));
        tParameterlist( 0 )( 0 ).set( "lagrange_pattern", std::string( "0" ));
        tParameterlist( 0 )( 0 ).set( "bspline_orders",   std::string( "1" ));
        tParameterlist( 0 )( 0 ).set( "bspline_pattern",  std::string( "0" ));

        tParameterlist( 0 )( 0 ).set( "lagrange_to_bspline", std::string("0") );

        tParameterlist( 0 )( 0 ).set( "truncate_bsplines",  1 );
        tParameterlist( 0 )( 0 ).set( "refinement_buffer",  3 );
        tParameterlist( 0 )( 0 ).set( "staircase_buffer",   3 );
        tParameterlist( 0 )( 0 ).set( "initial_refinement", 0 );

        tParameterlist( 0 )( 0 ).set( "use_multigrid",  0 );
        tParameterlist( 0 )( 0 ).set( "severity_level", 0 );

        tParameterlist( 0 )( 0 ).set( "adaptive_refinement_level", 1 );
    }

    void OPTParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_opt_problem_parameter_list();

        tParameterlist( 0 )( 0 ).set( "is_optimization_problem", false);
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
        tParameterlist( 0 )( 0 ).set( "exodus_output_XTK_ig_mesh", true );
        tParameterlist( 0 )( 0 ).set( "high_to_low_dbl_side_sets", true );
    }

    void GENParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 3 );
        tParameterlist( 0 ).resize( 1 );

        // Main GEN parameter list
        tParameterlist( 0 )( 0 ) = prm::create_gen_parameter_list();

        // init geometry counter
        uint tGeoCounter = 0;

        // Geometry parameter lists
        tParameterlist( 1 ).push_back( prm::create_user_defined_geometry_parameter_list() );
        tParameterlist( 1 )( tGeoCounter ).set( "field_function_name", "Dummy_LS");
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

        // create parameter list for property 1
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ) = prm::create_property_parameter_list();
        tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropDensity") );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      std::string("0.75") );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Const") );
        tPropCounter++;

        // create parameter list for property 2
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropConductivity") );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      std::string("2.1e-7") );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Const") );
        tPropCounter++;

        // create parameter list for property 3
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropHeatCapacity") );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      std::string("2.4") );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Const") );
        tPropCounter++;

        // create parameter list for property 3
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropLatentHeat") );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      std::string("175") );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Const") );
        tPropCounter++;

        // create parameter list for property 9
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropPCTemp") );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      std::string("314.5") );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Const") );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropPhaseState") );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      std::string("2.0") );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Const") );
        tPropCounter++;

        // create parameter list for property 11
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropPCconst") );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      std::string("3.0") );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Const") );
        tPropCounter++;

        // create parameter list for property 4
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropImposedTemp") );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Wall_Condition") );
        tPropCounter++;

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

        // create parameter list for constitutive model 1
        tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_name", std::string("CMDiffusionPhaseChange") );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::DIFF_LIN_ISO_PC ) );
        tParameterList( 1 )( tCMCounter ).set( "dof_dependencies",  std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( 1 )( tCMCounter ).set( "properties",
                std::string("PropConductivity,Conductivity;") +
                std::string("PropDensity,Density;") +
                std::string("PropHeatCapacity,HeatCapacity;") +
                std::string("PropLatentHeat,LatentHeat;") +
                std::string("PropPCTemp,PCTemp;") +
                std::string("PropPhaseState,PhaseStateFunction;") +
                std::string("PropPCconst,PhaseChangeConst") );
        tCMCounter++;

        //------------------------------------------------------------------------------

        // init SP counter
        uint tSPCounter = 0;

        // Ghost
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name",      std::string("SPGPTemp") );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type",      static_cast< uint >( fem::Stabilization_Type::GHOST_DISPL ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters",     std::string("0.01") );
        tParameterList( 2 )( tSPCounter ).set( "master_properties",       std::string("PropConductivity,Material") );
        tSPCounter++;

        // GGLS parameter
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name",      std::string("SPGGLSDiffusion") );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type",      static_cast< uint >( fem::Stabilization_Type::GGLS_DIFFUSION ) );
        tParameterList( 2 )( tSPCounter ).set( "master_dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( 2 )( tSPCounter ).set( "master_properties",
                std::string("PropConductivity , Conductivity;")         +
                std::string("PropDensity      , Density;")              +
                std::string("PropHeatCapacity , HeatCapacity;")        +
                std::string("PropLatentHeat   , LatentHeat;")          +
                std::string("PropPCTemp       , PCTemp;")              +
                std::string("PropPhaseState   , PhaseStateFunction;") +
                std::string("PropPCconst      , PhaseChangeConst")    );
        tSPCounter++;

        // Dirichlet SP
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name",      std::string("SPNitscheTemp") );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type",      static_cast< uint >( fem::Stabilization_Type::DIRICHLET_NITSCHE ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters",     std::string("1000.0") );
        tParameterList( 2 )( tSPCounter ).set( "master_properties",       std::string("PropConductivity,Material") );
        tSPCounter++;

        //------------------------------------------------------------------------------
        // init IWG counter
        uint tIWGCounter = 0;

        // Bulk
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   std::string("IWGDiffusionBulk") );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::SPATIALDIFF_BULK ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               std::string("TEMP") );
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    std::string("TEMP") );
        tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", std::string("CMDiffusionPhaseChange,Diffusion") );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   std::string("SPGGLSDiffusion,GGLSParam") );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             std::string("HMR_dummy_n_p1") );
        tIWGCounter++;

        // Dirichlet BC
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", 		          std::string("IWGOutletTemp") );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", 		          static_cast< uint >( fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               std::string("TEMP") );
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    std::string("TEMP") );
        tParameterList( 3 )( tIWGCounter ).set( "master_properties",	      std::string("PropImposedTemp,Dirichlet") );
        tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", std::string("CMDiffusionPhaseChange,Diffusion") );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   std::string("SPNitscheTemp,DirichletNitsche") );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",		      std::string("SideSet_4_n_p1") );
        tIWGCounter++;

        // Time Continuity
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   std::string("IWGTimeContinuityTemp") );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::TIME_CONTINUITY_DOF ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               std::string("TEMP") );
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    std::string("TEMP") );
        tParameterList( 3 )( tIWGCounter ).set( "master_properties",
                std::string("PropWeightCurrent,WeightCurrent;") +
                std::string("PropWeightPrevious,WeightPrevious;") +
                std::string("PropInitialCondition,InitialCondition") );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             std::string("HMR_dummy_n_p1") );
        tParameterList( 3 )( tIWGCounter ).set( "time_continuity",            true );
        tIWGCounter++;

        //------------------------------------------------------------------------------
        // init IQI counter
        uint tIQICounter = 0;

        // Temperature IQI
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name",                   std::string("IQIBulkTEMP") );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::DOF ) );
        tParameterList( 4 )( tIQICounter ).set( "IQI_output_type",            static_cast< uint >( vis::Output_Type::TEMP ) );
        tParameterList( 4 )( tIQICounter ).set( "master_dof_dependencies",    std::string("TEMP") );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index",      0 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names",             std::string("HMR_dummy_n_p1") );
        tIQICounter++;

        // Max Temperature IQI
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name",                   std::string("IQIMaxTemp") );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::MAX_DOF ) );
        tParameterList( 4 )( tIQICounter ).set( "IQI_output_type",            static_cast< uint >( vis::Output_Type::MAX_DOF ) );
        tParameterList( 4 )( tIQICounter ).set( "master_dof_dependencies",    std::string("TEMP") );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index",      0 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names",             std::string("HMR_dummy_n_p1") );
        tIQICounter++;

        // Latent Heat IQI
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name",                   std::string("IQILatentHeatAbsorption") );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::LATENT_HEAT_ABSORPTION ) );
        tParameterList( 4 )( tIQICounter ).set( "IQI_output_type",            static_cast< uint >( vis::Output_Type::LATENT_HEAT_ABSORPTION ) );
        tParameterList( 4 )( tIQICounter ).set( "master_dof_dependencies",    std::string("TEMP") );
        tParameterList( 4 )( tIQICounter ).set( "master_properties",
                std::string("PropDensity,Density;") +
                std::string("PropLatentHeat,LatentHeat;") +
                std::string("PropPCTemp,PCTemp;") +
                std::string("PropPhaseState,PhaseStateFunction;") +
                std::string("PropPCconst,PhaseChangeConst") );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names",             std::string("HMR_dummy_n_p1") );
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
        tParameterlist( 2 )( 0 ).set("NLA_rel_res_norm_drop", 5e-05 );
        tParameterlist( 2 )( 0 ).set("NLA_relaxation_parameter", 0.9 );
        tParameterlist( 2 )( 0 ).set("NLA_max_iter", 20 );

        tParameterlist( 3 )( 0 ) = moris::prm::create_nonlinear_solver_parameter_list();
        tParameterlist( 3 )( 0 ).set("NLA_DofTypes"      , std::string("TEMP") );

        tParameterlist( 4 )( 0 ) = moris::prm::create_time_solver_algorithm_parameter_list();
        tParameterlist( 4 )( 0 ).set("TSA_Num_Time_Steps",     tStep );
        tParameterlist( 4 )( 0 ).set("TSA_Time_Frame",         tTmax );

        tParameterlist( 5 )( 0 ) = moris::prm::create_time_solver_parameter_list();
        tParameterlist( 5 )( 0 ).set("TSA_DofTypes",            std::string("TEMP") );
        tParameterlist( 5 )( 0 ).set("TSA_Initialize_Sol_Vec",  std::string("TEMP,0.0") );
        tParameterlist( 5 )( 0 ).set("TSA_Output_Indices",      std::string("0") );
        tParameterlist( 5 )( 0 ).set("TSA_Output_Crteria",      std::string("Output_Criterion") );
        tParameterlist( 5 )( 0 ).set("TSA_time_level_per_type", std::string("TEMP,2") );

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
        tParameterlist( 0 )( 0 ).set( "File_Name"  , std::pair< std::string, std::string >( "./", "Comsol_conform.exo" ) );
        tParameterlist( 0 )( 0 ).set( "Mesh_Type"  , static_cast< uint >( vis::VIS_Mesh_Type::STANDARD ) );
        tParameterlist( 0 )( 0 ).set( "Set_Names"  , std::string( "HMR_dummy_n_p1" ) );
        tParameterlist( 0 )( 0 ).set( "Field_Names", std::string( "TEMP,MAX_DOF,LATENT_HEAT_ABSORPTION" ) );
        tParameterlist( 0 )( 0 ).set( "Field_Type" , std::string( "NODAL,GLOBAL,GLOBAL" ) );
        tParameterlist( 0 )( 0 ).set( "Output_Type", std::string( "TEMP,MAX_DOF,LATENT_HEAT_ABSORPTION" ) );
    }

    //------------------------------------------------------------------------------
}

//------------------------------------------------------------------------------
#ifdef  __cplusplus
}
#endif
