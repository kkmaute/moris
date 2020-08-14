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
#include "cl_HMR_Element.hpp"
#include "fn_equal_to.hpp"

#include "AztecOO.h"

//---------------------------------------------------------------

// global variable for interpolation order
extern uint gInterpolationOrder;

//---------------------------------------------------------------

#ifdef  __cplusplus
extern "C"
{
#endif
//------------------------------------------------------------------------------
namespace moris
{
    /* ------------------------------------------------------------------------ */
    // Mesh Set Information

    std::string tInnerPhase        = "HMR_dummy_n_p1,HMR_dummy_c_p1";
    std::string tOuterPhase        = "HMR_dummy_n_p0,HMR_dummy_c_p0";

    std::string tInterface         = "dbl_iside_p0_1_p1_0";

    std::string tOuterSurface      = "SideSet_1_n_p0,SideSet_2_n_p0,SideSet_3_n_p0,SideSet_4_n_p0,SideSet_5_n_p0,SideSet_6_n_p0";

    std::string tInnerPhaseGhost   = "ghost_p1";
    std::string tOuterPhaseGhost   = "ghost_p0";

    std::string tTotalDomain       = tInnerPhase + "," + tOuterPhase;

    /* ------------------------------------------------------------------------ */    
    // geometry parameters

    // general
    moris::real tCenterX = 0.0;
    moris::real tCenterY = 0.0;
    moris::real tCenterZ = 0.0;
    moris::real tRadius  = 0.5;

    /* ------------------------------------------------------------------------ */    
    // material parameters

    // capacity
    std::string tCapacityInner = "0.0"; 
    std::string tCapacityOuter = "0.0";

    // density
    std::string tDensityInner = "0.0"; 
    std::string tDensityOuter = "0.0";

    // conductivity
    std::string tConductivityInner = "1.0"; 
    std::string tConductivityOuter = "0.125";

    // body flux
    std::string tHeatLoadInner = "1.0";
    std::string tHeatLoadOuter = "0.0";

    /* ------------------------------------------------------------------------ */
    // HMR parameters

    std::string tNumElemsPerDim     = "8, 8, 8";
    std::string tDomainDims         = "2.0, 2.0, 2.0";
    std::string tDomainOffset       = "-1.0, -1.0, -1.0";
    std::string tDomainSidesets     = "1,2,3,4,5,6";

    int tRefineBuffer      = 1;
    int tInitialRefinement = 1;

    int tInterfaceRefinement = 0;

    /* ------------------------------------------------------------------------ */
    // Solver config

    moris::real tNLA_rel_res_norm_drop = 1.0e-08;
    moris::real tNLA_relaxation_parameter = 1.0;
    int tNLA_max_iter = 2;

    int tTSA_Num_Time_Steps = 1;
    moris::real tTSA_Time_Frame = 1.0e0;

    /* ------------------------------------------------------------------------ */
    // Minimum level set value
    moris::real tMinLevs = 1.0e-8;

    /* ------------------------------------------------------------------------ */
    // Minimum level set value
    bool tUseGhost = true;

    /* ------------------------------------------------------------------------ */
    // Output Config

    std::string tOutputFileName = "Heated_Sphere.exo";

    /* ------------------------------------------------------------------------ */

    // Level set function for diamond shaped wedge
    moris::real Inclusion(
            const moris::Matrix< DDRMat >     & aCoordinates,
            const moris::Cell< moris::real* > & aGeometryParameters )
    {        
        // distance from sphere center
        moris::real tDx = aCoordinates(0) - tCenterX;
        moris::real tDy = aCoordinates(1) - tCenterY;
        moris::real tDz = aCoordinates(2) - tCenterZ;

        // Compute Signed-Distance field 
        moris::real tVal = tRadius - std::sqrt( tDx*tDx + tDy*tDy + tDz*tDz);

        // clean return value to return non-zero value
        return std::abs(tVal) < tMinLevs ? tMinLevs : tVal;
    }

    /* ------------------------------------------------------------------------ */

    // Constant function for properties
    void Func_Const( moris::Matrix<
            moris::DDRMat >                                & aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
            moris::fem::Field_Interpolator_Manager         * aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    void Func_Exact_Temperature(
            moris::Matrix< moris::DDRMat >                 & aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
            moris::fem::Field_Interpolator_Manager         * aFIManager )
    {
        // distance from sphere center
        moris::real tDx = aFIManager->get_IP_geometry_interpolator()->valx()( 0 ) - tCenterX;
        moris::real tDy = aFIManager->get_IP_geometry_interpolator()->valx()( 1 ) - tCenterY;
        moris::real tDz = aFIManager->get_IP_geometry_interpolator()->valx()( 2 ) - tCenterZ;

        // Compute Signed-Distance field 
        moris::real distance = std::sqrt( tDx*tDx + tDy*tDy + tDz*tDz);

        if ( distance <= tRadius )
        {
            aPropMatrix = {{ 3.0/8.0-distance*distance/6.0 }};
        }
        else
        {
            aPropMatrix = {{ (1.0/3.0)*(1.0/distance-1.0) }};
        }
    }

    /* ------------------------------------------------------------------------ */

    void Func_Exact_TemperatureGradient(
            moris::Matrix< moris::DDRMat >                 & aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
            moris::fem::Field_Interpolator_Manager         * aFIManager )
    {
        // distance from sphere center
        moris::real tDx = aFIManager->get_IP_geometry_interpolator()->valx()( 0 ) - tCenterX;
        moris::real tDy = aFIManager->get_IP_geometry_interpolator()->valx()( 1 ) - tCenterY;
        moris::real tDz = aFIManager->get_IP_geometry_interpolator()->valx()( 2 ) - tCenterZ;

        // Compute Signed-Distance field
        moris::real distance = std::sqrt( tDx*tDx + tDy*tDy + tDz*tDz);

        moris::real ddistdx = tDx/distance;
        moris::real ddistdy = tDy/distance;
        moris::real ddistdz = tDz/distance;

        // set size for aPropMatrix
        aPropMatrix.set_size( 3, 1, 0.0 );

        // spatial gradients of analytic temperature distribution
        if ( distance <= tRadius )
        {
            aPropMatrix( 0, 0 ) = -2.0/6.0*distance*ddistdx;
            aPropMatrix( 1, 0 ) = -2.0/6.0*distance*ddistdy;
            aPropMatrix( 2, 0 ) = -2.0/6.0*distance*ddistdz;
        }
        else {
            aPropMatrix( 0, 0 ) = -(1.0/3.0)*(1.0/distance/distance*ddistdx);
            aPropMatrix( 1, 0 ) = -(1.0/3.0)*(1.0/distance/distance*ddistdy);
            aPropMatrix( 2, 0 ) = -(1.0/3.0)*(1.0/distance/distance*ddistdz);
        }
    }

    /* ------------------------------------------------------------------------ */

    bool Output_Criterion( moris::tsa::Time_Solver * aTimeSolver )
    {
        return true;
    }

    /* ------------------------------------------------------------------------ */

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

        tParameterlist( 0 )( 0 ).set( "number_of_elements_per_dimension", tNumElemsPerDim );
        tParameterlist( 0 )( 0 ).set( "domain_dimensions",                tDomainDims );
        tParameterlist( 0 )( 0 ).set( "domain_offset",                    tDomainOffset );
        tParameterlist( 0 )( 0 ).set( "domain_sidesets",                  tDomainSidesets);
        tParameterlist( 0 )( 0 ).set( "lagrange_output_meshes",           std::string("0"));

        tParameterlist( 0 )( 0 ).set( "lagrange_orders",  std::to_string( gInterpolationOrder ));
        tParameterlist( 0 )( 0 ).set( "lagrange_pattern", std::string( "0" ));
        tParameterlist( 0 )( 0 ).set( "bspline_orders",   std::to_string( gInterpolationOrder));
        tParameterlist( 0 )( 0 ).set( "bspline_pattern",  std::string( "0" ));

        tParameterlist( 0 )( 0 ).set( "lagrange_to_bspline", std::string("0") );

        tParameterlist( 0 )( 0 ).set( "truncate_bsplines",  1 );
        tParameterlist( 0 )( 0 ).set( "refinement_buffer",  tRefineBuffer );
        tParameterlist( 0 )( 0 ).set( "staircase_buffer",   tRefineBuffer );
        tParameterlist( 0 )( 0 ).set( "initial_refinement", tInitialRefinement );

        tParameterlist( 0 )( 0 ).set( "use_number_aura", 1);

        tParameterlist( 0 )( 0 ).set( "use_multigrid",  0 );
        tParameterlist( 0 )( 0 ).set( "severity_level", 0 );
//
//        tParameterlist( 0 )( 0 ).set( "lagrange_input_meshes", "0");
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
        tParameterlist( 0 )( 0 ).set( "ghost_stab",                tUseGhost );
        tParameterlist( 0 )( 0 ).set( "multigrid",                 false );
        tParameterlist( 0 )( 0 ).set( "verbose",                   true );
        tParameterlist( 0 )( 0 ).set( "print_enriched_ig_mesh",    false );
        tParameterlist( 0 )( 0 ).set( "exodus_output_XTK_ig_mesh", false );
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
        tParameterlist( 1 )( tGeoCounter ).set( "field_function_name",       "Inclusion" );
        tParameterlist( 1 )( tGeoCounter ).set( "number_of_refinements",      tInterfaceRefinement );
        tGeoCounter++;

    }

    void FEMParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterList )
    {
        // create a cell of cell of parameter list for fem
        tParameterList.resize( 5 );

        //------------------------------------------------------------------------------
        // init property counter
        uint tPropCounter = 0;

        // properties for inclusion

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ) = prm::create_property_parameter_list();
        tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropDensityInner") );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      tDensityInner );               
        tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Const") );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropCapacityInner") );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      tCapacityInner );                 
        tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Const") );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropConductivityInner") );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      tConductivityInner );              
        tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Const") );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropHeatLoadInner") );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      tHeatLoadInner );              
        tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Const") );
        tPropCounter++;

        // properties for outer material

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ) = prm::create_property_parameter_list();
        tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropDensityOuter") );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      tDensityOuter );               
        tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Const") );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropCapacityOuter") );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      tCapacityOuter );                 
        tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Const") );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropConductivityOuter") );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      tConductivityOuter );              
        tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Const") );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropHeatLoadOuter") );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      tHeatLoadOuter );              
        tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Const") );
        tPropCounter++;

        // time continuity weights        
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropWeightCurrent") );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      std::string("100.0") );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Const") );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropWeightPrevious") );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      std::string("100.0") );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Const") );
        tPropCounter++;

        // initial condition       
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropInitialCondition") );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Exact_Temperature") );
        tPropCounter++;

        // temperature at outer surface condition       
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropImposedTemperature") );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Exact_Temperature") );
        tPropCounter++;

        // create parameter list for property 4
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropExactTemperature") );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Exact_Temperature") );
        tPropCounter++;

        // create parameter list for property 5
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropExactTemperatureGradient") );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Exact_TemperatureGradient") );
        tPropCounter++;

        //------------------------------------------------------------------------------
        // init CM counter
        uint tCMCounter = 0;

        // create parameter list for constitutive model - Inclusion
        tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_name", std::string("CMDiffusionInner") );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::DIFF_LIN_ISO ) );
        tParameterList( 1 )( tCMCounter ).set( "dof_dependencies",  std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( 1 )( tCMCounter ).set( "properties",
                std::string("PropConductivityInner , Conductivity;") +
                std::string("PropDensityInner      , Density;")      +
                std::string("PropCapacityInner     , HeatCapacity") );
        tCMCounter++;

        // create parameter list for constitutive model - Outer Material
        tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_name", std::string("CMDiffusionOuter") );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::DIFF_LIN_ISO ) );
        tParameterList( 1 )( tCMCounter ).set( "dof_dependencies",  std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( 1 )( tCMCounter ).set( "properties",
                std::string("PropConductivityOuter , Conductivity;") +
                std::string("PropDensityOuter      , Density;")      +
                std::string("PropCapacityOuter     , HeatCapacity") );
        tCMCounter++;

        //------------------------------------------------------------------------------
        // init SP counter
        uint tSPCounter = 0;

        // create parameter list for ghost stabilization parameter for inclusion
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name",      std::string("SPGPTempInner") );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type",      static_cast< uint >( fem::Stabilization_Type::GHOST_DISPL ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters",     std::string("0.01") );
        tParameterList( 2 )( tSPCounter ).set( "master_properties",       std::string("PropConductivityInner,Material") );
        tSPCounter++;

        // create parameter list for ghost stabilization parameter for outer material
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name",      std::string("SPGPTempOuter") );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type",      static_cast< uint >( fem::Stabilization_Type::GHOST_DISPL ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters",     std::string("0.01") );
        tParameterList( 2 )( tSPCounter ).set( "master_properties",       std::string("PropConductivityOuter,Material") );
        tSPCounter++;

        // create parameter list for Nitsche stabilization parameter for inclusion-outer material interface
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name",  std::string("SPInterfaceNitsche") );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type",  static_cast< uint >( fem::Stabilization_Type::NITSCHE_INTERFACE ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", std::string("100.0") );
        tParameterList( 2 )( tSPCounter ).set( "master_properties",   std::string("PropConductivityInner,Material") );
        tParameterList( 2 )( tSPCounter ).set( "slave_properties",    std::string("PropConductivityOuter,Material") );
        tSPCounter++;

        // create parameter list for DBC on outer surface
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name",      std::string("SPNitscheTemp") );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type",      static_cast< uint >( fem::Stabilization_Type::DIRICHLET_NITSCHE ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters",     std::string("100.0") );
        tParameterList( 2 )( tSPCounter ).set( "master_properties",       std::string("PropConductivityOuter,Material") );
        tSPCounter++;

        //------------------------------------------------------------------------------
        // init IWG counter
        uint tIWGCounter = 0;

        // create IWG for inclusion - bulk diffusion
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   std::string("IWGDiffusionInnerBulk") );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::SPATIALDIFF_BULK ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               std::string("TEMP") );
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    std::string("TEMP") );
        tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", std::string("CMDiffusionInner,Diffusion") );
        tParameterList( 3 )( tIWGCounter ).set( "master_properties",          std::string("PropHeatLoadInner,Load") );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             tInnerPhase );
        tIWGCounter++;   

        // create IWG for outer material - bulk diffusion
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   std::string("IWGDiffusionOuterBulk") );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::SPATIALDIFF_BULK ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               std::string("TEMP") );
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    std::string("TEMP") );
        tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", std::string("CMDiffusionOuter,Diffusion") );
        tParameterList( 3 )( tIWGCounter ).set( "master_properties",          std::string("PropHeatLoadOuter,Load") );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             tOuterPhase );
        tIWGCounter++;   

        // create parameter list for outer boundary conditions
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   std::string("IWGOuterSurfaceTemp") );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               std::string("TEMP") );
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    std::string("TEMP") );
        tParameterList( 3 )( tIWGCounter ).set( "master_properties",          std::string("PropImposedTemperature,Dirichlet") );
        tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", std::string("CMDiffusionOuter,Diffusion") );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   std::string("SPNitscheTemp,DirichletNitsche") );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             tOuterSurface );
        tIWGCounter++;

        // create parameter list for interface conditions
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   std::string("IWGInterfaceInnerOuterTEMP") );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::SPATIALDIFF_INTERFACE ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               std::string("TEMP") );
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    std::string("TEMP") );
        tParameterList( 3 )( tIWGCounter ).set( "slave_dof_dependencies",     std::string("TEMP") );
        tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", std::string("CMDiffusionInner,Diffusion") );
        tParameterList( 3 )( tIWGCounter ).set( "slave_constitutive_models",  std::string("CMDiffusionOuter,Diffusion") );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   std::string("SPInterfaceNitsche,NitscheInterface") );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             tInterface );
        tIWGCounter++;

        if (tUseGhost)
        {
            // create IWG for outer material - ghost
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   std::string("IWGGPInnerTemp") );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::GHOST_NORMAL_FIELD ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               std::string("TEMP") );
            tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    std::string("TEMP") );
            tParameterList( 3 )( tIWGCounter ).set( "slave_dof_dependencies",     std::string("TEMP") );
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   std::string("SPGPTempInner,GhostSP") );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             tInnerPhaseGhost );
            tIWGCounter++;

            // create IWG for outer material - ghost
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   std::string("IWGGPInnerTemp") );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::GHOST_NORMAL_FIELD ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               std::string("TEMP") );
            tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    std::string("TEMP") );
            tParameterList( 3 )( tIWGCounter ).set( "slave_dof_dependencies",     std::string("TEMP") );
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   std::string("SPGPTempOuter,GhostSP") );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             tOuterPhaseGhost );
            tIWGCounter++;
        }

//        // create IWG for time continuity
//        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
//        tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   std::string("IWGTimeContinuityTemp") );
//        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::TIME_CONTINUITY_DOF ) );
//        tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               std::string("TEMP") );
//        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    std::string("TEMP") );
//        tParameterList( 3 )( tIWGCounter ).set( "master_properties",
//                std::string("PropWeightCurrent   ,WeightCurrent;")   +
//                std::string("PropWeightPrevious  ,WeightPrevious;")  +
//                std::string("PropInitialCondition,InitialCondition") );
//        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             tTotalDomain );
//        tParameterList( 3 )( tIWGCounter ).set( "time_continuity",            true );
//        tIWGCounter++;


        //------------------------------------------------------------------------------
        // init IQI counter
        uint tIQICounter = 0;

        // Nodal Temperature IQI
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name",                   std::string("IQIBulkTEMP") );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::DOF ) );
        tParameterList( 4 )( tIQICounter ).set( "IQI_output_type",            static_cast< uint >( vis::Output_Type::TEMP ) );
        tParameterList( 4 )( tIQICounter ).set( "master_dof_dependencies",    std::string("TEMP") );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index",      0 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names",             tTotalDomain );
        tIQICounter++; 

        // Nodal Analytic Temperature IQI
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name",                   std::string("IQIBulkTEMPAnalytic") );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::PROPERTY ) );
        tParameterList( 4 )( tIQICounter ).set( "IQI_output_type",            static_cast< uint >( vis::Output_Type::PROPERTY ) );
        tParameterList( 4 )( tIQICounter ).set( "master_properties",          std::string("PropExactTemperature,Property") );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names",             tTotalDomain );
        tIQICounter++;

        // L2 Error IQI
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name",                   std::string("IQIBulkL2Error") );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::L2_ERROR_ANALYTIC ) );
        tParameterList( 4 )( tIQICounter ).set( "IQI_output_type",            static_cast< uint >( vis::Output_Type::L2_ERROR_ANALYTIC ) );
        tParameterList( 4 )( tIQICounter ).set( "master_dof_dependencies",    std::string("TEMP") );
        tParameterList( 4 )( tIQICounter ).set( "master_properties",          std::string("PropExactTemperature,L2Check") );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names",             tTotalDomain );
        tIQICounter++;

        // H1 Semi-Norm Error IQI
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name",                   std::string("IQIBulkH1Error") );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::H1_ERROR_ANALYTIC ) );
        tParameterList( 4 )( tIQICounter ).set( "IQI_output_type",            static_cast< uint >( vis::Output_Type::H1_ERROR_ANALYTIC ) );
        tParameterList( 4 )( tIQICounter ).set( "master_dof_dependencies",    std::string("TEMP") );
        tParameterList( 4 )( tIQICounter ).set( "master_properties",          std::string("PropExactTemperatureGradient,H1Check") );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names",             tTotalDomain );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name",                   std::string("IQIVolume") );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::VOLUME ) );
        tParameterList( 4 )( tIQICounter ).set( "IQI_output_type",            static_cast< uint >( vis::Output_Type::VOLUME ) );
        tParameterList( 4 )( tIQICounter ).set( "master_dof_dependencies",    std::string("TEMP") );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names",             tTotalDomain );
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
        tParameterlist( 2 )( 0 ).set("NLA_rel_res_norm_drop",    tNLA_rel_res_norm_drop );
        tParameterlist( 2 )( 0 ).set("NLA_relaxation_parameter", tNLA_relaxation_parameter  );
        tParameterlist( 2 )( 0 ).set("NLA_max_iter",             tNLA_max_iter );

        tParameterlist( 3 )( 0 ) = moris::prm::create_nonlinear_solver_parameter_list();
        tParameterlist( 3 )( 0 ).set("NLA_DofTypes"      , std::string("TEMP") );

        tParameterlist( 4 )( 0 ) = moris::prm::create_time_solver_algorithm_parameter_list();
        tParameterlist( 4 )( 0 ).set("TSA_Num_Time_Steps", tTSA_Num_Time_Steps );
        tParameterlist( 4 )( 0 ).set("TSA_Time_Frame",     tTSA_Time_Frame );

        tParameterlist( 5 )( 0 ) = moris::prm::create_time_solver_parameter_list();
        tParameterlist( 5 )( 0 ).set("TSA_DofTypes",           std::string("TEMP") );
        tParameterlist( 5 )( 0 ).set("TSA_Initialize_Sol_Vec", std::string("TEMP,0.0") );
        tParameterlist( 5 )( 0 ).set("TSA_Output_Indices",     std::string("0") );
        tParameterlist( 5 )( 0 ).set("TSA_Output_Crteria",     std::string("Output_Criterion") );

        tParameterlist( 6 )( 0 ) = moris::prm::create_solver_warehouse_parameterlist();
    }

    void MSIParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_msi_parameter_list();
        tParameterlist( 0 )( 0 ).set("order_adofs_by_host",false);
    }

    void VISParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_vis_parameter_list();
        tParameterlist( 0 )( 0 ).set( "File_Name"  , std::pair< std::string, std::string >( "./", tOutputFileName ) );
        tParameterlist( 0 )( 0 ).set( "Mesh_Type"  , static_cast< uint >( vis::VIS_Mesh_Type::STANDARD ) );
        tParameterlist( 0 )( 0 ).set( "Set_Names"  , tTotalDomain );
        tParameterlist( 0 )( 0 ).set( "Field_Names", std::string("TEMP,TEMP_ANALYTIC,L2_ERROR_ANALYTIC,H1_ERROR_ANALYTIC,VOLUME") );
        tParameterlist( 0 )( 0 ).set( "Field_Type" , std::string("NODAL,NODAL,GLOBAL,GLOBAL,GLOBAL") );
        tParameterlist( 0 )( 0 ).set( "Output_Type", std::string("TEMP,PROPERTY,L2_ERROR_ANALYTIC,H1_ERROR_ANALYTIC,VOLUME") );
        tParameterlist( 0 )( 0 ).set( "Save_Frequency", 1 );
    }

    /* ------------------------------------------------------------------------ */
}

//------------------------------------------------------------------------------
#ifdef  __cplusplus
}
#endif

