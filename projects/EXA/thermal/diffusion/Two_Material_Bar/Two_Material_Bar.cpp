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
#include "cl_HMR_Element.hpp"
#include "fn_equal_to.hpp"

#include "AztecOO.h"

//---------------------------------------------------------------

// global variable for interpolation order
extern uint gInterpolationOrder;

// problem dimension: 2D or 3D
extern uint gDim;

// problem dimension: 2D or 3D
extern uint gTestCaseIndex;

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

    // Phase 1: back  - Material 1
    // Phase 2: front - Material 2

    std::string tPhase1        = "HMR_dummy_n_p1,HMR_dummy_c_p1";
    std::string tPhase2        = "HMR_dummy_n_p0,HMR_dummy_c_p0";

    std::string tInterface     = "dbl_iside_p0_1_p1_0";

    std::string tBackSurface   = "SideSet_4_n_p1";

    std::string tFrontSurface  = "SideSet_2_n_p0";

    std::string tPhase1Ghost   = "ghost_p1";
    std::string tPhase2Ghost   = "ghost_p0";

    std::string tTotalDomain   = tPhase1 + "," + tPhase2;

    /* ------------------------------------------------------------------------ */    
    // geometry parameters

    // dimensionality 2D or 3D
    uint gDim = 3;

    // general
    moris::real sL  = 10.0;     // total length
    moris::real sL1 =  4.1;
    moris::real sL2 = sL - sL1;

    /* ------------------------------------------------------------------------ */
    // boundary condition

    // prescribed temperature
    moris::real sTpre = 0.0;

    // flux at tip fac
    moris::real sP2 = 5.0;

    /* ------------------------------------------------------------------------ */    
    // material parameters

    // capacity
    std::string sCap1 = "0.0"; 
    std::string sCap2 = "0.0";

    // density
    std::string tDens1 = "0.0"; 
    std::string tDens2 = "0.0";

    // conductivity
    moris::real sK1 = 1.0; 
    moris::real sK2 = 0.125;

    // body flux
    moris::real  sQ1 = 1.0;
    moris::real  sQ2 = 2.0;

    /* ------------------------------------------------------------------------ */
    // HMR parameters

    std::string tNumElemsPerDim     = gDim == 2 ? "2,   1"                      : "2,   1,   1";
    std::string tDomainDims         = gDim == 2 ?  std::to_string(sL) + ", 1.0" : std::to_string(sL) + ", 1.0, 1.0";
    std::string tDomainOffset       = gDim == 2 ? "0.0,  0.0"                   : "0.0,  0.0, 0.0";
    std::string tDomainSidesets     = gDim == 2 ? "1,2,3,4"                     : "1,2,3,4,5,6";

    std::string tInterpolationOrder = "1";

    int tRefineBuffer      = 1;

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
    bool tUseGhost = false;

    /* ------------------------------------------------------------------------ */
    // Output Config

    std::string tOutputFileName = "Two_Material_Bar_" + std::to_string(gTestCaseIndex) + ".exo";

    /* ------------------------------------------------------------------------ */
    // Level set function for diamond shaped wedge

    moris::real Inclusion(
            const moris::Matrix< DDRMat >     & aCoordinates,
            const moris::Cell< moris::real* > & aGeometryParameters )
    {        
        // distance from sphere center
        moris::real tX = aCoordinates(0);

        // Compute Signed-Distance field 
        moris::real tVal = sL1 - tX;

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

    /* ------------------------------------------------------------------------ */
    // Exact temperature
    // see Matlab file in test src folder

    void Func_Exact_Temperature(
            moris::Matrix< moris::DDRMat >                 & aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
            moris::fem::Field_Interpolator_Manager         * aFIManager )
    {
        // distance from sphere center
        moris::real tX = aFIManager->get_IP_geometry_interpolator()->valx()( 0 );

        // Update geometry parameters
        sL2 = sL - sL1;

        if (tX<sL1)
        {
            aPropMatrix = {{ sTpre - (sQ1*tX*tX)/(2*sK1) + (tX*(sP2 + sL1*sQ1 + sL2*sQ2))/sK1 }};
        }
        else
        {
            aPropMatrix = {{ (2*sL1*sK2*sP2 - 2*sL1*sK1*sP2 + 2*sTpre*sK1*sK2 - sL1*sL1*sK1*sQ2 + sL1*sL1*sK2*sQ1 - 2*sL1*sL2*sK1*sQ2 + 2*sL1*sL2*sK2*sQ2)/(2*sK1*sK2) - (sQ2*tX*tX)/(2*sK2) + (tX*(sP2 + sL1*sQ2 + sL2*sQ2))/sK2 }};
        } 
    }

    /* ------------------------------------------------------------------------ */
    // Exact temperature gradients 
    // see Matlab file in test src folder

    void Func_Exact_TemperatureGradient(
            moris::Matrix< moris::DDRMat >                 & aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
            moris::fem::Field_Interpolator_Manager         * aFIManager )
    {
        // distance from sphere center
        moris::real tX = aFIManager->get_IP_geometry_interpolator()->valx()( 0 );

        // Update geometry parameters
        sL2 = sL - sL1;

        // set size for aPropMatrix
        aPropMatrix.set_size( gDim, 1, 0.0 );

        // spatial gradients of analytic temperature distribution
        if ( tX <= sL1 )
        {
            aPropMatrix( 0, 0 ) = (sP2 + sL1*sQ1 + sL2*sQ2)/sK1 - (sQ1*tX)/sK1;
        }
        else 
        {
            aPropMatrix( 0, 0 ) = (sP2 + sL1*sQ2 + sL2*sQ2)/sK2 - (sQ2*tX)/sK2;
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
        tParameterlist( 0 )( 0 ).set( "lagrange_output_meshes",           "0");

        tParameterlist( 0 )( 0 ).set( "lagrange_orders",  std::to_string(gInterpolationOrder) );
        tParameterlist( 0 )( 0 ).set( "lagrange_pattern",  "0" )  ;
        tParameterlist( 0 )( 0 ).set( "bspline_orders",   std::to_string(gInterpolationOrder) );
        tParameterlist( 0 )( 0 ).set( "bspline_pattern",   "0" )  ;

        tParameterlist( 0 )( 0 ).set( "lagrange_to_bspline", "0") ;

        tParameterlist( 0 )( 0 ).set( "truncate_bsplines",  1 );
        tParameterlist( 0 )( 0 ).set( "refinement_buffer",  tRefineBuffer );
        tParameterlist( 0 )( 0 ).set( "staircase_buffer",   tRefineBuffer );
        tParameterlist( 0 )( 0 ).set( "initial_refinement", "1" );
        tParameterlist( 0 )( 0 ).set( "initial_refinement_pattern", "0" );

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
        tParameterlist( 0 )( 0 ).set( "ghost_stab",                tUseGhost );
        tParameterlist( 0 )( 0 ).set( "multigrid",                 false );
        tParameterlist( 0 )( 0 ).set( "verbose",                   true );
        tParameterlist( 0 )( 0 ).set( "print_enriched_ig_mesh",    false );
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
        tParameterlist( 1 )( tGeoCounter ).set( "field_function_name",       "Inclusion" );
    }

    void FEMParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterList )
    {
        // create a cell of cell of parameter list for fem
        tParameterList.resize( 8 );

        //------------------------------------------------------------------------------
        // init property counter
        uint tPropCounter = 0;

        // properties for inclusion

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ) = prm::create_property_parameter_list();
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropDensity1") ;
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      tDens1 );               
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const") ;
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropCapacity1") ;
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      sCap1 );                 
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const") ;
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropConductivity1") ;
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      std::to_string(sK1) );              
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const") ;
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropHeatLoad1") ;
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      std::to_string(sQ1) );              
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const") ;
        tPropCounter++;

        // properties for 2 material

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ) = prm::create_property_parameter_list();
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropDensity2") ;
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      tDens2 );               
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const") ;
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropCapacity2") ;
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      sCap2 );                 
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const") ;
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropConductivity2") ;
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      std::to_string(sK2) );              
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const") ;
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropHeatLoad2") ;
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      std::to_string(sQ2) );              
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const") ;
        tPropCounter++;

        // surface flux
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropSurfaceFlux") ;
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      std::to_string(sP2) );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const") ;
        tPropCounter++;

        // temperature at back surface
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropImposedTemperature") ;
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Exact_Temperature") ;
        tPropCounter++;

        // time continuity weights        
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropWeightCurrent") ;
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      "100.0") ;
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const") ;
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropWeightPrevious") ;
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      "100.0") ;
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const") ;
        tPropCounter++;

        // initial condition       
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropInitialCondition") ;
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Exact_Temperature") ;
        tPropCounter++;

        // exact temperature
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropExactTemperature") ;
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Exact_Temperature") ;
        tPropCounter++;

        // exact temperature gradient
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropExactTemperatureGradient") ;
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Exact_TemperatureGradient") ;
        tPropCounter++;

        //------------------------------------------------------------------------------
        // init CM counter
        uint tCMCounter = 0;

        // create parameter list for constitutive model - Inclusion
        tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_name", "CMDiffusion1") ;
        tParameterList( 1 )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::DIFF_LIN_ISO ) );
        tParameterList( 1 )( tCMCounter ).set( "dof_dependencies",  std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( 1 )( tCMCounter ).set( "properties",
                "PropConductivity1 , Conductivity;"
                "PropDensity1      , Density;"
                "PropCapacity1     , HeatCapacity") ;
        tCMCounter++;

        // create parameter list for constitutive model - 2 Material
        tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_name", "CMDiffusion2") ;
        tParameterList( 1 )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::DIFF_LIN_ISO ) );
        tParameterList( 1 )( tCMCounter ).set( "dof_dependencies",  std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( 1 )( tCMCounter ).set( "properties",
                "PropConductivity2 , Conductivity;"
                "PropDensity2      , Density;"
                "PropCapacity2     , HeatCapacity") ;
        tCMCounter++;

        //------------------------------------------------------------------------------
        // init SP counter
        uint tSPCounter = 0;

        // create parameter list for Nitsche stabilization parameter for inclusion-2 material interface
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name",  "SPInterfaceNitsche") ;
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type",  static_cast< uint >( fem::Stabilization_Type::NITSCHE_INTERFACE ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", "100.0") ;
        tParameterList( 2 )( tSPCounter ).set( "master_properties",   "PropConductivity1,Material") ;
        tParameterList( 2 )( tSPCounter ).set( "slave_properties",    "PropConductivity2,Material") ;
        tSPCounter++;

        // create parameter list for DBC on back surface
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name",      "SPNitscheTemp") ;
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type",      static_cast< uint >( fem::Stabilization_Type::DIRICHLET_NITSCHE ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters",     "100.0") ;
        tParameterList( 2 )( tSPCounter ).set( "master_properties",       "PropConductivity1,Material") ;
        tSPCounter++;

        if (tUseGhost)
        {
            // create parameter list for ghost stabilization parameter for inclusion
            tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
            tParameterList( 2 )( tSPCounter ).set( "stabilization_name",      "SPGPTemp1") ;
            tParameterList( 2 )( tSPCounter ).set( "stabilization_type",      static_cast< uint >( fem::Stabilization_Type::GHOST_DISPL ) );
            tParameterList( 2 )( tSPCounter ).set( "function_parameters",     "0.01") ;
            tParameterList( 2 )( tSPCounter ).set( "master_properties",       "PropConductivity1,Material") ;
            tSPCounter++;

            // create parameter list for ghost stabilization parameter for 2 material
            tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
            tParameterList( 2 )( tSPCounter ).set( "stabilization_name",      "SPGPTemp2") ;
            tParameterList( 2 )( tSPCounter ).set( "stabilization_type",      static_cast< uint >( fem::Stabilization_Type::GHOST_DISPL ) );
            tParameterList( 2 )( tSPCounter ).set( "function_parameters",     "0.01") ;
            tParameterList( 2 )( tSPCounter ).set( "master_properties",       "PropConductivity2,Material") ;
            tSPCounter++;
        }

        //------------------------------------------------------------------------------
        // init IWG counter
        uint tIWGCounter = 0;

        // create IWG for inclusion - bulk diffusion
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWGDiffusion1Bulk") ;
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::SPATIALDIFF_BULK ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               "TEMP") ;
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    "TEMP") ;
        tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", "CMDiffusion1,Diffusion") ;
        tParameterList( 3 )( tIWGCounter ).set( "master_properties",          "PropHeatLoad1,Load") ;
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             tPhase1 );
        tIWGCounter++;   

        // create IWG for 2 material - bulk diffusion
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWGDiffusion2Bulk") ;
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::SPATIALDIFF_BULK ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               "TEMP") ;
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    "TEMP") ;
        tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", "CMDiffusion2,Diffusion") ;
        tParameterList( 3 )( tIWGCounter ).set( "master_properties",          "PropHeatLoad2,Load") ;
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             tPhase2 );
        tIWGCounter++;   

        // create IWG for Neumann boundary conditions
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWGInletFlux") ;
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::SPATIALDIFF_NEUMANN ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               "TEMP") ;
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    "TEMP") ;
        tParameterList( 3 )( tIWGCounter ).set( "master_properties",          "PropSurfaceFlux,Neumann") ;
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             tFrontSurface );
        tIWGCounter++;

        // create IWG for Dirichlet boundary conditions
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWG2SurfaceTemp") ;
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               "TEMP") ;
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    "TEMP") ;
        tParameterList( 3 )( tIWGCounter ).set( "master_properties",          "PropImposedTemperature,Dirichlet") ;
        tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", "CMDiffusion1,Diffusion") ;
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   "SPNitscheTemp,DirichletNitsche") ;
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             tBackSurface );
        tIWGCounter++;

        // create parameter list for interface conditions
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWGInterface12TEMP") ;
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               "TEMP");
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    "TEMP");
        tParameterList( 3 )( tIWGCounter ).set( "slave_dof_dependencies",     "TEMP");
        tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", "CMDiffusion1,Diffusion");
        tParameterList( 3 )( tIWGCounter ).set( "slave_constitutive_models",  "CMDiffusion2,Diffusion");
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   "SPInterfaceNitsche,NitscheInterface");
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             tInterface );
        tIWGCounter++;

        if (tUseGhost)
        {
            // create IWG for 2 material - ghost
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWGGP1Temp") ;
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::GHOST_NORMAL_FIELD ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               "TEMP") ;
            tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    "TEMP") ;
            tParameterList( 3 )( tIWGCounter ).set( "slave_dof_dependencies",     "TEMP") ;
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   "SPGPTemp1,GhostSP") ;
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             tPhase1Ghost );
            tIWGCounter++;

            // create IWG for 2 material - ghost
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWGGP1Temp") ;
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::GHOST_NORMAL_FIELD ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               "TEMP") ;
            tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    "TEMP") ;
            tParameterList( 3 )( tIWGCounter ).set( "slave_dof_dependencies",     "TEMP") ;
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   "SPGPTemp2,GhostSP") ;
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             tPhase2Ghost );
            tIWGCounter++;
        }

        //        // create IWG for time continuity
        //        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        //        tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWGTimeContinuityTemp") ;
        //        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::TIME_CONTINUITY_DOF ) );
        //        tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               "TEMP") ;
        //        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    "TEMP") ;
        //        tParameterList( 3 )( tIWGCounter ).set( "master_properties",
        //                "PropWeightCurrent   ,WeightCurrent;"
        //                "PropWeightPrevious  ,WeightPrevious;"
        //                "PropInitialCondition,InitialCondition") ;
        //        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             tTotalDomain );
        //        tParameterList( 3 )( tIWGCounter ).set( "time_continuity",            true );
        //        tIWGCounter++;


        //------------------------------------------------------------------------------
        // init IQI counter
        uint tIQICounter = 0;

        // Nodal Temperature IQI
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name",                   "IQIBulkTEMP") ;
        tParameterList( 4 )( tIQICounter ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::DOF ) );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity",               "TEMP");
        tParameterList( 4 )( tIQICounter ).set( "master_dof_dependencies",    "TEMP") ;
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index",      0 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names",             tTotalDomain );
        tIQICounter++; 

        // create parameter list for IQI 1
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name",                   "IQIBulkTEMPAnalytic") ;
        tParameterList( 4 )( tIQICounter ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::PROPERTY ) );
        tParameterList( 4 )( tIQICounter ).set( "master_properties",          "PropExactTemperature,Property") ;
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names",             tTotalDomain );
        tIQICounter++;

        // create parameter list for IQI 2
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name",                   "IQIBulkL2Error") ;
        tParameterList( 4 )( tIQICounter ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::L2_ERROR_ANALYTIC ) );
        tParameterList( 4 )( tIQICounter ).set( "master_dof_dependencies",    "TEMP") ;
        tParameterList( 4 )( tIQICounter ).set( "master_properties",          "PropExactTemperature,L2Check") ;
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names",             tTotalDomain );
        tIQICounter++;

        // create parameter list for IQI 3
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name",                   "IQIBulkH1Error") ;
        tParameterList( 4 )( tIQICounter ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::H1_ERROR_ANALYTIC ) );
        tParameterList( 4 )( tIQICounter ).set( "master_dof_dependencies",    "TEMP") ;
        tParameterList( 4 )( tIQICounter ).set( "master_properties",          "PropExactTemperatureGradient,H1Check") ;
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names",             tTotalDomain );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name",                   "IQIVolume") ;
        tParameterList( 4 )( tIQICounter ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::VOLUME ) );
        tParameterList( 4 )( tIQICounter ).set( "master_dof_dependencies",    "TEMP") ;
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names",             tTotalDomain );
        tIQICounter++;

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        tParameterList( 5 ).resize( 1 );
        tParameterList( 5 )( 0 ) = prm::create_computation_parameter_list();
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
        tParameterlist( 3 )( 0 ).set("NLA_DofTypes"      , "TEMP") ;

        tParameterlist( 4 )( 0 ) = moris::prm::create_time_solver_algorithm_parameter_list();
        tParameterlist( 4 )( 0 ).set("TSA_Num_Time_Steps", tTSA_Num_Time_Steps );
        tParameterlist( 4 )( 0 ).set("TSA_Time_Frame",     tTSA_Time_Frame );

        tParameterlist( 5 )( 0 ) = moris::prm::create_time_solver_parameter_list();
        tParameterlist( 5 )( 0 ).set("TSA_DofTypes",           "TEMP") ;
        tParameterlist( 5 )( 0 ).set("TSA_Initialize_Sol_Vec", "TEMP,0.0") ;
        tParameterlist( 5 )( 0 ).set("TSA_Output_Indices",     "0") ;
        tParameterlist( 5 )( 0 ).set("TSA_Output_Crteria",     "Output_Criterion") ;

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
        tParameterlist( 0 )( 0 ).set( "Field_Names", "TEMP,TEMP_ANALYTIC,L2_ERROR_ANALYTIC,H1_ERROR_ANALYTIC,VOLUME") ;
        tParameterlist( 0 )( 0 ).set( "Field_Type" , "NODAL,NODAL,GLOBAL,GLOBAL,GLOBAL") ;
        tParameterlist( 0 )( 0 ).set( "IQI_Names"  , "IQIBulkTEMP,IQIBulkTEMPAnalytic,IQIBulkL2Error,IQIBulkH1Error,IQIVolume") ;
        tParameterlist( 0 )( 0 ).set( "Save_Frequency", 1 );
    }

    void MORISGENERALParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {

    }

    /* ------------------------------------------------------------------------ */
}

//------------------------------------------------------------------------------
#ifdef  __cplusplus
}
#endif

