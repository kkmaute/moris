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

extern uint gInterpolationOrder;

#ifdef  __cplusplus
extern "C"
{
#endif
//------------------------------------------------------------------------------
namespace moris
{
    bool isGhost = false;
    
    std::string tBulk = "HMR_dummy_n_p0,HMR_dummy_n_p1,HMR_dummy_n_p2,HMR_dummy_n_p4,HMR_dummy_n_p8,HMR_dummy_c_p0,HMR_dummy_c_p1,HMR_dummy_c_p2,HMR_dummy_c_p4,HMR_dummy_c_p8";
    std::string tRightBC = "SideSet_2_c_p0,SideSet_2_c_p4,SideSet_2_c_p8,SideSet_2_n_p0,SideSet_2_n_p4,SideSet_2_n_p8";
    std::string tLeftBC = "SideSet_4_c_p1,SideSet_4_c_p2,SideSet_4_n_p1,SideSet_4_n_p2";
    std::string tInterface = "dbl_iside_p0_0_p1_1,dbl_iside_p0_0_p1_2,dbl_iside_p0_0_p1_8,dbl_iside_p0_0_p1_4,dbl_iside_p0_1_p1_2,dbl_iside_p0_1_p1_8";

    // Constant function for properties
    void Func_Const(
            moris::Matrix< moris::DDRMat >                 & aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
            moris::fem::Field_Interpolator_Manager         * aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    bool Output_Criterion( moris::tsa::Time_Solver * aTimeSolver )
    {
        return true;
    }

    moris::real Colors_1(
                const moris::Matrix< DDRMat >     & aCoordinates,
                const moris::Cell< moris::real* > & aGeometryParameters)
    {
        moris::real tVal = -1.0;
        if( aCoordinates(0) > 0.1117 && aCoordinates(1) <= -0.11)
        {
            tVal =1.0;
        }

        return tVal;
    }

    moris::real Colors_2(
                const moris::Matrix< DDRMat >     & aCoordinates,
                const moris::Cell< moris::real* > & aGeometryParameters)
    {
        moris::real tVal = -1.0;
        if( std::sqrt( std::pow( aCoordinates(0) -2.0, 2) + std::pow( aCoordinates(1) -2.0, 2) ) - 1.2 <= 0.0 &&
            !(aCoordinates(0) > 0.1117 && aCoordinates(1) <= -0.11 ) )
        {
            tVal =1.0;
        }

        return tVal;
    }

       moris::real Colors_3(
                const moris::Matrix< DDRMat >     & aCoordinates,
                const moris::Cell< moris::real* > & aGeometryParameters)
    {
        moris::real tVal = -1.0;
        if(  aCoordinates(0)-aCoordinates(1) +  1.5 <= 0.0 &&
            !( std::sqrt( std::pow( aCoordinates(0) -2.0, 2) + std::pow( aCoordinates(1) -2.0, 2) ) - 1.2 <= 0.0) &&
            !(aCoordinates(0) > 0.1117 && aCoordinates(1) <= -0.11 ))
        {
            tVal =1.0;
        }

        return tVal;
    }

       moris::real Colors_4(
                const moris::Matrix< DDRMat >     & aCoordinates,
                const moris::Cell< moris::real* > & aGeometryParameters)
    {
        moris::real tVal = -1.0;
        if( std::sqrt( std::pow( aCoordinates(0) +2.0, 2) + std::pow( aCoordinates(1) +2.0, 2) ) - 2.6 <= 0.0 &&
            !( aCoordinates(0)-aCoordinates(1) +  1.5 <= 0.0 ) &&
            !( std::sqrt( std::pow( aCoordinates(0) -2.0, 2) + std::pow( aCoordinates(1) -2.0, 2) ) - 1.2 <= 0.0) &&
            !( aCoordinates(0) > 0.1117 && aCoordinates(1) <= -0.11 ))
        {
            tVal =1.0;
        }

        return tVal;
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

        tParameterlist( 0 )( 0 ).set( "number_of_elements_per_dimension", "11,11");
        tParameterlist( 0 )( 0 ).set( "domain_dimensions",                "4,4");
        tParameterlist( 0 )( 0 ).set( "domain_offset",                    "-2.0,-2.0");
        tParameterlist( 0 )( 0 ).set( "domain_sidesets",                  "1,2,3,4");
        tParameterlist( 0 )( 0 ).set( "lagrange_output_meshes",           "0");

        tParameterlist( 0 )( 0 ).set( "lagrange_orders",  "1");
        tParameterlist( 0 )( 0 ).set( "lagrange_pattern", "0" );
        tParameterlist( 0 )( 0 ).set( "bspline_orders",   "1");
        tParameterlist( 0 )( 0 ).set( "bspline_pattern",  "0" );

        tParameterlist( 0 )( 0 ).set( "lagrange_to_bspline", "0") ;

        tParameterlist( 0 )( 0 ).set( "truncate_bsplines",  1 );
        tParameterlist( 0 )( 0 ).set( "refinement_buffer",  1 );
        tParameterlist( 0 )( 0 ).set( "staircase_buffer",   1 );
        tParameterlist( 0 )( 0 ).set( "initial_refinement", "0" );
        tParameterlist( 0 )( 0 ).set( "initial_refinement_pattern", "0" );

        tParameterlist( 0 )( 0 ).set( "use_number_aura",    1 );

        tParameterlist( 0 )( 0 ).set( "use_multigrid",  0 );
        tParameterlist( 0 )( 0 ).set( "severity_level", 0 );
    }

    void XTKParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_xtk_parameter_list();
        tParameterlist( 0 )( 0 ).set( "decompose",                   true );
        tParameterlist( 0 )( 0 ).set( "decomposition_type",          "conformal") ;
        tParameterlist( 0 )( 0 ).set( "enrich",                      true );
        tParameterlist( 0 )( 0 ).set( "basis_rank",                  "bspline") ;
        tParameterlist( 0 )( 0 ).set( "enrich_mesh_indices",         "0") ;
        tParameterlist( 0 )( 0 ).set( "ghost_stab",                  isGhost );
        tParameterlist( 0 )( 0 ).set( "multigrid",                   false );
        tParameterlist( 0 )( 0 ).set( "verbose",                     true );
        tParameterlist( 0 )( 0 ).set( "print_enriched_ig_mesh",      true );
        tParameterlist( 0 )( 0 ).set( "exodus_output_XTK_ig_mesh",   true );
    }

    void GENParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 3 );
        tParameterlist( 0 ).resize( 1 );
        tParameterlist( 1 ).resize( 4 );

        tParameterlist( 0 )( 0 ) = prm::create_gen_parameter_list();
        tParameterlist( 0 )( 0 ).set( "intersection_mode", "LEVEL_SET");
        tParameterlist( 0 )( 0 ).set( "isocontour_tolerance", 1E-12);

        // Geometry parameter lists
        tParameterlist( 1 )( 0 ) = prm::create_user_defined_geometry_parameter_list();
        tParameterlist( 1 )( 0 ).set( "field_function_name", "Colors_1");
        tParameterlist( 1 )( 0 ).set( "number_of_refinements", "2");
        tParameterlist( 1 )( 0 ).set( "refinement_mesh_index", "0");
        tParameterlist( 1 )( 0 ).set( "discretization_mesh_index", -1);
        tParameterlist( 1 )( 0 ).set( "multilinear_intersections", false);

                // Geometry parameter lists
        tParameterlist( 1 )( 1 ) = prm::create_user_defined_geometry_parameter_list();
        tParameterlist( 1 )( 1 ).set( "field_function_name", "Colors_2");
        tParameterlist( 1 )( 1 ).set( "number_of_refinements", "2");
        tParameterlist( 1 )( 1 ).set( "refinement_mesh_index", "0");
        tParameterlist( 1 )( 1 ).set( "discretization_mesh_index", -1);
        tParameterlist( 1 )( 1 ).set( "multilinear_intersections", false);

                // Geometry parameter lists
        tParameterlist( 1 )( 2 ) = prm::create_user_defined_geometry_parameter_list();
        tParameterlist( 1 )( 2 ).set( "field_function_name", "Colors_3");
        tParameterlist( 1 )( 2 ).set( "number_of_refinements", "2");
        tParameterlist( 1 )( 2 ).set( "refinement_mesh_index", "0");
        tParameterlist( 1 )( 2 ).set( "discretization_mesh_index", -1);
        tParameterlist( 1 )( 2 ).set( "multilinear_intersections", false);

                // Geometry parameter lists
        tParameterlist( 1 )( 3 ) = prm::create_user_defined_geometry_parameter_list();
        tParameterlist( 1 )( 3 ).set( "field_function_name", "Colors_4");
        tParameterlist( 1 )( 3 ).set( "number_of_refinements", "2");
        tParameterlist( 1 )( 3 ).set( "refinement_mesh_index", "0");
        tParameterlist( 1 )( 3 ).set( "discretization_mesh_index", -1);
        tParameterlist( 1 )( 3 ).set( "multilinear_intersections", false);


    }

    void FEMParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterList )
    {
        // create a cell of cell of parameter list for fem
        tParameterList.resize( 8 );

        //------------------------------------------------------------------------------
        // fill the property part of the parameter list

        // init property counter
        uint tPropCounter = 0;

        // create parameter list for property 1
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropDensity") ;
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      "1.0") ;
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const") ;
        tPropCounter++;

        // create parameter list for property 2
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropCapacity") ;
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      "1.0") ;
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const") ;
        tPropCounter++;

        // create parameter list for property 3
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropConductivity") ;
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      "0.00005") ;
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const") ;
        tPropCounter++;

        // create parameter list for property 5
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropInletTemp") ;
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      "2.0") ;
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const") ;
        tPropCounter++;

        // create parameter list for property 6
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropSideFlux") ;
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      "2.0") ;
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const") ;
        tPropCounter++;

        //------------------------------------------------------------------------------
        // fill the constitutive model part of the parameter list

        // init CM counter
        uint tCMCounter = 0;

        // create parameter list for constitutive model 2
        tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_name", "CMDiffusion") ;
        tParameterList( 1 )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::DIFF_LIN_ISO ) );
        tParameterList( 1 )( tCMCounter ).set( "dof_dependencies",  std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( 1 )( tCMCounter ).set( "properties",
                "PropConductivity,Conductivity;"
                "PropDensity     ,Density;"
                "PropCapacity    ,HeatCapacity") ;
        tCMCounter++;

                // create parameter list for constitutive model 2
        tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_name", "CMDiffusion_2") ;
        tParameterList( 1 )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::DIFF_LIN_ISO ) );
        tParameterList( 1 )( tCMCounter ).set( "dof_dependencies",  std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( 1 )( tCMCounter ).set( "properties",
                "PropConductivity,Conductivity;"
                "PropDensity     ,Density;"
                "PropCapacity    ,HeatCapacity") ;
        tCMCounter++;

        //------------------------------------------------------------------------------
        // fill the stabilization parameter part of the parameter list

        // init SP counter
        uint tSPCounter = 0;


        // create parameter list for stabilization parameter 2
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name",      "SPDirichletNitscheT") ;
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type",      static_cast< uint >( fem::Stabilization_Type::DIRICHLET_NITSCHE ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters",     "100.0") ;
        tParameterList( 2 )( tSPCounter ).set( "master_properties",       "PropConductivity,Material") ;
        tSPCounter++;

        // create parameter list for stabilization parameter 8
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name",      "SPGPTemp") ;
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type",      static_cast< uint >( fem::Stabilization_Type::GHOST_DISPL ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters",     "0.005") ;
        tParameterList( 2 )( tSPCounter ).set( "master_properties",       "PropConductivity,Material") ;
        tSPCounter++;

        // Temperature - Shell - PCM
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name",  std::string("SPInterfaceNitsche") );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type",  static_cast< uint >( fem::Stabilization_Type::NITSCHE_INTERFACE ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", std::string("100.0") );
        tParameterList( 2 )( tSPCounter ).set( "master_properties",   std::string("PropConductivity,Material") );
        tParameterList( 2 )( tSPCounter ).set( "slave_properties",    std::string("PropConductivity,Material") );
        tSPCounter++;

        //------------------------------------------------------------------------------
        // fill the IWG part of the parameter list

        // init IWG counter
        uint tIWGCounter = 0;

         // create parameter list for IWG 3
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWGDiffusionBulk") ;
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::SPATIALDIFF_BULK ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               "TEMP") ;
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    "TEMP") ;
        tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", "CMDiffusion,Diffusion") ;
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             tBulk) ;
        tIWGCounter++;

        // create parameter list for IWG 11
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWGInletTemp") ;
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::SPATIALDIFF_DIRICHLET_SYMMETRIC_NITSCHE ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               "TEMP") ;
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    "TEMP") ;
        tParameterList( 3 )( tIWGCounter ).set( "master_properties",          "PropInletTemp,Dirichlet") ;
        tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", "CMDiffusion,Diffusion") ;
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   "SPDirichletNitscheT,DirichletNitsche") ;
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             tRightBC) ;
        tIWGCounter++;


        // create parameter list for IWG 11
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWGCylinderFluxTemp") ;
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::SPATIALDIFF_NEUMANN ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               "TEMP") ;
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    "TEMP") ;
        tParameterList( 3 )( tIWGCounter ).set( "master_properties",          "PropSideFlux,Neumann") ;
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             tLeftBC) ;
        tIWGCounter++;


        // Temperature - Shell - Shell
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   std::string("IWGInterfaceShellShellTEMP") );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               std::string("TEMP") );
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    std::string("TEMP") );
        tParameterList( 3 )( tIWGCounter ).set( "slave_dof_dependencies",     std::string("TEMP") );
        tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", std::string("CMDiffusion,Diffusion") );
        tParameterList( 3 )( tIWGCounter ).set( "slave_constitutive_models",  std::string("CMDiffusion_2,Diffusion") );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",
                std::string("SPInterfaceNitsche     ,NitscheInterface") );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             tInterface );
        tIWGCounter++;



        if (isGhost)
        {
            // create parameter list for IWG 16
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWGGPTemp") ;
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::GHOST_NORMAL_FIELD ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               "TEMP") ;
            tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    "TEMP") ;
            tParameterList( 3 )( tIWGCounter ).set( "slave_dof_dependencies",     "TEMP") ;
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   "SPGPTemp,GhostSP") ;
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             "ghost_p160") ;
            tIWGCounter++;
        }

        //------------------------------------------------------------------------------
        // fill the IQI part of the parameter list

        // init IQI counter
        uint tIQICounter = 0;

        // create parameter list for IQI 3
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name",                   "IQIBulkTEMP") ;
        tParameterList( 4 )( tIQICounter ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::DOF ) );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity",               "TEMP");
        tParameterList( 4 )( tIQICounter ).set( "master_dof_dependencies",    "TEMP") ;
        //tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index",      0 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names",             tBulk) ;
        tIQICounter++;

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        tParameterList( 5 ).resize( 1 );
        tParameterList( 5 )( 0 ) = prm::create_computation_parameter_list();
    }

    void SOLParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 7 );
        for( uint Ik = 0; Ik < 7; Ik ++ )
        {
            tParameterlist( Ik ).resize( 1 );
        }

        tParameterlist( 0 )( 0 ) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::BELOS_IMPL );
        tParameterlist( 0 )( 0 ).set( "ifpack_prec_type", "ILU");

        tParameterlist( 1 )( 0 ) = moris::prm::create_linear_solver_parameter_list();

        tParameterlist( 2 )( 0 ) = moris::prm::create_nonlinear_algorithm_parameter_list();
        tParameterlist( 2 )( 0 ).set("NLA_rel_res_norm_drop",    1e-04 );
        tParameterlist( 2 )( 0 ).set("NLA_relaxation_parameter", 1.0  );
        tParameterlist( 2 )( 0 ).set("NLA_max_iter", 20 );

        tParameterlist( 3 )( 0 ) = moris::prm::create_nonlinear_solver_parameter_list();
        tParameterlist( 3 )( 0 ).set("NLA_DofTypes"      , "TEMP") ;

        tParameterlist( 4 )( 0 ) = moris::prm::create_time_solver_algorithm_parameter_list();

        tParameterlist( 5 )( 0 ) = moris::prm::create_time_solver_parameter_list();
        tParameterlist( 5 )( 0 ).set("TSA_DofTypes"       , "TEMP") ;
        tParameterlist( 5 )( 0 ).set("TSA_Initialize_Sol_Vec" , "TEMP,0.0") ;
        tParameterlist( 5 )( 0 ).set("TSA_Output_Indices" , "0") ;
        tParameterlist( 5 )( 0 ).set("TSA_Output_Crteria" , "Output_Criterion") ;

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
        tParameterlist( 0 )( 0 ).set( "File_Name"     , std::pair< std::string, std::string >( "./", "Grain_Structure.exo" ) );
        tParameterlist( 0 )( 0 ).set( "Mesh_Type"     , static_cast< uint >( vis::VIS_Mesh_Type::STANDARD ) );
        tParameterlist( 0 )( 0 ).set( "Set_Names"     ,  tBulk ) ;
        tParameterlist( 0 )( 0 ).set( "Field_Names"   ,  "TEMP" ) ;
        tParameterlist( 0 )( 0 ).set( "Field_Type"    ,  "NODAL" ) ;
        tParameterlist( 0 )( 0 ).set( "IQI_Names"     ,  "IQIBulkTEMP" ) ;
        tParameterlist( 0 )( 0 ).set( "Save_Frequency", 1 );
    }

    void MORISGENERALParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {

    }

    //------------------------------------------------------------------------------
}

//------------------------------------------------------------------------------
#ifdef  __cplusplus
}
#endif
