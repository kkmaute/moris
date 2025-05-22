// standard library
#include <string>
#include <iostream>

// Type definitions
#include "moris_typedefs.hpp"
#include "parameters.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

// modules
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_MSI_Equation_Object.hpp"
#include "cl_TSA_Time_Solver.hpp"
#include "cl_DLA_Solver_Interface.hpp"
#include "cl_DLA_Linear_Solver_Aztec.hpp"

// parameter lists
#include "fn_PRM_OPT_Parameters.hpp"
#include "fn_PRM_HMR_Parameters.hpp"
#include "fn_PRM_XTK_Parameters.hpp"
#include "fn_PRM_GEN_Parameters.hpp"
#include "fn_PRM_FEM_Parameters.hpp"
#include "fn_PRM_MSI_Parameters.hpp"
#include "fn_PRM_SOL_Parameters.hpp"
#include "fn_PRM_VIS_Parameters.hpp"

// LinAlg Operators
#include "fn_equal_to.hpp"
#include "fn_dot.hpp"
#include "fn_norm.hpp"

// Parsing Tools
#include "cl_Ascii.hpp"
#include "fn_Parsing_Tools.hpp"
#include "fn_assert.hpp"
#include "fn_stringify_matrix.hpp"


#ifdef __cplusplus
extern "C" {

#endif
//------------------------------------------------------------------------------
namespace moris
{
    // define the discretization
    uint tNumRefinements = 0;
    bool tIsEnriched = true;
    
    // Ghost
    std::string tGhostPenalty = "0.01";
    bool tUseGhost = i-0000;
    
    // shift beams in X-Direction
    real tBeamThickness = 1.4;
    real tBeamDistance  = 2.0;
    real tBeamOffSet    = i-0001 + ( tBeamThickness / 2.0 );

    // lower boundary (y-coordinate) where the beams meet the base plate
    real tBeamLowerBound = 1.6;
    
    // names of the mesh sets used
    std::string tBulkSet      = "HMR_dummy_n_p0,HMR_dummy_c_p0";
    std::string tGhostSet     = "ghost_p0";
    std::string tNeumannSet   = "SideSet_3_n_p0,SideSet_3_c_p0";
    std::string tDirichletSet = "SideSet_1_n_p0,SideSet_1_c_p0";

    // polynomial of the background interpolation
    std::string tPolyOrder = "2";
    
    // intersection computation mode
    bool tUseMultiLinear = false;

    /* ------------------------------------------------------------------------ */
    // geometry parameters & LS functions

    moris::real Func_Phi_Gutter(
            const moris::Matrix< DDRMat >     & aCoordinates,
            const moris::Vector< moris::real* > & aGeometryParameters )
    {
        // get coordinates
        real tX = aCoordinates( 0 );

        // Level-set function
        return std::sin( M_PI * ( tX - tBeamOffSet ) );
    }
    
    real tXMin0 = tBeamOffSet - 1.0 * tBeamDistance - 0.5 * tBeamThickness;
    real tXMax0 = tBeamOffSet - 1.0 * tBeamDistance + 0.5 * tBeamThickness;
    real tXMin1 = tBeamOffSet + 0.0 * tBeamDistance - 0.5 * tBeamThickness;
    real tXMax1 = tBeamOffSet + 0.0 * tBeamDistance + 0.5 * tBeamThickness;
    real tXMin2 = tBeamOffSet + 1.0 * tBeamDistance - 0.5 * tBeamThickness;
    real tXMax2 = tBeamOffSet + 1.0 * tBeamDistance + 0.5 * tBeamThickness;
    real tXMin3 = tBeamOffSet + 2.0 * tBeamDistance - 0.5 * tBeamThickness;
    real tXMax3 = tBeamOffSet + 2.0 * tBeamDistance + 0.5 * tBeamThickness;
    real tXMin4 = tBeamOffSet + 3.0 * tBeamDistance - 0.5 * tBeamThickness;
    real tXMax4 = tBeamOffSet + 3.0 * tBeamDistance + 0.5 * tBeamThickness;

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
        aParameterLists.set( "number_of_elements_per_dimension", 8, 12 );
        aParameterLists.set( "domain_dimensions", 8.0, 12.0 );
        aParameterLists.set( "domain_offset", 0.0, 0.0 );

        aParameterLists.set( "lagrange_output_meshes", "0" );

        aParameterLists.set( "lagrange_orders", tPolyOrder );
        aParameterLists.set( "lagrange_pattern", "0" );

        aParameterLists.set( "bspline_orders", tPolyOrder );
        aParameterLists.set( "bspline_pattern", "0" );
        aParameterLists.set( "lagrange_to_bspline", "0" );
        
        aParameterLists.set( "pattern_initial_refinement", tNumRefinements );
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
        aParameterLists.set( "ghost_stab", tUseGhost );
        aParameterLists.set( "visualize_ghost", tUseGhost );
        aParameterLists.set( "print_enriched_ig_mesh", false );
        aParameterLists.set( "exodus_output_XTK_ig_mesh", true );
        aParameterLists.set( "high_to_low_dbl_side_sets", true );
        aParameterLists.set( "only_generate_xtk_temp", false );
        aParameterLists.set( "output_cut_ig_mesh", true );
        aParameterLists.set( "verbose", true );
        aParameterLists.set( "verbose_level", moris::uint( 1 ) );

        if( !tIsEnriched )
        {
            aParameterLists.set( "unenriched_mesh_indices", "0" );
        }
    }

    //------------------------------------------------------------------------------

    void
    GENParameterList( Module_Parameter_Lists& aParameterLists )
    {
        // 3-beam version phase table
        Matrix< DDUMat > tPhaseMap( 128, 1, 1 );
        tPhaseMap(   0 ) = 0;
        tPhaseMap(  32 ) = 0;
        tPhaseMap(  48 ) = 0;
        tPhaseMap(  56 ) = 0;
        tPhaseMap(  60 ) = 0;
        tPhaseMap(  62 ) = 0;
        tPhaseMap(  63 ) = 0;
        tPhaseMap(  96 ) = 0;
        tPhaseMap( 120 ) = 0;
        tPhaseMap( 126 ) = 0;

        // // old (4-beam version) phase table
        // Matrix< DDUMat > tPhaseMap( 512, 1, 1 );
        // tPhaseMap(   0 ) = 0;
        // tPhaseMap( 128 ) = 0;
        // tPhaseMap( 192 ) = 0;
        // tPhaseMap( 224 ) = 0;
        // tPhaseMap( 240 ) = 0;
        // tPhaseMap( 248 ) = 0;
        // tPhaseMap( 252 ) = 0;
        // tPhaseMap( 254 ) = 0;
        // tPhaseMap( 255 ) = 0;
        // tPhaseMap( 384 ) = 0;
        // tPhaseMap( 480 ) = 0;
        // tPhaseMap( 504 ) = 0;
        // tPhaseMap( 510 ) = 0;

        // // new (5-beam version) phase table
        // Matrix< DDUMat > tPhaseMap( 2048, 1, 1 );
        // tPhaseMap( 1536 ) = 0;
        // tPhaseMap(  512 ) = 0;
        // tPhaseMap(  768 ) = 0;
        // tPhaseMap(  896 ) = 0;
        // tPhaseMap( 1920 ) = 0;
        // tPhaseMap(  896 ) = 0;
        // tPhaseMap(  960 ) = 0;
        // tPhaseMap(  992 ) = 0;
        // tPhaseMap( 2016 ) = 0;
        // tPhaseMap( 1008 ) = 0;
        // tPhaseMap( 1016 ) = 0;
        // tPhaseMap( 2040 ) = 0;
        // tPhaseMap( 1020 ) = 0;
        // tPhaseMap( 1022 ) = 0;
        // tPhaseMap( 2046 ) = 0;
        // tPhaseMap(    0 ) = 0;
        // tPhaseMap( 1023 ) = 0;
        
        std::string tPhaseMapString = moris::ios::stringify( tPhaseMap );
        aParameterLists.set( "phase_table", tPhaseMapString );
        aParameterLists.set( "print_phase_table", false );


        // Bottom Plane
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
        aParameterLists.set( "center_x", 0.0 );
        aParameterLists.set( "center_y", tBeamLowerBound );
        aParameterLists.set( "normal_x", 0.0 );
        aParameterLists.set( "normal_y", 1.0 );
        aParameterLists.set( "use_multilinear_interpolation", tUseMultiLinear );
        
        // Beam 0
        // aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
        // aParameterLists.set( "center_x", tXMin0 );
        // aParameterLists.set( "center_y", 0.0 );
        // aParameterLists.set( "normal_x", 1.0 );
        // aParameterLists.set( "normal_y", 0.0 );
        // aParameterLists.set( "use_multilinear_interpolation", tUseMultiLinear );

        // aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
        // aParameterLists.set( "center_x", tXMax0 );
        // aParameterLists.set( "center_y", 0.0 );
        // aParameterLists.set( "normal_x", 1.0 );
        // aParameterLists.set( "normal_y", 0.0 );
        // aParameterLists.set( "use_multilinear_interpolation", tUseMultiLinear );
        
        // Beam 1
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
        aParameterLists.set( "center_x", tXMin1 );
        aParameterLists.set( "center_y", 0.0 );
        aParameterLists.set( "normal_x", 1.0 );
        aParameterLists.set( "normal_y", 0.0 );
        aParameterLists.set( "use_multilinear_interpolation", tUseMultiLinear );
        
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
        aParameterLists.set( "center_x", tXMax1 );
        aParameterLists.set( "center_y", 0.0 );
        aParameterLists.set( "normal_x", 1.0 );
        aParameterLists.set( "normal_y", 0.0 );
        aParameterLists.set( "use_multilinear_interpolation", tUseMultiLinear );
        
        // Beam 2
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
        aParameterLists.set( "center_x", tXMin2 );
        aParameterLists.set( "center_y", 0.0 );
        aParameterLists.set( "normal_x", 1.0 );
        aParameterLists.set( "normal_y", 0.0 );
        aParameterLists.set( "use_multilinear_interpolation", tUseMultiLinear );
        
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
        aParameterLists.set( "center_x", tXMax2 );
        aParameterLists.set( "center_y", 0.0 );
        aParameterLists.set( "normal_x", 1.0 );
        aParameterLists.set( "normal_y", 0.0 );
        aParameterLists.set( "use_multilinear_interpolation", tUseMultiLinear );
        
        // Beam 3
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
        aParameterLists.set( "center_x", tXMin3 );
        aParameterLists.set( "center_y", 0.0 );
        aParameterLists.set( "normal_x", 1.0 );
        aParameterLists.set( "normal_y", 0.0 );
        aParameterLists.set( "use_multilinear_interpolation", tUseMultiLinear );
        
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
        aParameterLists.set( "center_x", tXMax3 );
        aParameterLists.set( "center_y", 0.0 );
        aParameterLists.set( "normal_x", 1.0 );
        aParameterLists.set( "normal_y", 0.0 );
        aParameterLists.set( "use_multilinear_interpolation", tUseMultiLinear );
        
        // // Beam 4
        // aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
        // aParameterLists.set( "center_x", tXMin4 );
        // aParameterLists.set( "center_y", 0.0 );
        // aParameterLists.set( "normal_x", 1.0 );
        // aParameterLists.set( "normal_y", 0.0 );
        // aParameterLists.set( "use_multilinear_interpolation", tUseMultiLinear );
        // 
        // aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
        // aParameterLists.set( "center_x", tXMax4 );
        // aParameterLists.set( "center_y", 0.0 );
        // aParameterLists.set( "normal_x", 1.0 );
        // aParameterLists.set( "normal_y", 0.0 );
        // aParameterLists.set( "use_multilinear_interpolation", tUseMultiLinear );
    }

    //------------------------------------------------------------------------------

    void
    FEMParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.hack_for_legacy_fem();

        //------------------------------------------------------------------------------
        // Properties

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropLoad" );
        aParameterLists.set( "function_parameters", "1.0;0.0" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropUnit" );
        aParameterLists.set( "function_parameters", "1.0" );
        
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropYoungs" );
        aParameterLists.set( "function_parameters", "1000.0" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropPoisson" );
        aParameterLists.set( "function_parameters", "0.3" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDirichlet" );
        aParameterLists.set( "function_parameters", "0.0;0.0" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropNeumann" );
        aParameterLists.set( "function_parameters", "1.0;0.0" );

        //------------------------------------------------------------------------------
        // Constitutive Models (CMs)

        // create parameter list for constitutive model 1
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMStrucLinIso1" );
        aParameterLists.set( "constitutive_type", (uint) fem::Constitutive_Type::STRUC_LIN_ISO );
        aParameterLists.set( "model_type", (uint) fem::Model_Type::PLANE_STRESS );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
        aParameterLists.set( "properties", "PropYoungs,YoungsModulus; PropPoisson,PoissonRatio" );

        //------------------------------------------------------------------------------
        // Stabilization Parameters (SPs)

        // create parameter list for stabilization parameter 1
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPNitsche" );
        aParameterLists.set( "stabilization_type", (uint) fem::Stabilization_Type::DIRICHLET_NITSCHE );
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "leader_properties", "PropYoungs,Material" );
        
        // create parameter list for ghost stabilization parameter for outer material
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPGhost" );
        aParameterLists.set( "stabilization_type", (uint) fem::Stabilization_Type::GHOST_DISPL );
        aParameterLists.set( "function_parameters", "0.01" );
        aParameterLists.set( "leader_properties", "PropYoungs,Material" );

        //------------------------------------------------------------------------------
        // IWGs

        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWG_Bulk" );
        aParameterLists.set( "IWG_type", (uint) fem::IWG_Type::STRUC_LINEAR_BULK );
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        aParameterLists.set( "leader_properties", "PropUnit,Thickness;PropLoad,Load" );
        aParameterLists.set( "mesh_set_names", tBulkSet );

        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWG_Dirichlet" );
        aParameterLists.set( "IWG_type", (uint) fem::IWG_Type::STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_properties", "PropDirichlet,Dirichlet;PropUnit,Thickness" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        aParameterLists.set( "stabilization_parameters", "SPNitsche,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", tDirichletSet );

        // aParameterLists( FEM::IWG ).add_parameter_list();
        // aParameterLists.set( "IWG_name", "IWG_Traction" );
        // aParameterLists.set( "IWG_type", (uint) fem::IWG_Type::STRUC_LINEAR_NEUMANN );
        // aParameterLists.set( "dof_residual", "UX,UY" );
        // aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        // aParameterLists.set( "leader_properties", "PropNeumann,Traction;PropUnit,Thickness" );
        // aParameterLists.set( "mesh_set_names", tNeumannSet );
        
        if ( tUseGhost )
        {
            // create IWG for outer material - ghost
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "IWG_Ghost" );
            aParameterLists.set( "IWG_type", (uint) fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists.set( "dof_residual", "UX,UY" );
            aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
            aParameterLists.set( "follower_dof_dependencies", "UX,UY" );
            aParameterLists.set( "stabilization_parameters", "SPGhost,GhostSP" );
            aParameterLists.set( "mesh_set_names", tGhostSet );
        }

        //------------------------------------------------------------------------------
        // IQIs

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkDISPX" );
        aParameterLists.set( "IQI_type", (uint) fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "mesh_set_names", tBulkSet );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkDISPY" );
        aParameterLists.set( "IQI_type", (uint) fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "vectorial_field_index", 1 );
        aParameterLists.set( "mesh_set_names", tBulkSet );

        // nodal von-mises stresses for shell
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIVMStress" );
        aParameterLists.set( "IQI_type", (uint) fem::IQI_Type::VON_MISES_STRESS );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        aParameterLists.set( "mesh_set_names", tBulkSet );

        // nodal stresses for shell
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIXStress" );
        aParameterLists.set( "IQI_type", (uint) fem::IQI_Type::NORMAL_STRESS );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "mesh_set_names", tBulkSet );

        // nodal von-mises stresses for shell
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIYStress" );
        aParameterLists.set( "IQI_type", (uint) fem::IQI_Type::NORMAL_STRESS );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        aParameterLists.set( "vectorial_field_index", 1 );
        aParameterLists.set( "mesh_set_names", tBulkSet );

        // nodal von-mises stresses for shell
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIOOPStress" );
        aParameterLists.set( "IQI_type", (uint) fem::IQI_Type::NORMAL_STRESS );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        aParameterLists.set( "vectorial_field_index", 2 );
        aParameterLists.set( "mesh_set_names", tBulkSet );

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        aParameterLists( FEM::COMPUTATION );
    }

    //------------------------------------------------------------------------------

    void
    SOLParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::AMESOS_IMPL );

        aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list();

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "NLA_rel_res_norm_drop", 1e-09 );
        aParameterLists.set( "NLA_relaxation_parameter", 1.0 );
        aParameterLists.set( "NLA_max_iter", 5 );
        aParameterLists.set( "NLA_combined_res_jac_assembly", false );

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "NLA_DofTypes", "UX,UY" );

        aParameterLists( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list();

        aParameterLists( SOL::TIME_SOLVERS ).add_parameter_list();
        aParameterLists.set( "TSA_DofTypes", "UX,UY" );
        aParameterLists.set( "TSA_Output_Indices", "0" );

        aParameterLists( SOL::SOLVER_WAREHOUSE );
        aParameterLists.set( "SOL_save_operator_to_matlab", "SOE" );
        
        aParameterLists( SOL::PRECONDITIONERS ).add_parameter_list(  sol::PreconditionerType::NONE );
    }

    //------------------------------------------------------------------------------

    void
    MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
    }

    //------------------------------------------------------------------------------

    void
    VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "File_Name", std::pair< std::string, std::string >( "./", "Multi_Beam.exo" ) );
        aParameterLists.set( "Mesh_Type", (uint) vis::VIS_Mesh_Type::STANDARD );
        aParameterLists.set( "Set_Names", tBulkSet );
        aParameterLists.set( "Field_Names", "UX,UY,STRESSVM,STRESSX,STRESSY,STRESSOOP" );
        aParameterLists.set( "Field_Type", "NODAL,NODAL,NODAL,NODAL,NODAL,NODAL" );
        aParameterLists.set( "IQI_Names", "IQIBulkDISPX,IQIBulkDISPY,IQIVMStress,IQIXStress,IQIYStress,IQIOOPStress" );
        aParameterLists.set( "Save_Frequency", 1 );
    }

    //------------------------------------------------------------------------------

    void
    MORISGENERALParameterList( moris::Vector< moris::Vector< Parameter_List > >& tParameterList )
    {
    }

    //------------------------------------------------------------------------------
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
