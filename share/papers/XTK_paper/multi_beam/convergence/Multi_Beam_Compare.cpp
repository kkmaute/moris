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
    uint tNumRefinements = i-0001;
    uint tNumLagRefinements = i-0002;
    bool tIsEnriched = i-0000;
    bool tUseGhost = false;
    
    // shift beams in X-Direction
    real tBeamThickness = 1.0;
    real tBeamDistance  = 2.0;
    real tBeamOffSet    = 0.6 + ( tBeamThickness / 2.0 );
    real tRadius        = ( tBeamDistance - tBeamThickness ) / 2.0;

    // lower boundary (y-coordinate) where the beams meet the base plate
    real tBeamLowerBound = 1.6;
    
    // names of the mesh sets used
    std::string tBulkSet      = "HMR_dummy_n_p0,HMR_dummy_c_p0";
    std::string tGhostSet     = "ghost_p0";
    std::string tNeumannSet   = "SideSet_3_n_p0,SideSet_3_c_p0";
    std::string tDirichletSet = "SideSet_1_n_p0,SideSet_1_c_p0";
    
    // intersection computation mode
    bool tUseMultiLinear = false;

    /* ------------------------------------------------------------------------ */
    // geometry parameters for LS functions
    
    real tXMin0 = tBeamOffSet + 0.0 * tBeamDistance - 0.5 * tBeamThickness;
    real tXMax0 = tBeamOffSet + 0.0 * tBeamDistance + 0.5 * tBeamThickness;
    real tXMin1 = tBeamOffSet + 1.0 * tBeamDistance - 0.5 * tBeamThickness;
    real tXMax1 = tBeamOffSet + 1.0 * tBeamDistance + 0.5 * tBeamThickness;
    real tXMin2 = tBeamOffSet + 2.0 * tBeamDistance - 0.5 * tBeamThickness;
    real tXMax2 = tBeamOffSet + 2.0 * tBeamDistance + 0.5 * tBeamThickness;
    real tXMin3 = tBeamOffSet + 3.0 * tBeamDistance - 0.5 * tBeamThickness;
    real tXMax3 = tBeamOffSet + 3.0 * tBeamDistance + 0.5 * tBeamThickness;
    
    real tXc00 = tXMin0 - tRadius;
    real tXc01 = 0.5 * ( tXMax0 + tXMin1 );
    real tXc12 = 0.5 * ( tXMax1 + tXMin2 );
    real tXc23 = 0.5 * ( tXMax2 + tXMin3 );
    real tXc33 = tXMax3 + tRadius;

    /* ------------------------------------------------------------------------ */
    // Read reference solution field and write as functions

    void
    Func_Field_Ref_UX(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix = aFIManager->get_field_interpolators_for_type( mtk::Field_Type::FIELD_1 )->val();
    }
    
    void
    Func_Field_Ref_UY(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix = aFIManager->get_field_interpolators_for_type( mtk::Field_Type::FIELD_2 )->val();
    }

    void
    Func_Field_Ref_U(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        real tUX = aFIManager->get_field_interpolators_for_type( mtk::Field_Type::FIELD_1 )->val()( 0 );
        real tUY = aFIManager->get_field_interpolators_for_type( mtk::Field_Type::FIELD_2 )->val()( 0 );
        
        aPropMatrix = { { tUX }, { tUY } };
    }
    
    void
    Func_Field_Ref_Grad_UXX(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix = aFIManager->get_field_interpolators_for_type( mtk::Field_Type::FIELD_3 )->val();
    }
    
    void
    Func_Field_Ref_Grad_UXY(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix = aFIManager->get_field_interpolators_for_type( mtk::Field_Type::FIELD_4 )->val();
    }
    
    void
    Func_Field_Ref_Grad_UYX(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix = aFIManager->get_field_interpolators_for_type( mtk::Field_Type::FIELD_5 )->val();
    }
    
    void
    Func_Field_Ref_Grad_UYY(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix = aFIManager->get_field_interpolators_for_type( mtk::Field_Type::FIELD_6 )->val();
    }
    
    void
    Func_Field_Ref_Grad_U(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        real tdUxdx = aFIManager->get_field_interpolators_for_type( mtk::Field_Type::FIELD_3 )->val()( 0 );
        real tdUxdy = aFIManager->get_field_interpolators_for_type( mtk::Field_Type::FIELD_4 )->val()( 0 );
        real tdUydx = aFIManager->get_field_interpolators_for_type( mtk::Field_Type::FIELD_5 )->val()( 0 );
        real tdUydy = aFIManager->get_field_interpolators_for_type( mtk::Field_Type::FIELD_6 )->val()( 0 );
        
        aPropMatrix = { { tdUxdx, tdUxdy }, { tdUydx, tdUydy } };
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
        aParameterLists.set( "number_of_elements_per_dimension", 8, 12 );
        aParameterLists.set( "domain_dimensions", 8.0, 12.0 );
        aParameterLists.set( "domain_offset", 0.0, 0.0 );

        aParameterLists.set( "lagrange_output_meshes", "0" );

        aParameterLists.set( "lagrange_orders", "2" );
        aParameterLists.set( "lagrange_pattern", "1" );

        aParameterLists.set( "bspline_orders", "2" );
        aParameterLists.set( "bspline_pattern", "0" );
        aParameterLists.set( "lagrange_to_bspline", "0" );
        
        aParameterLists.set( "pattern_initial_refinement", tNumRefinements, tNumLagRefinements );
                                                   // pattern 0: B-spline // pattern 1: Lagrange
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
        // aParameterLists.set( "verbose", true );
        // aParameterLists.set( "verbose_level", moris::uint( 1 ) );

        if( !tIsEnriched )
        {
            aParameterLists.set( "unenriched_mesh_indices", "0" );
        }
    }

    //------------------------------------------------------------------------------

    uint get_phase_index( const Bitset< 14 >& aGeometrySigns )
    {        
        // inside any of the circles (geoms. 9-13) -> void (1)
        if ( !aGeometrySigns.test( 9 ) || 
             !aGeometrySigns.test( 10 )|| 
             !aGeometrySigns.test( 11 )|| 
             !aGeometrySigns.test( 12 )|| 
             !aGeometrySigns.test( 13 ) )
        {
            return 1;
        }
        
        // outside the circles but below the bottom plate -> solid (0)
        if ( !aGeometrySigns.test( 0 ) )
        {
            return 0;
        }
        
        // inside beam 0
        if ( aGeometrySigns.test( 1 ) && !aGeometrySigns.test( 2 ) )
        {
            return 0;
        }
        
        // inside beam 1
        if ( aGeometrySigns.test( 3 ) && !aGeometrySigns.test( 4 ) )
        {
            return 0;
        }
        
        // inside beam 2
        if ( aGeometrySigns.test( 5 ) && !aGeometrySigns.test( 6 ) )
        {
            return 0;
        }
        
        // inside beam 3
        if ( aGeometrySigns.test( 7 ) && !aGeometrySigns.test( 8 ) )
        {
            return 0;
        }

        // if it's not inside one of the beams it must be void (1)
        return 1;
    }

    //------------------------------------------------------------------------------

    void
    GENParameterList( Module_Parameter_Lists& aParameterLists )
    {
        // phase function
        aParameterLists.set( "number_of_phases", 2 );
        aParameterLists.set( "phase_function_name", "get_phase_index" );

        // Bottom Plane
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
        aParameterLists.set( "center_x", 0.0 );
        aParameterLists.set( "center_y", tBeamLowerBound );
        aParameterLists.set( "normal_x", 0.0 );
        aParameterLists.set( "normal_y", 1.0 );
        aParameterLists.set( "use_multilinear_interpolation", tUseMultiLinear );
        
        // Beam 0
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
        aParameterLists.set( "center_x", tXMin0 );
        aParameterLists.set( "center_y", 0.0 );
        aParameterLists.set( "normal_x", 1.0 );
        aParameterLists.set( "normal_y", 0.0 );
        aParameterLists.set( "use_multilinear_interpolation", tUseMultiLinear );

        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
        aParameterLists.set( "center_x", tXMax0 );
        aParameterLists.set( "center_y", 0.0 );
        aParameterLists.set( "normal_x", 1.0 );
        aParameterLists.set( "normal_y", 0.0 );
        aParameterLists.set( "use_multilinear_interpolation", tUseMultiLinear );
        
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
        aParameterLists.set( "use_multilinear_interpolation", tUseMultiLinear );;
        
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
        
        // Cicle cutouts
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::CIRCLE ) );
        aParameterLists.set( "center_x", tXc00 );
        aParameterLists.set( "center_y", tBeamLowerBound );
        aParameterLists.set( "radius", tRadius );
        aParameterLists.set( "use_multilinear_interpolation", tUseMultiLinear );
        
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::CIRCLE ) );
        aParameterLists.set( "center_x", tXc01 );
        aParameterLists.set( "center_y", tBeamLowerBound );
        aParameterLists.set( "radius", tRadius );
        aParameterLists.set( "use_multilinear_interpolation", tUseMultiLinear );
        
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::CIRCLE ) );
        aParameterLists.set( "center_x", tXc12 );
        aParameterLists.set( "center_y", tBeamLowerBound );
        aParameterLists.set( "radius", tRadius );
        aParameterLists.set( "use_multilinear_interpolation", tUseMultiLinear );
        
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::CIRCLE ) );
        aParameterLists.set( "center_x", tXc23 );
        aParameterLists.set( "center_y", tBeamLowerBound );
        aParameterLists.set( "radius", tRadius );
        aParameterLists.set( "use_multilinear_interpolation", tUseMultiLinear );
        
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::CIRCLE ) );
        aParameterLists.set( "center_x", tXc33 );
        aParameterLists.set( "center_y", tBeamLowerBound );
        aParameterLists.set( "radius", tRadius );
        aParameterLists.set( "use_multilinear_interpolation", tUseMultiLinear );
    }

    //------------------------------------------------------------------------------

    void
    FEMParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.hack_for_legacy_fem();

        //------------------------------------------------------------------------------
        // Properties

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropLoad" );
        aParameterLists.set( "function_parameters", "1.0;0.0;0.0" );

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
        
        // Field properties for reference displacement
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropFieldRefUX" );
        aParameterLists.set( "value_function", "Func_Field_Ref_UX" );
        aParameterLists.set( "field_dependencies", "FIELD_1" );
        
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropFieldRefUY" );
        aParameterLists.set( "value_function", "Func_Field_Ref_UY" );
        aParameterLists.set( "field_dependencies", "FIELD_2" );
        
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropFieldRefU" );
        aParameterLists.set( "value_function", "Func_Field_Ref_U" );
        aParameterLists.set( "field_dependencies", "FIELD_1;FIELD_2" );
        
        // Field properties for reference displacement gradient
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropFieldRefGradUXX" );
        aParameterLists.set( "value_function", "Func_Field_Ref_Grad_UXX" );
        aParameterLists.set( "field_dependencies", "FIELD_3" );
        
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropFieldRefGradUXX" );
        aParameterLists.set( "value_function", "Func_Field_Ref_Grad_UXY" );
        aParameterLists.set( "field_dependencies", "FIELD_4" );
        
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropFieldRefGradUXX" );
        aParameterLists.set( "value_function", "Func_Field_Ref_Grad_UYX" );
        aParameterLists.set( "field_dependencies", "FIELD_5" );
        
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropFieldRefGradUXX" );
        aParameterLists.set( "value_function", "Func_Field_Ref_Grad_UYY" );
        aParameterLists.set( "field_dependencies", "FIELD_6" );
        
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropFieldRefGradU" );
        aParameterLists.set( "value_function", "Func_Field_Ref_Grad_U" );
        aParameterLists.set( "field_dependencies", "FIELD_3;FIELD_4;FIELD_5;FIELD_6" );

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

        /* aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWG_Traction" );
        aParameterLists.set( "IWG_type", (uint) fem::IWG_Type::STRUC_LINEAR_NEUMANN );
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_properties", "PropNeumann,Traction;PropUnit,Thickness" );
        aParameterLists.set( "mesh_set_names", tNeumannSet );
        */
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
        // init IQI counter
        uint tIQICounter = 0;

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQI_Bulk_Displ_X" );
        aParameterLists.set( "IQI_type", (uint) fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "mesh_set_names", tBulkSet );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQI_Bulk_Displ_Y" );
        aParameterLists.set( "IQI_type", (uint) fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "vectorial_field_index", 1 );
        aParameterLists.set( "mesh_set_names", tBulkSet );

        // nodal von-mises stresses for shell
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQI_VM_Stress" );
        aParameterLists.set( "IQI_type", (uint) fem::IQI_Type::VON_MISES_STRESS );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        aParameterLists.set( "mesh_set_names", tBulkSet );
        tIQICounter++;

        // nodal stresses for shell
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQI_X_Stress" );
        aParameterLists.set( "IQI_type", (uint) fem::IQI_Type::NORMAL_STRESS );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "mesh_set_names", tBulkSet );

        // nodal von-mises stresses for shell
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQI_Y_Stress" );
        aParameterLists.set( "IQI_type", (uint) fem::IQI_Type::NORMAL_STRESS );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        aParameterLists.set( "vectorial_field_index", 1 );
        aParameterLists.set( "mesh_set_names", tBulkSet );

        // nodal out-of-plane stresses for shell
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQI_OOP_Stress" );
        aParameterLists.set( "IQI_type", (uint) fem::IQI_Type::NORMAL_STRESS );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        aParameterLists.set( "vectorial_field_index", 2 );
        aParameterLists.set( "mesh_set_names", tBulkSet );
        
        // IQIs computing the L2 and H1s errors
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQI_L2_Error" );
        aParameterLists.set( "IQI_type", (uint) fem::IQI_Type::H1_ERROR );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "dof_quantity", "UX,UY" );
        aParameterLists.set( "leader_properties", "PropFieldRefU,L2_Reference" );
        aParameterLists.set( "function_parameters", "1.0/0.0" );
        aParameterLists.set( "mesh_set_names", tBulkSet );
        
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQI_H1s_Error" );
        aParameterLists.set( "IQI_type", (uint) fem::IQI_Type::H1_ERROR );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "dof_quantity", "UX,UY" );
        aParameterLists.set( "leader_properties", "PropFieldRefGradU,H1S_Reference" );
        aParameterLists.set( "function_parameters", "0.0/1.0" );
        aParameterLists.set( "mesh_set_names", tBulkSet );
        
        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        aParameterLists( FEM::COMPUTATION );
        
        //------------------------------------------------------------------------------
        // FIELDS
        
        // Input_fields
        aParameterLists( FEM::FIELDS ).add_parameter_list();
        aParameterLists.set( "field_name", "Field_Nodal_UX" );
        aParameterLists.set( "field_entity_type", "NODAL" );
        aParameterLists.set( "field_type", "FIELD_1" );
        aParameterLists.set( "field_create_from_file", "../../Reference_Solution/Field_Displ_X.hdf5" );
        
        aParameterLists( FEM::FIELDS ).add_parameter_list();
        aParameterLists.set( "field_name", "Field_Nodal_UY" );
        aParameterLists.set( "field_entity_type", "NODAL" );
        aParameterLists.set( "field_type", "FIELD_2" );
        aParameterLists.set( "field_create_from_file", "../../Reference_Solution/Field_Displ_Y.hdf5" );
        
        aParameterLists( FEM::FIELDS ).add_parameter_list();
        aParameterLists.set( "field_name", "Field_Nodal_Grad_UXX" );
        aParameterLists.set( "field_entity_type", "NODAL" );
        aParameterLists.set( "field_type", "FIELD_3" );
        aParameterLists.set( "field_create_from_file", "../../Reference_Solution/Field_Displ_Grad_XX.hdf5" );

        aParameterLists( FEM::FIELDS ).add_parameter_list();
        aParameterLists.set( "field_name", "Field_Nodal_Grad_UXY" );
        aParameterLists.set( "field_entity_type", "NODAL" );
        aParameterLists.set( "field_type", "FIELD_4" );
        aParameterLists.set( "field_create_from_file", "../../Reference_Solution/Field_Displ_Grad_XY.hdf5" );

        aParameterLists( FEM::FIELDS ).add_parameter_list();
        aParameterLists.set( "field_name", "Field_Nodal_Grad_UYX" );
        aParameterLists.set( "field_entity_type", "NODAL" );
        aParameterLists.set( "field_type", "FIELD_5" );
        aParameterLists.set( "field_create_from_file", "../../Reference_Solution/Field_Displ_Grad_YX.hdf5" );

        aParameterLists( FEM::FIELDS ).add_parameter_list();
        aParameterLists.set( "field_name", "Field_Nodal_Grad_UYY" );
        aParameterLists.set( "field_entity_type", "NODAL" );
        aParameterLists.set( "field_type", "FIELD_6" );
        aParameterLists.set( "field_create_from_file", "../../Reference_Solution/Field_Displ_Grad_YY.hdf5" );
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
        aParameterLists.set( "NLA_combined_res_jac_assembly", true );

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "NLA_DofTypes", "UX,UY" );

        aParameterLists( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list();

        aParameterLists( SOL::TIME_SOLVERS ).add_parameter_list();
        aParameterLists.set( "TSA_DofTypes", "UX,UY" );
        aParameterLists.set( "TSA_Output_Indices", "0" );

        aParameterLists( SOL::SOLVER_WAREHOUSE );
        // aParameterLists.set( "SOL_save_operator_to_matlab", "SOE" );
        // aParameterLists.set( "SOL_save_final_sol_vec_to_file", "Coarse_Solution_Vector.hdf5" );
        aParameterLists.set( "SOL_load_sol_vec_from_file", "../Coarse_Solution/Coarse_Solution_Vector.hdf5" );
        
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
        aParameterLists.set( "File_Name", std::pair< std::string, std::string >( "./", "Multi_Beam_Compare.exo" ) );
        aParameterLists.set( "Mesh_Type", (uint) vis::VIS_Mesh_Type::STANDARD );
        aParameterLists.set( "Set_Names", tBulkSet );
        aParameterLists.set( "Field_Names", "UX,UY,L2_Local,L2_Global,H1s_Local,H1s_Global" );
        aParameterLists.set( "Field_Type", "NODAL,NODAL,NODAL,GLOBAL,NODAL,GLOBAL" );
        aParameterLists.set( "IQI_Names", "IQI_Bulk_Displ_X,IQI_Bulk_Displ_Y,IQI_L2_Error,IQI_L2_Error,IQI_H1s_Error,IQI_H1s_Error" );
        aParameterLists.set( "Save_Frequency", 1 );
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
