/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Homogenization_3D.cpp
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
#include "fn_PRM_MIG_Parameters.hpp"
#include "cl_HMR_Element.hpp"
#include "fn_equal_to.hpp"

#include "AztecOO.h"

//---------------------------------------------------------------

#ifdef __cplusplus
extern "C" {
#endif
//------------------------------------------------------------------------------
namespace moris
{
    /* ------------------------------------------------------------------------ */
    // Mesh Set Information

    std::string tInnerPhase = "HMR_dummy_n_p1,HMR_dummy_c_p1";
    std::string tOuterPhase = "HMR_dummy_n_p0,HMR_dummy_c_p0";

    std::string tInterface = "dbl_iside_p0_1_p1_0";

    std::string tOuterSurface = "SideSet_1_n_p0,SideSet_2_n_p0,SideSet_3_n_p0,SideSet_4_n_p0,SideSet_5_n_p0,SideSet_6_n_p0";

    std::string tInnerPhaseGhost = "ghost_p1";
    std::string tOuterPhaseGhost = "ghost_p0";

    std::string tTotalDomain = tInnerPhase + "," + tOuterPhase;

    std::string tPeriodicSidePairs = "4,2;1,3;5,6";

    std::string tDBC = "SideSet_1_n_p0,SideSet_4_n_p0,SideSet_6_n_p0";

    real tDirichletLength = 0.1;
    /* ------------------------------------------------------------------------ */
    // geometry parameters

    moris::real tCenterX1 = 0.5;
    moris::real tCenterY1 = 0.5;
    moris::real tCenterZ1 = 0.5;
    moris::real tRadius   = 0.25;

    /* ------------------------------------------------------------------------ */
    // material parameters

    // Density, YoungsModulus, PoissonRatio
    std::string tDensityInner = "0.0";
    std::string tDensityOuter = "0.0";
    std::string tEmodInner    = "2.0";
    std::string tEmodOuter    = "1.0";
    std::string tPoisInner    = "0.1";
    std::string tPoisOuter    = "0.1";

    /* ------------------------------------------------------------------------ */
    // HMR parameters

    std::string tNumElemsPerDim = "5, 5 ,5 ";
    std::string tDomainDims     = "1.0, 1.0, 1.0";
    std::string tDomainOffset   = "0, 0, 0";
    std::string tDomainSidesets = "1,2,3,4, 5, 6";

    std::string tEigenStrain = "1.0;0.0;0.0;0.0;0.0;0.0";

    int tRefineBuffer = 1;

    /* ------------------------------------------------------------------------ */
    // Minimum level set value

    moris::real tMinLevs = 1.0e-8;

    /* ------------------------------------------------------------------------ */
    // Minimum level set value

    bool tUseGhost = true;

    /* ------------------------------------------------------------------------ */
    // Output Config

    std::string tOutputFileName = "Homogenization_3D.exo";

    /* ------------------------------------------------------------------------ */

    // Level set function for diamond shaped wedge
    moris::real
    Inclusion(
            const moris::Matrix< DDRMat >& aCoordinates,
            const Vector< real >&     aGeometryParameters )
    {
        // distance from sphere center
        moris::real tDx1 = aCoordinates( 0 ) - tCenterX1;
        moris::real tDy1 = aCoordinates( 1 ) - tCenterY1;
        moris::real tDz1 = aCoordinates( 2 ) - tCenterZ1;

        // Compute Signed-Distance field
        moris::real tVal = tRadius - std::sqrt( tDx1 * tDx1 + tDy1 * tDy1 + tDz1 * tDz1 );

        return ( std::abs( tVal ) < tMinLevs ) ? tMinLevs : tVal;
    }

    /* ------------------------------------------------------------------------ */
    // Constant function for properties

    void
    Func_Const( moris::Matrix< moris::DDRMat >&            aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    /* ------------------------------------------------------------------------ */

    // Dirichlet function for properties
    void
    Func_Dirichlet( moris::Matrix< moris::DDRMat >&        aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        // Coordinates of pont
        moris::Matrix< moris::DDRMat > tXp = aFIManager->get_IP_geometry_interpolator()->valx();

        aPropMatrix.set_size( 3, 3, 0.0 );

        if ( ( tXp( 0 ) < tDirichletLength ) and ( tXp( 1 ) < tDirichletLength ) and ( tXp( 2 ) < tDirichletLength ) )
        {
            aPropMatrix( 0, 0 ) = 1.0;
            aPropMatrix( 1, 1 ) = 1.0;
            aPropMatrix( 2, 2 ) = 1.0;
        }
    }

    /* ------------------------------------------------------------------------ */

    bool
    Output_Criterion( moris::tsa::Time_Solver* aTimeSolver )
    {
        return true;
    }

    /* ------------------------------------------------------------------------ */

    void
    OPTParameterList( Vector< Vector< ParameterList > >& aParameterLists )
    {
        aParameterLists( 0 ).push_back( prm::create_opt_problem_parameter_list() );

        aParameterLists.set( "is_optimization_problem", false );
    }

    void
    HMRParameterList( Vector< Vector< ParameterList > >& aParameterLists )
    {
        aParameterLists( 0 ).push_back( prm::create_hmr_parameter_list() );

        aParameterLists.set( "number_of_elements_per_dimension", tNumElemsPerDim );
        aParameterLists.set( "domain_dimensions", tDomainDims );
        aParameterLists.set( "domain_offset", tDomainOffset );
        aParameterLists.set( "domain_sidesets", tDomainSidesets );
        aParameterLists.set( "lagrange_output_meshes", "0" );

        aParameterLists.set( "lagrange_orders", std::to_string( 1 ) );
        aParameterLists.set( "lagrange_pattern", "0" );
        aParameterLists.set( "bspline_orders", std::to_string( 1 ) );
        aParameterLists.set( "bspline_pattern", "0" );

        aParameterLists.set( "lagrange_to_bspline", "0" );

        aParameterLists.set( "truncate_bsplines", 1 );
        aParameterLists.set( "refinement_buffer", tRefineBuffer );
        aParameterLists.set( "staircase_buffer", tRefineBuffer );
        aParameterLists.set( "initial_refinement", "0" );
        aParameterLists.set( "initial_refinement_pattern", "0" );

        aParameterLists.set( "use_number_aura", 1 );

        aParameterLists.set( "use_multigrid", 0 );
        aParameterLists.set( "severity_level", 0 );
    }

    void
    XTKParameterList( Vector< Vector< ParameterList > >& aParameterLists )
    {
        aParameterLists( 0 ).push_back( prm::create_xtk_parameter_list() );
        aParameterLists.set( "decompose", true );
        aParameterLists.set( "decomposition_type", "conformal" );
        aParameterLists.set( "enrich", true );
        aParameterLists.set( "basis_rank", "bspline" );
        aParameterLists.set( "enrich_mesh_indices", "0" );
        aParameterLists.set( "ghost_stab", tUseGhost );
        aParameterLists.set( "multigrid", false );
        aParameterLists.set( "verbose", true );
        aParameterLists.set( "print_enriched_ig_mesh", false );
        aParameterLists.set( "exodus_output_XTK_ig_mesh", true );
        aParameterLists.set( "high_to_low_dbl_side_sets", true );
    }

    void
    GENParameterList( Vector< Vector< ParameterList > >& aParameterLists )
    {
        aParameterLists( 0 ).push_back( prm::create_gen_parameter_list() );

        // Geometry parameter lists
        aParameterLists( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists.set( "field_function_name", "Inclusion" );
    }

    void
    MIGParameterList( Vector< Vector< ParameterList > >& aParameterLists )
    {

        // Main GEN parameter list
        aParameterLists( 0 ).push_back( prm::create_mig_parameter_list() );

        aParameterLists.set( "periodic_side_set_pair", tPeriodicSidePairs );
    }

    void
    FEMParameterList( Vector< Vector< ParameterList > >& aParameterLists )
    {
        // create a cell of cell of parameter list for fem

        //------------------------------------------------------------------------------

        // properties for inclusion

        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropDensityInner" );
        aParameterLists.set( "function_parameters", tDensityInner );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropYoungInner" );
        aParameterLists.set( "function_parameters", tEmodInner );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropPoisInner" );
        aParameterLists.set( "function_parameters", tPoisInner );
        aParameterLists.set( "value_function", "Func_Const" );

        // properties of boundary conditions
        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropDirichlet" );
        aParameterLists.set( "function_parameters", "0.0;0.0;0.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // properties for outer material
        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropDensityOuter" );
        aParameterLists.set( "function_parameters", tDensityOuter );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropYoungOuter" );
        aParameterLists.set( "function_parameters", tEmodOuter );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropPoisOuter" );
        aParameterLists.set( "function_parameters", tPoisOuter );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropEigenStrain" );
        aParameterLists.set( "function_parameters", tEigenStrain );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropSelect" );
        aParameterLists.set( "function_parameters", "0.0;0.0;0.0" );
        aParameterLists.set( "value_function", "Func_Dirichlet" );

        // time continuity weights

        //------------------------------------------------------------------------------

        // create parameter list for constitutive model 1
        aParameterLists( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        aParameterLists.set( "constitutive_name", "CMStrucLinIsoInner" );
        aParameterLists.set( "model_type",  fem::Model_Type::FULL ) ;
        aParameterLists.set( "constitutive_type",  fem::Constitutive_Type::STRUC_LIN_ISO ) ;
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY,UZ", "Displacement" ) );
        aParameterLists.set( "properties", "PropYoungInner,YoungsModulus;PropPoisInner,PoissonRatio;PropEigenStrain,EigenStrain" );

        // create parameter list for constitutive model 2
        aParameterLists( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        aParameterLists.set( "constitutive_name", "CMStrucLinIsoOuter" );
        aParameterLists.set( "model_type",  fem::Model_Type::FULL ) ;
        aParameterLists.set( "constitutive_type",  fem::Constitutive_Type::STRUC_LIN_ISO ) ;
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY,UZ", "Displacement" ) );
        aParameterLists.set( "properties", "PropYoungOuter,YoungsModulus;PropPoisOuter,PoissonRatio;PropEigenStrain,EigenStrain" );

        // create parameter list for constitutive model 2
        aParameterLists( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        aParameterLists.set( "constitutive_name", "CMStrucLinIsoPeriodicOuter" );
        aParameterLists.set( "model_type",  fem::Model_Type::FULL ) ;
        aParameterLists.set( "constitutive_type",  fem::Constitutive_Type::STRUC_LIN_ISO ) ;
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY,UZ", "Displacement" ) );
        aParameterLists.set( "properties", "PropYoungOuter,YoungsModulus;PropPoisOuter,PoissonRatio;PropEigenStrain,EigenStrain" );

        // create parameter list for constitutive model 1
        aParameterLists( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        aParameterLists.set( "constitutive_name", "CMStrucLinIsoPeriodicInner" );
        aParameterLists.set( "model_type",  fem::Model_Type::FULL ) ;
        aParameterLists.set( "constitutive_type",  fem::Constitutive_Type::STRUC_LIN_ISO ) ;
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY,UZ", "Displacement" ) );
        aParameterLists.set( "properties", "PropYoungInner,YoungsModulus;PropPoisInner,PoissonRatio;PropEigenStrain,EigenStrain" );

        //------------------------------------------------------------------------------

        // create parameter list for ghost stabilization parameter for outer material
        aParameterLists( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists.set( "stabilization_name", "SPGhostInner" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::GHOST_DISPL ) ;
        aParameterLists.set( "function_parameters", "0.01" );
        aParameterLists.set( "leader_properties", "PropYoungInner,Material" );

        // create parameter list for ghost stabilization parameter for outer material
        aParameterLists( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists.set( "stabilization_name", "SPGhostOuter" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::GHOST_DISPL ) ;
        aParameterLists.set( "function_parameters", "0.01" );
        aParameterLists.set( "leader_properties", "PropYoungOuter,Material" );

        // create parameter list for stabilization parameter 1
        aParameterLists( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists.set( "stabilization_name", "SPNitsche" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::DIRICHLET_NITSCHE ) ;
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "leader_properties", "PropYoungOuter,Material" );

        // create parameter list for Nitsche stabilization parameter for inclusion-outer material interface
        aParameterLists( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists.set( "stabilization_name", "SPInterfaceNitsche" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::NITSCHE_INTERFACE ) ;
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "leader_properties", "PropYoungInner,Material" );
        aParameterLists.set( "follower_properties", "PropYoungOuter,Material" );

        // create parameter list for Nitsche stabilization parameter for inclusion-outer material interface
        aParameterLists( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists.set( "stabilization_name", "SPPeriodicNitscheP00" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::NITSCHE_INTERFACE ) ;
        aParameterLists.set( "function_parameters", "10000.0" );
        aParameterLists.set( "leader_properties", "PropYoungOuter,Material" );
        aParameterLists.set( "follower_properties", "PropYoungOuter,Material" );

        // create parameter list for Nitsche stabilization parameter for inclusion-outer material interface
        aParameterLists( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists.set( "stabilization_name", "SPPeriodicNitscheP11" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::NITSCHE_INTERFACE ) ;
        aParameterLists.set( "function_parameters", "1000.0" );
        aParameterLists.set( "leader_properties", "PropYoungInner,Material" );
        aParameterLists.set( "follower_properties", "PropYoungInner,Material" );

        //------------------------------------------------------------------------------
        // create IWG for inclusion - bulk diffusion
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGBulkInner" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_BULK ) ;
        aParameterLists.set( "dof_residual", "UX,UY,UZ" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY,UZ" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIsoInner,ElastLinIso" );
        aParameterLists.set( "mesh_set_names", tInnerPhase );

        // create IWG for inclusion - bulk diffusion
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGBulkOuter" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_BULK ) ;
        aParameterLists.set( "dof_residual", "UX,UY,UZ" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY,UZ" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIsoOuter,ElastLinIso" );
        aParameterLists.set( "mesh_set_names", tOuterPhase );

        // create parameter list for IWG 2
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGDirichlet" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE ) ;
        aParameterLists.set( "dof_residual", "UX,UY,UZ" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY,UZ" );
        aParameterLists.set( "leader_properties", "PropDirichlet,Dirichlet;PropSelect,Select" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIsoOuter,ElastLinIso" );
        aParameterLists.set( "stabilization_parameters", "SPNitsche,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", tDBC );

        // create parameter list for interface conditions
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGInterfaceInnerOuter" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_INTERFACE_UNSYMMETRIC_NITSCHE ) ;
        aParameterLists.set( "dof_residual", "UX,UY,UZ" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY,UZ" );
        aParameterLists.set( "follower_dof_dependencies", "UX,UY,UZ" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIsoInner,ElastLinIso" );
        aParameterLists.set( "follower_constitutive_models", "CMStrucLinIsoOuter,ElastLinIso" );
        aParameterLists.set( "stabilization_parameters", "SPInterfaceNitsche,NitscheInterface" );
        aParameterLists.set( "mesh_set_names", tInterface );

        if ( tUseGhost )
        {
            // create IWG for outer material - ghost
            aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
            aParameterLists.set( "IWG_name", "IWGGPInnerDisp" );
            aParameterLists.set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
            aParameterLists.set( "dof_residual", "UX,UY,UZ" );
            aParameterLists.set( "leader_dof_dependencies", "UX,UY,UZ" );
            aParameterLists.set( "follower_dof_dependencies", "UX,UY,UZ" );
            aParameterLists.set( "stabilization_parameters", "SPGhostInner,GhostSP" );
            aParameterLists.set( "mesh_set_names", tInnerPhaseGhost );

            // create IWG for outer material - ghost
            aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
            aParameterLists.set( "IWG_name", "IWGGPOuterDisp" );
            aParameterLists.set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
            aParameterLists.set( "dof_residual", "UX,UY,UZ" );
            aParameterLists.set( "leader_dof_dependencies", "UX,UY,UZ" );
            aParameterLists.set( "follower_dof_dependencies", "UX,UY,UZ" );
            aParameterLists.set( "stabilization_parameters", "SPGhostOuter,GhostSP" );
            aParameterLists.set( "mesh_set_names", tOuterPhaseGhost );
        }

        // periodic IWG
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGIPeriodic" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_INTERFACE_UNSYMMETRIC_NITSCHE ) ;
        aParameterLists.set( "dof_residual", "UX,UY,UZ" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY,UZ" );
        aParameterLists.set( "follower_dof_dependencies", "UX,UY,UZ" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIsoOuter,ElastLinIso" );
        aParameterLists.set( "follower_constitutive_models", "CMStrucLinIsoPeriodicOuter,ElastLinIso" );
        aParameterLists.set( "stabilization_parameters", "SPPeriodicNitscheP00,NitscheInterface" );
        aParameterLists.set( "mesh_set_names", "P00" );
        //------------------------------------------------------------------------------
        // init IQI counter

        uint tIQICounter = 0;

        aParameterLists( 4 ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists.set( "IQI_name", "IQIBulkDISPX" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists.set( "dof_quantity", "UX,UY,UZ" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY,UZ" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "mesh_set_names", tTotalDomain );

        aParameterLists( 4 ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists.set( "IQI_name", "IQIBulkDISPY" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists.set( "dof_quantity", "UX,UY,UZ" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY,UZ" );
        aParameterLists.set( "vectorial_field_index", 1 );
        aParameterLists.set( "mesh_set_names", tTotalDomain );

        aParameterLists( 4 ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists.set( "IQI_name", "IQIYoungsModulus1" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::PROPERTY ) ;
        aParameterLists.set( "leader_properties", "PropYoungOuter,Property" );
        aParameterLists.set( "mesh_set_names", tOuterPhase );

        aParameterLists( 4 ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists.set( "IQI_name", "IQIYoungsModulus2" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::PROPERTY ) ;
        aParameterLists.set( "leader_properties", "PropYoungInner,Property" );
        aParameterLists.set( "mesh_set_names", tInnerPhase );

        aParameterLists( 4 ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists.set( "IQI_name", "IQIHomInner" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::HOMOGENIZED_CONSTITUTIVE ) ;
        aParameterLists.set( "leader_dof_dependencies", "UX,UY,UZ" );
        aParameterLists.set( "leader_properties", "PropEigenStrain,EigenStrain" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIsoInner,Elast" );
        aParameterLists.set( "mesh_set_names", tInnerPhase );

        aParameterLists( 4 ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists.set( "IQI_name", "IQIHomOuter" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::HOMOGENIZED_CONSTITUTIVE ) ;
        aParameterLists.set( "leader_dof_dependencies", "UX,UY,UZ" );
        aParameterLists.set( "leader_properties", "PropEigenStrain,EigenStrain" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIsoOuter,Elast" );
        aParameterLists.set( "mesh_set_names", tOuterPhase );

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        aParameterLists( 5 ).push_back( prm::create_computation_parameter_list() );
    }

    void
    SOLParameterList( Vector< Vector< ParameterList > >& aParameterLists )
    {

        aParameterLists( 0 ).push_back( add_parameter_list( sol::SolverType::AMESOS_IMPL ) );
        // aParameterLists.set( "ifpack_prec_type", "ILU");

        aParameterLists( 1 ).push_back( moris::prm::create_linear_solver_parameter_list() );

        aParameterLists( 2 ).push_back( moris::prm::create_nonlinear_algorithm_parameter_list() );
        aParameterLists.set( "NLA_combined_res_jac_assembly", false );

        aParameterLists( 3 ).push_back( moris::prm::create_nonlinear_solver_parameter_list() );
        aParameterLists.set( "NLA_DofTypes", "UX,UY,UZ" );

        aParameterLists( 4 ).push_back( moris::prm::create_time_solver_algorithm_parameter_list() );
        //         aParameterLists.set("TSA_Num_Time_Steps",     1 );
        //         aParameterLists.set("TSA_Time_Frame",         1.0 );

        aParameterLists( 5 ).push_back( moris::prm::create_time_solver_parameter_list() );
        aParameterLists.set( "TSA_DofTypes", "UX,UY,UZ" );
        aParameterLists.set( "TSA_Output_Indices", "0" );
        aParameterLists.set( "TSA_Output_Criteria", "Output_Criterion" );

        aParameterLists( 6 ).push_back( moris::prm::create_solver_warehouse_parameterlist() );
        aParameterLists.set( "SOL_save_operator_to_matlab", "jp2.dat" );

        aParameterLists( 7 ).push_back(  sol::PreconditionerType::NONE );
    }

    void
    MSIParameterList( Vector< Vector< ParameterList > >& aParameterLists )
    {
        aParameterLists( 0 ).push_back( prm::create_msi_parameter_list() );
        aParameterLists.set( "order_adofs_by_host", false );
    }

    void
    VISParameterList( Vector< Vector< ParameterList > >& aParameterLists )
    {
        aParameterLists( 0 ).push_back( prm::create_vis_parameter_list() );
        aParameterLists.set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName ) );
        aParameterLists.set( "Mesh_Type",  vis::VIS_Mesh_Type::STANDARD ) ;
        aParameterLists.set( "Set_Names", tTotalDomain );
        aParameterLists.set( "Field_Names", "UX,UY,YoungsModulus1,YoungsModulus2,IQIIn,IQIOut" );
        aParameterLists.set( "Field_Type", "NODAL,NODAL,NODAL,NODAL,GLOBAL,GLOBAL" );
        aParameterLists.set( "IQI_Names", "IQIBulkDISPX,IQIBulkDISPY,IQIYoungsModulus1,IQIYoungsModulus2,IQIHomInner,IQIHomOuter" );
        aParameterLists.set( "Save_Frequency", 1 );
    }

    /* ------------------------------------------------------------------------ */

    void
    MORISGENERALParameterList( Vector< Vector< ParameterList > >& aParameterLists )
    {
    }
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
