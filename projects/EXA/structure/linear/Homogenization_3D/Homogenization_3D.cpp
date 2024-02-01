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
#include "fn_PRM_FEM_Parameters.hpp"
#include "fn_PRM_MSI_Parameters.hpp"
#include "fn_PRM_SOL_Parameters.hpp"
#include "fn_PRM_VIS_Parameters.hpp"
#include "fn_PRM_HMR_Parameters.hpp"
#include "fn_PRM_GEN_Parameters.hpp"
#include "fn_PRM_XTK_Parameters.hpp"
#include "fn_PRM_OPT_Parameters.hpp"
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
            const moris::Cell< real >&     aGeometryParameters )
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
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    /* ------------------------------------------------------------------------ */

    // Dirichlet function for properties
    void
    Func_Dirichlet( moris::Matrix< moris::DDRMat >&        aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
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
    OPTParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_opt_problem_parameter_list();

        tParameterlist( 0 )( 0 ).set( "is_optimization_problem", false );
    }

    void
    HMRParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_hmr_parameter_list();

        tParameterlist( 0 )( 0 ).set( "number_of_elements_per_dimension", tNumElemsPerDim );
        tParameterlist( 0 )( 0 ).set( "domain_dimensions", tDomainDims );
        tParameterlist( 0 )( 0 ).set( "domain_offset", tDomainOffset );
        tParameterlist( 0 )( 0 ).set( "domain_sidesets", tDomainSidesets );
        tParameterlist( 0 )( 0 ).set( "lagrange_output_meshes", "0" );

        tParameterlist( 0 )( 0 ).set( "lagrange_orders", std::to_string( 1 ) );
        tParameterlist( 0 )( 0 ).set( "lagrange_pattern", "0" );
        tParameterlist( 0 )( 0 ).set( "bspline_orders", std::to_string( 1 ) );
        tParameterlist( 0 )( 0 ).set( "bspline_pattern", "0" );

        tParameterlist( 0 )( 0 ).set( "lagrange_to_bspline", "0" );

        tParameterlist( 0 )( 0 ).set( "truncate_bsplines", 1 );
        tParameterlist( 0 )( 0 ).set( "refinement_buffer", tRefineBuffer );
        tParameterlist( 0 )( 0 ).set( "staircase_buffer", tRefineBuffer );
        tParameterlist( 0 )( 0 ).set( "initial_refinement", "0" );
        tParameterlist( 0 )( 0 ).set( "initial_refinement_pattern", "0" );

        tParameterlist( 0 )( 0 ).set( "use_number_aura", 1 );

        tParameterlist( 0 )( 0 ).set( "use_multigrid", 0 );
        tParameterlist( 0 )( 0 ).set( "severity_level", 0 );
    }

    void
    XTKParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_xtk_parameter_list();
        tParameterlist( 0 )( 0 ).set( "decompose", true );
        tParameterlist( 0 )( 0 ).set( "decomposition_type", "conformal" );
        tParameterlist( 0 )( 0 ).set( "enrich", true );
        tParameterlist( 0 )( 0 ).set( "basis_rank", "bspline" );
        tParameterlist( 0 )( 0 ).set( "enrich_mesh_indices", "0" );
        tParameterlist( 0 )( 0 ).set( "ghost_stab", tUseGhost );
        tParameterlist( 0 )( 0 ).set( "multigrid", false );
        tParameterlist( 0 )( 0 ).set( "verbose", true );
        tParameterlist( 0 )( 0 ).set( "print_enriched_ig_mesh", false );
        tParameterlist( 0 )( 0 ).set( "exodus_output_XTK_ig_mesh", true );
        tParameterlist( 0 )( 0 ).set( "high_to_low_dbl_side_sets", true );
    }

    void
    GENParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 3 );
        tParameterlist( 0 ).resize( 1 );

        // Main GEN parameter list
        tParameterlist( 0 )( 0 ) = prm::create_gen_parameter_list();

        // init geometry counter
        uint tGeoCounter = 0;

        // Geometry parameter lists
        tParameterlist( 1 ).push_back( prm::create_user_defined_geometry_parameter_list() );
        tParameterlist( 1 )( tGeoCounter ).set( "field_function_name", "Inclusion" );
    }

    void
    MIGParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        // Main GEN parameter list
        tParameterlist( 0 )( 0 ) = prm::create_mig_parameter_list();

        tParameterlist( 0 )( 0 ).set( "periodic_side_set_pair", tPeriodicSidePairs );
    }

    void
    FEMParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterList )
    {
        // create a cell of cell of parameter list for fem
        tParameterList.resize( 8 );

        //------------------------------------------------------------------------------
        // init property counter
        uint tPropCounter = 0;

        // properties for inclusion

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ) = prm::create_property_parameter_list();
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropDensityInner" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tDensityInner );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropYoungInner" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tEmodInner );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropPoisInner" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tPoisInner );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // properties of boundary conditions
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropDirichlet" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "0.0;0.0;0.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // properties for outer material
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ) = prm::create_property_parameter_list();
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropDensityOuter" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tDensityOuter );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropYoungOuter" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tEmodOuter );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropPoisOuter" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tPoisOuter );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropEigenStrain" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tEigenStrain );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropSelect" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "0.0;0.0;0.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Dirichlet" );
        tPropCounter++;

        // time continuity weights

        //------------------------------------------------------------------------------

        // init CM counter
        uint tCMCounter = 0;

        // create parameter list for constitutive model 1
        tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_name", "CMStrucLinIsoInner" );
        tParameterList( 1 )( tCMCounter ).set( "model_type", static_cast< uint >( fem::Model_Type::FULL ) );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::STRUC_LIN_ISO ) );
        tParameterList( 1 )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY,UZ", "Displacement" ) );
        tParameterList( 1 )( tCMCounter ).set( "properties", "PropYoungInner,YoungsModulus;PropPoisInner,PoissonRatio;PropEigenStrain,EigenStrain" );
        tCMCounter++;

        // create parameter list for constitutive model 2
        tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_name", "CMStrucLinIsoOuter" );
        tParameterList( 1 )( tCMCounter ).set( "model_type", static_cast< uint >( fem::Model_Type::FULL ) );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::STRUC_LIN_ISO ) );
        tParameterList( 1 )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY,UZ", "Displacement" ) );
        tParameterList( 1 )( tCMCounter ).set( "properties", "PropYoungOuter,YoungsModulus;PropPoisOuter,PoissonRatio;PropEigenStrain,EigenStrain" );
        tCMCounter++;

        // create parameter list for constitutive model 2
        tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_name", "CMStrucLinIsoPeriodicOuter" );
        tParameterList( 1 )( tCMCounter ).set( "model_type", static_cast< uint >( fem::Model_Type::FULL ) );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::STRUC_LIN_ISO ) );
        tParameterList( 1 )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY,UZ", "Displacement" ) );
        tParameterList( 1 )( tCMCounter ).set( "properties", "PropYoungOuter,YoungsModulus;PropPoisOuter,PoissonRatio;PropEigenStrain,EigenStrain" );
        tCMCounter++;

        // create parameter list for constitutive model 1
        tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_name", "CMStrucLinIsoPeriodicInner" );
        tParameterList( 1 )( tCMCounter ).set( "model_type", static_cast< uint >( fem::Model_Type::FULL ) );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::STRUC_LIN_ISO ) );
        tParameterList( 1 )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY,UZ", "Displacement" ) );
        tParameterList( 1 )( tCMCounter ).set( "properties", "PropYoungInner,YoungsModulus;PropPoisInner,PoissonRatio;PropEigenStrain,EigenStrain" );
        tCMCounter++;

        //------------------------------------------------------------------------------
        // init SP counter
        uint tSPCounter = 0;

        // create parameter list for ghost stabilization parameter for outer material
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", "SPGhostInner" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::GHOST_DISPL ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", "0.01" );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", "PropYoungInner,Material" );
        tSPCounter++;

        // create parameter list for ghost stabilization parameter for outer material
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", "SPGhostOuter" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::GHOST_DISPL ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", "0.01" );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", "PropYoungOuter,Material" );
        tSPCounter++;

        // create parameter list for stabilization parameter 1
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", "SPNitsche" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::DIRICHLET_NITSCHE ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", "100.0" );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", "PropYoungOuter,Material" );
        tSPCounter++;

        // create parameter list for Nitsche stabilization parameter for inclusion-outer material interface
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", "SPInterfaceNitsche" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::NITSCHE_INTERFACE ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", "100.0" );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", "PropYoungInner,Material" );
        tParameterList( 2 )( tSPCounter ).set( "follower_properties", "PropYoungOuter,Material" );
        tSPCounter++;

        // create parameter list for Nitsche stabilization parameter for inclusion-outer material interface
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", "SPPeriodicNitscheP00" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::NITSCHE_INTERFACE ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", "10000.0" );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", "PropYoungOuter,Material" );
        tParameterList( 2 )( tSPCounter ).set( "follower_properties", "PropYoungOuter,Material" );
        tSPCounter++;

        // create parameter list for Nitsche stabilization parameter for inclusion-outer material interface
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", "SPPeriodicNitscheP11" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::NITSCHE_INTERFACE ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", "1000.0" );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", "PropYoungInner,Material" );
        tParameterList( 2 )( tSPCounter ).set( "follower_properties", "PropYoungInner,Material" );
        tSPCounter++;

        //------------------------------------------------------------------------------
        // init IWG counter
        uint tIWGCounter = 0;

        // create IWG for inclusion - bulk diffusion
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGBulkInner" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_BULK ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "UX,UY,UZ" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY,UZ" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMStrucLinIsoInner,ElastLinIso" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tInnerPhase );
        tIWGCounter++;

        // create IWG for inclusion - bulk diffusion
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGBulkOuter" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_BULK ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "UX,UY,UZ" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY,UZ" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMStrucLinIsoOuter,ElastLinIso" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tOuterPhase );
        tIWGCounter++;

        // create parameter list for IWG 2
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGDirichlet" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "UX,UY,UZ" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY,UZ" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties", "PropDirichlet,Dirichlet;PropSelect,Select" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMStrucLinIsoOuter,ElastLinIso" );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPNitsche,DirichletNitsche" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tDBC );
        tIWGCounter++;

        // create parameter list for interface conditions
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGInterfaceInnerOuter" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_INTERFACE_UNSYMMETRIC_NITSCHE ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "UX,UY,UZ" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY,UZ" );
        tParameterList( 3 )( tIWGCounter ).set( "follower_dof_dependencies", "UX,UY,UZ" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMStrucLinIsoInner,ElastLinIso" );
        tParameterList( 3 )( tIWGCounter ).set( "follower_constitutive_models", "CMStrucLinIsoOuter,ElastLinIso" );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPInterfaceNitsche,NitscheInterface" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tInterface );
        tIWGCounter++;

        if ( tUseGhost )
        {
            // create IWG for outer material - ghost
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGGPInnerDisp" );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::GHOST_NORMAL_FIELD ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "UX,UY,UZ" );
            tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY,UZ" );
            tParameterList( 3 )( tIWGCounter ).set( "follower_dof_dependencies", "UX,UY,UZ" );
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPGhostInner,GhostSP" );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tInnerPhaseGhost );
            tIWGCounter++;

            // create IWG for outer material - ghost
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGGPOuterDisp" );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::GHOST_NORMAL_FIELD ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "UX,UY,UZ" );
            tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY,UZ" );
            tParameterList( 3 )( tIWGCounter ).set( "follower_dof_dependencies", "UX,UY,UZ" );
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPGhostOuter,GhostSP" );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tOuterPhaseGhost );
            tIWGCounter++;
        }

        // periodic IWG
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGIPeriodic" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_INTERFACE_UNSYMMETRIC_NITSCHE ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "UX,UY,UZ" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY,UZ" );
        tParameterList( 3 )( tIWGCounter ).set( "follower_dof_dependencies", "UX,UY,UZ" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMStrucLinIsoOuter,ElastLinIso" );
        tParameterList( 3 )( tIWGCounter ).set( "follower_constitutive_models", "CMStrucLinIsoPeriodicOuter,ElastLinIso" );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPPeriodicNitscheP00,NitscheInterface" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", "P00" );
        tIWGCounter++;
        //------------------------------------------------------------------------------
        // init IQI counter

        uint tIQICounter = 0;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkDISPX" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::DOF ) );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "UX,UY,UZ" );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "UX,UY,UZ" );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 0 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tTotalDomain );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkDISPY" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::DOF ) );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "UX,UY,UZ" );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "UX,UY,UZ" );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 1 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tTotalDomain );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIYoungsModulus1" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::PROPERTY ) );
        tParameterList( 4 )( tIQICounter ).set( "leader_properties", "PropYoungOuter,Property" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tOuterPhase );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIYoungsModulus2" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::PROPERTY ) );
        tParameterList( 4 )( tIQICounter ).set( "leader_properties", "PropYoungInner,Property" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tInnerPhase );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIHomInner" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::HOMOGENIZED_CONSTITUTIVE ) );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "UX,UY,UZ" );
        tParameterList( 4 )( tIQICounter ).set( "leader_properties", "PropEigenStrain,EigenStrain" );
        tParameterList( 4 )( tIQICounter ).set( "leader_constitutive_models", "CMStrucLinIsoInner,Elast" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tInnerPhase );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIHomOuter" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::HOMOGENIZED_CONSTITUTIVE ) );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "UX,UY,UZ" );
        tParameterList( 4 )( tIQICounter ).set( "leader_properties", "PropEigenStrain,EigenStrain" );
        tParameterList( 4 )( tIQICounter ).set( "leader_constitutive_models", "CMStrucLinIsoOuter,Elast" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tOuterPhase );
        tIQICounter++;

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        tParameterList( 5 ).resize( 1 );
        tParameterList( 5 )( 0 ) = prm::create_computation_parameter_list();
    }

    void
    SOLParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 8 );

        tParameterlist( 0 ).push_back( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL ) );
        // tParameterlist( 0 )( 0 ).set( "ifpack_prec_type", "ILU");

        tParameterlist( 1 ).push_back( moris::prm::create_linear_solver_parameter_list() );

        tParameterlist( 2 ).push_back( moris::prm::create_nonlinear_algorithm_parameter_list() );
        tParameterlist( 2 )( 0 ).set( "NLA_combined_res_jac_assembly", false );

        tParameterlist( 3 ).push_back( moris::prm::create_nonlinear_solver_parameter_list() );
        tParameterlist( 3 )( 0 ).set( "NLA_DofTypes", "UX,UY,UZ" );

        tParameterlist( 4 ).push_back( moris::prm::create_time_solver_algorithm_parameter_list() );
        //         tParameterlist( 4 )( 0 ).set("TSA_Num_Time_Steps",     1 );
        //         tParameterlist( 4 )( 0 ).set("TSA_Time_Frame",         1.0 );

        tParameterlist( 5 ).push_back( moris::prm::create_time_solver_parameter_list() );
        tParameterlist( 5 )( 0 ).set( "TSA_DofTypes", "UX,UY,UZ" );
        tParameterlist( 5 )( 0 ).set( "TSA_Output_Indices", "0" );
        tParameterlist( 5 )( 0 ).set( "TSA_Output_Criteria", "Output_Criterion" );

        tParameterlist( 6 ).push_back( moris::prm::create_solver_warehouse_parameterlist() );
        tParameterlist( 6 )( 0 ).set( "SOL_save_operator_to_matlab", "jp2.dat" );

        tParameterlist( 7 ).push_back( moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::NONE ) );
    }

    void
    MSIParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_msi_parameter_list();
        tParameterlist( 0 )( 0 ).set( "order_adofs_by_host", false );
    }

    void
    VISParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_vis_parameter_list();
        tParameterlist( 0 )( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName ) );
        tParameterlist( 0 )( 0 ).set( "Mesh_Type", static_cast< uint >( vis::VIS_Mesh_Type::STANDARD ) );
        tParameterlist( 0 )( 0 ).set( "Set_Names", tTotalDomain );
        tParameterlist( 0 )( 0 ).set( "Field_Names", "UX,UY,YoungsModulus1,YoungsModulus2,IQIIn,IQIOut" );
        tParameterlist( 0 )( 0 ).set( "Field_Type", "NODAL,NODAL,NODAL,NODAL,GLOBAL,GLOBAL" );
        tParameterlist( 0 )( 0 ).set( "IQI_Names", "IQIBulkDISPX,IQIBulkDISPY,IQIYoungsModulus1,IQIYoungsModulus2,IQIHomInner,IQIHomOuter" );
        tParameterlist( 0 )( 0 ).set( "Save_Frequency", 1 );
    }

    /* ------------------------------------------------------------------------ */

    void
    MORISGENERALParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
    }
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif