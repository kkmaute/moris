/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * linear_contact_phasebased.cpp
 *
 */

#include "cl_Bitset.hpp"
#include "cl_DLA_Linear_Solver_Aztec.hpp"
#include "cl_DLA_Solver_Interface.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_MSI_Equation_Object.hpp"
#include "cl_Matrix.hpp"
#include "cl_TSA_Time_Solver.hpp"
#include "parameters.hpp"
#include "fn_equal_to.hpp"
#include "fn_norm.hpp"
#include "linalg_typedefs.hpp"
#include "moris_typedefs.hpp"
#include <iostream>
#include <string>

#include "AztecOO.h"
//---------------------------------------------------------------

#define F2STR( func ) #func

using moris::fem::IWG_Type;

//---------------------------------------------------------------

// global variable for interpolation order
extern uint gInterpolationOrder;

// global variable for case index
extern uint gCaseIndex;

//---------------------------------------------------------------

#ifdef __cplusplus
extern "C" {
#endif
namespace moris
{
    /* -------------------------------------- Geometry Setting -------------------------------------- */
    real tHeight = 1.0;
    real tWidth  = 1.0;

    uint tNumElems            = 11;
    uint tInterfaceRefinement = 0;
    uint tNumElemsX           = tNumElems;
    uint tNumElemsY           = tNumElems;

    real tInitialGap = 0.01;
    // real tPenetration = 0.05 * 0.85;
    real tPenetration          = 0.02;
    real tVerticalDisplacement = ( tInitialGap + tPenetration ) / 2;    // both bodies move

    /* --------------------------------------- Top Geometry ----------------------------------------- */
    bool tTopIsPlane         = false;
    bool tTopHasSymmetryWest = true;
    bool tTopHasSymmetryEast = true;

    real tInterfaceTopAngle = 0 * M_PI / 180.0;
    real tTopYShift         = tInitialGap / 2;
    real tTopXShift         = 0.0;
    real tTopParabolaFactor = 0.5;

    /* --------------------------------------- Bottom Geometry -------------------------------------- */
    bool tBottomIsPlane         = false;
    bool tBottomHasSymmetryWest = true;
    bool tBottomHasSymmetryEast = true;

    real tInterfaceBottomAngle = 0 * M_PI / 180.0;
    real tBottomYShift         = -tInitialGap / 2;
    real tBottomXShift         = 0.0;
    real tBottomParabolaFactor = -0.5;

    /* ------------------------------------- Material Parameters ------------------------------------ */
    real tTopEmod    = 1000.0;
    real tBottomEmod = 1000.0;
    real tTopPois    = 0.3;
    real tBottomPois = 0.3;

    /* --------------------------------------- Domain Setting --------------------------------------- */
    real tDomainOffsetX = -0.5 * tWidth;
    real tDomainOffsetY = -0.5 * tHeight;

    std::string tInterpolationOrder     = std::to_string( gInterpolationOrder );
    std::string tDomainTop              = "HMR_dummy_n_p1,HMR_dummy_c_p1";
    std::string tDomainBottom           = "HMR_dummy_n_p2,HMR_dummy_c_p2";
    std::string tContactInterfaceTop    = "ncss|iside_b0_1_b1_0|iside_b0_2_b1_0";
    std::string tContactInterfaceBottom = "ncss|iside_b0_2_b1_0|iside_b0_1_b1_0";
    std::string tContactInterface       = tContactInterfaceBottom + "," + tContactInterfaceTop;
    std::string tDomain                 = tDomainTop + "," + tDomainBottom;

    /* ------------------------------------------ File I/O ------------------------------------------ */
    std::string tOutputFileName =
            "Parabolic_Indenter_NonLinear_Case_" + std::to_string( gCaseIndex );

    /* ------------------------------------------- Solver ------------------------------------------- */
    int tMaxIterations = 200;

    real tRelResNormDrop     = 1e-8;
    real tRelaxation         = 0.6;
    real tFDPerturbationSize = 1e-8;

    fem::Perturbation_Type tFDPerturbationStrategy = fem::Perturbation_Type::RELATIVE;
    fem::FDScheme_Type     tFDScheme               = fem::FDScheme_Type::POINT_3_CENTRAL;

    int  tLoadControlSteps    = 20;
    real tLoadControlFactor   = 0.01;
    real tLoadControlRelRes   = 1e0;
    real tLoadControlExponent = 1.0;

    real tRemapResidualChange     = -1e-2;
    int  tRemapLoadStepFrequency  = 1;
    int  tRemapIterationFrequency = 1;
    auto tRaytracingStrategy      = sol::SolverRaytracingStrategy::EveryNthIteration;
    real tMaxNegativeRayLength    = -0.08;
    real tMaxPositiveRayLength    = 0.5;

    mtk::Integration_Order tNonconformalIntegrationOrder = mtk::Integration_Order::BAR_2;

    /* ------------------------------------------- Contact ------------------------------------------ */
    std::string tContactType          = "mlika";
    std::string tContactBias          = "unsymmetric";
    real        tContactStabilization = 30.00;

    /* -------------------------------------- Control Variables ------------------------------------- */
    bool tOnlyGenerateMesh   = false;
    bool tUseGhost           = true;
    bool tTopHasDirichlet    = true;
    bool tBottomHasDirichlet = true;

    /* ---------------------------------------------------------------------------------------------- */
    /*                                         Field Functions                                        */
    /* ---------------------------------------------------------------------------------------------- */

    /* --------------------------------------- Phase Indexing --------------------------------------- */
    uint Phase_Index_Split( const Bitset< 2 > &aGeometrySigns )
    {
        if ( aGeometrySigns.test( 1 ) )
        {
            if ( aGeometrySigns.test( 0 ) )
            {
                return 1;    // upper domain
            }
            return 0;    // interface
        }
        return 2;    // bottom
    }

    /* ----------------------------------- Geometry Field Function ---------------------------------- */

    real tMinimumLevelSetValue = 1.0e-8;

    // Level set function for a parabola that splits the domain into a top and a bottom part.
    real Parabola( const Matrix< DDRMat > &aCoordinates,
            const Vector< real >          &aGeometryParameters )
    {

        real tXShift = aGeometryParameters( 0 );
        real tYShift = aGeometryParameters( 1 );
        real tFactor = aGeometryParameters( 2 );

        real x = aCoordinates( 0 ) - tXShift;
        real y = aCoordinates( 1 ) - tYShift;

        // Compute Signed-Distance field
        real tVal = y - ( tFactor * x * x );

        return ( std::abs( tVal ) < tMinimumLevelSetValue ) ? tMinimumLevelSetValue : tVal;
    }

    /* ----------------------------------- Property Field Function ---------------------------------- */

    void Func_Const( Matrix< DDRMat >       &aPropMatrix,
            Vector< Matrix< DDRMat > >      &aParameters,
            fem::Field_Interpolator_Manager *aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    void Func_Dirichlet( Matrix< DDRMat >   &aPropMatrix,
            Vector< Matrix< DDRMat > >      &aParameters,
            fem::Field_Interpolator_Manager *aFIManager )
    {
        real tLoadFactor = gLogger.get_action_data( "NonLinearAlgorithm", "Newton", "Solve", "LoadFactor" );
        aPropMatrix      = aParameters( 0 ) * tLoadFactor;
    }

    void Func_Linear_Along_X( Matrix< DDRMat > &aPropMatrix,
            Vector< Matrix< DDRMat > >         &aParameters,
            fem::Field_Interpolator_Manager    *aFIManager )
    {

        real tLoadFactor = gLogger.get_action_data( "NonLinearAlgorithm", "Newton", "Solve", "LoadFactor" );
        auto tXp         = aFIManager->get_IP_geometry_interpolator()->valx();
        auto a           = aParameters( 0 )( 0 );
        auto b           = aParameters( 0 )( 1 );

        aPropMatrix = {
            { 0.0 },
            { a * tXp( 0 ) + b },
        };
        aPropMatrix = aPropMatrix * tLoadFactor;
    }

    void Func_Select_X( Matrix< DDRMat >    &aPropMatrix,
            Vector< Matrix< DDRMat > >      &aParameters,
            fem::Field_Interpolator_Manager *aFIManager )
    {
        aPropMatrix.set_size( 2, 2, 0.0 );
        aPropMatrix( 0, 0 ) = 1.0;
    }

    void Func_Body_Load( Matrix< DDRMat >   &aPropMatrix,
            Vector< Matrix< DDRMat > >      &aParameters,
            fem::Field_Interpolator_Manager *aFIManager )
    {

        auto tXp    = aFIManager->get_IG_geometry_interpolator()->valx();
        auto param  = aParameters( 0 )( 0 );
        aPropMatrix = {
            { 0.0 },
            { param },
        };
    }

    bool Output_Criterion( tsa::Time_Solver *aTimeSolver ) { return true; }

    /* ---------------------------------------------------------------------------------------------- */
    /*                                           ### OPT ###                                          */
    /* ---------------------------------------------------------------------------------------------- */
    void OPTParameterList( Module_Parameter_Lists &aParameterLists )
    {
        aParameterLists.set( "is_optimization_problem", false );
    }

    /* ---------------------------------------------------------------------------------------------- */
    /*                                           ### HMR ###                                          */
    /* ---------------------------------------------------------------------------------------------- */
    void HMRParameterList( Module_Parameter_Lists &aParameterLists )
    {
        /* --------------------------------------- Domain Settings ------------------------------------ */
        aParameterLists.set( "number_of_elements_per_dimension", tNumElemsX, tNumElemsY );
        aParameterLists.set( "domain_dimensions", tWidth, tHeight );
        aParameterLists.set( "domain_offset", tDomainOffsetX, tDomainOffsetY );

        /* ---------------------------------------- Lagrange Mesh ------------------------------------- */
        aParameterLists.set( "lagrange_pattern", "0" );
        aParameterLists.set( "lagrange_orders", tInterpolationOrder );
        aParameterLists.set( "lagrange_output_meshes", "0" );

        /* -------------------------------------- B-Spline Mesh --------------------------------------- */
        aParameterLists.set( "bspline_pattern", "0" );
        aParameterLists.set( "bspline_orders", tInterpolationOrder );
        aParameterLists.set( "lagrange_to_bspline", "0" );
        aParameterLists.set( "truncate_bsplines", true );

        /* ----------------------------------------- Refinement --------------------------------------- */
        aParameterLists.set( "refinement_buffer", 2 );
        aParameterLists.set( "staircase_buffer", 2 );
    }

    /* ---------------------------------------------------------------------------------------------- */
    /*                                           ### XTK ###                                          */
    /* ---------------------------------------------------------------------------------------------- */
    void XTKParameterList( Module_Parameter_Lists &aParameterLists )
    {
        aParameterLists.set( "decompose", true );
        aParameterLists.set( "decomposition_type", "conformal" );
        aParameterLists.set( "enrich_mesh_indices", "0" );
        aParameterLists.set( "ghost_stab", tUseGhost );
        aParameterLists.set( "multigrid", false );
        aParameterLists.set( "verbose", true );
        aParameterLists.set( "print_enriched_ig_mesh", true );
        aParameterLists.set( "exodus_output_XTK_ig_mesh", true );
        aParameterLists.set( "high_to_low_dbl_side_sets", true );
        aParameterLists.set( "only_generate_xtk_temp", tOnlyGenerateMesh );
    }

    /* ---------------------------------------------------------------------------------------------- */
    /*                                           ### GEN ###                                          */
    /* ---------------------------------------------------------------------------------------------- */
    void GENParameterList( Module_Parameter_Lists &aParameterLists )
    {
        /* -------------------------------------------------------------------------------------------- */
        /*                                        GEN Parameter                                         */
        /* -------------------------------------------------------------------------------------------- */
        aParameterLists.set( "number_of_phases", 3 );
        aParameterLists.set( "phase_function_name", F2STR( Phase_Index_Split ) );
        /* -------------------------------------------------------------------------------------------- */
        /*                                      Geometry Parameter                                      */
        /* -------------------------------------------------------------------------------------------- */

        /* ----------------------------------------- Top Plane ---------------------------------------- */
        if ( tTopIsPlane )
        {
            // calculate the components of the rotated normal
            real tNormalX = std::cos( tInterfaceTopAngle + M_PI_2 );
            real tNormalY = std::sin( tInterfaceTopAngle + M_PI_2 );

            aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
            aParameterLists.set( "center_x", tTopXShift );
            aParameterLists.set( "center_y", tTopYShift );
            aParameterLists.set( "normal_x", tNormalX );
            aParameterLists.set( "normal_y", tNormalY );
            aParameterLists.set( "number_of_refinements", tInterfaceRefinement );
            aParameterLists.set( "refinement_mesh_index", 0 );
        }
        else
        {
            // use a parabola for the top side
            aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
            aParameterLists.set( "field_function_name", "Parabola" );
            aParameterLists.set( "number_of_refinements", tInterfaceRefinement );
            aParameterLists.set( "refinement_mesh_index", 0 );
            aParameterLists.insert( "variable_1", tTopXShift );
            aParameterLists.insert( "variable_2", tTopYShift );
            aParameterLists.insert( "variable_3", tTopParabolaFactor );
        }

        /* --------------------------------------- Bottom Plane --------------------------------------- */
        if ( tBottomIsPlane )
        {
            // calculate the components of the rotated normal
            real tNormalX = std::cos( tInterfaceBottomAngle + M_PI_2 );
            real tNormalY = std::sin( tInterfaceBottomAngle + M_PI_2 );

            aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
            aParameterLists.set( "center_x", tBottomXShift );
            aParameterLists.set( "center_y", tBottomYShift );
            aParameterLists.set( "normal_x", tNormalX );
            aParameterLists.set( "normal_y", tNormalY );
            aParameterLists.set( "number_of_refinements", tInterfaceRefinement );
            aParameterLists.set( "refinement_mesh_index", 0 );
        }
        else
        {
            // use a parabola for the bottom side
            aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
            aParameterLists.set( "field_function_name", "Parabola" );
            aParameterLists.set( "number_of_refinements", tInterfaceRefinement );
            aParameterLists.set( "refinement_mesh_index", 0 );
            aParameterLists.insert( "variable_1", tBottomXShift );
            aParameterLists.insert( "variable_2", tBottomYShift );
            aParameterLists.insert( "variable_3", tBottomParabolaFactor );
        }
    }

    /* ---------------------------------------------------------------------------------------------- */
    /*                                           ### FEM ###                                          */
    /* ---------------------------------------------------------------------------------------------- */
    void FEMParameterList( Module_Parameter_Lists &aParameterLists )
    {
        /* -------------------------------------------------------------------------------------------- */
        /*                                            Phases                                            */
        /* -------------------------------------------------------------------------------------------- */
        Parameter_List pl = prm::create_phase_parameter_list();
        pl.set( "phase_name", "PhaseVoid" );
        pl.set( "phase_indices", "0" );
        aParameterLists( FEM::PHASES ).add_parameter_list( pl );

        pl = prm::create_phase_parameter_list();
        pl.set( "phase_name", "PhaseTop" );
        pl.set( "phase_indices", "1" );
        aParameterLists( FEM::PHASES ).add_parameter_list( pl );

        pl = prm::create_phase_parameter_list();
        pl.set( "phase_name", "PhaseBottom" );
        pl.set( "phase_indices", "2" );
        aParameterLists( FEM::PHASES ).add_parameter_list( pl );

        pl = prm::create_phase_parameter_list();
        pl.set( "phase_name", "PhaseAll" );
        pl.set( "phase_indices", "1,2" );
        aParameterLists( FEM::PHASES ).add_parameter_list( pl );

        /* -------------------------------------------------------------------------------------------- */
        /*                                           Properties                                         */
        /* -------------------------------------------------------------------------------------------- */

        /* -------------------------------------- Youngs Modulus -------------------------------------- */
        pl = prm::create_property_parameter_list();
        pl.set( "property_name", "PropYoungsTop" );
        pl.set( "function_parameters", std::to_string( tTopEmod ) );
        pl.set( "value_function", F2STR( Func_Const ) );
        aParameterLists( FEM::PROPERTIES ).add_parameter_list( pl );

        pl = prm::create_property_parameter_list();
        pl.set( "property_name", "PropYoungsBottom" );
        pl.set( "function_parameters", std::to_string( tBottomEmod ) );
        pl.set( "value_function", F2STR( Func_Const ) );
        aParameterLists( FEM::PROPERTIES ).add_parameter_list( pl );

        /* ------------------------------------------ Poisson ----------------------------------------- */
        pl = prm::create_property_parameter_list();
        pl.set( "property_name", "PropPoissonTop" );
        pl.set( "function_parameters", std::to_string( tTopPois ) );
        pl.set( "value_function", F2STR( Func_Const ) );
        aParameterLists( FEM::PROPERTIES ).add_parameter_list( pl );

        pl = prm::create_property_parameter_list();
        pl.set( "property_name", "PropPoissonBottom" );
        pl.set( "function_parameters", std::to_string( tBottomPois ) );
        pl.set( "value_function", F2STR( Func_Const ) );
        aParameterLists( FEM::PROPERTIES ).add_parameter_list( pl );

        /* --------------------------------------- Dirichlet Top -------------------------------------- */
        if ( tTopHasDirichlet )
        {
            pl = prm::create_property_parameter_list();
            pl.set( "property_name", "PropDirichletTop" );
            pl.set( "function_parameters", "0.0;" + std::to_string( -tVerticalDisplacement ) );
            pl.set( "value_function", F2STR( Func_Dirichlet ) );
            aParameterLists( FEM::PROPERTIES ).add_parameter_list( pl );
        }

        /* ------------------------------------- Dirichlet Bottom ------------------------------------- */
        if ( tBottomHasDirichlet )
        {
            pl = prm::create_property_parameter_list();
            pl.set( "property_name", "PropDirichletBottom" );
            pl.set( "function_parameters", "0.0;" + std::to_string( tVerticalDisplacement ) );
            pl.set( "value_function", F2STR( Func_Dirichlet ) );
            aParameterLists( FEM::PROPERTIES ).add_parameter_list( pl );
        }

        pl = prm::create_property_parameter_list();
        pl.set( "property_name", "PropFixed" );
        pl.set( "function_parameters", "0.0;0.0" );
        pl.set( "value_function", F2STR( Func_Const ) );
        aParameterLists( FEM::PROPERTIES ).add_parameter_list( pl );

        /* -------------------------------------- Symmetry Plane -------------------------------------- */
        pl = prm::create_property_parameter_list();
        pl.set( "property_name", "PropSelectX" );
        pl.set( "function_parameters", "0.0;0.0" );
        pl.set( "value_function", F2STR( Func_Select_X ) );
        aParameterLists( FEM::PROPERTIES ).add_parameter_list( pl );

        /* -------------------------------------------------------------------------------------------- */
        /*                                       Constitutive Model                                     */
        /* -------------------------------------------------------------------------------------------- */

        /* ----------------------------------- Elastic Material Top ----------------------------------- */
        pl = prm::create_constitutive_model_parameter_list();
        pl.set( "constitutive_name", "MaterialTop" );
        pl.set( "phase_name", "PhaseTop" );
        pl.set( "constitutive_type",
                (uint)fem::Constitutive_Type::STRUC_NON_LIN_ISO_SAINT_VENANT_KIRCHHOFF );
        pl.set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
        pl.set( "properties",
                "PropYoungsTop,YoungsModulus;"
                "PropPoissonTop,PoissonRatio" );
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list( pl );

        /* -------------------------------- Elastic Material Bottom ----------------------------------- */
        pl = prm::create_constitutive_model_parameter_list();
        pl.set( "constitutive_name", "MaterialBottom" );
        pl.set( "phase_name", "PhaseBottom" );
        pl.set( "constitutive_type",
                (uint)fem::Constitutive_Type::STRUC_NON_LIN_ISO_SAINT_VENANT_KIRCHHOFF );
        pl.set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
        pl.set( "properties",
                "PropYoungsBottom,YoungsModulus;"
                "PropPoissonBottom,PoissonRatio" );
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list( pl );

        /* -------------------------------------------------------------------------------------------- */
        /*                                     Stabilization Parameter                                  */
        /* -------------------------------------------------------------------------------------------- */

        pl = prm::create_stabilization_parameter_parameter_list();
        pl.set( "stabilization_name", "SPDirichletNitscheTop" );
        pl.set( "stabilization_type", (uint)fem::Stabilization_Type::DIRICHLET_NITSCHE );
        pl.set( "function_parameters", "100.0" );
        pl.set( "leader_phase_name", "PhaseTop" );
        pl.set( "leader_properties", "PropYoungsTop,Material" );
        aParameterLists( FEM::STABILIZATION ).add_parameter_list( pl );

        pl = prm::create_stabilization_parameter_parameter_list();
        pl.set( "stabilization_name", "SPDirichletNitscheBottom" );
        pl.set( "stabilization_type", (uint)fem::Stabilization_Type::DIRICHLET_NITSCHE );
        pl.set( "function_parameters", "100.0" );
        pl.set( "leader_phase_name", "PhaseBottom" );
        pl.set( "leader_properties", "PropYoungsBottom,Material" );
        aParameterLists( FEM::STABILIZATION ).add_parameter_list( pl );

        pl = prm::create_stabilization_parameter_parameter_list();
        pl.set( "stabilization_name", "SPContactInterfaceTop" );
        pl.set( "stabilization_type", (uint)fem::Stabilization_Type::DIRICHLET_NITSCHE );
        pl.set( "function_parameters", std::to_string( tContactStabilization ) + ";1.0" );
        pl.set( "leader_phase_name", "PhaseTop" );
        pl.set( "leader_properties", "PropYoungsTop,Material" );
        aParameterLists( FEM::STABILIZATION ).add_parameter_list( pl );

        pl = prm::create_stabilization_parameter_parameter_list();
        pl.set( "stabilization_name", "SPContactInterfaceBottom" );
        pl.set( "stabilization_type", (uint)fem::Stabilization_Type::DIRICHLET_NITSCHE );
        pl.set( "function_parameters", std::to_string( tContactStabilization ) + ";1.0" );
        pl.set( "leader_phase_name", "PhaseBottom" );
        pl.set( "leader_properties", "PropYoungsBottom,Material" );
        aParameterLists( FEM::STABILIZATION ).add_parameter_list( pl );

        // create ghost penalty displacement Matrix
        pl = prm::create_stabilization_parameter_parameter_list();
        pl.set( "stabilization_name", "SPGPDisplTop" );
        pl.set( "stabilization_type", (uint)fem::Stabilization_Type::GHOST_DISPL );
        pl.set( "function_parameters", "0.05" );
        pl.set( "leader_phase_name", "PhaseTop" );
        pl.set( "follower_phase_name", "PhaseTop" );
        pl.set( "leader_properties", "PropYoungsTop,Material" );
        aParameterLists( FEM::STABILIZATION ).add_parameter_list( pl );

        // create ghost penalty  displacement Fiber
        pl = prm::create_stabilization_parameter_parameter_list();
        pl.set( "stabilization_name", "SPGPDisplBottom" );
        pl.set( "stabilization_type", (uint)fem::Stabilization_Type::GHOST_DISPL );
        pl.set( "function_parameters", "0.05" );
        pl.set( "leader_phase_name", "PhaseBottom" );
        pl.set( "follower_phase_name", "PhaseBottom" );
        pl.set( "leader_properties", "PropYoungsBottom,Material" );
        aParameterLists( FEM::STABILIZATION ).add_parameter_list( pl );

        /* -------------------------------------------------------------------------------------------- */
        /*                                               IWG                                            */
        /* -------------------------------------------------------------------------------------------- */

        /* -------------------------------------- Body Material Top ----------------------------------- */
        pl = prm::create_IWG_parameter_list();
        pl.set( "IWG_name", "IWGBulkStructTop" );
        pl.set( "IWG_type", (uint)fem::IWG_Type::STRUC_NON_LINEAR_BULK_PF );
        pl.set( "IWG_bulk_type", (uint)fem::Element_Type::BULK );
        pl.set( "dof_residual", "UX,UY" );
        pl.set( "leader_phase_name", "PhaseTop" );
        pl.set( "leader_dof_dependencies", "UX,UY" );
        pl.set( "leader_constitutive_models", "MaterialTop,ElastLinIso" );
        aParameterLists( FEM::IWG ).add_parameter_list( pl );

        /* ----------------------------- Body Material Bottom (Left or Whole) ------------------------- */
        pl = prm::create_IWG_parameter_list();
        pl.set( "IWG_name", "IWGBulkStructBottom" );
        pl.set( "IWG_type", (uint)fem::IWG_Type::STRUC_NON_LINEAR_BULK_PF );
        pl.set( "IWG_bulk_type", (uint)fem::Element_Type::BULK );
        pl.set( "dof_residual", "UX,UY" );
        pl.set( "leader_phase_name", "PhaseBottom" );
        pl.set( "leader_dof_dependencies", "UX,UY" );
        pl.set( "leader_constitutive_models", "MaterialBottom,ElastLinIso" );
        aParameterLists( FEM::IWG ).add_parameter_list( pl );

        /* ---------------------------------------- Dirichlet Top ------------------------------------- */
        if ( tTopHasDirichlet )
        {
            pl = prm::create_IWG_parameter_list();
            pl.set( "IWG_name", "IWGDirichletTop" );
            pl.set( "IWG_type", (uint)fem::IWG_Type::STRUC_NON_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE_PF );
            pl.set( "IWG_bulk_type", (uint)fem::Element_Type::SIDESET );
            pl.set( "dof_residual", "UX,UY" );
            pl.set( "leader_phase_name", "PhaseTop" );
            pl.set( "side_ordinals", "3" );
            pl.set( "leader_dof_dependencies", "UX,UY" );
            pl.set( "leader_properties", "PropDirichletTop,Dirichlet" );
            pl.set( "leader_constitutive_models", "MaterialTop,ElastLinIso" );
            pl.set( "stabilization_parameters", "SPDirichletNitscheTop,DirichletNitsche" );
            aParameterLists( FEM::IWG ).add_parameter_list( pl );
        }

        /* -------------------------------------- Dirichlet Bottom ------------------------------------ */
        if ( tBottomHasDirichlet )
        {
            pl = prm::create_IWG_parameter_list();
            pl.set( "IWG_name", "IWGDirichletBottom" );
            pl.set( "IWG_type", (uint)fem::IWG_Type::STRUC_NON_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE_PF );
            pl.set( "IWG_bulk_type", (uint)fem::Element_Type::SIDESET );
            pl.set( "dof_residual", "UX,UY" );
            pl.set( "leader_phase_name", "PhaseBottom" );
            pl.set( "side_ordinals", "1" );
            pl.set( "leader_dof_dependencies", "UX,UY" );
            pl.set( "leader_properties", "PropDirichletBottom,Dirichlet" );
            pl.set( "leader_constitutive_models", "MaterialBottom,ElastLinIso" );
            pl.set( "stabilization_parameters", "SPDirichletNitscheBottom,DirichletNitsche" );
            aParameterLists( FEM::IWG ).add_parameter_list( pl );
        }

        /* ---------------------------------------- Symmetry ------------------------------------------ */
        if ( tTopHasSymmetryWest )
        {
            pl = prm::create_IWG_parameter_list();
            pl.set( "IWG_name", "IWGSymmetryTopWest" );
            pl.set( "IWG_type", (uint)fem::IWG_Type::STRUC_NON_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE_PF );
            pl.set( "IWG_bulk_type", (uint)fem::Element_Type::SIDESET );
            pl.set( "dof_residual", "UX,UY" );
            pl.set( "leader_phase_name", "PhaseTop" );
            pl.set( "side_ordinals", "4" );
            pl.set( "leader_dof_dependencies", "UX,UY" );
            pl.set( "leader_properties", "PropFixed,Dirichlet;PropSelectX,Select" );
            pl.set( "leader_constitutive_models", "MaterialTop,ElastLinIso" );
            pl.set( "stabilization_parameters", "SPDirichletNitscheTop,DirichletNitsche" );
            aParameterLists( FEM::IWG ).add_parameter_list( pl );
        }

        if ( tTopHasSymmetryEast )
        {
            pl = prm::create_IWG_parameter_list();
            pl.set( "IWG_name", "IWGSymmetryTopEast" );
            pl.set( "IWG_type", (uint)fem::IWG_Type::STRUC_NON_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE_PF );
            pl.set( "IWG_bulk_type", (uint)fem::Element_Type::SIDESET );
            pl.set( "dof_residual", "UX,UY" );
            pl.set( "leader_phase_name", "PhaseTop" );
            pl.set( "side_ordinals", "2" );
            pl.set( "leader_dof_dependencies", "UX,UY" );
            pl.set( "leader_properties", "PropFixed,Dirichlet;PropSelectX,Select" );
            pl.set( "leader_constitutive_models", "MaterialTop,ElastLinIso" );
            pl.set( "stabilization_parameters", "SPDirichletNitscheTop,DirichletNitsche" );
            aParameterLists( FEM::IWG ).add_parameter_list( pl );
        }

        if ( tBottomHasSymmetryWest )
        {
            pl = prm::create_IWG_parameter_list();
            pl.set( "IWG_name", "IWGSymmetryBottomWest" );
            pl.set( "IWG_type", (uint)fem::IWG_Type::STRUC_NON_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE_PF );
            pl.set( "IWG_bulk_type", (uint)fem::Element_Type::SIDESET );
            pl.set( "dof_residual", "UX,UY" );
            pl.set( "leader_phase_name", "PhaseBottom" );
            pl.set( "side_ordinals", "4" );
            pl.set( "leader_dof_dependencies", "UX,UY" );
            pl.set( "leader_properties", "PropFixed,Dirichlet;PropSelectX,Select" );
            pl.set( "leader_constitutive_models", "MaterialBottom,ElastLinIso" );
            pl.set( "stabilization_parameters", "SPDirichletNitscheBottom,DirichletNitsche" );
            aParameterLists( FEM::IWG ).add_parameter_list( pl );
        }

        if ( tBottomHasSymmetryEast )
        {
            pl = prm::create_IWG_parameter_list();
            pl.set( "IWG_name", "IWGSymmetryBottomEast" );
            pl.set( "IWG_type", (uint)fem::IWG_Type::STRUC_NON_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE_PF );
            pl.set( "IWG_bulk_type", (uint)fem::Element_Type::SIDESET );
            pl.set( "dof_residual", "UX,UY" );
            pl.set( "leader_phase_name", "PhaseBottom" );
            pl.set( "side_ordinals", "2" );
            pl.set( "leader_dof_dependencies", "UX,UY" );
            pl.set( "leader_properties", "PropFixed,Dirichlet;PropSelectX,Select" );
            pl.set( "leader_constitutive_models", "MaterialBottom,ElastLinIso" );
            pl.set( "stabilization_parameters", "SPDirichletNitscheBottom,DirichletNitsche" );
            aParameterLists( FEM::IWG ).add_parameter_list( pl );
        }

        /* ------------------------------------- Contact Interface ------------------------------------ */
        // determine the requested contact type
        IWG_Type tContactIWGType = IWG_Type::END_IWG_TYPE;
        if ( tContactType == "mlika" )
        {
            if ( tContactBias == "symmetric" )
                tContactIWGType = IWG_Type::STRUC_NONLINEAR_CONTACT_MLIKA_UNBIASED_SYMMETRIC;
            if ( tContactBias == "unsymmetric" )
                tContactIWGType = IWG_Type::STRUC_NONLINEAR_CONTACT_MLIKA_UNBIASED_UNSYMMETRIC;
            if ( tContactBias == "neutral" )
                tContactIWGType = IWG_Type::STRUC_NONLINEAR_CONTACT_MLIKA_UNBIASED_NEUTRAL;
        }
        if ( tContactType == "small" )
        {
            if ( tContactBias == "symmetric" )
                tContactIWGType = IWG_Type::STRUC_LINEAR_CONTACT_NORMAL_SYMMETRIC_NITSCHE_UNBIASED;
            if ( tContactBias == "unsymmetric" )
                tContactIWGType = IWG_Type::STRUC_LINEAR_CONTACT_NORMAL_UNSYMMETRIC_NITSCHE_UNBIASED;
            if ( tContactBias == "neutral" )
                tContactIWGType = IWG_Type::STRUC_LINEAR_CONTACT_NORMAL_NEUTRAL_NITSCHE_UNBIASED;
        }

        pl = prm::create_IWG_parameter_list();
        pl.set( "IWG_name", "IWGContactInterfaceTop" );
        pl.set( "IWG_type", (uint)tContactIWGType );
        pl.set( "IWG_bulk_type", (uint)fem::Element_Type::NONCONFORMAL_SIDESET );
        pl.set( "dof_residual", "UX,UY" );
        pl.set( "leader_phase_name", "PhaseTop" );
        pl.set( "neighbor_phases", "PhaseVoid" );
        pl.set( "follower_phase_name", "PhaseBottom" );
        pl.set( "leader_dof_dependencies", "UX,UY" );
        pl.set( "follower_dof_dependencies", "UX,UY" );
        pl.set( "leader_constitutive_models", "MaterialTop,ElastLinIso" );
        pl.set( "follower_constitutive_models", "MaterialBottom,ElastLinIso" );
        pl.set( "stabilization_parameters", "SPContactInterfaceTop,NitscheInterface" );
        pl.set( "analytical_jacobian", false );
        aParameterLists( FEM::IWG ).add_parameter_list( pl );

        pl = prm::create_IWG_parameter_list();
        pl.set( "IWG_name", "IWGContactInterfaceBottom" );
        pl.set( "IWG_type", (uint)tContactIWGType );
        pl.set( "IWG_bulk_type", (uint)fem::Element_Type::NONCONFORMAL_SIDESET );
        pl.set( "dof_residual", "UX,UY" );
        pl.set( "leader_phase_name", "PhaseBottom" );
        pl.set( "neighbor_phases", "PhaseVoid" );
        pl.set( "follower_phase_name", "PhaseTop" );
        pl.set( "leader_dof_dependencies", "UX,UY" );
        pl.set( "follower_dof_dependencies", "UX,UY" );
        pl.set( "leader_constitutive_models", "MaterialBottom,ElastLinIso" );
        pl.set( "follower_constitutive_models", "MaterialTop,ElastLinIso" );
        pl.set( "stabilization_parameters", "SPContactInterfaceBottom,NitscheInterface" );
        pl.set( "analytical_jacobian", false );
        aParameterLists( FEM::IWG ).add_parameter_list( pl );

        if ( tUseGhost )
        {
            // ghost  displacement
            pl = prm::create_IWG_parameter_list();
            pl.set( "IWG_name", "IWGGPDisplTop" );
            pl.set( "IWG_type", (uint)fem::IWG_Type::GHOST_NORMAL_FIELD );
            pl.set( "IWG_bulk_type", (uint)fem::Element_Type::DOUBLE_SIDESET );
            pl.set( "dof_residual", "UX,UY" );
            pl.set( "leader_phase_name", "PhaseTop" );
            pl.set( "follower_phase_name", "PhaseTop" );
            pl.set( "leader_dof_dependencies", "UX,UY" );
            pl.set( "follower_dof_dependencies", "UX,UY" );
            pl.set( "stabilization_parameters", "SPGPDisplTop,GhostSP" );
            aParameterLists( FEM::IWG ).add_parameter_list( pl );

            // ghost  displacement
            pl = prm::create_IWG_parameter_list();
            pl.set( "IWG_name", "IWGGPDisplBottom" );
            pl.set( "IWG_type", (uint)fem::IWG_Type::GHOST_NORMAL_FIELD );
            pl.set( "IWG_bulk_type", (uint)fem::Element_Type::DOUBLE_SIDESET );
            pl.set( "dof_residual", "UX,UY" );
            pl.set( "leader_phase_name", "PhaseBottom" );
            pl.set( "follower_phase_name", "PhaseBottom" );
            pl.set( "leader_dof_dependencies", "UX,UY" );
            pl.set( "follower_dof_dependencies", "UX,UY" );
            pl.set( "stabilization_parameters", "SPGPDisplBottom,GhostSP" );
            aParameterLists( FEM::IWG ).add_parameter_list( pl );
        }

        /* -------------------------------------------------------------------------------------------- */
        /*                                               IQI                                            */
        /* -------------------------------------------------------------------------------------------- */

        Vector< std::pair< std::string, std::string > > tSides      = { { "Top", "Bottom" },
                 { "Bottom", "Top" } };
        Vector< std::pair< int, std::string > >         tComponents = { { 0, "X" }, { 1, "Y" } };
        for ( auto const &[ tLeaderSide, tFollowerSide ] : tSides )
        {
            for ( auto const &[ tDir, tAxis ] : tComponents )
            {
                pl = prm::create_IQI_parameter_list();
                pl.set( "IQI_name", "IQIStress" + tLeaderSide + tAxis );
                pl.set( "IQI_type", (uint)fem::IQI_Type::NORMAL_STRESS_CAUCHY );
                pl.set( "IQI_bulk_type", (uint)fem::Element_Type::BULK );
                pl.set( "leader_dof_dependencies", "UX,UY" );
                pl.set( "leader_constitutive_models", "Material" + tLeaderSide + ",ElastLinIso" );
                pl.set( "leader_phase_name", "Phase" + tLeaderSide );
                pl.set( "vectorial_field_index", tDir );
                aParameterLists( FEM::IQI ).add_parameter_list( pl );

                pl = prm::create_IQI_parameter_list();
                pl.set( "IQI_name", "IQIDisp" + tLeaderSide + tAxis );
                pl.set( "IQI_type", (uint)fem::IQI_Type::DOF );
                pl.set( "IQI_bulk_type", (uint)fem::Element_Type::BULK );
                pl.set( "dof_quantity", "UX,UY" );
                pl.set( "leader_dof_dependencies", "UX,UY" );
                pl.set( "leader_phase_name", "Phase" + tLeaderSide );
                pl.set( "vectorial_field_index", tDir );
            }

            pl = prm::create_IQI_parameter_list();
            pl.set( "IQI_name", "IQIVonMises" + tLeaderSide );
            pl.set( "IQI_type", (uint)fem::IQI_Type::VON_MISES_STRESS_CAUCHY );
            pl.set( "IQI_bulk_type", (uint)fem::Element_Type::BULK );
            pl.set( "leader_dof_dependencies", "UX,UY" );
            pl.set( "leader_constitutive_models", "Material" + tLeaderSide + ",ElastLinIso" );
            pl.set( "leader_phase_name", "Phase" + tLeaderSide );
            aParameterLists( FEM::IQI ).add_parameter_list( pl );

            /*pl = prm::create_IQI_parameter_list();
            pl.set("IQI_name", "IQIContactPressure" + tLeaderSide);
            pl.set("IQI_type", (uint)fem::IQI_Type::CONTACT_PRESSURE_CURRENT);
            pl.set("IQI_bulk_type", (uint)fem::Element_Type::NONCONFORMAL_SIDESET);
            pl.set("leader_dof_dependencies", "UX,UY");
            pl.set("leader_phase_name", "Phase" + tLeaderSide);
            pl.set("neighbor_phases", "PhaseVoid");
            pl.set("follower_phase_name", "Phase" + tFollowerSide);
            pl.set("leader_constitutive_models", "Material" + tLeaderSide + ",TractionCM");
            aParameterLists( FEM::IQI ).add_parameter_list( pl );
            */

            pl = prm::create_IQI_parameter_list();
            pl.set( "IQI_name", "IQIGap" + tLeaderSide );
            pl.set( "IQI_type", (uint)fem::IQI_Type::GAP );
            pl.set( "IQI_bulk_type", (uint)fem::Element_Type::NONCONFORMAL_SIDESET );
            pl.set( "leader_phase_name", "Phase" + tLeaderSide );
            pl.set( "neighbor_phases", "PhaseVoid" );
            pl.set( "follower_phase_name", "Phase" + tFollowerSide );
            aParameterLists( FEM::IQI ).add_parameter_list( pl );
        }

        /* -------------------------------------------------------------------------------------------- */
        /*                                   Computation Parameter List                                 */
        /* -------------------------------------------------------------------------------------------- */
        aParameterLists( FEM::COMPUTATION ).set( "is_analytical_forward", true );
        aParameterLists.set( "finite_difference_scheme_forward", (uint)( tFDScheme ) );
        aParameterLists.set( "finite_difference_perturbation_size_forward", tFDPerturbationSize );
        aParameterLists.set( "finite_difference_perturbation_strategy", (uint)tFDPerturbationStrategy );
        aParameterLists.set( "nonconformal_integration_order", static_cast< uint >( tNonconformalIntegrationOrder ) );
        aParameterLists.set( "nonconformal_max_negative_ray_length", tMaxNegativeRayLength );
        aParameterLists.set( "nonconformal_max_positive_ray_length", tMaxPositiveRayLength );
    }

    /* ---------------------------------------------------------------------------------------------- */
    /*                                           ### SOL ###                                          */
    /* ---------------------------------------------------------------------------------------------- */
    void SOLParameterList( Module_Parameter_Lists &aParameterLists )
    {
        /* -------------------------------------------------------------------------------------------- */
        /*                                        Linear Algorithm                                      */
        /* -------------------------------------------------------------------------------------------- */
        Parameter_List pl = prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL );
        // pl.set( "Solver_Type", "Amesos_Umfpack" );
        aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( pl );

        /* -------------------------------------------------------------------------------------------- */
        /*                                         Linear Solver                                        */
        /* -------------------------------------------------------------------------------------------- */
        pl = prm::create_linear_solver_parameter_list();
        aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list( pl );

        /* -------------------------------------------------------------------------------------------- */
        /*                                      Nonlinear Algorithm                                     */
        /* -------------------------------------------------------------------------------------------- */
        pl = prm::create_nonlinear_algorithm_parameter_list();
        pl.set( "NLA_combined_res_jac_assembly", true );
        pl.set( "NLA_Linear_solver", 0 );
        pl.set( "NLA_max_iter", tMaxIterations );
        pl.set( "NLA_rel_res_norm_drop", tRelResNormDrop );
        pl.set( "NLA_load_control_strategy", (uint)( sol::SolverLoadControlType::Linear ) );
        pl.set( "NLA_load_control_factor", tLoadControlFactor );
        pl.set( "NLA_load_control_steps", tLoadControlSteps );
        pl.set( "NLA_load_control_relres", tLoadControlRelRes );
        pl.set( "NLA_relaxation_parameter", tRelaxation );
        pl.set( "NLA_remap_strategy", (uint)tRaytracingStrategy );
        pl.set( "NLA_remap_load_step_frequency", tRemapLoadStepFrequency );
        pl.set( "NLA_remap_iteration_frequency", tRemapIterationFrequency );
        pl.set( "NLA_remap_residual_change_tolerance", tRemapResidualChange );
        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list( pl );

        /* -------------------------------------------------------------------------------------------- */
        /*                                       Nonlinear Solver                                       */
        /* -------------------------------------------------------------------------------------------- */
        pl = prm::create_nonlinear_solver_parameter_list();
        pl.set( "NLA_DofTypes", "UX,UY" );
        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list( pl );

        /* -------------------------------------------------------------------------------------------- */
        /*                                        Time Algorithm                                        */
        /* -------------------------------------------------------------------------------------------- */
        pl = prm::create_time_solver_algorithm_parameter_list();
        pl.set( "TSA_Num_Time_Steps", 1 );
        pl.set( "TSA_Time_Frame", 100.0 );
        aParameterLists( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list( pl );

        /* -------------------------------------------------------------------------------------------- */
        /*                                         Time Solver                                          */
        /* -------------------------------------------------------------------------------------------- */
        pl = prm::create_time_solver_parameter_list();
        pl.set( "TSA_DofTypes", "UX,UY" );
        pl.set( "TSA_Output_Indices", "0" );
        pl.set( "TSA_Output_Criteria", F2STR( Output_Criterion ) );
        aParameterLists( SOL::TIME_SOLVERS ).add_parameter_list( pl );

        /* -------------------------------------------------------------------------------------------- */
        /*                                        Solver Warehouse                                      */
        /* -------------------------------------------------------------------------------------------- */

        /* -------------------------------------------------------------------------------------------- */
        /*                                      Preconditioner List                                     */
        /* -------------------------------------------------------------------------------------------- */
        pl = prm::create_preconditioner_parameter_list( sol::PreconditionerType::NONE );
        aParameterLists( SOL::PRECONDITIONERS ).add_parameter_list( pl );
    }

    void MSIParameterList( Module_Parameter_Lists &aParameterLists )
    {
        /* -------------------------------------------------------------------------------------------- */
        /*                                       MSI Parameter List                                     */
        /* -------------------------------------------------------------------------------------------- */
        aParameterLists.set( "UX", 0 );
        aParameterLists.set( "UY", 0 );
    }

    /* ---------------------------------------------------------------------------------------------- */
    /*                                           ### VIS ###                                          */
    /* ---------------------------------------------------------------------------------------------- */
    void VISParameterList( Module_Parameter_Lists &aParameterLists )
    {
        /* -------------------------------------------------------------------------------------------- */
        /*                                          VIS Parameter                                       */
        /* -------------------------------------------------------------------------------------------- */
        aParameterLists.set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName ) );
        aParameterLists.set( "Mesh_Type", (uint)vis::VIS_Mesh_Type::STANDARD );

        // std::string tSetNames = tDomain;
        std::string tSetNames = tDomain + "," + tContactInterface;

        std::string tFieldNames = "";
        std::string tFieldTypes = "";
        std::string tIQINames   = "";

        // add stress components as STRESSX, STRESSY
        Vector< std::string > tSides      = { "Top", "Bottom" };
        Vector< std::string > tComponents = { "X", "Y" };

        for ( auto const &tSide : tSides )
        {
            for ( auto const &tComponent : tComponents )
            {
                tFieldNames += ",STRESSBulk" + tSide + tComponent;
                tFieldTypes += ",ELEMENTAL_AVG";    // ELEMENTAL_AVG
                tIQINames += ",IQIStress" + tSide + tComponent;

                // add the displacement components
                tFieldNames += ",U" + tSide + tComponent;
                tFieldTypes += ",NODAL";    // ELEMENTAL_AVG
                tIQINames += ",IQIDisp" + tSide + tComponent;
            }

            tFieldNames += ",STRESSVonMises" + tSide;
            tFieldTypes += ",ELEMENTAL_AVG";    // ELEMENTAL_AVG
            tIQINames += ",IQIVonMises" + tSide;

            // add the traction components
            /*tFieldNames += ",ContactPressure" + tSide;
            tFieldTypes += ",NODAL";
            tIQINames += ",IQIContactPressure" + tSide;
            */

            // add the gap
            tFieldNames += ",Gap" + tSide;
            tFieldTypes += ",NODAL";
            tIQINames += ",IQIGap" + tSide;
        }

        // remove first comma
        tFieldNames = tFieldNames.substr( 1 );
        tFieldTypes = tFieldTypes.substr( 1 );
        tIQINames   = tIQINames.substr( 1 );

        aParameterLists.set( "Set_Names", tSetNames );
        aParameterLists.set( "Field_Names", tFieldNames );
        aParameterLists.set( "Field_Type", tFieldTypes );
        aParameterLists.set( "IQI_Names", tIQINames );
        aParameterLists.set( "Save_Frequency", 1 );
        aParameterLists.set( "Time_Offset", 0.0 );
    }

    void MORISGENERALParameterList( Module_Parameter_Lists &aParameterLists ) {}
}    // namespace moris
#ifdef __cplusplus
}
#endif
