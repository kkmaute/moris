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
#include "fn_PRM_FEM_Parameters.hpp"
#include "fn_PRM_GEN_Parameters.hpp"
#include "fn_PRM_HMR_Parameters.hpp"
#include "fn_PRM_MSI_Parameters.hpp"
#include "fn_PRM_OPT_Parameters.hpp"
#include "fn_PRM_SOL_Parameters.hpp"
#include "fn_PRM_VIS_Parameters.hpp"
#include "fn_PRM_XTK_Parameters.hpp"
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
    moris::real tHeight = 1.0;
    moris::real tWidth  = 1.0;

    uint tNumElems  = 10;
    uint tNumElemsX = tNumElems;
    uint tNumElemsY = tNumElems;

    moris::real tInitialGap           = 0.001;
    moris::real tInitialGapMiddle     = 0.503 * tHeight;
    moris::real tPenetration          = 0.05;
    moris::real tVerticalDisplacement = -tInitialGap - tPenetration;

    /* --------------------------------------- Top Geometry ----------------------------------------- */
    bool tTopIsPlane         = false;
    bool tTopHasSymmetryWest = false;
    bool tTopHasSymmetryEast = false;

    moris::real tInterfaceTopAngle = 0 * M_PI / 180.0;
    moris::real tTopYShift         = tInitialGapMiddle + tInitialGap / 2;
    moris::real tTopXShift         = 0.5;
    moris::real tTopParabolaFactor = 1;

    /* --------------------------------------- Bottom Geometry -------------------------------------- */
    bool tBottomIsPlane         = true;
    bool tBottomHasSymmetryWest = false;
    bool tBottomHasSymmetryEast = false;

    moris::real tInterfaceBottomAngle = 0 * M_PI / 180.0;
    moris::real tBottomYShift         = tInitialGapMiddle - tInitialGap / 2;
    moris::real tBottomXShift         = 0.5;
    moris::real tBottomParabolaFactor = -0.7;

    /* ------------------------------------- Material Parameters ------------------------------------ */
    real tTopEmod    = 1000.0;
    real tBottomEmod = 1000.0;
    real tTopPois    = 0.0;
    real tBottomPois = 0.0;

    /* --------------------------------------- Domain Setting --------------------------------------- */
    std::string tDomainOffset           = "0.0,0.0";    //
    std::string tInterpolationOrder     = std::to_string( gInterpolationOrder );
    std::string tDomainTop              = "HMR_dummy_n_p1,HMR_dummy_c_p1";
    std::string tDomainBottom           = "HMR_dummy_n_p2,HMR_dummy_c_p2";
    std::string tContactInterfaceTop    = "ncss|iside_b0_1_b1_0|iside_b0_2_b1_0";
    std::string tContactInterfaceBottom = "ncss|iside_b0_2_b1_0|iside_b0_1_b1_0";
    std::string tContactInterface       = tContactInterfaceBottom + "," + tContactInterfaceTop;
    std::string tDomain                 = tDomainTop + "," + tDomainBottom;

    /* ------------------------------------------ File I/O ------------------------------------------ */
    std::string tOutputFileName =
            "Parabolic_Indenter_Case_" + std::to_string( gCaseIndex );

    /* ------------------------------------------- Solver ------------------------------------------- */
    int                    tMaxIterations          = 200;
    real                   tRelResNormDrop         = 1e-8;
    real                   tRelaxation             = 1.0;
    real                   tFDPerturbationSize     = 1e-8;
    fem::Perturbation_Type tFDPerturbationStrategy = fem::Perturbation_Type::ABSOLUTE;
    fem::FDScheme_Type     tFDScheme               = fem::FDScheme_Type::POINT_3_CENTRAL;

    int  tLoadControlSteps    = 15;
    real tLoadControlFactor   = 0.1;
    real tLoadControlRelRes   = 1e-2;
    real tLoadControlExponent = 1.0;

    real tRemapResidualChange     = -1e-2;
    int  tRemapLoadStepFrequency  = 1;
    int  tRemapIterationFrequency = 10;
    auto tRaytracingStrategy      = sol::SolverRaytracingStrategy::MixedNthLoadStepAndResidualChange;
    real tMaxNegativeRayLength    = -0.005;
    real tMaxPositiveRayLength    = 0.01;

    mtk::Integration_Order tNonconformalIntegrationOrder = mtk::Integration_Order::BAR_6;

    /* ------------------------------------------- Contact ------------------------------------------ */
    std::string tContactType          = "small";
    std::string tContactBias          = "unsymmetric";
    real        tContactStabilization = 10.0;

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

    moris::real tMinimumLevelSetValue = 1.0e-8;

    // Level set function for a parabola that splits the domain into a top and a bottom part.
    moris::real Parabola( const moris::Matrix< DDRMat > &aCoordinates,
            const moris::Vector< moris::real >          &aGeometryParameters )
    {

        moris::real tXShift = aGeometryParameters( 0 );
        moris::real tYShift = aGeometryParameters( 1 );
        moris::real tFactor = aGeometryParameters( 2 );

        moris::real x = aCoordinates( 0 ) - tXShift;
        moris::real y = aCoordinates( 1 ) - tYShift;

        // Compute Signed-Distance field
        moris::real tVal = y - ( tFactor * x * x );

        return ( std::abs( tVal ) < tMinimumLevelSetValue ) ? tMinimumLevelSetValue : tVal;
    }

    /* ----------------------------------- Property Field Function ---------------------------------- */

    void Func_Const( moris::Matrix< moris::DDRMat >         &aPropMatrix,
            moris::Vector< moris::Matrix< moris::DDRMat > > &aParameters,
            moris::fem::Field_Interpolator_Manager          *aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    void Func_Dirichlet( moris::Matrix< moris::DDRMat >     &aPropMatrix,
            moris::Vector< moris::Matrix< moris::DDRMat > > &aParameters,
            moris::fem::Field_Interpolator_Manager          *aFIManager )
    {
        real tLoadFactor = gLogger.get_action_data( "NonLinearAlgorithm", "Newton", "Solve", "LoadFactor" );
        aPropMatrix      = aParameters( 0 ) * tLoadFactor;
    }

    void Func_Linear_Along_X( moris::Matrix< moris::DDRMat > &aPropMatrix,
            moris::Vector< moris::Matrix< moris::DDRMat > >  &aParameters,
            moris::fem::Field_Interpolator_Manager           *aFIManager )
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

    void Func_Select_X( moris::Matrix< moris::DDRMat > &aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >   &aParameters,
            moris::fem::Field_Interpolator_Manager     *aFIManager )
    {
        aPropMatrix.set_size( 2, 2, 0.0 );
        aPropMatrix( 0, 0 ) = 1.0;
    }

    void Func_Body_Load( moris::Matrix< moris::DDRMat >     &aPropMatrix,
            moris::Vector< moris::Matrix< moris::DDRMat > > &aParameters,
            moris::fem::Field_Interpolator_Manager          *aFIManager )
    {

        auto tXp    = aFIManager->get_IG_geometry_interpolator()->valx();
        auto param  = aParameters( 0 )( 0 );
        aPropMatrix = {
            { 0.0 },
            { param },
        };
    }

    bool Output_Criterion( moris::tsa::Time_Solver *aTimeSolver ) { return true; }

    /* ---------------------------------------------------------------------------------------------- */
    /*                                           ### OPT ###                                          */
    /* ---------------------------------------------------------------------------------------------- */
    void OPTParameterList( moris::Vector< moris::Submodule_Parameter_Lists > &tParameterlist )
    {
        tParameterlist.resize( 1 );
        Parameter_List pl = prm::create_opt_problem_parameter_list();
        pl.set( "is_optimization_problem", false );
        tParameterlist( 0 ).add_parameter_list( pl );
    }

    /* ---------------------------------------------------------------------------------------------- */
    /*                                           ### HMR ###                                          */
    /* ---------------------------------------------------------------------------------------------- */
    void HMRParameterList( moris::Vector< moris::Submodule_Parameter_Lists > &tParameterlist )
    {
        tParameterlist.resize( 1 );
        auto pl = prm::create_hmr_parameter_list();

        /* --------------------------------------- Domain Settings ------------------------------------ */
        pl.set( "number_of_elements_per_dimension",
                std::to_string( tNumElemsX ) + "," + std::to_string( tNumElemsY ) );
        pl.set( "domain_dimensions", std::to_string( tWidth ) + "," + std::to_string( tHeight ) );
        pl.set( "domain_offset", tDomainOffset );
        pl.set( "domain_sidesets", "1,2,3,4" );

        /* ---------------------------------------- Lagrange Mesh ------------------------------------- */
        pl.set( "lagrange_pattern", "0" );
        pl.set( "lagrange_orders", tInterpolationOrder );
        pl.set( "lagrange_output_meshes", "0" );

        /* -------------------------------------- B-Spline Mesh --------------------------------------- */
        pl.set( "bspline_pattern", "0" );
        pl.set( "bspline_orders", tInterpolationOrder );
        pl.set( "lagrange_to_bspline", "0" );
        pl.set( "truncate_bsplines", 1 );

        /* ----------------------------------------- Refinement --------------------------------------- */
        pl.set( "refinement_buffer", 0 );
        pl.set( "staircase_buffer", 0 );
        pl.set( "initial_refinement", "0" );
        pl.set( "initial_refinement_pattern", "0" );

        /* ---------------------------------------- Miscellaneous ------------------------------------- */
        pl.set( "use_number_aura", 1 );
        pl.set( "use_multigrid", 0 );
        pl.set( "severity_level", 0 );
        tParameterlist( 0 ).add_parameter_list( pl );
    }

    /* ---------------------------------------------------------------------------------------------- */
    /*                                           ### XTK ###                                          */
    /* ---------------------------------------------------------------------------------------------- */
    void XTKParameterList( moris::Vector< moris::Submodule_Parameter_Lists > &tParameterlist )
    {
        tParameterlist.resize( 1 );
        auto pl = prm::create_xtk_parameter_list();
        pl.set( "decompose", true );
        pl.set( "decomposition_type", "conformal" );
        pl.set( "enrich", true );
        pl.set( "basis_rank", "bspline" );
        pl.set( "enrich_mesh_indices", "0" );
        pl.set( "ghost_stab", tUseGhost );
        pl.set( "multigrid", false );
        pl.set( "verbose", true );
        pl.set( "print_enriched_ig_mesh", true );
        pl.set( "exodus_output_XTK_ig_mesh", true );
        pl.set( "high_to_low_dbl_side_sets", true );
        pl.set( "only_generate_xtk_temp", tOnlyGenerateMesh );
        tParameterlist( 0 ).add_parameter_list( pl );
    }

    /* ---------------------------------------------------------------------------------------------- */
    /*                                           ### GEN ###                                          */
    /* ---------------------------------------------------------------------------------------------- */
    void GENParameterList( moris::Vector< moris::Submodule_Parameter_Lists > &tParameterlist )
    {
        tParameterlist.resize( 3 );
        /* -------------------------------------------------------------------------------------------- */
        /*                                        GEN Parameter                                         */
        /* -------------------------------------------------------------------------------------------- */
        Parameter_List pl = prm::create_gen_parameter_list();
        pl.set( "number_of_phases", 3 );
        pl.set( "phase_function_name", F2STR( Phase_Index_Split ) );
        tParameterlist( 0 ).add_parameter_list( pl );

        /* -------------------------------------------------------------------------------------------- */
        /*                                      Geometry Parameter                                      */
        /* -------------------------------------------------------------------------------------------- */

        /* ----------------------------------------- Top Plane ---------------------------------------- */
        if ( tTopIsPlane )
        {
            // calculate the components of the rotated normal
            real tNormalX = std::cos( tInterfaceTopAngle + M_PI_2 );
            real tNormalY = std::sin( tInterfaceTopAngle + M_PI_2 );

            pl = prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE );
            pl.set( "plx", tTopXShift );
            pl.set( "center_y", tTopYShift );
            pl.set( "normal_x", tNormalX );
            pl.set( "normal_y", tNormalY );

            pl.set( "number_of_refinements", 0 );
            pl.set( "refinement_mesh_index", 0 );
            tParameterlist( 1 ).add_parameter_list( pl );
        }
        else
        {
            // use a parabola for the top side
            pl = prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED );
            pl.set( "field_function_name", "Parabola" );
            pl.insert( "variable_1", tTopXShift );
            pl.insert( "variable_2", tTopYShift );
            pl.insert( "variable_3", tTopParabolaFactor );
            tParameterlist( 1 ).add_parameter_list( pl );
        }

        /* --------------------------------------- Bottom Plane --------------------------------------- */
        if ( tBottomIsPlane )
        {
            // calculate the components of the rotated normal
            real tNormalX = std::cos( tInterfaceBottomAngle + M_PI_2 );
            real tNormalY = std::sin( tInterfaceBottomAngle + M_PI_2 );

            pl = prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE );
            pl.set( "center_x", tBottomXShift );
            pl.set( "center_y", tBottomYShift );
            pl.set( "normal_x", tNormalX );
            pl.set( "normal_y", tNormalY );

            pl.set( "number_of_refinements", 0 );
            pl.set( "refinement_mesh_index", 0 );
            tParameterlist( 1 ).add_parameter_list( pl );
        }
        else
        {
            // use a parabola for the bottom side
            pl = prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED );
            pl.set( "field_function_name", "Parabola" );
            pl.insert( "variable_1", tBottomXShift );
            pl.insert( "variable_2", tBottomYShift );
            pl.insert( "variable_3", tBottomParabolaFactor );
            tParameterlist( 1 ).add_parameter_list( pl );
        }
    }

    /* ---------------------------------------------------------------------------------------------- */
    /*                                           ### FEM ###                                          */
    /* ---------------------------------------------------------------------------------------------- */
    void FEMParameterList( moris::Vector< moris::Submodule_Parameter_Lists > &tParameterList )
    {
        tParameterList.resize( 9 );
        uint tPropIndex  = 0;
        uint tCMIndex    = 1;
        uint tSPIndex    = 2;
        uint tIWGIndex   = 3;
        uint tIQIIndex   = 4;
        uint tFEMIndex   = 5;
        uint tPhaseIndex = 7;

        /* -------------------------------------------------------------------------------------------- */
        /*                                            Phases                                            */
        /* -------------------------------------------------------------------------------------------- */
        Parameter_List pl = prm::create_phase_parameter_list();
        pl.set( "phase_name", "PhaseVoid" );
        pl.set( "phase_indices", "0" );
        tParameterList( tPhaseIndex ).add_parameter_list( pl );

        pl = prm::create_phase_parameter_list();
        pl.set( "phase_name", "PhaseTop" );
        pl.set( "phase_indices", "1" );
        tParameterList( tPhaseIndex ).add_parameter_list( pl );

        pl = prm::create_phase_parameter_list();
        pl.set( "phase_name", "PhaseBottom" );
        pl.set( "phase_indices", "2" );
        tParameterList( tPhaseIndex ).add_parameter_list( pl );

        pl = prm::create_phase_parameter_list();
        pl.set( "phase_name", "PhaseAll" );
        pl.set( "phase_indices", "1,2" );
        tParameterList( tPhaseIndex ).add_parameter_list( pl );

        /* -------------------------------------------------------------------------------------------- */
        /*                                           Properties                                         */
        /* -------------------------------------------------------------------------------------------- */

        /* -------------------------------------- Youngs Modulus -------------------------------------- */
        pl = prm::create_property_parameter_list();
        pl.set( "property_name", "PropYoungsTop" );
        pl.set( "function_parameters", std::to_string( tTopEmod ) );
        pl.set( "value_function", F2STR( Func_Const ) );
        tParameterList( tPropIndex ).add_parameter_list( pl );

        pl = prm::create_property_parameter_list();
        pl.set( "property_name", "PropYoungsBottom" );
        pl.set( "function_parameters", std::to_string( tBottomEmod ) );
        pl.set( "value_function", F2STR( Func_Const ) );
        tParameterList( tPropIndex ).add_parameter_list( pl );

        /* ------------------------------------------ Poisson ----------------------------------------- */
        pl = prm::create_property_parameter_list();
        pl.set( "property_name", "PropPoissonTop" );
        pl.set( "function_parameters", std::to_string( tTopPois ) );
        pl.set( "value_function", F2STR( Func_Const ) );
        tParameterList( tPropIndex ).add_parameter_list( pl );

        pl = prm::create_property_parameter_list();
        pl.set( "property_name", "PropPoissonBottom" );
        pl.set( "function_parameters", std::to_string( tBottomPois ) );
        pl.set( "value_function", F2STR( Func_Const ) );
        tParameterList( tPropIndex ).add_parameter_list( pl );

        /* --------------------------------------- Dirichlet Top -------------------------------------- */
        if ( tTopHasDirichlet )
        {
            pl = prm::create_property_parameter_list();
            pl.set( "property_name", "PropDirichletTop" );
            pl.set( "function_parameters", "0.0;" + std::to_string( tVerticalDisplacement ) );
            pl.set( "value_function", F2STR( Func_Dirichlet ) );
            tParameterList( tPropIndex ).add_parameter_list( pl );
        }

        /* ------------------------------------- Dirichlet Bottom ------------------------------------- */
        if ( tBottomHasDirichlet )
        {
            pl = prm::create_property_parameter_list();
            pl.set( "property_name", "PropDirichletBottom" );
            pl.set( "function_parameters", "0.0;0.0" );
            pl.set( "value_function", F2STR( Func_Const ) );
            tParameterList( tPropIndex ).add_parameter_list( pl );
        }

        pl = prm::create_property_parameter_list();
        pl.set( "property_name", "PropFixed" );
        pl.set( "function_parameters", "0.0;0.0" );
        pl.set( "value_function", F2STR( Func_Const ) );
        tParameterList( tPropIndex ).add_parameter_list( pl );

        /* -------------------------------------- Symmetry Plane -------------------------------------- */
        pl = prm::create_property_parameter_list();
        pl.set( "property_name", "PropSelectX" );
        pl.set( "function_parameters", "0.0;0.0" );
        pl.set( "value_function", F2STR( Func_Select_X ) );
        tParameterList( tPropIndex ).add_parameter_list( pl );

        /* -------------------------------------------------------------------------------------------- */
        /*                                       Constitutive Model                                     */
        /* -------------------------------------------------------------------------------------------- */

        /* ----------------------------------- Elastic Material Top ----------------------------------- */
        pl = prm::create_constitutive_model_parameter_list();
        pl.set( "constitutive_name", "MaterialTop" );
        pl.set( "phase_name", "PhaseTop" );
        pl.set( "constitutive_type", (uint)fem::Constitutive_Type::STRUC_LIN_ISO );
        pl.set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
        pl.set( "properties",
                "PropYoungsTop,YoungsModulus;"
                "PropPoissonTop,PoissonRatio" );
        tParameterList( tCMIndex ).add_parameter_list( pl );

        /* -------------------------------- Elastic Material Bottom ----------------------------------- */
        pl = prm::create_constitutive_model_parameter_list();
        pl.set( "constitutive_name", "MaterialBottom" );
        pl.set( "phase_name", "PhaseBottom" );
        pl.set( "constitutive_type", (uint)fem::Constitutive_Type::STRUC_LIN_ISO );
        pl.set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
        pl.set( "properties",
                "PropYoungsBottom,YoungsModulus;"
                "PropPoissonBottom,PoissonRatio" );
        tParameterList( tCMIndex ).add_parameter_list( pl );

        /* -------------------------------------------------------------------------------------------- */
        /*                                     Stabilization Parameter                                  */
        /* -------------------------------------------------------------------------------------------- */

        pl = prm::create_stabilization_parameter_parameter_list();
        pl.set( "stabilization_name", "SPDirichletNitscheTop" );
        pl.set( "stabilization_type", (uint)fem::Stabilization_Type::DIRICHLET_NITSCHE );
        pl.set( "function_parameters", "100.0" );
        pl.set( "leader_phase_name", "PhaseTop" );
        pl.set( "leader_properties", "PropYoungsTop,Material" );
        tParameterList( tSPIndex ).add_parameter_list( pl );

        pl = prm::create_stabilization_parameter_parameter_list();
        pl.set( "stabilization_name", "SPDirichletNitscheBottom" );
        pl.set( "stabilization_type", (uint)fem::Stabilization_Type::DIRICHLET_NITSCHE );
        pl.set( "function_parameters", "100.0" );
        pl.set( "leader_phase_name", "PhaseBottom" );
        pl.set( "leader_properties", "PropYoungsBottom,Material" );
        tParameterList( tSPIndex ).add_parameter_list( pl );

        pl = prm::create_stabilization_parameter_parameter_list();
        pl.set( "stabilization_name", "SPContactInterfaceTop" );
        pl.set( "stabilization_type", (uint)fem::Stabilization_Type::DIRICHLET_NITSCHE );
        pl.set( "function_parameters", std::to_string( tContactStabilization ) + ";1.0" );
        pl.set( "leader_phase_name", "PhaseTop" );
        pl.set( "leader_properties", "PropYoungsTop,Material" );
        tParameterList( tSPIndex ).add_parameter_list( pl );

        pl = prm::create_stabilization_parameter_parameter_list();
        pl.set( "stabilization_name", "SPContactInterfaceBottom" );
        pl.set( "stabilization_type", (uint)fem::Stabilization_Type::DIRICHLET_NITSCHE );
        pl.set( "function_parameters", std::to_string( tContactStabilization ) + ";1.0" );
        pl.set( "leader_phase_name", "PhaseBottom" );
        pl.set( "leader_properties", "PropYoungsBottom,Material" );
        tParameterList( tSPIndex ).add_parameter_list( pl );

        // create ghost penalty displacement Matrix
        pl = prm::create_stabilization_parameter_parameter_list();
        pl.set( "stabilization_name", "SPGPDisplTop" );
        pl.set( "stabilization_type", (uint)fem::Stabilization_Type::GHOST_DISPL );
        pl.set( "function_parameters", "0.05" );
        pl.set( "leader_phase_name", "PhaseTop" );
        pl.set( "follower_phase_name", "PhaseTop" );
        pl.set( "leader_properties", "PropYoungsTop,Material" );
        tParameterList( tSPIndex ).add_parameter_list( pl );

        // create ghost penalty  displacement Fiber
        pl = prm::create_stabilization_parameter_parameter_list();
        pl.set( "stabilization_name", "SPGPDisplBottom" );
        pl.set( "stabilization_type", (uint)fem::Stabilization_Type::GHOST_DISPL );
        pl.set( "function_parameters", "0.05" );
        pl.set( "leader_phase_name", "PhaseBottom" );
        pl.set( "follower_phase_name", "PhaseBottom" );
        pl.set( "leader_properties", "PropYoungsBottom,Material" );
        tParameterList( tSPIndex ).add_parameter_list( pl );

        /* -------------------------------------------------------------------------------------------- */
        /*                                               IWG                                            */
        /* -------------------------------------------------------------------------------------------- */

        /* -------------------------------------- Body Material Top ----------------------------------- */
        pl = prm::create_IWG_parameter_list();
        pl.set( "IWG_name", "IWGBulkStructTop" );
        pl.set( "IWG_type", (uint)fem::IWG_Type::STRUC_LINEAR_BULK );
        pl.set( "IWG_bulk_type", (uint)fem::Element_Type::BULK );
        pl.set( "dof_residual", "UX,UY" );
        pl.set( "leader_phase_name", "PhaseTop" );
        pl.set( "leader_dof_dependencies", "UX,UY" );
        pl.set( "leader_constitutive_models", "MaterialTop,ElastLinIso" );
        tParameterList( tIWGIndex ).add_parameter_list( pl );

        /* ----------------------------- Body Material Bottom (Left or Whole) ------------------------- */
        pl = prm::create_IWG_parameter_list();
        pl.set( "IWG_name", "IWGBulkStructBottom" );
        pl.set( "IWG_type", (uint)fem::IWG_Type::STRUC_LINEAR_BULK );
        pl.set( "IWG_bulk_type", (uint)fem::Element_Type::BULK );
        pl.set( "dof_residual", "UX,UY" );
        pl.set( "leader_phase_name", "PhaseBottom" );
        pl.set( "leader_dof_dependencies", "UX,UY" );
        pl.set( "leader_constitutive_models", "MaterialBottom,ElastLinIso" );
        tParameterList( tIWGIndex ).add_parameter_list( pl );

        /* ---------------------------------------- Dirichlet Top ------------------------------------- */
        if ( tTopHasDirichlet )
        {
            pl = prm::create_IWG_parameter_list();
            pl.set( "IWG_name", "IWGDirichletTop" );
            pl.set( "IWG_type", (uint)fem::IWG_Type::STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE );
            pl.set( "IWG_bulk_type", (uint)fem::Element_Type::SIDESET );
            pl.set( "dof_residual", "UX,UY" );
            pl.set( "leader_phase_name", "PhaseTop" );
            pl.set( "side_ordinals", "3" );
            pl.set( "leader_dof_dependencies", "UX,UY" );
            pl.set( "leader_properties", "PropDirichletTop,Dirichlet" );
            pl.set( "leader_constitutive_models", "MaterialTop,ElastLinIso" );
            pl.set( "stabilization_parameters", "SPDirichletNitscheTop,DirichletNitsche" );
            tParameterList( tIWGIndex ).add_parameter_list( pl );
        }

        /* -------------------------------------- Dirichlet Bottom ------------------------------------ */
        if ( tBottomHasDirichlet )
        {
            pl = prm::create_IWG_parameter_list();
            pl.set( "IWG_name", "IWGDirichletBottom" );
            pl.set( "IWG_type", (uint)fem::IWG_Type::STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE );
            pl.set( "IWG_bulk_type", (uint)fem::Element_Type::SIDESET );
            pl.set( "dof_residual", "UX,UY" );
            pl.set( "leader_phase_name", "PhaseBottom" );
            pl.set( "side_ordinals", "1" );
            pl.set( "leader_dof_dependencies", "UX,UY" );
            pl.set( "leader_properties", "PropDirichletBottom,Dirichlet" );
            pl.set( "leader_constitutive_models", "MaterialBottom,ElastLinIso" );
            pl.set( "stabilization_parameters", "SPDirichletNitscheBottom,DirichletNitsche" );
            tParameterList( tIWGIndex ).add_parameter_list( pl );
        }

        /* ---------------------------------------- Symmetry ------------------------------------------ */
        if ( tTopHasSymmetryWest )
        {
            pl = prm::create_IWG_parameter_list();
            pl.set( "IWG_name", "IWGSymmetryTopWest" );
            pl.set( "IWG_type", (uint)fem::IWG_Type::STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE );
            pl.set( "IWG_bulk_type", (uint)fem::Element_Type::SIDESET );
            pl.set( "dof_residual", "UX,UY" );
            pl.set( "leader_phase_name", "PhaseTop" );
            pl.set( "side_ordinals", "4" );
            pl.set( "leader_dof_dependencies", "UX,UY" );
            pl.set( "leader_properties", "PropFixed,Dirichlet;PropSelectX,Select" );
            pl.set( "leader_constitutive_models", "MaterialTop,ElastLinIso" );
            pl.set( "stabilization_parameters", "SPDirichletNitscheTop,DirichletNitsche" );
            tParameterList( tIWGIndex ).add_parameter_list( pl );
        }

        if ( tTopHasSymmetryEast )
        {
            pl = prm::create_IWG_parameter_list();
            pl.set( "IWG_name", "IWGSymmetryTopEast" );
            pl.set( "IWG_type", (uint)fem::IWG_Type::STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE );
            pl.set( "IWG_bulk_type", (uint)fem::Element_Type::SIDESET );
            pl.set( "dof_residual", "UX,UY" );
            pl.set( "leader_phase_name", "PhaseTop" );
            pl.set( "side_ordinals", "2" );
            pl.set( "leader_dof_dependencies", "UX,UY" );
            pl.set( "leader_properties", "PropFixed,Dirichlet;PropSelectX,Select" );
            pl.set( "leader_constitutive_models", "MaterialTop,ElastLinIso" );
            pl.set( "stabilization_parameters", "SPDirichletNitscheTop,DirichletNitsche" );
            tParameterList( tIWGIndex ).add_parameter_list( pl );
        }

        if ( tBottomHasSymmetryWest )
        {
            pl = prm::create_IWG_parameter_list();
            pl.set( "IWG_name", "IWGSymmetryBottomWest" );
            pl.set( "IWG_type", (uint)fem::IWG_Type::STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE );
            pl.set( "IWG_bulk_type", (uint)fem::Element_Type::SIDESET );
            pl.set( "dof_residual", "UX,UY" );
            pl.set( "leader_phase_name", "PhaseBottom" );
            pl.set( "side_ordinals", "4" );
            pl.set( "leader_dof_dependencies", "UX,UY" );
            pl.set( "leader_properties", "PropFixed,Dirichlet;PropSelectX,Select" );
            pl.set( "leader_constitutive_models", "MaterialBottom,ElastLinIso" );
            pl.set( "stabilization_parameters", "SPDirichletNitscheBottom,DirichletNitsche" );
            tParameterList( tIWGIndex ).add_parameter_list( pl );
        }

        if ( tBottomHasSymmetryEast )
        {
            pl = prm::create_IWG_parameter_list();
            pl.set( "IWG_name", "IWGSymmetryBottomEast" );
            pl.set( "IWG_type", (uint)fem::IWG_Type::STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE );
            pl.set( "IWG_bulk_type", (uint)fem::Element_Type::SIDESET );
            pl.set( "dof_residual", "UX,UY" );
            pl.set( "leader_phase_name", "PhaseBottom" );
            pl.set( "side_ordinals", "2" );
            pl.set( "leader_dof_dependencies", "UX,UY" );
            pl.set( "leader_properties", "PropFixed,Dirichlet;PropSelectX,Select" );
            pl.set( "leader_constitutive_models", "MaterialBottom,ElastLinIso" );
            pl.set( "stabilization_parameters", "SPDirichletNitscheBottom,DirichletNitsche" );
            tParameterList( tIWGIndex ).add_parameter_list( pl );
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
        tParameterList( tIWGIndex ).add_parameter_list( pl );

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
        tParameterList( tIWGIndex ).add_parameter_list( pl );

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
            tParameterList( tIWGIndex ).add_parameter_list( pl );

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
            tParameterList( tIWGIndex ).add_parameter_list( pl );
        }

        /* -------------------------------------------------------------------------------------------- */
        /*                                               IQI                                            */
        /* -------------------------------------------------------------------------------------------- */

        pl = prm::create_IQI_parameter_list();
        pl.set( "IQI_name", "IQIDispX" );
        pl.set( "IQI_type", (uint)fem::IQI_Type::DOF );
        pl.set( "dof_quantity", "UX,UY" );
        pl.set( "leader_dof_dependencies", "UX,UY" );
        pl.set( "leader_phase_name", "PhaseAll" );
        pl.set( "vectorial_field_index", 0 );
        tParameterList( tIQIIndex ).add_parameter_list( pl );

        pl = prm::create_IQI_parameter_list();
        pl.set( "IQI_name", "IQIDispY" );
        pl.set( "IQI_type", (uint)fem::IQI_Type::DOF );
        pl.set( "dof_quantity", "UX,UY" );
        pl.set( "leader_dof_dependencies", "UX,UY" );
        pl.set( "leader_phase_name", "PhaseAll" );
        pl.set( "vectorial_field_index", 1 );
        tParameterList( tIQIIndex ).add_parameter_list( pl );

        moris::Vector< std::string >                   tSides      = { "Top", "Bottom" };
        moris::Vector< std::pair< int, std::string > > tComponents = { { 0, "X" }, { 1, "Y" } };
        for ( auto const &tSide : tSides )
        {
            for ( auto const &[ tDir, tAxis ] : tComponents )
            {
                pl = prm::create_IQI_parameter_list();
                pl.set( "IQI_name", "IQIStress" + tSide + tAxis );
                pl.set( "IQI_type", (uint)fem::IQI_Type::NORMAL_STRESS );
                pl.set( "leader_dof_dependencies", "UX,UY" );
                pl.set( "leader_constitutive_models", "Material" + tSide + ",ElastLinIso" );
                pl.set( "leader_phase_name", "Phase" + tSide );
                pl.set( "vectorial_field_index", tDir );
                tParameterList( tIQIIndex ).add_parameter_list( pl );
            }

            pl = prm::create_IQI_parameter_list();
            pl.set( "IQI_name", "IQIContactPressure" + tSide );
            pl.set( "IQI_type", (uint)fem::IQI_Type::CONTACT_PRESSURE_REFERENCE );
            pl.set( "leader_dof_dependencies", "UX,UY" );
            pl.set( "leader_phase_name", "Phase" + tSide );
            pl.set( "leader_constitutive_models", "Material" + tSide + ",TractionCM" );
            tParameterList( tIQIIndex ).add_parameter_list( pl );
        }

        // /* ---------------------------------------- Ray Length ----------------------------------------
        pl = prm::create_IQI_parameter_list();
        pl.set( "IQI_name", "IQIGap" );
        pl.set( "IQI_type", (uint)fem::IQI_Type::GAP );
        pl.set( "leader_phase_name", "PhaseBottom" );
        pl.set( "follower_phase_name", "PhaseTop" );
        tParameterList( tIQIIndex ).add_parameter_list( pl );

        /* -------------------------------------------------------------------------------------------- */
        /*                                   Computation Parameter List                                 */
        /* -------------------------------------------------------------------------------------------- */
        pl = prm::create_computation_parameter_list();
        pl.set( "is_analytical_forward", true );
        pl.set( "finite_difference_scheme_forward", (uint)( tFDScheme ) );
        pl.set( "finite_difference_perturbation_size_forward", tFDPerturbationSize );
        pl.set( "finite_difference_perturbation_strategy", (uint)tFDPerturbationStrategy );
        pl.set( "nonconformal_integration_order", static_cast< uint >( tNonconformalIntegrationOrder ) );
        pl.set( "nonconformal_max_negative_ray_length", tMaxNegativeRayLength );
        pl.set( "nonconformal_max_positive_ray_length", tMaxPositiveRayLength );
        tParameterList( tFEMIndex ).add_parameter_list( pl );
    }

    /* ----------------------------------------------------------------------------------------------
     */
    /*                                           ### SOL ### */
    /* ----------------------------------------------------------------------------------------------
     */
    void SOLParameterList( moris::Vector< moris::Submodule_Parameter_Lists > &tParameterlist )
    {
        tParameterlist.resize( 8 );

        /* --------------------------------------------------------------------------------------------
         */
        /*                                        Linear Algorithm */
        /* --------------------------------------------------------------------------------------------
         */
        Parameter_List pl = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL );
        pl.set( "Solver_Type", "Amesos_Umfpack" );
        tParameterlist( 0 ).add_parameter_list( pl );

        /* --------------------------------------------------------------------------------------------
         */
        /*                                         Linear Solver */
        /* --------------------------------------------------------------------------------------------
         */
        pl = moris::prm::create_linear_solver_parameter_list();
        tParameterlist( 1 ).add_parameter_list( pl );

        /* --------------------------------------------------------------------------------------------
         */
        /*                                      Nonlinear Algorithm */
        /* --------------------------------------------------------------------------------------------
         */
        pl = moris::prm::create_nonlinear_algorithm_parameter_list();
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
        tParameterlist( 2 ).add_parameter_list( pl );

        /* --------------------------------------------------------------------------------------------
         */
        /*                                       Nonlinear Solver */
        /* --------------------------------------------------------------------------------------------
         */
        pl = moris::prm::create_nonlinear_solver_parameter_list();
        pl.set( "NLA_DofTypes", "UX,UY" );
        tParameterlist( 3 ).add_parameter_list( pl );

        /* --------------------------------------------------------------------------------------------
         */
        /*                                        Time Algorithm */
        /* --------------------------------------------------------------------------------------------
         */
        pl = moris::prm::create_time_solver_algorithm_parameter_list();
        pl.set( "TSA_Num_Time_Steps", 1 );
        pl.set( "TSA_Time_Frame", 10.0 );
        tParameterlist( 4 ).add_parameter_list( pl );

        /* --------------------------------------------------------------------------------------------
         */
        /*                                         Time Solver */
        /* --------------------------------------------------------------------------------------------
         */
        pl = moris::prm::create_time_solver_parameter_list();
        pl.set( "TSA_DofTypes", "UX,UY" );
        pl.set( "TSA_Output_Indices", "0" );
        pl.set( "TSA_Output_Criteria", F2STR( Output_Criterion ) );
        tParameterlist( 5 ).add_parameter_list( pl );

        /* --------------------------------------------------------------------------------------------
         */
        /*                                        Solver Warehouse */
        /* --------------------------------------------------------------------------------------------
         */
        pl = moris::prm::create_solver_warehouse_parameterlist();
        // pl.set("SOL_save_operator_to_matlab", "jacobian");
        tParameterlist( 6 ).add_parameter_list( pl );

        /* --------------------------------------------------------------------------------------------
         */
        /*                                      Preconditioner List */
        /* --------------------------------------------------------------------------------------------
         */
        pl = moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::NONE );
        tParameterlist( 7 ).add_parameter_list( pl );
    }

    void MSIParameterList( moris::Vector< moris::Submodule_Parameter_Lists > &tParameterlist )
    {
        tParameterlist.resize( 1 );

        /* --------------------------------------------------------------------------------------------
         */
        /*                                       MSI Parameter List */
        /* --------------------------------------------------------------------------------------------
         */
        Parameter_List pl = prm::create_msi_parameter_list();
        pl.set( "UX", 0 );
        pl.set( "UY", 0 );
        tParameterlist( 0 ).add_parameter_list( pl );
    }

    /* ----------------------------------------------------------------------------------------------
     */
    /*                                           ### VIS ### */
    /* ----------------------------------------------------------------------------------------------
     */
    void VISParameterList( moris::Vector< moris::Submodule_Parameter_Lists > &tParameterlist )
    {
        tParameterlist.resize( 1 );

        /* --------------------------------------------------------------------------------------------
         */
        /*                                          VIS Parameter */
        /* --------------------------------------------------------------------------------------------
         */
        Parameter_List pl = prm::create_vis_parameter_list();
        pl.set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName ) );
        pl.set( "Mesh_Type", (uint)vis::VIS_Mesh_Type::STANDARD );

        // std::string tSetNames = tDomain;
        std::string tSetNames   = tDomain + "," + tContactInterface;
        std::string tFieldNames = "GAP";
        std::string tFieldTypes = "FACETED_AVG";
        std::string tIQINames   = "IQIGap";

        // add the displacement components
        tFieldNames += ",UX,UY";
        tFieldTypes += ",NODAL,NODAL";
        tIQINames += ",IQIDispX,IQIDispY";

        // add stress components as STRESSX, STRESSY
        moris::Vector< std::string >
                                     tSides      = { "Top", "Bottom" };
        moris::Vector< std::string > tComponents = { "X", "Y" };

        for ( auto const &tSide : tSides )
        {
            for ( auto const &tComponent : tComponents )
            {

                // add the normal stress components
                // tFieldNames += ",STRESS" + tSide + tComponent;
                // tFieldTypes += ",ELEMENTAL_AVG"; // ELEMENTAL_AVG
                // tIQINames += ",IQIStress"  + tSide + tComponent;

                tFieldNames += ",STRESSBulk" + tSide + tComponent;
                tFieldTypes += ",ELEMENTAL_AVG";    // ELEMENTAL_AVG
                tIQINames += ",IQIStress" + tSide + tComponent;
            }

            // add the traction components
            tFieldNames += ",ContactPressure" + tSide;
            tFieldTypes += ",FACETED_INT";
            tIQINames += ",IQIContactPressure" + tSide;
        }

        for ( auto const &tComponent : tComponents )
        {
            // add the normal vector components
            tFieldNames += ",NORMAL" + tComponent;
            tFieldTypes += ",FACETED_AVG";
            tIQINames += ",IQINormal" + tComponent;
        }

        pl.set( "Set_Names", tSetNames );
        pl.set( "Field_Names", tFieldNames );
        pl.set( "Field_Type", tFieldTypes );
        pl.set( "IQI_Names", tIQINames );
        pl.set( "Save_Frequency", 1 );
        pl.set( "Time_Offset", 0.1 );
        tParameterlist( 0 ).add_parameter_list( pl );
    }

    void MORISGENERALParameterList( moris::Vector< moris::Submodule_Parameter_Lists > &tParameterlist ) {}
}    // namespace moris
#ifdef __cplusplus
}
#endif
