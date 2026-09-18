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
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_Matrix.hpp"
#include "cl_TSA_Time_Solver.hpp"
#include "parameters.hpp"
#include "linalg_typedefs.hpp"
#include "moris_typedefs.hpp"
#include <string>

//---------------------------------------------------------------

#define F2STR( func ) #func

using moris::fem::IWG_Type;

//---------------------------------------------------------------


//---------------------------------------------------------------

#ifdef __cplusplus
extern "C" {
#endif
namespace moris
{
    /* -------------------------------------- Geometry Setting -------------------------------------- */
    //real tDomainHeight = 1.0; // depends on the bottom body height, gap and radius of the top body
    real tDomainWidth  = 260.0;

    uint tNumElems  = 16;
    uint tNumElemsX = tNumElems;
    uint tNumElemsY = tNumElems;

    uint tInterfaceRefinement = 0;

    real tInitialGap  = 5.0;
    real tPenetration = 70.0; // Mlika: 70
    real tVerticalDisplacement = ( tInitialGap + tPenetration );    // only the upper body moves

    /* --------------------------------------- Top Geometry ----------------------------------------- */
    real tTopYShift         = tInitialGap;
    real tTopXShift         = 0.0;

    real tTopThickness      = 5.0; // thickness of the top body
    const real tTopRadius   = 100; // outer radius of upper body

    /* --------------------------------------- Bottom Geometry -------------------------------------- */
    real tBottomHeight = 50.0;

    real tDomainHeight = tBottomHeight + tInitialGap + tTopRadius;

    /* --------------------------------------- Material & BC Parameters ----------------------------- */
    real tTopEmodInner = 1.0e+5;
    real tTopEmodOuter = 1.0e+3;
    real tBottomEmod   = 300.0;
    real tTopPois      = 0.3;
    real tBottomPois   = 0.3;

    /* --------------------------------------- Domain Setting --------------------------------------- */
    real tDomainOffsetX = -0.5 * tDomainWidth;
    real tDomainOffsetY = -tBottomHeight;

    std::string tInterpolationOrder     = "2";

    std::string tDomainTopInner         = "HMR_dummy_n_p3,HMR_dummy_c_p3";
    std::string tDomainTopOuter         = "HMR_dummy_n_p2,HMR_dummy_c_p2";
    std::string tDomainBottom           = "HMR_dummy_n_p1,HMR_dummy_c_p1";
    std::string tDomain                 = tDomainTopInner + "," + tDomainTopOuter + "," + tDomainBottom;

    std::string tContactInterfaceTop    = "ncss|iside_b0_2_b1_0|iside_b0_1_b1_0";
    std::string tContactInterfaceBottom = "ncss|iside_b0_1_b1_0|iside_b0_2_b1_0";
    std::string tContactInterface       = tContactInterfaceBottom + "," + tContactInterfaceTop;

    /* ------------------------------------------ File I/O ------------------------------------------ */
    std::string tOutputFileName = "Elastic_Half_Ring";
    std::string tSoFile         = tOutputFileName + ".so";
    std::string tHdf5File       = tOutputFileName + ".hdf5";

    /* ------------------------------------------- Solver ------------------------------------------- */
    int tMaxIterations   = 200;
    real tRelResNormDrop = 1e-8;
    real tRelaxation     = 0.6;

    int tNewtonMaxIter = 20;
    real tNewtonRelRes = 1.0e-12;

    int  tLoadControlSteps    = 20;
    real tLoadControlFactor   = 0.01; //1.0/tLoadControlSteps;
    real tLoadControlRelRes   = 1e0;
    real tLoadControlExponent = 1.0;

    real tRemapResidualChange     = -1e-5;
    int  tRemapLoadStepFrequency  = 1;
    int  tRemapIterationFrequency = 1;
    auto tRaytracingStrategy      = sol::SolverRaytracingStrategy::EveryNthIteration;
    real tMaxNegativeRayLength    = -8.0;
    real tMaxPositiveRayLength    = 50.0;

    mtk::Integration_Order tNonconformalIntegrationOrder = mtk::Integration_Order::BAR_4;

    /* ------------------------------------------- Contact ------------------------------------------ */
    std::string tContactType          = "mlika";
    std::string tContactBias          = "unsymmetric";
    std::string tContactStabilization = "1.00/0.0";

    /* -------------------------------------- Control Variables ------------------------------------- */
    bool tOnlyGenerateMesh   = false;
    bool tUseGhost           = true;

    /* ---------------------------------------------------------------------------------------------- */
    /*                                     User-Defined Functions                                     */
    /* ---------------------------------------------------------------------------------------------- */

    void Func_Dirichlet( Matrix< DDRMat >   &aPropMatrix,
            Vector< Matrix< DDRMat > >      &aParameters,
            fem::Field_Interpolator_Manager *aFIManager )
    {
        real tLoadFactor = gLogger.get_action_data("NonLinearAlgorithm", "NLBGS", "Solve", "LoadFactor");
        aPropMatrix      = aParameters( 0 ) * tLoadFactor;
    }

    void Func_Select_X( Matrix< DDRMat >    &aPropMatrix,
            Vector< Matrix< DDRMat > >      &aParameters,
            fem::Field_Interpolator_Manager *aFIManager )
    {
        aPropMatrix.set_size( 2, 2, 0.0 );
        aPropMatrix( 0, 0 ) = 1.0;
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
        aParameterLists.set( "domain_dimensions", tDomainWidth, tDomainHeight );
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
        aParameterLists.set( "number_of_phases", 5 );

        Matrix< DDUMat > tPhaseMap( 16, 1, 4 ); // num rows, num cols, initial value
                                                // 4 is default = top void
        tPhaseMap(  9 )  = 3;	           	// upper body inner shell
        tPhaseMap( 13 )  = 2;	           	// upper body outer shell
        tPhaseMap( 14 )  = 1;	           	// lower body
        tPhaseMap( 15 )  = 0;	           	// interface

        aParameterLists.set( "phase_table", moris::ios::stringify( tPhaseMap ) );
        aParameterLists.set( "print_phase_table", true );
        aParameterLists.set( "output_mesh_file","gen.exo");

        /* -------------------------------------------------------------------------------------------- */
        /*                                      Geometry Parameter                                      */
        /* -------------------------------------------------------------------------------------------- */

        /* ----------------------------------------- Upper Body --------------------------------------- */

        /* ----------------------------------------- Top Plane ---------------------------------------- */

        // inner circle for upper body
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::CIRCLE ) );
        aParameterLists.set( "center_x", tTopXShift );
        aParameterLists.set( "center_y", tTopYShift + tTopRadius );
        aParameterLists.set( "radius", tTopRadius-2.0*tTopThickness );
        aParameterLists.set( "number_of_refinements", tInterfaceRefinement );
        aParameterLists.set( "refinement_mesh_index", 0 );

        // mid circle for upper body
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::CIRCLE ) );
        aParameterLists.set( "center_x", tTopXShift );
        aParameterLists.set( "center_y", tTopYShift + tTopRadius );
        aParameterLists.set( "radius", tTopRadius-tTopThickness );
        aParameterLists.set( "number_of_refinements", tInterfaceRefinement );
        aParameterLists.set( "refinement_mesh_index", 0 );

        // outer circle for upper body
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::CIRCLE ) );
        aParameterLists.set( "center_x", tTopXShift );
        aParameterLists.set( "center_y", tTopYShift + tTopRadius );
        aParameterLists.set( "radius", tTopRadius );
        aParameterLists.set( "number_of_refinements", tInterfaceRefinement );
        aParameterLists.set( "refinement_mesh_index", 0 );

        /* ----------------------------------------- Lower Body --------------------------------------- */

        /* --------------------------------------- Upper Plane ---------------------------------------- */
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
        aParameterLists.set( "center_x", 0.0 );
        aParameterLists.set( "center_y", 0.0 );
        aParameterLists.set( "normal_x", 0.0 );
        aParameterLists.set( "normal_y", 1.0 );
        aParameterLists.set( "number_of_refinements", tInterfaceRefinement );
        aParameterLists.set( "refinement_mesh_index", 0 );
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
        pl.set( "phase_name", "PhaseTopOuter" );
        pl.set( "phase_indices", "2" );
        aParameterLists( FEM::PHASES ).add_parameter_list( pl );

        pl = prm::create_phase_parameter_list();
        pl.set( "phase_name", "PhaseTopInner" );
        pl.set( "phase_indices", "3" );
        aParameterLists( FEM::PHASES ).add_parameter_list( pl );

        pl = prm::create_phase_parameter_list();
        pl.set( "phase_name", "PhaseBottom" );
        pl.set( "phase_indices", "1" );
        aParameterLists( FEM::PHASES ).add_parameter_list( pl );

        pl = prm::create_phase_parameter_list();
        pl.set( "phase_name", "PhaseTop" );
        pl.set( "phase_indices", "2,3" );
        aParameterLists( FEM::PHASES ).add_parameter_list( pl );

        pl = prm::create_phase_parameter_list();
        pl.set( "phase_name", "PhaseVoidUpper" );
        pl.set( "phase_indices", "4" );
        aParameterLists( FEM::PHASES ).add_parameter_list( pl );

        /* -------------------------------------------------------------------------------------------- */
        /*                                           Properties                                         */
        /* -------------------------------------------------------------------------------------------- */

        /* -------------------------------------- Youngs Modulus -------------------------------------- */
        pl = prm::create_property_parameter_list();
        pl.set( "property_name", "PropYoungsTopInner" );
        pl.set( "function_parameters", std::to_string( tTopEmodInner ) );
        aParameterLists( FEM::PROPERTIES ).add_parameter_list( pl );

        pl = prm::create_property_parameter_list();
        pl.set( "property_name", "PropYoungsTopOuter" );
        pl.set( "function_parameters", std::to_string( tTopEmodOuter ) );
        aParameterLists( FEM::PROPERTIES ).add_parameter_list( pl );

        pl = prm::create_property_parameter_list();
        pl.set( "property_name", "PropYoungsBottom" );
        pl.set( "function_parameters", std::to_string( tBottomEmod ) );
        aParameterLists( FEM::PROPERTIES ).add_parameter_list( pl );

        /* ------------------------------------------ Poisson ----------------------------------------- */
        pl = prm::create_property_parameter_list();
        pl.set( "property_name", "PropPoissonTop" );
        pl.set( "function_parameters", std::to_string( tTopPois ) );
        aParameterLists( FEM::PROPERTIES ).add_parameter_list( pl );

        pl = prm::create_property_parameter_list();
        pl.set( "property_name", "PropPoissonBottom" );
        pl.set( "function_parameters", std::to_string( tBottomPois ) );
        aParameterLists( FEM::PROPERTIES ).add_parameter_list( pl );

        /* --------------------------------------- Dirichlet Top -------------------------------------- */
        pl = prm::create_property_parameter_list();
        pl.set( "property_name", "PropDirichletTop" );
        pl.set( "function_parameters", "0.0;" + std::to_string( -tVerticalDisplacement ) );
        pl.set( "value_function", F2STR( Func_Dirichlet ) );
        aParameterLists( FEM::PROPERTIES ).add_parameter_list( pl );

        /* ------------------------------------- Dirichlet Bottom ------------------------------------- */
        pl = prm::create_property_parameter_list();
        pl.set( "property_name", "PropFixed" );
        pl.set( "function_parameters", "0.0;0.0" );
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
        pl.set( "constitutive_name", "MaterialTopOuter" );
        pl.set( "phase_name", "PhaseTopOuter" );
        pl.set( "constitutive_type",
                (uint)fem::Constitutive_Type::STRUC_NON_LIN_ISO_COMPRESSIBLE_NEO_HOOKEAN_WRIGGERS );
        pl.set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
        pl.set( "properties",
                "PropYoungsTopOuter,YoungsModulus;"
                "PropPoissonTop,PoissonRatio" );
        pl.set( "model_type",  fem::Model_Type::PLANE_STRAIN ) ;
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list( pl );

        pl = prm::create_constitutive_model_parameter_list();
        pl.set( "constitutive_name", "MaterialTopInner" );
        pl.set( "phase_name", "PhaseTopInner" );
        pl.set( "constitutive_type",
                (uint)fem::Constitutive_Type::STRUC_NON_LIN_ISO_COMPRESSIBLE_NEO_HOOKEAN_WRIGGERS );
        pl.set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
        pl.set( "properties",
                "PropYoungsTopInner,YoungsModulus;"
                "PropPoissonTop,PoissonRatio" );
        pl.set( "model_type",  fem::Model_Type::PLANE_STRAIN ) ;
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list( pl );

        /* -------------------------------- Elastic Material Bottom ----------------------------------- */
        pl = prm::create_constitutive_model_parameter_list();
        pl.set( "constitutive_name", "MaterialBottom" );
        pl.set( "phase_name", "PhaseBottom" );
        pl.set( "constitutive_type",
                (uint)fem::Constitutive_Type::STRUC_NON_LIN_ISO_COMPRESSIBLE_NEO_HOOKEAN_WRIGGERS );
        pl.set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
        pl.set( "properties",
                "PropYoungsBottom,YoungsModulus;"
                "PropPoissonBottom,PoissonRatio" );
        pl.set( "model_type",  fem::Model_Type::PLANE_STRAIN ) ;
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list( pl );

        /* -------------------------------------------------------------------------------------------- */
        /*                                     Stabilization Parameter                                  */
        /* -------------------------------------------------------------------------------------------- */

        pl = prm::create_stabilization_parameter_parameter_list();
        pl.set( "stabilization_name", "SPDirichletNitscheTop" );
        pl.set( "stabilization_type", (uint)fem::Stabilization_Type::DIRICHLET_NITSCHE );
        pl.set( "function_parameters", "100.0" );
        pl.set( "leader_phase_name", "PhaseTopOuter" );
        pl.set( "leader_properties", "PropYoungsTopOuter,Material" );
        aParameterLists( FEM::STABILIZATION ).add_parameter_list( pl );

        pl = prm::create_stabilization_parameter_parameter_list();
        pl.set( "stabilization_name", "SPDirichletNitscheTopInner" );
        pl.set( "stabilization_type", (uint)fem::Stabilization_Type::DIRICHLET_NITSCHE );
        pl.set( "function_parameters", "100.0" );
        pl.set( "leader_phase_name", "PhaseTopInner" );
        pl.set( "leader_properties", "PropYoungsTopInner,Material" );
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
        pl.set( "function_parameters", tContactStabilization );
        pl.set( "leader_phase_name", "PhaseTopOuter" );
        pl.set( "leader_properties", "PropYoungsTopOuter,Material" );
        aParameterLists( FEM::STABILIZATION ).add_parameter_list( pl );

        pl = prm::create_stabilization_parameter_parameter_list();
        pl.set( "stabilization_name", "SPContactInterfaceBottom" );
        pl.set( "stabilization_type", (uint)fem::Stabilization_Type::DIRICHLET_NITSCHE );
        pl.set( "function_parameters", tContactStabilization );
        pl.set( "leader_phase_name", "PhaseBottom" );
        pl.set( "leader_properties", "PropYoungsBottom,Material" );
        aParameterLists( FEM::STABILIZATION ).add_parameter_list( pl );

        // create ghost penalty displacement Matrix
        pl = prm::create_stabilization_parameter_parameter_list();
        pl.set( "stabilization_name", "SPGPDisplTop" );
        pl.set( "stabilization_type", (uint)fem::Stabilization_Type::GHOST_DISPL );
        pl.set( "function_parameters", "0.05" );
        pl.set( "leader_phase_name", "PhaseTopOuter" );
        pl.set( "follower_phase_name", "PhaseTopOuter" );
        pl.set( "leader_properties", "PropYoungsTopOuter,Material" );
        aParameterLists( FEM::STABILIZATION ).add_parameter_list( pl );

        pl = prm::create_stabilization_parameter_parameter_list();
        pl.set( "stabilization_name", "SPGPDisplTopInner" );
        pl.set( "stabilization_type", (uint)fem::Stabilization_Type::GHOST_DISPL );
        pl.set( "function_parameters", "0.05" );
        pl.set( "leader_phase_name", "PhaseTopInner" );
        pl.set( "follower_phase_name", "PhaseTopInner" );
        pl.set( "leader_properties", "PropYoungsTopInner,Material" );
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


        // displacement interface, tied contact between outer und inner shell
        pl = prm::create_stabilization_parameter_parameter_list();
        pl.set( "stabilization_name", "SP_Nitsche_Displ_Interface" );
        pl.set( "leader_phase_name", "PhaseTopOuter" );
        pl.set( "follower_phase_name", "PhaseTopInner" );
        pl.set( "stabilization_type", (uint)fem::Stabilization_Type::NITSCHE_INTERFACE );
        pl.set( "function_parameters", "50.0" );
        pl.set( "leader_properties", "PropYoungsTopOuter,Material" );
        pl.set( "follower_properties", "PropYoungsTopInner,Material" );
        aParameterLists( FEM::STABILIZATION ).add_parameter_list( pl );

        /* -------------------------------------------------------------------------------------------- */
        /*                                               IWG                                            */
        /* -------------------------------------------------------------------------------------------- */

        /* -------------------------------------- Body Material Top ----------------------------------- */
        pl = prm::create_IWG_parameter_list();
        pl.set( "IWG_name", "IWGBulkStructTopOuter" );
        pl.set( "IWG_type", (uint)fem::IWG_Type::STRUC_NON_LINEAR_BULK_PF );
        pl.set( "IWG_bulk_type", (uint)fem::Element_Type::BULK );
        pl.set( "dof_residual", "UX,UY" );
        pl.set( "leader_phase_name", "PhaseTopOuter" );
        pl.set( "leader_constitutive_models", "MaterialTopOuter,ElastLinIso" );
        aParameterLists( FEM::IWG ).add_parameter_list( pl );

        pl = prm::create_IWG_parameter_list();
        pl.set( "IWG_name", "IWGBulkStructTopInner" );
        pl.set( "IWG_type", (uint)fem::IWG_Type::STRUC_NON_LINEAR_BULK_PF );
        pl.set( "IWG_bulk_type", (uint)fem::Element_Type::BULK );
        pl.set( "dof_residual", "UX,UY" );
        pl.set( "leader_phase_name", "PhaseTopInner" );
        pl.set( "leader_constitutive_models", "MaterialTopInner,ElastLinIso" );
        aParameterLists( FEM::IWG ).add_parameter_list( pl );

        // displacement interface, fiber-filler -> tied contact between the two materials
        pl = prm::create_IWG_parameter_list();
        pl.set( "IWG_name", "IWGBulkInterfaceTop" );
        pl.set( "IWG_type", (uint)fem::IWG_Type::STRUC_NON_LINEAR_INTERFACE_UNSYMMETRIC_NITSCHE_PF );
        pl.set( "IWG_bulk_type", (uint)fem::Element_Type::DOUBLE_SIDESET );
        pl.set( "leader_phase_name", "PhaseTopOuter" );
        pl.set( "follower_phase_name", "PhaseTopInner" );
        pl.set( "dof_residual", "UX,UY" );
        pl.set( "leader_dof_dependencies", "UX,UY" );
        pl.set( "follower_dof_dependencies", "UX,UY" );
        pl.set( "leader_constitutive_models", "MaterialTopOuter,ElastLinIso" );
        pl.set( "follower_constitutive_models", "MaterialTopInner,ElastLinIso" );
        pl.set( "stabilization_parameters", "SP_Nitsche_Displ_Interface,NitscheInterface" );
        aParameterLists( FEM::IWG ).add_parameter_list( pl );

        /* ------------------------------------ Body Material Bottom  --------------------------------- */
        pl = prm::create_IWG_parameter_list();
        pl.set( "IWG_name", "IWGBulkStructBottom" );
        pl.set( "IWG_type", (uint)fem::IWG_Type::STRUC_NON_LINEAR_BULK_PF );
        pl.set( "IWG_bulk_type", (uint)fem::Element_Type::BULK );
        pl.set( "dof_residual", "UX,UY" );
        pl.set( "leader_phase_name", "PhaseBottom" );
        pl.set( "leader_constitutive_models", "MaterialBottom,ElastLinIso" );
        aParameterLists( FEM::IWG ).add_parameter_list( pl );

        /* ---------------------------------------- Dirichlet Top ------------------------------------- */
        pl = prm::create_IWG_parameter_list();
        pl.set( "IWG_name", "IWGDirichletTop" );
        pl.set( "IWG_type", (uint)fem::IWG_Type::STRUC_NON_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE_PF );
        pl.set( "IWG_bulk_type", (uint)fem::Element_Type::SIDESET );
        pl.set( "dof_residual", "UX,UY" );
        pl.set( "leader_phase_name", "PhaseTopOuter" );
        pl.set( "side_ordinals", "3" );
        //pl.set( "neighbor_phases", "PhaseVoidUpper" );
        pl.set( "leader_properties", "PropDirichletTop,Dirichlet" );
        pl.set( "leader_constitutive_models", "MaterialTopOuter,ElastLinIso" );
        pl.set( "stabilization_parameters", "SPDirichletNitscheTop,DirichletNitsche" );
        aParameterLists( FEM::IWG ).add_parameter_list( pl );

        pl = prm::create_IWG_parameter_list();
        pl.set( "IWG_name", "IWGDirichletTop" );
        pl.set( "IWG_type", (uint)fem::IWG_Type::STRUC_NON_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE_PF );
        pl.set( "IWG_bulk_type", (uint)fem::Element_Type::SIDESET );
        pl.set( "dof_residual", "UX,UY" );
        pl.set( "leader_phase_name", "PhaseTopInner" );
        pl.set( "side_ordinals", "3" );
        //pl.set( "neighbor_phases", "PhaseVoidUpper" );
        pl.set( "leader_properties", "PropDirichletTop,Dirichlet" );
        pl.set( "leader_constitutive_models", "MaterialTopInner,ElastLinIso" );
        pl.set( "stabilization_parameters", "SPDirichletNitscheTopInner,DirichletNitsche" );
        aParameterLists( FEM::IWG ).add_parameter_list( pl );

        /* -------------------------------------- Dirichlet Bottom ------------------------------------ */
        pl = prm::create_IWG_parameter_list();
        pl.set( "IWG_name", "IWGDirichletBottom" );
        pl.set( "IWG_type", (uint)fem::IWG_Type::STRUC_NON_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE_PF );
        pl.set( "IWG_bulk_type", (uint)fem::Element_Type::SIDESET );
        pl.set( "dof_residual", "UX,UY" );
        pl.set( "leader_phase_name", "PhaseBottom" );
        pl.set( "side_ordinals", "1" );
        pl.set( "leader_properties", "PropFixed,Dirichlet" );
        pl.set( "leader_constitutive_models", "MaterialBottom,ElastLinIso" );
        pl.set( "stabilization_parameters", "SPDirichletNitscheBottom,DirichletNitsche" );
        aParameterLists( FEM::IWG ).add_parameter_list( pl );


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
        pl.set( "leader_phase_name", "PhaseTopOuter" );
        pl.set( "neighbor_phases", "PhaseVoid" );
        pl.set( "follower_phase_name", "PhaseBottom" );
        pl.set( "leader_constitutive_models", "MaterialTopOuter,ElastLinIso" );
        pl.set( "follower_constitutive_models", "MaterialBottom,ElastLinIso" );
        pl.set( "stabilization_parameters", "SPContactInterfaceTop,NitscheInterface" );
        aParameterLists( FEM::IWG ).add_parameter_list( pl );

        pl = prm::create_IWG_parameter_list();
        pl.set( "IWG_name", "IWGContactInterfaceBottom" );
        pl.set( "IWG_type", (uint)tContactIWGType );
        pl.set( "IWG_bulk_type", (uint)fem::Element_Type::NONCONFORMAL_SIDESET );
        pl.set( "dof_residual", "UX,UY" );
        pl.set( "leader_phase_name", "PhaseBottom" );
        pl.set( "neighbor_phases", "PhaseVoid" );
        pl.set( "follower_phase_name", "PhaseTopOuter" );
        pl.set( "leader_constitutive_models", "MaterialBottom,ElastLinIso" );
        pl.set( "follower_constitutive_models", "MaterialTopOuter,ElastLinIso" );
        pl.set( "stabilization_parameters", "SPContactInterfaceBottom,NitscheInterface" );
        aParameterLists( FEM::IWG ).add_parameter_list( pl );

        if ( tUseGhost )
        {
            // ghost  displacement
            pl = prm::create_IWG_parameter_list();
            pl.set( "IWG_name", "IWGGPDisplTop" );
            pl.set( "IWG_type", (uint)fem::IWG_Type::GHOST_NORMAL_FIELD );
            pl.set( "IWG_bulk_type", (uint)fem::Element_Type::DOUBLE_SIDESET );
            pl.set( "dof_residual", "UX,UY" );
            pl.set( "leader_phase_name", "PhaseTopOuter" );
            pl.set( "follower_phase_name", "PhaseTopOuter" );
            pl.set( "stabilization_parameters", "SPGPDisplTop,GhostSP" );
            aParameterLists( FEM::IWG ).add_parameter_list( pl );

            // ghost  displacement
            pl = prm::create_IWG_parameter_list();
            pl.set( "IWG_name", "IWGGPDisplTopInner" );
            pl.set( "IWG_type", (uint)fem::IWG_Type::GHOST_NORMAL_FIELD );
            pl.set( "IWG_bulk_type", (uint)fem::Element_Type::DOUBLE_SIDESET );
            pl.set( "dof_residual", "UX,UY" );
            pl.set( "leader_phase_name", "PhaseTopInner" );
            pl.set( "follower_phase_name", "PhaseTopInner" );
            pl.set( "stabilization_parameters", "SPGPDisplTopInner,GhostSP" );
            aParameterLists( FEM::IWG ).add_parameter_list( pl );

            // ghost  displacement
            pl = prm::create_IWG_parameter_list();
            pl.set( "IWG_name", "IWGGPDisplBottom" );
            pl.set( "IWG_type", (uint)fem::IWG_Type::GHOST_NORMAL_FIELD );
            pl.set( "IWG_bulk_type", (uint)fem::Element_Type::DOUBLE_SIDESET );
            pl.set( "dof_residual", "UX,UY" );
            pl.set( "leader_phase_name", "PhaseBottom" );
            pl.set( "follower_phase_name", "PhaseBottom" );
            pl.set( "stabilization_parameters", "SPGPDisplBottom,GhostSP" );
            aParameterLists( FEM::IWG ).add_parameter_list( pl );
        }

        /* -------------------------------------------------------------------------------------------- */
        /*                                               IQI                                            */
        /* -------------------------------------------------------------------------------------------- */

        //Vector< std::pair< std::string, std::string > > tSides      = { { "TopOuter", "Bottom" },
        //                                                                { "Bottom", "TopOuter" },
        //                                                                { "TopInner", "Bottom" },
        //                                                                { "Bottom", "TopInner" } };
        Vector< std::pair< std::string, std::string > > tSides      = { { "Top", "Bottom" },
                                                                        { "Bottom", "Top" } };
        Vector< std::pair< int, std::string > >         tComponents = { { 0, "X" }, { 1, "Y" } };
        for ( auto const &[ tLeaderSide, tFollowerSide ] : tSides )
        {
            for ( auto const &[ tDir, tAxis ] : tComponents )
            {
            //    pl = prm::create_IQI_parameter_list();
            //    pl.set( "IQI_name", "IQIStress" + tLeaderSide + tAxis );
            //    pl.set( "IQI_type", (uint)fem::IQI_Type::NORMAL_STRESS_CAUCHY );
            //    pl.set( "IQI_bulk_type", (uint)fem::Element_Type::BULK );
            //    pl.set( "leader_constitutive_models", "Material" + tLeaderSide + ",ElastLinIso" );
            //    pl.set( "leader_phase_name", "Phase" + tLeaderSide );
            //    pl.set( "vectorial_field_index", tDir );
            //    aParameterLists( FEM::IQI ).add_parameter_list( pl );

                pl = prm::create_IQI_parameter_list();
                pl.set( "IQI_name", "IQIDisp" + tLeaderSide + tAxis );
                pl.set( "IQI_type", (uint)fem::IQI_Type::DOF );
                pl.set( "IQI_bulk_type", (uint)fem::Element_Type::BULK );
                pl.set( "dof_quantity", "UX,UY" );
                pl.set( "leader_phase_name", "Phase" + tLeaderSide );
                pl.set( "vectorial_field_index", tDir );
                aParameterLists( FEM::IQI ).add_parameter_list( pl );
            }

            //pl = prm::create_IQI_parameter_list();
            //pl.set( "IQI_name", "IQIVonMises" + tLeaderSide );
            //pl.set( "IQI_type", (uint)fem::IQI_Type::VON_MISES_STRESS_CAUCHY );
            //pl.set( "IQI_bulk_type", (uint)fem::Element_Type::BULK );
            //pl.set( "leader_constitutive_models", "Material" + tLeaderSide + ",ElastLinIso" );
            //pl.set( "leader_phase_name", "Phase" + tLeaderSide );
            //aParameterLists( FEM::IQI ).add_parameter_list( pl );

            /*pl = prm::create_IQI_parameter_list();
            pl.set( "IQI_name", "IQIGap" + tLeaderSide );
            pl.set( "IQI_type", (uint)fem::IQI_Type::GAP );
            pl.set( "IQI_bulk_type", (uint)fem::Element_Type::NONCONFORMAL_SIDESET );
            pl.set( "leader_phase_name", "Phase" + tLeaderSide );
            pl.set( "neighbor_phases", "PhaseVoid" );
            pl.set( "follower_phase_name", "Phase" + tFollowerSide );
            aParameterLists( FEM::IQI ).add_parameter_list( pl );
            */
        }

        // /* ---------------------------------------- Ray Length ----------------------------------------
        //pl = prm::create_IQI_parameter_list();
        //pl.set( "IQI_name", "IQIGap" );
        //pl.set( "IQI_type", (uint)fem::IQI_Type::GAP );
        //pl.set( "leader_phase_name", "PhaseBottom" );
        //pl.set( "follower_phase_name", "PhaseTopOuter" );
        //aParameterLists( FEM::IQI ).add_parameter_list( pl );

        /* -------------------------------------------------------------------------------------------- */
        /*                                   Computation Parameter List                                 */
        /* -------------------------------------------------------------------------------------------- */
        aParameterLists( FEM::COMPUTATION ).set( "is_analytical_forward", true );

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
        Parameter_List pl = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL );
        // pl.set( "Solver_Type", "Amesos_Umfpack" );
        aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( pl );

        /* -------------------------------------------------------------------------------------------- */
        /*                                         Linear Solver                                        */
        /* -------------------------------------------------------------------------------------------- */
        pl = moris::prm::create_linear_solver_parameter_list();
        aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list( pl );

        /* -------------------------------------------------------------------------------------------- */
        /*                                      Nonlinear Algorithm                                     */
        /* -------------------------------------------------------------------------------------------- */
        pl = moris::prm::create_nonlinear_algorithm_parameter_list();
        pl.set("NLA_Linear_solver",                   0);
        pl.set("NLA_max_iter",                        tNewtonMaxIter );
        pl.set("NLA_rel_res_norm_drop",               tNewtonRelRes);
        pl.set("NLA_time_offset",                     0.0);
        //pl.set("NLA_relaxation_parameter",            tRelaxation);
        pl.set("NLA_relaxation_strategy", 	    sol::SolverRelaxationType::InvResNormAdaptive );
        pl.set("NLA_relaxation_parameter",	    1.0 );
        pl.set("NLA_relaxation_damping",  	    0.5 );
        pl.set("NLA_remap_strategy",                  (uint)tRaytracingStrategy);
        pl.set("NLA_remap_load_step_frequency",       tRemapLoadStepFrequency);
        pl.set("NLA_remap_iteration_frequency",       tRemapIterationFrequency);
        pl.set("NLA_remap_residual_change_tolerance", tRemapResidualChange);
        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list( pl );

        pl = prm::create_nonlinear_algorithm_parameter_list();
        pl.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NLBGS_SOLVER );
        pl.set( "NLA_Linear_solver",         0 );
        pl.set( "NLA_rel_res_norm_drop",     tRelResNormDrop );
        pl.set( "NLA_relaxation_parameter",  1.0 );
        pl.set( "NLA_max_iter",              tMaxIterations );
        pl.set( "NLA_ref_iter",              1 );
        pl.set( "NLA_time_offset",           1.0);
        pl.set( "NLA_load_control_strategy", sol::SolverLoadControlType::Linear );
        pl.set( "NLA_load_control_factor",   tLoadControlFactor );
        pl.set( "NLA_load_control_steps",    tLoadControlSteps );
        pl.set( "NLA_load_control_relres",   tLoadControlRelRes );
        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list( pl );

        /* -------------------------------------------------------------------------------------------- */
        /*                                       Nonlinear Solver                                       */
        /* -------------------------------------------------------------------------------------------- */
        pl = moris::prm::create_nonlinear_solver_parameter_list();
        pl.set( "NLA_DofTypes", "UX,UY" );
        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list( pl );

        pl = prm::create_nonlinear_solver_parameter_list();
        pl.set( "NLA_Nonlinear_solver_algorithms", "1" );
        pl.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NLBGS_SOLVER );
        pl.set( "NLA_Sub_Nonlinear_Solver", "0" );
        pl.set( "NLA_DofTypes", "UX,UY" );
        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list( pl );

        /* -------------------------------------------------------------------------------------------- */
        /*                                        Time Algorithm                                        */
        /* -------------------------------------------------------------------------------------------- */
        pl = moris::prm::create_time_solver_algorithm_parameter_list();
        pl.set( "TSA_Num_Time_Steps", 1 );
        pl.set( "TSA_Time_Frame", 10.0 ); // nonlinear 100.0
        pl.set( "TSA_Nonlinear_Solver", 1 );  
        aParameterLists( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list( pl );

        /* -------------------------------------------------------------------------------------------- */
        /*                                         Time Solver                                          */
        /* -------------------------------------------------------------------------------------------- */
        pl = moris::prm::create_time_solver_parameter_list();
        pl.set( "TSA_DofTypes", "UX,UY" );
        pl.set( "TSA_Output_Indices", "0" );
        pl.set( "TSA_Output_Criteria", F2STR( Output_Criterion ) );
        aParameterLists( SOL::TIME_SOLVERS ).add_parameter_list( pl );

        /* -------------------------------------------------------------------------------------------- */
        /*                                        Solver Warehouse                                      */
        /* -------------------------------------------------------------------------------------------- */
        // pl.set("SOL_save_operator_to_matlab", "jacobian");

        /* -------------------------------------------------------------------------------------------- */
        /*                                      Preconditioner List                                     */
        /* -------------------------------------------------------------------------------------------- */
        pl = moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::NONE );
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
                tIQINames   += ",IQIStress" + tSide + tComponent;

                // add the displacement components
                tFieldNames += ",U" + tSide + tComponent;
                tFieldTypes += ",NODAL";    // ELEMENTAL_AVG
                tIQINames   += ",IQIDisp" + tSide + tComponent;
            }

            tFieldNames += ",STRESSVonMises" + tSide;
            tFieldTypes += ",NODAL";    // ELEMENTAL_AVG
            tIQINames   += ",IQIVonMises" + tSide;

            // add the traction components
            /*tFieldNames += ",ContactPressure" + tSide;
            tFieldTypes += ",NODAL";
            tIQINames += ",IQIContactPressure" + tSide;
            */

            // add the gap
            /*tFieldNames += ",Gap" + tSide;
            tFieldTypes += ",NODAL";
            tIQINames += ",IQIGap" + tSide;
	    */
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
        aParameterLists.set( "Time_Offset", 1.0 );
    }

    void MORISGENERALParameterList( Module_Parameter_Lists &aParameterLists ) {}

}    // namespace moris
#ifdef __cplusplus
}
#endif
