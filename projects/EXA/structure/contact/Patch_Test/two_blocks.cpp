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


//---------------------------------------------------------------

#ifdef __cplusplus
extern "C" {
#endif
namespace moris
{
    /* -------------------------------------- Geometry Setting -------------------------------------- */
    real tDomainHeight = 0.065;
    real tDomainWidth  = 0.05;
    
    real tBoxHeight = 0.01;
    real tBoxWidth  = 0.05;

    real tUpperCenter = 0.04;

    uint tNumElemsX = 6;
    uint tNumElemsY = 6;

    /* --------------------------------------- Material & BC Parameters ----------------------------- */

    real tEmod = 2e11;
    real tPois = 0.3;
    
    /* --------------------------------------- Domain Setting --------------------------------------- */
    real tDomainOffsetX = -0.5 * tDomainWidth;
    real tDomainOffsetY =  0.0;

    std::string tInterpolationOrder     = "2";
    std::string tDomainTop              = "HMR_dummy_n_p2,HMR_dummy_c_p2";
    std::string tDomainBottom           = "HMR_dummy_n_p1,HMR_dummy_c_p1";
    std::string tContactInterfaceTop    = "ncss|iside_b0_2_b1_0|iside_b0_1_b1_0";
    std::string tContactInterfaceBottom = "ncss|iside_b0_1_b1_0|iside_b0_2_b1_0";
    std::string tContactInterface       = tContactInterfaceBottom + "," + tContactInterfaceTop;
    std::string tDomain                 = tDomainTop + "," + tDomainBottom;

    /* ------------------------------------------ File I/O ------------------------------------------ */
    std::string tOutputFileName = "two_blocks";
    std::string tSoFile    = "two_blocks.so";
    std::string tHdf5File  = "two_blocks.hdf5";
    
    /* ------------------------------------------- Solver ------------------------------------------- */
    int tMaxIterations = 180;
    real tRelResNormDrop = 1e-6;
    real tRelaxation = 1.0;
    
    int tNewtonMaxIter = 15;
    real tNewtonRelRes=1.0e-12;
    
    int  tLoadControlSteps = 40;
    real tLoadControlFactor = 1.0/tLoadControlSteps;
    real tLoadControlRelRes = 1e0;
    real tLoadControlExponent = 1.0;
    
    real tT1 = 1.0-10.0/tLoadControlSteps;
    real tT2 = 1.0-2.0/tLoadControlSteps;
    real tT3 = 1.0-1.0/tLoadControlSteps;
    
    real tRemapResidualChange = -1e-2;
    int  tRemapLoadStepFrequency = 1;
    int  tRemapIterationFrequency = 1;
    auto tRaytracingStrategy = sol::SolverRaytracingStrategy::EveryNthIteration;
    real tMaxNegativeRayLength = -0.005;
    real tMaxPositiveRayLength = 0.005;
    
    mtk::Integration_Order tNonconformalIntegrationOrder = mtk::Integration_Order::BAR_4;

    /* ------------------------------------------- Contact ------------------------------------------ */
    //std::string tContactType = "mlika_frieder";
    std::string tContactType = "mlika";
    std::string tContactBias = "neutral";
    bool tUseAnalyticalJacbian= true;
    std::string tConsistentProjection = "1/0";
    std::string tContactStabilization = "100.0/0.0";   


    /* -------------------------------------- Control Variables ------------------------------------- */
    bool tOnlyGenerateMesh   = false;
    bool tUseGhost           = true;
    
    /* ---------------------------------------------------------------------------------------------- */
    /*                                         Field Functions                                        */
    /* ---------------------------------------------------------------------------------------------- */

    /* --------------------------------------- Phase Indexing --------------------------------------- */
    uint Phase_Index_Split( const Bitset< 4 > &aGeometrySigns )
    {
        if ( ! aGeometrySigns.test( 2 ) )
        {
           return 1;    // bottom domain
        }
    
        if ( ! aGeometrySigns.test( 3 ) )
        {
           return 0;    // gap  (void)
        }
           
    if ( ! aGeometrySigns.test( 0 ) )
        {
           return 0;    // upper left domain (void)
        }

    if ( ! aGeometrySigns.test( 1 ) )
        {
           return 2;    // top domain
        }

        return 3;    // upper right domain (void)

      }

    /* ----------------------------------- Property Field Function ---------------------------------- */

    void Func_Const( Matrix< DDRMat >       &aPropMatrix,
            Vector< Matrix< DDRMat > >      &aParameters,
            fem::Field_Interpolator_Manager *aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    void Func_Dirichlet_IFC( 
            Matrix< DDRMat >                &aPropMatrix,
            Vector< Matrix< DDRMat > >      &aParameters,
            fem::Field_Interpolator_Manager *aFIManager )
    {
         real tTime = gLogger.get_action_data("NonLinearAlgorithm", "NLBGS", "Solve", "LoadFactor");

        const Matrix< DDRMat > tX = trans(aFIManager->get_IP_geometry_interpolator()->valx());
        const Matrix< DDRMat > tC   = { {0.0},{tUpperCenter }};
        const Matrix< DDRMat > tEy  = { {0.0},{ 1.0 }};
        
        const real tAngle  = std::acos(0) * std::min(1.0, tTime/tT1);
	const real tYShift = std::min(1.0,std::max(0.0,( tTime - tT1 )/(tT2-tT1))) *(tC(1)-1.5000*tBoxHeight+0.0001);
	
        const Matrix< DDRMat > tRot = { {std::cos(tAngle), -std::sin(tAngle)},{std::sin(tAngle), std::cos(tAngle) }};
        
        const Matrix< DDRMat > tNewX = tC + tRot*(tX-tC) - tYShift*tEy;
    
        aPropMatrix = tNewX - tX;
	
	const Matrix< DDRMat > tRefPoint = {{tBoxHeight/2.0},{tDomainHeight}};
	
	if ( norm(tX-tRefPoint)<1e-3)
	{
	    Matrix<DDRMat> tDisp = aFIManager->get_field_interpolators_for_type( MSI::Dof_Type::UX )->val();
	    
	    Matrix<DDRMat> tCurPos = tX + tDisp;
	
	    fprintf(stdout," xyx time = %f  x-pos = %f  y-pos = %f  x-new = %f  y-new = %f  x-cur = %f  y-cur = %f  u = %f  v = %f yshift = %f (%f,%f)\n",
	        tTime,tX(0),tX(1),tNewX(0),tNewX(1),tCurPos(0),tCurPos(1),aPropMatrix(0),aPropMatrix(1),
		tYShift, std::min(1.0,std::max(0.0,( tTime - tT1 )/(tT2-tT1))),tC(1)-1.5*tBoxHeight);
	}
    }
    
    void Func_Select_IFC( 
            Matrix< DDRMat >                &aPropMatrix,
            Vector< Matrix< DDRMat > >      &aParameters,
            fem::Field_Interpolator_Manager *aFIManager )
    {
        real tTime = gLogger.get_action_data("NonLinearAlgorithm", "NLBGS", "Solve", "LoadFactor");

        aPropMatrix.set_size( 2, 2, 0.0 );
 
        aPropMatrix( 0, 0 ) = tTime>tT3 ? 0.0: 1.0;
        aPropMatrix( 1, 1 ) = tTime>tT3 ? 0.0: 1.0;
	
        const Matrix< DDRMat > tX = trans(aFIManager->get_IP_geometry_interpolator()->valx());

	const Matrix< DDRMat > tRefPoint = {{tBoxHeight/2.0},{tDomainHeight}};
	
	if ( norm(tX-tRefPoint)<1e-3)
	{
             fprintf(stdout," xyx ifc time = %f  x-pos = %f  y-pos = %f  x-prescribed = %f  y-prescribed = %f\n",
	         tTime,tX(0),tX(1),aPropMatrix( 0, 0 ),aPropMatrix( 1, 1 ));
	}
    }

    void Func_Select_SS3( 
            Matrix< DDRMat >                &aPropMatrix,
            Vector< Matrix< DDRMat > >      &aParameters,
            fem::Field_Interpolator_Manager *aFIManager )
    {
        real tTime = gLogger.get_action_data("NonLinearAlgorithm", "NLBGS", "Solve", "LoadFactor");

        aPropMatrix.set_size( 2, 2, 0.0 );
 
        aPropMatrix( 0, 0 ) =  1.0;
        aPropMatrix( 1, 1 ) =  tTime>tT3 ? 0.0: 1.0;
	
        const Matrix< DDRMat > tX = trans(aFIManager->get_IP_geometry_interpolator()->valx());

	const Matrix< DDRMat > tRefPoint = {{-tBoxHeight/2.0},{tDomainHeight}};
	
	if ( norm(tX-tRefPoint)<1e-3)
	{
             fprintf(stdout," xyx ss3 time = %f  x-pos = %f  y-pos = %f  x-prescribed = %f  y-prescribed = %f\n",
	          tTime,tX(0),tX(1),aPropMatrix( 0, 0 ),aPropMatrix( 1, 1 ));
	}
    }

    void Func_SurfaceLoad_IFC( 
            Matrix< DDRMat >                &aPropMatrix,
            Vector< Matrix< DDRMat > >      &aParameters,
            fem::Field_Interpolator_Manager *aFIManager )
    {
        real tTime = gLogger.get_action_data("NonLinearAlgorithm", "NLBGS", "Solve", "LoadFactor");

        aPropMatrix.set_size( 2, 1, 0.0 );
	
        aPropMatrix(1) = tTime>tT3 ? -1.0e4 : 0.0;
	
	const Matrix< DDRMat > tX = trans(aFIManager->get_IP_geometry_interpolator()->valx());

	const Matrix< DDRMat > tRefPoint = {{tBoxHeight/2.0},{tDomainHeight}};
	
	if ( norm(tX-tRefPoint)<1e-3)
	{
             fprintf(stdout," xyx time = %f  x-pos = %f  y-pos = %f  load-x %f  load-y = %f\n",tTime,tX(0),tX(1),aPropMatrix( 0 ),aPropMatrix( 1 ));
	}
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
        aParameterLists.set( "refinement_buffer", 1 );
        aParameterLists.set( "staircase_buffer", 0 );
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
        aParameterLists.set( "print_enriched_ig_mesh", false );
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
        aParameterLists.set( "number_of_phases", 4 );
        aParameterLists.set( "phase_function_name", F2STR( Phase_Index_Split ) );
        aParameterLists.set( "output_mesh_file","gen.exo");
        /* -------------------------------------------------------------------------------------------- */
        /*                                      Geometry Parameter                                      */
        /* -------------------------------------------------------------------------------------------- */
        auto pl = prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE );
        pl.set( "center_x", -tBoxHeight/2.0 );
        pl.set( "center_y", 0.0);
        pl.set( "normal_x", 1.0);
        pl.set( "normal_y", 0.0);     
        pl.set( "use_multilinear_interpolation", true );
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list(pl);

        pl = prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE );
        pl.set( "center_x", tBoxHeight/2.0 );
        pl.set( "center_y", 0.0);
        pl.set( "normal_x", 1.0);
        pl.set( "normal_y", 0.0);     
        pl.set( "use_multilinear_interpolation", true );
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list(pl);
    
        pl = prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE );
        pl.set( "center_x", 0.0 );
        pl.set( "center_y", tBoxHeight);
        pl.set( "normal_x", 0.0);
        pl.set( "normal_y", 1.0);     
        pl.set( "use_multilinear_interpolation", true );
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list(pl);

        pl = prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE );
        pl.set( "center_x", 0.0 );
        pl.set( "center_y", tUpperCenter - 0.5*tBoxWidth);
        pl.set( "normal_x", 0.0);
        pl.set( "normal_y", 1.0);     
        pl.set( "use_multilinear_interpolation", true );
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list(pl);
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
        pl.set( "phase_indices", "2" );
        aParameterLists( FEM::PHASES ).add_parameter_list( pl );

        pl = prm::create_phase_parameter_list();
        pl.set( "phase_name", "PhaseBottom" );
        pl.set( "phase_indices", "1" );
        aParameterLists( FEM::PHASES ).add_parameter_list( pl );
 
        pl = prm::create_phase_parameter_list();
        pl.set( "phase_name", "PhaseTopRight" );
        pl.set( "phase_indices", "3" );
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
        pl.set( "property_name", "PropYoungsModulus" );
        pl.set( "function_parameters", std::to_string( tEmod ) );
        pl.set( "value_function", F2STR( Func_Const ) );
        aParameterLists( FEM::PROPERTIES ).add_parameter_list( pl );

        /* ------------------------------------------ Poisson ----------------------------------------- */
        pl = prm::create_property_parameter_list();
        pl.set( "property_name", "PropPoisson" );
        pl.set( "function_parameters", std::to_string( tPois ) );
        pl.set( "value_function", F2STR( Func_Const ) );
        aParameterLists( FEM::PROPERTIES ).add_parameter_list( pl );

        /* --------------------------------------- Dirichlet BCs -------------------------------------- */
        pl = prm::create_property_parameter_list();
        pl.set( "property_name", "PropDirichletFixed" );
        pl.set( "function_parameters", "0.0;0.0" );
        pl.set( "value_function", F2STR( Func_Const ) );
        aParameterLists( FEM::PROPERTIES ).add_parameter_list( pl );

        pl = prm::create_property_parameter_list();
        pl.set( "property_name", "PropDirichletIFC" );
        pl.set( "value_function", F2STR( Func_Dirichlet_IFC ) );
        aParameterLists( FEM::PROPERTIES ).add_parameter_list( pl );
 
        pl = prm::create_property_parameter_list();
        pl.set( "property_name", "PropSelectX" );
        pl.set( "function_parameters", "1.0,0.0;0.0,0.0" );
        aParameterLists( FEM::PROPERTIES ).add_parameter_list( pl );

        pl = prm::create_property_parameter_list();
        pl.set( "property_name", "PropSelectY" );
        pl.set( "function_parameters", "0.0,0.0;0.0,1.0" );
        aParameterLists( FEM::PROPERTIES ).add_parameter_list( pl );

        pl = prm::create_property_parameter_list();
        pl.set( "property_name", "PropSelectIFC" );
        pl.set( "value_function", F2STR( Func_Select_IFC ) );
        aParameterLists( FEM::PROPERTIES ).add_parameter_list( pl );

        pl = prm::create_property_parameter_list();
        pl.set( "property_name", "PropSelectSS3" );
        pl.set( "value_function", F2STR( Func_Select_SS3 ) );
        aParameterLists( FEM::PROPERTIES ).add_parameter_list( pl );

        pl = prm::create_property_parameter_list();
        pl.set( "property_name", "PropSurfaceLoadIFC" );
        pl.set( "value_function", F2STR( Func_SurfaceLoad_IFC ) );
        aParameterLists( FEM::PROPERTIES ).add_parameter_list( pl );

        /* -------------------------------------------------------------------------------------------- */
        /*                                       Constitutive Model                                     */
        /* -------------------------------------------------------------------------------------------- */

        /* ----------------------------------- Elastic Material Top ----------------------------------- */
        pl = prm::create_constitutive_model_parameter_list();
        pl.set( "constitutive_name", "MaterialTop" );
        pl.set( "phase_name", "PhaseTop" );
        pl.set( "constitutive_type",
                (uint)fem::Constitutive_Type::STRUC_NON_LIN_ISO_COMPRESSIBLE_NEO_HOOKEAN_WRIGGERS );
        pl.set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
        pl.set( "properties",
                "PropYoungsModulus,YoungsModulus;"
                "PropPoisson,PoissonRatio" );
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
                "PropYoungsModulus,YoungsModulus;"
                "PropPoisson,PoissonRatio" );
        pl.set( "model_type",  fem::Model_Type::PLANE_STRAIN ) ;
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list( pl );

        /* -------------------------------------------------------------------------------------------- */
        /*                                     Stabilization Parameter                                  */
        /* -------------------------------------------------------------------------------------------- */

        pl = prm::create_stabilization_parameter_parameter_list();
        pl.set( "stabilization_name", "SPDirichletNitscheTop" );
        pl.set( "stabilization_type", (uint)fem::Stabilization_Type::DIRICHLET_NITSCHE );
        pl.set( "function_parameters", "100.0" );
        pl.set( "leader_phase_name", "PhaseTop" );
        pl.set( "leader_properties", "PropYoungsModulus,Material" );
        aParameterLists( FEM::STABILIZATION ).add_parameter_list( pl );

        pl = prm::create_stabilization_parameter_parameter_list();
        pl.set( "stabilization_name", "SPDirichletNitscheBottom" );
        pl.set( "stabilization_type", (uint)fem::Stabilization_Type::DIRICHLET_NITSCHE );
        pl.set( "function_parameters", "100.0" );
        pl.set( "leader_phase_name", "PhaseBottom" );
        pl.set( "leader_properties", "PropYoungsModulus,Material" );
        aParameterLists( FEM::STABILIZATION ).add_parameter_list( pl );

        pl = prm::create_stabilization_parameter_parameter_list();
        pl.set( "stabilization_name", "SPContactInterfaceTop" );
        pl.set( "stabilization_type", (uint)fem::Stabilization_Type::DIRICHLET_NITSCHE );
        pl.set( "function_parameters", tContactStabilization );
        pl.set( "leader_phase_name", "PhaseTop" );
        pl.set( "leader_properties", "PropYoungsModulus,Material" );
        aParameterLists( FEM::STABILIZATION ).add_parameter_list( pl );

        pl = prm::create_stabilization_parameter_parameter_list();
        pl.set( "stabilization_name", "SPContactInterfaceBottom" );
        pl.set( "stabilization_type", (uint)fem::Stabilization_Type::DIRICHLET_NITSCHE );
        pl.set( "function_parameters", tContactStabilization );
        pl.set( "leader_phase_name", "PhaseBottom" );
        pl.set( "leader_properties", "PropYoungsModulus,Material" );
        aParameterLists( FEM::STABILIZATION ).add_parameter_list( pl );

        // create ghost penalty displacement Matrix
        pl = prm::create_stabilization_parameter_parameter_list();
        pl.set( "stabilization_name", "SPGPDisplTop" );
        pl.set( "stabilization_type", (uint)fem::Stabilization_Type::GHOST_DISPL );
        pl.set( "function_parameters", "0.05" );
        pl.set( "leader_phase_name", "PhaseTop" );
        pl.set( "follower_phase_name", "PhaseTop" );
        pl.set( "leader_properties", "PropYoungsModulus,Material" );
        aParameterLists( FEM::STABILIZATION ).add_parameter_list( pl );

        // create ghost penalty  displacement Fiber
        pl = prm::create_stabilization_parameter_parameter_list();
        pl.set( "stabilization_name", "SPGPDisplBottom" );
        pl.set( "stabilization_type", (uint)fem::Stabilization_Type::GHOST_DISPL );
        pl.set( "function_parameters", "0.05" );
        pl.set( "leader_phase_name", "PhaseBottom" );
        pl.set( "follower_phase_name", "PhaseBottom" );
        pl.set( "leader_properties", "PropYoungsModulus,Material" );
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
        pl.set( "leader_constitutive_models", "MaterialTop,ElastLinIso" );
        aParameterLists( FEM::IWG ).add_parameter_list( pl );
 
        /* ----------------------------- Body Material Bottom (Left or Whole) ------------------------- */
        pl = prm::create_IWG_parameter_list();
        pl.set( "IWG_name", "IWGBulkStructBottom" );
        pl.set( "IWG_type", (uint)fem::IWG_Type::STRUC_NON_LINEAR_BULK_PF );
        pl.set( "IWG_bulk_type", (uint)fem::Element_Type::BULK );
        pl.set( "dof_residual", "UX,UY" );
        pl.set( "leader_phase_name", "PhaseBottom" );
        pl.set( "leader_constitutive_models", "MaterialBottom,ElastLinIso" );
        aParameterLists( FEM::IWG ).add_parameter_list( pl );

        /* ---------------------------------------- Neumann Top ---------------------------- */
        pl = prm::create_IWG_parameter_list();
        pl.set( "IWG_name", "IWGNeumannIFC" );
        pl.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_NEUMANN );
        pl.set( "IWG_bulk_type", (uint)fem::Element_Type::SIDESET );
        pl.set( "dof_residual", "UX,UY" );
        pl.set( "leader_phase_name", "PhaseTop" );
	pl.set( "neighbor_phases", "PhaseTopRight" );
        pl.set( "leader_properties", "PropSurfaceLoadIFC,Traction" );
        aParameterLists( FEM::IWG ).add_parameter_list( pl );
        
        /* ---------------------------------------- Dirichlet Bottomm --------------------------------- */
        pl = prm::create_IWG_parameter_list();
        pl.set( "IWG_name", "IWGDirichletBottomY" );
        pl.set( "IWG_type", (uint)fem::IWG_Type::STRUC_NON_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE_PF );
        pl.set( "IWG_bulk_type", (uint)fem::Element_Type::SIDESET );
        pl.set( "dof_residual", "UX,UY" );
        pl.set( "leader_phase_name", "PhaseBottom" );
        pl.set( "side_ordinals", "1" );
        pl.set( "leader_properties", "PropDirichletFixed,Dirichlet;PropSelectY,Select" );
        pl.set( "leader_constitutive_models", "MaterialBottom,ElastLinIso" );
        pl.set( "stabilization_parameters", "SPDirichletNitscheTop,DirichletNitsche" );
        aParameterLists( FEM::IWG ).add_parameter_list( pl );

        pl = prm::create_IWG_parameter_list();
        pl.set( "IWG_name", "IWGDirichletBottomX" );
        pl.set( "IWG_type", (uint)fem::IWG_Type::STRUC_NON_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE_PF );
        pl.set( "IWG_bulk_type", (uint)fem::Element_Type::SIDESET );
        pl.set( "dof_residual", "UX,UY" );
        pl.set( "leader_phase_name", "PhaseBottom" );
        pl.set( "side_ordinals", "4" );
        pl.set( "leader_properties", "PropDirichletFixed,Dirichlet;PropSelectX,Select" );
        pl.set( "leader_constitutive_models", "MaterialBottom,ElastLinIso" );
        pl.set( "stabilization_parameters", "SPDirichletNitscheBottom,DirichletNitsche" );
        aParameterLists( FEM::IWG ).add_parameter_list( pl );
	
        /* ---------------------------------------- Dirichlet Top ---------------------------- */
        pl = prm::create_IWG_parameter_list();
        pl.set( "IWG_name", "IWGDirichletTopIFC" );
        pl.set( "IWG_type", (uint)fem::IWG_Type::STRUC_NON_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE_PF );
        pl.set( "IWG_bulk_type", (uint)fem::Element_Type::SIDESET );
        pl.set( "dof_residual", "UX,UY" );
        pl.set( "leader_phase_name", "PhaseTop" );
	pl.set( "neighbor_phases", "PhaseTopRight" );
        pl.set( "leader_properties", "PropDirichletIFC,Dirichlet;PropSelectIFC,Select" );
        pl.set( "leader_constitutive_models", "MaterialTop,ElastLinIso" );
        pl.set( "stabilization_parameters", "SPDirichletNitscheTop,DirichletNitsche" );
        aParameterLists( FEM::IWG ).add_parameter_list( pl );

        pl = prm::create_IWG_parameter_list();
        pl.set( "IWG_name", "IWGDirichletTopSS3" );
        pl.set( "IWG_type", (uint)fem::IWG_Type::STRUC_NON_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE_PF );
        pl.set( "IWG_bulk_type", (uint)fem::Element_Type::SIDESET );
        pl.set( "dof_residual", "UX,UY" );
        pl.set( "leader_phase_name", "PhaseTop" );
        pl.set( "side_ordinals", "3" );
        pl.set( "leader_properties", "PropDirichletIFC,Dirichlet;PropSelectSS3,Select" );
        pl.set( "leader_constitutive_models", "MaterialTop,ElastLinIso" );
        pl.set( "stabilization_parameters", "SPDirichletNitscheTop,DirichletNitsche" );
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
        pl.set( "leader_phase_name", "PhaseTop" );
        pl.set( "neighbor_phases", "PhaseVoid" );
        pl.set( "follower_phase_name", "PhaseBottom" );
        pl.set( "leader_constitutive_models", "MaterialTop,ElastLinIso" );
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
        pl.set( "follower_phase_name", "PhaseTop" );
        pl.set( "leader_constitutive_models", "MaterialBottom,ElastLinIso" );
        pl.set( "follower_constitutive_models", "MaterialTop,ElastLinIso" );
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
            pl.set( "leader_phase_name", "PhaseTop" );
            pl.set( "follower_phase_name", "PhaseTop" );
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
                pl.set( "leader_constitutive_models", "Material" + tLeaderSide + ",ElastLinIso" );
                pl.set( "leader_phase_name", "Phase" + tLeaderSide );
                pl.set( "vectorial_field_index", tDir );
                aParameterLists( FEM::IQI ).add_parameter_list( pl );

                pl = prm::create_IQI_parameter_list();
                pl.set( "IQI_name", "IQIDisp" + tLeaderSide + tAxis );
                pl.set( "IQI_type", (uint)fem::IQI_Type::DOF );
                pl.set( "IQI_bulk_type", (uint)fem::Element_Type::BULK );
                pl.set( "dof_quantity", "UX,UY" );
                pl.set( "leader_phase_name", "Phase" + tLeaderSide );
                pl.set( "vectorial_field_index", tDir );
                aParameterLists( FEM::IQI ).add_parameter_list( pl );
            }

            pl = prm::create_IQI_parameter_list();
            pl.set( "IQI_name", "IQIVonMises" + tLeaderSide );
            pl.set( "IQI_type", (uint)fem::IQI_Type::VON_MISES_STRESS_CAUCHY );
            pl.set( "IQI_bulk_type", (uint)fem::Element_Type::BULK );
            pl.set( "leader_constitutive_models", "Material" + tLeaderSide + ",ElastLinIso" );
            pl.set( "leader_phase_name", "Phase" + tLeaderSide );
            aParameterLists( FEM::IQI ).add_parameter_list( pl );

            /*pl = prm::create_IQI_parameter_list();
            pl.set("IQI_name", "IQIContactPressure" + tLeaderSide);
            pl.set("IQI_type", (uint)fem::IQI_Type::CONTACT_PRESSURE_CURRENT);
            pl.set("IQI_bulk_type", (uint)fem::Element_Type::NONCONFORMAL_SIDESET);
            pl.set("leader_phase_name", "Phase" + tLeaderSide);
            pl.set("neighbor_phases", "PhaseVoid");
            pl.set("follower_phase_name", "Phase" + tFollowerSide);
            pl.set("leader_constitutive_models", "Material" + tLeaderSide + ",TractionCM");
            aParameterLists( FEM::IQI ).add_parameter_list( pl );
            */

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
    void SOLParameterList(Module_Parameter_Lists &aParameterLists) 
    {
      /* -------------------------------------------------------------------------------------------- */
      /*                    Linear Algorithm                      */
      /* -------------------------------------------------------------------------------------------- */
      Parameter_List pl = moris::prm::create_linear_algorithm_parameter_list(sol::SolverType::AMESOS_IMPL);
      // pl.set("Solver_Type", "Amesos_Umfpack");
      aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( pl );
   
      /* -------------------------------------------------------------------------------------------- */
      /*                     Linear Solver                        */
      /* -------------------------------------------------------------------------------------------- */
      pl = moris::prm::create_linear_solver_parameter_list();
      aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list( pl );
      
      /* -------------------------------------------------------------------------------------------- */
      /*                      Nonlinear Algorithm                      */
      /* -------------------------------------------------------------------------------------------- */
      pl = moris::prm::create_nonlinear_algorithm_parameter_list();
      pl.set("NLA_Linear_solver",                   0);
      pl.set("NLA_max_iter",                        tNewtonMaxIter );
      pl.set("NLA_rel_res_norm_drop",               tNewtonRelRes);
      pl.set("NLA_time_offset",                     1.0);
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
      /*                       Nonlinear Solver                       */
      /* -------------------------------------------------------------------------------------------- */
      pl = moris::prm::create_nonlinear_solver_parameter_list();
      pl.set("NLA_DofTypes", "UX,UY");
      aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list( pl );
   
      pl = prm::create_nonlinear_solver_parameter_list();
      pl.set( "NLA_Nonlinear_solver_algorithms", "1" );
      pl.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NLBGS_SOLVER );
      pl.set( "NLA_Sub_Nonlinear_Solver", "0" );
      pl.set( "NLA_DofTypes", "UX,UY" );
      aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list( pl );
   
      /* -------------------------------------------------------------------------------------------- */
      /*                    Time Algorithm                        */
      /* -------------------------------------------------------------------------------------------- */
      pl = moris::prm::create_time_solver_algorithm_parameter_list();
      pl.set("TSA_Num_Time_Steps", 1);
      pl.set("TSA_Time_Frame", 100.0);
      pl.set( "TSA_Nonlinear_Solver", 1 );  
      aParameterLists( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list( pl );
       
      /* -------------------------------------------------------------------------------------------- */
      /*                     Time Solver                          */
      /* -------------------------------------------------------------------------------------------- */
      pl = moris::prm::create_time_solver_parameter_list();
      pl.set("TSA_DofTypes", "UX,UY");
      pl.set("TSA_Output_Indices", "0");
      pl.set("TSA_Output_Criteria", F2STR(Output_Criterion));
      aParameterLists( SOL::TIME_SOLVERS ).add_parameter_list( pl );
   
      /* -------------------------------------------------------------------------------------------- */
      /*                    Solver Warehouse                      */
      /* -------------------------------------------------------------------------------------------- */
      pl = moris::prm::create_solver_warehouse_parameterlist();
      // pl.set("SOL_save_operator_to_matlab", "jacobian");
      aParameterLists( SOL::SOLVER_WAREHOUSE ).add_parameter_list( pl );
   
      /* -------------------------------------------------------------------------------------------- */
      /*                      Preconditioner List                      */
      /* -------------------------------------------------------------------------------------------- */
      pl = moris::prm::create_preconditioner_parameter_list(sol::PreconditionerType::NONE);
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
            tFieldTypes += ",NODAL";    // ELEMENTAL_AVG
            tIQINames += ",IQIVonMises" + tSide;

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
        aParameterLists.set( "Time_Offset", 10.0 );
    }

    void MORISGENERALParameterList( Module_Parameter_Lists &aParameterLists ) {}
}    // namespace moris
#ifdef __cplusplus
}
#endif
