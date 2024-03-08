/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Channel_with_Four_Cylinders_Static_Temp_Only.cpp
 *
 */

#include <string>
#include <iostream>
#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "cl_Bitset.hpp"
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

extern uint gTestCaseIndex;

#ifdef __cplusplus
extern "C" {
#endif
//------------------------------------------------------------------------------
namespace moris
{
    bool isGhost = true;

    // Minimum value to be returned by level set function
    moris::real tMinLSvalue = 1.0e-6;

    // Geometry Parameters
    moris::real tPlaneBottom = 0.0;  /* y bottom plane (m) */
    moris::real tPlaneTop    = 0.41; /* y top plane    (m) */
    moris::real tPlaneLeft   = 0.0;  /* x left plane   (m) */
    moris::real tPlaneRight  = 2.2;  /* x right plane  (m) */

    moris::real tCylinderCenterX = 0.2;
    moris::real tCylinderCenterY = 0.2;
    moris::real tCylinderRadius  = 0.05;
    moris::real tCylinderOffset  = 0.10;

    // Constant function for properties
    void
    Func_Const(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    // Inlet velocity function
    void
    Func_Inlet_U(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix.set_size( 2, 1, 0.0 );
        real tY          = aFIManager->get_IP_geometry_interpolator()->valx()( 1 );
        aPropMatrix( 0 ) = 4.0 * aParameters( 0 )( 0 ) * tY * ( 0.41 - tY ) / ( std::pow( 0.41, 2.0 ) );
    }

    bool
    Output_Criterion( moris::tsa::Time_Solver* aTimeSolver )
    {
        return true;
    }

    moris::real Func_Bottom_Plane(
            const moris::Matrix< DDRMat >& aCoordinates,
            const Vector< real >&     aGeometryParameters )
    {
        moris::real tValue = aCoordinates( 1 ) - tPlaneBottom;

        return std::abs( tValue ) < tMinLSvalue ? tMinLSvalue : tValue;
    }

    moris::real Func_Top_Plane(
            const moris::Matrix< DDRMat >& aCoordinates,
            const Vector< real >&     aGeometryParameters )
    {
        moris::real tValue = aCoordinates( 1 ) - tPlaneTop;

        return std::abs( tValue ) < tMinLSvalue ? tMinLSvalue : tValue;
    }

    moris::real Func_Left_Plane(
            const moris::Matrix< DDRMat >& aCoordinates,
            const Vector< real >&     aGeometryParameters )
    {
        moris::real tValue = aCoordinates( 0 ) - tPlaneLeft;

        return std::abs( tValue ) < tMinLSvalue ? tMinLSvalue : tValue;
    }

    moris::real Func_Right_Plane(
            const moris::Matrix< DDRMat >& aCoordinates,
            const Vector< real >&     aGeometryParameters )
    {
        moris::real tValue = aCoordinates( 0 ) - tPlaneRight;

        return std::abs( tValue ) < tMinLSvalue ? tMinLSvalue : tValue;
    }

    moris::real Func_Cylinder(
            const moris::Matrix< DDRMat >& aCoordinates,
            const Vector< real >&     aGeometryParameters )
    {
        moris::real tValue = tCylinderRadius - std::pow( std::pow( aCoordinates( 0 ) - 0.4, 2.0 ) + std::pow( aCoordinates( 1 ) - ( 0.2 - tCylinderOffset ), 2.0 ), 0.5 );

        return std::abs( tValue ) < tMinLSvalue ? tMinLSvalue : tValue;
    }

    moris::real Func_Cylinder2(
            const moris::Matrix< DDRMat >& aCoordinates,
            const Vector< real >&     aGeometryParameters )
    {
        moris::real tValue = tCylinderRadius - std::pow( std::pow( aCoordinates( 0 ) - 0.8, 2.0 ) + std::pow( aCoordinates( 1 ) - ( 0.2 + tCylinderOffset ), 2.0 ), 0.5 );

        return std::abs( tValue ) < tMinLSvalue ? tMinLSvalue : tValue;
    }

    moris::real Func_Cylinder3(
            const moris::Matrix< DDRMat >& aCoordinates,
            const Vector< real >&     aGeometryParameters )
    {
        moris::real tValue = tCylinderRadius - std::pow( std::pow( aCoordinates( 0 ) - 1.2, 2.0 ) + std::pow( aCoordinates( 1 ) - ( 0.2 - tCylinderOffset ), 2.0 ), 0.5 );

        return std::abs( tValue ) < tMinLSvalue ? tMinLSvalue : tValue;
    }

    moris::real Func_Cylinder4(
            const moris::Matrix< DDRMat >& aCoordinates,
            const Vector< real >&     aGeometryParameters )
    {
        moris::real tValue = tCylinderRadius - std::pow( std::pow( aCoordinates( 0 ) - 1.6, 2.0 ) + std::pow( aCoordinates( 1 ) - ( 0.2 + tCylinderOffset ), 2.0 ), 0.5 );

        return std::abs( tValue ) < tMinLSvalue ? tMinLSvalue : tValue;
    }

    void
    OPTParameterList( Vector< Vector< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_opt_problem_parameter_list();

        tParameterlist( 0 )( 0 ).set( "is_optimization_problem", false );
    }

    void
    HMRParameterList( Vector< Vector< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_hmr_parameter_list();

        tParameterlist( 0 )( 0 ).set( "number_of_elements_per_dimension", "176,88" );
        tParameterlist( 0 )( 0 ).set( "processor_decomposition_method", 1 );
        tParameterlist( 0 )( 0 ).set( "processor_dimensions", "2,1" );
        tParameterlist( 0 )( 0 ).set( "domain_dimensions", "4,2" );
        tParameterlist( 0 )( 0 ).set( "domain_offset", "-1.24,-0.86" );
        tParameterlist( 0 )( 0 ).set( "domain_sidesets", "1,2,3,4" );
        tParameterlist( 0 )( 0 ).set( "lagrange_output_meshes", "0" );

        tParameterlist( 0 )( 0 ).set( "lagrange_orders", std::to_string( gInterpolationOrder ) );
        tParameterlist( 0 )( 0 ).set( "lagrange_pattern", "0" );
        tParameterlist( 0 )( 0 ).set( "bspline_orders", std::to_string( gInterpolationOrder ) );
        tParameterlist( 0 )( 0 ).set( "bspline_pattern", "0" );

        tParameterlist( 0 )( 0 ).set( "lagrange_to_bspline", "0" );

        tParameterlist( 0 )( 0 ).set( "truncate_bsplines", 1 );
        tParameterlist( 0 )( 0 ).set( "refinement_buffer", (int)gInterpolationOrder );
        tParameterlist( 0 )( 0 ).set( "staircase_buffer", 1 );
        tParameterlist( 0 )( 0 ).set( "initial_refinement", "0" );
        tParameterlist( 0 )( 0 ).set( "initial_refinement_pattern", "0" );

        tParameterlist( 0 )( 0 ).set( "use_number_aura", 1 );

        tParameterlist( 0 )( 0 ).set( "use_multigrid", 0 );
        tParameterlist( 0 )( 0 ).set( "severity_level", 0 );
    }

    void
    XTKParameterList( Vector< Vector< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_xtk_parameter_list();
        tParameterlist( 0 )( 0 ).set( "decompose", true );
        tParameterlist( 0 )( 0 ).set( "decomposition_type", "conformal" );
        tParameterlist( 0 )( 0 ).set( "enrich", true );
        tParameterlist( 0 )( 0 ).set( "basis_rank", "bspline" );
        tParameterlist( 0 )( 0 ).set( "enrich_mesh_indices", "0" );
        tParameterlist( 0 )( 0 ).set( "ghost_stab", isGhost );
        tParameterlist( 0 )( 0 ).set( "multigrid", false );
        tParameterlist( 0 )( 0 ).set( "verbose", true );
        tParameterlist( 0 )( 0 ).set( "print_enriched_ig_mesh", false );
        tParameterlist( 0 )( 0 ).set( "exodus_output_XTK_ig_mesh", true );
    }

    uint
    get_phase_index( const Bitset< 8 >& aGeometrySigns )
    {
        // phase table
        // 0 - fluid
        // 1 - void inlet
        // 2 - cylinder
        // 3 - void walls
        // 4 - void outlet (not used)

        // in between lower and upper planes
        if ( aGeometrySigns.test( 0 ) && !aGeometrySigns.test( 1 ) )
        {
            // void inlet
            if ( !aGeometrySigns.test( 2 ) )
            {
                return 1;
            }
            // void outlet
            if ( aGeometrySigns.test( 3 ) )
            {
                return 4;
            }
            // cylinder
            if (                                                               //
                    aGeometrySigns.test( 4 ) || aGeometrySigns.test( 5 ) ||    //
                    aGeometrySigns.test( 6 ) || aGeometrySigns.test( 7 ) )
            {
                return 2;
            }
            // fluid
            return 0;
        }

        // void walls
        return 3;
    }

    void
    GENParameterList( Vector< Vector< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 3 );
        tParameterlist( 0 ).resize( 1 );
        tParameterlist( 1 ).resize( 8 );

        tParameterlist( 0 )( 0 ) = prm::create_gen_parameter_list();
        tParameterlist( 0 )( 0 ).set( "number_of_phases", 5 );
        tParameterlist( 0 )( 0 ).set( "phase_function_name", "get_phase_index" );
        tParameterlist( 0 )( 0 ).set( "output_mesh_file", "GEN_Channel_with_Four_Cylinders_Static_Temp_Only.exo" );

        // Geometry parameter lists
        tParameterlist( 1 )( 0 ) = prm::create_user_defined_geometry_parameter_list();
        tParameterlist( 1 )( 0 ).set( "field_function_name", "Func_Bottom_Plane" );
        tParameterlist( 1 )( 0 ).set( "number_of_refinements", "1" );
        tParameterlist( 1 )( 0 ).set( "refinement_mesh_index", "0" );

        tParameterlist( 1 )( 1 ) = prm::create_user_defined_geometry_parameter_list();
        tParameterlist( 1 )( 1 ).set( "field_function_name", "Func_Top_Plane" );
        tParameterlist( 1 )( 1 ).set( "number_of_refinements", "1" );
        tParameterlist( 1 )( 1 ).set( "refinement_mesh_index", "0" );

        tParameterlist( 1 )( 2 ) = prm::create_user_defined_geometry_parameter_list();
        tParameterlist( 1 )( 2 ).set( "field_function_name", "Func_Left_Plane" );
        tParameterlist( 1 )( 2 ).set( "number_of_refinements", "1" );
        tParameterlist( 1 )( 2 ).set( "refinement_mesh_index", "0" );

        tParameterlist( 1 )( 3 ) = prm::create_user_defined_geometry_parameter_list();
        tParameterlist( 1 )( 3 ).set( "field_function_name", "Func_Right_Plane" );
        tParameterlist( 1 )( 3 ).set( "number_of_refinements", "1" );
        tParameterlist( 1 )( 3 ).set( "refinement_mesh_index", "0" );

        tParameterlist( 1 )( 4 ) = prm::create_user_defined_geometry_parameter_list();
        tParameterlist( 1 )( 4 ).set( "field_function_name", "Func_Cylinder" );
        tParameterlist( 1 )( 4 ).set( "number_of_refinements", "1" );
        tParameterlist( 1 )( 4 ).set( "refinement_mesh_index", "0" );

        tParameterlist( 1 )( 5 ) = prm::create_user_defined_geometry_parameter_list();
        tParameterlist( 1 )( 5 ).set( "field_function_name", "Func_Cylinder2" );
        tParameterlist( 1 )( 5 ).set( "number_of_refinements", "1" );
        tParameterlist( 1 )( 5 ).set( "refinement_mesh_index", "0" );

        tParameterlist( 1 )( 6 ) = prm::create_user_defined_geometry_parameter_list();
        tParameterlist( 1 )( 6 ).set( "field_function_name", "Func_Cylinder3" );
        tParameterlist( 1 )( 6 ).set( "number_of_refinements", "1" );
        tParameterlist( 1 )( 6 ).set( "refinement_mesh_index", "0" );

        tParameterlist( 1 )( 7 ) = prm::create_user_defined_geometry_parameter_list();
        tParameterlist( 1 )( 7 ).set( "field_function_name", "Func_Cylinder4" );
        tParameterlist( 1 )( 7 ).set( "number_of_refinements", "1" );
        tParameterlist( 1 )( 7 ).set( "refinement_mesh_index", "0" );
    }

    void
    FEMParameterList( Vector< Vector< ParameterList > >& tParameterList )
    {
        // create a cell of cell of parameter list for fem
        tParameterList.resize( 9 );
        uint tPropIndex = 0;
        uint tCMIndex   = 1;
        uint tSPIndex   = 2;
        uint tIWGIndex  = 3;
        uint tIQIIndex  = 4;
        uint tFEMIndex  = 5;

        uint tPhaseIndex = 7;

        //------------------------------------------------------------------------------
        // phase info
        uint tPhaseCounter = 0;

        tParameterList( tPhaseIndex ).push_back( prm::create_phase_parameter_list() );
        tParameterList( tPhaseIndex )( tPhaseCounter ).set( "phase_name", "PhaseFluid" );
        tParameterList( tPhaseIndex )( tPhaseCounter ).set( "phase_indices", "0" );
        tPhaseCounter++;

        tParameterList( tPhaseIndex ).push_back( prm::create_phase_parameter_list() );
        tParameterList( tPhaseIndex )( tPhaseCounter ).set( "phase_name", "PhaseInlet" );
        tParameterList( tPhaseIndex )( tPhaseCounter ).set( "phase_indices", "1" );
        tPhaseCounter++;

        tParameterList( tPhaseIndex ).push_back( prm::create_phase_parameter_list() );
        tParameterList( tPhaseIndex )( tPhaseCounter ).set( "phase_name", "PhaseCylinder" );
        tParameterList( tPhaseIndex )( tPhaseCounter ).set( "phase_indices", "2" );
        tPhaseCounter++;

        tParameterList( tPhaseIndex ).push_back( prm::create_phase_parameter_list() );
        tParameterList( tPhaseIndex )( tPhaseCounter ).set( "phase_name", "PhaseWall" );
        tParameterList( tPhaseIndex )( tPhaseCounter ).set( "phase_indices", "3" );
        tPhaseCounter++;

        //------------------------------------------------------------------------------
        // fill the property part of the parameter list

        // init property counter
        uint tPropCounter = 0;

        // create parameter list for property 1
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropDensity" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", "1.0" );
        tPropCounter++;

        // create parameter list for property 2
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropCapacity" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", "1.0" );
        tPropCounter++;

        // create parameter list for property 3
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropConductivity" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", "0.00005" );
        tPropCounter++;

        // create parameter list for property 5
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropInletTemp" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", "0.0" );
        tPropCounter++;

        // create parameter list for property 6
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropSideFlux" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", "2.0" );
        tPropCounter++;

        //------------------------------------------------------------------------------
        // fill the constitutive model part of the parameter list

        // init CM counter
        uint tCMCounter = 0;

        // create parameter list for constitutive model 2
        tParameterList( tCMIndex ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( tCMIndex )( tCMCounter ).set( "constitutive_name", "CMDiffusion" );
        tParameterList( tCMIndex )( tCMCounter ).set( "phase_name", "PhaseFluid" );
        tParameterList( tCMIndex )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::DIFF_LIN_ISO ) );
        tParameterList( tCMIndex )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( tCMIndex )( tCMCounter ).set( "properties",
                "PropConductivity,Conductivity;"
                "PropDensity     ,Density;"
                "PropCapacity    ,HeatCapacity" );
        tCMCounter++;

        //------------------------------------------------------------------------------
        // fill the stabilization parameter part of the parameter list

        // init SP counter
        uint tSPCounter = 0;

        // create parameter list for stabilization parameter 2
        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name", "SPDirichletNitscheT" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::DIRICHLET_NITSCHE ) );
        tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters", "100.0" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_properties", "PropConductivity,Material" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_phase_name", "PhaseFluid" );
        tSPCounter++;

        // create parameter list for stabilization parameter 8
        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name", "SPGPTemp" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::GHOST_DISPL ) );
        tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters", "0.005" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_properties", "PropConductivity,Material" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tSPIndex )( tSPCounter ).set( "follower_phase_name", "PhaseFluid" );
        tSPCounter++;

        //------------------------------------------------------------------------------
        // fill the IWG part of the parameter list

        // init IWG counter
        uint tIWGCounter = 0;

        // create parameter list for IWG 3
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGDiffusionBulk" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::SPATIALDIFF_BULK ) );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "TEMP" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", "CMDiffusion,Diffusion" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFluid" );
        tIWGCounter++;

        // create parameter list for IWG 11
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGInletTemp" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::SPATIALDIFF_DIRICHLET_SYMMETRIC_NITSCHE ) );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "TEMP" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_properties", "PropInletTemp,Dirichlet" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", "CMDiffusion,Diffusion" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SPDirichletNitscheT,DirichletNitsche" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", static_cast< uint >( fem::Element_Type::SIDESET ) );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "neighbor_phases", "PhaseInlet" );
        tIWGCounter++;

        // create parameter list for IWG 11
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGCylinderFluxTemp" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::SPATIALDIFF_NEUMANN ) );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "TEMP" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_properties", "PropSideFlux,Neumann" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", static_cast< uint >( fem::Element_Type::SIDESET ) );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "neighbor_phases", "PhaseCylinder" );
        tIWGCounter++;

        if ( isGhost )
        {
            // create parameter list for IWG 16
            tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGGPTemp" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::GHOST_NORMAL_FIELD ) );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "TEMP" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_dof_dependencies", "TEMP" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "follower_dof_dependencies", "TEMP" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SPGPTemp,GhostSP" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", static_cast< uint >( fem::Element_Type::DOUBLE_SIDESET ) );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFluid" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "follower_phase_name", "PhaseFluid" );
            tIWGCounter++;
        }

        //------------------------------------------------------------------------------
        // fill the IQI part of the parameter list

        // init IQI counter
        uint tIQICounter = 0;

        // create parameter list for IQI 3
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIBulkTEMP" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::DOF ) );
        tParameterList( tIQIIndex )( tIQICounter ).set( "dof_quantity", "TEMP" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index", 0 );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseFluid" );
        tIQICounter++;

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        tParameterList( tFEMIndex ).resize( 1 );
        tParameterList( tFEMIndex )( 0 ) = prm::create_computation_parameter_list();
    }

    void
    SOLParameterList( Vector< Vector< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 8 );
        for ( uint Ik = 0; Ik < 8; Ik++ )
        {
            tParameterlist( Ik ).resize( 1 );
        }

        if ( gTestCaseIndex == 0 )
        {
            tParameterlist( 0 )( 0 ) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::PETSC );
        }
        else
        {
            tParameterlist( 0 )( 0 ) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL );
        }

        tParameterlist( 1 )( 0 ) = moris::prm::create_linear_solver_parameter_list();

        tParameterlist( 2 )( 0 ) = moris::prm::create_nonlinear_algorithm_parameter_list();
        tParameterlist( 2 )( 0 ).set( "NLA_rel_res_norm_drop", 1e-04 );
        tParameterlist( 2 )( 0 ).set( "NLA_relaxation_parameter", 1.0 );
        tParameterlist( 2 )( 0 ).set( "NLA_max_iter", 20 );

        tParameterlist( 3 )( 0 ) = moris::prm::create_nonlinear_solver_parameter_list();
        tParameterlist( 3 )( 0 ).set( "NLA_DofTypes", "TEMP" );

        tParameterlist( 4 )( 0 ) = moris::prm::create_time_solver_algorithm_parameter_list();

        tParameterlist( 5 )( 0 ) = moris::prm::create_time_solver_parameter_list();
        tParameterlist( 5 )( 0 ).set( "TSA_DofTypes", "TEMP" );
        tParameterlist( 5 )( 0 ).set( "TSA_Initialize_Sol_Vec", "TEMP,0.0" );
        tParameterlist( 5 )( 0 ).set( "TSA_Output_Indices", "0" );
        tParameterlist( 5 )( 0 ).set( "TSA_Output_Criteria", "Output_Criterion" );

        tParameterlist( 6 )( 0 ) = moris::prm::create_solver_warehouse_parameterlist();
        if ( gTestCaseIndex == 0 )
        {
            tParameterlist( 6 )( 0 ).set( "SOL_TPL_Type", static_cast< uint >( sol::MapType::Petsc ) );
        }

        tParameterlist( 7 )( 0 ) = moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::NONE );
    }

    void
    MSIParameterList( Vector< Vector< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_msi_parameter_list();
        tParameterlist( 0 )( 0 ).set( "order_adofs_by_host", false );
    }

    void
    VISParameterList( Vector< Vector< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_vis_parameter_list();
        tParameterlist( 0 )( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", "Channel_with_Four_Cylinders_Static_Temp_Only.exo" ) );
        tParameterlist( 0 )( 0 ).set( "Mesh_Type", static_cast< uint >( vis::VIS_Mesh_Type::STANDARD ) );
        tParameterlist( 0 )( 0 ).set( "Set_Names", "HMR_dummy_n_p0,HMR_dummy_c_p0" );
        tParameterlist( 0 )( 0 ).set( "Field_Names", "TEMP,IQIBulkTEMP" );
        tParameterlist( 0 )( 0 ).set( "Field_Type", "NODAL,GLOBAL" );
        tParameterlist( 0 )( 0 ).set( "IQI_Names", "IQIBulkTEMP,IQIBulkTEMP" );
        tParameterlist( 0 )( 0 ).set( "Save_Frequency", 1 );
    }

    void
    MORISGENERALParameterList( Vector< Vector< ParameterList > >& tParameterlist )
    {
    }

    //------------------------------------------------------------------------------
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
