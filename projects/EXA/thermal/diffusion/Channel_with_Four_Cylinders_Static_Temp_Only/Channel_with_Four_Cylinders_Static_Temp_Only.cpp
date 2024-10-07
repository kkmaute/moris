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
#include "parameters.hpp"
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
    OPTParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_opt_problem_parameter_list() );

        aParameterLists( 0 ).set( "is_optimization_problem", false );
    }

    void
    HMRParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_hmr_parameter_list() );

        aParameterLists( 0 ).set( "number_of_elements_per_dimension", "176,88" );
        aParameterLists( 0 ).set( "processor_decomposition_method", 1 );
        aParameterLists( 0 ).set( "processor_dimensions", "2,1" );
        aParameterLists( 0 ).set( "domain_dimensions", "4,2" );
        aParameterLists( 0 ).set( "domain_offset", "-1.24,-0.86" );
        aParameterLists( 0 ).set( "domain_sidesets", "1,2,3,4" );
        aParameterLists( 0 ).set( "lagrange_output_meshes", "0" );

        aParameterLists( 0 ).set( "lagrange_orders", std::to_string( gInterpolationOrder ) );
        aParameterLists( 0 ).set( "lagrange_pattern", "0" );
        aParameterLists( 0 ).set( "bspline_orders", std::to_string( gInterpolationOrder ) );
        aParameterLists( 0 ).set( "bspline_pattern", "0" );

        aParameterLists( 0 ).set( "lagrange_to_bspline", "0" );

        aParameterLists( 0 ).set( "truncate_bsplines", 1 );
        aParameterLists( 0 ).set( "refinement_buffer", (int)gInterpolationOrder );
        aParameterLists( 0 ).set( "staircase_buffer", 1 );
        aParameterLists( 0 ).set( "initial_refinement", "0" );
        aParameterLists( 0 ).set( "initial_refinement_pattern", "0" );

        aParameterLists( 0 ).set( "use_number_aura", 1 );

        aParameterLists( 0 ).set( "use_multigrid", 0 );
        aParameterLists( 0 ).set( "severity_level", 0 );
    }

    void
    XTKParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_xtk_parameter_list() );
        aParameterLists( 0 ).set( "decompose", true );
        aParameterLists( 0 ).set( "decomposition_type", "conformal" );
        aParameterLists( 0 ).set( "enrich", true );
        aParameterLists( 0 ).set( "basis_rank", "bspline" );
        aParameterLists( 0 ).set( "enrich_mesh_indices", "0" );
        aParameterLists( 0 ).set( "ghost_stab", isGhost );
        aParameterLists( 0 ).set( "multigrid", false );
        aParameterLists( 0 ).set( "verbose", true );
        aParameterLists( 0 ).set( "print_enriched_ig_mesh", false );
        aParameterLists( 0 ).set( "exodus_output_XTK_ig_mesh", true );
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
    GENParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists( 0 ).add_parameter_list( prm::create_gen_parameter_list() );
        aParameterLists( 0 ).set( "number_of_phases", 5 );
        aParameterLists( 0 ).set( "phase_function_name", "get_phase_index" );
        aParameterLists( 0 ).set( "output_mesh_file", "GEN_Channel_with_Four_Cylinders_Static_Temp_Only.exo" );

        // Geometry parameter lists
        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists( 1 ).set( "field_function_name", "Func_Bottom_Plane" );
        aParameterLists( 1 ).set( "number_of_refinements", 1 );
        aParameterLists( 1 ).set( "refinement_mesh_index", 0 );

        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists( 1 ).set( "field_function_name", "Func_Top_Plane" );
        aParameterLists( 1 ).set( "number_of_refinements", 1 );
        aParameterLists( 1 ).set( "refinement_mesh_index", 0 );

        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists( 1 ).set( "field_function_name", "Func_Left_Plane" );
        aParameterLists( 1 ).set( "number_of_refinements", 1 );
        aParameterLists( 1 ).set( "refinement_mesh_index", 0 );

        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists( 1 ).set( "field_function_name", "Func_Right_Plane" );
        aParameterLists( 1 ).set( "number_of_refinements", 1 );
        aParameterLists( 1 ).set( "refinement_mesh_index", 0 );

        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists( 1 ).set( "field_function_name", "Func_Cylinder" );
        aParameterLists( 1 ).set( "number_of_refinements", 1 );
        aParameterLists( 1 ).set( "refinement_mesh_index", 0 );

        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists( 1 ).set( "field_function_name", "Func_Cylinder2" );
        aParameterLists( 1 ).set( "number_of_refinements", 1 );
        aParameterLists( 1 ).set( "refinement_mesh_index", 0 );

        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists( 1 ).set( "field_function_name", "Func_Cylinder3" );
        aParameterLists( 1 ).set( "number_of_refinements", 1 );
        aParameterLists( 1 ).set( "refinement_mesh_index", 0 );

        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists( 1 ).set( "field_function_name", "Func_Cylinder4" );
        aParameterLists( 1 ).set( "number_of_refinements", 1 );
        aParameterLists( 1 ).set( "refinement_mesh_index", 0 );
    }

    void
    FEMParameterList( Module_Parameter_Lists& aParameterLists )
    {
        // create a cell of cell of parameter list for fem
        uint tPropIndex = 0;
        uint tCMIndex   = 1;
        uint tSPIndex   = 2;
        uint tIWGIndex  = 3;
        uint tIQIIndex  = 4;
        uint tFEMIndex  = 5;

        uint tPhaseIndex = 7;

        //------------------------------------------------------------------------------

        aParameterLists( tPhaseIndex ).add_parameter_list( prm::create_phase_parameter_list() );
        aParameterLists( tPhaseIndex ).set( "phase_name", "PhaseFluid" );
        aParameterLists( tPhaseIndex ).set( "phase_indices", "0" );

        aParameterLists( tPhaseIndex ).add_parameter_list( prm::create_phase_parameter_list() );
        aParameterLists( tPhaseIndex ).set( "phase_name", "PhaseInlet" );
        aParameterLists( tPhaseIndex ).set( "phase_indices", "1" );

        aParameterLists( tPhaseIndex ).add_parameter_list( prm::create_phase_parameter_list() );
        aParameterLists( tPhaseIndex ).set( "phase_name", "PhaseCylinder" );
        aParameterLists( tPhaseIndex ).set( "phase_indices", "2" );

        aParameterLists( tPhaseIndex ).add_parameter_list( prm::create_phase_parameter_list() );
        aParameterLists( tPhaseIndex ).set( "phase_name", "PhaseWall" );
        aParameterLists( tPhaseIndex ).set( "phase_indices", "3" );

        //------------------------------------------------------------------------------
        // fill the property part of the parameter list

        // create parameter list for property 1
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropDensity" );
        aParameterLists( tPropIndex ).set( "function_parameters", "1.0" );

        // create parameter list for property 2
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropCapacity" );
        aParameterLists( tPropIndex ).set( "function_parameters", "1.0" );

        // create parameter list for property 3
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropConductivity" );
        aParameterLists( tPropIndex ).set( "function_parameters", "0.00005" );

        // create parameter list for property 5
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropInletTemp" );
        aParameterLists( tPropIndex ).set( "function_parameters", "0.0" );

        // create parameter list for property 6
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropSideFlux" );
        aParameterLists( tPropIndex ).set( "function_parameters", "2.0" );

        //------------------------------------------------------------------------------
        // fill the constitutive model part of the parameter list

        // create parameter list for constitutive model 2
        aParameterLists( tCMIndex ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists( tCMIndex ).set( "constitutive_name", "CMDiffusion" );
        aParameterLists( tCMIndex ).set( "phase_name", "PhaseFluid" );
        aParameterLists( tCMIndex ).set( "constitutive_type",  fem::Constitutive_Type::DIFF_LIN_ISO ) ;
        aParameterLists( tCMIndex ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists( tCMIndex ).set( "properties",
                "PropConductivity,Conductivity;"
                "PropDensity     ,Density;"
                "PropCapacity    ,HeatCapacity" );

        //------------------------------------------------------------------------------
        // fill the stabilization parameter part of the parameter list

        // create parameter list for stabilization parameter 2
        aParameterLists( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( tSPIndex ).set( "stabilization_name", "SPDirichletNitscheT" );
        aParameterLists( tSPIndex ).set( "stabilization_type",  fem::Stabilization_Type::DIRICHLET_NITSCHE ) ;
        aParameterLists( tSPIndex ).set( "function_parameters", "100.0" );
        aParameterLists( tSPIndex ).set( "leader_properties", "PropConductivity,Material" );
        aParameterLists( tSPIndex ).set( "leader_phase_name", "PhaseFluid" );

        // create parameter list for stabilization parameter 8
        aParameterLists( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( tSPIndex ).set( "stabilization_name", "SPGPTemp" );
        aParameterLists( tSPIndex ).set( "stabilization_type",  fem::Stabilization_Type::GHOST_DISPL ) ;
        aParameterLists( tSPIndex ).set( "function_parameters", "0.005" );
        aParameterLists( tSPIndex ).set( "leader_properties", "PropConductivity,Material" );
        aParameterLists( tSPIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tSPIndex ).set( "follower_phase_name", "PhaseFluid" );

        //------------------------------------------------------------------------------
        // fill the IWG part of the parameter list

        // create parameter list for IWG 3
        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGDiffusionBulk" );
        aParameterLists( tIWGIndex ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_BULK ) ;
        aParameterLists( tIWGIndex ).set( "dof_residual", "TEMP" );
        aParameterLists( tIWGIndex ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( tIWGIndex ).set( "leader_constitutive_models", "CMDiffusion,Diffusion" );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );

        // create parameter list for IWG 11
        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGInletTemp" );
        aParameterLists( tIWGIndex ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        aParameterLists( tIWGIndex ).set( "dof_residual", "TEMP" );
        aParameterLists( tIWGIndex ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( tIWGIndex ).set( "leader_properties", "PropInletTemp,Dirichlet" );
        aParameterLists( tIWGIndex ).set( "leader_constitutive_models", "CMDiffusion,Diffusion" );
        aParameterLists( tIWGIndex ).set( "stabilization_parameters", "SPDirichletNitscheT,DirichletNitsche" );
        aParameterLists( tIWGIndex ).set( "IWG_bulk_type",  fem::Element_Type::SIDESET ) ;
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIWGIndex ).set( "neighbor_phases", "PhaseInlet" );

        // create parameter list for IWG 11
        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGCylinderFluxTemp" );
        aParameterLists( tIWGIndex ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_NEUMANN ) ;
        aParameterLists( tIWGIndex ).set( "dof_residual", "TEMP" );
        aParameterLists( tIWGIndex ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( tIWGIndex ).set( "leader_properties", "PropSideFlux,Neumann" );
        aParameterLists( tIWGIndex ).set( "IWG_bulk_type",  fem::Element_Type::SIDESET ) ;
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIWGIndex ).set( "neighbor_phases", "PhaseCylinder" );

        if ( isGhost )
        {
            // create parameter list for IWG 16
            aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists( tIWGIndex ).set( "IWG_name", "IWGGPTemp" );
            aParameterLists( tIWGIndex ).set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
            aParameterLists( tIWGIndex ).set( "dof_residual", "TEMP" );
            aParameterLists( tIWGIndex ).set( "leader_dof_dependencies", "TEMP" );
            aParameterLists( tIWGIndex ).set( "follower_dof_dependencies", "TEMP" );
            aParameterLists( tIWGIndex ).set( "stabilization_parameters", "SPGPTemp,GhostSP" );
            aParameterLists( tIWGIndex ).set( "IWG_bulk_type",  fem::Element_Type::DOUBLE_SIDESET ) ;
            aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
            aParameterLists( tIWGIndex ).set( "follower_phase_name", "PhaseFluid" );
            }

        //------------------------------------------------------------------------------
        // fill the IQI part of the parameter list

        // create parameter list for IQI 3
        aParameterLists( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIBulkTEMP" );
        aParameterLists( tIQIIndex ).set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists( tIQIIndex ).set( "dof_quantity", "TEMP" );
        aParameterLists( tIQIIndex ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( tIQIIndex ).set( "vectorial_field_index", 0 );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        aParameterLists( tFEMIndex ).add_parameter_list( prm::create_computation_parameter_list() );
    }

    void
    SOLParameterList( Module_Parameter_Lists& aParameterLists )
    {

        if ( gTestCaseIndex == 0 )
        {
            aParameterLists( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::PETSC ) );
        }
        else
        {
            aParameterLists( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL ) );
        }

        aParameterLists( 1 ).add_parameter_list( moris::prm::create_linear_solver_parameter_list() );

        aParameterLists( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
        aParameterLists( 2 ).set( "NLA_rel_res_norm_drop", 1e-04 );
        aParameterLists( 2 ).set( "NLA_relaxation_parameter", 1.0 );
        aParameterLists( 2 ).set( "NLA_max_iter", 20 );

        aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
        aParameterLists( 3 ).set( "NLA_DofTypes", "TEMP" );

        aParameterLists( 4 ).add_parameter_list( moris::prm::create_time_solver_algorithm_parameter_list() );

        aParameterLists( 5 ).add_parameter_list( moris::prm::create_time_solver_parameter_list() );
        aParameterLists( 5 ).set( "TSA_DofTypes", "TEMP" );
        aParameterLists( 5 ).set( "TSA_Initialize_Sol_Vec", "TEMP,0.0" );
        aParameterLists( 5 ).set( "TSA_Output_Indices", "0" );
        aParameterLists( 5 ).set( "TSA_Output_Criteria", "Output_Criterion" );

        aParameterLists( 6 ).add_parameter_list( moris::prm::create_solver_warehouse_parameterlist() );
        if ( gTestCaseIndex == 0 )
        {
            aParameterLists( 6 ).set( "SOL_TPL_Type",  sol::MapType::Petsc ) ;
        }

        aParameterLists( 7 ).add_parameter_list( moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::NONE ) );
    }

    void
    MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_msi_parameter_list() );
        aParameterLists( 0 ).set( "order_adofs_by_host", false );
    }

    void
    VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_vis_parameter_list() );
        aParameterLists( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", "Channel_with_Four_Cylinders_Static_Temp_Only.exo" ) );
        aParameterLists( 0 ).set( "Mesh_Type",  vis::VIS_Mesh_Type::STANDARD ) ;
        aParameterLists( 0 ).set( "Set_Names", "HMR_dummy_n_p0,HMR_dummy_c_p0" );
        aParameterLists( 0 ).set( "Field_Names", "TEMP,IQIBulkTEMP" );
        aParameterLists( 0 ).set( "Field_Type", "NODAL,GLOBAL" );
        aParameterLists( 0 ).set( "IQI_Names", "IQIBulkTEMP,IQIBulkTEMP" );
        aParameterLists( 0 ).set( "Save_Frequency", 1 );
    }

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
