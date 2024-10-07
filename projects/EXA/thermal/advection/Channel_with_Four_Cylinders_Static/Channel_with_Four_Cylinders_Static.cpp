/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Channel_with_Four_Cylinders_Static.cpp
 *
 */

#include <string>
#include <iostream>
#include "moris_typedefs.hpp"
#include "cl_Bitset.hpp"
#include "cl_Matrix.hpp"
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

#ifdef __cplusplus
extern "C" {
#endif
//------------------------------------------------------------------------------
namespace moris
{
    // Files
    std::string tFileName      = "Channel_with_Four_Cylinders_Static";
    std::string tExoFile       = tFileName + ".exo";
    std::string tGENOutputFile = tFileName + "_GEN.exo";

    // use ghost
    bool isGhost = true;

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
        aParameterLists( 0 ).set( "refinement_buffer", 1 );
        aParameterLists( 0 ).set( "staircase_buffer", 1 );
        aParameterLists( 0 ).set( "initial_refinement", "1" );
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
            if (                                                                 //
                    !aGeometrySigns.test( 4 ) || !aGeometrySigns.test( 5 ) ||    //
                    !aGeometrySigns.test( 6 ) || !aGeometrySigns.test( 7 ) )
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
        aParameterLists( 0 ).set( "output_mesh_file", tGENOutputFile );

        aParameterLists( 0 ).set( "number_of_phases", 5 );
        aParameterLists( 0 ).set( "phase_function_name", "get_phase_index" );

        // Bottom plane
        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
        aParameterLists( 1 ).set( "center_x", 0.0 );
        aParameterLists( 1 ).set( "center_y", 0.0 );
        aParameterLists( 1 ).set( "normal_x", 0.0 );
        aParameterLists( 1 ).set( "normal_y", 1.0 );
        aParameterLists( 1 ).set( "number_of_refinements", 1 );
        aParameterLists( 1 ).set( "refinement_mesh_index", 0 );

        // Top plane
        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
        aParameterLists( 1 ).set( "center_x", 0.0 );
        aParameterLists( 1 ).set( "center_y", 0.41 );
        aParameterLists( 1 ).set( "normal_x", 0.0 );
        aParameterLists( 1 ).set( "normal_y", 1.0 );
        aParameterLists( 1 ).set( "number_of_refinements", 1 );
        aParameterLists( 1 ).set( "refinement_mesh_index", 0 );

        // Left plane
        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
        aParameterLists( 1 ).set( "center_x", 0.0 );
        aParameterLists( 1 ).set( "center_y", 0.0 );
        aParameterLists( 1 ).set( "normal_x", 1.0 );
        aParameterLists( 1 ).set( "normal_y", 0.0 );
        aParameterLists( 1 ).set( "number_of_refinements", 1 );
        aParameterLists( 1 ).set( "refinement_mesh_index", 0 );

        // Right plane
        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
        aParameterLists( 1 ).set( "center_x", 2.2 );
        aParameterLists( 1 ).set( "center_y", 0.0 );
        aParameterLists( 1 ).set( "normal_x", 1.0 );
        aParameterLists( 1 ).set( "normal_y", 0.0 );
        aParameterLists( 1 ).set( "number_of_refinements", 1 );
        aParameterLists( 1 ).set( "refinement_mesh_index", 0 );

        // Cylinder 1
        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::CIRCLE ) );
        aParameterLists( 1 ).set( "center_x", 0.4 );
        aParameterLists( 1 ).set( "center_y", 0.1 );
        aParameterLists( 1 ).set( "radius", 0.05 );
        aParameterLists( 1 ).set( "number_of_refinements", 1 );
        aParameterLists( 1 ).set( "refinement_mesh_index", 0 );

        // Cylinder 2
        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::CIRCLE ) );
        aParameterLists( 1 ).set( "center_x", 0.8 );
        aParameterLists( 1 ).set( "center_y", 0.3 );
        aParameterLists( 1 ).set( "radius", 0.05 );
        aParameterLists( 1 ).set( "number_of_refinements", 1 );
        aParameterLists( 1 ).set( "refinement_mesh_index", 0 );

        // Cylinder 3
        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::CIRCLE ) );
        aParameterLists( 1 ).set( "center_x", 1.2 );
        aParameterLists( 1 ).set( "center_y", 0.1 );
        aParameterLists( 1 ).set( "radius", 0.05 );
        aParameterLists( 1 ).set( "number_of_refinements", 1 );
        aParameterLists( 1 ).set( "refinement_mesh_index", 0 );

        // Cylinder 4
        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::CIRCLE ) );
        aParameterLists( 1 ).set( "center_x", 1.6 );
        aParameterLists( 1 ).set( "center_y", 0.3 );
        aParameterLists( 1 ).set( "radius", 0.05 );
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

        // fluid properties ------------------------------------------------------------
        // viscosity
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropViscosity" );
        aParameterLists( tPropIndex ).set( "function_parameters", "0.001" );

        // density
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropDensity" );
        aParameterLists( tPropIndex ).set( "function_parameters", "1.0" );

        // capacity
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropCapacity" );
        aParameterLists( tPropIndex ).set( "function_parameters", "1.0" );

        // conductivity
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropConductivity" );
        aParameterLists( tPropIndex ).set( "function_parameters", "0.00005" );

        // BC properties ---------------------------------------------------------------
        // inlet velocity
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropDirichletInletU" );
        aParameterLists( tPropIndex ).set( "function_parameters", "1.5" );
        aParameterLists( tPropIndex ).set( "value_function", "Func_Inlet_U" );

        // zero velocity
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropDirichletZeroU" );
        aParameterLists( tPropIndex ).set( "function_parameters", "0.0;0.0" );

        // inlet temperature
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropInletTemp" );
        aParameterLists( tPropIndex ).set( "function_parameters", "0.0" );

        // side flux on cylinders
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropSideFlux" );
        aParameterLists( tPropIndex ).set( "function_parameters", "2.0" );

        //------------------------------------------------------------------------------
        // fill the constitutive model part of the parameter list

        // Fluid CM --------------------------------------------------------------------
        // NS incompressible
        aParameterLists( tCMIndex ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists( tCMIndex ).set( "constitutive_name", "CMFluid" );
        aParameterLists( tCMIndex ).set( "phase_name", "PhaseFluid" );
        aParameterLists( tCMIndex ).set( "constitutive_type",  fem::Constitutive_Type::FLUID_INCOMPRESSIBLE ) ;
        aParameterLists( tCMIndex ).set( "dof_dependencies", std::pair< std::string, std::string >( "VX,VY;P", "Velocity,Pressure" ) );
        aParameterLists( tCMIndex ).set( "properties",
                "PropViscosity,Viscosity;"
                "PropDensity  ,Density" );

        // diffusion
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

        // SUPG/PSPG fluid incompressible (note: use time step in tau although static; for testing purposes only)
        aParameterLists( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( tSPIndex ).set( "stabilization_name", "SPIncFlow" );
        aParameterLists( tSPIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tSPIndex ).set( "stabilization_type",  fem::Stabilization_Type::INCOMPRESSIBLE_FLOW ) ;
        aParameterLists( tSPIndex ).set( "function_parameters", "36.0/1.0" );
        aParameterLists( tSPIndex ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY;P", "Velocity,Pressure" ) );
        aParameterLists( tSPIndex ).set( "leader_properties",
                "PropViscosity,Viscosity;"
                "PropDensity,Density" );

        // SUPG advection (note: use time step in tau although static; for testing purposes only)
        aParameterLists( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( tSPIndex ).set( "stabilization_name", "SPSUPGTemp" );
        aParameterLists( tSPIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tSPIndex ).set( "stabilization_type",  fem::Stabilization_Type::SUPG_ADVECTION ) ;
        aParameterLists( tSPIndex ).set( "function_parameters", "1.0" );
        aParameterLists( tSPIndex ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
        aParameterLists( tSPIndex ).set( "leader_properties", "PropConductivity,Conductivity" );

        // Dirichlet Nitsche for velocity
        aParameterLists( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( tSPIndex ).set( "stabilization_name", "SPDirichletNitscheU" );
        aParameterLists( tSPIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tSPIndex ).set( "stabilization_type",  fem::Stabilization_Type::VELOCITY_DIRICHLET_NITSCHE ) ;
        aParameterLists( tSPIndex ).set( "function_parameters", "100.0/1.0" );
        aParameterLists( tSPIndex ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
        aParameterLists( tSPIndex ).set( "leader_properties",
                "PropViscosity,Viscosity;"
                "PropDensity,Density" );

        // Dirichhlet Nitsche for temperature
        aParameterLists( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( tSPIndex ).set( "stabilization_name", "SPDirichletNitscheT" );
        aParameterLists( tSPIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tSPIndex ).set( "stabilization_type",  fem::Stabilization_Type::DIRICHLET_NITSCHE ) ;
        aParameterLists( tSPIndex ).set( "function_parameters", "100.0" );
        aParameterLists( tSPIndex ).set( "leader_properties", "PropConductivity,Material" );

        if ( isGhost )
        {
            // ghost viscous
            aParameterLists( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
            aParameterLists( tSPIndex ).set( "stabilization_name", "SPGPViscosity" );
            aParameterLists( tSPIndex ).set( "leader_phase_name", "PhaseFluid" );
            aParameterLists( tSPIndex ).set( "stabilization_type",  fem::Stabilization_Type::VISCOUS_GHOST ) ;
            aParameterLists( tSPIndex ).set( "function_parameters", "0.05" );
            aParameterLists( tSPIndex ).set( "leader_properties", "PropViscosity,Viscosity" );

            // ghost convective
            aParameterLists( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
            aParameterLists( tSPIndex ).set( "stabilization_name", "SPGPVelocity" );
            aParameterLists( tSPIndex ).set( "leader_phase_name", "PhaseFluid" );
            aParameterLists( tSPIndex ).set( "stabilization_type",  fem::Stabilization_Type::CONVECTIVE_GHOST ) ;
            aParameterLists( tSPIndex ).set( "function_parameters", "0.05" );
            aParameterLists( tSPIndex ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
            aParameterLists( tSPIndex ).set( "leader_properties", "PropDensity,Density" );

            // ghost pressure
            aParameterLists( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
            aParameterLists( tSPIndex ).set( "stabilization_name", "SPGPPressure" );
            aParameterLists( tSPIndex ).set( "leader_phase_name", "PhaseFluid" );
            aParameterLists( tSPIndex ).set( "stabilization_type",  fem::Stabilization_Type::PRESSURE_GHOST ) ;
            aParameterLists( tSPIndex ).set( "function_parameters", "0.005/1.0" );
            aParameterLists( tSPIndex ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
            aParameterLists( tSPIndex ).set( "leader_properties",
                    "PropViscosity,Viscosity;"
                    "PropDensity,Density" );

            // ghost temperature
            aParameterLists( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
            aParameterLists( tSPIndex ).set( "stabilization_name", "SPGPTemp" );
            aParameterLists( tSPIndex ).set( "leader_phase_name", "PhaseFluid" );
            aParameterLists( tSPIndex ).set( "follower_phase_name", "PhaseFluid" );
            aParameterLists( tSPIndex ).set( "stabilization_type",  fem::Stabilization_Type::GHOST_DISPL ) ;
            aParameterLists( tSPIndex ).set( "function_parameters", "0.005" );
            aParameterLists( tSPIndex ).set( "leader_properties", "PropConductivity,Material" );
            aParameterLists( tSPIndex ).set( "follower_properties", "PropConductivity,Material" );
            }

        //------------------------------------------------------------------------------
        // fill the IWG part of the parameter list

        // bulk NS incompressible for velocity
        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGVelocityBulk" );
        aParameterLists( tIWGIndex ).set( "IWG_bulk_type",  fem::Element_Type::BULK ) ;
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIWGIndex ).set( "IWG_type",  fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_BULK ) ;
        aParameterLists( tIWGIndex ).set( "dof_residual", "VX,VY" );
        aParameterLists( tIWGIndex ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists( tIWGIndex ).set( "stabilization_parameters", "SPIncFlow,IncompressibleFlow" );

        // bulk NS incompressible for pressure
        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGPressureBulk" );
        aParameterLists( tIWGIndex ).set( "IWG_bulk_type",  fem::Element_Type::BULK ) ;
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIWGIndex ).set( "IWG_type",  fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_BULK ) ;
        aParameterLists( tIWGIndex ).set( "dof_residual", "P" );
        aParameterLists( tIWGIndex ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists( tIWGIndex ).set( "stabilization_parameters", "SPIncFlow,IncompressibleFlow" );

        // bulk diffusion
        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGDiffusionBulk" );
        aParameterLists( tIWGIndex ).set( "IWG_bulk_type",  fem::Element_Type::BULK ) ;
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIWGIndex ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_BULK ) ;
        aParameterLists( tIWGIndex ).set( "dof_residual", "TEMP" );
        aParameterLists( tIWGIndex ).set( "leader_constitutive_models", "CMDiffusion,Diffusion" );

        // bulk advection
        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGAdvectionBulk" );
        aParameterLists( tIWGIndex ).set( "IWG_bulk_type",  fem::Element_Type::BULK ) ;
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIWGIndex ).set( "IWG_type",  fem::IWG_Type::ADVECTION_BULK ) ;
        aParameterLists( tIWGIndex ).set( "dof_residual", "TEMP" );
        aParameterLists( tIWGIndex ).set( "leader_constitutive_models", "CMDiffusion,Diffusion" );
        aParameterLists( tIWGIndex ).set( "stabilization_parameters", "SPSUPGTemp,SUPG" );

        // inlet velocity for velocity
        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGInletVelocity" );
        aParameterLists( tIWGIndex ).set( "IWG_bulk_type",  fem::Element_Type::SIDESET ) ;
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIWGIndex ).set( "IWG_type",  fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        aParameterLists( tIWGIndex ).set( "dof_residual", "VX,VY" );
        aParameterLists( tIWGIndex ).set( "leader_properties", "PropDirichletInletU,Dirichlet" );
        aParameterLists( tIWGIndex ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists( tIWGIndex ).set( "stabilization_parameters", "SPDirichletNitscheU,DirichletNitsche" );
        aParameterLists( tIWGIndex ).set( "neighbor_phases", "PhaseInlet" );

        // inlet velocity for pressure
        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGInletPressure" );
        aParameterLists( tIWGIndex ).set( "IWG_bulk_type",  fem::Element_Type::SIDESET ) ;
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIWGIndex ).set( "IWG_type",  fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_UNSYMMETRIC_NITSCHE ) ;
        aParameterLists( tIWGIndex ).set( "dof_residual", "P" );
        aParameterLists( tIWGIndex ).set( "leader_properties", "PropDirichletInletU,Dirichlet" );
        aParameterLists( tIWGIndex ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists( tIWGIndex ).set( "neighbor_phases", "PhaseInlet" );

        // inlet temperature
        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGInletTemp" );
        aParameterLists( tIWGIndex ).set( "IWG_bulk_type",  fem::Element_Type::SIDESET ) ;
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIWGIndex ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        aParameterLists( tIWGIndex ).set( "dof_residual", "TEMP" );
        aParameterLists( tIWGIndex ).set( "leader_properties", "PropInletTemp,Dirichlet" );
        aParameterLists( tIWGIndex ).set( "leader_constitutive_models", "CMDiffusion,Diffusion" );
        aParameterLists( tIWGIndex ).set( "stabilization_parameters", "SPDirichletNitscheT,DirichletNitsche" );
        aParameterLists( tIWGIndex ).set( "neighbor_phases", "PhaseInlet" );

        // side flux on cylinder
        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGCylinderFluxTemp" );
        aParameterLists( tIWGIndex ).set( "IWG_bulk_type",  fem::Element_Type::SIDESET ) ;
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIWGIndex ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_NEUMANN ) ;
        aParameterLists( tIWGIndex ).set( "dof_residual", "TEMP" );
        aParameterLists( tIWGIndex ).set( "leader_properties", "PropSideFlux,Neumann" );
        aParameterLists( tIWGIndex ).set( "neighbor_phases", "PhaseCylinder" );

        // zero velocity on wall and cylinders for velocity
        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGZeroVelocity" );
        aParameterLists( tIWGIndex ).set( "IWG_bulk_type",  fem::Element_Type::SIDESET ) ;
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIWGIndex ).set( "IWG_type",  fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        aParameterLists( tIWGIndex ).set( "dof_residual", "VX,VY" );
        aParameterLists( tIWGIndex ).set( "leader_properties", "PropDirichletZeroU,Dirichlet" );
        aParameterLists( tIWGIndex ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists( tIWGIndex ).set( "stabilization_parameters", "SPDirichletNitscheU,DirichletNitsche" );
        aParameterLists( tIWGIndex ).set( "neighbor_phases", "PhaseCylinder,PhaseWall" );

        // zero velocity on wall and cylinders for pressure
        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGZeroPressure" );
        aParameterLists( tIWGIndex ).set( "IWG_bulk_type",  fem::Element_Type::SIDESET ) ;
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIWGIndex ).set( "IWG_type",  fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_UNSYMMETRIC_NITSCHE ) ;
        aParameterLists( tIWGIndex ).set( "dof_residual", "P" );
        aParameterLists( tIWGIndex ).set( "leader_properties", "PropDirichletZeroU,Dirichlet" );
        aParameterLists( tIWGIndex ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists( tIWGIndex ).set( "neighbor_phases", "PhaseCylinder,PhaseWall" );

        if ( isGhost )
        {
            // ghost viscous
            aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists( tIWGIndex ).set( "IWG_name", "IWGGPViscous" );
            aParameterLists( tIWGIndex ).set( "IWG_bulk_type",  fem::Element_Type::DOUBLE_SIDESET ) ;
            aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
            aParameterLists( tIWGIndex ).set( "follower_phase_name", "PhaseFluid" );
            aParameterLists( tIWGIndex ).set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
            aParameterLists( tIWGIndex ).set( "dof_residual", "VX,VY" );
            aParameterLists( tIWGIndex ).set( "stabilization_parameters", "SPGPViscosity,GhostSP" );

            // ghost convective
            aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists( tIWGIndex ).set( "IWG_name", "IWGGPConvective" );
            aParameterLists( tIWGIndex ).set( "IWG_bulk_type",  fem::Element_Type::DOUBLE_SIDESET ) ;
            aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
            aParameterLists( tIWGIndex ).set( "follower_phase_name", "PhaseFluid" );
            aParameterLists( tIWGIndex ).set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
            aParameterLists( tIWGIndex ).set( "dof_residual", "VX,VY" );
            aParameterLists( tIWGIndex ).set( "stabilization_parameters", "SPGPVelocity,GhostSP" );

            // ghost pressure
            aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists( tIWGIndex ).set( "IWG_name", "IWGGPPressure" );
            aParameterLists( tIWGIndex ).set( "IWG_bulk_type",  fem::Element_Type::DOUBLE_SIDESET ) ;
            aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
            aParameterLists( tIWGIndex ).set( "follower_phase_name", "PhaseFluid" );
            aParameterLists( tIWGIndex ).set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
            aParameterLists( tIWGIndex ).set( "dof_residual", "P" );
            aParameterLists( tIWGIndex ).set( "stabilization_parameters", "SPGPPressure,GhostSP" );

            // ghost temperature
            aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists( tIWGIndex ).set( "IWG_name", "IWGGPTemp" );
            aParameterLists( tIWGIndex ).set( "IWG_bulk_type",  fem::Element_Type::DOUBLE_SIDESET ) ;
            aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
            aParameterLists( tIWGIndex ).set( "follower_phase_name", "PhaseFluid" );
            aParameterLists( tIWGIndex ).set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
            aParameterLists( tIWGIndex ).set( "dof_residual", "TEMP" );
            aParameterLists( tIWGIndex ).set( "stabilization_parameters", "SPGPTemp,GhostSP" );
            }

        //------------------------------------------------------------------------------
        // fill the IQI part of the parameter list

        // dof VX
        aParameterLists( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIBulkVX" );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIQIIndex ).set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists( tIQIIndex ).set( "dof_quantity", "VX,VY" );
        aParameterLists( tIQIIndex ).set( "vectorial_field_index", 0 );

        // dof VY
        aParameterLists( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIBulkVY" );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIQIIndex ).set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists( tIQIIndex ).set( "dof_quantity", "VX,VY" );
        aParameterLists( tIQIIndex ).set( "vectorial_field_index", 1 );

        // dof P
        aParameterLists( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIBulkP" );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIQIIndex ).set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists( tIQIIndex ).set( "dof_quantity", "P" );
        aParameterLists( tIQIIndex ).set( "vectorial_field_index", 0 );

        // dof TEMP
        aParameterLists( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIBulkTEMP" );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIQIIndex ).set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists( tIQIIndex ).set( "dof_quantity", "TEMP" );
        aParameterLists( tIQIIndex ).set( "vectorial_field_index", 0 );

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        aParameterLists( tFEMIndex ).add_parameter_list( prm::create_computation_parameter_list() );
        aParameterLists( tFEMIndex ).set( "print_physics_model", true );
    }

    void
    SOLParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL ) );

        aParameterLists( 1 ).add_parameter_list( moris::prm::create_linear_solver_parameter_list() );

        aParameterLists( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
        aParameterLists( 2 ).set( "NLA_combined_res_jac_assembly", true );
        aParameterLists( 2 ).set( "NLA_rel_res_norm_drop", 1e-04 );
        aParameterLists( 2 ).set( "NLA_relaxation_parameter", 0.80 );
        aParameterLists( 2 ).set( "NLA_max_iter", 20 );

        aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
        aParameterLists( 3 ).set( "NLA_DofTypes", "VX,VY;P;TEMP" );

        aParameterLists( 4 ).add_parameter_list( moris::prm::create_time_solver_algorithm_parameter_list() );

        aParameterLists( 5 ).add_parameter_list( moris::prm::create_time_solver_parameter_list() );
        aParameterLists( 5 ).set( "TSA_DofTypes", "VX,VY,P,TEMP" );
        aParameterLists( 5 ).set( "TSA_Initialize_Sol_Vec", "VX,1.0;VY,0.0;P,0.0;TEMP,0.0" );
        aParameterLists( 5 ).set( "TSA_Output_Indices", "0" );
        aParameterLists( 5 ).set( "TSA_Output_Criteria", "Output_Criterion" );

        aParameterLists( 6 ).add_parameter_list( moris::prm::create_solver_warehouse_parameterlist() );

        aParameterLists( 7 ).add_parameter_list( moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::NONE ) );
    }

    void
    MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_msi_parameter_list() );
    }

    void
    VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_vis_parameter_list() );
        aParameterLists( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", tExoFile ) );
        aParameterLists( 0 ).set( "Mesh_Type",  vis::VIS_Mesh_Type::STANDARD ) ;
        aParameterLists( 0 ).set( "Set_Names", "HMR_dummy_n_p0,HMR_dummy_c_p0" );
        aParameterLists( 0 ).set( "Field_Names",
                "VX,VY,P,TEMP,"
                "IQIBulkVX,IQIBulkVY,IQIBulkP,IQIBulkTEMP" );
        aParameterLists( 0 ).set( "Field_Type",
                "NODAL,NODAL,NODAL,NODAL,"
                "GLOBAL,GLOBAL,GLOBAL,GLOBAL" );
        aParameterLists( 0 ).set( "IQI_Names",
                "IQIBulkVX,IQIBulkVY,IQIBulkP,IQIBulkTEMP,"
                "IQIBulkVX,IQIBulkVY,IQIBulkP,IQIBulkTEMP" );
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
