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
    OPTParameterList( Vector< Submodule_Parameter_Lists >& tParameterlist )
    {
        tParameterlist.resize( 1 );

        tParameterlist( 0 ).add_parameter_list( prm::create_opt_problem_parameter_list() );
        tParameterlist( 0 ).set( "is_optimization_problem", false );
    }

    void
    HMRParameterList( Vector< Submodule_Parameter_Lists >& tParameterlist )
    {
        tParameterlist.resize( 1 );

        tParameterlist( 0 ).add_parameter_list( prm::create_hmr_parameter_list() );

        tParameterlist( 0 ).set( "number_of_elements_per_dimension", "176,88" );
        tParameterlist( 0 ).set( "domain_dimensions", "4,2" );
        tParameterlist( 0 ).set( "domain_offset", "-1.24,-0.86" );
        tParameterlist( 0 ).set( "domain_sidesets", "1,2,3,4" );
        tParameterlist( 0 ).set( "lagrange_output_meshes", "0" );

        tParameterlist( 0 ).set( "lagrange_orders", std::to_string( gInterpolationOrder ) );
        tParameterlist( 0 ).set( "lagrange_pattern", "0" );
        tParameterlist( 0 ).set( "bspline_orders", std::to_string( gInterpolationOrder ) );
        tParameterlist( 0 ).set( "bspline_pattern", "0" );

        tParameterlist( 0 ).set( "lagrange_to_bspline", "0" );

        tParameterlist( 0 ).set( "truncate_bsplines", 1 );
        tParameterlist( 0 ).set( "refinement_buffer", 1 );
        tParameterlist( 0 ).set( "staircase_buffer", 1 );
        tParameterlist( 0 ).set( "initial_refinement", "1" );
        tParameterlist( 0 ).set( "initial_refinement_pattern", "0" );

        tParameterlist( 0 ).set( "use_number_aura", 1 );

        tParameterlist( 0 ).set( "use_multigrid", 0 );
        tParameterlist( 0 ).set( "severity_level", 0 );
    }

    void
    XTKParameterList( Vector< Submodule_Parameter_Lists >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).add_parameter_list( prm::create_xtk_parameter_list() );
        tParameterlist( 0 ).set( "decompose", true );
        tParameterlist( 0 ).set( "decomposition_type", "conformal" );
        tParameterlist( 0 ).set( "enrich", true );
        tParameterlist( 0 ).set( "basis_rank", "bspline" );
        tParameterlist( 0 ).set( "enrich_mesh_indices", "0" );
        tParameterlist( 0 ).set( "ghost_stab", isGhost );
        tParameterlist( 0 ).set( "multigrid", false );
        tParameterlist( 0 ).set( "verbose", true );
        tParameterlist( 0 ).set( "print_enriched_ig_mesh", false );
        tParameterlist( 0 ).set( "exodus_output_XTK_ig_mesh", true );
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
    GENParameterList( Vector< Submodule_Parameter_Lists >& tParameterlist )
    {
        tParameterlist.resize( 3 );
        tParameterlist( 0 ).add_parameter_list( prm::create_gen_parameter_list() );
        tParameterlist( 0 ).set( "output_mesh_file", tGENOutputFile );

        tParameterlist( 0 ).set( "number_of_phases", 5 );
        tParameterlist( 0 ).set( "phase_function_name", "get_phase_index" );

        // Bottom plane
        tParameterlist( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
        tParameterlist( 1 ).set( "center_x", 0.0 );
        tParameterlist( 1 ).set( "center_y", 0.0 );
        tParameterlist( 1 ).set( "normal_x", 0.0 );
        tParameterlist( 1 ).set( "normal_y", 1.0 );
        tParameterlist( 1 ).set( "number_of_refinements", 1 );
        tParameterlist( 1 ).set( "refinement_mesh_index", 0 );

        // Top plane
        tParameterlist( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
        tParameterlist( 1 ).set( "center_x", 0.0 );
        tParameterlist( 1 ).set( "center_y", 0.41 );
        tParameterlist( 1 ).set( "normal_x", 0.0 );
        tParameterlist( 1 ).set( "normal_y", 1.0 );
        tParameterlist( 1 ).set( "number_of_refinements", 1 );
        tParameterlist( 1 ).set( "refinement_mesh_index", 0 );

        // Left plane
        tParameterlist( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
        tParameterlist( 1 ).set( "center_x", 0.0 );
        tParameterlist( 1 ).set( "center_y", 0.0 );
        tParameterlist( 1 ).set( "normal_x", 1.0 );
        tParameterlist( 1 ).set( "normal_y", 0.0 );
        tParameterlist( 1 ).set( "number_of_refinements", 1 );
        tParameterlist( 1 ).set( "refinement_mesh_index", 0 );

        // Right plane
        tParameterlist( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
        tParameterlist( 1 ).set( "center_x", 2.2 );
        tParameterlist( 1 ).set( "center_y", 0.0 );
        tParameterlist( 1 ).set( "normal_x", 1.0 );
        tParameterlist( 1 ).set( "normal_y", 0.0 );
        tParameterlist( 1 ).set( "number_of_refinements", 1 );
        tParameterlist( 1 ).set( "refinement_mesh_index", 0 );

        // Cylinder 1
        tParameterlist( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::CIRCLE ) );
        tParameterlist( 1 ).set( "center_x", 0.4 );
        tParameterlist( 1 ).set( "center_y", 0.1 );
        tParameterlist( 1 ).set( "radius", 0.05 );
        tParameterlist( 1 ).set( "number_of_refinements", 1 );
        tParameterlist( 1 ).set( "refinement_mesh_index", 0 );

        // Cylinder 2
        tParameterlist( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::CIRCLE ) );
        tParameterlist( 1 ).set( "center_x", 0.8 );
        tParameterlist( 1 ).set( "center_y", 0.3 );
        tParameterlist( 1 ).set( "radius", 0.05 );
        tParameterlist( 1 ).set( "number_of_refinements", 1 );
        tParameterlist( 1 ).set( "refinement_mesh_index", 0 );

        // Cylinder 3
        tParameterlist( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::CIRCLE ) );
        tParameterlist( 1 ).set( "center_x", 1.2 );
        tParameterlist( 1 ).set( "center_y", 0.1 );
        tParameterlist( 1 ).set( "radius", 0.05 );
        tParameterlist( 1 ).set( "number_of_refinements", 1 );
        tParameterlist( 1 ).set( "refinement_mesh_index", 0 );

        // Cylinder 4
        tParameterlist( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::CIRCLE ) );
        tParameterlist( 1 ).set( "center_x", 1.6 );
        tParameterlist( 1 ).set( "center_y", 0.3 );
        tParameterlist( 1 ).set( "radius", 0.05 );
        tParameterlist( 1 ).set( "number_of_refinements", 1 );
        tParameterlist( 1 ).set( "refinement_mesh_index", 0 );
    }

    void
    FEMParameterList( Vector< Submodule_Parameter_Lists >& tParameterList )
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

        tParameterList( tPhaseIndex ).add_parameter_list( prm::create_phase_parameter_list() );
        tParameterList( tPhaseIndex ).set( "phase_name", "PhaseFluid" );
        tParameterList( tPhaseIndex ).set( "phase_indices", "0" );

        tParameterList( tPhaseIndex ).add_parameter_list( prm::create_phase_parameter_list() );
        tParameterList( tPhaseIndex ).set( "phase_name", "PhaseInlet" );
        tParameterList( tPhaseIndex ).set( "phase_indices", "1" );

        tParameterList( tPhaseIndex ).add_parameter_list( prm::create_phase_parameter_list() );
        tParameterList( tPhaseIndex ).set( "phase_name", "PhaseCylinder" );
        tParameterList( tPhaseIndex ).set( "phase_indices", "2" );

        tParameterList( tPhaseIndex ).add_parameter_list( prm::create_phase_parameter_list() );
        tParameterList( tPhaseIndex ).set( "phase_name", "PhaseWall" );
        tParameterList( tPhaseIndex ).set( "phase_indices", "3" );

        //------------------------------------------------------------------------------
        // fill the property part of the parameter list


        // fluid properties ------------------------------------------------------------
        // viscosity
        tParameterList( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( tPropIndex ).set( "property_name", "PropViscosity" );
        tParameterList( tPropIndex ).set( "function_parameters", "0.001" );

        // density
        tParameterList( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( tPropIndex ).set( "property_name", "PropDensity" );
        tParameterList( tPropIndex ).set( "function_parameters", "1.0" );

        // capacity
        tParameterList( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( tPropIndex ).set( "property_name", "PropCapacity" );
        tParameterList( tPropIndex ).set( "function_parameters", "1.0" );

        // conductivity
        tParameterList( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( tPropIndex ).set( "property_name", "PropConductivity" );
        tParameterList( tPropIndex ).set( "function_parameters", "0.00005" );

        // BC properties ---------------------------------------------------------------
        // inlet velocity
        tParameterList( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( tPropIndex ).set( "property_name", "PropDirichletInletU" );
        tParameterList( tPropIndex ).set( "function_parameters", "1.5" );
        tParameterList( tPropIndex ).set( "value_function", "Func_Inlet_U" );

        // zero velocity
        tParameterList( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( tPropIndex ).set( "property_name", "PropDirichletZeroU" );
        tParameterList( tPropIndex ).set( "function_parameters", "0.0;0.0" );

        // inlet temperature
        tParameterList( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( tPropIndex ).set( "property_name", "PropInletTemp" );
        tParameterList( tPropIndex ).set( "function_parameters", "0.0" );

        // side flux on cylinders
        tParameterList( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( tPropIndex ).set( "property_name", "PropSideFlux" );
        tParameterList( tPropIndex ).set( "function_parameters", "2.0" );

        //------------------------------------------------------------------------------
        // fill the constitutive model part of the parameter list


        // Fluid CM --------------------------------------------------------------------
        // NS incompressible
        tParameterList( tCMIndex ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        tParameterList( tCMIndex ).set( "constitutive_name", "CMFluid" );
        tParameterList( tCMIndex ).set( "phase_name", "PhaseFluid" );
        tParameterList( tCMIndex ).set( "constitutive_type",  fem::Constitutive_Type::FLUID_INCOMPRESSIBLE ) ;
        tParameterList( tCMIndex ).set( "dof_dependencies", std::pair< std::string, std::string >( "VX,VY;P", "Velocity,Pressure" ) );
        tParameterList( tCMIndex ).set( "properties",
                "PropViscosity,Viscosity;"
                "PropDensity  ,Density" );

        // diffusion
        tParameterList( tCMIndex ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        tParameterList( tCMIndex ).set( "constitutive_name", "CMDiffusion" );
        tParameterList( tCMIndex ).set( "phase_name", "PhaseFluid" );
        tParameterList( tCMIndex ).set( "constitutive_type",  fem::Constitutive_Type::DIFF_LIN_ISO ) ;
        tParameterList( tCMIndex ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( tCMIndex ).set( "properties",
                "PropConductivity,Conductivity;"
                "PropDensity     ,Density;"
                "PropCapacity    ,HeatCapacity" );

        //------------------------------------------------------------------------------
        // fill the stabilization parameter part of the parameter list

        // SUPG/PSPG fluid incompressible (note: use time step in tau although static; for testing purposes only)
        tParameterList( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex ).set( "stabilization_name", "SPIncFlow" );
        tParameterList( tSPIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tSPIndex ).set( "stabilization_type",  fem::Stabilization_Type::INCOMPRESSIBLE_FLOW ) ;
        tParameterList( tSPIndex ).set( "function_parameters", "36.0/1.0" );
        tParameterList( tSPIndex ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY;P", "Velocity,Pressure" ) );
        tParameterList( tSPIndex ).set( "leader_properties",
                "PropViscosity,Viscosity;"
                "PropDensity,Density" );

        // SUPG advection (note: use time step in tau although static; for testing purposes only)
        tParameterList( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex ).set( "stabilization_name", "SPSUPGTemp" );
        tParameterList( tSPIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tSPIndex ).set( "stabilization_type",  fem::Stabilization_Type::SUPG_ADVECTION ) ;
        tParameterList( tSPIndex ).set( "function_parameters", "1.0" );
        tParameterList( tSPIndex ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
        tParameterList( tSPIndex ).set( "leader_properties", "PropConductivity,Conductivity" );

        // Dirichlet Nitsche for velocity
        tParameterList( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex ).set( "stabilization_name", "SPDirichletNitscheU" );
        tParameterList( tSPIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tSPIndex ).set( "stabilization_type",  fem::Stabilization_Type::VELOCITY_DIRICHLET_NITSCHE ) ;
        tParameterList( tSPIndex ).set( "function_parameters", "100.0/1.0" );
        tParameterList( tSPIndex ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
        tParameterList( tSPIndex ).set( "leader_properties",
                "PropViscosity,Viscosity;"
                "PropDensity,Density" );

        // Dirichhlet Nitsche for temperature
        tParameterList( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex ).set( "stabilization_name", "SPDirichletNitscheT" );
        tParameterList( tSPIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tSPIndex ).set( "stabilization_type",  fem::Stabilization_Type::DIRICHLET_NITSCHE ) ;
        tParameterList( tSPIndex ).set( "function_parameters", "100.0" );
        tParameterList( tSPIndex ).set( "leader_properties", "PropConductivity,Material" );

        if ( isGhost )
        {
            // ghost viscous
            tParameterList( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
            tParameterList( tSPIndex ).set( "stabilization_name", "SPGPViscosity" );
            tParameterList( tSPIndex ).set( "leader_phase_name", "PhaseFluid" );
            tParameterList( tSPIndex ).set( "stabilization_type",  fem::Stabilization_Type::VISCOUS_GHOST ) ;
            tParameterList( tSPIndex ).set( "function_parameters", "0.05" );
            tParameterList( tSPIndex ).set( "leader_properties", "PropViscosity,Viscosity" );

            // ghost convective
            tParameterList( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
            tParameterList( tSPIndex ).set( "stabilization_name", "SPGPVelocity" );
            tParameterList( tSPIndex ).set( "leader_phase_name", "PhaseFluid" );
            tParameterList( tSPIndex ).set( "stabilization_type",  fem::Stabilization_Type::CONVECTIVE_GHOST ) ;
            tParameterList( tSPIndex ).set( "function_parameters", "0.05" );
            tParameterList( tSPIndex ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
            tParameterList( tSPIndex ).set( "leader_properties", "PropDensity,Density" );

            // ghost pressure
            tParameterList( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
            tParameterList( tSPIndex ).set( "stabilization_name", "SPGPPressure" );
            tParameterList( tSPIndex ).set( "leader_phase_name", "PhaseFluid" );
            tParameterList( tSPIndex ).set( "stabilization_type",  fem::Stabilization_Type::PRESSURE_GHOST ) ;
            tParameterList( tSPIndex ).set( "function_parameters", "0.005/1.0" );
            tParameterList( tSPIndex ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
            tParameterList( tSPIndex ).set( "leader_properties",
                    "PropViscosity,Viscosity;"
                    "PropDensity,Density" );

            // ghost temperature
            tParameterList( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
            tParameterList( tSPIndex ).set( "stabilization_name", "SPGPTemp" );
            tParameterList( tSPIndex ).set( "leader_phase_name", "PhaseFluid" );
            tParameterList( tSPIndex ).set( "follower_phase_name", "PhaseFluid" );
            tParameterList( tSPIndex ).set( "stabilization_type",  fem::Stabilization_Type::GHOST_DISPL ) ;
            tParameterList( tSPIndex ).set( "function_parameters", "0.005" );
            tParameterList( tSPIndex ).set( "leader_properties", "PropConductivity,Material" );
            tParameterList( tSPIndex ).set( "follower_properties", "PropConductivity,Material" );
            }

        //------------------------------------------------------------------------------
        // fill the IWG part of the parameter list

        // bulk NS incompressible for velocity
        tParameterList( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex ).set( "IWG_name", "IWGVelocityBulk" );
        tParameterList( tIWGIndex ).set( "IWG_bulk_type",  fem::Element_Type::BULK ) ;
        tParameterList( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex ).set( "IWG_type",  fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_BULK ) ;
        tParameterList( tIWGIndex ).set( "dof_residual", "VX,VY" );
        tParameterList( tIWGIndex ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        tParameterList( tIWGIndex ).set( "stabilization_parameters", "SPIncFlow,IncompressibleFlow" );

        // bulk NS incompressible for pressure
        tParameterList( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex ).set( "IWG_name", "IWGPressureBulk" );
        tParameterList( tIWGIndex ).set( "IWG_bulk_type",  fem::Element_Type::BULK ) ;
        tParameterList( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex ).set( "IWG_type",  fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_BULK ) ;
        tParameterList( tIWGIndex ).set( "dof_residual", "P" );
        tParameterList( tIWGIndex ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        tParameterList( tIWGIndex ).set( "stabilization_parameters", "SPIncFlow,IncompressibleFlow" );

        // bulk diffusion
        tParameterList( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex ).set( "IWG_name", "IWGDiffusionBulk" );
        tParameterList( tIWGIndex ).set( "IWG_bulk_type",  fem::Element_Type::BULK ) ;
        tParameterList( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_BULK ) ;
        tParameterList( tIWGIndex ).set( "dof_residual", "TEMP" );
        tParameterList( tIWGIndex ).set( "leader_constitutive_models", "CMDiffusion,Diffusion" );

        // bulk advection
        tParameterList( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex ).set( "IWG_name", "IWGAdvectionBulk" );
        tParameterList( tIWGIndex ).set( "IWG_bulk_type",  fem::Element_Type::BULK ) ;
        tParameterList( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex ).set( "IWG_type",  fem::IWG_Type::ADVECTION_BULK ) ;
        tParameterList( tIWGIndex ).set( "dof_residual", "TEMP" );
        tParameterList( tIWGIndex ).set( "leader_constitutive_models", "CMDiffusion,Diffusion" );
        tParameterList( tIWGIndex ).set( "stabilization_parameters", "SPSUPGTemp,SUPG" );

        // inlet velocity for velocity
        tParameterList( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex ).set( "IWG_name", "IWGInletVelocity" );
        tParameterList( tIWGIndex ).set( "IWG_bulk_type",  fem::Element_Type::SIDESET ) ;
        tParameterList( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex ).set( "IWG_type",  fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        tParameterList( tIWGIndex ).set( "dof_residual", "VX,VY" );
        tParameterList( tIWGIndex ).set( "leader_properties", "PropDirichletInletU,Dirichlet" );
        tParameterList( tIWGIndex ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        tParameterList( tIWGIndex ).set( "stabilization_parameters", "SPDirichletNitscheU,DirichletNitsche" );
        tParameterList( tIWGIndex ).set( "neighbor_phases", "PhaseInlet" );

        // inlet velocity for pressure
        tParameterList( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex ).set( "IWG_name", "IWGInletPressure" );
        tParameterList( tIWGIndex ).set( "IWG_bulk_type",  fem::Element_Type::SIDESET ) ;
        tParameterList( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex ).set( "IWG_type",  fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_UNSYMMETRIC_NITSCHE ) ;
        tParameterList( tIWGIndex ).set( "dof_residual", "P" );
        tParameterList( tIWGIndex ).set( "leader_properties", "PropDirichletInletU,Dirichlet" );
        tParameterList( tIWGIndex ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        tParameterList( tIWGIndex ).set( "neighbor_phases", "PhaseInlet" );

        // inlet temperature
        tParameterList( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex ).set( "IWG_name", "IWGInletTemp" );
        tParameterList( tIWGIndex ).set( "IWG_bulk_type",  fem::Element_Type::SIDESET ) ;
        tParameterList( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        tParameterList( tIWGIndex ).set( "dof_residual", "TEMP" );
        tParameterList( tIWGIndex ).set( "leader_properties", "PropInletTemp,Dirichlet" );
        tParameterList( tIWGIndex ).set( "leader_constitutive_models", "CMDiffusion,Diffusion" );
        tParameterList( tIWGIndex ).set( "stabilization_parameters", "SPDirichletNitscheT,DirichletNitsche" );
        tParameterList( tIWGIndex ).set( "neighbor_phases", "PhaseInlet" );

        // side flux on cylinder
        tParameterList( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex ).set( "IWG_name", "IWGCylinderFluxTemp" );
        tParameterList( tIWGIndex ).set( "IWG_bulk_type",  fem::Element_Type::SIDESET ) ;
        tParameterList( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_NEUMANN ) ;
        tParameterList( tIWGIndex ).set( "dof_residual", "TEMP" );
        tParameterList( tIWGIndex ).set( "leader_properties", "PropSideFlux,Neumann" );
        tParameterList( tIWGIndex ).set( "neighbor_phases", "PhaseCylinder" );

        // zero velocity on wall and cylinders for velocity
        tParameterList( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex ).set( "IWG_name", "IWGZeroVelocity" );
        tParameterList( tIWGIndex ).set( "IWG_bulk_type",  fem::Element_Type::SIDESET ) ;
        tParameterList( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex ).set( "IWG_type",  fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        tParameterList( tIWGIndex ).set( "dof_residual", "VX,VY" );
        tParameterList( tIWGIndex ).set( "leader_properties", "PropDirichletZeroU,Dirichlet" );
        tParameterList( tIWGIndex ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        tParameterList( tIWGIndex ).set( "stabilization_parameters", "SPDirichletNitscheU,DirichletNitsche" );
        tParameterList( tIWGIndex ).set( "neighbor_phases", "PhaseCylinder,PhaseWall" );

        // zero velocity on wall and cylinders for pressure
        tParameterList( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex ).set( "IWG_name", "IWGZeroPressure" );
        tParameterList( tIWGIndex ).set( "IWG_bulk_type",  fem::Element_Type::SIDESET ) ;
        tParameterList( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex ).set( "IWG_type",  fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_UNSYMMETRIC_NITSCHE ) ;
        tParameterList( tIWGIndex ).set( "dof_residual", "P" );
        tParameterList( tIWGIndex ).set( "leader_properties", "PropDirichletZeroU,Dirichlet" );
        tParameterList( tIWGIndex ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        tParameterList( tIWGIndex ).set( "neighbor_phases", "PhaseCylinder,PhaseWall" );

        if ( isGhost )
        {
            // ghost viscous
            tParameterList( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
            tParameterList( tIWGIndex ).set( "IWG_name", "IWGGPViscous" );
            tParameterList( tIWGIndex ).set( "IWG_bulk_type",  fem::Element_Type::DOUBLE_SIDESET ) ;
            tParameterList( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
            tParameterList( tIWGIndex ).set( "follower_phase_name", "PhaseFluid" );
            tParameterList( tIWGIndex ).set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
            tParameterList( tIWGIndex ).set( "dof_residual", "VX,VY" );
            tParameterList( tIWGIndex ).set( "stabilization_parameters", "SPGPViscosity,GhostSP" );

            // ghost convective
            tParameterList( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
            tParameterList( tIWGIndex ).set( "IWG_name", "IWGGPConvective" );
            tParameterList( tIWGIndex ).set( "IWG_bulk_type",  fem::Element_Type::DOUBLE_SIDESET ) ;
            tParameterList( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
            tParameterList( tIWGIndex ).set( "follower_phase_name", "PhaseFluid" );
            tParameterList( tIWGIndex ).set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
            tParameterList( tIWGIndex ).set( "dof_residual", "VX,VY" );
            tParameterList( tIWGIndex ).set( "stabilization_parameters", "SPGPVelocity,GhostSP" );

            // ghost pressure
            tParameterList( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
            tParameterList( tIWGIndex ).set( "IWG_name", "IWGGPPressure" );
            tParameterList( tIWGIndex ).set( "IWG_bulk_type",  fem::Element_Type::DOUBLE_SIDESET ) ;
            tParameterList( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
            tParameterList( tIWGIndex ).set( "follower_phase_name", "PhaseFluid" );
            tParameterList( tIWGIndex ).set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
            tParameterList( tIWGIndex ).set( "dof_residual", "P" );
            tParameterList( tIWGIndex ).set( "stabilization_parameters", "SPGPPressure,GhostSP" );

            // ghost temperature
            tParameterList( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
            tParameterList( tIWGIndex ).set( "IWG_name", "IWGGPTemp" );
            tParameterList( tIWGIndex ).set( "IWG_bulk_type",  fem::Element_Type::DOUBLE_SIDESET ) ;
            tParameterList( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
            tParameterList( tIWGIndex ).set( "follower_phase_name", "PhaseFluid" );
            tParameterList( tIWGIndex ).set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
            tParameterList( tIWGIndex ).set( "dof_residual", "TEMP" );
            tParameterList( tIWGIndex ).set( "stabilization_parameters", "SPGPTemp,GhostSP" );
            }

        //------------------------------------------------------------------------------
        // fill the IQI part of the parameter list

        // dof VX
        tParameterList( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex ).set( "IQI_name", "IQIBulkVX" );
        tParameterList( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIQIIndex ).set( "IQI_type",  fem::IQI_Type::DOF ) ;
        tParameterList( tIQIIndex ).set( "dof_quantity", "VX,VY" );
        tParameterList( tIQIIndex ).set( "vectorial_field_index", 0 );

        // dof VY
        tParameterList( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex ).set( "IQI_name", "IQIBulkVY" );
        tParameterList( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIQIIndex ).set( "IQI_type",  fem::IQI_Type::DOF ) ;
        tParameterList( tIQIIndex ).set( "dof_quantity", "VX,VY" );
        tParameterList( tIQIIndex ).set( "vectorial_field_index", 1 );

        // dof P
        tParameterList( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex ).set( "IQI_name", "IQIBulkP" );
        tParameterList( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIQIIndex ).set( "IQI_type",  fem::IQI_Type::DOF ) ;
        tParameterList( tIQIIndex ).set( "dof_quantity", "P" );
        tParameterList( tIQIIndex ).set( "vectorial_field_index", 0 );

        // dof TEMP
        tParameterList( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex ).set( "IQI_name", "IQIBulkTEMP" );
        tParameterList( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIQIIndex ).set( "IQI_type",  fem::IQI_Type::DOF ) ;
        tParameterList( tIQIIndex ).set( "dof_quantity", "TEMP" );
        tParameterList( tIQIIndex ).set( "vectorial_field_index", 0 );

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        tParameterList( tFEMIndex ).add_parameter_list( prm::create_computation_parameter_list() );
        tParameterList( tFEMIndex ).set( "print_physics_model", true );
    }

    void
    SOLParameterList( Vector< Submodule_Parameter_Lists >& tParameterlist )
    {
        tParameterlist.resize( 8 );

        tParameterlist( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL ) );

        tParameterlist( 1 ).add_parameter_list( moris::prm::create_linear_solver_parameter_list() );

        tParameterlist( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
        tParameterlist( 2 ).set( "NLA_combined_res_jac_assembly", true );
        tParameterlist( 2 ).set( "NLA_rel_res_norm_drop", 1e-04 );
        tParameterlist( 2 ).set( "NLA_relaxation_parameter", 0.80 );
        tParameterlist( 2 ).set( "NLA_max_iter", 20 );

        tParameterlist( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
        tParameterlist( 3 ).set( "NLA_DofTypes", "VX,VY;P;TEMP" );

        tParameterlist( 4 ).add_parameter_list( moris::prm::create_time_solver_algorithm_parameter_list() );

        tParameterlist( 5 ).add_parameter_list( moris::prm::create_time_solver_parameter_list() );
        tParameterlist( 5 ).set( "TSA_DofTypes", "VX,VY,P,TEMP" );
        tParameterlist( 5 ).set( "TSA_Initialize_Sol_Vec", "VX,1.0;VY,0.0;P,0.0;TEMP,0.0" );
        tParameterlist( 5 ).set( "TSA_Output_Indices", "0" );
        tParameterlist( 5 ).set( "TSA_Output_Criteria", "Output_Criterion" );

        tParameterlist( 6 ).add_parameter_list( moris::prm::create_solver_warehouse_parameterlist() );

        tParameterlist( 7 ).add_parameter_list( moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::NONE ) );
    }

    void
    MSIParameterList( Vector< Submodule_Parameter_Lists >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).add_parameter_list( prm::create_msi_parameter_list() );
    }

    void
    VISParameterList( Vector< Submodule_Parameter_Lists >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).add_parameter_list( prm::create_vis_parameter_list() );
        tParameterlist( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", tExoFile ) );
        tParameterlist( 0 ).set( "Mesh_Type",  vis::VIS_Mesh_Type::STANDARD ) ;
        tParameterlist( 0 ).set( "Set_Names", "HMR_dummy_n_p0,HMR_dummy_c_p0" );
        tParameterlist( 0 ).set( "Field_Names",
                "VX,VY,P,TEMP,"
                "IQIBulkVX,IQIBulkVY,IQIBulkP,IQIBulkTEMP" );
        tParameterlist( 0 ).set( "Field_Type",
                "NODAL,NODAL,NODAL,NODAL,"
                "GLOBAL,GLOBAL,GLOBAL,GLOBAL" );
        tParameterlist( 0 ).set( "IQI_Names",
                "IQIBulkVX,IQIBulkVY,IQIBulkP,IQIBulkTEMP,"
                "IQIBulkVX,IQIBulkVY,IQIBulkP,IQIBulkTEMP" );
        tParameterlist( 0 ).set( "Save_Frequency", 1 );
    }

    void
    MORISGENERALParameterList( Vector< Submodule_Parameter_Lists >& tParameterlist )
    {
    }

    //------------------------------------------------------------------------------
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
