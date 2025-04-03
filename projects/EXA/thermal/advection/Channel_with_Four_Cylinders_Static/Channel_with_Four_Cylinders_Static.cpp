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
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
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

        aParameterLists.set( "is_optimization_problem", false );
    }

    void
    HMRParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists.set( "number_of_elements_per_dimension", 176, 88 );
        aParameterLists.set( "domain_dimensions", 4.0, 2.0 );
        aParameterLists.set( "domain_offset", -1.24, -0.86 );
        aParameterLists.set( "lagrange_output_meshes", "0" );

        aParameterLists.set( "lagrange_orders", std::to_string( gInterpolationOrder ) );
        aParameterLists.set( "lagrange_pattern", "0" );
        aParameterLists.set( "bspline_orders", std::to_string( gInterpolationOrder ) );
        aParameterLists.set( "bspline_pattern", "0" );

        aParameterLists.set( "lagrange_to_bspline", "0" );

        aParameterLists.set( "refinement_buffer", 1 );
        aParameterLists.set( "staircase_buffer", 1 );
        aParameterLists.set( "pattern_initial_refinement", 1 );
    }

    void
    XTKParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "decompose", true );
        aParameterLists.set( "decomposition_type", "conformal" );
        aParameterLists.set( "enrich_mesh_indices", "0" );
        aParameterLists.set( "ghost_stab", isGhost );
        aParameterLists.set( "multigrid", false );
        aParameterLists.set( "verbose", true );
        aParameterLists.set( "print_enriched_ig_mesh", false );
        aParameterLists.set( "exodus_output_XTK_ig_mesh", true );
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
        aParameterLists.set( "output_mesh_file", tGENOutputFile );

        aParameterLists.set( "number_of_phases", 5 );
        aParameterLists.set( "phase_function_name", "get_phase_index" );

        // Bottom plane
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( gen::Field_Type::LINE );
        aParameterLists.set( "center_x", 0.0 );
        aParameterLists.set( "center_y", 0.0 );
        aParameterLists.set( "normal_x", 0.0 );
        aParameterLists.set( "normal_y", 1.0 );
        aParameterLists.set( "number_of_refinements", 1 );
        aParameterLists.set( "refinement_mesh_index", 0 );

        // Top plane
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( gen::Field_Type::LINE );
        aParameterLists.set( "center_x", 0.0 );
        aParameterLists.set( "center_y", 0.41 );
        aParameterLists.set( "normal_x", 0.0 );
        aParameterLists.set( "normal_y", 1.0 );
        aParameterLists.set( "number_of_refinements", 1 );
        aParameterLists.set( "refinement_mesh_index", 0 );

        // Left plane
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( gen::Field_Type::LINE );
        aParameterLists.set( "center_x", 0.0 );
        aParameterLists.set( "center_y", 0.0 );
        aParameterLists.set( "normal_x", 1.0 );
        aParameterLists.set( "normal_y", 0.0 );
        aParameterLists.set( "number_of_refinements", 1 );
        aParameterLists.set( "refinement_mesh_index", 0 );

        // Right plane
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( gen::Field_Type::LINE );
        aParameterLists.set( "center_x", 2.2 );
        aParameterLists.set( "center_y", 0.0 );
        aParameterLists.set( "normal_x", 1.0 );
        aParameterLists.set( "normal_y", 0.0 );
        aParameterLists.set( "number_of_refinements", 1 );
        aParameterLists.set( "refinement_mesh_index", 0 );

        // Cylinder 1
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( gen::Field_Type::CIRCLE );
        aParameterLists.set( "center_x", 0.4 );
        aParameterLists.set( "center_y", 0.1 );
        aParameterLists.set( "radius", 0.05 );
        aParameterLists.set( "number_of_refinements", 1 );
        aParameterLists.set( "refinement_mesh_index", 0 );

        // Cylinder 2
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( gen::Field_Type::CIRCLE );
        aParameterLists.set( "center_x", 0.8 );
        aParameterLists.set( "center_y", 0.3 );
        aParameterLists.set( "radius", 0.05 );
        aParameterLists.set( "number_of_refinements", 1 );
        aParameterLists.set( "refinement_mesh_index", 0 );

        // Cylinder 3
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( gen::Field_Type::CIRCLE );
        aParameterLists.set( "center_x", 1.2 );
        aParameterLists.set( "center_y", 0.1 );
        aParameterLists.set( "radius", 0.05 );
        aParameterLists.set( "number_of_refinements", 1 );
        aParameterLists.set( "refinement_mesh_index", 0 );

        // Cylinder 4
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( gen::Field_Type::CIRCLE );
        aParameterLists.set( "center_x", 1.6 );
        aParameterLists.set( "center_y", 0.3 );
        aParameterLists.set( "radius", 0.05 );
        aParameterLists.set( "number_of_refinements", 1 );
        aParameterLists.set( "refinement_mesh_index", 0 );
    }

    void
    FEMParameterList( Module_Parameter_Lists& aParameterLists )
    {

        //------------------------------------------------------------------------------

        aParameterLists( FEM::PHASES ).add_parameter_list();
        aParameterLists.set( "phase_name", "PhaseFluid" );
        aParameterLists.set( "phase_indices", "0" );

        aParameterLists( FEM::PHASES ).add_parameter_list();
        aParameterLists.set( "phase_name", "PhaseInlet" );
        aParameterLists.set( "phase_indices", "1" );

        aParameterLists( FEM::PHASES ).add_parameter_list();
        aParameterLists.set( "phase_name", "PhaseCylinder" );
        aParameterLists.set( "phase_indices", "2" );

        aParameterLists( FEM::PHASES ).add_parameter_list();
        aParameterLists.set( "phase_name", "PhaseWall" );
        aParameterLists.set( "phase_indices", "3" );

        //------------------------------------------------------------------------------
        // fill the property part of the parameter list

        // fluid properties ------------------------------------------------------------
        // viscosity
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropViscosity" );
        aParameterLists.set( "function_parameters", "0.001" );

        // density
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDensity" );
        aParameterLists.set( "function_parameters", "1.0" );

        // capacity
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropCapacity" );
        aParameterLists.set( "function_parameters", "1.0" );

        // conductivity
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropConductivity" );
        aParameterLists.set( "function_parameters", "0.00005" );

        // BC properties ---------------------------------------------------------------
        // inlet velocity
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDirichletInletU" );
        aParameterLists.set( "function_parameters", "1.5" );
        aParameterLists.set( "value_function", "Func_Inlet_U" );

        // zero velocity
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDirichletZeroU" );
        aParameterLists.set( "function_parameters", "0.0;0.0" );

        // inlet temperature
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropInletTemp" );
        aParameterLists.set( "function_parameters", "0.0" );

        // side flux on cylinders
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropSideFlux" );
        aParameterLists.set( "function_parameters", "2.0" );

        //------------------------------------------------------------------------------
        // fill the constitutive model part of the parameter list

        // Fluid CM --------------------------------------------------------------------
        // NS incompressible
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMFluid" );
        aParameterLists.set( "phase_name", "PhaseFluid" );
        aParameterLists.set( "constitutive_type", fem::Constitutive_Type::FLUID_INCOMPRESSIBLE );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "VX,VY;P", "Velocity,Pressure" ) );
        aParameterLists.set( "properties",
                "PropViscosity,Viscosity;"
                "PropDensity  ,Density" );

        // diffusion
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMDiffusion" );
        aParameterLists.set( "phase_name", "PhaseFluid" );
        aParameterLists.set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists.set( "properties",
                "PropConductivity,Conductivity;"
                "PropDensity     ,Density;"
                "PropCapacity    ,HeatCapacity" );

        //------------------------------------------------------------------------------
        // fill the stabilization parameter part of the parameter list

        // SUPG/PSPG fluid incompressible (note: use time step in tau although static; for testing purposes only)
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPIncFlow" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::INCOMPRESSIBLE_FLOW );
        aParameterLists.set( "function_parameters", "36.0/1.0" );
        aParameterLists.set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY;P", "Velocity,Pressure" ) );
        aParameterLists.set( "leader_properties",
                "PropViscosity,Viscosity;"
                "PropDensity,Density" );

        // SUPG advection (note: use time step in tau although static; for testing purposes only)
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPSUPGTemp" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::SUPG_ADVECTION );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
        aParameterLists.set( "leader_properties", "PropConductivity,Conductivity" );

        // Dirichlet Nitsche for velocity
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPDirichletNitscheU" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::VELOCITY_DIRICHLET_NITSCHE );
        aParameterLists.set( "function_parameters", "100.0/1.0" );
        aParameterLists.set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
        aParameterLists.set( "leader_properties",
                "PropViscosity,Viscosity;"
                "PropDensity,Density" );

        // Dirichhlet Nitsche for temperature
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPDirichletNitscheT" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "leader_properties", "PropConductivity,Material" );

        if ( isGhost )
        {
            // ghost viscous
            aParameterLists( FEM::STABILIZATION ).add_parameter_list();
            aParameterLists.set( "stabilization_name", "SPGPViscosity" );
            aParameterLists.set( "leader_phase_name", "PhaseFluid" );
            aParameterLists.set( "stabilization_type", fem::Stabilization_Type::VISCOUS_GHOST );
            aParameterLists.set( "function_parameters", "0.05" );
            aParameterLists.set( "leader_properties", "PropViscosity,Viscosity" );

            // ghost convective
            aParameterLists( FEM::STABILIZATION ).add_parameter_list();
            aParameterLists.set( "stabilization_name", "SPGPVelocity" );
            aParameterLists.set( "leader_phase_name", "PhaseFluid" );
            aParameterLists.set( "stabilization_type", fem::Stabilization_Type::CONVECTIVE_GHOST );
            aParameterLists.set( "function_parameters", "0.05" );
            aParameterLists.set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
            aParameterLists.set( "leader_properties", "PropDensity,Density" );

            // ghost pressure
            aParameterLists( FEM::STABILIZATION ).add_parameter_list();
            aParameterLists.set( "stabilization_name", "SPGPPressure" );
            aParameterLists.set( "leader_phase_name", "PhaseFluid" );
            aParameterLists.set( "stabilization_type", fem::Stabilization_Type::PRESSURE_GHOST );
            aParameterLists.set( "function_parameters", "0.005/1.0" );
            aParameterLists.set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
            aParameterLists.set( "leader_properties",
                    "PropViscosity,Viscosity;"
                    "PropDensity,Density" );

            // ghost temperature
            aParameterLists( FEM::STABILIZATION ).add_parameter_list();
            aParameterLists.set( "stabilization_name", "SPGPTemp" );
            aParameterLists.set( "leader_phase_name", "PhaseFluid" );
            aParameterLists.set( "follower_phase_name", "PhaseFluid" );
            aParameterLists.set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
            aParameterLists.set( "function_parameters", "0.005" );
            aParameterLists.set( "leader_properties", "PropConductivity,Material" );
            aParameterLists.set( "follower_properties", "PropConductivity,Material" );
        }

        //------------------------------------------------------------------------------
        // fill the IWG part of the parameter list

        // bulk NS incompressible for velocity
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGVelocityBulk" );
        aParameterLists.set( "IWG_bulk_type", fem::Element_Type::BULK );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_BULK );
        aParameterLists.set( "dof_residual", "VX,VY" );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists.set( "stabilization_parameters", "SPIncFlow,IncompressibleFlow" );

        // bulk NS incompressible for pressure
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGPressureBulk" );
        aParameterLists.set( "IWG_bulk_type", fem::Element_Type::BULK );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_BULK );
        aParameterLists.set( "dof_residual", "P" );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists.set( "stabilization_parameters", "SPIncFlow,IncompressibleFlow" );

        // bulk diffusion
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGDiffusionBulk" );
        aParameterLists.set( "IWG_bulk_type", fem::Element_Type::BULK );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusion,Diffusion" );

        // bulk advection
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGAdvectionBulk" );
        aParameterLists.set( "IWG_bulk_type", fem::Element_Type::BULK );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::ADVECTION_BULK );
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusion,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPSUPGTemp,SUPG" );

        // inlet velocity for velocity
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGInletVelocity" );
        aParameterLists.set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_SYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "VX,VY" );
        aParameterLists.set( "leader_properties", "PropDirichletInletU,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists.set( "stabilization_parameters", "SPDirichletNitscheU,DirichletNitsche" );
        aParameterLists.set( "neighbor_phases", "PhaseInlet" );

        // inlet velocity for pressure
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGInletPressure" );
        aParameterLists.set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "P" );
        aParameterLists.set( "leader_properties", "PropDirichletInletU,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists.set( "neighbor_phases", "PhaseInlet" );

        // inlet temperature
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGInletTemp" );
        aParameterLists.set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_SYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_properties", "PropInletTemp,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusion,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPDirichletNitscheT,DirichletNitsche" );
        aParameterLists.set( "neighbor_phases", "PhaseInlet" );

        // side flux on cylinder
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGCylinderFluxTemp" );
        aParameterLists.set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_NEUMANN );
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_properties", "PropSideFlux,Neumann" );
        aParameterLists.set( "neighbor_phases", "PhaseCylinder" );

        // zero velocity on wall and cylinders for velocity
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGZeroVelocity" );
        aParameterLists.set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_SYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "VX,VY" );
        aParameterLists.set( "leader_properties", "PropDirichletZeroU,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists.set( "stabilization_parameters", "SPDirichletNitscheU,DirichletNitsche" );
        aParameterLists.set( "neighbor_phases", "PhaseCylinder,PhaseWall" );

        // zero velocity on wall and cylinders for pressure
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGZeroPressure" );
        aParameterLists.set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "P" );
        aParameterLists.set( "leader_properties", "PropDirichletZeroU,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists.set( "neighbor_phases", "PhaseCylinder,PhaseWall" );

        if ( isGhost )
        {
            // ghost viscous
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "IWGGPViscous" );
            aParameterLists.set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            aParameterLists.set( "leader_phase_name", "PhaseFluid" );
            aParameterLists.set( "follower_phase_name", "PhaseFluid" );
            aParameterLists.set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists.set( "dof_residual", "VX,VY" );
            aParameterLists.set( "stabilization_parameters", "SPGPViscosity,GhostSP" );

            // ghost convective
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "IWGGPConvective" );
            aParameterLists.set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            aParameterLists.set( "leader_phase_name", "PhaseFluid" );
            aParameterLists.set( "follower_phase_name", "PhaseFluid" );
            aParameterLists.set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists.set( "dof_residual", "VX,VY" );
            aParameterLists.set( "stabilization_parameters", "SPGPVelocity,GhostSP" );

            // ghost pressure
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "IWGGPPressure" );
            aParameterLists.set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            aParameterLists.set( "leader_phase_name", "PhaseFluid" );
            aParameterLists.set( "follower_phase_name", "PhaseFluid" );
            aParameterLists.set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists.set( "dof_residual", "P" );
            aParameterLists.set( "stabilization_parameters", "SPGPPressure,GhostSP" );

            // ghost temperature
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "IWGGPTemp" );
            aParameterLists.set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            aParameterLists.set( "leader_phase_name", "PhaseFluid" );
            aParameterLists.set( "follower_phase_name", "PhaseFluid" );
            aParameterLists.set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists.set( "dof_residual", "TEMP" );
            aParameterLists.set( "stabilization_parameters", "SPGPTemp,GhostSP" );
        }

        //------------------------------------------------------------------------------
        // fill the IQI part of the parameter list

        // dof VX
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkVX" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "VX,VY" );
        aParameterLists.set( "vectorial_field_index", 0 );

        // dof VY
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkVY" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "VX,VY" );
        aParameterLists.set( "vectorial_field_index", 1 );

        // dof P
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkP" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "P" );
        aParameterLists.set( "vectorial_field_index", 0 );

        // dof TEMP
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkTEMP" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "TEMP" );
        aParameterLists.set( "vectorial_field_index", 0 );

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        aParameterLists( FEM::COMPUTATION );
        aParameterLists.set( "print_physics_model", true );
    }

    void
    SOLParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::AMESOS_IMPL );

        aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list();

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "NLA_combined_res_jac_assembly", true );
        aParameterLists.set( "NLA_rel_res_norm_drop", 1e-04 );
        aParameterLists.set( "NLA_relaxation_parameter", 0.80 );
        aParameterLists.set( "NLA_max_iter", 20 );

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "NLA_DofTypes", "VX,VY;P;TEMP" );

        aParameterLists( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list();

        aParameterLists( SOL::TIME_SOLVERS ).add_parameter_list();
        aParameterLists.set( "TSA_DofTypes", "VX,VY,P,TEMP" );
        aParameterLists.set( "TSA_Initialize_Sol_Vec", "VX,1.0;VY,0.0;P,0.0;TEMP,0.0" );
        aParameterLists.set( "TSA_Output_Indices", "0" );
        aParameterLists.set( "TSA_Output_Criteria", "Output_Criterion" );

        aParameterLists( SOL::PRECONDITIONERS ).add_parameter_list( sol::PreconditionerType::NONE );
    }

    void
    MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
    }

    void
    VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "File_Name", std::pair< std::string, std::string >( "./", tExoFile ) );
        aParameterLists.set( "Mesh_Type", vis::VIS_Mesh_Type::STANDARD );
        aParameterLists.set( "Set_Names", "HMR_dummy_n_p0,HMR_dummy_c_p0" );
        aParameterLists.set( "Field_Names",
                "VX,VY,P,TEMP,"
                "IQIBulkVX,IQIBulkVY,IQIBulkP,IQIBulkTEMP" );
        aParameterLists.set( "Field_Type",
                "NODAL,NODAL,NODAL,NODAL,"
                "GLOBAL,GLOBAL,GLOBAL,GLOBAL" );
        aParameterLists.set( "IQI_Names",
                "IQIBulkVX,IQIBulkVY,IQIBulkP,IQIBulkTEMP,"
                "IQIBulkVX,IQIBulkVY,IQIBulkP,IQIBulkTEMP" );
        aParameterLists.set( "Save_Frequency", 1 );
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
