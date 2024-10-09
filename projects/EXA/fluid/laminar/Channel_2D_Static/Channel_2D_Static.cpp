/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Channel_2D_Static.cpp
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
#include "parameters.hpp"
#include "fn_equal_to.hpp"

#include "AztecOO.h"

extern bool gInletVelocityBCFlag;
extern bool gInletPressureBCFlag;

#ifdef  __cplusplus
extern "C"
{
#endif
    //------------------------------------------------------------------------------
    namespace moris
    {

        bool sUseGhost          = true;

        std::string sFluid      = "HMR_dummy_c_p10,HMR_dummy_n_p10";
        std::string sFluidGhost = "ghost_p10";
        std::string sInlet      = "iside_b0_10_b1_8";
        std::string sWalls      = "iside_b0_10_b1_2,iside_b0_10_b1_14";

        moris::real tPlaneBottom = -0.5;                  /* y bottom plane (m) */
        moris::real tPlaneTop    =  0.5;                  /* y top plane    (m) */
        moris::real tPlaneLeft   = -2.5;                  /* x left plane   (m) */
        moris::real tPlaneRight  =  2.5;                  /* x right plane  (m) */
        moris::real tChannelRadius = ( tPlaneTop - tPlaneBottom ) / 2.0; /* channel radius  (m) */

        //Material Parameters
        moris::real tFluidDensity   = 1.0; /* Fluid density   () */
        moris::real tFluidViscosity = 1.0; /* Fluid viscosity () */

        // Boundary Conditions
        moris::real tInletPressure  = 20.0;  /* Inlet pressure  () */
        moris::real tGammaNitsche   = 50.0;  /* Penalty for Dirichlet BC */
        moris::real tGammaGPmu      = 0.05;  /* Penalty for ghost viscosity */
        moris::real tGammaGPu       = 0.05;  /* Penalty for ghost velocity */
        moris::real tGammaGPp       = 0.005; /* Penalty for ghost pressure */

        // Constant function for properties
        void Func_Const(
                moris::Matrix< moris::DDRMat >                 & aPropMatrix,
                Vector< moris::Matrix< moris::DDRMat > >  & aParameters,
                moris::fem::Field_Interpolator_Manager         * aFIManager )
        {
            aPropMatrix = aParameters( 0 );
        }

        // Inlet velocity function
        void Func_Inlet_U(
                moris::Matrix< moris::DDRMat >                 & aPropMatrix,
                Vector< moris::Matrix< moris::DDRMat > >  & aParameters,
                moris::fem::Field_Interpolator_Manager         * aFIManager )
        {
            // unpack parameters
            real tRadiusChannel = aParameters( 0 )( 0 );
            real tYChannel      = aParameters( 1 )( 0 );

            // get position in space
            real tY = aFIManager->get_IP_geometry_interpolator()->valx()( 1 );

            // set size for aPropMatrix
            aPropMatrix.set_size( 2, 1, 0.0 );

            // velocity along x direction
            aPropMatrix( 0 ) = - ( tY - ( tYChannel + tRadiusChannel ) ) *
                    ( tY - ( tYChannel - tRadiusChannel ) ) /
                    ( 2.0 * std::pow( tRadiusChannel, 2.0 ) );
        }

        // Output criterion function
        bool Output_Criterion( moris::tsa::Time_Solver * aTimeSolver )
        {
            return true;
        }

        moris::real Func_Plane(
                const moris::Matrix< DDRMat >     & aCoordinates,
                const Vector< real > & aGeometryParameters )
        {
            moris::real tXNormal = aGeometryParameters( 0 );
            moris::real tYNormal = aGeometryParameters( 1 );

            moris::real tXCenter = aGeometryParameters( 2 );
            moris::real tYCenter = aGeometryParameters( 3 );

            moris::real aReturnValue =
                    tXNormal * ( aCoordinates( 0 ) - tXCenter ) + tYNormal * ( aCoordinates( 1 ) - tYCenter );
            return aReturnValue;
        }

        void OPTParameterList( Module_Parameter_Lists & aParameterLists )
        {

                aParameterLists.set( "is_optimization_problem", false);
        }

        void HMRParameterList( Module_Parameter_Lists & aParameterLists )
        {

                aParameterLists.set( "number_of_elements_per_dimension", "200,200");
            aParameterLists.set( "processor_decomposition_method",   1 );
            aParameterLists.set( "processor_dimensions",             "2,1");
            aParameterLists.set( "domain_dimensions",                "10.0,10.0");
            aParameterLists.set( "domain_offset",                    "-4.63,-4.74");
            aParameterLists.set( "domain_sidesets",                  "1,2,3,4");
            aParameterLists.set( "lagrange_output_meshes",           "0");

            aParameterLists.set( "lagrange_orders",  "1" );
            aParameterLists.set( "lagrange_pattern", "0" );
            aParameterLists.set( "bspline_orders",   "1" );
            aParameterLists.set( "bspline_pattern",  "0" );

            aParameterLists.set( "lagrange_to_bspline", "0") ;

            aParameterLists.set( "truncate_bsplines",  1 );
            aParameterLists.set( "refinement_buffer",  3 );
            aParameterLists.set( "staircase_buffer",   3 );
            aParameterLists.set( "initial_refinement", "0" );
            aParameterLists.set( "initial_refinement_pattern", "0" );

            aParameterLists.set( "use_number_aura", 1 );

            aParameterLists.set( "use_multigrid",  0 );
            aParameterLists.set( "severity_level", 0 );

            aParameterLists.set( "adaptive_refinement_level", 0 );
        }

        void XTKParameterList( Module_Parameter_Lists & aParameterLists )
        {
            aParameterLists.set( "decompose",                 true );
            aParameterLists.set( "decomposition_type",        "conformal") ;
            aParameterLists.set( "enrich",                    true );
            aParameterLists.set( "basis_rank",                "bspline") ;
            aParameterLists.set( "enrich_mesh_indices",       "0") ;
            aParameterLists.set( "ghost_stab",                true );
            aParameterLists.set( "multigrid",                 false );
            aParameterLists.set( "print_enriched_ig_mesh",    false );
            aParameterLists.set( "exodus_output_XTK_ig_mesh", false );
        }

        void GENParameterList( Module_Parameter_Lists & aParameterLists )
        {

            // Main GEN parameter list
                // Bottom plane
            aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
            aParameterLists.set( "field_function_name", "Func_Plane");
            aParameterLists( 1 ).insert( "variable_1", 0.0 );
            aParameterLists( 1 ).insert( "variable_2", 1.0 );
            aParameterLists( 1 ).insert( "variable_3", 0.0 );
            aParameterLists( 1 ).insert( "variable_4", tPlaneBottom );

            // Top plane
            aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
            aParameterLists.set( "field_function_name", "Func_Plane");
            aParameterLists( 1 ).insert( "variable_1", 0.0 );
            aParameterLists( 1 ).insert( "variable_2", 1.0 );
            aParameterLists( 1 ).insert( "variable_3", 0.0 );
            aParameterLists( 1 ).insert( "variable_4", tPlaneTop );

            // Left plane
            aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
            aParameterLists.set( "field_function_name", "Func_Plane");
            aParameterLists( 1 ).insert( "variable_1", 1.0 );
            aParameterLists( 1 ).insert( "variable_2", 0.0 );
            aParameterLists( 1 ).insert( "variable_3", tPlaneLeft );
            aParameterLists( 1 ).insert( "variable_4", 0.0 );

            // Right plane
            aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
            aParameterLists.set( "field_function_name", "Func_Plane");
            aParameterLists( 1 ).insert( "variable_1", 1.0 );
            aParameterLists( 1 ).insert( "variable_2", 0.0 );
            aParameterLists( 1 ).insert( "variable_3", tPlaneRight );
            aParameterLists( 1 ).insert( "variable_4", 0.0 );
        }

        void FEMParameterList( Module_Parameter_Lists & aParameterLists )
        {
            aParameterLists.hack_for_legacy_fem();
            //------------------------------------------------------------------------------
            // fill the property part of the parameter list

            // create viscosity property
            aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
            aParameterLists.set( "property_name",            "PropViscosity") ;
            aParameterLists.set( "function_parameters",      std::to_string(tFluidViscosity) );
            aParameterLists.set( "value_function",           "Func_Const") ;

            // create density property
            aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
            aParameterLists.set( "property_name",            "PropDensity") ;
            aParameterLists.set( "function_parameters",      std::to_string(tFluidDensity) );
            aParameterLists.set( "value_function",           "Func_Const") ;

            // create zero velocity property
            aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
            aParameterLists.set( "property_name",            "PropZeroU") ;
            aParameterLists.set( "function_parameters",      "0.0;0.0") ;
            aParameterLists.set( "value_function",           "Func_Const") ;

            // create init velocity property
            aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
            aParameterLists.set( "property_name",            "PropInletU") ;
            aParameterLists.set( "value_function",           "Func_Inlet_U") ;
            aParameterLists.set( "function_parameters",      std::to_string(tChannelRadius) + "/0.0" );

            // create pressure property
            aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
            aParameterLists.set( "property_name",            "InletPressure") ;
            aParameterLists.set( "function_parameters",      std::to_string(tInletPressure) );
            aParameterLists.set( "value_function",           "Func_Const") ;

            // create total pressure property
            aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
            aParameterLists.set( "property_name",            "PropInletUpwind") ;
            aParameterLists.set( "function_parameters",      "1.0") ;
            aParameterLists.set( "value_function",           "Func_Const") ;

            //------------------------------------------------------------------------------
            // fill the constitutive model part of the parameter list

            // create fluid constitutive model
            aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
            aParameterLists.set( "constitutive_name", "CMFluid") ;
            aParameterLists.set( "constitutive_type",  fem::Constitutive_Type::FLUID_INCOMPRESSIBLE ) ;
            aParameterLists.set( "dof_dependencies",  std::pair< std::string, std::string >( "VX,VY;P", "Velocity,Pressure" ) );
            aParameterLists.set( "properties",        "PropViscosity,Viscosity;PropDensity,Density") ;

            //------------------------------------------------------------------------------
            // fill the stabilization parameter part of the parameter list

            // create SUPG stabilization parameter for Navier-Stokes
            aParameterLists( FEM::STABILIZATION ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
            aParameterLists.set( "stabilization_name",      "SPSUPGNS") ;
            aParameterLists.set( "stabilization_type",       fem::Stabilization_Type::INCOMPRESSIBLE_FLOW ) ;
            aParameterLists.set( "function_parameters",     "36.0") ;
            aParameterLists.set( "leader_properties",       "PropViscosity,Viscosity;PropDensity,Density") ;
            aParameterLists.set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY;P", "Velocity,Pressure" ) );

            // create Nitsche stabilization parameter for velocity
            aParameterLists( FEM::STABILIZATION ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
            aParameterLists.set( "stabilization_name",      "SPNitscheU") ;
            aParameterLists.set( "stabilization_type",       fem::Stabilization_Type::VELOCITY_DIRICHLET_NITSCHE ) ;
            aParameterLists.set( "function_parameters",     std::to_string(tGammaNitsche) +"/1.0" );
            aParameterLists.set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
            aParameterLists.set( "leader_properties",       "PropViscosity,Viscosity;PropDensity,Density") ;

            // create Ghost stabilization parameter for viscous velocity
            aParameterLists( FEM::STABILIZATION ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
            aParameterLists.set( "stabilization_name",  "SPGPViscous") ;
            aParameterLists.set( "stabilization_type",   fem::Stabilization_Type::VISCOUS_GHOST ) ;
            aParameterLists.set( "function_parameters", std::to_string(tGammaGPmu) );
            aParameterLists.set( "leader_properties",   "PropViscosity,Viscosity") ;

            // create Ghost stabilization parameter for velocity
            aParameterLists( FEM::STABILIZATION ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
            aParameterLists.set( "stabilization_name",      "SPGPVelocity") ;
            aParameterLists.set( "stabilization_type",       fem::Stabilization_Type::CONVECTIVE_GHOST ) ;
            aParameterLists.set( "function_parameters",     std::to_string(tGammaGPu) );
            aParameterLists.set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
            aParameterLists.set( "leader_properties",       "PropDensity,Density") ;

            // create Ghost stabilization parameter for pressure
            aParameterLists( FEM::STABILIZATION ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
            aParameterLists.set( "stabilization_name",      "SPGPPressure") ;
            aParameterLists.set( "stabilization_type",       fem::Stabilization_Type::PRESSURE_GHOST ) ;
            aParameterLists.set( "function_parameters",     std::to_string(tGammaGPp) +"/1.0" );
            aParameterLists.set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
            aParameterLists.set( "leader_properties",       "PropViscosity,Viscosity;PropDensity,Density") ;

            //------------------------------------------------------------------------------
            // fill the IWG part of the parameter list

            // create incompressible NS velocity bulk IWG
            aParameterLists( FEM::IWG ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists.set( "IWG_name",                   "IWGVelocityBulk") ;
            aParameterLists.set( "IWG_type",                    fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_BULK ) ;
            aParameterLists.set( "dof_residual",               "VX,VY") ;
            aParameterLists.set( "leader_dof_dependencies",    "VX,VY;P") ;
            aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid") ;
            aParameterLists.set( "stabilization_parameters",   "SPSUPGNS,IncompressibleFlow") ;
            aParameterLists.set( "mesh_set_names",             sFluid );

            // create incompressible NS pressure bulk IWG
            aParameterLists( FEM::IWG ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists.set( "IWG_name",                   "IWGPressureBulk") ;
            aParameterLists.set( "IWG_type",                    fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_BULK ) ;
            aParameterLists.set( "dof_residual",               "P") ;
            aParameterLists.set( "leader_dof_dependencies",    "VX,VY;P") ;
            aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid") ;
            aParameterLists.set( "stabilization_parameters",   "SPSUPGNS,IncompressibleFlow") ;
            aParameterLists.set( "mesh_set_names",             sFluid );

            if (gInletPressureBCFlag)
            {
                // create inlet total pressure BC
                aParameterLists( FEM::IWG ).add_parameter_list( prm::create_IWG_parameter_list() );
                aParameterLists.set( "IWG_name",                   "IWGInletImposedPressure") ;
                aParameterLists.set( "IWG_type",                    fem::IWG_Type::INCOMPRESSIBLE_NS_IMPOSED_PRESSURE ) ;
                aParameterLists.set( "dof_residual",               "VX,VY") ;
                aParameterLists.set( "leader_dof_dependencies",    "VX,VY;P") ;
                aParameterLists.set( "leader_properties",          "InletPressure,Pressure") ;
                aParameterLists.set( "mesh_set_names",             sInlet );
                    }

            if (gInletVelocityBCFlag)
            {
                // create incompressible NS velocity Dirichlet IWG for inlet
                aParameterLists( FEM::IWG ).add_parameter_list( prm::create_IWG_parameter_list() );
                aParameterLists.set( "IWG_name",                   "IWGInletVelocity") ;
                aParameterLists.set( "IWG_type",                    fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_SYMMETRIC_NITSCHE ) ;
                aParameterLists.set( "dof_residual",               "VX,VY") ;
                aParameterLists.set( "leader_dof_dependencies",    "VX,VY;P") ;
                aParameterLists.set( "leader_properties",          "PropInletU,Dirichlet") ;
                aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid") ;
                aParameterLists.set( "stabilization_parameters",   "SPNitscheU,DirichletNitsche") ;
                aParameterLists.set( "mesh_set_names",             sInlet );

                // create incompressible NS pressure Dirichlet IWG for inlet
                aParameterLists( FEM::IWG ).add_parameter_list( prm::create_IWG_parameter_list() );
                aParameterLists.set( "IWG_name",                   "IWGInletPressure") ;
                aParameterLists.set( "IWG_type",                    fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_SYMMETRIC_NITSCHE ) ;
                aParameterLists.set( "dof_residual",               "P") ;
                aParameterLists.set( "leader_dof_dependencies",    "VX,VY;P") ;
                aParameterLists.set( "leader_properties",          "PropInletU,Dirichlet") ;
                aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid") ;
                aParameterLists.set( "mesh_set_names",             sInlet );
                    }

            // create incompressible NS velocity Dirichlet IWG for walls
            aParameterLists( FEM::IWG ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists.set( "IWG_name",                   "IWGZeroVelocity") ;
            aParameterLists.set( "IWG_type",                    fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_SYMMETRIC_NITSCHE ) ;
            aParameterLists.set( "dof_residual",               "VX,VY") ;
            aParameterLists.set( "leader_dof_dependencies",    "VX,VY;P") ;
            aParameterLists.set( "leader_properties",          "PropZeroU,Dirichlet") ;
            aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid") ;
            aParameterLists.set( "stabilization_parameters",   "SPNitscheU,DirichletNitsche") ;
            aParameterLists.set( "mesh_set_names",             sWalls );

            // create incompressible NS pressure Dirichlet IWG for walls
            aParameterLists( FEM::IWG ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists.set( "IWG_name",                   "IWGZeroPressure") ;
            aParameterLists.set( "IWG_type",                    fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_SYMMETRIC_NITSCHE ) ;
            aParameterLists.set( "dof_residual",               "P") ;
            aParameterLists.set( "leader_dof_dependencies",    "VX,VY;P") ;
            aParameterLists.set( "leader_properties",          "PropZeroU,Dirichlet") ;
            aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid") ;
            aParameterLists.set( "mesh_set_names",             sWalls );

            if( sUseGhost )
            {
                // create Ghost stabilization viscous IWG
                aParameterLists( FEM::IWG ).add_parameter_list( prm::create_IWG_parameter_list() );
                aParameterLists.set( "IWG_name",                   "IWGGPViscous") ;
                aParameterLists.set( "IWG_type",                    fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
                aParameterLists.set( "dof_residual",               "VX,VY") ;
                aParameterLists.set( "leader_dof_dependencies",    "VX,VY;P") ;
                aParameterLists.set( "follower_dof_dependencies",     "VX,VY;P") ;
                aParameterLists.set( "stabilization_parameters",   "SPGPViscous,GhostSP") ;
                aParameterLists.set( "mesh_set_names",             sFluidGhost );

                // create Ghost stabilization convective IWG
                aParameterLists( FEM::IWG ).add_parameter_list( prm::create_IWG_parameter_list() );
                aParameterLists.set( "IWG_name",                   "IWGGPConvective") ;
                aParameterLists.set( "IWG_type",                    fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
                aParameterLists.set( "dof_residual",               "VX,VY") ;
                aParameterLists.set( "leader_dof_dependencies",    "VX,VY;P") ;
                aParameterLists.set( "follower_dof_dependencies",     "VX,VY;P") ;
                aParameterLists.set( "stabilization_parameters",   "SPGPVelocity,GhostSP") ;
                aParameterLists.set( "mesh_set_names",             sFluidGhost );

                // create Ghost stabilization pressure IWG
                aParameterLists( FEM::IWG ).add_parameter_list( prm::create_IWG_parameter_list() );
                aParameterLists.set( "IWG_name",                   "IWGGPPressure") ;
                aParameterLists.set( "IWG_type",                    fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
                aParameterLists.set( "dof_residual",               "P") ;
                aParameterLists.set( "leader_dof_dependencies",    "VX,VY;P") ;
                aParameterLists.set( "follower_dof_dependencies",     "VX,VY;P") ;
                aParameterLists.set( "stabilization_parameters",   "SPGPPressure,GhostSP") ;
                aParameterLists.set( "mesh_set_names",             sFluidGhost );
                    }

            //------------------------------------------------------------------------------
            // fill the IQI part of the parameter list

            // create parameter list for IQI 1
            aParameterLists( FEM::IQI ).add_parameter_list( prm::create_IQI_parameter_list() );
            aParameterLists.set( "IQI_name",                   "IQIBulkVX") ;
            aParameterLists.set( "IQI_type",                    fem::IQI_Type::DOF ) ;
            aParameterLists.set( "dof_quantity",               "VX,VY");
            aParameterLists.set( "leader_dof_dependencies",    "VX,VY") ;
            aParameterLists.set( "vectorial_field_index",      0 );
            aParameterLists.set( "mesh_set_names",             sFluid );

            // create parameter list for IQI 2
            aParameterLists( FEM::IQI ).add_parameter_list( prm::create_IQI_parameter_list() );
            aParameterLists.set( "IQI_name",                   "IQIBulkVY") ;
            aParameterLists.set( "IQI_type",                    fem::IQI_Type::DOF ) ;
            aParameterLists.set( "dof_quantity",               "VX,VY");
            aParameterLists.set( "leader_dof_dependencies",    "VX,VY") ;
            aParameterLists.set( "vectorial_field_index",      1 );
            aParameterLists.set( "mesh_set_names",             sFluid );

            // create parameter list for IQI 3
            aParameterLists( FEM::IQI ).add_parameter_list( prm::create_IQI_parameter_list() );
            aParameterLists.set( "IQI_name",                   "IQIBulkP") ;
            aParameterLists.set( "IQI_type",                    fem::IQI_Type::DOF ) ;
            aParameterLists.set( "dof_quantity",               "P");
            aParameterLists.set( "leader_dof_dependencies",    "P") ;
            aParameterLists.set( "vectorial_field_index",      0 );
            aParameterLists.set( "mesh_set_names",             sFluid );

            //------------------------------------------------------------------------------
            // fill the computation part of the parameter list
            aParameterLists( FEM::COMPUTATION );
        }

        void SOLParameterList( Module_Parameter_Lists & aParameterLists )
        {

            aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL ) );

            aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list( moris::prm::create_linear_solver_parameter_list() );

            aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
            aParameterLists.set("NLA_rel_res_norm_drop",    1e-06 );
            aParameterLists.set("NLA_relaxation_parameter", 1.0 );
            aParameterLists.set("NLA_max_iter",             50 );

            aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
            aParameterLists.set("NLA_DofTypes", "VX,VY;P") ;

            aParameterLists( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list( moris::prm::create_time_solver_algorithm_parameter_list() );

            aParameterLists( SOL::TIME_SOLVERS ).add_parameter_list( moris::prm::create_time_solver_parameter_list() );
            aParameterLists.set("TSA_DofTypes",            "VX,VY;P") ;
            aParameterLists.set("TSA_Initialize_Sol_Vec",  "VX,0.5;VY,0.0;P,0.0") ;
            aParameterLists.set("TSA_Output_Indices",      "0") ;
            aParameterLists.set("TSA_Output_Criteria",      "Output_Criterion") ;

                aParameterLists( SOL::PRECONDITIONERS ).add_parameter_list( moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::NONE ) );
        }

        void MSIParameterList( Module_Parameter_Lists & aParameterLists )
        {

            }

        void VISParameterList( Module_Parameter_Lists & aParameterLists )
        {

                aParameterLists.set( "File_Name"  , std::pair< std::string, std::string >( "./", "Channel_2D_Static.exo" ) );
            aParameterLists.set( "Mesh_Type"  ,  vis::VIS_Mesh_Type::STANDARD ) ;
            aParameterLists.set( "Set_Names"  , sFluid );
            aParameterLists.set( "Field_Names",  "VX,VY,P" ) ;
            aParameterLists.set( "Field_Type" ,  "NODAL,NODAL,NODAL" ) ;
            aParameterLists.set( "IQI_Names"  ,  "IQIBulkVX,IQIBulkVY,IQIBulkP" ) ;
            aParameterLists.set( "Save_Frequency",1);
        }

        void MORISGENERALParameterList( Module_Parameter_Lists & aParameterLists )
        {

        }

        //------------------------------------------------------------------------------
    }

    //------------------------------------------------------------------------------
#ifdef  __cplusplus
}
#endif

