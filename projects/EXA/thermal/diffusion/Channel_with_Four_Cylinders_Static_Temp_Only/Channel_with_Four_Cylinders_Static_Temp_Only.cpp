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
        aParameterLists.set( "is_optimization_problem", false );
    }

    void
    HMRParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "number_of_elements_per_dimension", 176, 88 );
        aParameterLists.set( "processor_decomposition_method", 1 );
        aParameterLists.set( "processor_dimensions", 2, 1 );
        aParameterLists.set( "domain_dimensions", 4.0, 2.0 );
        aParameterLists.set( "domain_offset", -1.24, -0.86 );
        aParameterLists.set( "lagrange_output_meshes", "0" );

        aParameterLists.set( "lagrange_orders", std::to_string( gInterpolationOrder ) );
        aParameterLists.set( "lagrange_pattern", "0" );
        aParameterLists.set( "bspline_orders", std::to_string( gInterpolationOrder ) );
        aParameterLists.set( "bspline_pattern", "0" );

        aParameterLists.set( "lagrange_to_bspline", "0" );

        aParameterLists.set( "refinement_buffer", (int)gInterpolationOrder );
        aParameterLists.set( "staircase_buffer", 1 );
        aParameterLists.set( "initial_refinement", "0" );
        aParameterLists.set( "initial_refinement_pattern", "0" );

        aParameterLists.set( "use_number_aura", 1 );

        aParameterLists.set( "severity_level", 0 );
    }

    void
    XTKParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "decompose", true );
        aParameterLists.set( "decomposition_type", "conformal" );
        aParameterLists.set( "enrich", true );
        aParameterLists.set( "basis_rank", "bspline" );
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

        aParameterLists.set( "number_of_phases", 5 );
        aParameterLists.set( "phase_function_name", "get_phase_index" );
        aParameterLists.set( "output_mesh_file", "GEN_Channel_with_Four_Cylinders_Static_Temp_Only.exo" );

        // Geometry parameter lists
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists.set( "field_function_name", "Func_Bottom_Plane" );
        aParameterLists.set( "number_of_refinements", 1 );
        aParameterLists.set( "refinement_mesh_index", 0 );

        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists.set( "field_function_name", "Func_Top_Plane" );
        aParameterLists.set( "number_of_refinements", 1 );
        aParameterLists.set( "refinement_mesh_index", 0 );

        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists.set( "field_function_name", "Func_Left_Plane" );
        aParameterLists.set( "number_of_refinements", 1 );
        aParameterLists.set( "refinement_mesh_index", 0 );

        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists.set( "field_function_name", "Func_Right_Plane" );
        aParameterLists.set( "number_of_refinements", 1 );
        aParameterLists.set( "refinement_mesh_index", 0 );

        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists.set( "field_function_name", "Func_Cylinder" );
        aParameterLists.set( "number_of_refinements", 1 );
        aParameterLists.set( "refinement_mesh_index", 0 );

        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists.set( "field_function_name", "Func_Cylinder2" );
        aParameterLists.set( "number_of_refinements", 1 );
        aParameterLists.set( "refinement_mesh_index", 0 );

        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists.set( "field_function_name", "Func_Cylinder3" );
        aParameterLists.set( "number_of_refinements", 1 );
        aParameterLists.set( "refinement_mesh_index", 0 );

        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists.set( "field_function_name", "Func_Cylinder4" );
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

        // create parameter list for property 1
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDensity" );
        aParameterLists.set( "function_parameters", "1.0" );

        // create parameter list for property 2
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropCapacity" );
        aParameterLists.set( "function_parameters", "1.0" );

        // create parameter list for property 3
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropConductivity" );
        aParameterLists.set( "function_parameters", "0.00005" );

        // create parameter list for property 5
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropInletTemp" );
        aParameterLists.set( "function_parameters", "0.0" );

        // create parameter list for property 6
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropSideFlux" );
        aParameterLists.set( "function_parameters", "2.0" );

        //------------------------------------------------------------------------------
        // fill the constitutive model part of the parameter list

        // create parameter list for constitutive model 2
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMDiffusion" );
        aParameterLists.set( "phase_name", "PhaseFluid" );
        aParameterLists.set( "constitutive_type",  fem::Constitutive_Type::DIFF_LIN_ISO ) ;
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists.set( "properties",
                "PropConductivity,Conductivity;"
                "PropDensity     ,Density;"
                "PropCapacity    ,HeatCapacity" );

        //------------------------------------------------------------------------------
        // fill the stabilization parameter part of the parameter list

        // create parameter list for stabilization parameter 2
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPDirichletNitscheT" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::DIRICHLET_NITSCHE ) ;
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "leader_properties", "PropConductivity,Material" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );

        // create parameter list for stabilization parameter 8
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPGPTemp" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::GHOST_DISPL ) ;
        aParameterLists.set( "function_parameters", "0.005" );
        aParameterLists.set( "leader_properties", "PropConductivity,Material" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "follower_phase_name", "PhaseFluid" );

        //------------------------------------------------------------------------------
        // fill the IWG part of the parameter list

        // create parameter list for IWG 3
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGDiffusionBulk" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_BULK ) ;
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusion,Diffusion" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );

        // create parameter list for IWG 11
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGInletTemp" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "leader_properties", "PropInletTemp,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusion,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPDirichletNitscheT,DirichletNitsche" );
        aParameterLists.set( "IWG_bulk_type",  fem::Element_Type::SIDESET ) ;
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "neighbor_phases", "PhaseInlet" );

        // create parameter list for IWG 11
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGCylinderFluxTemp" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_NEUMANN ) ;
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "leader_properties", "PropSideFlux,Neumann" );
        aParameterLists.set( "IWG_bulk_type",  fem::Element_Type::SIDESET ) ;
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "neighbor_phases", "PhaseCylinder" );

        if ( isGhost )
        {
            // create parameter list for IWG 16
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "IWGGPTemp" );
            aParameterLists.set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
            aParameterLists.set( "dof_residual", "TEMP" );
            aParameterLists.set( "leader_dof_dependencies", "TEMP" );
            aParameterLists.set( "follower_dof_dependencies", "TEMP" );
            aParameterLists.set( "stabilization_parameters", "SPGPTemp,GhostSP" );
            aParameterLists.set( "IWG_bulk_type",  fem::Element_Type::DOUBLE_SIDESET ) ;
            aParameterLists.set( "leader_phase_name", "PhaseFluid" );
            aParameterLists.set( "follower_phase_name", "PhaseFluid" );
            }

        //------------------------------------------------------------------------------
        // fill the IQI part of the parameter list

        // create parameter list for IQI 3
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkTEMP" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists.set( "dof_quantity", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        aParameterLists( FEM::COMPUTATION );
    }

    void
    SOLParameterList( Module_Parameter_Lists& aParameterLists )
    {

        if ( gTestCaseIndex == 0 )
        {
            aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::PETSC );
        }
        else
        {
            aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::AMESOS_IMPL );
        }

        aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list();

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "NLA_rel_res_norm_drop", 1e-04 );
        aParameterLists.set( "NLA_relaxation_parameter", 1.0 );
        aParameterLists.set( "NLA_max_iter", 20 );

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "NLA_DofTypes", "TEMP" );

        aParameterLists( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list();

        aParameterLists( SOL::TIME_SOLVERS ).add_parameter_list();
        aParameterLists.set( "TSA_DofTypes", "TEMP" );
        aParameterLists.set( "TSA_Initialize_Sol_Vec", "TEMP,0.0" );
        aParameterLists.set( "TSA_Output_Indices", "0" );
        aParameterLists.set( "TSA_Output_Criteria", "Output_Criterion" );

        aParameterLists( SOL::SOLVER_WAREHOUSE );
        if ( gTestCaseIndex == 0 )
        {
            aParameterLists.set( "SOL_TPL_Type",  sol::MapType::Petsc ) ;
        }

        aParameterLists( SOL::PRECONDITIONERS ).add_parameter_list(  sol::PreconditionerType::NONE );
    }

    void
    MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "order_adofs_by_host", false );
    }

    void
    VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "File_Name", std::pair< std::string, std::string >( "./", "Channel_with_Four_Cylinders_Static_Temp_Only.exo" ) );
        aParameterLists.set( "Mesh_Type",  vis::VIS_Mesh_Type::STANDARD ) ;
        aParameterLists.set( "Set_Names", "HMR_dummy_n_p0,HMR_dummy_c_p0" );
        aParameterLists.set( "Field_Names", "TEMP,IQIBulkTEMP" );
        aParameterLists.set( "Field_Type", "NODAL,GLOBAL" );
        aParameterLists.set( "IQI_Names", "IQIBulkTEMP,IQIBulkTEMP" );
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
