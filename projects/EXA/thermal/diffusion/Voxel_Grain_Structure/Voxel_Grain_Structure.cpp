/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Voxel_Grain_Structure.cpp
 *
 */

#include <string>
#include <iostream>

#include "paths.hpp"

#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "GEN_Data_Types.hpp"
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
    bool isGhost     = false;
    bool isMultigrid = false;

    uint tNumBulkPhases     = 10;
    uint tNumBulkPhasesHalf = std::floor( tNumBulkPhases / 2 );

    // Minimum value to be returned by level set function
    moris::real tMinLSvalue = 1.0e-6;

    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------

    std::string tBulk_1       = "";
    std::string tBulk_2       = "";
    std::string tBulk         = "";
    std::string tRightBC_1    = "";
    std::string tRightBC_2    = "";
    std::string tLeftBC       = "";
    std::string tInterface_11 = "";
    std::string tInterface_22 = "";
    std::string tInterface_21 = "";
    std::string tInterface_12 = "";
    std::string tGhost_1      = "";
    std::string tGhost_2      = "";

    void
    Name_func()
    {
        tBulk_1 = "";
        tBulk_2 = "";

        for ( uint Ik = 1; Ik < tNumBulkPhasesHalf + 1; Ik++ )
        {
            tBulk_1  = tBulk_1 + "HMR_dummy_n_p" + std::to_string( Ik ) + "," + "HMR_dummy_c_p" + std::to_string( Ik ) + ",";
            tGhost_1 = tGhost_1 + "ghost_p" + std::to_string( Ik ) + ",";
        }

        for ( uint Ik = tNumBulkPhasesHalf + 1; Ik < tNumBulkPhases + 1; Ik++ )
        {
            tBulk_2  = tBulk_2 + "HMR_dummy_n_p" + std::to_string( Ik ) + "," + "HMR_dummy_c_p" + std::to_string( Ik ) + ",";
            tGhost_2 = tGhost_2 + "ghost_p" + std::to_string( Ik ) + ",";
        }

        tBulk_1.pop_back();
        tBulk_2.pop_back();

        tGhost_1.pop_back();
        tGhost_2.pop_back();

        tBulk = tBulk_1 + "," + tBulk_2;

        //--------------------------------------------------------------------------------------------------------

        tRightBC_1 = "";

        for ( uint Ik = 1; Ik < tNumBulkPhasesHalf + 1; Ik++ )
        {
            tRightBC_1 = tRightBC_1 + "SideSet_2_c_p" + std::to_string( Ik ) + "," + "SideSet_2_n_p" + std::to_string( Ik ) + ",";
        }
        for ( uint Ik = tNumBulkPhasesHalf + 1; Ik < tNumBulkPhases + 1; Ik++ )
        {
            tRightBC_2 = tRightBC_2 + "SideSet_2_c_p" + std::to_string( Ik ) + "," + "SideSet_2_n_p" + std::to_string( Ik ) + ",";
        }
        tRightBC_1.pop_back();
        tRightBC_2.pop_back();

        //--------------------------------------------------------------------------------------------------------

        tLeftBC = "";

        for ( uint Ik = 1; Ik < tNumBulkPhases + 1; Ik++ )
        {
            tLeftBC = tLeftBC + "SideSet_4_c_p" + std::to_string( Ik ) + "," + "SideSet_4_n_p" + std::to_string( Ik ) + ",";
        }
        tLeftBC.pop_back();
        //--------------------------------------------------------------------------------------------------------

        tInterface_11 = "";
        tInterface_22 = "";
        tInterface_21 = "";
        tInterface_12 = "";

        for ( uint Ik = 1; Ik < tNumBulkPhases + 1; Ik++ )
        {
            for ( uint Ii = Ik; Ii < tNumBulkPhases + 1; Ii++ )
            {
                if ( Ik != Ii )
                {
                    if ( Ik < tNumBulkPhasesHalf + 1 )
                    {
                        if ( Ii < tNumBulkPhasesHalf + 1 )
                        {
                            tInterface_11 = tInterface_11 + "dbl_iside_p0_" + std::to_string( Ik ) + "_p1_" + std::to_string( Ii ) + ",";
                        }
                        else
                        {
                            tInterface_12 = tInterface_12 + "dbl_iside_p0_" + std::to_string( Ik ) + "_p1_" + std::to_string( Ii ) + ",";
                        }
                    }
                    else
                    {
                        if ( Ii < tNumBulkPhasesHalf + 1 )
                        {
                            std::cout << "Ik: " << Ik << ""
                                      << " Ii: " << Ii << '\n';
                            tInterface_21 = tInterface_21 + "dbl_iside_p0_" + std::to_string( Ik ) + "_p1_" + std::to_string( Ii ) + ",";
                            std::cout << tInterface_21 << Ii << '\n';
                        }
                        else
                        {
                            tInterface_22 = tInterface_22 + "dbl_iside_p0_" + std::to_string( Ik ) + "_p1_" + std::to_string( Ii ) + ",";
                        }
                    }
                }
            }
        }
        tInterface_11.pop_back();
        tInterface_22.pop_back();
        tInterface_12.pop_back();
    }

    //--------------------------------------------------------------------------------------------------------

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

    uint
    get_phase_index( const moris::gen::Geometry_Bitset& aGeometrySigns )
    {
        uint tPhaseIndex = 0;

        for ( uint tGeometryIndex = 1; tGeometryIndex < tNumBulkPhases + 1; tGeometryIndex++ )
        {
            if ( aGeometrySigns.test( tGeometryIndex ) )
            {
                tPhaseIndex = tGeometryIndex;
                break;
            }
        }

        return tPhaseIndex;
    }

    void
    OPTParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "is_optimization_problem", false );
    }

    void
    HMRParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "number_of_elements_per_dimension", 39, 39, 39 );
        aParameterLists.set( "domain_dimensions", "4.0,4.0,4.0" );
        aParameterLists.set( "domain_offset", "-2.0,-2.0,-2.0" );
        aParameterLists.set( "lagrange_output_meshes", "0" );

        aParameterLists.set( "lagrange_orders", "1" );
        aParameterLists.set( "lagrange_pattern", "0" );
        aParameterLists.set( "bspline_orders", "1" );
        aParameterLists.set( "bspline_pattern", "0" );

        aParameterLists.set( "lagrange_to_bspline", "0" );

        aParameterLists.set( "truncate_bsplines", 1 );
        aParameterLists.set( "refinement_buffer", 1 );
        aParameterLists.set( "staircase_buffer", 1 );
        aParameterLists.set( "initial_refinement", "0" );
        aParameterLists.set( "initial_refinement_pattern", "0" );

        aParameterLists.set( "use_number_aura", 1 );

        aParameterLists.set( "use_multigrid", 0 );
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
        aParameterLists.set( "multigrid", isMultigrid );
        aParameterLists.set( "verbose", true );
        aParameterLists.set( "print_enriched_ig_mesh", false );
        aParameterLists.set( "exodus_output_XTK_ig_mesh", true );
        aParameterLists.set( "write_cluster_measures_to_exo", false );
    }

    void
    GENParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "number_of_phases", (sint)( tNumBulkPhases + 1 ) );
        aParameterLists.set( "phase_function_name", "get_phase_index" );
        aParameterLists.set( "output_mesh_file", "Voxel_Grain_Structure_Gen.exo" );

        // Geometry parameter lists
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_voxel_geometry_parameter_list() );

        std::string tMORISROOT = moris::get_base_moris_dir();
        std::string tPath      = tMORISROOT + "projects/EXA/thermal/diffusion/Voxel_Grain_Structure/Voxel_file";

        aParameterLists.set( "voxel_field_file", tPath );
        aParameterLists.set( "domain_dimensions", 4.0, 4.0, 4.0 );
        aParameterLists.set( "domain_offset", -2.0, -2.0, -2.0 );
    }

    void
    FEMParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.hack_for_legacy_fem();
        Name_func();
        // create a cell of cell of parameter list for fem

        //------------------------------------------------------------------------------
        // fill the property part of the parameter list

        // create parameter list for property 1
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDensity" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 2
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropCapacity" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 3
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropConductivity_1" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 3
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropConductivity_2" );
        aParameterLists.set( "function_parameters", "10.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 5
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropInletTemp" );
        aParameterLists.set( "function_parameters", "2.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 6
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropSideFlux" );
        aParameterLists.set( "function_parameters", "2.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        //------------------------------------------------------------------------------
        // fill the constitutive model part of the parameter list

        // create parameter list for constitutive model 2
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMDiffusion_1" );
        aParameterLists.set( "constitutive_type",  fem::Constitutive_Type::DIFF_LIN_ISO ) ;
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists.set( "properties",
                "PropConductivity_1,Conductivity;"
                "PropDensity     ,Density;"
                "PropCapacity    ,HeatCapacity" );

        // create parameter list for constitutive model 2
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMDiffusion_2" );
        aParameterLists.set( "constitutive_type",  fem::Constitutive_Type::DIFF_LIN_ISO ) ;
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists.set( "properties",
                "PropConductivity_1,Conductivity;"
                "PropDensity     ,Density;"
                "PropCapacity    ,HeatCapacity" );

        // create parameter list for constitutive model 2
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMDiffusion_3" );
        aParameterLists.set( "constitutive_type",  fem::Constitutive_Type::DIFF_LIN_ISO ) ;
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists.set( "properties",
                "PropConductivity_2,Conductivity;"
                "PropDensity     ,Density;"
                "PropCapacity    ,HeatCapacity" );

        // create parameter list for constitutive model 2
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMDiffusion_4" );
        aParameterLists.set( "constitutive_type",  fem::Constitutive_Type::DIFF_LIN_ISO ) ;
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists.set( "properties",
                "PropConductivity_2,Conductivity;"
                "PropDensity     ,Density;"
                "PropCapacity    ,HeatCapacity" );

        //------------------------------------------------------------------------------
        // fill the stabilization parameter part of the parameter list

        // create parameter list for stabilization parameter 2
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPDirichletNitscheT_1" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::DIRICHLET_NITSCHE ) ;
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "leader_properties", "PropConductivity_1,Material" );

        // create parameter list for stabilization parameter 2
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPDirichletNitscheT_2" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::DIRICHLET_NITSCHE ) ;
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "leader_properties", "PropConductivity_2,Material" );

        // create parameter list for stabilization parameter 8
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPGPTemp_1" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::GHOST_DISPL ) ;
        aParameterLists.set( "function_parameters", "0.005" );
        aParameterLists.set( "leader_properties", "PropConductivity_1,Material" );

        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPGPTemp_2" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::GHOST_DISPL ) ;
        aParameterLists.set( "function_parameters", "0.005" );
        aParameterLists.set( "leader_properties", "PropConductivity_2,Material" );

        // Temperature - Shell - PCM
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", std::string( "SPInterfaceNitsche_11" ) );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::NITSCHE_INTERFACE ) ;
        aParameterLists.set( "function_parameters", std::string( "100.0" ) );
        aParameterLists.set( "leader_properties", std::string( "PropConductivity_1,Material" ) );
        aParameterLists.set( "follower_properties", std::string( "PropConductivity_1,Material" ) );

        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", std::string( "SPInterfaceNitsche_12" ) );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::NITSCHE_INTERFACE ) ;
        aParameterLists.set( "function_parameters", std::string( "100.0" ) );
        aParameterLists.set( "leader_properties", std::string( "PropConductivity_1,Material" ) );
        aParameterLists.set( "follower_properties", std::string( "PropConductivity_2,Material" ) );

        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", std::string( "SPInterfaceNitsche_22" ) );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::NITSCHE_INTERFACE ) ;
        aParameterLists.set( "function_parameters", std::string( "100.0" ) );
        aParameterLists.set( "leader_properties", std::string( "PropConductivity_2,Material" ) );
        aParameterLists.set( "follower_properties", std::string( "PropConductivity_2,Material" ) );

        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", std::string( "SPInterfaceNitsche_21" ) );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::NITSCHE_INTERFACE ) ;
        aParameterLists.set( "function_parameters", std::string( "100.0" ) );
        aParameterLists.set( "leader_properties", std::string( "PropConductivity_2,Material" ) );
        aParameterLists.set( "follower_properties", std::string( "PropConductivity_1,Material" ) );

        //------------------------------------------------------------------------------
        // fill the IWG part of the parameter list

        // create parameter list for IWG 3
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGDiffusionBulk_1" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_BULK ) ;
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusion_1,Diffusion" );
        aParameterLists.set( "mesh_set_names", tBulk_1 );

        // create parameter list for IWG 3
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGDiffusionBulk_2" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_BULK ) ;
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusion_3,Diffusion" );
        aParameterLists.set( "mesh_set_names", tBulk_2 );

        // create parameter list for IWG 11
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGInletTemp_1" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "leader_properties", "PropInletTemp,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusion_1,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPDirichletNitscheT_1,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", tRightBC_1 );

        // create parameter list for IWG 11
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGInletTemp_2" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "leader_properties", "PropInletTemp,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusion_3,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPDirichletNitscheT_2,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", tRightBC_2 );

        // create parameter list for IWG 11
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGCylinderFluxTemp" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_NEUMANN ) ;
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "leader_properties", "PropSideFlux,Neumann" );
        aParameterLists.set( "mesh_set_names", tLeftBC );

        // Temperature - Shell - Shell
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", std::string( "IWGInterfaceShellShellTEMP_1" ) );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_INTERFACE_UNSYMMETRIC_NITSCHE ) ;
        aParameterLists.set( "dof_residual", std::string( "TEMP" ) );
        aParameterLists.set( "leader_dof_dependencies", std::string( "TEMP" ) );
        aParameterLists.set( "follower_dof_dependencies", std::string( "TEMP" ) );
        aParameterLists.set( "leader_constitutive_models", std::string( "CMDiffusion_1,Diffusion" ) );
        aParameterLists.set( "follower_constitutive_models", std::string( "CMDiffusion_2,Diffusion" ) );
        aParameterLists.set( "stabilization_parameters", std::string( "SPInterfaceNitsche_11     ,NitscheInterface" ) );
        aParameterLists.set( "mesh_set_names", tInterface_11 );

        // Temperature - Shell - Shell
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", std::string( "IWGInterfaceShellShellTEMP_2" ) );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_INTERFACE_UNSYMMETRIC_NITSCHE ) ;
        aParameterLists.set( "dof_residual", std::string( "TEMP" ) );
        aParameterLists.set( "leader_dof_dependencies", std::string( "TEMP" ) );
        aParameterLists.set( "follower_dof_dependencies", std::string( "TEMP" ) );
        aParameterLists.set( "leader_constitutive_models", std::string( "CMDiffusion_1,Diffusion" ) );
        aParameterLists.set( "follower_constitutive_models", std::string( "CMDiffusion_3,Diffusion" ) );
        aParameterLists.set( "stabilization_parameters", std::string( "SPInterfaceNitsche_12     ,NitscheInterface" ) );
        aParameterLists.set( "mesh_set_names", tInterface_12 );

        // Temperature - Shell - Shell
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", std::string( "IWGInterfaceShellShellTEMP_3" ) );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_INTERFACE_UNSYMMETRIC_NITSCHE ) ;
        aParameterLists.set( "dof_residual", std::string( "TEMP" ) );
        aParameterLists.set( "leader_dof_dependencies", std::string( "TEMP" ) );
        aParameterLists.set( "follower_dof_dependencies", std::string( "TEMP" ) );
        aParameterLists.set( "leader_constitutive_models", std::string( "CMDiffusion_3,Diffusion" ) );
        aParameterLists.set( "follower_constitutive_models", std::string( "CMDiffusion_4,Diffusion" ) );
        aParameterLists.set( "stabilization_parameters", std::string( "SPInterfaceNitsche_22     ,NitscheInterface" ) );
        aParameterLists.set( "mesh_set_names", tInterface_22 );

        // Temperature - Shell - Shell
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", std::string( "IWGInterfaceShellShellTEMP_4" ) );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_INTERFACE_UNSYMMETRIC_NITSCHE ) ;
        aParameterLists.set( "dof_residual", std::string( "TEMP" ) );
        aParameterLists.set( "leader_dof_dependencies", std::string( "TEMP" ) );
        aParameterLists.set( "follower_dof_dependencies", std::string( "TEMP" ) );
        aParameterLists.set( "leader_constitutive_models", std::string( "CMDiffusion_3,Diffusion" ) );
        aParameterLists.set( "follower_constitutive_models", std::string( "CMDiffusion_1,Diffusion" ) );
        aParameterLists.set( "stabilization_parameters", std::string( "SPInterfaceNitsche_21     ,NitscheInterface" ) );
        aParameterLists.set( "mesh_set_names", tInterface_21 );

        if ( isGhost )
        {
            // create parameter list for IWG 16
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "GPenTemp_1" );
            aParameterLists.set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
            aParameterLists.set( "dof_residual", "TEMP" );
            aParameterLists.set( "leader_dof_dependencies", "TEMP" );
            aParameterLists.set( "follower_dof_dependencies", "TEMP" );
            aParameterLists.set( "stabilization_parameters", "SPGPTemp_1,GhostSP" );
            aParameterLists.set( "mesh_set_names", tGhost_1 );

            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "GPenTemp_2" );
            aParameterLists.set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
            aParameterLists.set( "dof_residual", "TEMP" );
            aParameterLists.set( "leader_dof_dependencies", "TEMP" );
            aParameterLists.set( "follower_dof_dependencies", "TEMP" );
            aParameterLists.set( "stabilization_parameters", "SPGPTemp_2,GhostSP" );
            aParameterLists.set( "mesh_set_names", tGhost_2 );
            }

        //------------------------------------------------------------------------------
        // fill the IQI part of the parameter list

        // create parameter list for IQI 3
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkTEMP" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists.set( "dof_quantity", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "mesh_set_names", tBulk );

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        aParameterLists( FEM::COMPUTATION );
    }

    void
    SOLParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::AMESOS_IMPL );
#ifdef MORIS_USE_MUMPS
        aParameterLists.set( "Solver_Type", "Amesos_Mumps" );
#else
        aParameterLists.set( "Solver_Type", "Amesos_Superludist" );
#endif

        aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list();

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "NLA_rel_res_norm_drop", 1e-04 );
        aParameterLists.set( "NLA_relaxation_parameter", 1.0 );
        aParameterLists.set( "NLA_max_iter", 1 );

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "NLA_DofTypes", "TEMP" );

        aParameterLists( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list();

        aParameterLists( SOL::TIME_SOLVERS ).add_parameter_list();
        aParameterLists.set( "TSA_DofTypes", "TEMP" );
        aParameterLists.set( "TSA_Initialize_Sol_Vec", "TEMP,0.0" );
        aParameterLists.set( "TSA_Output_Indices", "0" );
        aParameterLists.set( "TSA_Output_Criteria", "Output_Criterion" );

        aParameterLists( SOL::PRECONDITIONERS ).add_parameter_list(  sol::PreconditionerType::NONE );
    }

    void
    MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "order_adofs_by_host", false );
        aParameterLists.set( "multigrid", isMultigrid );
    }

    void
    VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "File_Name", std::pair< std::string, std::string >( "./", "Voxel_Grain_Structure.exo" ) );
        aParameterLists.set( "Mesh_Type",  vis::VIS_Mesh_Type::STANDARD ) ;
        aParameterLists.set( "Set_Names", tBulk );
        aParameterLists.set( "Field_Names", "TEMP" );
        aParameterLists.set( "Field_Type", "NODAL" );
        aParameterLists.set( "IQI_Names", "IQIBulkTEMP" );
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
