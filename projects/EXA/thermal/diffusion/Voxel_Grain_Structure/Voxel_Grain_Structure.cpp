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
        aParameterLists( 0 ).add_parameter_list( prm::create_opt_problem_parameter_list() );

        aParameterLists( 0 ).set( "is_optimization_problem", false );
    }

    void
    HMRParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_hmr_parameter_list() );

        aParameterLists( 0 ).set( "number_of_elements_per_dimension", "39,39,39" );
        aParameterLists( 0 ).set( "domain_dimensions", "4.0,4.0,4.0" );
        aParameterLists( 0 ).set( "domain_offset", "-2.0,-2.0,-2.0" );
        aParameterLists( 0 ).set( "domain_sidesets", "1,2,3,4,5,6" );
        aParameterLists( 0 ).set( "lagrange_output_meshes", "0" );

        aParameterLists( 0 ).set( "lagrange_orders", "1" );
        aParameterLists( 0 ).set( "lagrange_pattern", "0" );
        aParameterLists( 0 ).set( "bspline_orders", "1" );
        aParameterLists( 0 ).set( "bspline_pattern", "0" );

        aParameterLists( 0 ).set( "lagrange_to_bspline", "0" );

        aParameterLists( 0 ).set( "truncate_bsplines", 1 );
        aParameterLists( 0 ).set( "refinement_buffer", 1 );
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
        aParameterLists( 0 ).set( "multigrid", isMultigrid );
        aParameterLists( 0 ).set( "verbose", true );
        aParameterLists( 0 ).set( "print_enriched_ig_mesh", false );
        aParameterLists( 0 ).set( "exodus_output_XTK_ig_mesh", true );
        aParameterLists( 0 ).set( "write_cluster_measures_to_exo", false );
    }

    void
    GENParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_gen_parameter_list() );
        aParameterLists( 0 ).set( "number_of_phases", (sint)( tNumBulkPhases + 1 ) );
        aParameterLists( 0 ).set( "phase_function_name", "get_phase_index" );
        aParameterLists( 0 ).set( "output_mesh_file", "Voxel_Grain_Structure_Gen.exo" );

        // Geometry parameter lists
        aParameterLists( 1 ).add_parameter_list( prm::create_voxel_geometry_parameter_list() );

        std::string tMORISROOT = moris::get_base_moris_dir();
        std::string tPath      = tMORISROOT + "projects/EXA/thermal/diffusion/Voxel_Grain_Structure/Voxel_file";

        aParameterLists( 1 ).set( "voxel_field_file", tPath );
        aParameterLists( 1 ).set( "domain_dimensions", 4.0, 4.0, 4.0 );
        aParameterLists( 1 ).set( "domain_offset", -2.0, -2.0, -2.0 );
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
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropDensity" );
        aParameterLists( 0 ).set( "function_parameters", "1.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // create parameter list for property 2
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropCapacity" );
        aParameterLists( 0 ).set( "function_parameters", "1.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // create parameter list for property 3
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropConductivity_1" );
        aParameterLists( 0 ).set( "function_parameters", "1.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // create parameter list for property 3
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropConductivity_2" );
        aParameterLists( 0 ).set( "function_parameters", "10.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // create parameter list for property 5
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropInletTemp" );
        aParameterLists( 0 ).set( "function_parameters", "2.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // create parameter list for property 6
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropSideFlux" );
        aParameterLists( 0 ).set( "function_parameters", "2.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        //------------------------------------------------------------------------------
        // fill the constitutive model part of the parameter list

        // create parameter list for constitutive model 2
        aParameterLists( 1 ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists( 1 ).set( "constitutive_name", "CMDiffusion_1" );
        aParameterLists( 1 ).set( "constitutive_type",  fem::Constitutive_Type::DIFF_LIN_ISO ) ;
        aParameterLists( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists( 1 ).set( "properties",
                "PropConductivity_1,Conductivity;"
                "PropDensity     ,Density;"
                "PropCapacity    ,HeatCapacity" );

        // create parameter list for constitutive model 2
        aParameterLists( 1 ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists( 1 ).set( "constitutive_name", "CMDiffusion_2" );
        aParameterLists( 1 ).set( "constitutive_type",  fem::Constitutive_Type::DIFF_LIN_ISO ) ;
        aParameterLists( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists( 1 ).set( "properties",
                "PropConductivity_1,Conductivity;"
                "PropDensity     ,Density;"
                "PropCapacity    ,HeatCapacity" );

        // create parameter list for constitutive model 2
        aParameterLists( 1 ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists( 1 ).set( "constitutive_name", "CMDiffusion_3" );
        aParameterLists( 1 ).set( "constitutive_type",  fem::Constitutive_Type::DIFF_LIN_ISO ) ;
        aParameterLists( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists( 1 ).set( "properties",
                "PropConductivity_2,Conductivity;"
                "PropDensity     ,Density;"
                "PropCapacity    ,HeatCapacity" );

        // create parameter list for constitutive model 2
        aParameterLists( 1 ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists( 1 ).set( "constitutive_name", "CMDiffusion_4" );
        aParameterLists( 1 ).set( "constitutive_type",  fem::Constitutive_Type::DIFF_LIN_ISO ) ;
        aParameterLists( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists( 1 ).set( "properties",
                "PropConductivity_2,Conductivity;"
                "PropDensity     ,Density;"
                "PropCapacity    ,HeatCapacity" );

        //------------------------------------------------------------------------------
        // fill the stabilization parameter part of the parameter list

        // create parameter list for stabilization parameter 2
        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPDirichletNitscheT_1" );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::DIRICHLET_NITSCHE ) ;
        aParameterLists( 2 ).set( "function_parameters", "100.0" );
        aParameterLists( 2 ).set( "leader_properties", "PropConductivity_1,Material" );

        // create parameter list for stabilization parameter 2
        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPDirichletNitscheT_2" );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::DIRICHLET_NITSCHE ) ;
        aParameterLists( 2 ).set( "function_parameters", "100.0" );
        aParameterLists( 2 ).set( "leader_properties", "PropConductivity_2,Material" );

        // create parameter list for stabilization parameter 8
        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPGPTemp_1" );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::GHOST_DISPL ) ;
        aParameterLists( 2 ).set( "function_parameters", "0.005" );
        aParameterLists( 2 ).set( "leader_properties", "PropConductivity_1,Material" );

        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPGPTemp_2" );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::GHOST_DISPL ) ;
        aParameterLists( 2 ).set( "function_parameters", "0.005" );
        aParameterLists( 2 ).set( "leader_properties", "PropConductivity_2,Material" );

        // Temperature - Shell - PCM
        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", std::string( "SPInterfaceNitsche_11" ) );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::NITSCHE_INTERFACE ) ;
        aParameterLists( 2 ).set( "function_parameters", std::string( "100.0" ) );
        aParameterLists( 2 ).set( "leader_properties", std::string( "PropConductivity_1,Material" ) );
        aParameterLists( 2 ).set( "follower_properties", std::string( "PropConductivity_1,Material" ) );

        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", std::string( "SPInterfaceNitsche_12" ) );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::NITSCHE_INTERFACE ) ;
        aParameterLists( 2 ).set( "function_parameters", std::string( "100.0" ) );
        aParameterLists( 2 ).set( "leader_properties", std::string( "PropConductivity_1,Material" ) );
        aParameterLists( 2 ).set( "follower_properties", std::string( "PropConductivity_2,Material" ) );

        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", std::string( "SPInterfaceNitsche_22" ) );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::NITSCHE_INTERFACE ) ;
        aParameterLists( 2 ).set( "function_parameters", std::string( "100.0" ) );
        aParameterLists( 2 ).set( "leader_properties", std::string( "PropConductivity_2,Material" ) );
        aParameterLists( 2 ).set( "follower_properties", std::string( "PropConductivity_2,Material" ) );

        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", std::string( "SPInterfaceNitsche_21" ) );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::NITSCHE_INTERFACE ) ;
        aParameterLists( 2 ).set( "function_parameters", std::string( "100.0" ) );
        aParameterLists( 2 ).set( "leader_properties", std::string( "PropConductivity_2,Material" ) );
        aParameterLists( 2 ).set( "follower_properties", std::string( "PropConductivity_1,Material" ) );

        //------------------------------------------------------------------------------
        // fill the IWG part of the parameter list

        // create parameter list for IWG 3
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGDiffusionBulk_1" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_BULK ) ;
        aParameterLists( 3 ).set( "dof_residual", "TEMP" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMDiffusion_1,Diffusion" );
        aParameterLists( 3 ).set( "mesh_set_names", tBulk_1 );

        // create parameter list for IWG 3
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGDiffusionBulk_2" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_BULK ) ;
        aParameterLists( 3 ).set( "dof_residual", "TEMP" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMDiffusion_3,Diffusion" );
        aParameterLists( 3 ).set( "mesh_set_names", tBulk_2 );

        // create parameter list for IWG 11
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGInletTemp_1" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        aParameterLists( 3 ).set( "dof_residual", "TEMP" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 3 ).set( "leader_properties", "PropInletTemp,Dirichlet" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMDiffusion_1,Diffusion" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPDirichletNitscheT_1,DirichletNitsche" );
        aParameterLists( 3 ).set( "mesh_set_names", tRightBC_1 );

        // create parameter list for IWG 11
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGInletTemp_2" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        aParameterLists( 3 ).set( "dof_residual", "TEMP" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 3 ).set( "leader_properties", "PropInletTemp,Dirichlet" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMDiffusion_3,Diffusion" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPDirichletNitscheT_2,DirichletNitsche" );
        aParameterLists( 3 ).set( "mesh_set_names", tRightBC_2 );

        // create parameter list for IWG 11
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGCylinderFluxTemp" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_NEUMANN ) ;
        aParameterLists( 3 ).set( "dof_residual", "TEMP" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 3 ).set( "leader_properties", "PropSideFlux,Neumann" );
        aParameterLists( 3 ).set( "mesh_set_names", tLeftBC );

        // Temperature - Shell - Shell
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", std::string( "IWGInterfaceShellShellTEMP_1" ) );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_INTERFACE_UNSYMMETRIC_NITSCHE ) ;
        aParameterLists( 3 ).set( "dof_residual", std::string( "TEMP" ) );
        aParameterLists( 3 ).set( "leader_dof_dependencies", std::string( "TEMP" ) );
        aParameterLists( 3 ).set( "follower_dof_dependencies", std::string( "TEMP" ) );
        aParameterLists( 3 ).set( "leader_constitutive_models", std::string( "CMDiffusion_1,Diffusion" ) );
        aParameterLists( 3 ).set( "follower_constitutive_models", std::string( "CMDiffusion_2,Diffusion" ) );
        aParameterLists( 3 ).set( "stabilization_parameters", std::string( "SPInterfaceNitsche_11     ,NitscheInterface" ) );
        aParameterLists( 3 ).set( "mesh_set_names", tInterface_11 );

        // Temperature - Shell - Shell
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", std::string( "IWGInterfaceShellShellTEMP_2" ) );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_INTERFACE_UNSYMMETRIC_NITSCHE ) ;
        aParameterLists( 3 ).set( "dof_residual", std::string( "TEMP" ) );
        aParameterLists( 3 ).set( "leader_dof_dependencies", std::string( "TEMP" ) );
        aParameterLists( 3 ).set( "follower_dof_dependencies", std::string( "TEMP" ) );
        aParameterLists( 3 ).set( "leader_constitutive_models", std::string( "CMDiffusion_1,Diffusion" ) );
        aParameterLists( 3 ).set( "follower_constitutive_models", std::string( "CMDiffusion_3,Diffusion" ) );
        aParameterLists( 3 ).set( "stabilization_parameters", std::string( "SPInterfaceNitsche_12     ,NitscheInterface" ) );
        aParameterLists( 3 ).set( "mesh_set_names", tInterface_12 );

        // Temperature - Shell - Shell
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", std::string( "IWGInterfaceShellShellTEMP_3" ) );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_INTERFACE_UNSYMMETRIC_NITSCHE ) ;
        aParameterLists( 3 ).set( "dof_residual", std::string( "TEMP" ) );
        aParameterLists( 3 ).set( "leader_dof_dependencies", std::string( "TEMP" ) );
        aParameterLists( 3 ).set( "follower_dof_dependencies", std::string( "TEMP" ) );
        aParameterLists( 3 ).set( "leader_constitutive_models", std::string( "CMDiffusion_3,Diffusion" ) );
        aParameterLists( 3 ).set( "follower_constitutive_models", std::string( "CMDiffusion_4,Diffusion" ) );
        aParameterLists( 3 ).set( "stabilization_parameters", std::string( "SPInterfaceNitsche_22     ,NitscheInterface" ) );
        aParameterLists( 3 ).set( "mesh_set_names", tInterface_22 );

        // Temperature - Shell - Shell
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", std::string( "IWGInterfaceShellShellTEMP_4" ) );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_INTERFACE_UNSYMMETRIC_NITSCHE ) ;
        aParameterLists( 3 ).set( "dof_residual", std::string( "TEMP" ) );
        aParameterLists( 3 ).set( "leader_dof_dependencies", std::string( "TEMP" ) );
        aParameterLists( 3 ).set( "follower_dof_dependencies", std::string( "TEMP" ) );
        aParameterLists( 3 ).set( "leader_constitutive_models", std::string( "CMDiffusion_3,Diffusion" ) );
        aParameterLists( 3 ).set( "follower_constitutive_models", std::string( "CMDiffusion_1,Diffusion" ) );
        aParameterLists( 3 ).set( "stabilization_parameters", std::string( "SPInterfaceNitsche_21     ,NitscheInterface" ) );
        aParameterLists( 3 ).set( "mesh_set_names", tInterface_21 );

        if ( isGhost )
        {
            // create parameter list for IWG 16
            aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists( 3 ).set( "IWG_name", "GPenTemp_1" );
            aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
            aParameterLists( 3 ).set( "dof_residual", "TEMP" );
            aParameterLists( 3 ).set( "leader_dof_dependencies", "TEMP" );
            aParameterLists( 3 ).set( "follower_dof_dependencies", "TEMP" );
            aParameterLists( 3 ).set( "stabilization_parameters", "SPGPTemp_1,GhostSP" );
            aParameterLists( 3 ).set( "mesh_set_names", tGhost_1 );

            aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists( 3 ).set( "IWG_name", "GPenTemp_2" );
            aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
            aParameterLists( 3 ).set( "dof_residual", "TEMP" );
            aParameterLists( 3 ).set( "leader_dof_dependencies", "TEMP" );
            aParameterLists( 3 ).set( "follower_dof_dependencies", "TEMP" );
            aParameterLists( 3 ).set( "stabilization_parameters", "SPGPTemp_2,GhostSP" );
            aParameterLists( 3 ).set( "mesh_set_names", tGhost_2 );
            }

        //------------------------------------------------------------------------------
        // fill the IQI part of the parameter list

        // create parameter list for IQI 3
        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkTEMP" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists( 4 ).set( "dof_quantity", "TEMP" );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 4 ).set( "mesh_set_names", tBulk );

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        aParameterLists( 5 ).add_parameter_list( prm::create_computation_parameter_list() );
    }

    void
    SOLParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL ) );
#ifdef MORIS_USE_MUMPS
        aParameterLists( 0 ).set( "Solver_Type", "Amesos_Mumps" );
#else
        aParameterLists( 0 ).set( "Solver_Type", "Amesos_Superludist" );
#endif

        aParameterLists( 1 ).add_parameter_list( moris::prm::create_linear_solver_parameter_list() );

        aParameterLists( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
        aParameterLists( 2 ).set( "NLA_rel_res_norm_drop", 1e-04 );
        aParameterLists( 2 ).set( "NLA_relaxation_parameter", 1.0 );
        aParameterLists( 2 ).set( "NLA_max_iter", 1 );

        aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
        aParameterLists( 3 ).set( "NLA_DofTypes", "TEMP" );

        aParameterLists( 4 ).add_parameter_list( moris::prm::create_time_solver_algorithm_parameter_list() );

        aParameterLists( 5 ).add_parameter_list( moris::prm::create_time_solver_parameter_list() );
        aParameterLists( 5 ).set( "TSA_DofTypes", "TEMP" );
        aParameterLists( 5 ).set( "TSA_Initialize_Sol_Vec", "TEMP,0.0" );
        aParameterLists( 5 ).set( "TSA_Output_Indices", "0" );
        aParameterLists( 5 ).set( "TSA_Output_Criteria", "Output_Criterion" );

        aParameterLists( 6 ).add_parameter_list( moris::prm::create_solver_warehouse_parameterlist() );

        aParameterLists( 7 ).add_parameter_list( moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::NONE ) );
    }

    void
    MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_msi_parameter_list() );
        aParameterLists( 0 ).set( "order_adofs_by_host", false );
        aParameterLists( 0 ).set( "multigrid", isMultigrid );
    }

    void
    VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_vis_parameter_list() );
        aParameterLists( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", "Voxel_Grain_Structure.exo" ) );
        aParameterLists( 0 ).set( "Mesh_Type",  vis::VIS_Mesh_Type::STANDARD ) ;
        aParameterLists( 0 ).set( "Set_Names", tBulk );
        aParameterLists( 0 ).set( "Field_Names", "TEMP" );
        aParameterLists( 0 ).set( "Field_Type", "NODAL" );
        aParameterLists( 0 ).set( "IQI_Names", "IQIBulkTEMP" );
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
